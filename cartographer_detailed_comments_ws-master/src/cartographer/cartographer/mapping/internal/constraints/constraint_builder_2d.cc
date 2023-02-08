/*
 * Copyright 2016 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "cartographer/mapping/internal/constraints/constraint_builder_2d.h"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#include "Eigen/Eigenvalues"
#include "absl/memory/memory.h"
#include "cartographer/common/math.h"
#include "cartographer/common/thread_pool.h"
#include "cartographer/mapping/proto/scan_matching/ceres_scan_matcher_options_2d.pb.h"
#include "cartographer/mapping/proto/scan_matching/fast_correlative_scan_matcher_options_2d.pb.h"
#include "cartographer/metrics/counter.h"
#include "cartographer/metrics/gauge.h"
#include "cartographer/metrics/histogram.h"
#include "cartographer/transform/transform.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {
namespace constraints {

static auto* kConstraintsSearchedMetric = metrics::Counter::Null();
static auto* kConstraintsFoundMetric = metrics::Counter::Null();
static auto* kGlobalConstraintsSearchedMetric = metrics::Counter::Null();
static auto* kGlobalConstraintsFoundMetric = metrics::Counter::Null();
static auto* kQueueLengthMetric = metrics::Gauge::Null();
static auto* kConstraintScoresMetric = metrics::Histogram::Null();
static auto* kGlobalConstraintScoresMetric = metrics::Histogram::Null();
static auto* kNumSubmapScanMatchersMetric = metrics::Gauge::Null();

// 返回submap的原点在local坐标系下的二维坐标
transform::Rigid2d ComputeSubmapPose(const Submap2D& submap) {      //logic:由pose_graph.cc 56行调用
  return transform::Project2D(submap.local_pose());
}

/**
 * @brief 构造函数
 * 
 * @param[in] options 约束构造器的配置参数
 * @param[in] thread_pool map_builder中构造的线程池
 */
ConstraintBuilder2D::ConstraintBuilder2D(                                  //logic:由pose_graph_2d.cc 65 行调用
    const constraints::proto::ConstraintBuilderOptions& options,
    common::ThreadPoolInterface* const thread_pool)
    : options_(options),
      thread_pool_(thread_pool),
      finish_node_task_(absl::make_unique<common::Task>()),
      when_done_task_(absl::make_unique<common::Task>()),
      ceres_scan_matcher_(options.ceres_scan_matcher_options()) {}

ConstraintBuilder2D::~ConstraintBuilder2D() {
  absl::MutexLock locker(&mutex_);
  CHECK_EQ(finish_node_task_->GetState(), common::Task::NEW);
  CHECK_EQ(when_done_task_->GetState(), common::Task::NEW);
  CHECK_EQ(constraints_.size(), 0) << "WhenDone() was not called";
  CHECK_EQ(num_started_nodes_, num_finished_nodes_);
  CHECK(when_done_ == nullptr);
}

/**
 * @brief 进行局部搜索窗口的约束计算(对局部子图进行回环检测)
 * 
 * @param[in] submap_id submap的id
 * @param[in] submap 单个submap
 * @param[in] node_id 节点的id
 * @param[in] constant_data 节点的数据
 * @param[in] initial_relative_pose 约束的初值
 */
void ConstraintBuilder2D::MaybeAddConstraint(             //logic:由pose_graph_2d.cc 427行  调用
    const SubmapId& submap_id, const Submap2D* const submap, 
    const NodeId& node_id, const TrajectoryNode::Data* const constant_data,
    const transform::Rigid2d& initial_relative_pose) {
  // 超过范围的不进行约束的计算
  if (initial_relative_pose.translation().norm() >//jc:max_constraint_distance 超过这个范围就不进行回环检测，因为要和所有已经完成的子图进行匹配，所以计算量大
      options_.max_constraint_distance()) { // param: max_constraint_distance  //jc:pose_graph.lua 22 行最大的检测距离 单位 米
    return;
  }
  // 根据参数配置添加约束的频率
  if (!per_submap_sampler_
           .emplace(std::piecewise_construct, std::forward_as_tuple(submap_id),
                    std::forward_as_tuple(options_.sampling_ratio()))
           .first->second.Pulse()) {
    return;
  }

  absl::MutexLock locker(&mutex_);
  // 当when_done_正在处理任务时调用本函数, 报个警告
  if (when_done_) {
    LOG(WARNING)
        << "MaybeAddConstraint was called while WhenDone was scheduled.";
  }

  // 在队列中新建一个指向Constraint数据的指针
  constraints_.emplace_back();
  kQueueLengthMetric->Set(constraints_.size());
  auto* const constraint = &constraints_.back();
  
  // 为子图新建一个匹配器
  const auto* scan_matcher =
      DispatchScanMatcherConstruction(submap_id, submap->grid());   //logic:调用本文件226行  ，里面构建多分辨率地图

  // 生成个计算约束的任务
  auto constraint_task = absl::make_unique<common::Task>();
  constraint_task->SetWorkItem([=]() LOCKS_EXCLUDED(mutex_) {
    ComputeConstraint(submap_id, submap, node_id, false, /* match_full_submap */     //logic:调用本文件269行
                      constant_data, initial_relative_pose, *scan_matcher,
                      constraint);
  });

  // 等匹配器之后初始化才能进行约束的计算
   //jc:添加本任务依赖的任务scan_matcher->creation_task_handle
  constraint_task->AddDependency(scan_matcher->creation_task_handle);//jc:调用task.cc89行   //jc:这里要使用多分辨率地图必须等多分辨率地图建完，
  // 将计算约束这个任务放入线程池等待执行 //jc:scan_matcher->creation_task_handle就是构造多分辨率地图
  auto constraint_task_handle =
      thread_pool_->Schedule(std::move(constraint_task));   //jc:Schedule直接加入线程池
  // 将计算约束这个任务 添加到 finish_node_task_的依赖项中
  finish_node_task_->AddDependency(constraint_task_handle);  //jc:AddDependency先添加依赖再进入线程池
}

/**
 * @brief 进行全局搜索窗口的约束计算(对整体子图进行回环检测)
 * 
 * @param[in] submap_id submap的id
 * @param[in] submap 单个submap
 * @param[in] node_id 节点的id
 * @param[in] constant_data 节点的数据
 */
void ConstraintBuilder2D::MaybeAddGlobalConstraint(              //logic:由pose_graph_2d.cc 152行 以及433行 调用    
    const SubmapId& submap_id, const Submap2D* const submap,
    const NodeId& node_id, const TrajectoryNode::Data* const constant_data) {
  absl::MutexLock locker(&mutex_);
  if (when_done_) {
    LOG(WARNING)
        << "MaybeAddGlobalConstraint was called while WhenDone was scheduled.";
  }
  
  // note: 对整体子图进行回环检测时没有距离的限制

  constraints_.emplace_back();
  kQueueLengthMetric->Set(constraints_.size());
  auto* const constraint = &constraints_.back();
  // 为子图新建一个匹配器
  const auto* scan_matcher =
      DispatchScanMatcherConstruction(submap_id, submap->grid()); //jc:生成多分辨率地图 ，调用本文件226行
  auto constraint_task = absl::make_unique<common::Task>();
  // 生成个计算全局约束的任务
  constraint_task->SetWorkItem([=]() LOCKS_EXCLUDED(mutex_) {
    ComputeConstraint(submap_id, submap, node_id, true, /* match_full_submap */
                      constant_data, transform::Rigid2d::Identity(),
                      *scan_matcher, constraint);
  });
  constraint_task->AddDependency(scan_matcher->creation_task_handle);
  auto constraint_task_handle =
      thread_pool_->Schedule(std::move(constraint_task));
  finish_node_task_->AddDependency(constraint_task_handle);
}

// 告诉ConstraintBuilder2D的对象, 刚刚完成了一个节点的约束的计算
void ConstraintBuilder2D::NotifyEndOfNode() {        //logic:由pose_graph_2d.cc 561行 调用
  absl::MutexLock locker(&mutex_);
  CHECK(finish_node_task_ != nullptr);
  
  // 生成个任务: 将num_finished_nodes_自加, 记录完成约束计算节点的总个数 //jc:记录完成约束计算,计算这个节点约束的总个数ComputeConstraint
  finish_node_task_->SetWorkItem([this] {
    absl::MutexLock locker(&mutex_);
    ++num_finished_nodes_;                   //jc:记录完成约束计算节点的总个数
  });

  // 将这个任务传入线程池中等待执行, 由于之前添加了依赖, 所以finish_node_task_一定会比计算约束更晚完成
  auto finish_node_task_handle =
      thread_pool_->Schedule(std::move(finish_node_task_));

  // move之后finish_node_task_就没有指向的地址了, 所以这里要重新初始化
  finish_node_task_ = absl::make_unique<common::Task>();
  // 设置when_done_task_依赖finish_node_task_handle
  when_done_task_->AddDependency(finish_node_task_handle); //jc:when_done_task_依赖于finish_node_task_handle
  ++num_started_nodes_;
}

// 约束计算完成之后执行一下回调函数
void ConstraintBuilder2D::WhenDone(
    const std::function<void(const ConstraintBuilder2D::Result&)>& callback) {
  absl::MutexLock locker(&mutex_);
  CHECK(when_done_ == nullptr);

  // TODO(gaschler): Consider using just std::function, it can also be empty.
  // 将回调函数赋值给when_done_
  when_done_ = absl::make_unique<std::function<void(const Result&)>>(callback);
  CHECK(when_done_task_ != nullptr);

  // 生成 执行when_done_的任务
  when_done_task_->SetWorkItem([this] { RunWhenDoneCallback(); });
  // 将任务放入线程池中等待执行
  thread_pool_->Schedule(std::move(when_done_task_));

  // when_done_task_的重新初始化
  when_done_task_ = absl::make_unique<common::Task>();
}

// 为每个子图新建一个匹配器
const ConstraintBuilder2D::SubmapScanMatcher*
ConstraintBuilder2D::DispatchScanMatcherConstruction(const SubmapId& submap_id,   //logic:由本文件168行调用
                                                     const Grid2D* const grid) {
  CHECK(grid);
  // 如果匹配器里已经存在, 则直接返回对应id的匹配器
  if (submap_scan_matchers_.count(submap_id) != 0) {
    return &submap_scan_matchers_.at(submap_id);
  }
  // submap_scan_matchers_新增加一个 key
  auto& submap_scan_matcher = submap_scan_matchers_[submap_id];
  kNumSubmapScanMatchersMetric->Set(submap_scan_matchers_.size());
  // 保存栅格地图的指针
  submap_scan_matcher.grid = grid;

  auto& scan_matcher_options = options_.fast_correlative_scan_matcher_options();
  auto scan_matcher_task = absl::make_unique<common::Task>();
  // 生成一个将初始化匹配器的任务, 初始化时会计算多分辨率地图, 比较耗时
  scan_matcher_task->SetWorkItem(                                             //logic:调用task.cc  37行
      [&submap_scan_matcher, &scan_matcher_options]() {
        // 进行匹配器的初始化, 与多分辨率地图的创建
        submap_scan_matcher.fast_correlative_scan_matcher =
            absl::make_unique<scan_matching::FastCorrelativeScanMatcher2D>(  //logic:调用fast_correlative_scan_matcher_2d.cc  273 行
                *submap_scan_matcher.grid, scan_matcher_options);  //jc:这里一个地图就变成7张地图了，生成了多分辨率地图，其实栅格数一样，只不过分辨率低的多个栅格代表一个
      });                                                            //jc:*submap_scan_matcher.grid 为栅格地图
  // 将初始化匹配器的任务放入线程池中, 并且将任务的智能指针保存起来
  submap_scan_matcher.creation_task_handle =
      thread_pool_->Schedule(std::move(scan_matcher_task));

  return &submap_scan_matchers_.at(submap_id);
}

/**
 * @brief 计算节点和子图之间的一个约束(回环检测)
 *        用基于分支定界算法的匹配器进行粗匹配,然后用ceres进行精匹配
 * 
 * @param[in] submap_id submap的id
 * @param[in] submap 地图数据
 * @param[in] node_id 节点id
 * @param[in] match_full_submap 是局部匹配还是全子图匹配
 * @param[in] constant_data 节点数据
 * @param[in] initial_relative_pose 约束的初值
 * @param[in] submap_scan_matcher 匹配器
 * @param[out] constraint 计算出的约束
 */
void ConstraintBuilder2D::ComputeConstraint(                              //logic:有本文件129 行调用 pose_graph_2d.cc的ComputeConstraint调用MaybeAddConstraint，MaybeAddConstraint在调用ComputeConstraint
    const SubmapId& submap_id, const Submap2D* const submap,
    const NodeId& node_id, bool match_full_submap,                     //jc:match_full_submap是否匹配全子图
    const TrajectoryNode::Data* const constant_data,
    const transform::Rigid2d& initial_relative_pose,
    const SubmapScanMatcher& submap_scan_matcher,
    std::unique_ptr<ConstraintBuilder2D::Constraint>* constraint) {
  CHECK(submap_scan_matcher.fast_correlative_scan_matcher);

  // Step:1 得到节点在local frame下的坐标
  const transform::Rigid2d initial_pose =                      //jc:tracking_frame就是节点
      ComputeSubmapPose(*submap) * initial_relative_pose;  //jc:initial_relative_pose在pose_graph_2d.cc 422行，代表子图间相对变换，ComputeSubmapPose(*submap)为子图在local坐标系下的坐标

  // The 'constraint_transform' (submap i <- node j) is computed from:
  // - a 'filtered_gravity_aligned_point_cloud' in node j,
  // - the initial guess 'initial_pose' for (map <- node j),
  // - the result 'pose_estimate' of Match() (map <- node j).
  // - the ComputeSubmapPose() (map <- submap i)

  float score = 0.;
  transform::Rigid2d pose_estimate = transform::Rigid2d::Identity();

  // Compute 'pose_estimate' in three stages:
  // 1. Fast estimate using the fast correlative scan matcher.
  // 2. Prune if the score is too low.
  // 3. Refine.
  // param: global_localization_min_score 对整体子图进行回环检测时的最低分数阈值
  // param: min_score 对局部子图进行回环检测时的最低分数阈值

  // Step:2 使用基于分支定界算法的匹配器进行粗匹配
  if (match_full_submap) {
    // 节点与全地图进行匹配
    kGlobalConstraintsSearchedMetric->Increment();
    if (submap_scan_matcher.fast_correlative_scan_matcher->MatchFullSubmap(   //jc:基于分支定界算法的粗匹配
            constant_data->filtered_gravity_aligned_point_cloud,
            options_.global_localization_min_score(), &score, &pose_estimate)) {
      CHECK_GT(score, options_.global_localization_min_score());
      CHECK_GE(node_id.trajectory_id, 0);
      CHECK_GE(submap_id.trajectory_id, 0);
      kGlobalConstraintsFoundMetric->Increment();
      kGlobalConstraintScoresMetric->Observe(score);
    } else {
      // 计算失败了就退出
      return;
    }
  } 
  else {
    // 节点与局部地图进行匹配
    kConstraintsSearchedMetric->Increment();
    if (submap_scan_matcher.fast_correlative_scan_matcher->Match( //jc:基于分支定界算法的匹配  调用fast_correlative_scan_matcher_2d.cc 294 行
            initial_pose, constant_data->filtered_gravity_aligned_point_cloud,
            options_.min_score(), &score, &pose_estimate)) {
      // We've reported a successful local match.
      CHECK_GT(score, options_.min_score());
      kConstraintsFoundMetric->Increment();
      kConstraintScoresMetric->Observe(score);
    } else {
      return;
    }
  }
  
  {
    absl::MutexLock locker(&mutex_);
    score_histogram_.Add(score);
  }

  // Use the CSM estimate as both the initial and previous pose. This has the
  // effect that, in the absence of better information, we prefer the original
  // CSM estimate.

  // Step:3 使用ceres进行精匹配, 就是前端扫描匹配使用的函数
  ceres::Solver::Summary unused_summary;
  ceres_scan_matcher_.Match(pose_estimate.translation(), pose_estimate,                      //jc:基于ceres算法的精匹配
                            constant_data->filtered_gravity_aligned_point_cloud,
                            *submap_scan_matcher.grid, &pose_estimate,
                            &unused_summary);

  //jc:这里是计算子图间约束
  // Step:4 获取节点到submap坐标系原点间的坐标变换 //jc:子图原点到tracking_frame(tracking_frame为carto跟踪的坐标系)的坐标变换；节点就是指tracking_frame
  // pose_estimate 是 节点在 loacl frame 下的坐标 //jc:子图内约束和子图间都是子图原点到节点（节点==pose_estimate==tracking_frame）的坐标变换
  const transform::Rigid2d constraint_transform =     //jc:子图原点到匹配之后的local坐标系下的pose_estimate，得到相对坐标变换
      ComputeSubmapPose(*submap).inverse() * pose_estimate; //jc:子图在local坐标系下的inverse * pose_estimate在local坐标系的坐标得到子图到节点在Local坐标系下的坐标变换，即子图间约束，看箭头图

  // Step:5 返回计算后的约束
  constraint->reset(new Constraint{submap_id,
                                   node_id,
                                   {transform::Embed3D(constraint_transform),
                                    options_.loop_closure_translation_weight(),   //jc:平移旋转的权重，不用改
                                    options_.loop_closure_rotation_weight()},      //jc:平移旋转的权重，不用改
                                   Constraint::INTER_SUBMAP});

  // log相关
  if (options_.log_matches()) {
    std::ostringstream info;
    info << "Node " << node_id << " with "
         << constant_data->filtered_gravity_aligned_point_cloud.size()
         << " points on submap " << submap_id << std::fixed;
    if (match_full_submap) {
      info << " matches";
    } else {
      const transform::Rigid2d difference =
          initial_pose.inverse() * pose_estimate;
      info << " differs by translation " << std::setprecision(2) // c++11: std::setprecision(2) 保留2个小数点
           << difference.translation().norm() << " rotation "
           << std::setprecision(3) << std::abs(difference.normalized_angle());
    }
    info << " with score " << std::setprecision(1) << 100. * score << "%.";
    LOG(INFO) << info.str();
  }
}

// 将临时保存的所有约束数据传入回调函数, 并执行回调函数
void ConstraintBuilder2D::RunWhenDoneCallback() {
  Result result;
  std::unique_ptr<std::function<void(const Result&)>> callback;
  {
    absl::MutexLock locker(&mutex_);
    CHECK(when_done_ != nullptr);  //jc:这个when_done在205行由call_back赋值

    // 将计算完的约束进行保存
    for (const std::unique_ptr<Constraint>& constraint : constraints_) {
      if (constraint == nullptr) continue;
      result.push_back(*constraint);
    }

    if (options_.log_matches()) {
      LOG(INFO) << constraints_.size() << " computations resulted in "
                << result.size() << " additional constraints.";
      LOG(INFO) << "Score histogram:\n" << score_histogram_.ToString(10);
    }

    // 这些约束已经保存过了, 就可以删掉了
    constraints_.clear();

    callback = std::move(when_done_);
    when_done_.reset();
    kQueueLengthMetric->Set(constraints_.size());
  }
  // 执行回调函数 HandleWorkQueue
  (*callback)(result); //jc:这里的callback就是when_done，所以就是传入的callback也就是handleworkqueue 
}

// 获取完成约束计算节点的总个数
int ConstraintBuilder2D::GetNumFinishedNodes() {
  absl::MutexLock locker(&mutex_);
  return num_finished_nodes_;
}

// 删除指定submap_id的匹配器
void ConstraintBuilder2D::DeleteScanMatcher(const SubmapId& submap_id) {
  absl::MutexLock locker(&mutex_);
  if (when_done_) {
    LOG(WARNING)
        << "DeleteScanMatcher was called while WhenDone was scheduled.";
  }
  submap_scan_matchers_.erase(submap_id);
  per_submap_sampler_.erase(submap_id);
  kNumSubmapScanMatchersMetric->Set(submap_scan_matchers_.size());
}

void ConstraintBuilder2D::RegisterMetrics(metrics::FamilyFactory* factory) {
  auto* counts = factory->NewCounterFamily(
      "mapping_constraints_constraint_builder_2d_constraints",
      "Constraints computed");
  kConstraintsSearchedMetric =
      counts->Add({{"search_region", "local"}, {"matcher", "searched"}});
  kConstraintsFoundMetric =
      counts->Add({{"search_region", "local"}, {"matcher", "found"}});
  kGlobalConstraintsSearchedMetric =
      counts->Add({{"search_region", "global"}, {"matcher", "searched"}});
  kGlobalConstraintsFoundMetric =
      counts->Add({{"search_region", "global"}, {"matcher", "found"}});
  auto* queue_length = factory->NewGaugeFamily(
      "mapping_constraints_constraint_builder_2d_queue_length", "Queue length");
  kQueueLengthMetric = queue_length->Add({});
  auto boundaries = metrics::Histogram::FixedWidth(0.05, 20);
  auto* scores = factory->NewHistogramFamily(
      "mapping_constraints_constraint_builder_2d_scores",
      "Constraint scores built", boundaries);
  kConstraintScoresMetric = scores->Add({{"search_region", "local"}});
  kGlobalConstraintScoresMetric = scores->Add({{"search_region", "global"}});
  auto* num_matchers = factory->NewGaugeFamily(
      "mapping_constraints_constraint_builder_2d_num_submap_scan_matchers",
      "Current number of constructed submap scan matchers");
  kNumSubmapScanMatchersMetric = num_matchers->Add({});
}

}  // namespace constraints
}  // namespace mapping
}  // namespace cartographer
