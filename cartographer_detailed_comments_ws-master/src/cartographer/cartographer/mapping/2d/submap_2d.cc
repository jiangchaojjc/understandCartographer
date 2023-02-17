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

#include "cartographer/mapping/2d/submap_2d.h"

#include <cinttypes>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>

#include "Eigen/Geometry"
#include "absl/memory/memory.h"
#include "cartographer/common/port.h"
#include "cartographer/mapping/2d/probability_grid_range_data_inserter_2d.h"
#include "cartographer/mapping/internal/2d/tsdf_range_data_inserter_2d.h"
#include "cartographer/mapping/range_data_inserter_interface.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

// submap2d的参数配置配置
proto::SubmapsOptions2D CreateSubmapsOptions2D(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::SubmapsOptions2D options;
  options.set_num_range_data(
      parameter_dictionary->GetNonNegativeInt("num_range_data"));
  *options.mutable_grid_options_2d() = CreateGridOptions2D(
      parameter_dictionary->GetDictionary("grid_options_2d").get());
  *options.mutable_range_data_inserter_options() =
      CreateRangeDataInserterOptions(
          parameter_dictionary->GetDictionary("range_data_inserter").get());

  bool valid_range_data_inserter_grid_combination = false;
  // 地图类型
  const proto::GridOptions2D_GridType& grid_type =
      options.grid_options_2d().grid_type();
  // 将scan写成地图的方式
  const proto::RangeDataInserterOptions_RangeDataInserterType&
      range_data_inserter_type =
          options.range_data_inserter_options().range_data_inserter_type();
  if (grid_type == proto::GridOptions2D::PROBABILITY_GRID &&
      range_data_inserter_type ==
          proto::RangeDataInserterOptions::PROBABILITY_GRID_INSERTER_2D) {
    valid_range_data_inserter_grid_combination = true;
  }
  if (grid_type == proto::GridOptions2D::TSDF &&
      range_data_inserter_type ==
          proto::RangeDataInserterOptions::TSDF_INSERTER_2D) {
    valid_range_data_inserter_grid_combination = true;
  }
  CHECK(valid_range_data_inserter_grid_combination)
      << "Invalid combination grid_type " << grid_type
      << " with range_data_inserter_type " << range_data_inserter_type;
  CHECK_GT(options.num_range_data(), 0);
  return options;
}

/**
 * @brief 构造函数
 * 
 * @param[in] origin Submap2D的原点,保存在Submap类里
 * @param[in] grid 地图数据的指针
 * @param[in] conversion_tables 地图数据的转换表
 */
Submap2D::Submap2D(const Eigen::Vector2f& origin, std::unique_ptr<Grid2D> grid,
                   ValueConversionTables* conversion_tables)
    : Submap(transform::Rigid3d::Translation(
          Eigen::Vector3d(origin.x(), origin.y(), 0.))),
      conversion_tables_(conversion_tables) {
  grid_ = std::move(grid);
}

// 根据proto::Submap格式的数据生成Submap2D
Submap2D::Submap2D(const proto::Submap2D& proto,
                   ValueConversionTables* conversion_tables)
    : Submap(transform::ToRigid3(proto.local_pose())),
      conversion_tables_(conversion_tables) {
  if (proto.has_grid()) {
    if (proto.grid().has_probability_grid_2d()) {
      grid_ =
          absl::make_unique<ProbabilityGrid>(proto.grid(), conversion_tables_);
    } else if (proto.grid().has_tsdf_2d()) {
      grid_ = absl::make_unique<TSDF2D>(proto.grid(), conversion_tables_);
    } else {
      LOG(FATAL) << "proto::Submap2D has grid with unknown type.";
    }
  }
  set_num_range_data(proto.num_range_data());
  set_insertion_finished(proto.finished());
}

// 根据mapping::Submap2D生成proto::Submap格式的数据
proto::Submap Submap2D::ToProto(const bool include_grid_data) const {
  proto::Submap proto;
  auto* const submap_2d = proto.mutable_submap_2d();
  *submap_2d->mutable_local_pose() = transform::ToProto(local_pose());
  submap_2d->set_num_range_data(num_range_data());
  submap_2d->set_finished(insertion_finished());
  if (include_grid_data) {
    CHECK(grid_);
    *submap_2d->mutable_grid() = grid_->ToProto();
  }
  return proto;
}

// 根据proto::Submap格式的数据更新地图
void Submap2D::UpdateFromProto(const proto::Submap& proto) {
  CHECK(proto.has_submap_2d());
  const auto& submap_2d = proto.submap_2d();
  set_num_range_data(submap_2d.num_range_data());
  set_insertion_finished(submap_2d.finished());
  if (proto.submap_2d().has_grid()) {
    if (proto.submap_2d().grid().has_probability_grid_2d()) {
      grid_ = absl::make_unique<ProbabilityGrid>(proto.submap_2d().grid(),
                                                 conversion_tables_);
    } else if (proto.submap_2d().grid().has_tsdf_2d()) {
      grid_ = absl::make_unique<TSDF2D>(proto.submap_2d().grid(),
                                        conversion_tables_);
    } else {
      LOG(FATAL) << "proto::Submap2D has grid with unknown type.";
    }
  }
}

/**
 * @brief 将地图进行压缩, 放入response
 * 
 * @param[out] response 压缩后的地图数据
 */
void Submap2D::ToResponseProto(
    const transform::Rigid3d&,
    proto::SubmapQuery::Response* const response) const {
  if (!grid_) return;
  response->set_submap_version(num_range_data());
  // note: const在*后边, 指针指向的地址不能变,而内存单元中的内容可变
  proto::SubmapQuery::Response::SubmapTexture* const texture =
      response->add_textures();
  // 填充压缩后的数据
  grid()->DrawToSubmapTexture(texture, local_pose());
}

// 将雷达数据写到栅格地图中
void Submap2D::InsertRangeData(  //logic:由本文件206行调用
    const sensor::RangeData& range_data,
    const RangeDataInserterInterface* range_data_inserter) {
  CHECK(grid_);
  CHECK(!insertion_finished());
  // 将雷达数据写到栅格地图中                                 //jc:insert里面会更新地图栅格值
  range_data_inserter->Insert(range_data, grid_.get());   //jc:这个指针在activesubmap2s构造的时候 181行 根据CreateRangeDataInserter创建PROBABILITY_GRID_INSERTER_2D
  // 插入到地图中的雷达数据的个数加1                         //jc:所以range_data_inserter指的是ProbabilityGridRangeDataInserter2D的对象，
                                                          //logic:调用probability_grid_range_data_insert_2d.cc 152行
  set_num_range_data(num_range_data() + 1);
}

// 将子图标记为完成状态
void Submap2D::Finish() {
  CHECK(grid_);
  CHECK(!insertion_finished());
  grid_ = grid_->ComputeCroppedGrid(); //jc:在这里设置栅格地图概率情况
  // 将子图标记为完成状态
  set_insertion_finished(true);
}

/********** ActiveSubmaps2D *****************/

// ActiveSubmaps2D构造函数
ActiveSubmaps2D::ActiveSubmaps2D(const proto::SubmapsOptions2D& options)     //logic:由local_trajectory_builder_2d.cc48行调用并赋值
    : options_(options), range_data_inserter_(CreateRangeDataInserter()) {}   //logic:调用本文件212行

// 返回指向 Submap2D 的 shared_ptr指针 的vector
std::vector<std::shared_ptr<const Submap2D>> ActiveSubmaps2D::submaps() const {
  return std::vector<std::shared_ptr<const Submap2D>>(submaps_.begin(),
                                                      submaps_.end());
}

// 将点云数据写入到submap中 //jc:写入数据到地图的操作从这里开始。
//jc:初始时候，没有submap,调用addsubmap新建一个子图，然后再不断的调用insertrangedata,直到submaps_.back().num.. == 90 又addSubmap
//jc:又会新建一个子图，继续不断调用insertRangeData,直到直到 submaps_.front()->num_range_data()== 180，将第一个子图标记成完成状态，
//jc:下次再调用InsertRangeData,submaps_.back().num.. == 90,删掉第一个子图，新建一个新的子图。重复如此
     //jc:submaps_.front()->Finish(); 将第一个子图标记成完成状态
std::vector<std::shared_ptr<const Submap2D>> ActiveSubmaps2D::InsertRangeData( //logic:由local_trajectory_builder_2d.cc 399行调用/ 不断的调用这个函数
    const sensor::RangeData& range_data) {                                          
  // 如果第二个子图插入节点的数据等于num_range_data时,就新建个子图                 
  // 因为这时第一个子图应该已经处于完成状态了
  if (submaps_.empty() ||                                                            //jc: options_.num_range_data() = 90
      submaps_.back()->num_range_data() == options_.num_range_data()) {         //jc:当submaps_.front()->num_range_data() == 180 时 submaps_.back()->num_range_data() ==90
    AddSubmap(range_data.origin.head<2>());                                     //jc:将扫描匹配之后的数据的坐标原点取出来作为地图的原点（local坐标系下）  到submap_2d.cc的278行
  }                                                                                //jc:如果前端完成更新了，后端又不会对子图进行更改，所以出现叠图
  // 将一帧雷达数据同时写入两个子图中
  for (auto& submap : submaps_) {                                                     //jc:里面由submap1和submap2
    submap->InsertRangeData(range_data, range_data_inserter_.get());                        //jc: 插入一次数据之后num_range_data 就会+1
  }                                                                                         //logic:调用本文件158行
  // 第一个子图的节点数量等于2倍的num_range_data时,第二个子图节点数量应该等于num_range_data
  if (submaps_.front()->num_range_data() == 2 * options_.num_range_data()) {                       //jc: 直到 submaps_.front()->num_range_data()== 180
    submaps_.front()->Finish();                                                                //jc:将第一个子图标记成完成状态
  }
  return submaps();
}

// 创建地图数据写入器
std::unique_ptr<RangeDataInserterInterface>
ActiveSubmaps2D::CreateRangeDataInserter() {   //logic:由本文件182行  构造函数调用
  switch (options_.range_data_inserter_options().range_data_inserter_type()) {  //jc:根据trajectory_builder_2d.lua 104 行的变量PROBABILITY_GRID_INSERTER_2D决定
    // 概率栅格地图的写入器
    case proto::RangeDataInserterOptions::PROBABILITY_GRID_INSERTER_2D:
      return absl::make_unique<ProbabilityGridRangeDataInserter2D>(  //logic:调用probability_grid_range_data_inserter_2d.cc 136 行 ，创建这个对象
          options_.range_data_inserter_options()
              .probability_grid_range_data_inserter_options_2d());  
    // tsdf地图的写入器 
    case proto::RangeDataInserterOptions::TSDF_INSERTER_2D:
      return absl::make_unique<TSDFRangeDataInserter2D>(
          options_.range_data_inserter_options()
              .tsdf_range_data_inserter_options_2d());
    default:
      LOG(FATAL) << "Unknown RangeDataInserterType.";
  }
}

// 以当前雷达原点为地图原件创建地图
std::unique_ptr<GridInterface> ActiveSubmaps2D::CreateGrid(  //logic:由本文件282行调用
    const Eigen::Vector2f& origin) {
  // 地图初始大小,100个栅格
  constexpr int kInitialSubmapSize = 100;
  float resolution = options_.grid_options_2d().resolution(); // param: grid_options_2d.resolution
  switch (options_.grid_options_2d().grid_type()) {  //jc:trajectory_builder_2d.lua 100行的参数
    // 概率栅格地图
    case proto::GridOptions2D::PROBABILITY_GRID:
      return absl::make_unique<ProbabilityGrid>(   //jc:新建一个probabilityGrid的对象
          MapLimits(resolution,
                    // 左上角坐标为坐标系的最大值, origin位于地图的中间
                    origin.cast<double>() + 0.5 * kInitialSubmapSize *   //jc:这里计算完之后origin就位于地图的中间，左上角坐标最大，物理坐标
                                                resolution *
                                                Eigen::Vector2d::Ones(), //jc:物理坐标变为二维的向量
                    CellLimits(kInitialSubmapSize, kInitialSubmapSize)),
          &conversion_tables_);
    // tsdf地图
    case proto::GridOptions2D::TSDF:
      return absl::make_unique<TSDF2D>(
          MapLimits(resolution,
                    origin.cast<double>() + 0.5 * kInitialSubmapSize *
                                                resolution *
                                                Eigen::Vector2d::Ones(),
                    CellLimits(kInitialSubmapSize, kInitialSubmapSize)),
          options_.range_data_inserter_options()
              .tsdf_range_data_inserter_options_2d()
              .truncation_distance(),               // 0.3
          options_.range_data_inserter_options()
              .tsdf_range_data_inserter_options_2d()
              .maximum_weight(),                    // 10.0
          &conversion_tables_);
    default:
      LOG(FATAL) << "Unknown GridType.";
  }
}

// 新增一个子图,根据子图个数判断是否删掉第一个子图
void ActiveSubmaps2D::AddSubmap(const Eigen::Vector2f& origin) {  
  // 调用AddSubmap时第一个子图一定是完成状态,所以子图数为2时就可以删掉第一个子图了
  if (submaps_.size() >= 2) {   //jc:当->num_range_data()第一次==180时，已经加入了两个submap，此时删掉第一个指针，
                            //jc:新建submap3的子图，submap3的num_range_data = 1,submap2的num_range_data =91，不断循环
    // This will crop the finished Submap before inserting a new Submap to
    // reduce peak memory usage a bit.
    CHECK(submaps_.front()->insertion_finished());
    // 删掉第一个子图的指针
    submaps_.erase(submaps_.begin());  //jc:删掉submap第一个shared_ptr的指针
  }
  // 新建一个子图, 并保存指向新子图的智能指针
  submaps_.push_back(absl::make_unique<Submap2D>(
      origin,                                         //jc:整个前端的origin都是在local坐标系下，也就是track_frame
      std::unique_ptr<Grid2D>(
          static_cast<Grid2D*>(CreateGrid(origin).release())),
      &conversion_tables_));
}

}  // namespace mapping
}  // namespace cartographer
