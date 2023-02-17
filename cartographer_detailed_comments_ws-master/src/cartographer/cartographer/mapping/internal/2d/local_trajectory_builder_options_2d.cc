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

#include "cartographer/mapping/internal/2d/local_trajectory_builder_options_2d.h"

#include "cartographer/mapping/2d/submap_2d.h"
#include "cartographer/mapping/internal/2d/scan_matching/ceres_scan_matcher_2d.h"
#include "cartographer/mapping/internal/motion_filter.h"
#include "cartographer/mapping/internal/scan_matching/real_time_correlative_scan_matcher.h"
#include "cartographer/mapping/pose_extrapolator_interface.h"
#include "cartographer/sensor/internal/voxel_filter.h"

namespace cartographer {
namespace mapping {

proto::LocalTrajectoryBuilderOptions2D CreateLocalTrajectoryBuilderOptions2D( //logic:由trajectory_builder_insterface 59 行调用
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::LocalTrajectoryBuilderOptions2D options;
  options.set_min_range(parameter_dictionary->GetDouble("min_range"));
  options.set_max_range(parameter_dictionary->GetDouble("max_range"));
  options.set_min_z(parameter_dictionary->GetDouble("min_z"));
  options.set_max_z(parameter_dictionary->GetDouble("max_z"));
  options.set_missing_data_ray_length(
      parameter_dictionary->GetDouble("missing_data_ray_length"));
  options.set_num_accumulated_range_data(
      parameter_dictionary->GetInt("num_accumulated_range_data"));
  options.set_voxel_filter_size(
      parameter_dictionary->GetDouble("voxel_filter_size"));
  options.set_use_online_correlative_scan_matching(
      parameter_dictionary->GetBool("use_online_correlative_scan_matching"));
  *options.mutable_adaptive_voxel_filter_options() =
      sensor::CreateAdaptiveVoxelFilterOptions(
          parameter_dictionary->GetDictionary("adaptive_voxel_filter").get());
  *options.mutable_loop_closure_adaptive_voxel_filter_options() =
      sensor::CreateAdaptiveVoxelFilterOptions(
          parameter_dictionary
              ->GetDictionary("loop_closure_adaptive_voxel_filter")
              .get());
  *options.mutable_real_time_correlative_scan_matcher_options() =
      mapping::scan_matching::CreateRealTimeCorrelativeScanMatcherOptions(
          parameter_dictionary
              ->GetDictionary("real_time_correlative_scan_matcher")
              .get());
  *options.mutable_ceres_scan_matcher_options() =
      mapping::scan_matching::CreateCeresScanMatcherOptions2D(
          parameter_dictionary->GetDictionary("ceres_scan_matcher").get());
  *options.mutable_motion_filter_options() = mapping::CreateMotionFilterOptions(
      parameter_dictionary->GetDictionary("motion_filter").get());
  *options.mutable_pose_extrapolator_options() = CreatePoseExtrapolatorOptions( //jc:这里调用了位子推测器的函数pose_extrapolator_interface.cc 68行
      parameter_dictionary->GetDictionary("pose_extrapolator").get());   //jc:pose_extrapolator这个字段根据lua文件中中括号前面的字段确定
  options.set_imu_gravity_time_constant(
      parameter_dictionary->GetDouble("imu_gravity_time_constant"));
  *options.mutable_submaps_options() = CreateSubmapsOptions2D(
      parameter_dictionary->GetDictionary("submaps").get());
  options.set_use_imu_data(parameter_dictionary->GetBool("use_imu_data"));
  return options;
}

}  // namespace mapping
}  // namespace cartographer
