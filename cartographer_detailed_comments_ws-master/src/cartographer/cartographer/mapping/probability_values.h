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

#ifndef CARTOGRAPHER_MAPPING_PROBABILITY_VALUES_H_
#define CARTOGRAPHER_MAPPING_PROBABILITY_VALUES_H_

#include <cmath>
#include <vector>

#include "cartographer/common/math.h"
#include "cartographer/common/port.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

namespace {
// 为了避免浮点运算, 将[0~0.9]的浮点数转成[1~32767]之间的值
inline uint16 BoundedFloatToValue(const float float_value,
                                  const float lower_bound,
                                  const float upper_bound) {
  const int value =
      common::RoundToInt(
          (common::Clamp(float_value, lower_bound, upper_bound) - lower_bound) *
          (32766.f / (upper_bound - lower_bound))) +
      1;
  // DCHECK for performance.
  DCHECK_GE(value, 1);
  DCHECK_LE(value, 32767);
  return value;
}

}  // namespace

// 通过概率计算Odd值 论文里的 odds(p)函数
inline float Odds(float probability) {   //jc:将probability转换到odd
  return probability / (1.f - probability);
}

// 通过Odd值计算概率值 论文里的 odds^-1 函数
inline float ProbabilityFromOdds(const float odds) { // //jc:将odd转换到probability
  return odds / (odds + 1.f);
}

// probability与CorrespondenceCost的关系, CorrespondenceCost代表free的概率  //jc:probability代表占用的概率，CorrespondenceCost代表free的概率
inline float ProbabilityToCorrespondenceCost(const float probability) {
  return 1.f - probability;
}

inline float CorrespondenceCostToProbability(const float correspondence_cost) {
  return 1.f - correspondence_cost;
}

constexpr float kMinProbability = 0.1f;                         // 0.1
constexpr float kMaxProbability = 1.f - kMinProbability;        // 0.9
constexpr float kMinCorrespondenceCost = 1.f - kMaxProbability; // 0.1
constexpr float kMaxCorrespondenceCost = 1.f - kMinProbability; // 0.9

// Clamps probability to be in the range [kMinProbability, kMaxProbability].
// 对数据进行上下界的限定
inline float ClampProbability(const float probability) {
  return common::Clamp(probability, kMinProbability, kMaxProbability);
}

// Clamps correspondece cost to be in the range [kMinCorrespondenceCost,
// kMaxCorrespondenceCost].
// 对数据进行上下界的限定
inline float ClampCorrespondenceCost(const float correspondence_cost) {
  return common::Clamp(correspondence_cost, kMinCorrespondenceCost,
                       kMaxCorrespondenceCost);
}

constexpr uint16 kUnknownProbabilityValue = 0; // 0
constexpr uint16 kUnknownCorrespondenceValue = kUnknownProbabilityValue; // 0
constexpr uint16 kUpdateMarker = 1u << 15; // 32768

// Converts a correspondence_cost to a uint16 in the [1, 32767] range.
// 将浮点数correspondence_cost转成[1, 32767]范围内的 uint16 整数
inline uint16 CorrespondenceCostToValue(const float correspondence_cost) {
  return BoundedFloatToValue(correspondence_cost, kMinCorrespondenceCost,
                             kMaxCorrespondenceCost);
}

// Converts a probability to a uint16 in the [1, 32767] range.
// 将浮点数probability转成[1, 32767]范围内的 uint16 整数
inline uint16 ProbabilityToValue(const float probability) {
  return BoundedFloatToValue(probability, kMinProbability, kMaxProbability);
}


// c++11: extern c风格
extern const std::vector<float>* const kValueToProbability;   //jc:extern说明这个值是从其他文件中获取来的，栅格占用概率
extern const std::vector<float>* const kValueToCorrespondenceCost; //jc:栅格是空闲的概率 在probability_value.cc80行进行初始化  ，这里.cc初始化，.h申明再调用

// Converts a uint16 (which may or may not have the update marker set) to a
// probability in the range [kMinProbability, kMaxProbability].
inline float ValueToProbability(const uint16 value) {
  return (*kValueToProbability)[value];
}

// Converts a uint16 (which may or may not have the update marker set) to a
// correspondence cost in the range [kMinCorrespondenceCost,
// kMaxCorrespondenceCost].
inline float ValueToCorrespondenceCost(const uint16 value) {
  return (*kValueToCorrespondenceCost)[value];
}

// 测试用的函数
inline uint16 ProbabilityValueToCorrespondenceCostValue(
    uint16 probability_value) {
  if (probability_value == kUnknownProbabilityValue) {
    return kUnknownCorrespondenceValue;
  }
  bool update_carry = false;
  if (probability_value > kUpdateMarker) {
    probability_value -= kUpdateMarker;
    update_carry = true;
  }
  uint16 result = CorrespondenceCostToValue(
      ProbabilityToCorrespondenceCost(ValueToProbability(probability_value)));
  if (update_carry) result += kUpdateMarker;
  return result;
}

// 测试用的函数
inline uint16 CorrespondenceCostValueToProbabilityValue(
    uint16 correspondence_cost_value) {
  if (correspondence_cost_value == kUnknownCorrespondenceValue)
    return kUnknownProbabilityValue;
  bool update_carry = false;
  if (correspondence_cost_value > kUpdateMarker) {
    correspondence_cost_value -= kUpdateMarker;
    update_carry = true;
  }
  uint16 result = ProbabilityToValue(CorrespondenceCostToProbability(
      ValueToCorrespondenceCost(correspondence_cost_value)));
  if (update_carry) result += kUpdateMarker;
  return result;
}

std::vector<uint16> ComputeLookupTableToApplyOdds(float odds);
std::vector<uint16> ComputeLookupTableToApplyCorrespondenceCostOdds(float odds);

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_PROBABILITY_VALUES_H_
