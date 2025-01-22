// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   TrendCheck.h
/// \author Andrea Ferrero
/// \brief  Generic checker for trending graphs
///

#ifndef QC_MODULE_COMMON_TRENDCHECK_H
#define QC_MODULE_COMMON_TRENDCHECK_H

#include "QualityControl/CheckInterface.h"

#include <optional>
#include <array>
#include <unordered_map>

class TGraph;
class TObject;

namespace o2::ctp
{
  class CTPRateFetcher;
  class CTPRunScalers;
}

using namespace o2::quality_control::core;

namespace o2::quality_control_modules::common
{

class CheckerThresholdsConfig;

/// \brief  Generic checker for trending graphs
///
/// \author Andrea Ferrero
class TrendCheck : public o2::quality_control::checker::CheckInterface
{
 public:
  /// Default constructor
  TrendCheck() = default;
  /// Destructor
  ~TrendCheck() override = default;

  void configure() override;
  Quality check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap) override;
  void beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult = Quality::Null) override;
  std::string getAcceptedType() override;

  void startOfActivity(const Activity& activity) override;
  void endOfActivity(const Activity& activity) override;

  ClassDefOverride(TrendCheck, 1);

 private:
  enum ThresholdsMode {
    ExpectedRange,
    DeviationFromMean,
    StdDeviation
  };

  std::array<std::optional<std::pair<double, double>>, 2> getThresholds(std::string key, TGraph* graph, double rate);
  void getGraphsFromObject(TObject* object, std::vector<TGraph*>& graphs);

  void updateScalers(uint64_t timestamp);
  double getRateForTimestamp(uint64_t timestamp);
  double getRateForMO(MonitorObject* mo);

  Activity mActivity;
  ThresholdsMode mTrendCheckMode{ ExpectedRange };
  int mNPointsForAverage{ 0 };
  std::pair<float, float> mQualityLabelPosition{ 0.12, 0.8 };
  std::pair<float, float> mQualityLabelSize{ 0.5, 0.07 };
  std::shared_ptr<CheckerThresholdsConfig> mThresholds;
  std::unordered_map<std::string, std::vector<std::pair<double, std::pair<double, double>>>> mAverageTrend;
  std::unordered_map<std::string, std::vector<std::pair<double, std::pair<double, double>>>> mThresholdsTrendBad;
  std::unordered_map<std::string, std::vector<std::pair<double, std::pair<double, double>>>> mThresholdsTrendMedium;
  std::unordered_map<std::string, Quality> mQualities;

  /// \brief CCDB URL used to initialize the rate fetcher
  /// FIXME: this should be taken from the global QC configuration key "qc.config.conditionDB.url"?
  std::string mCTPConfigURL{ "https://alice-ccdb.cern.ch" };

  /// \brief URL and path for accessing the CTP scalers object
  std::string mCTPScalerURL{ "https://alice-ccdb.cern.ch" };
  std::string mCTPScalerPath{ "CTP/Calib/Scalers" };
  // in online we must access the temporary scaler objects in the QCDB
  //std::string mCTPScalerURL{ "http://ali-qcdb.cern.ch:8083" };
  //std::string mCTPScalerPath{ "qc/CTP/Scalers" };

  /// \brief trigger source for the collision rate computation
  /// To be set differently for pp and Pb-Pb
  std::string mCTPScalerSourceName{ "T0VTX" };
  //std::string mCTPScalerSourceName{ "ZNC-hadronic" };

  /// \brief minimum interval between updates of the CTP scalers object
  /// Set to zero to fetch the object only once
  int mScalerUpdateIntervalMs{ 0 };
  /// \brief time stamp of the most recent update of the scalers object
  uint64_t mLastScalerUpdateTimestamp{ 0 };

  /// helper object to fetch the collision rate for a given run number and time stamp
  std::shared_ptr<o2::ctp::CTPRateFetcher> mRateFetcher;
  o2::ctp::CTPRunScalers* mScalers{ nullptr };
};

} // namespace o2::quality_control_modules::common

#endif // QC_MODULE_COMMON_TRENDCHECK_H
