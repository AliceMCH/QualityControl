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
/// \file    ReferenceComparatorTask.xx
/// \author  Andrea Ferrero
/// \brief   Post-processing task that compares a given set of plots with reference ones
///

#include "Common/ReferenceComparatorTask.h"
#include "Common/ReferenceComparatorPlot.h"
#include "Common/Utils.h"
#include "QualityControl/ReferenceUtils.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/ActivityHelpers.h"
#include "QualityControl/CcdbDatabase.h"
#include "DataFormatsCTP/CTPRateFetcher.h"
#include "DataFormatsParameters/GRPECSObject.h"
// ROOT
#include <TClass.h>
#include <TH1.h>

using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control::core;
using namespace o2::quality_control;
using namespace o2::quality_control::repository;

namespace o2::quality_control_modules::common
{

//_________________________________________________________________________________________
// Helper function for retrieving a MonitorObject from the QCDB, in the form of a std::pair<std::shared_ptr<MonitorObject>, bool>
// A non-null MO is returned in the first element of the pair if the MO is found in the QCDB
// The second element of the pair is set to true if the MO has a time stamp more recent than a user-supplied threshold

static std::pair<std::shared_ptr<MonitorObject>, bool> getMO(repository::DatabaseInterface& qcdb, const std::string& fullPath, Trigger trigger, long notOlderThan)
{
  // find the time-stamp of the most recent object matching the current activity
  // if ignoreActivity is true the activity matching criteria are not applied
  auto objectTimestamp = trigger.timestamp;
  const auto filterMetadata = activity_helpers::asDatabaseMetadata(trigger.activity, false);
  const auto objectValidity = qcdb.getLatestObjectValidity(trigger.activity.mProvenance + "/" + fullPath, filterMetadata);
  if (objectValidity.isValid()) {
    objectTimestamp = objectValidity.getMax() - 1;
  } else {
    ILOG(Warning, Devel) << "Could not find the object '" << fullPath << "' for activity " << trigger.activity << ENDM;
    return { nullptr, false };
  }

  auto [success, path, name] = o2::quality_control::core::RepoPathUtils::splitObjectPath(fullPath);
  if (!success) {
    return { nullptr, false };
  }
  // retrieve QO from CCDB - do not associate to trigger activity if ignoreActivity is true
  auto qo = qcdb.retrieveMO(path, name, objectTimestamp, trigger.activity);
  if (!qo) {
    return { nullptr, false };
  }

  long elapsed = static_cast<long>(trigger.timestamp) - objectTimestamp;
  // check if the object is not older than a given number of milliseconds
  if (elapsed > notOlderThan) {
    ILOG(Warning, Devel) << "Object '" << fullPath << "' for activity " << trigger.activity << " is too old: " << elapsed << " > " << notOlderThan << ENDM;
    return { qo, false };
  }

  return { qo, true };
}

//_________________________________________________________________________________________

void ReferenceComparatorTask::configure(const boost::property_tree::ptree& config)
{
  mConfig = ReferenceComparatorTaskConfig(getID(), config);
}

//_________________________________________________________________________________________

static std::string getCustomParameter(const o2::quality_control::core::CustomParameters& customParameters, const std::string& key, const Activity& activity, const std::string& defaultValue)
{
  std::string value;
  auto valueOptional = customParameters.atOptional(key, activity);
  if (valueOptional.has_value()) {
    value = valueOptional.value();
  } else {
    value = customParameters.atOptional(key).value_or(defaultValue);
  }

  return value;
}

//_________________________________________________________________________________________

void ReferenceComparatorTask::initialize(quality_control::postprocessing::Trigger trigger, framework::ServiceRegistryRef services)
{
  // reset all existing objects
  mPlotNames.clear();
  mReferencePlots.clear();
  mHistograms.clear();

  auto& qcdb = services.get<repository::DatabaseInterface>();
  mNotOlderThan = getFromExtendedConfig<int>(trigger.activity, mCustomParameters, "notOlderThan", 120);
  mReferenceRun = getFromExtendedConfig<int>(trigger.activity, mCustomParameters, "referenceRun", 0);
  mIgnorePeriodForReference = getFromExtendedConfig<bool>(trigger.activity, mCustomParameters, "ignorePeriodForReference", true);
  mIgnorePassForReference = getFromExtendedConfig<bool>(trigger.activity, mCustomParameters, "ignorePassForReference", true);

  ILOG(Info, Devel) << "Reference run set to '" << mReferenceRun << "' for activity " << trigger.activity << ENDM;

  auto referenceRunsList = o2::utils::Str::tokenize(getCustomParameter(mCustomParameters, "referenceRuns", trigger.activity, ""),
                                                    ';', false, true);

  for (const auto& referenceRunStr : referenceRunsList) {
    auto rateAndRun = o2::utils::Str::tokenize(referenceRunStr, ':', false, true);
    if (rateAndRun.size() == 2) {
      double rate = std::stod(rateAndRun[0]);
      size_t run = std::stoi(rateAndRun[1]);
      mReferenceRunForRate.insert(std::make_pair(rate, run));
    }
  }

  auto referenceActivity = trigger.activity;
  referenceActivity.mId = mReferenceRun;
  if (mIgnorePeriodForReference) {
    referenceActivity.mPeriodName = "";
  }
  if (mIgnorePassForReference) {
    referenceActivity.mPassName = "";
  }

  // initialise (new run)
  auto& ccdbManager = o2::ccdb::BasicCCDBManager::instance();
  ccdbManager.setURL("https://alice-ccdb.cern.ch");

  auto database = std::make_unique<CcdbDatabase>();
  database->connect("https://alice-ccdb.cern.ch", "", "", "");


  for (const auto& [rate, runNumber] : mReferenceRunForRate) {
    o2::ctp::CTPRateFetcher fetcher;
    fetcher.setupRun(runNumber, &ccdbManager, 1732165391957-1 /*repository::DatabaseInterface::Timestamp::Latest*/, false);

    // get scalers
    std::map<string, string> metadata;
    metadata["runNumber"] = std::to_string(runNumber);
    //auto grpECS = ccdbManager.getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", 1698896865746-1 /*repository::DatabaseInterface::Timestamp::Latest*/, metadata);
    o2::ctp::CTPRunScalers* ctpscalers = ccdbManager.getSpecific<ctp::CTPRunScalers>("CTP/Calib/Scalers", 1732179095000-1 /*repository::DatabaseInterface::Timestamp::Latest*/, metadata);
    //o2::ctp::CTPRunScalers* ctpscalers = ccdbManager.getSpecific<ctp::CTPRunScalers>("CTP/Calib/Scalers", repository::DatabaseInterface::Timestamp::Latest, metadata);
    ctpscalers->convertRawToO2();
    std::cout << "Before updateScalers()" << std::endl;
    fetcher.updateScalers(*ctpscalers);
    std::cout << "After updateScalers()" << std::endl;

    unsigned long timeStart;
    unsigned long timeEnd;
    std::tie(timeStart, timeEnd) = ctpscalers->getTimeLimit();
    std::cout << "After getTimeLimit()" << std::endl;
    unsigned long timeStep = 30000; // 30 seconds

    //string sourceName = "T0VTX"; // this will be provided by the user
    string sourceName = "ZNC-hadronic";
    for (unsigned long timestamp = timeStart; timestamp <= timeEnd; timestamp += timeStep) {
      auto rate = fetcher.fetchNoPuCorr(&ccdbManager, timestamp, runNumber, sourceName);
      std::cout << "Run " << runNumber << " timestamp " << timestamp << " source " << sourceName << " rate " << rate << std::endl;
    }
  }

  // load and initialize the input groups
  for (auto group : mConfig.dataGroups) {
    auto groupName = group.name;
    auto& plotVec = mPlotNames[groupName];
    for (auto path : group.objects) {
      auto fullPath = group.inputPath + "/" + path;
      auto fullRefPath = group.referencePath + "/" + path;
      auto fullOutPath = group.outputPath + "/" + path;

      // retrieve the reference MO
      auto referencePlot = o2::quality_control::checker::getReferencePlot(&qcdb, fullRefPath, referenceActivity);
      if (!referencePlot) {
        ILOG(Warning, Support) << "Could not load reference plot for object \"" << fullRefPath << "\" and activity " << referenceActivity << ENDM;
        continue;
      }

      // extract the reference histogram
      TH1* referenceHistogram = dynamic_cast<TH1*>(referencePlot->getObject());
      if (!referenceHistogram) {
        continue;
      }
      std::cout << "Loaded reference plot for object \"" << fullRefPath << "\" and activity " << referenceActivity
          << " - integral = " << referenceHistogram->Integral() << std::endl;

      // store the reference MO
      mReferencePlots[fullPath] = referencePlot;

      // fill an array with the full paths of the plots associated to this group
      plotVec.push_back(fullPath);

      // create and store the plotter object
      mHistograms[fullPath] = std::make_shared<ReferenceComparatorPlot>(referenceHistogram, mReferenceRun, fullOutPath,
                                                                        group.normalizeReference,
                                                                        group.drawRatioOnly,
                                                                        group.legendHeight,
                                                                        group.drawOption1D,
                                                                        group.drawOption2D);
      auto* outObject = mHistograms[fullPath]->getMainCanvas();
      // publish the object created by the plotter
      if (outObject) {
        getObjectsManager()->startPublishing(outObject);
      }
    }
  }
}

//_________________________________________________________________________________________

std::shared_ptr<ReferenceComparatorPlot> ReferenceComparatorTask::getComparatorPlot(std::string plotName)
{
  // check if a corresponding output plot was initialized
  auto iter = mHistograms.find(plotName);
  if (iter == mHistograms.end()) {
    return {};
  }
  return iter->second;
}

//_________________________________________________________________________________________

void ReferenceComparatorTask::update(quality_control::postprocessing::Trigger trigger, framework::ServiceRegistryRef services)
{
  auto& qcdb = services.get<repository::DatabaseInterface>();

  for (auto& [groupName, plotVec] : mPlotNames) {
    for (auto& plotName : plotVec) {
      // get object for current timestamp - age limit is converted to milliseconds
      auto object = getMO(qcdb, plotName, trigger, mNotOlderThan * 1000);

      // skip objects that are not found or too old
      if (!object.first || !object.second) {
        continue;
      }

      // only process objects inheriting from TH1
      auto* histogram = dynamic_cast<TH1*>(object.first->getObject());
      if (!histogram) {
        continue;
      }

      // check if a corresponding output plot was initialized
      auto iter = mHistograms.find(plotName);
      if (iter == mHistograms.end()) {
        continue;
      }

      // update the plot ratios and the histograms with superimposed reference
      iter->second->update(histogram);
    }
  }
}

//_________________________________________________________________________________________

void ReferenceComparatorTask::finalize(quality_control::postprocessing::Trigger t, framework::ServiceRegistryRef)
{
}

} // namespace o2::quality_control_modules::common
