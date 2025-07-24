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
/// \file    DEQualityTask.cxx
/// \author  Andrea Ferrero andrea.ferrero@cern.ch
/// \brief   Post-processing of the MCH pre-clusters
/// \since   21/06/2022
///

#include "MCH/DEQualityTask.h"
#include "MCH/Helpers.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/ActivityHelpers.h"
#include "QualityControl/ReferenceUtils.h"
#include "Common/Utils.h"
#include "CCDB/CcdbObjectInfo.h"
#include "CommonUtils/MemFileHelper.h"
#include "TObjString.h"

#include <boost/property_tree/ptree.hpp>
#include <gsl/span>

namespace o2::aod::evsel
{
// Detector quality flags
enum DetectorQualityFlags {
  kCPVBad = 0,
  kEMCBad,
  kFDDBad,
  kFT0Bad,
  kFV0Bad,
  kHMPBad,
  kITSBad,
  kMCHBad,
  kMFTBad,
  kMIDBad,
  kPHSBad,
  kTOFBad,
  kTPCBadTracking = 0,
  kTPCBadPID,
  kTRDBad,
  kZDCBad,
  kNDetectorQualityFlags
};

class BunchCrossingTable
{
 public:
  bool detectorQuality_bit(int bit) const { return true; }
};

class CollisionTable
{
 public:
  bool detectorQuality_bit(int bit) const { return true; }
};

template <typename T>
concept HasDetectorQuality = requires(T a, int bit) {
  //{a.detectorQuality_bit(bit) } -> std::same_as<bool>;
  { a.detectorQuality_bit(bit) } -> std::convertible_to<bool>;
};

class QualityChecker
{
 public:
  QualityChecker() = default;
  QualityChecker(std::initializer_list<DetectorQualityFlags> bitsTocheck) : mBitsToCheck(bitsTocheck.size())
  {
    std::copy(bitsTocheck.begin(), bitsTocheck.end(), mBitsToCheck.begin());
  }

  QualityChecker(const std::vector<DetectorQualityFlags> bitsTocheck) : mBitsToCheck(bitsTocheck) {}

  // template <class TABLE>
  // bool checkTable(const TABLE& table)
  bool checkTable(const HasDetectorQuality auto& table)
  {
    if (mBitsToCheck.empty()) {
      throw std::out_of_range("QualityChecker with empty DetectorQualityFlags bits vector");
    }

    for (auto bit : mBitsToCheck) {
      if (table.detectorQuality_bit(bit)) {
        return false;
      }
    }
    return true;
  }

  // template <class TABLE>
  // bool operator ()(const TABLE& table)
  bool operator()(const HasDetectorQuality auto& table)
  {
    return checkTable(table);
  }

 private:
  std::vector<DetectorQualityFlags> mBitsToCheck;
};

class QualitySelectionFactory
{
 public:
  static QualityChecker create(const std::string& label)
  {
    if (label == "CBT") {
      return { kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID };
    } else if (label == "CBT_hadronPID") {
      return { kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kTOFBad };
    } else if (label == "CBT_electronPID") {
      return { kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kTRDBad };
    } else if (label == "CBT_calo") {
      return { kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kEMCBad };
    } else if (label == "CBT_muon") {
      return { kFT0Bad, kITSBad, kTPCBadTracking, kMCHBad, kMIDBad };
    } else if (label == "CBT_muon_glo") {
      return { kFT0Bad, kITSBad, kTPCBadTracking, kMCHBad, kMFTBad, kMIDBad };
    }

    return {};
  }
};

} // namespace o2::aod::evsel

using namespace o2::quality_control::core;
using namespace o2::quality_control::repository;
using namespace o2::quality_control_modules::common;
using namespace o2::quality_control_modules::muonchambers;

template <typename T>
o2::ccdb::CcdbObjectInfo createCcdbInfo(const T& object, uint64_t timeStamp, std::string_view reason)
{
  auto clName = o2::utils::MemFileHelper::getClassName(object);
  auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
  std::map<std::string, std::string> md;
  md["upload-reason"] = reason;
  constexpr auto fiveDays = 5 * o2::ccdb::CcdbObjectInfo::DAY;
  return o2::ccdb::CcdbObjectInfo("MCH/Calib/BadDE", clName, flName, md, timeStamp, timeStamp + fiveDays);
}

/*
std::string toCSV(const std::set<int>& deIds)
{
  std::stringstream csv;
  csv << fmt::format("deid\n");

  for (auto deId : deIds) {
    csv << fmt::format("{}\n", deId);
  }

  return csv.str();
}
*/
//_________________________________________________________________________________________
// Helper function for retrieving a MonitorObject from the QCDB, in the form of a std::pair<std::shared_ptr<MonitorObject>, bool>
// A non-null MO is returned in the first element of the pair if the MO is found in the QCDB
// The second element of the pair is set to true if the MO has a time stamp more recent than a user-supplied threshold

static std::pair<std::shared_ptr<MonitorObject>, bool> getMO(DatabaseInterface& qcdb, const std::string& fullPath, Trigger trigger, long notOlderThan)
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
  // retrieve QO from CCDB
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

void DEQualityTask::configure(const boost::property_tree::ptree& config)
{
  // input plots
  if (const auto& inputs = config.get_child_optional("qc.postprocessing." + getID() + ".inputs"); inputs.has_value()) {
    for (const auto& input : inputs.value()) {
      mPlotPaths.push_back(input.second.get_value<std::string>());
    }
  }
}

//_________________________________________________________________________________________

void DEQualityTask::initialize(Trigger t, framework::ServiceRegistryRef services)
{
  mCCDBpath = getFromExtendedConfig<std::string>(t.activity, mCustomParameters, "CCDBpath", mCCDBpath);
  mAPI.init(mCCDBpath);

  mObjectPath = getFromExtendedConfig<std::string>(t.activity, mCustomParameters, "objectPath", mObjectPath);

  auto detectorQuality = o2::aod::evsel::QualitySelectionFactory::create("CBT_muon");
  o2::aod::evsel::CollisionTable table;
  bool good = detectorQuality(table);
}

//_________________________________________________________________________________________

void DEQualityTask::update(Trigger trigger, framework::ServiceRegistryRef services)
{
  auto& qcdb = services.get<repository::DatabaseInterface>();

  std::set<int> badDEs;

  // fill vector of bad DEs
  for (auto plotPath : mPlotPaths) {
    auto [mo, success] = getMO(qcdb, plotPath, trigger, 600000);
    if (!success || !mo) {
      continue;
    }

    TH2F* hist = dynamic_cast<TH2F*>(mo->getObject());
    if (!hist) {
      ILOG(Warning, Devel) << "Could not cast the object '" << plotPath << "' to TH2F" << ENDM;
      continue;
    }

    if (hist->GetXaxis()->GetNbins() != getNumDE()) {
      ILOG(Warning, Devel) << "Wrong number of bins for object '" << plotPath << "': "
                           << hist->GetXaxis()->GetNbins() << " while " << getNumDE() << " were expected" << ENDM;
      continue;
    }

    for (int index = 0; index < hist->GetXaxis()->GetNbins(); index++) {
      bool isGood = (hist->GetBinContent(index + 1, 1) != 0);
      if (isGood) {
        continue;
      }

      int deId = getDEFromIndex(index);
      badDEs.insert(deId);

      ILOG(Info, Devel) << "Bad detection element DE" << deId << " found in \'" << plotPath << "\'" << ENDM;
    }
  }

  // check if the list of bad DEs has changed
  bool changed = false;
  if (mPreviousBadDEs.has_value()) {
    for (auto deId : badDEs) {
      if (!mPreviousBadDEs.value().count(deId) > 0) {
        changed = true;
        break;
      }
    }
  }

  if (changed) {
    std::string sBadDEs = "deid\n";
    for (auto deId : badDEs) {
      sBadDEs += fmt::format("{}\n", deId);
    }
    TObjString badDEsObject(sBadDEs.c_str());

    // time stamp
    auto uploadTS = o2::ccdb::getCurrentTimestamp();

    // CCDB object info
    auto clName = o2::utils::MemFileHelper::getClassName(badDEsObject);
    auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
    constexpr auto fiveDays = 5 * o2::ccdb::CcdbObjectInfo::DAY;
    auto info = o2::ccdb::CcdbObjectInfo(mObjectPath, clName, flName, std::map<std::string, std::string>(), uploadTS, uploadTS + fiveDays);

    // CCDB object image
    auto image = o2::ccdb::CcdbApi::createObjectImage(&badDEsObject, &info);

    ILOG(Info, Devel) << "Storing object " << info.getPath()
                      << " of type " << info.getObjectType()
                      << " / " << info.getFileName()
                      << " and size " << image->size()
                      << " bytes, valid for " << info.getStartValidityTimestamp()
                      << " : " << info.getEndValidityTimestamp() << ENDM;

    int res = mAPI.storeAsBinaryFile(image->data(), image->size(), info.getFileName(), info.getObjectType(), info.getPath(), info.getMetaData(), info.getStartValidityTimestamp(), info.getEndValidityTimestamp());
    if (res) {
      LOGP(error, "uploading to {} / {} failed for [{}:{}]", mAPI.getURL(), info.getPath(), info.getStartValidityTimestamp(), info.getEndValidityTimestamp());
    } else {
      mPreviousBadDEs = badDEs;
    }
  }
}

//_________________________________________________________________________________________

void DEQualityTask::finalize(Trigger t, framework::ServiceRegistryRef)
{
}
