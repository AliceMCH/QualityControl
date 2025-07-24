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

#include "MUONCommon/TracksTask.h"

#include "MUONCommon/Helpers.h"
#include "QualityControl/ObjectsManager.h"
#include <DataFormatsMCH/Cluster.h>
#include <DataFormatsMCH/Digit.h>
#include <DataFormatsMCH/ROFRecord.h>
#include <DataFormatsMCH/TrackMCH.h>
#include <DetectorsBase/GeometryManager.h>
#include "QualityControl/QcInfoLogger.h"
#include <Framework/DataRefUtils.h>
#include <Framework/InputRecord.h>
#include <Framework/TimingInfo.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <ReconstructionDataFormats/TrackMCHMID.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <DataFormatsZDC/RecEventFlat.h>
#include <CommonConstants/LHCConstants.h>
#include <gsl/span>

#include <fstream>

namespace o2::quality_control_modules::muon
{

TracksTask::TracksTask()
{
}

TracksTask::~TracksTask() = default;

GID::mask_t adaptSource(GID::mask_t src)
{
  if (src[GID::Source::MFTMCHMID] == 1) {
    src.reset(GID::Source::MFTMCHMID); // does not exist
    src.set(GID::Source::MFTMCH);
    // ensure we request the individual tracks as we use their information in the plotter
    src.set(GID::Source::MFT);
    src.set(GID::Source::MCH);
    src.set(GID::Source::MID);
  }
  if (src[GID::Source::MFTMCH] == 1) {
    // ensure we request the individual tracks as we use their information in the plotter
    src.set(GID::Source::MFT);
    src.set(GID::Source::MCH);
  }
  if (src[GID::Source::MCHMID] == 1) {
    // ensure we request the individual tracks as we use their information in the plotter
    src.set(GID::Source::MCH);
    src.set(GID::Source::MID);
  }
  return src;
}

void TracksTask::initialize(o2::framework::InitContext& /*ic*/)
{
  ILOG(Debug, Devel) << "initialize TracksTask" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  mBackgroundITS.clear();
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108626425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108626425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108657381));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108657381));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108695964));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 108700781));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108722471));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108722471));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108734899));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108734899));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108756602));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108765870));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108765885));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 108778286));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108799981));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 108821668));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 108826126));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108826139));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 108852641));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 108852658));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 108852658));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108886778));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108886778));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108896055));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108896055));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 108930147));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108956294));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108977996));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 108995244));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109016929));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109021403));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109047915));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109047915));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109199778));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109233912));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109233912));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109255595));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109260035));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109286586));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109308283));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109320699));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109329978));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109329978));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109385784));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109472552));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109472564));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109472564));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109477005));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109477005));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109477020));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109498699));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109498699));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109525240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109525240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109542085));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109542092));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109542092));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109559349));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109559349));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109672277));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109715658));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109759065));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109824138));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109824151));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109824151));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109841379));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109850656));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109884786));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109910913));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109915773));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109915773));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109932620));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109932620));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109932621));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109932621));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109937443));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109937443));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109980859));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110014974));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110019393));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110210215));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110241190));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110241190));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110275303));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110275303));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110327993));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110349678));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110349678));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110393070));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110427188));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110427188));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110436464));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110453322));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110453322));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110501565));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110535679));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110561792));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110561792));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110566655));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110588332));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110610043));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110631714));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110631714));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110631722));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110687524));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110687538));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110691975));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110740212));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110740212));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110796014));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110796021));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110827003));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110843858));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110843858));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110913781));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110913781));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110935486));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110947883));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110952333));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110952333));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110952333));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111000556));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111000556));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111022252));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111034684));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111034684));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111060800));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111060828));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111143156));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111147609));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111164860));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111164860));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111217525));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111229946));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111260904));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111277770));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111277770));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111295036));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111295036));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111304299));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111304312));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111304312));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111360098));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111360098));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111364549));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111364549));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111364561));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111364561));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111364570));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111369391));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111381800));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111386253));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111386253));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111391072));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111391072));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111391097));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111391097));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111446893));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111468582));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111468582));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111468584));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111473053));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111473053));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111538129));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111538129));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111608036));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111651439));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111733377));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111733384));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111750648));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111750648));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111755085));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111755085));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111759924));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111815725));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111820186));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111820186));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111863557));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111885266));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111885266));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111933493));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111955188));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111976889));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111976889));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112010985));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112010985));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112015425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112015425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112037123));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112037123));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112037126));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112037126));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112063669));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112063669));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112119464));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112128744));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112128744));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112162851));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112210709));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112232397));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112232411));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112232411));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112237216));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112237216));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112297494));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112297494));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112389094));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112410803));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112423219));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112423219));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112466607));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112471053));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112488298));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112497582));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112497582));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112536158));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112536158));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112596771));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112596771));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112601230));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112601230));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112622914));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112622914));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112671158));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112683552));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112683552));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112705249));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112736240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112736240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112801318));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112801318));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112844704));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112844729));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112900518));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112900518));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 112900519));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112922212));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112922212));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112965625));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112965625));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112970074));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113035146));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113035146));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113039986));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113039986));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113039989));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113039989));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113052393));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113074088));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113083381));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113083381));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 113117496));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113126758));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113126758));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113139193));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113139193));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113160878));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113182562));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113182562));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 113208713));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113235240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113235240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113269348));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113382300));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113382300));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 113403982));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113408803));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113408803));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113452198));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113452198));

  ILOG(Info, Support) << "loading sources" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  auto srcFixed = mSrc;
  // For track type selection
  if (auto param = mCustomParameters.find("GID"); param != mCustomParameters.end()) {
    ILOG(Info, Devel) << "Custom parameter - GID (= sources by user): " << param->second << ENDM;
    ILOG(Info, Devel) << "Allowed sources           = " << mAllowedSources << " " << GID::getSourcesNames(mAllowedSources) << ENDM;
    auto requested = GID::getSourcesMask(param->second);
    ILOG(Info, Devel) << "Requested Sources         = " << requested << " " << GID::getSourcesNames(requested) << ENDM;
    mSrc = mAllowedSources & requested;
    srcFixed = adaptSource(mSrc);
    ILOG(Info, Devel) << "Allowed requested sources = " << mSrc << " " << GID::getSourcesNames(mSrc) << ENDM;
    ILOG(Info, Devel) << "Sources for data request  = " << srcFixed << " " << GID::getSourcesNames(srcFixed) << ENDM;
  }

  ILOG(Info, Support) << "Will do DataRequest for " << GID::getSourcesNames(srcFixed) << ENDM;
  if (srcFixed[GID::Source::MFTMCHMID] == 1) {
    srcFixed.reset(GID::Source::MFTMCHMID);
    srcFixed.set(GID::Source::MFTMCH);
  }
  mDataRequest = std::make_shared<o2::globaltracking::DataRequest>();
  mDataRequest->requestTracks(srcFixed, false);
}

void TracksTask::createTrackHistos(const Activity& activity)
{
  bool fullHistos = getConfigurationParameter<bool>(mCustomParameters, "fullHistos", false, activity);

  double maxTracksPerTF = getConfigurationParameter<double>(mCustomParameters, "maxTracksPerTF", 400, activity);
  double cutRAbsMin = getConfigurationParameter<double>(mCustomParameters, "cutRAbsMin", 17.6, activity);
  double cutRAbsMax = getConfigurationParameter<double>(mCustomParameters, "cutRAbsMax", 89.5, activity);
  double cutEtaMin = getConfigurationParameter<double>(mCustomParameters, "cutEtaMin", -4.0, activity);
  double cutEtaMax = getConfigurationParameter<double>(mCustomParameters, "cutEtaMax", -2.5, activity);
  double cutPtMin = getConfigurationParameter<double>(mCustomParameters, "cutPtMin", 0.5, activity);
  double cutPMin = getConfigurationParameter<double>(mCustomParameters, "cutPMin", 0, activity);
  double cutChi2Min = getConfigurationParameter<double>(mCustomParameters, "cutChi2Min", 0, activity);
  double cutChi2Max = getConfigurationParameter<double>(mCustomParameters, "cutChi2Max", 1000, activity);
  double nSigmaPDCA = getConfigurationParameter<double>(mCustomParameters, "nSigmaPDCA", 6, activity);
  double matchScoreMaxMFT = getConfigurationParameter<double>(mCustomParameters, "matchScoreMaxMFT", 1000, activity);
  double matchChi2MaxMFT = getConfigurationParameter<double>(mCustomParameters, "matchChi2MaxMFT", 1000, activity);
  double diMuonTimeCut = getConfigurationParameter<double>(mCustomParameters, "diMuonTimeCut", 100, activity) / 1000;

  int etaBins = getConfigurationParameter<int>(mCustomParameters, "etaBins", 200, activity);
  int phiBins = getConfigurationParameter<int>(mCustomParameters, "phiBins", 180, activity);
  int ptBins = getConfigurationParameter<int>(mCustomParameters, "ptBins", 300, activity);

  //======================================
  // Track plotters without cuts

  auto createPlotter = [&](GID::Source source, std::string path) {
    if (mSrc[source] == 1) {
      ILOG(Info, Devel) << "Creating plotter for path " << path << ENDM;
      mTrackPlotters[source] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, source, path, fullHistos);
      mTrackPlotters[source]->publish(getObjectsManager());
      mTrackPlottersAll.push_back(mTrackPlotters[source].get());
    }
  };

  createPlotter(GID::Source::MCH, "");
  createPlotter(GID::Source::MCHMID, "MCH-MID/");
  createPlotter(GID::Source::MFTMCH, "MFT-MCH/");
  createPlotter(GID::Source::MFTMCHMID, "MFT-MCH-MID/");

  //======================================
  // Track plotters with cuts

  std::vector<MuonCutFunc> muonCuts{
    // Rabs cut
    [cutRAbsMin, cutRAbsMax](const MuonTrack& t) { return ((t.getRAbs() >= cutRAbsMin) && (t.getRAbs() <= cutRAbsMax)); },
    // Eta cut
    [cutEtaMin, cutEtaMax](const MuonTrack& t) { return ((t.getMuonMomentumAtVertexMCH().eta() >= cutEtaMin) && (t.getMuonMomentumAtVertexMCH().eta() <= cutEtaMax)); },
    // Pt cut
    [cutPtMin](const MuonTrack& t) { return ((t.getMuonMomentumAtVertexMCH().Pt() >= cutPtMin)); },
    // P cut
    [cutPMin](const MuonTrack& t) { return ((std::fabs(t.getMuonMomentumMCH().P()) >= cutPMin)); },
    // pDCA cut
    [nSigmaPDCA](const MuonTrack& t) {
      static const double sigmaPDCA23 = 80.;
      static const double sigmaPDCA310 = 54.;
      static const double relPRes = 0.0004;
      static const double slopeRes = 0.0005;

      double thetaAbs = TMath::ATan(t.getRAbs() / 505.) * TMath::RadToDeg();

      double pUncorr = t.getTrackParamMCH().p();
      double p = t.getMuonMomentumAtVertexMCH().P();

      double pDCA = pUncorr * t.getDCAMCH();
      double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
      double nrp = nSigmaPDCA * relPRes * p;
      double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
      double slopeResEffect = 535. * slopeRes * p;
      double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
      if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
        return false;
      }

      return true;
    },
    // MFT-MCH match score
    [matchScoreMaxMFT](const MuonTrack& t) {
      if (t.hasMFT() && t.hasMCH() && t.getMatchInfoFwd().getMFTMCHMatchingScore() > matchScoreMaxMFT)
        return false;
      return true;
    },
    // MFT-MCH chi2 score
    [matchChi2MaxMFT](const MuonTrack& t) {
      if (t.hasMFT() && t.hasMCH() && t.getMatchInfoFwd().getMFTMCHMatchingChi2() > matchChi2MaxMFT)
        return false;
      return true;
    },
    // MCH chi2 cut
    [cutChi2Min, cutChi2Max](const MuonTrack& t) { return ((t.getChi2OverNDFMCH() >= cutChi2Min) && (t.getChi2OverNDFMCH() <= cutChi2Max)); }
  };

  // ZDC background selection/rejection
  std::vector<MuonCutFunc> zdcBdgSelection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& zdcIR : mBackgroundZDC) {
        auto diffIR = muonIR.toLong() - zdcIR.toLong() - 31;
        if (std::abs(diffIR) <= 5) {
          std::cout << "[PIPPO] ZDC background found at " << zdcIR << " (" << muonIR << ")" << std::endl;
          return true;
        }
      }
      return false;
    }
  };

  std::vector<MuonCutFunc> zdcBdgRejection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& zdcIR : mBackgroundZDC) {
        auto diffIR = muonIR.toLong() - zdcIR.toLong() - 31;
        if (std::abs(diffIR) <= 5) {
          return false;
        }
      }
      return true;
    }
  };

  // ITS background selection/rejection
  std::vector<MuonCutFunc> itsBdgSelection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& itsIR : mBackgroundITS) {
        auto diffIR = muonIR.toLong() - itsIR.toLong();
        if (diffIR >= 0 && diffIR <= 594) {
          std::cout << "[PIPPO] ITS background found at " << itsIR << " (" << muonIR << ")" << std::endl;
          return true;
        }
      }
      return false;
    }
  };

  std::vector<MuonCutFunc> itsBdgRejection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& itsIR : mBackgroundITS) {
        auto diffIR = muonIR.toLong() - itsIR.toLong();
        if (diffIR >= 0 && diffIR <= 594) {
          std::cout << "[PIPPO] ITS background rejected " << std::endl;
          return false;
        }
      }
      return true;
    }
  };

  std::vector<MuonCutFunc> bdgSelection{ zdcBdgSelection[0], itsBdgSelection[0] };
  std::vector<MuonCutFunc> bdgRejection{ zdcBdgRejection[0], itsBdgRejection[0] };

  std::vector<MuonCutFunc> zdcOnlyBdgSelection{ zdcBdgSelection[0], itsBdgRejection[0] };
  std::vector<MuonCutFunc> itsOnlyBdgSelection{ zdcBdgRejection[0], itsBdgSelection[0] };

  std::vector<DiMuonCutFunc> diMuonCuts{
    // cut on time difference between the two muon tracks
    [diMuonTimeCut](const MuonTrack& t1, const MuonTrack& t2) { return (std::abs(t1.getTime().getTimeStamp() - t2.getTime().getTimeStamp()) < diMuonTimeCut); }
  };

  auto createPlotterWithCuts = [&](std::map<GID::Source, std::unique_ptr<TrackPlotter>>& trackPlotters,
                                   GID::Source source,
                                   std::string path,
                                   std::vector<MuonCutFunc>& cuts) {
    if (mSrc[source] == 1) {
      ILOG(Info, Devel) << "Creating plotter with cuts for path " << path << ENDM;
      trackPlotters[source] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, source, path, fullHistos);
      trackPlotters[source]->setMuonCuts(cuts);
      trackPlotters[source]->setDiMuonCuts(diMuonCuts);
      trackPlotters[source]->publish(getObjectsManager());
      mTrackPlottersAll.push_back(trackPlotters[source].get());
      ILOG(Info, Devel) << "Plotter with cuts for path " << path << " created" << ENDM;
    }
  };

  createPlotterWithCuts(mTrackPlottersWithCuts, GID::Source::MCH, "WithCuts/", muonCuts);
  createPlotterWithCuts(mTrackPlottersWithCuts, GID::Source::MCHMID, "MCH-MID/WithCuts/", muonCuts);
  createPlotterWithCuts(mTrackPlottersWithCuts, GID::Source::MFTMCH, "MFT-MCH/WithCuts/", muonCuts);
  createPlotterWithCuts(mTrackPlottersWithCuts, GID::Source::MFTMCHMID, "MFT-MCH-MID/WithCuts/", muonCuts);

  std::vector<MuonCutFunc> muonFT0OrACuts {
    [&](const MuonTrack& track) {
      if (track.getIR().bc >= 250 && track.getIR().bc <= 355) return true;
      return false;
    }
  };
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrA, GID::Source::MCH, "WithBgdFT0OrA/", muonFT0OrACuts);
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrA, GID::Source::MCHMID, "MCH-MID/WithBgdFT0OrA/", muonFT0OrACuts);
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrA, GID::Source::MFTMCH, "MFT-MCH/WithBgdFT0OrA/", muonFT0OrACuts);
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrA, GID::Source::MFTMCHMID, "MFT-MCH-MID/WithBgdFT0OrA/", muonFT0OrACuts);

  muonFT0OrACuts.insert(muonFT0OrACuts.end(), muonCuts.begin(), muonCuts.end());
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrA, GID::Source::MCH, "WithCutsAndBgdFT0OrA/", muonFT0OrACuts);
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrA, GID::Source::MCHMID, "MCH-MID/WithCutsAndBgdFT0OrA/", muonFT0OrACuts);
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrA, GID::Source::MFTMCH, "MFT-MCH/WithCutsAndBgdFT0OrA/", muonFT0OrACuts);
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrA, GID::Source::MFTMCHMID, "MFT-MCH-MID/WithCutsAndBgdFT0OrA/", muonFT0OrACuts);

  std::vector<MuonCutFunc> muonFT0OrCCuts {
    [&](const MuonTrack& track) {
      if (track.getIR().bc >= 2120 && track.getIR().bc <= 2155) return true;
      if (track.getIR().bc >= 2888 && track.getIR().bc < 3006) return true;
      return false;
    }
  };
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrC, GID::Source::MCH, "WithBgdFT0OrC/", muonFT0OrCCuts);
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrC, GID::Source::MCHMID, "MCH-MID/WithBgdFT0OrC/", muonFT0OrCCuts);
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrC, GID::Source::MFTMCH, "MFT-MCH/WithBgdFT0OrC/", muonFT0OrCCuts);
  createPlotterWithCuts(mTrackPlottersWithBgdFT0OrC, GID::Source::MFTMCHMID, "MFT-MCH-MID/WithBgdFT0OrC/", muonFT0OrCCuts);

  muonFT0OrCCuts.insert(muonFT0OrCCuts.end(), muonCuts.begin(), muonCuts.end());
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrC, GID::Source::MCH, "WithCutsAndBgdFT0OrC/", muonFT0OrCCuts);
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrC, GID::Source::MCHMID, "MCH-MID/WithCutsAndBgdFT0OrC/", muonFT0OrCCuts);
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrC, GID::Source::MFTMCH, "MFT-MCH/WithCutsAndBgdFT0OrC/", muonFT0OrCCuts);
  createPlotterWithCuts(mTrackPlottersWithCutsAndBgdFT0OrC, GID::Source::MFTMCHMID, "MFT-MCH-MID/WithCutsAndBgdFT0OrC/", muonFT0OrCCuts);

  std::vector<MuonCutFunc> muonNoBgdCuts {
    [&](const MuonTrack& track) {
      if (track.getIR().bc > 405 && track.getIR().bc < 2000) return true;
      return false;
    }
  };
  createPlotterWithCuts(mTrackPlottersNoBgd, GID::Source::MCH, "NoBgd/", muonNoBgdCuts);
  createPlotterWithCuts(mTrackPlottersNoBgd, GID::Source::MCHMID, "MCH-MID/NoBgd/", muonNoBgdCuts);
  createPlotterWithCuts(mTrackPlottersNoBgd, GID::Source::MFTMCH, "MFT-MCH/NoBgd/", muonNoBgdCuts);
  createPlotterWithCuts(mTrackPlottersNoBgd, GID::Source::MFTMCHMID, "MFT-MCH-MID/NoBgd/", muonNoBgdCuts);

  muonNoBgdCuts.insert(muonNoBgdCuts.end(), muonCuts.begin(), muonCuts.end());
  createPlotterWithCuts(mTrackPlottersWithCutsNoBgd, GID::Source::MCH, "WithCutsNoBgd/", muonNoBgdCuts);
  createPlotterWithCuts(mTrackPlottersWithCutsNoBgd, GID::Source::MCHMID, "MCH-MID/WithCutsNoBgd/", muonNoBgdCuts);
  createPlotterWithCuts(mTrackPlottersWithCutsNoBgd, GID::Source::MFTMCH, "MFT-MCH/WithCutsNoBgd/", muonNoBgdCuts);
  createPlotterWithCuts(mTrackPlottersWithCutsNoBgd, GID::Source::MFTMCHMID, "MFT-MCH-MID/WithCutsNoBgd/", muonNoBgdCuts);

  mTrackPlottersBgdZDC[0] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdZDC/", fullHistos);
  mTrackPlottersBgdZDC[0]->setMuonCuts(zdcBdgSelection);
  mTrackPlottersBgdZDC[0]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[0]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[1] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "NoBgdZDC/", fullHistos);
  mTrackPlottersBgdZDC[1]->setMuonCuts(zdcBdgRejection);
  mTrackPlottersBgdZDC[1]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[1]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[2] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdITS/", fullHistos);
  mTrackPlottersBgdZDC[2]->setMuonCuts(itsBdgSelection);
  mTrackPlottersBgdZDC[2]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[2]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[3] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "NoBgdITS/", fullHistos);
  mTrackPlottersBgdZDC[3]->setMuonCuts(itsBdgRejection);
  mTrackPlottersBgdZDC[3]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[3]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[4] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdITSZDC/", fullHistos);
  mTrackPlottersBgdZDC[4]->setMuonCuts(bdgSelection);
  mTrackPlottersBgdZDC[4]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[4]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[5] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "NoBgdITSZDC/", fullHistos);
  mTrackPlottersBgdZDC[5]->setMuonCuts(bdgRejection);
  mTrackPlottersBgdZDC[5]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[5]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[6] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdZDC-NoBgdITS/", fullHistos);
  mTrackPlottersBgdZDC[6]->setMuonCuts(zdcOnlyBdgSelection);
  mTrackPlottersBgdZDC[6]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[6]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[7] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdITS-NoBgdZDC/", fullHistos);
  mTrackPlottersBgdZDC[7]->setMuonCuts(itsOnlyBdgSelection);
  mTrackPlottersBgdZDC[7]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[7]->publish(getObjectsManager());

  mBcZDC = std::make_unique<TH1F>("BcZDC", "BcZDC;bc", o2::constants::lhc::LHCMaxBunches, 0, o2::constants::lhc::LHCMaxBunches);
  getObjectsManager()->startPublishing(mBcZDC.get());

  mDCAvsBcZDC = std::make_unique<TH2F>("DCAvsBcZDC", "DCAvsBcZDC;bc;DCA (cm)", 200, -100, 100, 100, 0, 100);
  getObjectsManager()->startPublishing(mDCAvsBcZDC.get());

  mBgdZDCTrackMult = std::make_unique<TH1F>("BgdZDCTrackMult", "Tracks multiplicity - ZDC background;# of tracks;", 100, 0, 100);
  getObjectsManager()->startPublishing(mBgdZDCTrackMult.get());
}

void TracksTask::removeTrackHistos()
{
  ILOG(Debug, Devel) << "Un-publishing objects" << ENDM;
  for (auto p : mTrackPlottersAll) {
    p->unpublish(getObjectsManager());
  }
  /*
  for (auto& p : mTrackPlotters) {
    p.second->unpublish(getObjectsManager());
  }
  for (auto& p : mTrackPlottersWithCuts) {
    p.second->unpublish(getObjectsManager());
  }
  */
  ILOG(Debug, Devel) << "Destroying objects" << ENDM;
  mTrackPlotters.clear();
  mTrackPlottersWithCuts.clear();
  mTrackPlottersAll.clear();
}

void TracksTask::startOfActivity(const Activity& activity)
{
  ILOG(Debug, Devel) << "startOfActivity : " << activity << ENDM;
  createTrackHistos(activity);
}

void TracksTask::startOfCycle()
{
  ILOG(Debug, Devel) << "startOfCycle" << ENDM;
}

bool TracksTask::assertInputs(o2::framework::ProcessingContext& ctx)
{
  if (!ctx.inputs().isValid("trackMCH")) {
    ILOG(Info, Support) << "no mch tracks available on input" << ENDM;
    return false;
  }
  if (!ctx.inputs().isValid("trackMCHROF")) {
    ILOG(Info, Support) << "no mch track rofs available on input" << ENDM;
    return false;
  }
  if (!ctx.inputs().isValid("trackMCHTRACKCLUSTERS")) {
    ILOG(Info, Support) << "no mch track clusters available on input" << ENDM;
    return false;
  }
  if (mSrc[GID::Source::MCHMID] == 1) {
    if (!ctx.inputs().isValid("matchMCHMID")) {
      ILOG(Info, Support) << "no muon (mch+mid) track available on input" << ENDM;
      return false;
    }
    if (!ctx.inputs().isValid("trackMID")) {
      ILOG(Info, Support) << "no mid track available on input" << ENDM;
      return false;
    }
  }
  if (mSrc[GID::Source::MFTMCH] == 1) {
    if (!ctx.inputs().isValid("fwdtracks")) {
      ILOG(Info, Support) << "no muon (mch+mft) track available on input" << ENDM;
      return false;
    }
  }
  if (mSrc[GID::Source::MFTMCHMID] == 1) {
    if (!ctx.inputs().isValid("fwdtracks")) {
      ILOG(Info, Support) << "no muon (mch+mft+mid) track available on input" << ENDM;
      return false;
    }
  }
  return true;
}

void TracksTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  ILOG(Debug, Devel) << "Debug: MonitorData" << ENDM;

  int firstTForbit = ctx.services().get<o2::framework::TimingInfo>().firstTForbit;
  ILOG(Debug, Devel) << "Debug: firstTForbit=" << firstTForbit << ENDM;

  if (!assertInputs(ctx)) {
    return;
  }

  ILOG(Debug, Devel) << "Debug: Asserted inputs" << ENDM;

  mRecoCont.collectData(ctx, *mDataRequest.get());
/*
  // ZDC data
  auto bcrecZDC = ctx.inputs().get<gsl::span<o2::zdc::BCRecData>>("zdc-bcrec");
  auto energyZDC = ctx.inputs().get<gsl::span<o2::zdc::ZDCEnergy>>("zdc-energyrec");
  auto tdcZDC = ctx.inputs().get<gsl::span<o2::zdc::ZDCTDCData>>("zdc-tdcrec");
  auto infoZDC = ctx.inputs().get<gsl::span<uint16_t>>("zdc-inforec");

  o2::zdc::RecEventFlat ev;
  ev.init(bcrecZDC, energyZDC, tdcZDC, infoZDC);

  mBackgroundZDC.clear();
  while (ev.next()) {
    int32_t itdc = o2::zdc::TDCZNAC;
    int nhit = ev.NtdcV(itdc);

    for (int32_t ipos = 0; ipos < nhit; ipos++) {
      double mytdc = o2::zdc::FTDCVal * ev.TDCVal[itdc][ipos];

      if (mytdc > 5.7 && mytdc < 8.7) {
        // Backgroud event found here!
        std::cout << "[TOTO] ZDC background event found: " << ev.ir << std::endl;
        mBackgroundZDC.emplace_back(ev.ir);
        mBcZDC->Fill(ev.ir.bc);
        for (const auto& itsIR : mBackgroundITS) {
          if (itsIR.orbit == ev.ir.orbit) {
            std::cout << "ITS+ZDC background event in same orbit: ZDC [" << ev.ir << "]  ITS [ " << itsIR << "]  BC difference: " << (ev.ir.toLong() - itsIR.toLong()) << std::endl;
          }
        }
      }
    }
  }
*/
  ILOG(Debug, Devel) << "Debug: Collected data" << ENDM;

  for (auto p : mTrackPlottersAll) {
    if (p) {
      p->setFirstTForbit(firstTForbit);
    }
  }
  /*
  for (auto& p : mTrackPlotters) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersWithCuts) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersWithBgd) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersNoBgd) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersWithCutsAndBgd) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersWithCutsNoBgd) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  */

  if (mSrc[GID::MCH] == 1) {
    ILOG(Debug, Devel) << "Debug: MCH requested" << ENDM;
    if (mRecoCont.isTrackSourceLoaded(GID::MCH)) {
      ILOG(Debug, Devel) << "Debug: MCH source loaded" << ENDM;
      mTrackPlotters[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrA[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrC[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersNoBgd[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrA[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrC[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsNoBgd[GID::MCH]->fillHistograms(mRecoCont);
    }
  }
  if (mSrc[GID::MCHMID] == 1) {
    ILOG(Debug, Devel) << "Debug: MCHMID requested" << ENDM;
    if (mRecoCont.isMatchSourceLoaded(GID::MCHMID)) {
      ILOG(Debug, Devel) << "Debug: MCHMID source loaded" << ENDM;
      mTrackPlotters[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrA[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrC[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersNoBgd[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrA[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrC[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsNoBgd[GID::MCHMID]->fillHistograms(mRecoCont);
    }
  }
  if (mSrc[GID::MFTMCH] == 1) {
    ILOG(Debug, Devel) << "Debug: MFTMCH requested" << ENDM;
    if (mRecoCont.isTrackSourceLoaded(GID::MFTMCH)) {
      mTrackPlotters[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrA[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrC[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersNoBgd[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrA[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrC[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsNoBgd[GID::MFTMCH]->fillHistograms(mRecoCont);
    }
  }
  if (mSrc[GID::MFTMCHMID] == 1) {
    ILOG(Debug, Devel) << "Debug: MFTMCHMID requested" << ENDM;
    if (mRecoCont.isTrackSourceLoaded(GID::MFTMCH)) {
      mTrackPlotters[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrA[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithBgdFT0OrC[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersNoBgd[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrA[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsAndBgdFT0OrC[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCutsNoBgd[GID::MFTMCHMID]->fillHistograms(mRecoCont);
    }
  }
}

void TracksTask::endOfCycle()
{
  ILOG(Debug, Devel) << "endOfCycle" << ENDM;
  for (auto p : mTrackPlottersAll) {
    p->endOfCycle();
  }
  /*
  for (auto& p : mTrackPlotters) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersWithCuts) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersWithBgd) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersNoBgd) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersWithCutsAndBgd) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersWithCutsNoBgd) {
    p.second->endOfCycle();
  }*/
}

void TracksTask::endOfActivity(const Activity& /*activity*/)
{
  ILOG(Debug, Devel) << "endOfActivity" << ENDM;
  removeTrackHistos();
}

void TracksTask::reset()
{
  ILOG(Debug, Devel) << "reset" << ENDM;
  for (auto p : mTrackPlottersAll) {
    p->reset();
  }
  /*
  for (auto& p : mTrackPlotters) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersWithCuts) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersWithBgd) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersNoBgd) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersWithCutsAndBgd) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersWithCutsNoBgd) {
    p.second->reset();
  }*/
}

} // namespace o2::quality_control_modules::muon
