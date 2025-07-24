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

#include "MUONCommon/TrackPlotter.h"
#include "MUONCommon/Helpers.h"

#include "DetectorsBase/GRPGeomHelper.h"
#include <DataFormatsITSMFT/ROFRecord.h>
#include <DataFormatsMFT/TrackMFT.h>
#include <DataFormatsMCH/TrackMCH.h>
#include <ReconstructionDataFormats/TrackMCHMID.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <MCHTracking/TrackExtrap.h>
#include <DetectorsBase/GeometryManager.h>
#include <Framework/DataRefUtils.h>
#include <Framework/InputRecord.h>
#include <CommonConstants/LHCConstants.h>
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include <TH2F.h>
#include <gsl/span>
#include <set>

namespace
{

static void setXAxisLabels(TProfile* h)
{
  TAxis* axis = h->GetXaxis();
  for (int i = 1; i <= 10; i++) {
    auto label = fmt::format("CH{}", i);
    axis->SetBinLabel(i, label.c_str());
  }
}

template <class HIST>
static void Fill(std::unique_ptr<HIST>& hist, double x)
{
  if (!hist) {
    return;
  }
  hist->Fill(x);
}

template <class HIST>
static void Fill(std::unique_ptr<HIST>& hist, double x, double y)
{
  if (!hist) {
    return;
  }
  hist->Fill(x, y);
}

template <>
void Fill<TH1DRatio>(std::unique_ptr<TH1DRatio>& hist, double x)
{
  if (!hist) {
    return;
  }
  hist->getNum()->Fill(x);
}

template <>
void Fill<TH2DRatio>(std::unique_ptr<TH2DRatio>& hist, double x, double y)
{
  if (!hist) {
    return;
  }
  hist->getNum()->Fill(x, y);
}

} // namespace

using namespace o2::dataformats;

namespace o2::quality_control_modules::muon
{

TrackPlotter::TrackPlotter(int maxTracksPerTF,
                           int etaBins, int phiBins, int ptBins,
                           GID::Source source,
                           std::string path,
                           bool fullHistos)
  : mSrc(source),
    mPath(path),
    mFullHistos(fullHistos)
{
  createTrackHistos(maxTracksPerTF, etaBins, phiBins, ptBins);
  createTrackPairHistos();
}

template <>
std::unique_ptr<TH1DRatio> TrackPlotter::createHisto(const char* name, const char* title,
                                                     int nbins, double xmin, double xmax,
                                                     bool optional,
                                                     bool statBox,
                                                     const char* drawOptions,
                                                     const char* displayHints)
{
  if (optional && !mFullHistos) {
    return nullptr;
  }
  std::string fullTitle = GID::getSourceName(mSrc) + " " + title;
  auto h = std::make_unique<TH1DRatio>(name, fullTitle.c_str(), nbins, xmin, xmax, true);
  if (!statBox) {
    h->SetStats(0);
  }
  histograms().emplace_back(HistInfo{ h.get(), drawOptions, displayHints });
  return h;
}

template <>
std::unique_ptr<TH1DRatio> TrackPlotter::createHisto(const char* name, const char* title,
                                                     int nbins, double* xbins,
                                                     bool optional,
                                                     bool statBox,
                                                     const char* drawOptions,
                                                     const char* displayHints)
{
  if (optional && !mFullHistos) {
    return nullptr;
  }
  std::string fullTitle = GID::getSourceName(mSrc) + " " + title;
  auto h = std::make_unique<TH1DRatio>(name, fullTitle.c_str(), nbins, xbins, true);
  if (!statBox) {
    h->SetStats(0);
  }
  histograms().emplace_back(HistInfo{ h.get(), drawOptions, displayHints });
  return h;
}

template <>
std::unique_ptr<TH2DRatio> TrackPlotter::createHisto(const char* name, const char* title,
                                                     int nbins, double xmin, double xmax,
                                                     int nbinsy, double ymin, double ymax,
                                                     bool optional,
                                                     bool statBox,
                                                     const char* drawOptions,
                                                     const char* displayHints)
{
  if (optional && !mFullHistos) {
    return nullptr;
  }
  std::string fullTitle = GID::getSourceName(mSrc) + " " + title;
  auto h = std::make_unique<TH2DRatio>(name, fullTitle.c_str(), nbins, xmin, xmax, nbinsy, ymin, ymax, true);
  if (!statBox) {
    h->SetStats(0);
  }
  histograms().emplace_back(HistInfo{ h.get(), drawOptions, displayHints });
  return h;
}

void TrackPlotter::createTrackHistos(int maxTracksPerTF, int etaBins, int phiBins, int ptBins)
{
  int nbins = 100;
  auto xLogBins = makeLogBinning(1, maxTracksPerTF, nbins);
  mNofTracksPerTF[0] = createHisto<TH1D>(TString::Format("%sPositive/TracksPerTF", mPath.c_str()), "Number of tracks per TimeFrame (+);Number of tracks per TF", nbins, xLogBins.data(), true, false, "hist logx");
  mNofTracksPerTF[1] = createHisto<TH1D>(TString::Format("%sNegative/TracksPerTF", mPath.c_str()), "Number of tracks per TimeFrame (-);Number of tracks per TF", nbins, xLogBins.data(), true, false, "hist logx");
  mNofTracksPerTF[2] = createHisto<TH1D>(TString::Format("%sTracksPerTF", mPath.c_str()), "Number of tracks per TimeFrame;Number of tracks per TF", nbins, xLogBins.data(), false, false, "hist logx");

  mTrackChi2OverNDF[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackMCHChi2OverNDF", mPath.c_str()), "Track #chi^{2}/ndf (MCH +);#chi^{2}/ndf;entries/s", 500, 0, 50, true, false, "hist");
  mTrackChi2OverNDF[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackMCHChi2OverNDF", mPath.c_str()), "Track #chi^{2}/ndf (MCH -);#chi^{2}/ndf;entries/s", 500, 0, 50, true, false, "hist");
  mTrackChi2OverNDF[2] = createHisto<TH1DRatio>(TString::Format("%sTrackMCHChi2OverNDF", mPath.c_str()), "Track #chi^{2}/ndf (MCH);#chi^{2}/ndf;entries/s", 500, 0, 50, false, false, "hist");

  mTrackDCA[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackDCA", mPath.c_str()), "Track DCA (+);DCA (cm);entries/s", 500, 0, 500, true, false, "hist");
  mTrackDCA[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackDCA", mPath.c_str()), "Track DCA (-);DCA (cm);entries/s", 500, 0, 500, true, false, "hist");
  mTrackDCA[2] = createHisto<TH1DRatio>(TString::Format("%sTrackDCA", mPath.c_str()), "Track DCA;DCA (cm);entries/s", 500, 0, 500, false, false, "hist");

  mTrackPDCA[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackPDCA", mPath.c_str()), "Track p#timesDCA (+);p#timesDCA (GeVcm/c);entries/s", 5000, 0, 5000, true, false, "hist");
  mTrackPDCA[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackPDCA", mPath.c_str()), "Track p#timesDCA (-);p#timesDCA (GeVcm/c);entries/s", 5000, 0, 5000, true, false, "hist");
  mTrackPDCA[2] = createHisto<TH1DRatio>(TString::Format("%sTrackPDCA", mPath.c_str()), "Track p#timesDCA;p#timesDCA (GeVcm/c);entries/s", 5000, 0, 5000, false, false, "hist");

  mTrackRAbs[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackRAbs", mPath.c_str()), "Track R_{abs} (+);R_{abs} (cm);entries/s", 1000, 0, 100, true, false, "hist");
  mTrackRAbs[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackRAbs", mPath.c_str()), "Track R_{abs} (-);R_{abs} (cm);entries/s", 1000, 0, 100, true, false, "hist");
  mTrackRAbs[2] = createHisto<TH1DRatio>(TString::Format("%sTrackRAbs", mPath.c_str()), "Track R_{abs};R_{abs} (cm);entries/s", 1000, 0, 100, false, false, "hist");

  mTrackEta[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackEta", mPath.c_str()), "Track #eta (+);#eta;entries/s", etaBins, -4.5, -2, true, false, "hist");
  mTrackEta[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackEta", mPath.c_str()), "Track #eta (-);#eta;entries/s", etaBins, -4.5, -2, true, false, "hist");
  mTrackEta[2] = createHisto<TH1DRatio>(TString::Format("%sTrackEta", mPath.c_str()), "Track #eta;#eta;entries/s", etaBins, -4.5, -2, false, false, "hist");

  mTrackPhi[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackPhi", mPath.c_str()), "Track #phi (+);#phi (deg);entries/s", phiBins, -180, 180, true, false, "hist");
  mTrackPhi[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackPhi", mPath.c_str()), "Track #phi (-);#phi (deg);entries/s", phiBins, -180, 180, true, false, "hist");
  mTrackPhi[2] = createHisto<TH1DRatio>(TString::Format("%sTrackPhi", mPath.c_str()), "Track #phi;#phi (deg);entries/s", phiBins, -180, 180, false, false, "hist");

  mTrackEtaPhi[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackEtaPhi", mPath.c_str()), "Track #phi vs #eta (+);#eta;#phi", etaBins / 5, -4.5, -2, phiBins / 5, -180, 180, true, false, "colz");
  mTrackEtaPhi[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackEtaPhi", mPath.c_str()), "Track #phi vs #eta (-);#eta;#phi", etaBins / 5, -4.5, -2, phiBins / 5, -180, 180, true, false, "colz");
  mTrackEtaPhi[2] = createHisto<TH2DRatio>(TString::Format("%sTrackEtaPhi", mPath.c_str()), "Track #phi vs #eta;#eta;#phi", etaBins / 5, -4.5, -2, phiBins / 5, -180, 180, false, false, "colz");

  mTrackPt[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackPt", mPath.c_str()), "Track p_{T} (+);p_{T} (GeV/c);entries/s", ptBins, 0, 30, true, false, "hist logy");
  mTrackPt[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackPt", mPath.c_str()), "Track p_{T} (-);p_{T} (GeV/c);entries/s", ptBins, 0, 30, true, false, "hist logy");
  mTrackPt[2] = createHisto<TH1DRatio>(TString::Format("%sTrackPt", mPath.c_str()), "Track p_{T};p_{T} (GeV/c);entries/s", ptBins, 0, 30, false, false, "hist logy");

  mTrackP[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackP", mPath.c_str()), "Track momentum (+);p (GeV/c);entries/s", 2000, 0, 1000, true, false, "hist logy");
  mTrackP[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackP", mPath.c_str()), "Track momentum (-);p (GeV/c);entries/s", 2000, 0, 1000, true, false, "hist logy");
  mTrackP[2] = createHisto<TH1DRatio>(TString::Format("%sTrackP", mPath.c_str()), "Track p_{T};p (GeV/c);entries/s", 2000, 0, 1000, false, false, "hist logy");

  mTrackQOverPt = createHisto<TH1DRatio>(TString::Format("%sTrackQOverPt", mPath.c_str()), "Track q/p_{T};q/p_{T} (GeV/c)^{-1};entries/s", 200, -10, 10, false, false, "hist logy");

  mTrackEtaPt[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackEtaPt", mPath.c_str()), "Track p_{T} vs #eta (+);#eta;p_{T} (GeV/c)", etaBins / 5, -4.5, -2, ptBins / 5, 0, 30, true, false, "colz");
  mTrackEtaPt[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackEtaPt", mPath.c_str()), "Track p_{T} vs #eta (-);#eta;p_{T} (GeV/c)", etaBins / 5, -4.5, -2, ptBins / 5, 0, 30, true, false, "colz");
  mTrackEtaPt[2] = createHisto<TH2DRatio>(TString::Format("%sTrackEtaPt", mPath.c_str()), "Track p_{T} vs #eta;#eta;p_{T} (GeV/c)", etaBins / 5, -4.5, -2, ptBins / 5, 0, 30, true, false, "colz");

  mTrackPhiPt[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackPhiPt", mPath.c_str()), "Track p_{T} vs #phi (+);#phi;p_{T} (GeV/c)", phiBins / 5, -180, 180, ptBins / 5, 0, 30, true, false, "colz");
  mTrackPhiPt[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackPhiPt", mPath.c_str()), "Track p_{T} vs #phi (-);#phi;p_{T} (GeV/c)", phiBins / 5, -180, 180, ptBins / 5, 0, 30, true, false, "colz");
  mTrackPhiPt[2] = createHisto<TH2DRatio>(TString::Format("%sTrackPhiPt", mPath.c_str()), "Track p_{T} vs #phi;#phi;p_{T} (GeV/c)", phiBins / 5, -180, 180, ptBins / 5, 0, 30, true, false, "colz");

  mTrackBC = createHisto<TH1DRatio>(TString::Format("%sTrackBC", mPath.c_str()), "Track BC;BC;entries/s", o2::constants::lhc::LHCMaxBunches, 0, o2::constants::lhc::LHCMaxBunches, false, false, "hist");

  if (mSrc == GID::MFTMCH || mSrc == GID::MFTMCHMID) {
    mMatchScoreMFTMCH = createHisto<TH1D>(TString::Format("%sMatchScoreMFTMCH", mPath.c_str()), "Match Score MFT-MCH;score", 1000, 0, 100, false, false, "hist");
    mMatchChi2MFTMCH = createHisto<TH1D>(TString::Format("%sMatchChi2MFTMCH", mPath.c_str()), "Match #chi^{2} MFT-MCH;#chi^{2}", 1000, 0, 100, false, false, "hist");

    mMatchNMFTCandidates = createHisto<TH1D>(TString::Format("%sMatchNMFTCandidates", mPath.c_str()), "MFT Candidates;candidates", 1000, 0, 1000, false, false, "hist");
    mMatchNMFTCandidatesVsROFSize = createHisto<TH2F>(TString::Format("%sMatchNMFTCandidatesVsROFSize", mPath.c_str()), "MFT Candidates vs. ROF size;candidates;ROF size", 1000, 0, 1000, 1000, 0, 1000, false, false, "colz");

    mMatchDRMFTMCH = createHisto<TH1D>(TString::Format("%sMatchDRMFTMCH", mPath.c_str()), "MFT-MCH dR;dR (cm)", 1000, -50, 50, false, false, "hist");
    mMatchDPhiMFTMCH = createHisto<TH1D>(TString::Format("%sMatchDPhiMFTMCH", mPath.c_str()), "MFT-MCH #Delta#phi;#Delta#phi (rad)", 1000, -0.5, 0.5, false, false, "hist");
    mMatchDXMFTMCH = createHisto<TH1D>(TString::Format("%sMatchDXMFTMCH", mPath.c_str()), "MFT-MCH dX;dX (cm)", 1000, -50, 50, false, false, "hist");
    mMatchDYMFTMCH = createHisto<TH1D>(TString::Format("%sMatchDYMFTMCH", mPath.c_str()), "MFT-MCH dY;dY (cm)", 1000, -50, 50, false, false, "hist");

    mMatchDXMFTMCHVsX = createHisto<TH2F>(TString::Format("%sMatchDXMFTMCHVsX", mPath.c_str()), "MFT-MCH dX vs. X;X (cm);dX (cm)", 30, -15, 15, 200, -10, 10, false, false, "colz");
    mMatchDXMFTMCHVsY = createHisto<TH2F>(TString::Format("%sMatchDXMFTMCHVsY", mPath.c_str()), "MFT-MCH dX vs. Y;Y (cm);dX (cm)", 30, -15, 15, 200, -10, 10, false, false, "colz");
    mMatchDYMFTMCHVsX = createHisto<TH2F>(TString::Format("%sMatchDYMFTMCHVsX", mPath.c_str()), "MFT-MCH dY vs. X;X (cm);dY (cm)", 30, -15, 15, 200, -10, 10, false, false, "colz");
    mMatchDYMFTMCHVsY = createHisto<TH2F>(TString::Format("%sMatchDYMFTMCHVsY", mPath.c_str()), "MFT-MCH dY vs. Y;Y (cm);dY (cm)", 30, -15, 15, 200, -10, 10, false, false, "colz");

    mAbsDXMFTMCHVsX = createHisto<TH2F>(TString::Format("%sAbsDXMFTMCHVsX", mPath.c_str()), "MFT-MCH dX vs. X at absorber;X (cm);dX (cm)", 50, -100, 100, 200, -10, 10, false, false, "colz");
    mAbsDXMFTMCHVsY = createHisto<TH2F>(TString::Format("%sAbsDXMFTMCHVsY", mPath.c_str()), "MFT-MCH dX vs. Y at absorber;Y (cm);dX (cm)", 50, -100, 100, 200, -10, 10, false, false, "colz");
    mAbsDYMFTMCHVsX = createHisto<TH2F>(TString::Format("%sAbsDYMFTMCHVsX", mPath.c_str()), "MFT-MCH dY vs. X at absorber;X (cm);dY (cm)", 50, -100, 100, 200, -10, 10, false, false, "colz");
    mAbsDYMFTMCHVsY = createHisto<TH2F>(TString::Format("%sAbsDYMFTMCHVsY", mPath.c_str()), "MFT-MCH dY vs. Y at absorber;Y (cm);dY (cm)", 50, -100, 100, 200, -10, 10, false, false, "colz");

    mVertexXVsMatchingXMFT = createHisto<TH2F>(TString::Format("%sVertexXVsMatchingXMFT", mPath.c_str()), "MFT X_{vertex} vs. X_{matching};X_{matching} (cm);X_{vertex} (cm)", 60, -15, 15, 400, -10, 10, false, false, "colz");
    mVertexYVsMatchingXMFT = createHisto<TH2F>(TString::Format("%sVertexYVsMatchingXMFT", mPath.c_str()), "MFT Y_{vertex} vs. X_{matching};X_{matching} (cm);Y_{vertex} (cm)", 60, -15, 15, 400, -10, 10, false, false, "colz");
    mVertexXVsMatchingYMFT = createHisto<TH2F>(TString::Format("%sVertexXVsMatchingYMFT", mPath.c_str()), "MFT X_{vertex} vs. Y_{matching};Y_{matching} (cm);X_{vertex} (cm)", 60, -15, 15, 400, -10, 10, false, false, "colz");
    mVertexYVsMatchingYMFT = createHisto<TH2F>(TString::Format("%sVertexYVsMatchingYMFT", mPath.c_str()), "MFT Y_{vertex} vs. Y_{matching};Y_{matching} (cm);Y_{vertex} (cm)", 60, -15, 15, 400, -10, 10, false, false, "colz");
    mVertexXVsMatchingXMCH = createHisto<TH2F>(TString::Format("%sVertexXVsXMCH", mPath.c_str()), "MCH X_{vertex} vs. X_;X (cm);X_{vertex} (cm)", 60, -120, 120, 400, -10, 10, false, false, "colz");
    mVertexYVsMatchingXMCH = createHisto<TH2F>(TString::Format("%sVertexYVsXMCH", mPath.c_str()), "MCH Y_{vertex} vs. X;X (cm);Y_{vertex} (cm)", 60, -120, 120, 400, -10, 10, false, false, "colz");
    mVertexXVsMatchingYMCH = createHisto<TH2F>(TString::Format("%sVertexXVsYMCH", mPath.c_str()), "MCH X_{vertex} vs. Y;Y (cm);X_{vertex} (cm)", 60, -120, 120, 400, -10, 10, false, false, "colz");
    mVertexYVsMatchingYMCH = createHisto<TH2F>(TString::Format("%sVertexYVsYMCH", mPath.c_str()), "MCH Y_{vertex} vs. Y;Y (cm);Y_{vertex} (cm)", 60, -120, 120, 400, -10, 10, false, false, "colz");

    mTrackEtaGlobal[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackEtaGlobal", mPath.c_str()), "Global track #eta (+);#eta;entries/s", etaBins, -4.5, -2, true, false, "hist");
    mTrackEtaGlobal[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackEtaGlobal", mPath.c_str()), "Global track #eta (-);#eta;entries/s", etaBins, -4.5, -2, true, false, "hist");
    mTrackEtaGlobal[2] = createHisto<TH1DRatio>(TString::Format("%sTrackEtaGlobal", mPath.c_str()), "Global track #eta;#eta;entries/s", etaBins, -4.5, -2, false, false, "hist");

    mTrackPhiGlobal[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackPhiGlobal", mPath.c_str()), "Global track #phi (+);#phi (deg);entries/s", phiBins, -180, 180, true, false, "hist");
    mTrackPhiGlobal[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackPhiGlobal", mPath.c_str()), "Global track #phi (-);#phi (deg);entries/s", phiBins, -180, 180, true, false, "hist");
    mTrackPhiGlobal[2] = createHisto<TH1DRatio>(TString::Format("%sTrackPhiGlobal", mPath.c_str()), "Global track #phi;#phi (deg);entries/s", phiBins, -180, 180, false, false, "hist");

    mTrackPtGlobal[0] = createHisto<TH1DRatio>(TString::Format("%sPositive/TrackPtGlobal", mPath.c_str()), "Global track p_{T} (+);p_{T} (GeV/c);entries/s", ptBins, 0, 30, true, false, "hist logy");
    mTrackPtGlobal[1] = createHisto<TH1DRatio>(TString::Format("%sNegative/TrackPtGlobal", mPath.c_str()), "Global track p_{T} (-);p_{T} (GeV/c);entries/s", ptBins, 0, 30, true, false, "hist logy");
    mTrackPtGlobal[2] = createHisto<TH1DRatio>(TString::Format("%sTrackPtGlobal", mPath.c_str()), "Global track p_{T};p_{T} (GeV/c);entries/s", ptBins, 0, 30, false, false, "hist logy");

    mTrackQOverPtGlobal = createHisto<TH1DRatio>(TString::Format("%sTrackQOverPtGlobal", mPath.c_str()), "Global track q/p_{T};q/p_{T} (GeV/c)^{-1};entries/s", 200, -10, 10, false, false, "hist logy");

    mTrackEtaPhiGlobal[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackEtaPhiGlobal", mPath.c_str()), "Global track #phi vs #eta (+);#eta;#phi", etaBins / 5, -4.5, -2, phiBins / 5, -180, 180, true, false, "colz");
    mTrackEtaPhiGlobal[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackEtaPhiGlobal", mPath.c_str()), "Global track #phi vs #eta (-);#eta;#phi", etaBins / 5, -4.5, -2, phiBins / 5, -180, 180, true, false, "colz");
    mTrackEtaPhiGlobal[2] = createHisto<TH2DRatio>(TString::Format("%sTrackEtaPhiGlobal", mPath.c_str()), "Global track #phi vs #eta;#eta;#phi", etaBins / 5, -4.5, -2, phiBins / 5, -180, 180, false, false, "colz");

    mTrackEtaPtGlobal[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackEtaPtGlobal", mPath.c_str()), "Global track p_{T} vs #eta (+);#eta;p_{T} (GeV/c)", etaBins / 5, -4.5, -2, ptBins / 5, 0, 30, true, false, "colz");
    mTrackEtaPtGlobal[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackEtaPtGlobal", mPath.c_str()), "Global track p_{T} vs #eta (-);#eta;p_{T} (GeV/c)", etaBins / 5, -4.5, -2, ptBins / 5, 0, 30, true, false, "colz");
    mTrackEtaPtGlobal[2] = createHisto<TH2DRatio>(TString::Format("%sTrackEtaPtGlobal", mPath.c_str()), "Global track p_{T} vs #eta;#eta;p_{T} (GeV/c)", etaBins / 5, -4.5, -2, ptBins / 5, 0, 30, true, false, "colz");

    mTrackPhiPtGlobal[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackPhiPtGlobal", mPath.c_str()), "Global track p_{T} vs #phi (+);#phi;p_{T} (GeV/c)", phiBins / 5, -180, 180, ptBins / 5, 0, 30, true, false, "colz");
    mTrackPhiPtGlobal[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackPhiPtGlobal", mPath.c_str()), "Global track p_{T} vs #phi (-);#phi;p_{T} (GeV/c)", phiBins / 5, -180, 180, ptBins / 5, 0, 30, true, false, "colz");
    mTrackPhiPtGlobal[2] = createHisto<TH2DRatio>(TString::Format("%sTrackPhiPtGlobal", mPath.c_str()), "Global track p_{T} vs #phi;#phi;p_{T} (GeV/c)", phiBins / 5, -180, 180, ptBins / 5, 0, 30, true, false, "colz");

    mTrackEtaCorr[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackEtaCorr", mPath.c_str()), "Track #eta - GLO vs MCH (+);#eta^{MCH};#eta^{GLO}", etaBins / 5, -4.5, -2, etaBins / 5, -4.5, -2, true, false, "colz");
    mTrackEtaCorr[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackEtaCorr", mPath.c_str()), "Track #eta - GLO vs MCH (-);#eta^{MCH};#eta^{GLO}", etaBins / 5, -4.5, -2, etaBins / 5, -4.5, -2, true, false, "colz");
    mTrackEtaCorr[2] = createHisto<TH2DRatio>(TString::Format("%sTrackEtaCorr", mPath.c_str()), "Track #eta - GLO vs MCH;#eta^{MCH};#eta^{GLO}", etaBins / 5, -4.5, -2, etaBins / 5, -4.5, -2, true, false, "colz");

    mTrackDEtaVsEta[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackDEtaVsEta", mPath.c_str()), "Track #eta^{GLO}-#eta^{MCH} vs #eta^{MCH} (+);#eta^{MCH};#eta^{GLO}-#eta^{MCH}", etaBins / 5, -4.5, -2, 200, -1, 1, true, false, "colz");
    mTrackDEtaVsEta[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackDEtaVsEta", mPath.c_str()), "Track #eta^{GLO}-#eta^{MCH} vs #eta^{MCH} (-);#eta^{MCH};#eta^{GLO}-#eta^{MCH}", etaBins / 5, -4.5, -2, 200, -1, 1, true, false, "colz");
    mTrackDEtaVsEta[2] = createHisto<TH2DRatio>(TString::Format("%sTrackDEtaVsEta", mPath.c_str()), "Track #eta^{GLO}-#eta^{MCH} vs #eta^{MCH};#eta^{MCH};#eta^{GLO}-#eta^{MCH}", etaBins / 5, -4.5, -2, 200, -1, 1, false, false, "colz");

    mTrackPhiCorr[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackPhiCorr", mPath.c_str()), "Track #phi - GLO vs MCH (+);#phi^{MCH};#phi^{GLO}", phiBins / 5, -180, 180, phiBins / 5, -180, 180, true, false, "colz");
    mTrackPhiCorr[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackPhiCorr", mPath.c_str()), "Track #phi - GLO vs MCH (-);#phi^{MCH};#phi^{GLO}", phiBins / 5, -180, 180, phiBins / 5, -180, 180, true, false, "colz");
    mTrackPhiCorr[2] = createHisto<TH2DRatio>(TString::Format("%sTrackPhiCorr", mPath.c_str()), "Track #phi - GLO vs MCH;#phi^{MCH};#phi^{GLO}", phiBins / 5, -180, 180, phiBins / 5, -180, 180, true, false, "colz");

    mTrackDPhiVsPhi[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackDPhiVsPhi", mPath.c_str()), "Track #phi^{GLO}-#phi^{MCH} vs #phi^{MCH} (+);#phi^{MCH};#phi^{GLO}-#phi^{MCH}", phiBins / 5, -180, 180, 200, -100, 100, true, false, "colz");
    mTrackDPhiVsPhi[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackDPhiVsPhi", mPath.c_str()), "Track #phi^{GLO}-#phi^{MCH} vs #phi^{MCH} (-);#phi^{MCH};#phi^{GLO}-#phi^{MCH}", phiBins / 5, -180, 180, 200, -100, 100, true, false, "colz");
    mTrackDPhiVsPhi[2] = createHisto<TH2DRatio>(TString::Format("%sTrackDPhiVsPhi", mPath.c_str()), "Track #phi^{GLO}-#phi^{MCH} vs #phi^{MCH};#phi^{MCH};#phi^{GLO}-#phi^{MCH}", phiBins / 5, -180, 180, 200, -100, 100, false, false, "colz");

    mTrackPtCorr[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackPtCorr", mPath.c_str()), "Track p_{T} - GLO vs MCH (+);p_{T}^{MCH};p_{T}^{GLO}", ptBins / 5, 0, 30, ptBins / 5, 0, 30, true, false, "colz");
    mTrackPtCorr[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackPtCorr", mPath.c_str()), "Track p_{T} - GLO vs MCH (-);p_{T}^{MCH};p_{T}^{GLO}", ptBins / 5, 0, 30, ptBins / 5, 0, 30, true, false, "colz");
    mTrackPtCorr[2] = createHisto<TH2DRatio>(TString::Format("%sTrackPtCorr", mPath.c_str()), "Track p_{T} - GLO vs MCH;p_{T}^{MCH};p_{T}^{GLO}", ptBins / 5, 0, 30, ptBins / 5, 0, 30, true, false, "colz");

    mTrackDPtVsPt[0] = createHisto<TH2DRatio>(TString::Format("%sPositive/TrackDPtVsPt", mPath.c_str()), "Track p_{T}^{GLO}-p_{T}^{MCH} vs p_{T}^{MCH} (+);p_{T}^{MCH};p_{T}^{GLO}-p_{T}^{MCH}", ptBins / 5, 0, 30, 200, -10, 10, true, false, "colz");
    mTrackDPtVsPt[1] = createHisto<TH2DRatio>(TString::Format("%sNegative/TrackDPtVsPt", mPath.c_str()), "Track p_{T}^{GLO}-p_{T}^{MCH} vs p_{T}^{MCH} (-);p_{T}^{MCH};p_{T}^{GLO}-p_{T}^{MCH}", ptBins / 5, 0, 30, 200, -10, 10, true, false, "colz");
    mTrackDPtVsPt[2] = createHisto<TH2DRatio>(TString::Format("%sTrackDPtVsPt", mPath.c_str()), "Track p_{T}^{GLO}-p_{T}^{MCH} vs p_{T}^{MCH};p_{T}^{MCH};p_{T}^{GLO}-p_{T}^{MCH}", ptBins / 5, 0, 30, 200, -10, 10, false, false, "colz");

    mTrackPosAtMatchingPlaneMCH = createHisto<TH2DRatio>(TString::Format("%sTrackPosAtMatchingPlaneMCH", mPath.c_str()), "MCH Track position at MFT-MCH matching plane;X (cm);Y (cm)", 100, -50, 50, 100, -50, 50, false, false, "colz");
    mTrackPosAtMatchingPlaneMFT = createHisto<TH2DRatio>(TString::Format("%sTrackPosAtMatchingPlaneMFT", mPath.c_str()), "MFT Track position at MFT-MCH matching plane;X (cm);Y (cm)", 100, -50, 50, 100, -50, 50, false, false, "colz");
  }

  if (mSrc == GID::MCHMID || mSrc == GID::MFTMCHMID) {
    mMatchChi2MCHMID = createHisto<TH1D>(TString::Format("%sMatchChi2MCHMID", mPath.c_str()), "Match #chi^{2} MCH-MID;#chi^{2}", 1000, 0, 100, false, true, "hist");
    mTrackDT = createHisto<TH1D>(TString::Format("%sTrackDT", mPath.c_str()), "MCH-MID time correlation;ns", 4000, -500, 500, false, true, "hist");
  }

  mTrackPosAtVertex = createHisto<TH2DRatio>(TString::Format("%sTrackPosAtVertex", mPath.c_str()), "MCH Track position at vertex;X (cm);Y (cm)", 200, -200, 200, 200, -200, 200, false, false, "colz");
  mTrackPosAtAbsorber = createHisto<TH2DRatio>(TString::Format("%sTrackPosAtAbsorber", mPath.c_str()), "MCH Track position at absorber exit;X (cm);Y (cm)", 200, -200, 200, 200, -200, 200, false, false, "colz");
  mTrackPosAtMID = createHisto<TH2DRatio>(TString::Format("%sTrackPosAtMID", mPath.c_str()), "MCH Track position at MID entrance;X (cm);Y (cm)", 200, -200, 200, 200, -200, 200, false, false, "colz");

  mTrackRPhiAtVertex = createHisto<TH2DRatio>(TString::Format("%sTrackRPhiAtVertex", mPath.c_str()), "MCH Track (R, #phi) at vertex;#phi (deg);R (cm)", 90, -180, 180, 200, 0, 200, false, false, "colz");
  mTrackRPhiAtAbsorber = createHisto<TH2DRatio>(TString::Format("%sTrackRPhiAtAbsorber", mPath.c_str()), "MCH Track (R, #phi) at absorber exit;#phi (deg);R (cm)", 90, -180, 180, 200, 0, 200, false, false, "colz");
  mTrackRPhiAtMID = createHisto<TH2DRatio>(TString::Format("%sTrackRPhiAtMID", mPath.c_str()), "MCH Track (R, #phi) at MID entrance;#phi (deg);R (cm)", 90, -180, 180, 200, 0, 200, false, false, "colz");

  mTrackRVsThetaAtVertex = createHisto<TH2DRatio>(TString::Format("%sTrackRVsThetaAtVertex", mPath.c_str()), "MCH Track (R, #theta) at vertex;#theta (deg);R (cm)", 100, 0, 10, 200, 0, 200, false, false, "colz");
  mTrackRVsThetaAtAbsorber = createHisto<TH2DRatio>(TString::Format("%sTrackRVsThetaAtAbsorber", mPath.c_str()), "MCH Track (R, #theta) at absorber exit;#theta (deg);R (cm)", 100, 0, 10, 200, 0, 200, false, false, "colz");
  mTrackRVsThetaAtMID = createHisto<TH2DRatio>(TString::Format("%sTrackRVsThetaAtMID", mPath.c_str()), "MCH Track (R, #theta) at MID entrance;#theta (deg);R (cm)", 100, 0, 10, 200, 0, 200, false, false, "colz");

  mTrackRVsMomAtVertex = createHisto<TH2DRatio>(TString::Format("%sTrackRVsMomAtVertex", mPath.c_str()), "MCH Track (R, mom) at vertex;p (GeV/c);R (cm)", 200, 0, 1000, 200, 0, 200, false, false, "colz");
  mTrackRVsMomAtAbsorber = createHisto<TH2DRatio>(TString::Format("%sTrackRVsMomAtAbsorber", mPath.c_str()), "MCH Track (R, mom) at absorber exit;p (GeV/c);R (cm)", 200, 0, 1000, 200, 0, 200, false, false, "colz");
  mTrackRVsMomAtMID = createHisto<TH2DRatio>(TString::Format("%sTrackRVsMomAtMID", mPath.c_str()), "MCH Track (R, mom) at MID entrance;p (GeV/c);R (cm)", 200, 0, 1000, 200, 0, 200, false, false, "colz");

  mTrackXSlopeVsPhiAtVertex = createHisto<TH2DRatio>(TString::Format("%sTrackXSlopeVsPhiAtVertex", mPath.c_str()), "MCH Track (xslope, #phi) at vertex;#phi (deg);x_slope", 90, -180, 180, 200, -0.1, 0.1, false, false, "colz");
  mTrackXSlopeVsPhiAtAbsorber = createHisto<TH2DRatio>(TString::Format("%sTrackXSlopeVsPhiAtAbsorber", mPath.c_str()), "MCH Track (xslope, #phi) at absorber exit;#phi (deg);x_slope", 90, -180, 180, 200, -0.1, 0.1, false, false, "colz");
  mTrackXSlopeVsPhiAtMID = createHisto<TH2DRatio>(TString::Format("%sTrackXSlopeVsPhiAtMID", mPath.c_str()), "MCH Track (xslope, #phi) at MID entrance;#phi (deg);x_slope", 90, -180, 180, 200, -0.1, 0.1, false, false, "colz");

  mTrackYSlopeVsPhiAtVertex = createHisto<TH2DRatio>(TString::Format("%sTrackYSlopeVsPhiAtVertex", mPath.c_str()), "MCH Track (yslope, #phi) at vertex;#phi (deg);y_slope", 90, -180, 180, 200, -0.1, 0.1, false, false, "colz");
  mTrackYSlopeVsPhiAtAbsorber = createHisto<TH2DRatio>(TString::Format("%sTrackYSlopeVsPhiAtAbsorber", mPath.c_str()), "MCH Track (yslope, #phi) at absorber exit;#phi (deg);y_slope", 90, -180, 180, 200, -0.1, 0.1, false, false, "colz");
  mTrackYSlopeVsPhiAtMID = createHisto<TH2DRatio>(TString::Format("%sTrackYSlopeVsPhiAtMID", mPath.c_str()), "MCH Track (yslope, #phi) at MID entrance;#phi (deg);y_slope", 90, -180, 180, 200, -0.1, 0.1, false, false, "colz");

  mSigmaXvsP = createHisto<TH2D>(TString::Format("%sSigmaXvsP", mPath.c_str()), "MCH Track X resolution vs. P at matching plane;momentum (GeV/c);#sigma_{X} (cm)", 200, 0, 1000, 100, 0, 10, false, false, "colz");
  mSigmaYvsP = createHisto<TH2D>(TString::Format("%sSigmaYvsP", mPath.c_str()), "MCH Track Y resolution vs. P at matching plane;momentum (GeV/c);#sigma_{Y} (cm)", 200, 0, 1000, 100, 0, 10, false, false, "colz");
}

void TrackPlotter::createTrackPairHistos()
{
  mMinvFull = createHisto<TH1DRatio>(TString::Format("%sMinvFull", mPath.c_str()), "#mu^{+}#mu^{-} invariant mass;M_{#mu^{+}#mu^{-}} (GeV/c^{2})", 5000, 0, 100, false, true, "hist");
  mMinv = createHisto<TH1DRatio>(TString::Format("%sMinv", mPath.c_str()), "#mu^{+}#mu^{-} invariant mass;M_{#mu^{+}#mu^{-}} (GeV/c^{2})", 200, 1, 5, false, true, "hist");
  mMinvBgd = createHisto<TH1DRatio>(TString::Format("%sMinvBgd", mPath.c_str()), "#mu^{+}#mu^{-} inv. mass background;M_{#mu^{+}#mu^{-}} (GeV/c^{2})", 200, 1, 5, true, true, "hist");
  mDimuonDT = createHisto<TH1DRatio>(TString::Format("%sDimuonTimeDiff", mPath.c_str()), "#mu^{+}#mu^{-} time difference;ns", 4000, -2000, 2000, false, true, "hist");
}

void TrackPlotter::normalizePlot(TH1* hist)
{
  static constexpr double sOrbitLengthInSeconds = o2::constants::lhc::LHCOrbitMUS / 1000000;

  TH1DRatio* h1 = dynamic_cast<TH1DRatio*>(hist);
  if (h1) {
    h1->getDen()->Fill((Double_t)0, sOrbitLengthInSeconds * mNOrbitsPerTF);
  } else {
    TH2DRatio* h2 = dynamic_cast<TH2DRatio*>(hist);
    if (h2) {
      h2->getDen()->Fill((Double_t)0, (Double_t)0, sOrbitLengthInSeconds * mNOrbitsPerTF);
    }
  }
}

void TrackPlotter::fillTrackPairHistos(gsl::span<const std::pair<MuonTrack, bool>> tracks)
{
  if (tracks.size() > 1) {
    for (auto i = 0; i < tracks.size(); i++) {
      if (!(tracks[i].second)) {
        continue;
      }
      const MuonTrack& ti = tracks[i].first;
      auto pi = ti.getMuonMomentumAtVertex();

      for (auto j = i + 1; j < tracks.size(); j++) {
        if (!(tracks[j].second)) {
          continue;
        }
        const MuonTrack& tj = tracks[j].first;

        if (ti.getSign() == tj.getSign()) {
          continue;
        }

        auto dt = (ti.getIR() - tj.getIR()).bc2ns();
        auto dtMUS = ti.getTime().getTimeStamp() - tj.getTime().getTimeStamp();
        Fill(mDimuonDT, dtMUS * 1000);

        bool diMuonOK = true;
        for (auto& cut : mDiMuonCuts) {
          if (!cut(ti, tj)) {
            diMuonOK = false;
            break;
          }
        }
        // tracks are considered to be correlated if they are closer than 1 us in time
        auto pj = tj.getMuonMomentumAtVertex();
        auto p = pi + pj;
        if (diMuonOK) {
          Fill(mMinv, p.M());
          Fill(mMinvFull, p.M());
        }
        // the shape of the combinatorial background is derived by combining tracks
        // than belong to different orbits (more than 90 us apart)
        if (std::abs(dtMUS) > 1) {
          Fill(mMinvBgd, p.M());
        }
      }
    }
  }
}

void TrackPlotter::fillTrackHistos(const MuonTrack& track)
{
  int q = (track.getSign() < 0) ? 1 : 0;

  Fill(mTrackBC, track.getIR().bc);

  auto dtMUS = track.getTimeMID().getTimeStamp() - track.getTimeMCH().getTimeStamp();
  Fill(mTrackDT, dtMUS * 1000);

  if (mSrc == GID::MCH) {
    auto mchTrack = track.getTrackParamMCH();
    //std::cout << std::format("MCH:  X = {:6.2f}  Y = {:6.2f}  Z = {:6.2f}  P = {:12.6f}\n",
    //    mchTrack.getNonBendingCoor(), mchTrack.getBendingCoor(), mchTrack.getZ(), mchTrack.p());
  }

  switch (mSrc) {
    case GID::MCHMID: {
      Fill(mMatchChi2MCHMID, track.getMatchInfoFwd().getMIDMatchingChi2());
      break;
    }
    case GID::MFTMCHMID: {
      Fill(mMatchChi2MCHMID, track.getMatchInfoFwd().getMIDMatchingChi2());
    }
    case GID::MFTMCH: {
      Fill(mMatchScoreMFTMCH, track.getMatchInfoFwd().getMFTMCHMatchingScore());
      Fill(mMatchChi2MFTMCH, track.getMatchInfoFwd().getMFTMCHMatchingChi2());
      Fill(mMatchNMFTCandidates, track.getMatchInfoFwd().getNMFTCandidates());
      Fill(mMatchNMFTCandidatesVsROFSize, track.getMatchInfoFwd().getNMFTCandidates(), track.getRofMFT().getNEntries());

      Fill(mTrackPosAtMatchingPlaneMCH, track.getXMatchMCH(), track.getYMatchMCH());
      Fill(mTrackPosAtMatchingPlaneMFT, track.getXMatchMFT(), track.getYMatchMFT());

      double RMCH = TMath::Sqrt(track.getXMatchMCH() * track.getXMatchMCH() + track.getYMatchMCH() * track.getYMatchMCH());
      double RMFT = TMath::Sqrt(track.getXMatchMFT() * track.getXMatchMFT() + track.getYMatchMFT() * track.getYMatchMFT());
      double dR = RMFT - RMCH;
      Fill(mMatchDRMFTMCH, dR);

      Fill(mMatchDXMFTMCH, track.getXMatchMFT() - track.getXMatchMCH());
      Fill(mMatchDYMFTMCH, track.getYMatchMFT() - track.getYMatchMCH());

      Fill(mMatchDXMFTMCHVsX, track.getXMatchMFT(), track.getXMatchMFT() - track.getXMatchMCH());
      Fill(mMatchDXMFTMCHVsY, track.getYMatchMFT(), track.getXMatchMFT() - track.getXMatchMCH());
      Fill(mMatchDYMFTMCHVsX, track.getXMatchMFT(), track.getYMatchMFT() - track.getYMatchMCH());
      Fill(mMatchDYMFTMCHVsY, track.getYMatchMFT(), track.getYMatchMFT() - track.getYMatchMCH());

      Fill(mSigmaXvsP, std::fabs(track.getMuonMomentumMCH().P()), track.getSigmaXMatchMCH());
      Fill(mSigmaYvsP, std::fabs(track.getMuonMomentumMCH().P()), track.getSigmaYMatchMCH());

      Fill(mAbsDXMFTMCHVsX, track.getXAbsMFT(), track.getXAbsMFT() - track.getXAbsMCH());
      Fill(mAbsDXMFTMCHVsY, track.getYAbsMFT(), track.getXAbsMFT() - track.getXAbsMCH());
      Fill(mAbsDYMFTMCHVsX, track.getXAbsMFT(), track.getYAbsMFT() - track.getYAbsMCH());
      Fill(mAbsDYMFTMCHVsY, track.getYAbsMFT(), track.getYAbsMFT() - track.getYAbsMCH());

      Fill(mVertexXVsMatchingXMFT, track.getXMatchMFT(), track.getXVertexMFT());
      Fill(mVertexYVsMatchingXMFT, track.getXMatchMFT(), track.getYVertexMFT());
      Fill(mVertexXVsMatchingYMFT, track.getYMatchMFT(), track.getXVertexMFT());
      Fill(mVertexYVsMatchingYMFT, track.getYMatchMFT(), track.getYVertexMFT());

      Fill(mVertexXVsMatchingXMCH, track.getTrackParamMCH().getNonBendingCoor(), track.getXVertexMCH());
      Fill(mVertexYVsMatchingXMCH, track.getTrackParamMCH().getNonBendingCoor(), track.getYVertexMCH());
      Fill(mVertexXVsMatchingYMCH, track.getTrackParamMCH().getBendingCoor(), track.getXVertexMCH());
      Fill(mVertexYVsMatchingYMCH, track.getTrackParamMCH().getBendingCoor(), track.getYVertexMCH());

      //std::cout << std::format("QC MFT extrapolation at matching: {:0.2},{:0.2},{:0.2}\n", track.getXMatchMFT(), track.getYMatchMFT(), track.getYMatchMFT());
      //std::cout << std::format("QC MFT extrapolation at vertex:   {:0.2},{:0.2},{:0.2}\n", track.getXVertexMFT(), track.getYVertexMFT(), track.getZVertexMFT());

      break;
    }
    default:
      break;
  }

  double chi2 = track.getChi2OverNDF();
  Fill(mTrackChi2OverNDF[q], chi2);
  Fill(mTrackChi2OverNDF[2], chi2);

  double dca = track.getDCA();
  Fill(mTrackDCA[q], dca);
  Fill(mTrackDCA[2], dca);

  double pdca = track.getPDCAMCH();
  Fill(mTrackPDCA[q], pdca);
  Fill(mTrackPDCA[2], pdca);

  auto rAbs = track.getRAbs();
  Fill(mTrackRAbs[q], rAbs);
  Fill(mTrackRAbs[2], rAbs);

  // Kinematic distributions, both from MCH tracks parameters and from global tracks parameters if MFT is included
  double eta = track.getEta();
  double etaMCH = track.hasMCH() ? track.getEtaMCH() : 0;
  Fill(mTrackEta[q], etaMCH);
  Fill(mTrackEta[2], etaMCH);
  Fill(mTrackEtaGlobal[q], eta);
  Fill(mTrackEtaGlobal[2], eta);
  Fill(mTrackEtaCorr[q], etaMCH, eta);
  Fill(mTrackEtaCorr[2], etaMCH, eta);
  Fill(mTrackDEtaVsEta[q], etaMCH, eta - etaMCH);
  Fill(mTrackDEtaVsEta[2], etaMCH, eta - etaMCH);

  double phi = track.getPhi();
  double phiMCH = track.hasMCH() ? track.getPhiMCH() : 0;
  Fill(mTrackPhi[q], phiMCH);
  Fill(mTrackPhi[2], phiMCH);
  Fill(mTrackPhiGlobal[q], phi);
  Fill(mTrackPhiGlobal[2], phi);
  Fill(mTrackPhiCorr[q], phiMCH, phi);
  Fill(mTrackPhiCorr[2], phiMCH, phi);
  Fill(mTrackDPhiVsPhi[q], phiMCH, phi - phiMCH);
  Fill(mTrackDPhiVsPhi[2], phiMCH, phi - phiMCH);

  Fill(mTrackEtaPhi[q], etaMCH, phiMCH);
  Fill(mTrackEtaPhi[2], etaMCH, phiMCH);
  Fill(mTrackEtaPhiGlobal[q], eta, phi);
  Fill(mTrackEtaPhiGlobal[2], eta, phi);

  double pt = track.getPt();
  double ptMCH = track.hasMCH() ? track.getPtMCH() : 0;
  Fill(mTrackPt[q], ptMCH);
  Fill(mTrackPt[2], ptMCH);
  Fill(mTrackPtGlobal[q], pt);
  Fill(mTrackPtGlobal[2], pt);
  Fill(mTrackPtCorr[q], ptMCH, pt);
  Fill(mTrackPtCorr[2], ptMCH, pt);
  Fill(mTrackDPtVsPt[q], ptMCH, pt - ptMCH);
  Fill(mTrackDPtVsPt[2], ptMCH, pt - ptMCH);

  Fill(mTrackEtaPt[q], etaMCH, ptMCH);
  Fill(mTrackEtaPt[2], etaMCH, ptMCH);
  Fill(mTrackEtaPtGlobal[q], eta, pt);
  Fill(mTrackEtaPtGlobal[2], eta, pt);

  Fill(mTrackPhiPt[q], phiMCH, ptMCH);
  Fill(mTrackPhiPt[2], phiMCH, ptMCH);
  Fill(mTrackPhiPtGlobal[q], phi, pt);
  Fill(mTrackPhiPtGlobal[2], phi, pt);

  double pMCH = track.hasMCH() ? std::fabs(track.getMuonMomentumMCH().P()) : 0;
  Fill(mTrackP[q], pMCH);
  Fill(mTrackP[2], pMCH);

  if (track.getSign() != 0 && ptMCH != 0) {
    double qOverPtMCH = track.getSign() / ptMCH;
    Fill(mTrackQOverPt, qOverPtMCH);
  }

  if (track.getSign() != 0 && pt != 0) {
    double qOverPt = track.getSign() / pt;
    Fill(mTrackQOverPtGlobal, qOverPt);
  }

  o2::mch::TrackParam trackParamAtVertex;
  track.extrapToZMCH(trackParamAtVertex, 0);
  double xVtx = trackParamAtVertex.getNonBendingCoor();
  double yVtx = trackParamAtVertex.getBendingCoor();
  double RVtx = TMath::Sqrt(xVtx * xVtx + yVtx * yVtx);
  double xSlopeVtx = trackParamAtVertex.getNonBendingSlope();
  double ySlopeVtx = trackParamAtVertex.getBendingSlope();
  double phiVtx = TMath::ATan2(yVtx, xVtx) * 180 / TMath::Pi();
  double pVtx = trackParamAtVertex.p();
  double pxVtx = trackParamAtVertex.px();
  double pyVtx = trackParamAtVertex.py();
  double pzVtx = trackParamAtVertex.pz();
  double ptVtx = TMath::Sqrt(pxVtx * pxVtx + pyVtx * pyVtx);
  double thetaVtx = TMath::ATan2(ptVtx, -pzVtx) * 180 / TMath::Pi();
  Fill(mTrackPosAtVertex, xVtx, yVtx);
  Fill(mTrackRPhiAtVertex, phiVtx, RVtx);
  Fill(mTrackRVsThetaAtVertex, thetaVtx, RVtx);
  Fill(mTrackRVsMomAtVertex, pVtx, RVtx);
  Fill(mTrackXSlopeVsPhiAtVertex, phiVtx, xSlopeVtx);
  Fill(mTrackYSlopeVsPhiAtVertex, phiVtx, ySlopeVtx);

  o2::mch::TrackParam trackParamAtAbs = track.getTrackParamAtAbsMCH();
  //track.extrapToZMCH(trackParamAtAbs, o2::quality_control_modules::muon::MuonTrack::sAbsZEnd);
  double xAbs = trackParamAtAbs.getNonBendingCoor();
  double yAbs = trackParamAtAbs.getBendingCoor();
  double RAbs = TMath::Sqrt(xAbs * xAbs + yAbs * yAbs);
  double xSlopeAbs = trackParamAtVertex.getNonBendingSlope();
  double ySlopeAbs = trackParamAtVertex.getBendingSlope();
  double phiAbs = TMath::ATan2(yAbs, xAbs) * 180 / TMath::Pi();
  double pAbs = trackParamAtVertex.p();
  double pxAbs = trackParamAtVertex.px();
  double pyAbs = trackParamAtVertex.py();
  double pzAbs = trackParamAtVertex.pz();
  double ptAbs = TMath::Sqrt(pxAbs * pxAbs + pyAbs * pyAbs);
  double thetaAbs = TMath::ATan2(ptAbs, -pzAbs) * 180 / TMath::Pi();
  Fill(mTrackPosAtAbsorber, xAbs, yAbs);
  Fill(mTrackRPhiAtAbsorber, phiAbs, RAbs);
  Fill(mTrackRVsThetaAtAbsorber, thetaAbs, RAbs);
  Fill(mTrackRVsMomAtAbsorber, pAbs, RAbs);
  Fill(mTrackXSlopeVsPhiAtAbsorber, phiAbs, xSlopeAbs);
  Fill(mTrackYSlopeVsPhiAtAbsorber, phiAbs, ySlopeAbs);

  //Fill(mTrackPosAtAbsorber, track.getXAbsMCH(), track.getYAbsMCH());

  /*
  o2::mch::TrackParam trackParamAtMFT = track.getTrackParamMCH();
  float zMFT = sLastMFTPlaneZ;
  //track.extrapToZMCH(trackParamAtMFT, zMFT);
  double xMCH = 1000000;
  double yMCH = 1000000;
  if (o2::mch::TrackExtrap::extrapToVertexWithoutBranson(trackParamAtMFT, zMFT)) {
    xMCH = trackParamAtMFT.getNonBendingCoor();
    yMCH = trackParamAtMFT.getBendingCoor();
  }
  Fill(mTrackPosAtMFT, xMCH, yMCH);
  */

  o2::mch::TrackParam trackParamAtMID = track.getTrackParamAtMID();
  double xMID = trackParamAtMID.getNonBendingCoor();
  double yMID = trackParamAtMID.getBendingCoor();
  double RMID = TMath::Sqrt(xMID * xMID + yMID * yMID);
  double xSlopeMID = trackParamAtVertex.getNonBendingSlope();
  double ySlopeMID = trackParamAtVertex.getBendingSlope();
  double phiMID = TMath::ATan2(yMID, xMID) * 180 / TMath::Pi();
  double pMID = trackParamAtVertex.p();
  double pxMID = trackParamAtVertex.px();
  double pyMID = trackParamAtVertex.py();
  double pzMID = trackParamAtVertex.pz();
  double ptMID = TMath::Sqrt(pxMID * pxMID + pyMID * pyMID);
  double thetaMID = TMath::ATan2(ptMID, -pzMID) * 180 / TMath::Pi();
  Fill(mTrackPosAtMID, xMID, yMID);
  Fill(mTrackRPhiAtMID, phiMID, RMID);
  Fill(mTrackRVsThetaAtMID, thetaMID, RMID);
  Fill(mTrackRVsMomAtMID, pMID, RMID);
  Fill(mTrackXSlopeVsPhiAtMID, phiMID, xSlopeMID);
  Fill(mTrackYSlopeVsPhiAtMID, phiMID, ySlopeMID);

  //Fill(mTrackPosAtMID, track.getXMid(), track.getYMid());
}

void TrackPlotter::fillHistograms(const o2::globaltracking::RecoContainer& recoCont)
{
  static bool sFirst = true;

  if (sFirst) {
    o2::mch::TrackExtrap::setField();
    sFirst = false;
  }

  if (!o2::mch::TrackExtrap::isFieldON()) {
    o2::mch::TrackExtrap::setField();
  }

  if (mNOrbitsPerTF < 0) {
    mNOrbitsPerTF = o2::base::GRPGeomHelper::instance().getNHBFPerTF();
  }

  if (mBzMFT == 0) {
    auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (field) {
    double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
    mBzMFT = field->getBz(centerMFT);
    }
  }

  mMuonTracks.clear();

  if (mSrc == GID::MCH) {
    auto tracksMCH = recoCont.getMCHTracks();
    int trackID = 0;
    for (auto& t : tracksMCH) {
      mMuonTracks.emplace_back(std::make_pair<MuonTrack, bool>({ &t, trackID, recoCont, mFirstTForbit, mBzMFT }, true));
      trackID += 1;
    }
  }
  if (mSrc == GID::MFTMCH || mSrc == GID::MFTMCHMID) {
    auto tracksFwd = recoCont.getGlobalFwdTracks();
    for (auto& t : tracksFwd) {
      //MuonTrack mt(&t, recoCont, mFirstTForbit, mBzMFT);
      // skip tracks without MID if full matching is requested
      if (mSrc == GID::MFTMCHMID && t.getMIDTrackID() < 0) {
        continue;
      }
      mMuonTracks.emplace_back(std::make_pair<MuonTrack, bool>({ &t, recoCont, mFirstTForbit, mBzMFT }, true));
    }
  }
  if (mSrc == GID::MCHMID) {
    auto tracksMCHMID = recoCont.getMCHMIDMatches();
    for (auto& t : tracksMCHMID) {
      mMuonTracks.emplace_back(std::make_pair<MuonTrack, bool>({ &t, recoCont, mFirstTForbit, mBzMFT }, true));
    }
  }

  int nPos{ 0 };
  int nNeg{ 0 };
  int nTot{ 0 };

  for (auto& t : mMuonTracks) {
    bool ok = true;
    for (auto& cut : mMuonCuts) {
      if (!cut(t.first)) {
        ok = false;
        break;
      }
    }
    t.second = ok;
    if (!t.second) {
      continue;
    }

    nTot += 1;
    if (t.first.getSign() < 0) {
      nNeg += 1;
    } else {
      nPos += 1;
    }
  }

  Fill(mNofTracksPerTF[0], nPos);
  Fill(mNofTracksPerTF[1], nNeg);
  Fill(mNofTracksPerTF[2], nTot);

  for (const auto& mt : mMuonTracks) {
    if (!mt.second) {
      continue;
    }
    fillTrackHistos(mt.first);
  }

  fillTrackPairHistos(mMuonTracks);

  for (auto hinfo : histograms()) {
    TH1* h1 = dynamic_cast<TH1*>(hinfo.object);
    if (h1) {
      normalizePlot(h1);
    }
  }
}

void TrackPlotter::endOfCycle()
{
  for (auto hinfo : histograms()) {
    TH1DRatio* h1 = dynamic_cast<TH1DRatio*>(hinfo.object);
    if (h1) {
      h1->update();
    } else {
      TH2DRatio* h2 = dynamic_cast<TH2DRatio*>(hinfo.object);
      if (h2) {
        h2->update();
      }
    }
  }
}

} // namespace o2::quality_control_modules::muon
