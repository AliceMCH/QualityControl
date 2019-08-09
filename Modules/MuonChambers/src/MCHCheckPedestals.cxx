// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   MCHCheckPedestals.cxx
/// \author Andrea Ferrero
///

#include "MuonChambers/MCHCheckPedestals.h"

// ROOT
#include <fairlogger/Logger.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveText.h>

using namespace std;

namespace o2::quality_control_modules::muonchambers
{

MCHCheckPedestals::MCHCheckPedestals() : minMCHpedestal(50.f), maxMCHpedestal(100.f)
{
}

MCHCheckPedestals::~MCHCheckPedestals() {}

void MCHCheckPedestals::configure(std::string)
{
  // if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic) {
  //   minTOFrawTime = 150.; //ns
  //   maxTOFrawTime = 250.; //ns
  // }
}

Quality MCHCheckPedestals::check(const MonitorObject* mo)
{
  Quality result = Quality::Null;

  // const Double_t binWidthTOFrawTime = 2.44;

  // if ((histname.EndsWith("RawsTime")) || (histname.Contains("RawsTime") && suffixTrgCl)) {
  if (mo->getName().find("QcMuonChambers_Pedestals") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    if( !h ) return result;

    if (h->GetEntries() == 0) {
      result = Quality::Medium;
      // flag = AliQAv1::kWARNING;
    } else {
      int nbinsx = 6;//h->GetXaxis()->GetNbins();
      int nbinsy = h->GetYaxis()->GetNbins();
      int nbad = 0;
      for(int i = 1; i <= nbinsx; i++) {
        for(int j = 1; j <= nbinsy; j++) {
          Float_t ped = h->GetBinContent(i, j);
          if( ped < minMCHpedestal || ped > maxMCHpedestal ) nbad += 1;
        }
      }
      if( nbad < 1 ) result = Quality::Good;
      else result = Quality::Bad;
    }
  }
  return result;
}

std::string MCHCheckPedestals::getAcceptedType() { return "TH1"; }

void MCHCheckPedestals::beautify(MonitorObject* mo, Quality checkResult)
{
  if (mo->getName().find("QcMuonChambers_Pedestals") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetDrawOption("colz");
    TPaveText* msg = new TPaveText(0.1, 0.9, 0.9, 0.95, "NDC");
    h->GetListOfFunctions()->Add(msg);
    msg->SetName(Form("%s_msg", mo->GetName()));

    if (checkResult == Quality::Good) {
      msg->Clear();
      msg->AddText("All pedestals within limits: OK!!!");
      msg->SetFillColor(kGreen);
      //
      h->SetFillColor(kGreen);
    } else if (checkResult == Quality::Bad) {
      LOG(INFO) << "Quality::Bad, setting to red";
      //
      msg->Clear();
      msg->AddText("Call MCH on-call.");
      msg->SetFillColor(kRed);
      //
      h->SetFillColor(kRed);
    } else if (checkResult == Quality::Medium) {
      LOG(INFO) << "Quality::medium, setting to orange";
      //
      msg->Clear();
      msg->AddText("No entries. If MCH in the run, check MCH TWiki");
      msg->SetFillColor(kYellow);
      // text->Clear();
      // text->AddText(Form("Raw time peak/total integral = %5.2f%%", peakIntegral * 100. / totIntegral));
      // text->AddText(Form("Mean = %5.2f ns", timeMean));
      // text->AddText(Form("Allowed range: %3.0f-%3.0f ns", minTOFrawTime, maxTOFrawTime));
      // text->AddText("If multiple peaks, check filling scheme");
      // text->AddText("See TOF TWiki.");
      // text->SetFillColor(kYellow);
      //
      h->SetFillColor(kOrange);
    }
    h->SetLineColor(kBlack);
  }
}

} // namespace o2::quality_control_modules::tof