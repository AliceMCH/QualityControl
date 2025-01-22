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

#include "MUONCommon/Helpers.h"
#include <TLine.h>
#include <TAxis.h>
#include <TH1.h>
#include <TList.h>
#include <TPolyLine.h>
#include <regex>
#include "QualityControl/ObjectsManager.h"

#include <cmath>

namespace o2::quality_control_modules::muon
{
template <>
std::string getConfigurationParameter<std::string>(o2::quality_control::core::CustomParameters customParameters, std::string parName, const std::string defaultValue)
{
  std::string result = defaultValue;
  auto parOpt = customParameters.atOptional(parName);
  if (parOpt.has_value()) {
    result = parOpt.value();
  }
  return result;
}

template <>
std::string getConfigurationParameter<std::string>(o2::quality_control::core::CustomParameters customParameters, std::string parName, const std::string defaultValue, const o2::quality_control::core::Activity& activity)
{
  auto parOpt = customParameters.atOptional(parName, activity);
  if (parOpt.has_value()) {
    std::string result = parOpt.value();
    return result;
  }
  return getConfigurationParameter<std::string>(customParameters, parName, defaultValue);
}

template <>
bool getConfigurationParameter<bool>(o2::quality_control::core::CustomParameters customParameters, std::string parName, const bool defaultValue)
{
  bool result = defaultValue;
  auto parOpt = customParameters.atOptional(parName);
  if (parOpt.has_value()) {
    std::string value = parOpt.value();
    std::transform(value.begin(), value.end(), value.begin(), ::toupper);
    if (value == "TRUE" || value == "YES" || value == "1") {
      return true;
    }
    if (value == "FALSE" || value == "NO" || value == "0") {
      return false;
    }
    throw std::invalid_argument(std::string("error parsing boolean configurable parameter: key=") + parName + " value=" + value);
  }
  return result;
}

template <>
bool getConfigurationParameter<bool>(o2::quality_control::core::CustomParameters customParameters, std::string parName, const bool defaultValue, const o2::quality_control::core::Activity& activity)
{
  auto parOpt = customParameters.atOptional(parName, activity);
  if (parOpt.has_value()) {
    std::string value = parOpt.value();
    std::transform(value.begin(), value.end(), value.begin(), ::toupper);
    if (value == "TRUE" || value == "YES" || value == "1") {
      return true;
    }
    if (value == "FALSE" || value == "NO" || value == "0") {
      return false;
    }
    throw std::invalid_argument(std::string("error parsing boolean configurable parameter: key=") + parName + " value=" + value);
  }
  return getConfigurationParameter<bool>(customParameters, parName, defaultValue);
}

//_________________________________________________________________________________________

std::vector<double> makeLogBinning(double min, double max, int nbins)
{
  auto logMin = std::log10(min);
  auto logMax = std::log10(max);
  auto binWidth = (logMax - logMin) / nbins;
  std::vector<double> bins(nbins + 1);
  for (int i = 0; i <= nbins; i++) {
    bins[i] = std::pow(10, logMin + i * binWidth);
  }
  return bins;
}

//_________________________________________________________________________________________

TLine* addHorizontalLine(TH1& histo, double y,
                         int lineColor, int lineStyle,
                         int lineWidth)
{
  auto nbins = histo.GetXaxis()->GetNbins();

  TLine* line = new TLine(histo.GetBinLowEdge(1), y, histo.GetBinLowEdge(nbins) + histo.GetBinWidth(nbins), y);
  line->SetLineColor(lineColor);
  line->SetLineStyle(lineStyle);
  line->SetLineWidth(lineWidth);
  histo.GetListOfFunctions()->Add(line);
  return line;
}

void markBunchCrossing(TH1& histo,
                       gsl::span<int> bunchCrossings)
{
  for (auto b : bunchCrossings) {
    addVerticalLine(histo, b * 1.0, 1, 10, 1);
  }
}

TLine* addVerticalLine(TH1& histo, double x,
                       int lineColor, int lineStyle,
                       int lineWidth)
{
  double max = histo.GetBinContent(histo.GetMaximumBin());
  TLine* line = new TLine(x, histo.GetMinimum(),
                          x, max * 1.05);
  line->SetLineColor(lineColor);
  line->SetLineStyle(lineStyle);
  line->SetLineWidth(lineWidth);
  histo.GetListOfFunctions()->Add(line);
  return line;
}

/// Add a marker to an histogram at a given position
/// The marker is draw with a TPolyLine such that it scales nicely with the size of the pad
/// The default dimensions of the marker are
/// * horizontal: 1/20 of the X-axis range
/// * vertical: 1/10 of the histogram values range
/// Parameters:
/// * histo: the histogram to which the marker is associated
/// * x0, y0: coordinates of the tip of the marker
/// * color: ROOT index of the marker fill color
/// * markerSize: overall scaling factor for the marker dimensions
/// * logx, logy: wether the X or Y axis are in logarithmic scale
TPolyLine* _addMarker(TH1& histo, double x, double y, int markerColor, float markerSize, bool logx, bool logy)
{
  double x0 = x;
  double y0 = y;
  double xmin = logx ? std::log10(histo.GetXaxis()->GetXmin()) : histo.GetXaxis()->GetXmin();
  double xmax = logx ? std::log10(histo.GetXaxis()->GetXmax()) : histo.GetXaxis()->GetXmax();
  double xSize = (xmax - xmin) / 20;
  double ySize = (histo.GetMaximum() - histo.GetMinimum()) / 10;
  double x1 = logx ? std::pow(10, std::log10(x0) - xSize * markerSize / 2) : x0 - xSize * markerSize / 2;
  double x2 = logx ? std::pow(10, std::log10(x0) + xSize * markerSize / 2) : x0 + xSize * markerSize / 2;
  double xMarker[4] = { x0, x1, x2, x0 };
  double yMarker[4] = { y0, y0 + ySize * markerSize, y0 + ySize * markerSize, y0 };
  auto* m = new TPolyLine(4, xMarker, yMarker);
  m->SetNDC(kFALSE);
  m->SetFillColor(markerColor);
  m->SetOption("f");
  m->SetLineWidth(0);
  histo.GetListOfFunctions()->Add(m);
  return m;
}

// remove all elements of class c
// from histo->GetListOfFunctions()
void cleanup(TH1& histo, const char* classname)
{
  TList* elements = histo.GetListOfFunctions();
  TIter next(elements);
  TObject* obj;
  std::vector<TObject*> toBeRemoved;
  while ((obj = next())) {
    if (strcmp(obj->ClassName(), classname) == 0) {
      toBeRemoved.push_back(obj);
    }
  }
  for (auto o : toBeRemoved) {
    elements->Remove(o);
  }
}
} // namespace o2::quality_control_modules::muon
