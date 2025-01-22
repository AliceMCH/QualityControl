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

#ifndef QC_MODULE_MUON_COMMON_HELPERS_H
#define QC_MODULE_MUON_COMMON_HELPERS_H

#include "QualityControl/CustomParameters.h"
#include <gsl/span>
#include <sstream>

class TH1;
class TLine;
class TPolyLine;

namespace o2::quality_control_modules::muon
{
template <class T>
T getConfigurationParameter(o2::quality_control::core::CustomParameters customParameters, std::string parName, const T defaultValue)
{
  T result = defaultValue;
  auto parOpt = customParameters.atOptional(parName);
  if (parOpt.has_value()) {
    std::stringstream ss(parOpt.value());
    ss >> result;
  }
  return result;
}

template <class T>
T getConfigurationParameter(o2::quality_control::core::CustomParameters customParameters, std::string parName, const T defaultValue, const o2::quality_control::core::Activity& activity)
{
  auto parOpt = customParameters.atOptional(parName, activity);
  if (parOpt.has_value()) {
    T result;
    std::stringstream ss(parOpt.value());
    ss >> result;
    return result;
  }
  return getConfigurationParameter<T>(customParameters, parName, defaultValue);
}

// create an array of bins in a given range with log10 spacing
std::vector<double> makeLogBinning(double xmin, double xmax, int nbins);

// add lines to an histogram with specified color, style and width, using ROOT's TAttLine conventions
TLine* addHorizontalLine(TH1& histo, double y,
                         int lineColor = 1, int lineStyle = 10,
                         int lineWidth = 1);
TLine* addVerticalLine(TH1& histo, double x,
                       int lineColor = 1, int lineStyle = 10,
                       int lineWidth = 1);

/// Add a marker to an histogram at a given position
/// The marker is draw with a TPolyLine such that it scales nicely with the size of the pad
/// The default dimensions of the marker are
/// * horizontal: 1/20 of the X-axis range
/// * vertical: 1/10 of the histogram values range
/// Parameters:
/// * histo: the histogram to which the marker is added
/// * x, y: coordinates of the tip of the marker
/// * color: ROOT index of the marker fill color
/// * markerSize: overall scaling factor for the marker dimensions
/// * logx, logy: wether the X or Y axis are in logarithmic scale
TPolyLine* _addMarker(TH1& histo, double x, double y, int markerColor, float markerSize, bool logx, bool logy);

void cleanup(TH1& histo, const char* classname);
void markBunchCrossing(TH1& histo,
                       gsl::span<int> bunchCrossings);

} // namespace o2::quality_control_modules::muon

#endif
