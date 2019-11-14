///
/// \file   RawDataProcessor.h
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_MUONCHAMBERS_PHYSICSDATAPROCESSOR_H
#define QC_MODULE_MUONCHAMBERS_PHYSICSDATAPROCESSOR_H

#include "QualityControl/TaskInterface.h"
#include "MCH/MuonChambersMapping.h"
#include "MCH/MuonChambersDataDecoder.h"

class TH1F;
class TH2F;

using namespace o2::quality_control::core;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{


/// \brief Example Quality Control DPL Task
/// It is final because there is no reason to derive from it. Just remove it if needed.
/// \author Barthelemy von Haller
/// \author Piotr Konopka
class PhysicsDataProcessor /*final*/ : public TaskInterface // todo add back the "final" when doxygen is fixed
{
 public:
  /// \brief Constructor
  PhysicsDataProcessor();
  /// Destructor
  ~PhysicsDataProcessor() override;

  // Definition of the methods for the template method pattern
  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(Activity& activity) override;
  void startOfCycle() override;
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(Activity& activity) override;
  void reset() override;

 private:
  int count;
  MuonChambersDataDecoder mDecoder;
  uint64_t nhits[24][40][64];

  TH2F* mHistogramNhits[24];
  TH1F* mHistogramADCamplitude[24];
  std::map<int, TH1F*> mHistogramADCamplitudeDE;
  std::map<int, TH2F*> mHistogramNhitsDE;
  std::map<int, TH2F*> mHistogramNhitsHighAmplDE;
};

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

#endif // QC_MODULE_MUONCHAMBERS_PHYSICSDATAPROCESSOR_H