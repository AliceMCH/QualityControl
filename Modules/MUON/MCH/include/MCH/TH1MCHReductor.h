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
/// \file   TH1MCHReductor.h
/// \author Piotr Konopka, Sebastien Perrin
///
#ifndef QUALITYCONTROL_TH1MCHREDUCTOR_H
#define QUALITYCONTROL_TH1MCHREDUCTOR_H

#include "QualityControl/Reductor.h"

namespace o2::quality_control_modules::muonchambers
{

/// \brief A Reductor which obtains the most popular characteristics of TH1.
///
/// A Reductor which obtains the most popular characteristics of TH1.

class TH1MCHReductor : public quality_control::postprocessing::Reductor
{
 public:
  TH1MCHReductor() = default;
  ~TH1MCHReductor() = default;

  void* getBranchAddress() override;
  const char* getBranchLeafList() override;
  void update(TObject* obj) override;

 private:
  struct {
    union {
      struct {
        Double_t occ500;
        Double_t occ501;
        Double_t occ502;
        Double_t occ503;
        Double_t occ504;
        Double_t occ505;
        Double_t occ506;
        Double_t occ507;
        Double_t occ508;
        Double_t occ509;
        Double_t occ510;
        Double_t occ511;
        Double_t occ512;
        Double_t occ513;
        Double_t occ514;
        Double_t occ515;
        Double_t occ516;
        Double_t occ517;
        Double_t occ600;
        Double_t occ601;
        Double_t occ602;
        Double_t occ603;
        Double_t occ604;
        Double_t occ605;
        Double_t occ606;
        Double_t occ607;
        Double_t occ608;
        Double_t occ609;
        Double_t occ610;
        Double_t occ611;
        Double_t occ612;
        Double_t occ613;
        Double_t occ614;
        Double_t occ615;
        Double_t occ616;
        Double_t occ617;
        Double_t occ700;
        Double_t occ701;
        Double_t occ702;
        Double_t occ703;
        Double_t occ704;
        Double_t occ705;
        Double_t occ706;
        Double_t occ707;
        Double_t occ708;
        Double_t occ709;
        Double_t occ710;
        Double_t occ711;
        Double_t occ712;
        Double_t occ713;
        Double_t occ714;
        Double_t occ715;
        Double_t occ716;
        Double_t occ717;
        Double_t occ718;
        Double_t occ719;
        Double_t occ720;
        Double_t occ721;
        Double_t occ722;
        Double_t occ723;
        Double_t occ724;
        Double_t occ725;
        Double_t occ800;
        Double_t occ801;
        Double_t occ802;
        Double_t occ803;
        Double_t occ804;
        Double_t occ805;
        Double_t occ806;
        Double_t occ807;
        Double_t occ808;
        Double_t occ809;
        Double_t occ810;
        Double_t occ811;
        Double_t occ812;
        Double_t occ813;
        Double_t occ814;
        Double_t occ815;
        Double_t occ816;
        Double_t occ817;
        Double_t occ818;
        Double_t occ819;
        Double_t occ820;
        Double_t occ821;
        Double_t occ822;
        Double_t occ823;
        Double_t occ824;
        Double_t occ825;
        Double_t occ900;
        Double_t occ901;
        Double_t occ902;
        Double_t occ903;
        Double_t occ904;
        Double_t occ905;
        Double_t occ906;
        Double_t occ907;
        Double_t occ908;
        Double_t occ909;
        Double_t occ910;
        Double_t occ911;
        Double_t occ912;
        Double_t occ913;
        Double_t occ914;
        Double_t occ915;
        Double_t occ916;
        Double_t occ917;
        Double_t occ918;
        Double_t occ919;
        Double_t occ920;
        Double_t occ921;
        Double_t occ922;
        Double_t occ923;
        Double_t occ924;
        Double_t occ925;
        Double_t occ1000;
        Double_t occ1001;
        Double_t occ1002;
        Double_t occ1003;
        Double_t occ1004;
        Double_t occ1005;
        Double_t occ1006;
        Double_t occ1007;
        Double_t occ1008;
        Double_t occ1009;
        Double_t occ1010;
        Double_t occ1011;
        Double_t occ1012;
        Double_t occ1013;
        Double_t occ1014;
        Double_t occ1015;
        Double_t occ1016;
        Double_t occ1017;
        Double_t occ1018;
        Double_t occ1019;
        Double_t occ1020;
        Double_t occ1021;
        Double_t occ1022;
        Double_t occ1023;
        Double_t occ1024;
        Double_t occ1025;
      } named;
      Double_t indiv[140];
    } indiv_occs;
    union {
      struct {
        Double_t occ5I;
        Double_t occ5O;
        Double_t occ6I;
        Double_t occ6O;
        Double_t occ7I;
        Double_t occ7O;
        Double_t occ8I;
        Double_t occ8O;
        Double_t occ9I;
        Double_t occ9O;
        Double_t occ10I;
        Double_t occ10O;
      } named;
      Double_t halfch[12];
    } halfch_occs;
    Double_t entries;
  } mStats;
  int detCH5I[9] = { 500, 501, 502, 503, 504, 514, 515, 516, 517 };
  int detCH5O[9] = { 505, 506, 507, 508, 509, 510, 511, 512, 513 };
  int detCH6I[9] = { 600, 601, 602, 603, 604, 614, 615, 616, 617 };
  int detCH6O[9] = { 605, 606, 607, 608, 609, 610, 611, 612, 613 };
  int detCH7I[13] = { 700, 701, 702, 703, 704, 705, 706, 720, 721, 722, 723, 724, 725 };
  int detCH7O[13] = { 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719 };
  int detCH8I[13] = { 800, 801, 802, 803, 804, 805, 806, 820, 821, 822, 823, 824, 825 };
  int detCH8O[13] = { 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819 };
  int detCH9I[13] = { 900, 901, 902, 903, 904, 905, 906, 920, 921, 922, 923, 924, 925 };
  int detCH9O[13] = { 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919 };
  int detCH10I[13] = { 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1020, 1021, 1022, 1023, 1024, 1025 };
  int detCH10O[13] = { 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019 };
};

} // namespace o2::quality_control_modules::muonchambers

#endif //QUALITYCONTROL_TH1mchREDUCTOR_H