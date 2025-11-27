#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <set>

#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void SetStyle() {
  gStyle->SetLineScalePS(1.2);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTitleOffset(1.25, "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetPalette(1);
  // TColor::InvertPalette();
}

void setHist(TH2D* h)
{
  h->GetYaxis()->SetTickSize(0.01);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetMaxDigits(5);
  h->LabelsOption("av");
  h->SetContour(20);
}

void filter()
{
  double threshold = 0.8;
  SetStyle();

  auto* fdigits = new TFile("/afs/cern.ch/work/s/sgaretti/rejectList/PbPb/analyse/analyse/figures/digits.root", "r");

  std::map<int, std::vector<int>> deIdsPerChamber{};
  deIdsPerChamber[1] = {100, 101, 102, 103};
  deIdsPerChamber[2] = {200, 201, 202, 203};
  deIdsPerChamber[3] = {300, 301, 302, 303};
  deIdsPerChamber[4] = {400, 401, 402, 403};
  deIdsPerChamber[5] = {500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517};
  deIdsPerChamber[6] = {600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617};
  deIdsPerChamber[7] = {700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725};
  deIdsPerChamber[8] = {800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825};
  deIdsPerChamber[9] = {900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925};
  deIdsPerChamber[10] = {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025};

  gSystem->Exec("if [[ ! -d figures_digits ]]; then mkdir figures; fi");
  auto* fout = new TFile("./figures/filter.root", "recreate");
  auto* l_out = new TList();

  int nRuns = 0;

  std::map<std::string, TH2D*> mapTrendsPerDE{};

  std::fstream fs("banned_ds.list", std::fstream::out);

  for (const auto& [chId, deIds] : deIdsPerChamber) {
    for (auto deId : deIds) {
      for (auto side : {"B", "NB"}) {
        auto* hDigits = dynamic_cast<TH2D*>(fdigits->Get(Form("hTrendDSperDE%d%s", deId, side)));
        if (!hDigits) {
          std::cerr << "Warning: histogram not found for DE " << deId << " side " << side << std::endl;
          continue;
        }
        hDigits->SetName(Form("hTrendDSperDE%d", deId));
        auto* hFilter = dynamic_cast<TH2D*>(hDigits->Clone(Form("hTrendDSperDE%d%s_filter", deId, side)));
        if (nRuns == 0) nRuns = hDigits->GetNbinsX();
        for (int ir = 1; ir <= hDigits->GetNbinsX(); ++ir) {
          int run = std::stoi(hDigits->GetXaxis()->GetBinLabel(ir));
          for (int ids = 1; ids <= hDigits->GetNbinsY(); ++ids) {
            int ds = trunc(hDigits->GetYaxis()->GetBinLowEdge(ids));
            double ndig = hDigits->GetBinContent(ir, ids);
            bool isBanned = false;
            /* if (ndig < 0.5)                 { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig < 0.66 && deId == 300) { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig < 0.66 && deId == 301) { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig < 0.66 && deId == 302) { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig < 0.66 && deId == 400) { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig < 0.66 && deId == 401) { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig < 0.66 && deId == 403) { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }*/
            
            if (ndig < threshold)                 { hFilter->SetBinContent(ir, ids, 1.5);  isBanned = true; }
            if (ndig >= threshold)                { hFilter->SetBinContent(ir, ids, ndig); isBanned = false; }

            /* if (isBannedByRect)             { hFilter->SetBinContent(ir, ids, 2.);   isBanned = true; }
            if (isAllowedByRect)            { hFilter->SetBinContent(ir, ids, ndig); isBanned = false; } */
            if (isBanned) {
              fs << Form("%d %d\n", run, ds);
            }
          }
        }
        mapTrendsPerDE[hFilter->GetName()] = hFilter;
        l_out->Add(hFilter);
      }
    }
  }

  fs.close();

  // draw/store distributions

  auto* ctrend = new TCanvas("ctrend", "ctrend", 5000, 3000);
  ctrend->cd();

  auto* padLeft = new TPad("padLeft", "padLeft", 0., 0., 0.5, 1.);
  padLeft->SetLeftMargin(0.08);
  padLeft->SetRightMargin(0.);
  padLeft->SetTopMargin(0.07);
  padLeft->SetBottomMargin(0.11);
  padLeft->Draw();

  auto* padRight = new TPad("padRight", "padRight", 0.5, 0., 1., 1.);
  padRight->SetLeftMargin(0.08);
  padRight->SetRightMargin(0.15);
  padRight->SetRightMargin(0.005);
  padRight->SetTopMargin(0.07);
  padRight->SetBottomMargin(0.11);
  padRight->Draw();

  for (auto& [chId, deIds] : deIdsPerChamber) {
    for (auto deId : deIds) {
      auto* hTrendDSperDE_B = mapTrendsPerDE.at(Form("hTrendDSperDE%dB_filter", deId));
      auto* hTrendDSperDE_NB = mapTrendsPerDE.at(Form("hTrendDSperDE%dNB_filter", deId));
      setHist(hTrendDSperDE_B);
      setHist(hTrendDSperDE_NB);
      ctrend->cd();
      padLeft->cd();
      hTrendDSperDE_B->SetTitle(Form("DE%d Bending", deId));
      hTrendDSperDE_B->GetZaxis()->SetRangeUser(0., 2.);
      hTrendDSperDE_B->Draw("col");
      ctrend->cd();
      padRight->cd();
      hTrendDSperDE_NB->SetTitle(Form("DE%d Non-Bending", deId));
      hTrendDSperDE_NB->GetZaxis()->SetRangeUser(0., 2.);
      hTrendDSperDE_NB->Draw("col");
      ctrend->cd();
      if (deId == 100)
        ctrend->Print("./figures/trending_filter_ds_per_de.pdf(", Form("Title: DE%d", deId));
      else if (deId == 1025)
        ctrend->Print("./figures/trending_filter_ds_per_de.pdf)", Form("Title: DE%d", deId));
      else
        ctrend->Print("./figures/trending_filter_ds_per_de.pdf", Form("Title: DE%d", deId));
    }
  }

  fout->cd();
  l_out->Write();
}
