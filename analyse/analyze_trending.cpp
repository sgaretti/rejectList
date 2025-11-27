#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TList.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "CommonConstants/LHCConstants.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingInterface/CathodeSegmentationCInterface.h"
#include "MCHMappingInterface/Segmentation.h"

#include "MCHGeometryTransformer/Transformations.h"

#include "MCHConstants/DetectionElements.h"
#include "MCHGlobalMapping/ChannelCode.h"
#include "MCHGlobalMapping/DsIndex.h"
#include "MCHConditions/Plane.h"
#include "MCHGlobalMapping/Mapper.h"

#include "MCHRawElecMap/Mapper.h"

#include "DetectorsBase/GeometryManager.h"
#endif

std::map<int, TString> periods;

void getList(const std::string& fname, std::vector<std::string>& list)
{
  std::ifstream rfile;
  rfile.open(fname);
  std::string line;
  while (std::getline(rfile, line)) {
    std::istringstream iss(line);
    std::string l;
    iss >> l;
    list.push_back(l);
  }
  rfile.close();
}

void getList(const std::string& fname, std::vector<int>& list)
{
  std::ifstream rfile;
  rfile.open(fname);
  std::string line;
  while (std::getline(rfile, line)) {
    std::istringstream iss(line);
    std::string l;
    iss >> l;
    list.push_back(std::atoi(l.c_str()));
  }
  rfile.close();
}

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

void normalize(TH2D* h)
{
  auto* axisRuns = h->GetXaxis();
  int nRuns = h->GetNbinsX();
  int nDS = h->GetNbinsY();
  double firstEdge = axisRuns->GetBinLowEdge(1);
  //int binRuns1 = axisRuns->FindFixBin(firstEdge);
  int binRuns1 = axisRuns->FindFixBin("560229");
  int binRuns2 = axisRuns->FindFixBin("560402");
  int binRuns3 = axisRuns->FindFixBin("559544");
  int binRuns4 = axisRuns->FindFixBin("560163");

  //int nRuns = h->GetNbinsX();
  //int nDS = h->GetNbinsY();
printf("\n runs found, nRuns: %d \n", nRuns);
printf("binRuns1=%d, binRuns2=%d, binRuns3=%d, binRuns4=%d\n", binRuns1, binRuns2, binRuns3, binRuns4);
printf("Sample bin content: binX=1, binY=1 -> %f\n", h->GetBinContent(1,1));

  auto normalizePart = [&h, nDS, nRuns] (int bin1, int bin2)
  {
    // normalize by 'max'
//printf("\n STD norm by max start \n");
    std::vector<double> nMaxPerDS(nDS, -999.);
    for (int ids = 1; ids <= nDS; ++ids) {
//printf("\n Processing DS: %d \n", ids);
      std::set<double> nPerRun{};
//printf("\n Creating set  \n");
      for (int ir = bin1; ir <= bin2; ++ir) {
       //printf("\n ENTERING IN THE LOOP processing DS bin with value: %d \n", ir);
        double n = h->GetBinContent(ir, ids);
        //printf("\n Processing DS bin with value: %f \n", n);
        nPerRun.insert(n);
//printf("\n Processing DS bin: %d \n", ir);
      }
      auto it = nPerRun.rbegin();
      --it;
      --it;
      double maxSecond = *it;
      nMaxPerDS[ids-1] = maxSecond;
//printf("DS %d: maxSecond = %f\n", ids, maxSecond);
    }
//printf("\n 1st norm loop done \n");
    for (int ir = bin1; ir <= bin2; ++ir) {
      for (int ids = 1; ids <= nDS; ++ids) {
        double n = h->GetBinContent(ir, ids);
        double nMax = nMaxPerDS[ids - 1] > 0. ? nMaxPerDS[ids - 1] : 1.;
        double nNorm = n / nMax;
        h->SetBinContent(ir, ids, nNorm);
//printf("Normalized value DS=%d, run=%d: %f\n", ids, ir, nNorm);

      }
    }
    printf("\n 2nd norm loop done \n");
  };

  normalizePart(binRuns1, binRuns2);
  //normalizePart(1, nRuns);
  normalizePart(binRuns3, binRuns4);
 printf("\n STD norm done \n");
}

// 0 -> digits
// 1 -> clusters
void analyze_trending(int flag=0)
{ std::cout << "entering in the function" << endl;
  std::string type = "";
  std::string hname = "";
  if (flag == 0) {
    type = "digits";
    hname = "hDigits";
  }
  if (flag == 1) {
    type = "clusters";
    hname = "hClusters";
  }
  if (flag == 2) {
    type = "digits";
    hname = "hDigitsRat";
  }
  if (flag == 3) {
    type = "digits";
    hname = "hDigitsOrb";
  }

  SetStyle();

  // MCH mappings
  std::cout << "MCH mapping" << endl;
  std::map<int, std::vector<int>> deIdsPerChamber{};
  deIdsPerChamber[1] = {100, 101, 102, 103};
  deIdsPerChamber[2] = {200, 201, 202, 203};
  //deIdsPerChamber[3] = {300, 301, 302, 303};
  //deIdsPerChamber[4] = {400, 401, 402, 403};
  //deIdsPerChamber[5] = {500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517};
  //deIdsPerChamber[6] = {600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617};
  //deIdsPerChamber[7] = {700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725};
  //deIdsPerChamber[8] = {800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825};
  //deIdsPerChamber[9] = {900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925};
  //deIdsPerChamber[10] = {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025};

  // input file list
  std::vector<std::string> flist;
  std::vector<int> rlist;

  getList("input_hists.txt", flist);
  getList("runs.list", rlist);
 // std::cout << "run list acquired and size: " << rlist.size() << endl;

  auto nRuns = static_cast<int>(rlist.size());

  // draw/store distributions

  auto* fout = new TFile(Form("./figures/%s.root", type.c_str()), "recreate");
  auto* l_out = new TList();

  std::map<std::string, TH2D*> mapTrendsPerDE{};
  //printf("entering figures loop");

  for (auto& [chId, deIds] : deIdsPerChamber) {
    for (auto deId : deIds) {
	cout << "ID DE: " << deId << endl;
      for (auto planeId : {0, 1}) {
        auto plane = static_cast<o2::mch::dcs::Plane>(planeId);
        o2::mch::dcs::Cathode cath{deId, plane};
        auto dsIds = o2::mch::dcs::getDsIndices({cath});
        int nDS = dsIds.size();
        int first = *dsIds.begin();
        int last = *dsIds.rbegin();
        //printf("DE %d %s: nds=%d first=%d last=%d\n",
          //deId, planeId == 0 ? "bending" : "non-bending", nDS, first, last);
        fflush(stdout);
        auto* hTrendDSperDE = new TH2D(
          Form("hTrendDSperDE%d%s", deId, planeId == 0 ? "B" : "NB"),
          Form("Number of %s per dual sampa: DE%d %s;;", type.c_str(), deId, planeId == 0 ? "Bending" : "Non-Bending"),
          nRuns, 0, nRuns,
          nDS, first, last+1);
        l_out->Add(hTrendDSperDE);
        mapTrendsPerDE[hTrendDSperDE->GetName()] = hTrendDSperDE;
        for (int i = 0; i < nRuns; ++i) {
          int r = rlist[i];
          hTrendDSperDE->GetXaxis()->SetBinLabel(i+1, Form("%d", r));
        }
      }
    }
  }

  // collect clusters per run
//printf("entering in the cluster loop, nRuns: %d\n", nRuns);
for (int i = 0; i < nRuns; ++i) {
    auto r = rlist[i];
    const auto& fn = flist[i];
    printf("%s\n", fn.c_str());
    auto* file = new TFile(fn.c_str(), "r");

    TH1* hNPerDS = nullptr;
    if (flag == 0) hNPerDS = dynamic_cast<TH1D*>(file->Get(Form("%s", hname.c_str()))); // hDigits
    if (flag == 1) hNPerDS = dynamic_cast<TH1F*>(file->Get(Form("%s", hname.c_str()))); // hClusters
    if (flag == 2) hNPerDS = dynamic_cast<TH1F*>(file->Get(Form("%s", hname.c_str()))); // hDigitsRat
    if (flag == 3) hNPerDS = dynamic_cast<TH1D*>(file->Get(Form("%s", hname.c_str()))); // hDigitsOrb

    printf("flag: %d\n", flag);
    int index = 0;

    for (auto& [chId, deIds] : deIdsPerChamber) {
      std::cout << "first loop" << endl;
      for (auto deId : deIds) {
        for (auto planeId : {0, 1}) {
          auto plane = static_cast<o2::mch::dcs::Plane>(planeId);
          o2::mch::dcs::Cathode cath{deId, plane};
          auto dsIds = o2::mch::dcs::getDsIndices({cath});
          int nDS = dsIds.size();
          int first = *dsIds.begin();
          int last = *dsIds.rbegin();
          auto* hTrendDSperDE =
            mapTrendsPerDE.at(Form("hTrendDSperDE%d%s", deId, planeId == 0 ? "B" : "NB"));
          for (int ds = first; ds <= last; ++ds) {
            double n = hNPerDS->GetBinContent(ds+1);
            int dsLocBin = hTrendDSperDE->GetYaxis()->FindFixBin(ds);
            hTrendDSperDE->SetBinContent(i+1, dsLocBin, n);
          }
        }
      }
	index ++;
        printf("INDEX: %d\n", index);
    }
    file->Close();
    delete file;
    //exit(0);
  } 

/*for (int i = 0; i < nRuns; ++i) {
    auto r = rlist[i];
    const auto& fn = flist[i];
    printf("Processing run %d: %s\n", r, fn.c_str());

    // Apri il file
    TFile* file = TFile::Open(fn.c_str(), "READ");
    if (!file || file->IsZombie()) {
        printf("ERRORE: impossibile aprire %s\n", fn.c_str());
        delete file; // se non è nullptr
        continue;
    }

    // Carica l'istogramma
    TH1* hTmp = nullptr;
    if (flag == 0) hTmp = dynamic_cast<TH1D*>(file->Get(hname.c_str())); // hDigits
    if (flag == 1) hTmp = dynamic_cast<TH1F*>(file->Get(hname.c_str())); // hClusters
    if (flag == 2) hTmp = dynamic_cast<TH1F*>(file->Get(hname.c_str())); // hDigitsRat
    if (flag == 3) hTmp = dynamic_cast<TH1D*>(file->Get(hname.c_str())); // hDigitsOrb

    if (!hTmp) {
        printf("ERRORE: istogramma %s non trovato in %s\n", hname.c_str(), fn.c_str());
        file->Close();
        delete file;
        continue;
    }

    // Clona l'istogramma in memoria e scollega dal TFile
    TH1* hNPerDS = dynamic_cast<TH1*>(hTmp->Clone());
    hNPerDS->SetDirectory(nullptr);

    // Numero di bin
    int nBins = hNPerDS->GetNbinsX();

    // Loop su camere / DE / piani
    for (auto& [chId, deIds] : deIdsPerChamber) {
        for (auto deId : deIds) {
            for (auto planeId : {0, 1}) {
                auto plane = static_cast<o2::mch::dcs::Plane>(planeId);
                o2::mch::dcs::Cathode cath{deId, plane};
                auto dsIds = o2::mch::dcs::getDsIndices({cath});
                int first = *dsIds.begin();
                int last = *dsIds.rbegin();

                auto* hTrendDSperDE =
                    mapTrendsPerDE.at(Form("hTrendDSperDE%d%s", deId, planeId == 0 ? "B" : "NB"));

                for (int ds = first; ds <= last; ++ds) {
                    if (ds+1 > nBins) {
                        printf("ATTENZIONE: ds+1=%d > nBins=%d in %s\n", ds+1, nBins, fn.c_str());
                        continue;
                    }
                    double n = hNPerDS->GetBinContent(ds+1);
                    int dsLocBin = hTrendDSperDE->GetYaxis()->FindFixBin(ds);
                    hTrendDSperDE->SetBinContent(i+1, dsLocBin, n);
                }
            }
        }
    }

    // Chiudi e libera
    //printf("\n before closing the file \n");
    file->Close();
   //printf("\n after closing the file \n");
   delete hNPerDS;
   //printf("\n after deleting the DS histo \n"); 
   delete file;
   //printf("\n after deleting the file \n");
}*/



  // normalize distributions
printf("\n Normalizing per distribution \n");
  std::vector<double> nIntPerRun(nRuns, 0.);
  for (int deId : {100,101,102,103}) {
    auto* hDEB = mapTrendsPerDE.at(Form("hTrendDSperDE%dB", deId));
    auto* hDENB = mapTrendsPerDE.at(Form("hTrendDSperDE%dNB", deId));
    for (int ir = 1; ir <= nRuns; ++ir) {
      int run = rlist[ir-1];
      for (int ids = 1; ids <= hDEB->GetNbinsY(); ++ids) {
        double n_b = hDEB->GetBinContent(ir, ids);
        double n_nb = hDENB->GetBinContent(ir, ids);
        int ds = trunc(hDENB->GetYaxis()->GetBinLowEdge(ids));
        if (deId == 100 && (ds == 248 || ds == 408) && run == 559544)
          n_nb = 0;
        nIntPerRun[ir - 1] += n_b + n_nb;
      }
    }
  }

  // normalize by lumi
printf("\n Normalizing per lumi \n");
  for (auto& [hn, h] : mapTrendsPerDE) {
    for (int ir = 1; ir <= nRuns; ++ir) {
      for (int ids = 1; ids <= h->GetNbinsY(); ++ids) {
        double n = h->GetBinContent(ir, ids);
        double nNorm = n / nIntPerRun[ir - 1];
        h->SetBinContent(ir, ids, nNorm);
      }
    }
    //printf("\n before cloning histo_nn \n");
    auto* h_nn = (TH2D*) h->Clone(Form("%s_nn", h->GetName()));
    //h_nn->SetDirectory(0);
    //printf("\n after cloning the histo_nn \n");
    l_out->Add(h_nn);
  }
  printf("\n out of the lumi loop \n");
  for (auto& [hn, h] : mapTrendsPerDE)
    normalize(h);
    //normalize(h_nn);
printf("\n normalize(h) done \n");

  auto* ctrend = new TCanvas("ctrend", "ctrend", 5000, 3000);
  ctrend->cd();
printf("\n ctrend canvas created \n");

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
  printf("Starting final loop");
  for (auto& [chId, deIds] : deIdsPerChamber) {
    for (auto deId : deIds) {
      auto* hTrendDSperDE_B = mapTrendsPerDE.at(Form("hTrendDSperDE%dB", deId));
      auto* hTrendDSperDE_NB = mapTrendsPerDE.at(Form("hTrendDSperDE%dNB", deId));
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
      printf("We are at DE %d\n", deId);
      ctrend->cd();
      if (deId == 100)
        ctrend->Print(Form("./figures/%s/trending_ds_per_de.pdf(", type.c_str()), Form("Title: DE%d", deId));
      else if (deId == 1025)
        ctrend->Print(Form("./figures/%s/trending_ds_per_de.pdf)", type.c_str()), Form("Title: DE%d", deId));
      else
        ctrend->Print(Form("./figures/%s/trending_ds_per_de.pdf", type.c_str()), Form("Title: DE%d", deId));
    }
  }
  //std::vector<int> dsBending  = {2296, 2319};
  //std::vector<int> dsNonBend  = {1119, 1264};
  std::vector<int> dsBending  = {2296, 2319};
  std::vector<int> dsNonBend  = {758, 898};

  std::cout << "\n>>> Starting DS projections for all DEs...\n";

  for (auto& [chId, deIds] : deIdsPerChamber) {
    for (auto deId : deIds) {
      // Taking histograms bending
      auto* hTrendDSperDE_B = mapTrendsPerDE.at(Form("hTrendDSperDE%dB", deId));

      for (int ds : dsBending) {
        int binY_B  = hTrendDSperDE_B->GetYaxis()->FindBin(ds);

        auto* projB = hTrendDSperDE_B->ProjectionX(
            Form("proj_DE%d_DS%d_B", deId, ds), binY_B, binY_B);
        projB->SetTitle(Form("DE%d DS%d Bending;Run index;Counts", deId, ds));

        // saving in ROOT file
        fout->cd();
        projB->Write();
      }
    }
  }
for (auto& [chId, deIds] : deIdsPerChamber) {
    for (auto deId : deIds) {
      // Taking histograms bending 
      auto* hTrendDSperDE_NB = mapTrendsPerDE.at(Form("hTrendDSperDE%dNB", deId));

      for (int ds : dsNonBend) {
        int binY_NB  = hTrendDSperDE_NB->GetYaxis()->FindBin(ds);
 
        auto* projNB = hTrendDSperDE_NB->ProjectionX(
            Form("proj_DE%d_DS%d_NB", deId, ds), binY_NB, binY_NB);
        projNB->SetTitle(Form("DE%d DS%d Non Bending;Run index;Counts", deId, ds));

        // saving in ROOT file
        fout->cd();
        projNB->Write();
      }
    }
  }

  fout->cd();
  //l_out->Write();
  fout->Write();
  fout->Close();

  printf("\nAll projections saved successfully!\n");
  fout->cd();
TString pdfOut = "./figures/projections/all_projections.pdf";
auto* cProj = new TCanvas("cProj", "Projections", 5000, 3000);
bool firstPage = true;

for (auto& [chId, deIds] : deIdsPerChamber) {
  for (auto deId : deIds) {

    // ---- BENDING ----
    for (int ds : dsBending) {
      TString hname = Form("proj_DE%d_DS%d_B", deId, ds);
      auto* projB = (TH1D*) fout->Get(hname);
      if (!projB) {
        printf("Missing histogram %s\n", hname.Data());
        continue;
      }

      cProj->cd();
      projB->SetLineWidth(2);
      projB->Draw("hist");

      if (firstPage) {
        cProj->Print(pdfOut + "(", Form("Title: DE%d_DS%d_B", deId, ds)); // apre PDF multipagina
        firstPage = false;
      } else {
        cProj->Print(pdfOut, Form("Title: DE%d_DS%d_B", deId, ds));
      }
    }

    printf(" B projections for DE %d done.\n", deId);

    // ---- NON-BENDING ----
    for (int ds : dsNonBend) {
      TString hname = Form("proj_DE%d_DS%d_NB", deId, ds);
      auto* projNB = (TH1D*) fout->Get(hname);
      if (!projNB) {
        printf("Missing histogram %s\n", hname.Data());
        continue;
      }

      cProj->cd();
      projNB->SetLineWidth(2);
      projNB->Draw("hist");
      cProj->Print(pdfOut, Form("Title: DE%d_DS%d_NB", deId, ds));
    }

    printf("NB projections for DE %d done.\n", deId);
  }
}

// Chiudi il PDF multipagina
cProj->Print(pdfOut + ")", "Title: End");
printf("\n All projections successfully exported to %s\n", pdfOut.Data());
/*
auto* cProjB = new TCanvas("cProjB", "ProjectionsB", 5000, 3000);
auto* cProjNB = new TCanvas("cProjNB", "ProjectionsNB", 5000, 3000);
gSystem->Exec("mkdir -p ./figures/projections/");

for (auto& [chId, deIds] : deIdsPerChamber) {
  for (auto deId : deIds) {
    for (int ds : dsBending) {
      cProjB->cd();
      auto* projB = (TH1D*) fout->Get(Form("proj_DE%d_DS%d_B", deId, ds));
      if (!projB) continue;
      projB->Draw("hist");
      cProjB->Print(Form("./figures/projections/DE%d_DS%d_B.pdf", deId, ds));
    }
printf("B done successfully! \n");
    for (int ds : dsNonBend) {
      cProjNB->cd();
      auto* projNB = (TH1D*) fout->Get(Form("proj_DE%d_DS%d_NB", deId, ds));
      if (!projNB) continue;
      projNB->Draw("hist");
      cProjNB->Print(Form("./figures/projections/DE%d_DS%d_NB.pdf", deId, ds));
    }
printf("NB done successfully! \n");
  }
}*/

  //fout->cd();
  //l_out->Write();
  //fout->Write();
  fout->Close();
  //l_out->Write();
}
