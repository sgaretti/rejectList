#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TList.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGrid.h"

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "QualityControl/MonitorObject.h"
#include "QualityControl/MonitorObjectCollection.h"
#include "Common/TH1Ratio.h"
#endif


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


void extract_qc_hists(TString inputs="input_data.txt")
{
  printf("Input list: %s\n", inputs.Data());

  std::vector<std::string> files_list{};
  getList(std::string(inputs.Data()), files_list);

  TGrid::Connect("alien://");

  std::vector<TH1D*> vHistsDigitsRat{};
  std::vector<TH1F*> vHistsDigitsOrb{};
  std::vector<TH1F*> vHistsClusters{};

  int iFile = 0;
  auto nFiles = (int) files_list.size();
  for (auto &fname : files_list) {
    iFile++;
    printf("\n%d/%d -- File:%s\n", iFile, nFiles, fname.c_str()); fflush(stdout);
    auto* file = TFile::Open(fname.c_str());
    auto* dir = dynamic_cast<TDirectoryFile*>(file->GetDirectory("int"));
    auto* mhcDir = dynamic_cast<TDirectoryFile*>(dir->GetDirectory("MCH"));

    auto* clusters_coll = dynamic_cast<o2::quality_control::core::MonitorObjectCollection*>(mhcDir->FindObjectAny("Clusters"));
    if (clusters_coll != nullptr) {
      auto* clusters_obj = dynamic_cast<o2::quality_control::core::MonitorObject*>(clusters_coll->FindObject("ClustersPerDualSampa"));
      if (clusters_obj != nullptr) {
        auto* hClusters = dynamic_cast<TH1F*>(clusters_obj->getObject());
        hClusters->SetName(Form("hClusters_%d", iFile));
        hClusters->SetDirectory(nullptr);
        vHistsClusters.push_back(hClusters);
      }
    }

    auto* digits_coll = dynamic_cast<o2::quality_control::core::MonitorObjectCollection*>(mhcDir->FindObjectAny("Digits"));
    if (digits_coll != nullptr) {
      auto* digits_obj_rate = dynamic_cast<o2::quality_control::core::MonitorObject*>(digits_coll->FindObject("RateSignalPerDualSampa"));
      if (digits_obj_rate != nullptr) {
        auto* hDigitsQcRatio = dynamic_cast<o2::quality_control_modules::common::TH1DRatio*>(digits_obj_rate->getObject());
        auto* hDigitsRat = hDigitsQcRatio->getNum();
        hDigitsRat->SetName(Form("hDigitsRat_%d", iFile));
        hDigitsRat->SetTitle(Form("Number of digits per dual sampa"));
        hDigitsRat->SetDirectory(nullptr);
        vHistsDigitsRat.push_back(hDigitsRat);
      }

      auto* digits_obj_orbits = dynamic_cast<o2::quality_control::core::MonitorObject*>(digits_coll->FindObject("DigitSignalOrbit_Elec"));
      if (digits_obj_orbits != nullptr) {
        auto* hDigitOrbits = dynamic_cast<TH2F*>(digits_obj_orbits->getObject());
        auto* hDigitsOrb = (TH1F*)hDigitOrbits->ProjectionX();
        hDigitsOrb->SetName(Form("hDigitsOrb_%d", iFile));
        hDigitsOrb->SetTitle(Form("Number of digits per dual sampa"));
        hDigitsOrb->SetDirectory(nullptr);
        vHistsDigitsOrb.push_back(hDigitsOrb);
      }
    }

    file->Close();
  }

  auto* fout = new TFile("qc.root", "recreate");

  if (!vHistsDigitsRat.empty()) {
    auto* hDigitsRat = (TH1F*)vHistsDigitsRat[0]->Clone("hDigitsRat");
    if (vHistsDigitsRat.size() > 1) {
      for (int i = 1; i < vHistsDigitsRat.size(); ++i) { 
        hDigitsRat->Add(vHistsDigitsRat[i]);
      }
    }
    hDigitsRat->Write();
  }

  if (!vHistsDigitsOrb.empty()) {
    auto* hDigitsOrb = (TH1F*)vHistsDigitsOrb[0]->Clone("hDigitsOrb");
    if (vHistsDigitsOrb.size() > 1) {
      for (int i = 1; i < vHistsDigitsOrb.size(); ++i) { 
        hDigitsOrb->Add(vHistsDigitsOrb[i]);
      }
    }
    hDigitsOrb->Write();
  }

  if (!vHistsClusters.empty()) {
    auto* hClusters = (TH1F*)vHistsClusters[0]->Clone("hClusters");
    if (vHistsClusters.size() > 1) {
      for (int i = 1; i < vHistsClusters.size(); ++i) { 
        hClusters->Add(vHistsClusters[i]);
      }
    }
    hClusters->Write();
  }

  fout->Close();
}