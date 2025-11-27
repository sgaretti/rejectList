#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TColor.h"
#include "TROOT.h"
#include "TList.h"
#include "TStyle.h"
#include "TSystem.h"

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#endif

#include "fmt/format.h"

std::vector<std::string> parseString(const TString& rs)
{
  std::string s = rs.Data();
  std::stringstream ss(s);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings(begin, end);
  // std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
  return vstrings;
}

void upload_reject_lists()
{
  std::map<int, std::vector<int>> runRejectLists{};

  std::fstream fs("banned_ds.list", std::fstream::in);

  int run, ds;
  while (fs >> run >> ds) {
    runRejectLists[run].push_back(ds);
  }

  o2::ccdb::CcdbApi ccdb;
  std::string ccdbPath = "http://alice-ccdb.cern.ch";
  std::string ccdbPathUpload = "http://alice-ccdb.cern.ch/Users/n/nburmaso/2023_pbpb_apass5_v1";
  ccdb.init(ccdbPath.c_str());

  auto upload = [&](int runNumber)
  {
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb, runNumber);
    auto eor = soreor.second;
    auto sor = soreor.first;
    const auto& rejectList = runRejectLists.at(runNumber);
    std::string dsList = "";
    for (auto ds : rejectList) {
      dsList += std::to_string(ds);
      dsList += " ";
    }
    auto cmd = fmt::format("o2-mch-bad-channels-ccdb -u --ccdb {} --put --starttimestamp {} --endtimestamp {} --type RejectList --ds {}",
                           ccdbPathUpload, sor - 1000, eor + 10000, dsList);
    // fmt::print("{}\n", cmd.c_str());
    fmt::print("uploading: {} {} {}\n", run, eor, sor); fflush(stdout);
    system(fmt::format("{}", cmd).c_str());
  };

  for (const auto& item : runRejectLists) {
    int run = item.first;
    upload(run);
  }
}
