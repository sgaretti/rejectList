#!/bin/bash
# script for collecting input QC histograms per run per period

wd=$PWD

datadir="/alice/data/2024"

periods=()
for p in `cat periods.list`; do
  periods+=(${p})
done

pass="apass1_muon"

for period in ${periods[@]}; do
  echo $period
  runs_=$(alien.py ls ${datadir}/${period}/) # list of runs
  runs=()
  for r in ${runs_[@]}; do
    r="${r%/}" # remove trailing '/'
    runs+=(${r})
  done
  echo ${runs[@]}
  for r in ${runs[@]}; do
    alien.py find2 ${datadir}/${period}/${r}/${pass}/*/QC/QC.root > list.txt
    # alien.py find2 ${datadir}/${period}/${r}/${pass}/*/QC/001/QC.root > list.txt
    if [ ! -s list.txt ]; then
      # no qc / run is not in reco pass ...
      continue
    fi
    mkdir -p ${wd}/qc_hists/${datadir}/${period}/${r}/${pass}/
    rm ${wd}/qc_hists/${datadir}/${period}/${r}/${pass}/input_data.txt &> /dev/null
    for l in `cat list.txt`; do
      echo "alien://${l}" >> ${wd}/qc_hists/${datadir}/${period}/${r}/${pass}/input_data.txt
    done
    list="${wd}/qc_hists/${datadir}/${period}/${r}/${pass}/input_data.txt"
    root -l -b -q "extract_qc_hists.cpp+(\"${list}\")"
    mv qc.root ${wd}/qc_hists/${datadir}/${period}/${r}/${pass}/
  done
done
