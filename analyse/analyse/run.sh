#!/bin/bash

year="2024"
pass="apass1_muon"

# list of histograms
find $PWD/ -name qc.root | grep ${year} | grep ${pass} | grep -v 559803 | grep -v 559827 | grep -v 560164 | grep -v 560168 | grep -v 560169 | grep -v 560184 | grep -v 560223 | grep -v 560310 | grep -v 560313 | grep -v 560330 | grep -v 560331 | sort -V > input_hists.txt
echo "list of histos taken"

# extract list of runs
rm runs.list &> /dev/null
for line in `cat input_hists.txt`; do
  f="${pass}/qc.root"
  rest=${line#*${f}}
  echo "LINE: $line"
  pos=$(( ${#line} - ${#rest} - ${#f} ))
  echo "POS: $pos"
  start=$((pos-7))
  echo "START: $start"
  r=${line:start:6}
  echo $r >> runs.list
done

echo "run list done"
# run for digits
if [[ ! -e figures/digits ]]; then mkdir -p ./figures/digits/; fi
echo "running command digits"
root -q -e "gSystem->Load(\"${O2_ROOT}/lib/libO2MCHMappingImpl4.so\"); gSystem->Load(\"${O2_ROOT}/lib/libO2MCHConditions.so\")" -x 'analyze_trending.cpp++(3)' #&> digits.log
echo "digits done"
# run for clusters
#if [[ ! -e figures/clusters ]]; then mkdir -p ./figures/clusters/; fi
#root -b -q -e "gSystem->Load(\"${O2_ROOT}/lib/libO2MCHMappingImpl4.so\");gSystem->Load(\"${O2_ROOT}/lib/libO2MCHConditions.so\")" -x 'analyze_trending.cpp++(1)' &> clusters.log
#echo "clusters done"
