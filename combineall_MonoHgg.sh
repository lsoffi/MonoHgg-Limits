#!/bin/bash

# set -x

folder_old=$1 && shift
coupl=$1 && shift

# echo $coupl

args=""
while [[ -n $1 ]]; do
    case $1 in
	-M)
	    method=$2
	    args="$args $1 $2"
	    shift
	    ;;
	-n)
	    label=$2
	    shift
	    ;;
	--suffix)
	    suffix=$2
	    shift
	    ;;
	--hadd)
	    hadd="hadd"
	    ;;
	--dry-run)
	    dry="1"
	    ;;
	--cont)
	    cont="1"
	    ;;
	*)	    
	    args="$args $1"
	    ;;	    
    esac
    shift
done
shift
folder=${folder_old}${suffix}
cd $folder
echo $folder

libs="-L libdiphotonsUtils"
rootversion=$(root-config --version| tr '.' ' ')
#[[ $rootversion -gt 5 ]] && libs="-L libdiphotonsRooUtils"

for mA0 in "300"  "400" "500" "600" "700" "800"
  do
  for mZ in "600" "800" "1000" "1200" "1400" "1700" "2000" "2500"
    do # 
    echo $mZ,$mA0,$suffix
    label=MonoHgg${suffix}
    card=datacard_${folder}_${label}_sig_2HDM_mZP${mZ}_mA0${mA0}_13TeV.txt
#    card=DataCard_2HDM_mZP${mZ}_mA0${mA0}.txt
    echo $card
    binary=$(echo $card | sed 's%.txt$%.root%')
    signame=$(echo $card | sed 's%.*MonoHgg_%%; s%.txt%%')
    echo $binary
    echo $signame

    log=combine_log_${method}_${label}_mH125_mZP${mZ}.log
    set -x
    combine $libs $args -n "${label}_mZP${mZ}_mA0${mA0}" -m ${mZ} $card 2>&1 | tee $log
    set +x
    tail -5 $log
done
hadd -f higgsCombine${label}_mA0${mA0}.Combined.$method.2HDM.root higgsCombine${label}_mZP*_mA0${mA0}.${method}.mH*.root
done


