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
label=MonoHgg${suffix}_sig_ScalarZp
libs="-L libdiphotonsUtils"
rootversion=$(root-config --version| tr '.' ' ')
#[[ $rootversion -gt 5 ]] && libs="-L libdiphotonsRooUtils"

for mZ in "10" "15" "20" "50" "100" "200" "295" "300" "500" "1000" "1995" "2000" "10000" ; do # 
    echo $mZ
    for mDM in "1" "10" "50" "150" "500" "1000" ; do
	card=datacard_${folder}_${label}_mZP${mZ}_mChi${mDM}_13TeV.txt
	echo $card
	binary=$(echo $card | sed 's%.txt$%.root%')
	signame=$(echo $card | sed 's%.*MonoHgg_%%; s%.txt%%')
	echo $binary
	echo $signame
	log=combine_log_${method}_${label}_mH125_mZP${mZ}.log
	set -x
	combine $libs $args -n "${label}_mZP${mZ}_mChi${mDM}"  -m ${mDM} $card 2>&1 | tee $log
	set +x
	tail -5 $log
    done
hadd -f higgsCombine${label}Combined.$method.ScalarZp_mZP${mZ}.root higgsCombine${label}_mZP${mZ}_mChi*.Asymptotic.mH*.root
done


