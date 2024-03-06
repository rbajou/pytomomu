#!/bin/bash 

det=$1 #e.g izen
prefix=$2 #e.g IZEN_COLIS_MID
irun=$3 #e.g 1
indir=${HOME}/Projects/tomomu/${det}/mnt
outdir=${HOME}/Projects/tomomu/Data/${det}/run${irun} 
scriptdir=${HOME}/Projects/tomomu/pytomomu/macro
if [ ! -d $outdir ]
then 
    mkdir $outdir 
fi 
unix_timestamp=$(stat -c "%Y" ${indir}/${prefix}_run${irun}_analyse.root)
echo $unix_timestamp >  ${outdir}/timestamp0.txt
# date0=$(date -d "@$timestamp" "+%Y-%m-%d %H:%M:%S")
date0=$(date -d "@${unix_timestamp}" "+%Y-%m-%d %H:%M:%S")
echo "Start run: ${date0}"
rootfile=${outdir}/${prefix}_run${irun}_analyse1.root
if [ ! -f $rootfile ]
then  
    echo "${rootfile} don't exist"
    root -q "${scriptdir}/shorten_analysis1_${det}.C($irun)" 
fi
echo $outdir
ls -lh $outdir