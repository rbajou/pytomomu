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

list_file=$(ls -d ${indir}/*mon.root)

for f in ${list_file}
do 
        echo $f
        ts=$(stat -c "%Y" ${f})
        date=$(date -d "@${ts}" "+%Y-%m-%d %H:%M:%S")
        if [[ $f =~ ([0-9]{8})_([0-9]{2}H[0-9]{2}) ]]; then
            extracted_date="${BASH_REMATCH[1]}"
            extracted_time="${BASH_REMATCH[2]}"
            # echo "Extracted date: $extracted_date"
            # echo "Extracted time: $extracted_time"
            # Replace 'H' with ':' in the time component
            formatted_time="${extracted_time//H/:}"
            # Combine date and time for the 'date' command
            datetime="$extracted_date $formatted_time"
            # Convert the combined date and time to a Unix timestamp
            unix_timestamp=$(date -d "$datetime" +"%s")

            if [ $? -eq 0 ]; then
                echo "Unix timestamp: $unix_timestamp"
            else
                echo "Error converting date and time to Unix timestamp."
            fi
        else
            echo "Date and/or time not found in the filename."
        fi
done


# unix_timestamp=$(stat -c "%Y" ${indir}/${prefix}_run${irun}_analyse.root)
# echo $unix_timestamp >  ${outdir}/timestamp0.txt
# # date0=$(date -d "@$timestamp" "+%Y-%m-%d %H:%M:%S")
# date0=$(date -d "@${unix_timestamp}" "+%Y-%m-%d %H:%M:%S")
# echo "Start run: ${date0}"
# rootfile=${outdir}/${prefix}_run${irun}_analyse1.root
# if [ ! -f $rootfile ]
# then  
#     echo "${rootfile} don't exist"
#     root -q "${scriptdir}/shorten_analysis1_${det}.C($irun)" 
# fi
# echo $outdir
# ls -lh $outdir