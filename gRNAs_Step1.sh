#!/bin/bash

helpFunction()
{
   echo "Usage: $0 -i [] -o [] -t [all/mut/snp/nd]"
   echo -e "\t-i Input xlsx file on s3"
   echo -e "\t-o output s3 directory"
   echo -e "\t-t Description of what gRNAs to produce"
   echo ""
   exit 1 # Exit script after printing help
}

runner()
{
  if [ $ex -eq 0 ]
  then
    echo ""
  elif [ $ex -eq 2 ]
  then
    echo "Failed generating gRNAs for $1" >&2
    exit 1
  else
    echo "Unexpected error occured, $1" >&2
    exit 1
  fi
}

while getopts "i:o:t:" opt
do
   case "$opt" in
      i ) input_path="$OPTARG" ;;
      o ) output_path="$OPTARG" ;;
      t ) type="${OPTARG}" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$input_path" ] || [ -z "$output_path" ] || [ -z "$type" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

if [ -z $CODE_LOCATION ] ; then
   CODE_LOCATION=$(dirname $(readlink -f $0))
fi

case $type in

  all)
    runner Mut
    runner SNP
    runner ND
    ;;

  mut)
    runner Mut
    ;;

  snp)
    runner SNP
    ;;

  nd)
    runner ND
    ;;
  *)
    echo -n "unknown"
    ;;
esac

d=$(date +%d-%m-%Y)
filename="${input_path##*/}"
filename=${filename%".xlsx"}
outfilename=$output_path/${filename}_OUT_$d.xlsx

