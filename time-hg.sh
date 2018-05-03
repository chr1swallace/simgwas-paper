#!/usr/bin/bash

N=$1
NEFF=$2
NSIM=$3

FG="/rds/user/cew54/hpc-work/simgwas/input"
TMP1=$(mktemp /tmp/hg.XXXXXX)
TMP2=$(mktemp /tmp/hg.XXXXXX)
MAP="/home/cew54/share/Data/reference/1000GP_Phase3/genetic_map_chr21_combined_b37.txt"

case $NEFF in
    1)
	GSTR="9557710 1 1.5 2.25"
	;;	
    2)
	  GSTR="9557710 1 1.5 2.25 9429998 1 1.8 3.24"
	  ;;
    3)
	GSTR="9557710 1 1.5 2.25 9429998 1 1.8 3.24 9579507 1 1.2 1.44"
	;;
    4)
	GSTR="9557710 1 1.5 2.25 9429998 1 1.8 3.24 9579507 1 1.2 1.44 9828538 1 1.8 3.24"
	;;
    5)
	GSTR="9557710 1 1.5 2.25 9429998 1 1.8 3.24 9579507 1 1.2 1.44 9828538 1 1.8 3.24 9507863 1 1.2 1.44 9476405"
	;;
    6)
	GSTR="9557710 1 1.5 2.25 9429998 1 1.8 3.24 9579507 1 1.2 1.44 9828538 1 1.8 3.24 9507863 1 1.2 1.44 9476405 1 1.5 2.25"
	;;
esac

for i in `seq 1 $NSIM`; do
    /home/cew54/localc/bin/hapgen2 \
	-m $MAP -l ${FG}.leg -h ${FG}.hap -o $TMP1 \
	-dl $GSTR -n $N $N -no_haps_output
    /home/cew54/localc/bin/snptest -data ${TMP1}.controls.gen ${TMP1}.controls.sample \
				   ${TMP1}.cases.gen ${TMP1}.cases.sample \
				   -o $TMP2 -frequentist 1 -method score -pheno pheno
done
	 
rm ${TMP1}.*
rm ${TMP2}.*
