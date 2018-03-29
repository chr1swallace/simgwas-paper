#!/usr/local/bin/bash

## make genotype data
hapgen2 -m ./example/ex.map -l ./example/ex.leg -h ./example/ex.haps -o ./example/ex.out -dl 2190692 1 1.5 2.25 2190692 1 2 4 -n 1000 1000 -t ./example/ex.tags \
 -no_haps_output 

## create input for snptest
cp example/ex.out.controls.sample out.sample
sed '1,2d' example/ex.out.cases.sample >> out.sample
cat example/ex.out.controls.gen example/ex.out.cases.gen > out.gen

## run snptest
# snptest -data out.gen out.sample \
    snptest -data example/ex.out.controls.gen example/ex.out.controls.sample example/ex.out.cases.gen example/ex.out.cases.sample \
-o snptest.out \
-frequentist 1 \
-method score \
-pheno pheno

# filename1=""
# filename10="/nfs/trans_ethnic_sim/simulated_data/scenario3/scenario3c/setting3/chr3_IGF2BP2/reps_to_redo_again.txt"
# #for((filenum=1;$filenum<$totalfiles;filenum+=1))
# for filenum in `cat ${filename10}`
# do
# filename1="/nfs/trans_ethnic_sim/simulated_data/scenario3/scenario3c/setting3/chr3_IGF2BP2/CEU/CEU_causal_s3_rep_"$filenum
# filename2="/nfs/trans_ethnic_sim/simulated_data/scenario3/scenario3c/setting3/chr3_IGF2BP2/CEU/CEU_causal_s3_causals_rep_"$filenum".txt"
# cut -f 1 -d " " $filename2 > causal.txt
# read cpos < causal.txt
# /nfs/trans_ethnic_sim/hapgen2 -h /nfs/trans_ethnic_sim/1KG_RPs/CEU_1000G_phase1interim_jun2011_chr3_impute.hap \
#  -l /nfs/trans_ethnic_sim/1KG_RPs/ALL_1000G_phase1interim_jun2011_chr3_impute.legend \
#  -m /nfs/trans_ethnic_sim/1KG_RPs/genetic_map_chr3_combined_b37.txt \
#  -o $filename1".gz" \
#  -no_haps_output \
#  -dl $cpos 1 1.2 1.44 \
#  -int 184811528 186092827 \
#  -n 1000 1000 \
#  -Ne 20000 \
#  -t /nfs/trans_ethnic_sim/1KG_RPs/Human660W-Quad_v1_C-b37-chr3.strand
# done
