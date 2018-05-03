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
