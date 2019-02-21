# simgwas-paper
calculations run for our simGWAS manuscript 

> simGWAS: a fast method for simulation of large scale case-control GWAS summary statistics
> Mary D Fortune, Chris Wallace
> Bioinformatics 2018; https://doi.org/10.1093/bioinformatics/bty898
> preprint at bioRxiv 313023; doi: https://doi.org/10.1101/313023

Note that commands like `q.rb` and `qR.rb` are wrappers to run these commands on out local HPC.  This code is provided so that the details of what was done can be understood, but it is expected that commands to submit jobs to local HPCs would need to be customised.


## Timings for comparison to HAPGEN2 + SNPTEST2

### how does number of causal variants change things?
```{sh}
for i in `seq 1 6`; do
    q.rb -r -y 0-49 -j hg_cvs-$i ./time-hg.sh 2000 $i 1
    q.rb -r -y 0-49 -j meth_cvs-$i ./time-meth1.R --args N=2000 NCV=$i NSIM=1
done

jobs=""
for i in `seq 1 6`; do
    jobs="${jobs},hg_cvs-$i,meth_cvs-$i"
done

sacct --name $jobs --format=JobID,JobName%50,State,TotalCPU,UserCPU | grep _cvs |grep COMPLETED > timings/cvs.txt
```

### how does number of cases/controls change things?
```{sh}
for i in 1000 2000 4000 8000 16000 32000 64000; do
    q.rb -r -y 0-9 -j hg_smp-$i ./time-hg.sh $i 3 1
    q.rb -r -y 0-9 -j meth_smp-$i ./time-meth1.R --args N=$i NCV=3 NSIM=1
done


jobs=""
for i in 1000 2000 4000 8000 16000 32000 64000; do
    jobs="${jobs},hg_smp-$i,meth_smp-$i"
done

sacct --name $jobs --format=JobID,JobName%50,State,TotalCPU,UserCPU | grep _smp > timings/smp.txt

# i=16000
# jobs="${jobs},hg-smp-$i,meth-smp-$i"
# sacct --name $jobs --format=JobID,JobName%50,State,TotalCPU,UserCPU | grep smp >> timings/smp.txt
```

### how does number of repetitions change things?
```{sh}
for i in 1 2 4 8 16 32 64 128; do
    q.rb -r -y 0-9 -j hg_sim-$i ./time-hg.sh 2000 3 $i
    q.rb -r -y 0-9 -j meth_sim-$i ./time-meth1.R --args N=2000 NCV=3 NSIM=$i
done

jobs=""
for i in 1 2 4 8 16 32 64 128; do
    jobs="${jobs},hg_sim-$i,meth_sim-$i"
done
sacct --name $jobs --format=JobID,JobName%50,State,TotalCPU,UserCPU | grep _sim > timings/sim.txt

```

### plot results

```{sh}
./plot-timings.R

```

## Simulations to check accuracy

```{sh}
for i in `seq 1 5`; do
    qR.rb -y 0-19 -c 1 -r -t "1:00:00" runsims-1kg.R --args SPECIAL=$i NSIM=50 NCSE=1000 NCTL=1000
    qR.rb -y 0-99 -c 1 -r -t "1:00:00" runsims-1kg.R --args SPECIAL=$i NSIM=10 NCSE=5000 NCTL=5000
done
```

### forward sims to check at CVs
for i in `seq 1 5`; do
    qR.rb -c 1 -r -t "1:00:00" forwardsims-1kg.R --args SPECIAL=$i NSIM=1000 NCSE=1000 NCTL=1000
    qR.rb -c 1 -r -t "1:00:00" forwardsims-1kg.R --args SPECIAL=$i NSIM=1000 NCSE=5000 NCTL=5000
done

collate the results and plot

```{sh}
./summarise-1kg.R
```


# Parts of what is needed to run on UK10K samples

## prepare bcf files for UK10K samples
The data were downloaded in the format of hap-legend-sample.  Indexed bcf files were created using the script below:

```{sh}
for i in 5 `seq 10 22`; do
    echo $i
    zcat _EGAZ00001017893_${i}.UK10K_COHORT.REL-2012-06-02.beagle.anno.csq.shapeit.20140306.legend.gz | sed "s/chr${i}[^ ]*  /chr${i}:/" | sed 's/ /_/g' | gzip -c > chr${i}.legend.gz
    cp _EGAZ00001017893_${i}.UK10K_COHORT.REL-2012-06-02.beagle.anno.csq.shapeit.20140306.hap.gz chr${i}.hap.gz
    cp _EGAZ00001017893_${i}.UK10K_COHORT.REL-2012-06-02.beagle.anno.csq.shapeit.20140306.sample chr${i}.samples
    q.rb -r -c1 -t "02:00:00" "bcftools convert -o chr${i}.bcf.gz -Ob  --haplegendsample2vcf chr${i} && bcftools index chr${i}.bcf.gz"
done

rm chr${i}.hap.gz chr${i}.samples chr${i}.legend.gz

```

## precompute LD matrices for UK10K data
```{sh}
qR.rb -j ld -t "4:00:00" -r -y 1-22 -n chr ./uk10k-ld.R
```

## precompute standard error matrices for UK10K data

These two steps should have been combined into one.  But weren't, and no point to rewrite now.
```{sh}
qR.rb -j ld -t "4:00:00" -r -y 1-22 -n chr ./uk10k-ld.R
qR.rb -j ld -t "4:00:00" -r -y 1-22 -n chr ./uk10k-se.R
```

## run a whole chromosome
NB, not yet optimised to use the pre-calculated ld matrices yet
```{sh}
./run-chr.R
```
