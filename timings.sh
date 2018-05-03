## how does number of cvs change things?
for i in `seq 1 6`; do
    q.rb -r -y 0-49 -j hg_cvs-$i ./time-hg.sh 2000 $i 1
    q.rb -r -y 0-49 -j meth_cvs-$i ./time-meth1.R --args N=2000 NCV=$i NSIM=1
done

jobs=""
for i in `seq 1 6`; do
    jobs="${jobs},hg_cvs-$i,meth_cvs-$i"
done

sacct --name $jobs --format=JobID,JobName%50,State,TotalCPU,UserCPU | grep _cvs |grep COMPLETED > timings/cvs.txt


## how does number of cases/controls change things?
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



## how does number of repetitions change things?
for i in 1 2 4 8 16 32 64 128; do
    q.rb -r -y 0-9 -j hg_sim-$i ./time-hg.sh 2000 3 $i
    q.rb -r -y 0-9 -j meth_sim-$i ./time-meth1.R --args N=2000 NCV=3 NSIM=$i
done

jobs=""
for i in 1 2 4 8 16 32 64 128; do
    jobs="${jobs},hg_sim-$i,meth_sim-$i"
done
sacct --name $jobs --format=JobID,JobName%50,State,TotalCPU,UserCPU | grep _sim > timings/sim.txt




./plot-timings.R
