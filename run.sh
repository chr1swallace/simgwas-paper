
for i in `seq 7 8`; do
    # qR.rb -a TODD-SL3-CPU  -y 0-9 -c 1 -r -t "10:00:00" runsims-1kg.R --args NCV=$i NSIM=500 NCSE=3000
    qR.rb -a TODD-SL3-CPU   -y 0-19 -c 1 -r -t "1:00:00" runsims-1kg.R --args SPECIAL=$i NSIM=50 NCSE=1000 NCTL=1000
    qR.rb -a TODD-SL3-CPU   -y 0-99 -c 1 -r -t "1:00:00" runsims-1kg.R --args SPECIAL=$i NSIM=10 NCSE=5000 NCTL=5000
    # qR.rb -a TODD-SL3-CPU -y 0-9 -c 1 -r -t "10:00:00" runsims-1kg.R --args NCV=$i NSIM=500 NCSE=1000
done
