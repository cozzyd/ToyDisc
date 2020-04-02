#! /bin/sh

root -b -q -e '.L ToyDisc.C+' 


N=1e9 
jobsfile=JOBS
rm -f $jobsfile

for E in 1e18 3.16e17 1e17 3.16e16 1e16
do
  for forced in 0 1
  do
  echo "root -b -q ToyDisc.C+\($forced,$E,$N\) > ${forced}_$E.log" >> $jobsfile 
done 
done

parallel < $jobsfile 
