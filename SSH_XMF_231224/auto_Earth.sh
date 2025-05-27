#!/bin/bash

#----------------------------------------------------------------------#
#run switchs
END_TIME=300.1
OUT_ITV=5

#model switchs
list_RHO='0.05'

list_TMP='600 800 1000 1200 1500'

list_IR='0.5' #0.5 0.2

list_GM='1.0001'

list_LB='2' #5 10

#----------------------------------------------------------------------#
INI_FILE='p_Earth.ini'
maindir='data_Earth_240120_r10'

#run settings

sed -i "s|tstop .*|tstop ${END_TIME}|g" ${INI_FILE}
sed -i "s|OUTPUT_INTERVAL .*|OUTPUT_INTERVAL ${OUT_ITV}|g" ${INI_FILE}

#model settings
mkdir ${maindir}

for lr in ${list_RHO};do
for lt in ${list_TMP};do
for ir in ${list_IR};do
for gm in ${list_GM};do
for lb in ${list_LB};do
  #Set initial file
  dirname="T_B${lb}_RHO${lr}_T${lt}_IR${ir}_GM${gm}_8w"

  rho=$(echo "$lp*$list_RHOFAC" | bc -l)

  sed -i "s|output_dir .*|output_dir ${maindir}/${dirname}|g" ${INI_FILE}
  
  sed -i "s|LB_RHO .*|LB_RHO ${lr}|g" ${INI_FILE}
  sed -i "s|LB_TMP .*|LB_TMP ${lt}|g" ${INI_FILE}
  sed -i "s|H_ION_RATE .*|H_ION_RATE ${ir}|g" ${INI_FILE}
  sed -i "s|GAMMA .*|GAMMA ${gm}|g" ${INI_FILE}
  sed -i "s|LB_B .*|LB_B ${lb}|g" ${INI_FILE}
  

  #Make output directory
  mkdir ${maindir}/${dirname}

  #run
  mpirun -machinefile $PBS_O_WORKDIR/nodefile -np 128 ./pluto -i ${INI_FILE} -dec 8 16
   
done
done
done
done
done



echo -e "\n  Finished \n"
date "+  DATE: %Y-%m-%d%n  TIME: %H:%M:%S%n"








































