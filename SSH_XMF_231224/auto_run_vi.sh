#!/bin/bash

#----------------------------------------------------------------------#

list_LP='0.05'
#0.2 1.0
list_IR='0.99 0.5 0.2'
#
list_GM='1.0001'
#
list_LB='5'
#0.1 1 2 
#list_LIR='0.1 0.15 0.2'
list_PM=''
list_RP=''
list_OR=''

#----------------------------------------------------------------------#


for lp in ${list_LP};do
for ir in ${list_IR};do
for gm in ${list_GM};do
for lb in ${list_LB};do
  #Set initial file
  dirname="T_B${lb}_P${lp}_IR${ir}_GM${gm}_8w"
  maindir='data_ISO/data_vi_test_231025'

  sed -i "s|output_dir .*|output_dir ${maindir}/${dirname}|g" p_ISO.ini

  sed -i "s|LB_PRS .*|LB_PRS ${lp}|g" p_ISO.ini
  sed -i "s|LB_RHO .*|LB_RHO ${lp}|g" p_ISO.ini
  sed -i "s|H_ION_RATE .*|H_ION_RATE ${ir}|g" p_ISO.ini
  sed -i "s|GAMMA .*|GAMMA ${gm}|g" p_ISO.ini
  sed -i "s|LB_B .*|LB_B ${lb}|g" p_ISO.ini

  #Make output directory
  mkdir ${maindir}/${dirname}

  #run pluto
  mpirun -np 64 ./pluto -dec 8 8 -i p_ISO.ini
   
done
done
done
done
  
echo -e "\n  Finished \n"
date "+  DATE: %Y-%m-%d%n  TIME: %H:%M:%S%n"








































