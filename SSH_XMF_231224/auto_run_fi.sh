#!/bin/bash

#----------------------------------------------------------------------#

list_LB='0.05 0.1 1.0'
list_IR='0.99 0.5 0.2'
list_GM='1.0001 1.1 1.2' 
#list_LIR='0.1 0.15 0.2'

#----------------------------------------------------------------------#

cd c_test
make

for mp in ${list_MP};do
for rp in ${list_RP};do
  sed -i "s|R_planet .*|R_planet ${mp}|g" test_rn.ini
  sed -i "s|M_planet .*|M_planet ${rp}|g" test_rn.ini   
  ./main   
done
done
  
echo -e "\n  Finished \n"
date "+  DATE: %Y-%m-%d%n  TIME: %H:%M:%S%n"








































