#!/bin/bash

path=`pwd`
fold='5_results'

for i in 1 3
do
    cd ${path}/pack
    k=$(($i+2))
    # mv log ${k}log
    # mv e_v.dat ${k}e_v.dat
    # mv coef.dat ${k}coef.dat
    # mv PROFESS_KERNEL.dat ${k}PROFESS_KERNEL.dat
    # mv optcoef.dat ${k}optcoef.dat
    # mv opte_v.dat ${k}opte_v.dat
    # mv opt_PROFESS_KERNEL.dat ${k}opt_PROFESS_KERNEL.dat
    # mv origincoef.dat ${k}origincoef.dat
    # mv origine_v.dat ${k}origine_v.dat
    # mv origin_PROFESS_KERNEL.dat ${k}origin_PROFESS_KERNEL.dat 
    cp ${i}optcoef.dat ${path}/${fold}/${k}optcoef.dat
    cp ${i}origincoef.dat ${path}/${fold}/${k}origincoef.dat
    cp ${i}opte_v.dat ${path}/${fold}/${k}opte_v.dat
    cp ${i}origine_v.dat ${path}/${fold}/${k}origine_v.dat
    cp ${i}opt_PROFESS_KERNEL.dat ${path}/${fold}/${k}PROFESS_KERNEL.dat
    cp ${i}origin_PROFESS_KERNEL.dat ${path}/${fold}/${k}origin_PROFESS_KERNEL.dat
    cp ${i}log ${path}/${fold}/${k}log
done
# cd ${path}
# tar -zcvf pack.tar.gz pack/*
