#!/bash/bin

path=$PWD
for k in 1 4 7 10 13
do
    cd ${path}/${k}PROFESS_KERNEL/cd100
    cp Si.inpt cd100.inpt
    cp Si.ion cd100.ion
    sed -i "s/MAXI DEN 50/PRIN DEN/" cd100.inpt
    sed -i "s/ECUT\t900/GDEN 10.4/g" cd100.inpt
    cd ../ && pPROFESS cd100/cd100
done
    