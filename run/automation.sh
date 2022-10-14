#!/bin/bash
j=0
while read line 
do
    let j++
done < listofvalues.txt
#echo "line numbers "$j
for ind in {1..10}
do
    line=$(sed -n ${ind}'p' ./listofvalues.txt) 
    echo $line
    position=$(expr ${#line} - 3) # position of t_0
    t=${line:$position} # t_0
    echo t is $t
    timesteps=$(expr $t \* 4878) 
    echo timestep is $timesteps 
    sed -i "41s/.*/$line/" ./biozement/input.dat
    sed -i "13s/.*/max  $timesteps/" ./biozement/input.dat
    mpirun -n 1 program_files/bin/biozement/lb_ejah 
    echo $line > ./biozement/out/inputline.txt
    zip -r out-$ind.zip ./biozement/out
    rm -r biozement/out
done
