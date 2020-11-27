#!/bin/bash

# Bash script for correcting the itp files for circularized MSP

itp=$1

grep ' BB ' ${itp} > temp_BB

head -n 3 temp_BB > temp_BB_sel
tail -n 3 temp_BB >> temp_BB_sel

for i in `seq 6`; do sed -n ${i}p temp_BB_sel | gawk '{print $1}' > ResBB_${i}; done


#Constraint
c="      $(cat ResBB_1)     $(cat ResBB_6)      1   0.31000"

#Angles
a="    $(cat ResBB_5)     $(cat ResBB_6)     $(cat ResBB_1)      2     96   700"
a1="    $(cat ResBB_6)     $(cat ResBB_1)     $(cat ResBB_2)      2     96   700"

#Dihedrals
d=" $(cat ResBB_5) $(cat ResBB_6)     $(cat ResBB_1)     $(cat ResBB_2)      1   -120   400     1"
d1=" $(cat ResBB_4) $(cat ResBB_5)     $(cat ResBB_6)     $(cat ResBB_1)      1   -120   400     1"
d2=" $(cat ResBB_6) $(cat ResBB_1)     $(cat ResBB_2)     $(cat ResBB_3)      1   -120   400     1"


sed -i "/\[ con/a$c" ${itp} 

sed -i "/\[ ang/a$a" ${itp}
sed -i "/\[ ang/a$a1" ${itp}

sed -i "/\[ dih/a$d" ${itp}
sed -i "/\[ dih/a$d1" ${itp}
sed -i "/\[ dih/a$d2" ${itp}

rm temp_BB temp_BB_sel ResBB_[1-6]

