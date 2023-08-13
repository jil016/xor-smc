#!/bin/bash

set -x
# for i in 0 10 20
# do
# for j in 0 10 11 16 18 32 60 102 113 114 119 121 147 166 172 200 203 218 219 222 234 239 240 242
# do
#    ./sharpSAT -decot 50 -decow 100 -tmpdir . -cs 15000 ../shelter_CNFs/shelter_CNFs_250/src${i}_sink${j}.cnf > ../shelter_CNFs/shelter_CNFs_250/src${i}_sink${j}.out
#    done
# done

for i in 0 10 20
do
for j in 1 10 13 32 39 42 50 84 85 87 113 122 165 173 175 176 200 218 248 311 314 371
do
   ./sharpSAT -decot 50 -decow 100 -tmpdir . -cs 15000 ../shelter_CNFs/shelter_CNFs_388/src${i}_sink${j}.cnf > ../shelter_CNFs/shelter_CNFs_388/src${i}_sink${j}.out
   done
done