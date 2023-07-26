



#wget https://github.com/meelgroup/approxmc/releases/download/3.0/approxmc-linux-x64.gz


start=`date +%s%N`
./approxmc $1.cnf
end=`date +%s%N`
echo Execution time was `expr $end - $start` nanoseconds.