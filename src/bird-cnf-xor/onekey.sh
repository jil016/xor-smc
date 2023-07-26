#!/usr/bin/zsh
for i in {1..100}
do
    echo "Number: $i"
    ./bird-cnf-xor.sh $i > $i.log
done
