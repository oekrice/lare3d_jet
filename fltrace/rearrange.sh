#!/bin/bash

rm -r ./plots/*.png
scp -r pjb20205@login.hpc.strath.ac.uk:/users/pjb20205/lare3d_jet/fltrace/plots ./


counter=0
divcounter=0
while [ $counter -le 500 ]
    do
    countzeros=$(printf "%04d" $counter)
    divcountzeros=$(printf "%04d" $divcounter)
    fname=./plots/b$countzeros.png
    fnamenew=./plots/b$divcountzeros.png
    if [ -f $fname ]; then
        echo $fname
        echo $fnamenew
        ((divcounter++))
        mv $fname $fnamenew
    fi
    ((counter++))
    done
echo All done

ffmpeg -y -framerate 10 -i ./plots/b%04d.png -b:v 10M emerge.mp4
