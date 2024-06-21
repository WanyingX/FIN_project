#!/bin/bash

expt=$1
lib=/mnt/rds/genetics01/JinLab/xww/0326_data/insulation_score/
bed_ref=/mnt/rds/genetics01/JinLab/xww/0326_data/insulation_score/mm10.20kb.bed
cat $expt | cut -f1-3 | python2 $lib/to_pre.py $bed_ref | sed 's/bin//g' |  sort -k2,2d -k6,6d   > $expt.pre
java -jar ~/software/juicer/juicer_ml_tools.jar pre -d -r 25000 -k KR $expt.pre $expt.hic mm10
fanc insulation ./$expt.hic \
                ./$expt.insulation \
                -w 1000000 \
                -o bed
