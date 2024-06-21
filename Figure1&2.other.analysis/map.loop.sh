# This is the demo.

#!/bin/bash

ref=$1
mergefile=$2
k=$3
cut=$4
#ln -s ../wt.merged
#ln -s ../ko.merged

#1. lifted chromatin loop from one Hi-C platform to another Hi-C platform

cat gain.loops | awk 'NR>1' | sed 's/:/\t/g' | cut -f1-2 | python ~/Ren/WT_2days/2fold/frag.2.anchor.py $ref  > gain.loops.lifted 
cat lost.loops | awk 'NR>1' | sed 's/:/\t/g' | cut -f1-2 | python ~/Ren/WT_2days/2fold/frag.2.anchor.py $ref  > lost.loops.lifted
cat common.loops | awk 'NR>1' | sed 's/:/\t/g' | cut -f1-2 | python ~/Ren/WT_2days/2fold/frag.2.anchor.py $ref  > common.loops.lifted


#2. extract the loop strength of lifted loop from corresponding background
python ~/mESC_wapl/to_nora/merged.py wt ko gain.loops.lifted 
python ~/mESC_wapl/to_nora/merged.py wt ko lost.loops.lifted
python ~/mESC_wapl/to_nora/merged.py wt ko common.loops.lifted

#3. map lifted loop to its corresponding background. 

Rscript ~/test_1129/03222023_script/map.crop.frame.r $k $cut $mergefile

#4. create the inlet of gain loop
Rscript ~/degradation/new_CTCF/3.16.2021.comp/crop_frame.r 1 5 gain.loops.lifted.lifted.ratio 
