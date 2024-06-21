#!/bin/bash

#!/bin/bash

cat H3K4me3_H3K27ac H3K27ac_H3K27ac H3K4me3_H3K4me3  | cut -f1-2 | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > ../merged.EP
cat *CTCF* | cut -f1-2 | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > ../merged.CTCF

cat *H3K27me3* | cut -f1-2 | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > ../merged.H3K27me3

cat *H3K36me3* | cut -f1-2 | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > ../merged.H3K36me3
