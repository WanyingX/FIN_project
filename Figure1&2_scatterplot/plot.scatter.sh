# Step1. Define chromatin loops by using fold change
# Let user give regression slope

k=$1
maxs=$2
Rscript getloop.r $k merged.loops

# Step2. Plot scatterplot

Rscript heatscatter.r $k $maxs merged.loops 
