#!/bin/sh
echo "Executing baobabs.R script in the background"
Rscript --vanilla R/baobabs.R > baobabs.log 2>&1 &
echo "Check the progress with command 'tail -f baobabs.log'"
echo "Check the processor usage with command 'top'"
## End of script