#!/usr/bin/env bash

# Execute the command to run the MCMC
now=$(date +%s)
log="output_b5.log"
process_file="process_number_b5.log"
OperatingSystem=$(uname -s)
if [ "$OperatingSystem" = "Linux" ]; then
  nohup python blazar_run_mcmc.py > "$log" &
fi

if [ "$OperatingSystem" = "Darwin" ]; then
  nohup caffeinate python blazar_run_mcmc.py > "$log" &
fi

echo $! > "$process_file"
echo "MCMC command started"
