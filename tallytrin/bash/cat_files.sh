#!/bin/bash

# Define base string
base_command="cat seperate_samples.dir/*SampleNUM.fastq | gzip > UCL-MM-Samples1-12_SampleNUM.fastq.gz"

# Loop over 12 samples
for i in $(seq 1 12)
do
    # Replace 'NUM' with the current sample number
    command=${base_command//NUM/$i}
    
    # Execute the command
    echo "Running: $command"
    eval $command
done
