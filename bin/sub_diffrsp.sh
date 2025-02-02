#!/bin/bash

# Check if the correct number of parameters is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <start_time_seconds> <end_time_seconds>"
    exit 1
fi

# Input start and end times in seconds
start_time="$1"
end_time="$2"

# Calculate the total duration
total_duration=$((end_time - start_time))

# Check if the end time is after the start time
if [ "$total_duration" -le 0 ]; then
    echo "End time must be greater than start time."
    exit 1
fi

# Calculate the duration of each sub-interval
interval_duration=$((total_duration / 100))

# Loop to submit 100 jobs
for i in $(seq 0 99); do
    # Calculate start and stop time for each job
    sub_start_time=$((start_time + i * interval_duration))
    sub_end_time=$((start_time + (i + 1) * interval_duration))

    # Submit the sbatch job with start and end times as parameters
    sbatch sb_run.sh "$sub_start_time" "$sub_end_time"
done

echo "Submitted 100 jobs from $start_time to $end_time."
