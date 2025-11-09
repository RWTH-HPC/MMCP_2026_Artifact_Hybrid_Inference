#!/bin/bash

if [[ ${BASH_SOURCE[1]} == "" ]]
then
  echo "You should source this script instead of calling it, otherwise it won't work."
  exit 1
fi

# check if time is specified by the user
if [[ $1 != "" ]]
then
    runtime=$1
else
  # getting the full path of the script sourcing this one
  calling_script=$(readlink -f ${BASH_SOURCE[1]})
  # grep for a time in hh:mm:ss format which should be the maximum runtime of the job
  runtime=$(grep -o "[0-9]\{1,2\}:[0-9]\{1,2\}:[0-9]\{1,2\}" "$calling_script" | head -n 1)
  if [[ $runtime == "" ]]
  then
    echo "Could not detect runtime please specify it either in the hh:mm, hh:mm:ss or ssss format."
    exit 1
  fi
fi

if [[ $(echo "$runtime" | grep -o "[0-9]\{1,2\}:[0-9]\{1,2\}:[0-9]\{1,2\}") != "" ]]
then
  # runtime in hh:mm:ss format
  runtime_end=$(($(echo "$runtime" | awk -F: '{print ($1 * 3600) + ($2 * 60) + $3}') + $(date +%s)))
elif [[ $(echo "$runtime" | grep -o "[0-9]\{1,2\}:[0-9]\{1,2\}") != "" ]]
then
  # runtime in hh:mm format
  runtime_end=$(($(echo "$runtime" | awk -F: '{print ($1 * 3600) + ($2 * 60)}') + $(date +%s)))
elif [[ $(echo "$runtime" | grep -o "[0-9]") != "" ]]
then
  # runtime in seconds
  runtime_end=$(($runtime + $(date +%s)))
else
  echo "Please specify the runtime either the hh:mm, hh:mm:ss or ssss format."
  exit 1
fi

# export end time of the job
export ZFS_JOB_END_TIME=$runtime_end
echo "Activated autosave feature."
