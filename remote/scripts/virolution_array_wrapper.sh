#!/usr/bin/env bash
# virolution array wrapper for euler
#
# USAGE: beast_array_wrapper.sh <tasklist>

tasklist=$1
mapfile -t tasklist_lines < $tasklist
task=${tasklist_lines[$SLURM_ARRAY_TASK_ID]}

echo "virolution_array_wrapper.sh:"

echo "Running tasklist: $tasklist"
echo "Run task: $task"

source ./scripts/virolution_wrapper.sh $task

