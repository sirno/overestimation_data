#!/usr/bin/env bash
# beast array wrapper for euler
#
# USAGE: beast_array_wrapper.sh <tasklist> <offset>

tasklist=$1
tasklist_lines=(`cat $1`)
task_id=$(($SLURM_ARRAY_TASK_ID + $2))
task=${tasklist_lines[$task_id]}

seed=$(echo $task | rg -o bs_[0-9]+ | rg -o [0-9]+)

echo "Running tasklist: $tasklist"
echo "Run task id: $task_id"
echo "Run task: $task"
echo "Using beast seed: $seed"

source ./scripts/beast_wrapper.sh $task $seed

