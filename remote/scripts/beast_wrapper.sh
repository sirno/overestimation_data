#!/usr/bin/env bash
# beast wrapper for euler
#
# USAGE: beast_wrapper.sh <dir> [<seed>]

# setup environment
module load openjdk/14.0.2

# set formatting variables
bold=$(tput -T xterm bold)
normal=$(tput -T xterm sgr0)

# read command line arguments
DIR=$1
SEED=$2

# specify directories
LOCAL=${PWD}/${DIR}
CACHE=/scratch/${DIR}
OUTDIR=${LOCAL}

BEAST_FLAGS=(-working)

STORE_TREES=false

# set seed if it is supplied
if [ ! -z "$SEED" ] ; then
    BEAST_FLAGS+=(-seed $SEED)
fi

# create directory structure
echo "${bold}Creating cache directory: ${CACHE}${normal}"
mkdir -p $CACHE

# continue run if it already exists
echo "${bold}Checking if there is a preexisting run${normal}"
if [ -f $OUTDIR/run.xml.state ] ; then
    echo "Resume existing run"
    cp -r $OUTDIR/* $CACHE
    BEAST_FLAGS+=(-resume)
else
    echo "Initiate new run"
    cp $LOCAL/run.xml $CACHE
fi

# start beast
echo "${bold}Running beast${normal}"
BEAST=beast$(rg -o "version=\"([2-9]+)\.([6-7]+)\"" ${CACHE}/run.xml --replace \$1\$2)
srun $BEAST ${BEAST_FLAGS[*]} ${CACHE}/run.xml &
pid=$!

# create exit functions
function store_local_scratch {
    if [[ $STORE_TREES == "false" ]] ; then
        echo "Removing trees files"
        rm ${CACHE}/*.trees
    fi
    echo "Storing local scratch to $OUTDIR"
    mkdir -p ${OUTDIR}
    cp -r ${CACHE}/* ${OUTDIR}
}

function purge_local_scratch {
    echo "Purging local scrath: $CACHE"
    rm -r ${CACHE}
}

function fix_trees_file {
    if [ "$(tail -n1 $1)" != "End;" ]; then
        echo "Fixing trees file: $1"
        echo "End;" >> $1
    fi
}

function fix_trees_files {
    for trees in $(find $OUTDIR -name "*.trees"); do
        fix_trees_file $trees
    done
}

# specify signal handling
function handle_usr1 {
    echo "${bold}Received SIGUSR1!${normal}"
    store_local_scratch
    fix_trees_files
}

function handle_term {
    echo "${bold}Received SIGTERM!${normal}"
    kill $pid
    wait $pid
    store_local_scratch
    purge_local_scratch
    fix_trees_files
}

trap handle_usr1 USR1
trap handle_term TERM

# wait for beast to finish
# we need to loop here, because otherwise we cannot continue to run
# after handling USR1
while ps -p $pid &> /dev/null ; do
    wait $pid
    exit_code=$?
done

echo "${bold}BEAST2 finished!"

# if cache has not been flushed - flush it
if [ -d $CACHE ] ; then
    echo "${normal} Saving non-empty cache."
    store_local_scratch
    purge_local_scratch
fi

# copy stdout and stderr to scratch
#cp -r ${LOCAL}/* ${OUTDIR}

