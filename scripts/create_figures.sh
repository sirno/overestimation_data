#!/usr/bin/env bash

# Use strict mode
set -euo pipefail

flags=""

while test $# != 0
do
    case "$1" in
    -l|--load-store) flags="--load-store" ;;
    -s|--save-store) flags="--save-store" ;;
    esac
    shift
done

mkdir -p figures

for clock in strict_clock relaxed_clock; do
  # Create figures
  for rate in nmr hmr; do
    for pop in pop_tiny pop_small; do
      python notebooks/overestimations.py --selector ${clock}/${rate}/${pop} ${flags};
    done;
  done;

  # Create tables with highlighted values
  python notebooks/tables.py --glob-type ${clock};
done;

