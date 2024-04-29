#!/usr/bin/env python

import argparse
import os
import subprocess
import sys


def build_command(args, log_path):
    cmd = ["sbatch"]

    # set array arguments
    if args.array:
        cmd += [
            f"--array={args.array}",
        ]

    # set resource arguments
    cmd += [
        f"--time={args.time}",
        f"--cpus-per-task={args.cores}",
        f"--mem-per-cpu={args.mem}",
        f"--tmp={args.tmp}",
    ]

    # set output arguments
    cmd += [
        f"--job-name={args.directory}",
        f"--output={log_path}",
    ]

    # set job script
    cmd += [
        "scripts/virolution_wrapper.sh",
        args.directory,
        str(args.generations),
        str(args.compartments),
        args.sequence,
    ]

    # set replicates
    cmd += [
        str(args.n_replicates),
    ]

    return cmd


def main(args):
    """Main."""
    # check arguments
    if args.array and args.n_replicates > 1:
        print("Cannot run array with multiple replicates.")
        sys.exit(1)

    # setup multiple replicates
    if args.n_replicates > 1:
        args.array = f"[1-{args.n_replicates}]"

    # setup path
    log_path = "logs/virolution_wrapper"
    os.makedirs(log_path, exist_ok=True)
    if args.array:
        log_path = os.path.join(log_path, "job_%A_%a.out")
    else:
        log_path = os.path.join(log_path, "job_%j.out")

    # submit batch
    p = subprocess.run(build_command(args, log_path), stdout=sys.stdout, stderr=sys.stderr)
    sys.exit(p.returncode)


def entry():
    """Entry."""
    parser = argparse.ArgumentParser(description="Submit virolution jobs on cluster.")
    parser.add_argument("directory", type=str, help="Path to config directory.")
    parser.add_argument(
        "--generations",
        "-g",
        type=int,
        default=1000,
        help="Number of generations. Default: 1000.",
    )
    parser.add_argument(
        "--compartments",
        "-c",
        type=int,
        default=1,
        help="Number of compartments. Default: 1.",
    )
    parser.add_argument(
        "--sequence",
        "-s",
        type=str,
        help="Path to sequence fasta file.",
    )
    parser.add_argument(
        "--cores",
        "-n",
        type=int,
        default=1,
        help="Number of cores. Default: 1.",
    )
    parser.add_argument(
        "--time",
        "-t",
        type=str,
        default="04-00",
        help="Time limit. Default: 4 days.",
    )
    parser.add_argument(
        "--mem",
        "-m",
        type=str,
        default="8G",
        help="Memory limit. Default: 8G.",
    )
    parser.add_argument(
        "--tmp",
        "-T",
        type=str,
        default="20G",
        help="Local scratch memory limit. Default: 20G.",
    )
    parser.add_argument(
        "--array",
        "-a",
        type=str,
        default="",
        help="Run job array.",
    )
    parser.add_argument(
        "--n-replicates",
        "-r",
        type=int,
        default=1,
        help="Number of replicates to be submitted.",
    )

    args = parser.parse_args()

    main(args)


if __name__ == "__main__":
    entry()
