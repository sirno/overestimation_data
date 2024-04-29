#!/usr/bin/env python
"""Submit beast.

Usage: submit_beast.py <config>...
"""

import argparse
import datetime
import fnmatch
import itertools
import logging
import os
import shutil
import subprocess
from tempfile import TemporaryDirectory

import pandas as pd
import phynalysis.cli as phyn_cli
import yaml
from joblib import Parallel, delayed
from phynalysis.configs import BeastJobConfig
from phynalysis.remote import PhynalysisRemoteClient
from rich.logging import RichHandler
from tqdm import tqdm

USER = os.environ.get("SSH_USER")
HOST = os.environ.get("SSH_HOST")

MAX_JOB_ARRAY_SIZE = 10000

LOCAL_BASE = "out/"
TEMPLATE_BASE = "data/"
REMOTE_BASE = "phages/simulations"

remote_client = None


def init_remote_client():
    """Init remote client if it does not exist."""
    global remote_client
    if remote_client is None:
        logging.info("[dim]Creating remote client.")
        remote_client = PhynalysisRemoteClient(HOST, USER)


def get_config_path(config: BeastJobConfig):
    """Get config path.

    This path is the relative path to the `run.xml` file."""
    return os.path.join(LOCAL_BASE, config.get_config_path(), "run.xml")


def get_sample_path(config: BeastJobConfig):
    """Get path to the sample data.

    This is the relative path to the `samples.haplotypes.csv` file. It identifies the
    data that is used for the specific BEAST run.
    """
    return os.path.join(LOCAL_BASE, config.sample, "samples.haplotypes.csv")


def get_simulation_path(config: BeastJobConfig):
    """Get simulation path."""
    return os.path.join(LOCAL_BASE, config.sample)


def get_template_path(config: BeastJobConfig):
    """Get template path."""
    return os.path.join(TEMPLATE_BASE, config.template, "template.xml")


def get_local_config_path(config: BeastJobConfig):
    """Get local path."""
    return os.path.join(LOCAL_BASE, config.get_config_path())


def get_local_run_path(config: BeastJobConfig):
    """Get local path."""
    return os.path.join(LOCAL_BASE, config.get_run_path())


def get_remote_path(config: BeastJobConfig):
    """Get remote path."""
    return os.path.join(REMOTE_BASE, config.get_run_path())


def has_config(config: BeastJobConfig):
    """Check if config exists."""
    return os.path.exists(get_config_path(config))


def has_remote_config(config: BeastJobConfig):
    """Check if config exists on cluster."""
    return remote_client.exists(os.path.join(get_remote_path(config), "run.xml"))


def has_remote_log(config: BeastJobConfig):
    return remote_client.exists(os.path.join(get_remote_path(config), "run.log"))


def has_local_log(config: BeastJobConfig):
    return os.path.exists(os.path.join(get_local_run_path(config), "run.log"))


def has_sample(config: BeastJobConfig):
    """Check if sample exists."""
    return os.path.exists(get_sample_path(config))


def has_template(config: BeastJobConfig):
    """Check if template exists."""
    return os.path.exists(get_template_path(config))


def expand_configs(configs):
    """Expand virolution configurations."""
    return itertools.chain.from_iterable(config.expand_paths() for config in configs)


def check_configs(configs: list[BeastJobConfig]):
    """Evaluate which configs need to be submitted."""
    configs = list(configs)
    logging.warning(f"Checking {len(configs)} configs.")
    logging.debug(f"Checking configs: {configs}")

    missing_samples = set(filter(lambda config: not has_sample(config), configs))
    missing_templates = set(filter(lambda config: not has_template(config), configs))

    missing_config = set(filter(lambda config: not has_config(config), configs))

    logging.warning(f"Missing samples: {len(missing_samples)}")
    logging.debug(f"Missing samples: {missing_samples}")
    logging.warning(f"Missing templates: {len(missing_templates)}")
    logging.debug(f"Missing templates: {missing_templates}")
    logging.warning(f"Missing configs: {len(missing_config)}")
    logging.debug(f"Missing configs: {missing_config}")

    return missing_config - missing_samples - missing_templates


def generate_config(beast_config: BeastJobConfig, reference: str):
    """Generate config file."""
    config_path = get_config_path(beast_config)
    local_path = get_local_config_path(beast_config)

    sample_path = get_sample_path(beast_config)
    template_path = get_template_path(beast_config)

    simulation_path = get_simulation_path(beast_config)

    logging.info("[bold blue]Generating `run.xml`...")
    logging.info(config_path)
    os.makedirs(local_path, exist_ok=True)

    # load data
    data = pd.read_csv(sample_path)

    # data preprocessing: filter, sample
    logging.info(f"Running phyn filter with query `{beast_config.query}`.")
    data = phyn_cli.filter(data, beast_config.query + "and count > 0", False)

    logging.info(f"Running phyn sample with {beast_config.n_samples} samples.")
    data = phyn_cli.choose_random(data, beast_config.n_samples, warnings=False)

    data.to_csv(os.path.join(local_path, "samples.haplotypes.csv"), index=False)

    # compute maximum likelihood tree
    if beast_config.tree == "ml" and not has_config(
        os.path.join(local_path, "ml_tree.nwk")
    ):
        logging.debug("[bold blue]Generating maximum likelihood `tree.nwk`...")
        logging.debug("Running phyn convert to phylip.")
        phyn_convert_phylip = subprocess.Popen(
            [
                "phyn",
                "convert",
                os.path.join(local_path, "samples.haplotypes.csv"),
                "--format",
                "phylip",
                "--reference",
                reference,
                "--output",
                os.path.join(local_path, "run.phylip"),
            ],
            stdin=subprocess.PIPE,
        )

        if phyn_convert_phylip.wait():
            raise RuntimeError("Phyn convert to phylip failed.")

        phyml = subprocess.Popen(
            [
                "phyml",
                "-i",
                os.path.join(local_path, "run.phylip"),
            ],
            stdin=subprocess.PIPE,
        )

        if phyml.wait():
            raise RuntimeError("Phyml failed.")

        os.rename(
            os.path.join(local_path, "run.phylip_phyml_tree.txt"),
            os.path.join(local_path, "ml_tree.nwk"),
        )

    # create beast2 xml
    logging.debug(f"Running phyn convert with xml template `{template_path}`.")
    phyn_convert_xml = subprocess.Popen(
        [
            "phyn",
            "convert",
            os.path.join(local_path, "samples.haplotypes.csv"),
            "--format",
            "xml",
            "--reference",
            reference,
            "--template",
            f"xml={template_path}",
            "--output",
            os.path.join(local_path, "run.xml"),
        ],
        stdin=subprocess.PIPE,
    )

    if phyn_convert_xml.wait():
        raise RuntimeError("Phyn convert to xml failed.")

    # add non-random tree to xml
    if beast_config.tree != "random":
        logging.debug(f"[bold blue]Inserting `{beast_config.tree}` tree...")
        # map tree type to directory
        tree_dir = {
            "ml": local_path,
            "ancestry": simulation_path,
            "genealogy": simulation_path,
        }

        # check if tree type is valid
        if beast_config.tree not in tree_dir:
            raise RuntimeError(
                f"Unknown tree type `{beast_config.tree}`. Can be `random`, `ml` or `true`."
            )

        phyn_insert_tree = subprocess.Popen(
            [
                "phyn",
                "beast-xml",
                "insert-tree",
                os.path.join(local_path, "run.xml"),
                os.path.join(
                    tree_dir[beast_config.tree], f"{beast_config.tree}_tree.nwk"
                ),
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        if phyn_insert_tree.wait():
            # this error should does not block further execution
            # cleanup then return an error code
            shutil.copy(os.path.join(local_path, "run.xml"), "run.xml")
            os.remove(os.path.join(local_path, "run.xml"))
            logging.error(f"Error in: {local_path}")
            return 1

    # remove tree operators and likelihood
    if beast_config.tree == "genealogy" or beast_config.tree == "ancestry":
        with open(os.path.join(local_path, "run.xml"), "r+") as xml_file:
            content = xml_file.read()
            content = content.replace("tree-operators-->", "")
            content = content.replace("<!--/tree-operators/", "")
            content = content.replace("likelihood-->", "")
            content = content.replace("<!--/likelihood/", "")
            xml_file.seek(0)
            xml_file.write(content)
            xml_file.truncate()

    logging.debug(f"[bold green]Generated `run.xml` config in `{local_path}`.")

    return 0


def send_configs(configs: set[BeastJobConfig]):
    """Package and send configs to cluster.

    The function will create a directory structure of all run directories in a temporary
    directory. The directory structure will be compressed and sent to the cluster and
    extracted there.
    """

    def get_files(config):
        config_path = config.get_config_path()
        run_path = config.get_run_path()
        return list(
            (os.path.join(config_path, file_name), os.path.join(run_path, file_name))
            for file_name in os.listdir(get_local_config_path(config))
            if os.path.isfile(os.path.join(get_local_config_path(config), file_name))
        )

    with TemporaryDirectory() as tmpdir:
        # create copy list
        copy_list = list(map(get_files, configs))

        # create run tree directory
        os.makedirs(os.path.join(tmpdir, "configs"), exist_ok=True)
        for origin, target in itertools.chain.from_iterable(copy_list):
            target_path = os.path.join(tmpdir, "configs", target)
            os.makedirs(
                os.path.join(tmpdir, "configs", os.path.dirname(target_path)),
                exist_ok=True,
            )
            shutil.copy(
                os.path.join(LOCAL_BASE, origin),
                target_path,
            )

        # create tarball
        # create content list
        content_list = list(
            map(lambda x: x[1], itertools.chain.from_iterable(copy_list))
        )
        content_list_path = os.path.join(tmpdir, "content_list.txt")
        with open(content_list_path, "w") as f:
            f.write("\n".join(content_list) + "\n")

        # create tarball
        local_archive_path = os.path.join(tmpdir, "configs.tar.gz")
        compress_cmd = [
            "gtar",
            "-czf",
            local_archive_path,
            "--directory",
            os.path.join(tmpdir, "configs"),
            "--files-from",
            content_list_path,
        ]
        logging.info("Compressing configs...")
        logging.debug(f"[dim] {' '.join(compress_cmd)}")
        subprocess.run(
            compress_cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # send tarball
        remote_archive_path = os.path.join(REMOTE_BASE, ".cache/configs.tar.gz")
        logging.info("Sending configs...")
        remote_client.put(local_archive_path, remote_archive_path)

        # extract tarball
        extract_cmd = [
            "tar",
            "-xzf",
            remote_archive_path,
            "--directory",
            REMOTE_BASE,
        ]
        logging.info("Extracting configs...")
        logging.debug(f"[dim] {' '.join(extract_cmd)}")
        remote_client.execute(" ".join(extract_cmd), print_output=False)


def batch_configs(configs: set[BeastJobConfig]):
    """Batch submit configs using a slurm job array."""
    random_id = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    with TemporaryDirectory() as tmpdir:
        local_task_file = os.path.join(tmpdir, f"beast_{random_id}.txt")

        # create task file
        get_path = lambda config: config.get_run_path()
        with open(local_task_file, "w") as f:
            f.write("\n".join(map(get_path, configs)))

        # transfer task file
        remote_task_file = os.path.join(
            REMOTE_BASE,
            ".cache/",
            f"beast_{random_id}.txt",
        )
        logging.info("Transferring task file...")
        remote_client.put(local_task_file, remote_task_file)

        # check if job array is too large
        # if it is split it into smaller parts
        # and create a schedule file
        start = 0
        splits = []
        while start < len(configs):
            end = min(start + MAX_JOB_ARRAY_SIZE, len(configs)) - 1
            splits.append(f"{start}-{end}")
            start += MAX_JOB_ARRAY_SIZE

        if splits[1:]:
            local_schedule_file = "out/cache/schedule.txt"
            with open(local_schedule_file, "w") as f:
                f.write(f".cache/beast_{random_id}.txt\n")
                f.write("\n".join(splits[1:]))

        # submit job array
        beast_job_array_cmd = f"cd phages/simulations && just beast_array .cache/beast_{random_id}.txt {splits[0]}"
        logging.info("Submitting job array...")
        logging.debug(f"[dim] {beast_job_array_cmd}")
        remote_client.execute(beast_job_array_cmd)


def batch_schedule(args):
    """Batch scheduled configs."""
    if not os.path.exists("out/cache/schedule.txt"):
        logging.error("No schedule file found.")
        return

    # read schedule file
    with open("out/cache/schedule.txt", "r") as f:
        schedule = f.readlines()

    # submit job array
    init_remote_client()
    task_file = schedule[0].strip()
    splits = schedule[1:]
    split = list(map(int, splits[0].strip().split("-")))
    offset = split[0]
    end = split[1] - offset
    beast_job_array_cmd = (
        f"cd phages/simulations && just beast_array {task_file} 0-{end} {offset}"
    )
    logging.info("Submitting job array...")
    logging.debug(f"[dim] {beast_job_array_cmd}")

    if remote_client.execute(beast_job_array_cmd):
        logging.error("Job array submission failed.")
        return

    if not splits[1:]:
        logging.info("All jobs submitted.")
        os.remove("out/cache/schedule.txt")
        return

    # update schedule file
    logging.info("Updating schedule file.")
    with open("out/cache/schedule.txt", "w") as f:
        f.write(task_file + "\n")
        f.write("\n".join(splits[1:]))


def dispatch(args):
    """Dispatch beast jobs on cluster."""
    config = yaml.safe_load(open(args.runs, "r"))
    beast_configs = set(
        expand_configs(BeastJobConfig.from_dict(config) for config in config["beast"])
    )

    missing_configs = check_configs(beast_configs)
    existing_configs = beast_configs - missing_configs

    if args.jobs:
        n = max(args.jobs - len(existing_configs), 0)
        missing_configs = set(list(missing_configs)[:n])

    missing_configs_ls = list(missing_configs)
    missing_configs_df = pd.DataFrame(
        data={
            "config": missing_configs_ls,
            "path": map(lambda config: config.get_config_path(), missing_configs_ls),
        }
    )
    unique_missing_configs = missing_configs_df.drop_duplicates(subset="path")["config"]

    logging.warning("---- ---- ----")
    logging.warning(f"Generating {len(unique_missing_configs)} configs.")

    errors = Parallel(n_jobs=8)(
        map(
            lambda config: delayed(generate_config)(config, args.reference),
            # prevent tqdm from printing when there is nothing to generate
            tqdm(unique_missing_configs) if not unique_missing_configs.empty else [],
        )
    )

    if any(errors):
        logging.warning("Some configs could not be generated.")
        logging.warning("Abort execution.")
        error_samples = list(
            map(
                lambda x: x[1].sample,
                filter(lambda x: x[0], zip(errors, unique_missing_configs)),
            )
        )
        print(error_samples)
        exit(1)

    logging.warning("---- ---- ----")
    local_configs = set(filter(has_config, beast_configs))
    logging.warning(f"Found {len(local_configs)} local configs.")

    if args.no_submit:
        return

    if args.local:
        list(map(run_local, local_configs))
        return

    init_remote_client()

    finished_configs = set(filter(has_local_log, beast_configs))
    logging.warning(f"Found {len(finished_configs)} finished configs.")
    unfinished_configs = local_configs - finished_configs

    remote_configs = set(filter(has_remote_config, unfinished_configs))
    logging.warning(f"Found {len(remote_configs)} remote configs.")
    executing_configs = set(filter(has_remote_log, unfinished_configs))
    logging.warning(f"Found {len(executing_configs)} executing configs.")
    transmit_configs = unfinished_configs - remote_configs
    execute_configs = unfinished_configs - executing_configs - finished_configs

    logging.warning("---- ---- ----")

    if args.query:
        filter_func = lambda config: fnmatch.fnmatch(
            config.get_config_path(),
            args.query,
        )
        transmit_configs = set(filter(filter_func, transmit_configs))
        execute_configs = set(filter(filter_func, execute_configs))

    if args.jobs:
        transmit_configs = list(transmit_configs)[: args.jobs]
        execute_configs = list(execute_configs)[: args.jobs]

    logging.warning(f"Submitting {len(transmit_configs)} jobs.")
    logging.debug(f"Submitting configs: {transmit_configs}")
    send_configs(transmit_configs)

    if args.missing:
        # report missing jobs and exit
        execute_configs_paths = [config.get_config_path() for config in execute_configs]
        logging.warning(f"Executing configs: {execute_configs_paths}")
        return

    logging.warning(f"Executing {len(execute_configs)} jobs.")
    logging.debug(f"Executing configs: {execute_configs}")

    batch_configs(execute_configs)


def run_local(beast_config: BeastJobConfig):
    logging.info("[bold blue]Local mode.")

    local_path = get_local_config_path(beast_config)

    def execute(cmd):
        ps = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        print(f"Running `{cmd}`.")
        for stdout_line in iter(ps.stdout.readline, ""):
            yield stdout_line
        ps.stdout.close()
        return_code = ps.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    for line in execute(["beast", "-working", os.path.join(local_path, "run.xml")]):
        print(line, end="")

    return


def submit(args, beast_config: BeastJobConfig):
    """Submit beast config to cluster."""
    config_path = beast_config.get_run_path()
    local_path = get_local_config_path(beast_config)
    remote_path = get_remote_path(beast_config)

    logging.info(f"[bold blue]Submit `{config_path}`...")

    def execute(cmd):
        ps = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        print(f"Running `{cmd}`.")
        for stdout_line in iter(ps.stdout.readline, ""):
            yield stdout_line
        ps.stdout.close()
        return_code = ps.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    if args.local:
        logging.info("[bold blue]Local mode.")
        for line in execute(["beast", "-working", os.path.join(local_path, "run.xml")]):
            print(line, end="")
        return

    init_remote_client()

    if has_remote_config(beast_config):
        logging.info(
            "[bold red]Found existing remote config, already running on cluster."
        )
        return

    logging.info("[bold blue]Send config to host.")
    logging.debug(f"Creating remote directory `{remote_path}`.")
    if not args.dry_run:
        remote_client.mkdir(remote_path)

    logging.debug(f"Copying config to `{remote_path}`.")
    if not args.dry_run:
        for file_name in os.listdir(local_path):
            local_file_path = os.path.join(local_path, file_name)
            remote_file_path = os.path.join(remote_path, file_name)
            remote_client.put(local_file_path, remote_file_path)

    logging.info("[bold blue]Submit batch on host.")
    if not args.dry_run:
        remote_client.execute(
            f"cd phages/simulations && just beast {config_path} {beast_config.beast_seed}",
            print_output=False,
        )

    logging.info("[bold green]Successfully submitted batch on host.")


def entry():
    parser = argparse.ArgumentParser(description="Dispatch beast jobs on cluster.")

    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="command",
        required=True,
    )

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        help="Verbose output.",
        default=0,
    )
    common_parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force config.",
    )
    common_parser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        help="Dry run.",
    )
    common_parser.add_argument(
        "--no-submit",
        action="store_true",
        help="Do not submit any code anywhere.",
    )
    common_parser.add_argument(
        "-l",
        "--local",
        action="store_true",
        help="Run locally only.",
    )

    submit_parser = subparsers.add_parser(
        "submit",
        help="Submit beast on cluster.",
        parents=[common_parser],
    )
    submit_parser.add_argument("reference", type=str, help="Path to reference genome.")
    submit_parser.set_defaults(func=submit)

    dispatch_parser = subparsers.add_parser(
        "dispatch",
        help="Dispatch beast on cluster.",
        parents=[common_parser],
    )
    dispatch_parser.add_argument("runs", type=str, help="Path to run configuration.")
    dispatch_parser.add_argument(
        "reference",
        type=str,
        help="Path to reference genome.",
    )
    dispatch_parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        help="Number of jobs to dispatch.",
    )
    dispatch_parser.add_argument(
        "-q",
        "--query",
        type=str,
        help="Query to filter jobs.",
    )
    dispatch_parser.add_argument(
        "--missing",
        action="store_true",
        help="Report missing jobs.",
    )
    dispatch_parser.set_defaults(func=dispatch)

    schedule_parser = subparsers.add_parser(
        "schedule",
        help="Execute scheduled jobs.",
        parents=[common_parser],
    )
    schedule_parser.set_defaults(func=batch_schedule)

    args = parser.parse_args()

    log_level = [logging.WARNING, logging.INFO, logging.DEBUG][min(args.verbose, 2)]

    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        handlers=[
            RichHandler(
                markup=True,
                show_path=False,
                show_time=False,
                show_level=False,
            )
        ],
    )

    args.func(args)


if __name__ == "__main__":
    entry()
