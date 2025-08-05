import fnmatch
import hashlib
import os
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
from phynalysis.beast import read_beast_log
from phynalysis.configs import BeastJobConfig


# %%
def is_notebook() -> bool:
    try:
        from IPython.core.getipython import get_ipython

        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False  # Probably standard Python interpreter
    except ModuleNotFoundError:
        return False


# %%
def fast_scandir(dirname, exclude=None):
    """Recursively scan for run.log files in subdirectories.

    Args:
        dirname (str): The directory path to start scanning from.

    Returns:
        list: A list of paths to run.log files found within the directory and its subdirectories.
    """
    matches = []
    for root, dirnames, filenames in os.walk(dirname):
        for filename in fnmatch.filter(filenames, "run.log"):
            if exclude is not None and any(e in root for e in exclude):
                continue
            matches.append(os.path.join(root, filename))
    return matches


def get_config(path):
    config = BeastJobConfig.from_path(path)
    selection = config.sample.split("/")[-3]
    migration = config.sample.split("/")[-2].split("_")[1]
    pop_size = config.sample.split("/")[-4]

    return {
        "path": path,
        **BeastJobConfig.to_dict(config),
        "selection": selection,
        "migration_rate": float(migration),
        "pop_size": pop_size,
    }


key_errors = []
n_key_errors = 0
empty_data_errors = []


def get_migration_rates(row):
    try:
        trace = read_beast_log(row.path, burn_in=0.1)
        model_name = "/".join(row.template.split("/")[0:2])
        fields = MIGRATION_FIELDS[model_name]
        migration = trace[fields[0]].values
        hpi = az.stats.hdi(migration)
        ess = az.ess(migration).reshape(-1)
        n_steps = np.array(migration.shape[0], ndmin=1)
        if migration.shape[0] <= 100:
            print(f"Warning: {row.path} has only {migration.shape[0]} steps.")
        overestimation = np.mean(migration > row.migration_rate).reshape(-1)
        s = pd.Series(
            np.concatenate(
                (trace[fields].mean().values, hpi, ess, n_steps, overestimation),
                axis=0,
            ),
            index=["mean", "lower", "upper", "ess", "n_steps", "overestimation"],
        )
        return s
    except KeyError as e:
        print(e)
        global key_errors
        global n_key_errors
        key_errors.append(e)
        n_key_errors += 1
        return pd.Series(
            [np.nan] * 6,
            index=["mean", "lower", "upper", "ess", "n_steps", "overestimation"],
        )
    except pd.errors.EmptyDataError:
        global empty_data_errors
        empty_data_errors.append(row.path)
        os.remove(row.path)
        return pd.Series(
            [np.nan] * 6,
            index=["mean", "lower", "upper", "ess", "n_steps", "overestimation"],
        )
    except:
        print("Error in", row.path)
        raise


def generate_hash(sorted_strings):
    # Convert the list of strings to a single string
    combined_string = "".join(sorted_strings)

    # Create a hash object
    hasher = hashlib.sha256()

    # Update the hash object with he combined string
    hasher.update(combined_string.encode("utf-8"))

    # Get the hexadecimal representation of the hash
    hash_value = hasher.hexdigest()

    # Return the hash value
    return hash_value


def get_function_hash(func):
    # Get the bytecode of the function
    bytecode = func.__code__.co_code

    # Calculate the hash using the MD5 algorithm
    md5_hash = hashlib.md5(bytecode).hexdigest()

    return md5_hash


# %%
MIGRATION_FIELDS = {
    "beast/lemey": ["trait_clock_rate"],
    "beast/mascot": ["b_migration_constant.0_and_1"],
    "beast/bdmm": ["rateMatrix.1"],
}


def load_data(prefix, identifier, log_files, force=False, show_exponential=False):
    """This function loads the data from the log files and processes it.

    The output is a dataframe with following columns:
        - multiple columns with simulation and inference parameters
        - mean: the mean of the migration rate
        - lower: the lower bound of the 95% HPD interval
        - upper: the upper bound of the 95% HPD interval
        - ess: the effective sample size of the migration rate
        - n_steps: the number of steps in the trace
    """
    cache_path = Path(f"{prefix}/out/cache/{identifier}.csv")
    print(f"Attempt to load data from {cache_path.name}")

    # if the cache does not exist or force is True, process the data
    if force or not cache_path.exists():
        # process data and save it
        print(f"Processing {len(log_files)} files")
        configs_raw = (get_config(path) for path in log_files)
        configs = pd.DataFrame(configs_raw)
        inferred_rates = configs.progress_apply(get_migration_rates, axis=1)
        if empty_data_errors:
            print("--- --- ---")
            print("Empty data errors: ", len(empty_data_errors))
        configs[["mean", "lower", "upper", "ess", "n_steps", "overestimation"]] = (
            inferred_rates
        )
        configs.to_csv(f"{prefix}/out/cache/{identifier}.csv")
        print("--- --- ---")

    data = pd.read_csv(f"{prefix}/out/cache/{identifier}.csv", index_col=0)
    data["method"] = data.template.map(lambda x: x.split("/")[1])
    data["mutation_rate"] = data.template.map(
        lambda x: "2.16e-5" if "nmr" in x else "2.16e-4"
    )
    data["population_size"] = data.template.map(
        lambda x: 100 if "pop_tiny" in x else 1000
    )

    if not show_exponential:
        data = data.query("selection != 'exponential'")

    if data.pop_size.unique().shape[0] > 1:
        data = data.query("pop_size == 'tiny'")

    print("===================")
    print("LOADED DATA REPORT")
    print("===================")
    print(data.groupby(["template", "tree", "n_samples"]).size())

    return data


def write_data(data, prefix, identifier):
    data.to_csv(f"{prefix}/out/cache/{identifier}.csv")
