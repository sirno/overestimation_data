# %%
import argparse
import os
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

try:
    from IPython.display import Markdown
except ImportError:
    Markdown = lambda x: print(x)

from load_data import (
    MIGRATION_FIELDS,
    fast_scandir,
    generate_hash,
    get_function_hash,
    get_migration_rates,
    is_notebook,
    load_data,
)
from tqdm import tqdm

tqdm.pandas()


# %%
parser = argparse.ArgumentParser()
parser.add_argument(
    "--show_exponential",
    action="store_true",
    help="Show results for exponential selection",
)
parser.add_argument(
    "--selector",
    type=str,
    help="The selector to use",
)
parser.add_argument(
    "--save-store",
    action="store_true",
    help="Save data to store.",
)
parser.add_argument(
    "--load-store",
    action="store_true",
    help="Use data from store.",
)

if is_notebook():
    args = parser.parse_args(
        args=[
            "--selector",
            "relaxed_clock/hmr/pop_small",
            # "tests/no_geo_prior/pop_tiny",
        ]
    )
    prefix = ".."
else:
    args = parser.parse_args()
    prefix = "."
    plt.show = lambda: None

show_exponential = args.show_exponential
selector = args.selector
selector_label = selector.replace("/", "_")
force = False

trees = ["ancestry", "random"]
selections = ["neutral", "lognormal"]
methods = ["lemey", "mascot"]
migration_rates = [0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2]
# mutation_rates = ["2.16e-5", "2.16e-4"]

experiment_group = [
    "method",
    "mutation_rate",
    "population_size",
    "migration_rate",
    "selection",
    "tree",
]

if not is_notebook():
    print(f"## Selected: {selector}")
    print("--- --- ---")

Markdown(f"## Selected: {selector}")

base_paths = [f"{prefix}/out/beast/{method}/{selector}" for method in methods]
paths = [fast_scandir(base_path) for base_path in base_paths]
log_files = sorted(list(p for path_list in paths for p in path_list))

if args.save_store:
    identifier = selector_label
    force = True
elif args.load_store:
    identifier = selector_label
else:
    identifier = (
        get_function_hash(get_migration_rates) + "_" + generate_hash(log_files)[:12]
    )

data = load_data(
    prefix,
    identifier,
    log_files,
    force=force,
    show_exponential=show_exponential,
)


# %%
def get_query(method, tree, query="time_eq_1000"):
    return f"method == '{method}' and tree == '{tree}' and query == '{query}'"


def plot_migration_rates(template, data, ylim=None, xlim=None):
    # get extent of x axis
    xextent = data.migration_rate.min(), data.migration_rate.max()

    # create figure and axes
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # setup random tree plot
    axes[0].set_title("based on phylogenetic reconstruction")
    sns.lineplot(
        ax=axes[0],
        x="migration_rate",
        y="mean",
        data=data.query(
            f"template.str.startswith('{template}').values and query == 'time_eq_1000' and tree == 'random'",
        ),
        hue="selection",
        style="pop_size",
    )
    axes[0].plot(xextent, xextent, color="black", linestyle="--", alpha=0.7)
    axes[0].set_ylabel("inferred migration rate")

    # setup fixed tree plot
    axes[1].set_title("based on true genealogy")
    sns.lineplot(
        ax=axes[1],
        x="migration_rate",
        y="mean",
        data=data.query(
            f"template.str.startswith('{template}').values and query == 'time_eq_1000' and tree == 'ancestry'"
        ),
        hue="selection",
        style="pop_size",
    )
    axes[1].plot(xextent, xextent, color="black", linestyle="--", alpha=0.7)

    for ax in axes:
        ax.set_xlabel("migration rate")

    # Create a legend for the whole figure
    handles, labels = axes[0].get_legend_handles_labels()
    labels = [
        label.capitalize() if label != "lognormal" else "Selection" for label in labels
    ]
    handles[0] = Line2D([0], [0], color="black", linestyle="--", label="True migration")
    labels[0] = "Identity"
    fig.legend(
        # handles,
        # labels,
        handles[0:-2],
        labels[0:-2],
        loc="center right",
        bbox_to_anchor=(1, 0.5),
        borderaxespad=0.0,
    )

    # modify axes with common attributes
    for ax in axes:
        ax.legend().remove()
        if ylim is not None:
            ax.set_ylim(*ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)

    fig.tight_layout()
    fig.subplots_adjust(right=0.88)


def plot_migration_rates_all(data, alpha=0.1, xlim=None, ylim=None):
    xextent = data.migration_rate.min(), data.migration_rate.max()
    fig, axes = plt.subplots(2, 2, figsize=(12, 6), sharey=True)

    rows = ["lemey", "mascot"]
    cols = ["random", "ancestry"]

    row_labels = ["DTA", "MASCOT"]
    col_labels = ["based on phylogenetic reconstruction", "based on true genealogy"]

    for row_idx, row in enumerate(rows):
        for col_idx, col in enumerate(cols):
            sns.lineplot(
                ax=axes[row_idx, col_idx],
                x="migration_rate",
                y="mean",
                data=data.query(get_query(row, col)),
                hue="selection",
                style="pop_size",
            )
            axes[row_idx, col_idx].plot(
                xextent, xextent, color="black", linestyle="--", alpha=0.7
            )

    for ax in axes[:-1, :].flat:
        ax.set_xlabel("")
    for ax in axes[-1, :]:
        ax.set_xlabel("migration rate")
    for ax in axes[:, 0]:
        ax.set_ylabel("inferred migration rate")

    pad = 5
    for ax, col in zip(axes[0], col_labels):
        ax.annotate(
            col,
            xy=(0.5, 1),
            xytext=(0, pad),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="center",
            va="baseline",
        )

    for ax, row in zip(axes[:, 0], row_labels):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
        )

    # Create a legend for the whole figure
    handles, labels = axes[0, 0].get_legend_handles_labels()
    labels = [
        label.capitalize() if label != "lognormal" else "Selection" for label in labels
    ]
    handles[0] = Line2D([0], [0], color="black", linestyle="--", label="True migration")
    labels[0] = "Identity"
    fig.legend(
        # handles,
        # labels,
        handles[0:-2],
        labels[0:-2],
        loc="center right",
        bbox_to_anchor=(1, 0.5),
        borderaxespad=0.0,
    )

    for ax in axes.flat:
        # h, l = ax.get_legend_handles_labels()
        # ax.legend(h[:-2], l[:-2])
        ax.legend().remove()
        if ylim is not None:
            ax.set_ylim(*ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95, right=0.88)


# %% [markdown]
## Migration rates

# %%
plot_migration_rates_all(data, alpha=0.1, xlim=(0.0, 0.01), ylim=(0.0, 0.02))

os.makedirs(f"{prefix}/figures/migration_rates", exist_ok=True)
plt.savefig(
    f"{prefix}/figures/migration_rates/{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## Lemey

# %%
plot_migration_rates("beast/lemey", data, ylim=(0.0, 0.02), xlim=(0.0, 0.01))

plt.savefig(
    f"{prefix}/figures/migration_rates/lemey_{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## Mascot

# %%
plot_migration_rates("beast/mascot", data, ylim=(0.0, 0.02), xlim=(0.0, 0.01))

plt.savefig(
    f"{prefix}/figures/migration_rates/mascot_{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## BDMM

# %%
# plot_migration_rates("beast/bdmm", data, ylim=(0.0, 0.5))


# %%
# functions to plot the posterior distribution of the migration rate
# for each category of the data mean of means and HPD intervals are calculated
# the HPD intervals are then plotted as a shaded area around the mean


def plot_posterior_subplot(data, ax, alpha=0.1):
    plot_data = data.groupby(["migration_rate", "selection", "pop_size"])[
        ["mean", "lower", "upper"]
    ].mean()
    plot_data = plot_data.reset_index()

    xextent = data.migration_rate.min(), data.migration_rate.max()

    colors = sns.color_palette("colorblind", 3)
    colors_iter = iter(colors)

    for selection in plot_data.selection.unique():
        pld = plot_data.query(f"selection == '{selection}'")
        c = next(colors_iter)
        ax.plot(pld.migration_rate, pld["mean"], label=selection, color=c)
        ax.fill_between(
            pld.migration_rate, pld["lower"], pld["upper"], color=c, alpha=alpha
        )
        ax.plot(xextent, xextent, color="black", linestyle="--", alpha=0.5)


def plot_posterior(data, alpha=0.1, ylim=None):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    axes[0].set_title("based on phylogenetic reconstruction")
    plot_posterior_subplot(data.query("tree == 'random'"), axes[0], alpha=alpha)

    axes[1].set_title("based on true genealogy")
    plot_posterior_subplot(data.query("tree == 'ancestry'"), axes[1], alpha=alpha)

    if ylim is not None:
        axes[0].set_ylim(*ylim)
        axes[1].set_ylim(*ylim)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center right",
        bbox_to_anchor=(1, 0.5),
        borderaxespad=0.0,
    )

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95, right=0.88)


def plot_posterior_all(data, alpha=0.1, ylim=None):
    fig, axes = plt.subplots(2, 2, figsize=(12, 6), sharey=True)

    rows = ["lemey", "mascot"]
    cols = ["random", "ancestry"]

    row_labels = ["DTA", "MASCOT"]
    col_labels = ["based on phylogenetic reconstruction", "based on true genealogy"]

    for row_idx, row in enumerate(rows):
        for col_idx, col in enumerate(cols):
            plot_posterior_subplot(
                data.query(get_query(row, col)),
                axes[row_idx, col_idx],
                alpha=alpha,
            )

    for ax in axes[:-1, :].flat:
        ax.set_xlabel("")
    for ax in axes[-1, :]:
        ax.set_xlabel("migration rate")
    for ax in axes[:, 0]:
        ax.set_ylabel("posterior migration rate hdi")

    pad = 5
    for ax, col in zip(axes[0], col_labels):
        ax.annotate(
            col,
            xy=(0.5, 1),
            xytext=(0, pad),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="center",
            va="baseline",
        )

    for ax, row in zip(axes[:, 0], row_labels):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
        )

    if ylim is not None:
        for ax in axes.flat:
            ax.set_ylim(*ylim)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center right",
        bbox_to_anchor=(1, 0.5),
        borderaxespad=0.0,
    )

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95, right=0.88)


# %% [markdown]
## Migration rate posterior
plot_posterior_all(data, alpha=0.1, ylim=None)

os.makedirs(f"{prefix}/figures/migration_rate_hpd", exist_ok=True)
plt.savefig(
    f"{prefix}/figures/migration_rate_hpd/{selector_label}.pdf", bbox_inches="tight"
)

# %%

# %% [markdown]
## Lemey

# %%
template = "beast/lemey"
plot_posterior(
    data.query(f"template.str.startswith('{template}') and query == 'time_eq_1000'"),
    0.1,
    # ylim=(0.0, 0.2),
)

plt.savefig(
    f"{prefix}/figures/migration_rate_hpd/lemey_{selector_label}.pdf",
    bbox_inches="tight",
)

# %% [markdown]
## Mascot

# %%
template = "beast/mascot"
plot_posterior(
    data.query(f"template.str.startswith('{template}') and query == 'time_eq_1000'"),
    0.1,
    # ylim=(0.0, 0.2),
)

plt.savefig(
    f"{prefix}/figures/migration_rate_hpd/{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## BDMM

# %%
# template = "beast/bdmm"
# plot_posterior(
#     data.query(f"template.str.startswith('{template}') and query == 'time_eq_1000'"),
#     0.1,
#     # ylim=(0.0, 0.5),
# )

# %%
# plot the effective sample size of the migration rate by selection
sns.boxplot(data=data, x="template", y="ess", hue="selection")

# %%
# plot the total number of samples by tree
sns.boxplot(data=data, x="template", y="n_samples", hue="tree")

# %% [markdown]
## Statistical tests

### Wilcoxon signed-rank test of the mean inferred migration rate against the true migration rate

# %%
pvals = (
    data.groupby(
        [
            "method",
            "mutation_rate",
            "population_size",
            "selection",
            "tree",
            "migration_rate",
        ]
    )[["migration_rate", "mean"]]
    .apply(
        lambda group: 2
        * sp.stats.wilcoxon(group["mean"] - group["migration_rate"]).pvalue
        * (int(group["mean"].mean() >= group["migration_rate"].iloc[0]) - 0.5)
    )
    .rename("pval")
    .reset_index()
)

index = pd.MultiIndex.from_product(
    [
        methods,
        data.mutation_rate.unique(),
        data.population_size.unique(),
        selections,
        trees,
        migration_rates,
    ],
    names=[
        "method",
        "mutation_rate",
        "population_size",
        "selection",
        "tree",
        "migration_rate",
    ],
)

pvals = (
    pvals.set_index(
        [
            "method",
            "mutation_rate",
            "population_size",
            "selection",
            "tree",
            "migration_rate",
        ]
    )
    .reindex(index)
    .reset_index()
)

if not selector.startswith("tests"):
    os.makedirs(f"{prefix}/figures/pvals_mean", exist_ok=True)
    pvals.to_csv(f"{prefix}/figures/pvals_mean/{selector_label}.csv", index=False)

# %%
pvals


# %%
def plot_pvals(data, method):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    data = data.copy()
    data["pval"] = np.abs(data["pval"])
    axes[0].set_title("random tree")
    sns.scatterplot(
        data=data.query(f"method == '{method}' and tree == 'random'"),
        x="migration_rate",
        y="pval",
        hue="selection",
        ax=axes[0],
    )
    axes[1].set_title("fixed tree")
    sns.scatterplot(
        data=data.query(f"method == '{method}' and tree == 'ancestry'"),
        x="migration_rate",
        y="pval",
        hue="selection",
        ax=axes[1],
    )
    lower = min(data["pval"].min(), 1e-2)
    for ax in axes:
        ax.set_xlim(-0.0005, 0.0105)
        ax.fill_between(
            [-1, 1],
            [0.05, 0.05],
            [lower, lower],
            color="red",
            alpha=0.1,
            label="p < 0.05",
        )
        ax.set_xlabel("migration rate")
        ax.set_yscale("log")


# %%
os.makedirs(f"{prefix}/figures/misc", exist_ok=True)
plot_pvals(pvals, "lemey")
plt.savefig(
    f"{prefix}/figures/misc/lemey_pvals_{selector_label}.pdf", bbox_inches="tight"
)

# %%
plot_pvals(pvals, "mascot")
plt.savefig(
    f"{prefix}/figures/misc/mascot_pvals_{selector_label}.pdf", bbox_inches="tight"
)

# %%
# Visualize distribution of fitness effects empirically

"""yaml
distribution: !Exponential
  weights:
    beneficial: 0.29
    deleterious: 0.51
    lethal: 0.2
    neutral: 0.0
    lambda_beneficial: 0.03
    lambda_deleterious: 0.21
distribution: !Lognormal
    lethal: 0.045
    mu: -0.248
    sigma: 0.149
"""


def draw_exp(
    beneficial, deleterious, lethal, neutral, lambda_beneficial, lambda_deleterious
):
    random_number = random.random()
    if random_number < lethal:
        fitness = 0
    elif random_number < lethal + neutral:
        fitness = 1
    elif random_number < lethal + neutral + beneficial:
        fitness = 1 + np.random.exponential(lambda_beneficial)
    else:
        fitness = 1 - np.random.exponential(lambda_deleterious)
    if fitness < 0:
        fitness = 0

    return fitness


exp = [draw_exp(0.29, 0.51, 0.2, 0.0, 0.03, 0.21) for _ in range(10000)]

plt.hist(exp, bins=100, density=True)

plt.show()

# %%
# Draw MFEDs
x = np.linspace(0, 1, 200)

exp = (x < 1) * 0.51 * sp.stats.expon.pdf(1 - x, scale=0.21) + 0.29 * (
    x > 1
) * sp.stats.expon.pdf(x - 1, scale=0.03)

exp_beneficial = 0.29 * sp.stats.expon.pdf(x, scale=0.03)
exp_deleterious = 0.51 * sp.stats.expon.pdf(1 - x, scale=0.21)

plt.plot(1 + x, exp_beneficial, color="tab:blue", label="exponential", alpha=0.7)
plt.plot(x, exp_deleterious, color="tab:blue", alpha=0.7)
plt.plot(
    2 * x,
    sp.stats.lognorm.pdf(2 * x, s=0.149, loc=-0.248),
    color="tab:orange",
    label="lognormal",
    alpha=0.7,
)
plt.yticks([])
plt.xticks([0, 1, 2], ["0", "1", "2"])
plt.xlabel("fitness")
plt.vlines(1, 0, 10, color="tab:red", alpha=0.7)
plt.legend()
plt.show()
plt.savefig(f"{prefix}/figures/fitness_distributions.pdf", bbox_inches="tight")


# %%
def mean_squared_error(group):
    diff = group["mean"] - group["migration_rate"]
    return (diff**2).mean()


mse = (
    data.groupby(experiment_group)[experiment_group + ["mean"]]
    .apply(mean_squared_error)
    .rename("mse")
    .reset_index()
)
mse = mse.set_index(experiment_group).reset_index()

if not selector.startswith("tests"):
    os.makedirs(f"{prefix}/figures/mse_mean", exist_ok=True)
    mse.to_csv(f"{prefix}/figures/mse_mean/{selector_label}.csv", index=False)

g = sns.FacetGrid(mse, col="tree", row="selection", hue="method")
g.map(sns.scatterplot, "migration_rate", "mse", alpha=0.7)
g.add_legend()
# plt.xlim(0, 0.01)
# plt.ylim(0, 0.0001)


# %%
def selection_pairwise_difference(group):
    neutral = group[group["selection"] == "neutral"]["mean"].values
    selection = group[group["selection"] == "lognormal"]["mean"].values
    return selection - neutral


pairwise_differences = (
    (
        data.groupby(list(set(experiment_group) - {"selection"}))
        .apply(selection_pairwise_difference, include_groups=False)
        .rename("difference")
        .explode()
        .explode()
    )
    .to_frame()
    .reset_index()
)
pairwise_differences = pairwise_differences.query("migration_rate <= 1e-2")

pairwise_differences["scenario"] = (
    pairwise_differences[
        [
            "population_size",
            "mutation_rate",
            "migration_rate",
        ]
    ]
    .map(str)
    .agg("\n".join, axis=1)
)
pairwise_differences["method"] = (
    pairwise_differences[["method", "tree"]]
    .map(str)
    .agg(" ".join, axis=1)
    .str.replace("ancestry", "/ fixed tree")
    .str.replace("random", "/ sampled tree")
    .str.replace("lemey", "DTA")
    .str.replace("mascot", "MASCOT")
)

pairwise_differences.drop(["population_size", "mutation_rate", "tree"], axis=1)

fig = plt.figure(figsize=(12, 36))
sns.set_style("whitegrid")
sns.set_context("paper")
sns.boxplot(
    data=pairwise_differences,
    x="difference",
    y="scenario",
    hue="method",
    palette="pastel",
)
# Add in points to show each observation
sns.stripplot(
    pairwise_differences,
    x="difference",
    y="scenario",
    hue="method",
    palette="pastel",
    dodge=True,
    size=2,
    color=".1",
)


# %% [markdown]
## Test significance between selection scenarios


# %%
# Compute the Wilcoxon signed-rank test for pairwise differences of estimates between selection scenarios
def color_pvalues(val):
    """
    Colors elements in a DataFrame
    red if p-value < 0.01,
    orange if 0.01 <= p-value < 0.05,
    green if p-value >= 0.05.
    """
    val = abs(val)
    bgcolor = "white"
    if not np.isnan(val):
        bgcolor = "lightgreen"
    if val < 0.005:
        bgcolor = "gold"
    if val < 0.001:
        bgcolor = "tomato"

    color = "black" if val < 0.005 else "transparent"
    return f"background-color: {bgcolor}; color: {color};"


def selection_wilcoxon(group):
    neutral = group[group["selection"] == "neutral"]["mean"].values
    selection = group[group["selection"] == "lognormal"]["mean"].values
    min_length = min(len(neutral), len(selection))
    distance = selection[:min_length] - neutral[:min_length]
    pval = sp.stats.wilcoxon(distance, alternative="greater").pvalue
    if np.isnan(pval):
        print(
            "NaN:",
            group[experiment_group],
        )
    return pval


selection_pvals = (
    (
        data.groupby(list(set(experiment_group) - {"selection"}))[
            list(set(experiment_group) | {"mean", "selection"})
        ]
        .apply(selection_wilcoxon)
        .rename("pval")
    )
    .to_frame()
    .reset_index()
)

os.makedirs(f"{prefix}/figures/pvals_selection", exist_ok=True)
selection_pvals.to_csv(
    f"{prefix}/figures/pvals_selection/{selector_label}.csv", index=False
)

# %%
data.groupby(["selection", "tree", "migration_rate"])[
    "mean"
].mean().reset_index().set_index(["tree", "migration_rate"]).groupby("selection").diff()


# %%
# Compute the Wilcoxon signed-rank test for paired samples from distributions


def load_compute_pair(path1, path2, n_samples=1000, alternative="greater"):
    """Is the difference between the two distributions significant?

    By default the test is one-sided, testing if the distribution of the second sample
    is greater than the first.
    """
    trace1 = read_beast_log(path1)
    trace2 = read_beast_log(path2)

    trace1_n_samples = trace1.shape[0]
    if trace1_n_samples < n_samples:
        print(
            f"Warning: {path1} has less than {n_samples} samples: {trace1_n_samples}."
        )

    trace2_n_samples = trace2.shape[0]
    if trace2_n_samples < n_samples:
        print(
            f"Warning: {path2} has less than {n_samples} samples: {trace2_n_samples}."
        )

    model_name1 = "/".join(path1.split("/")[2:4])
    model_name2 = "/".join(path2.split("/")[2:4])

    n_samples = min(trace1_n_samples, trace2_n_samples, n_samples)

    sample1 = trace1[MIGRATION_FIELDS[model_name1]].sample(n_samples).values
    sample2 = trace2[MIGRATION_FIELDS[model_name2]].sample(n_samples).values

    diff = sample2 - sample1
    pval = sp.stats.wilcoxon(diff, alternative=alternative).pvalue

    return pval


def wilcoxon_distributed_samples(group):
    """Split group in neutral and selection scenarios and compute the Wilcoxon
    signed-rank test for all pairs of runs. Return the proportion of significant
    differences.
    """
    neutral = group[group["selection"] == "neutral"]["path"].values
    selection = group[group["selection"] == "lognormal"]["path"].values

    if len(neutral) != len(selection):
        print("Warning: different number of runs for.")

    significance = 0.05

    pvals = np.array(
        [
            load_compute_pair(neutral_path, selection_path)
            for neutral_path, selection_path in zip(selection, neutral)
        ]
    )

    return np.sum(np.array(pvals) < significance) / len(pvals)


# %%
if False:
    print("Distribution Test --- Selection")
    selection_pvals_dist = data.groupby(list(set(experiment_group) - {"selection"}))[
        ["path", "selection"]
    ].progress_apply(
        wilcoxon_distributed_samples,
    )

    os.makedirs(f"{prefix}/figures/selection_gt_pvals_dist", exist_ok=True)
    selection_pvals_dist.to_csv(
        f"{prefix}/figures/selection_gt_pvals_dist/{selector_label}.csv", index=False
    )

# %%
# Compute the Wilcoxon signed-rank test for distributed samples


def load_compute(path, rate, n_samples=1000, alternative="greater"):
    """Is the difference between the two distributions significant?

    By default the test is one-sided, testing if the distribution of the second sample
    is greater than the first.
    """
    trace = read_beast_log(path)
    trace_n_samples = trace.shape[0]
    if trace_n_samples < n_samples:
        print(f"Warning: {path} has less than {n_samples} samples: {trace_n_samples}.")

    model_name = "/".join(path.split("/")[2:4])

    n_samples = min(trace_n_samples, n_samples)
    sample = trace[MIGRATION_FIELDS[model_name]].sample(n_samples).values

    diff = sample - rate
    pval = sp.stats.wilcoxon(diff, alternative=alternative).pvalue

    return pval


def wilcoxon_distributed_samples_inference(group):
    """Split group in neutral and selection scenarios and compute the Wilcoxon
    signed-rank test for all pairs of runs. Return the proportion of significant
    differences.
    """
    entries = group["path"].values
    migration_rate = group["migration_rate"].values[0]

    significance = 0.05

    pvals = np.array([load_compute(path, migration_rate) for path in entries])

    return np.sum(np.array(pvals) < significance) / len(pvals)


# %%
if False:
    print("Distribution Test --- Inference")
    pvals_dist = data.groupby(experiment_group)[
        ["path", "migration_rate"]
    ].progress_apply(
        wilcoxon_distributed_samples_inference,
    )

    os.makedirs(f"{prefix}/figures/{selector_label}_pvals_dist", exist_ok=True)
    pvals_dist.to_csv(f"{prefix}/figures/pvals_dist/{selector_label}.csv", index=False)

# %%
data["method"] = data.method.map({"lemey": "DTA", "mascot": "MASCOT"})
sns.set_theme(font_scale=1.5, style="whitegrid")
grid = sns.FacetGrid(
    data.query("migration_rate <= 1e-2"),
    row="method",
    row_order=["DTA", "MASCOT"],
    col="tree",
    col_order=trees,
    height=6,
    aspect=1.5,
    margin_titles=True,
)
grid.map_dataframe(
    sns.violinplot,
    x="migration_rate",
    y="overestimation",
    hue="selection",
    palette=sns.color_palette(n_colors=2),
    alpha=0.5,
    split=True,
    gap=0.2,
    cut=0,
    inner="quart",
)
grid.map_dataframe(
    sns.stripplot,
    x="migration_rate",
    y="overestimation",
    hue="selection",
    palette=sns.color_palette(n_colors=2),
    alpha=0.65,
    dodge=True,
)
grid.add_legend()
grid.set(
    ylim=(0, 1),
    xlabel="migration rate",
    ylabel="posterior overestimation probability",
)
grid.set_titles(row_template="{row_name}")

grid_axes = grid.axes
grid_axes[0, 0].set_title("based on pylogenetic reconstruction")
grid_axes[0, 1].set_title("based on true genealogy")

os.makedirs(f"{prefix}/figures/posterior_overestimation_probability", exist_ok=True)
plt.savefig(
    f"{prefix}/figures/posterior_overestimation_probability/{selector_label}.pdf",
    bbox_inches="tight",
)

# %% [markdown]
## Gather statistics for posterior overestimation probability

# %%
overestimation_stats = (
    data.groupby(experiment_group)[["overestimation"]].agg("mean").reset_index()
)
overestimation_stats

# %%
os.makedirs(f"{prefix}/figures/overestimation_stats", exist_ok=True)
overestimation_stats.to_csv(
    f"{prefix}/figures/overestimation_stats/{selector_label}.csv", index=False
)


# %%
def load_and_sample_pair(path1, path2, n_samples=1000):
    """Is the difference between the two distributions significant?

    By default the test is one-sided, testing if the distribution of the second sample
    is greater than the first.
    """
    trace1 = read_beast_log(path1)
    trace2 = read_beast_log(path2)

    trace1_n_samples = trace1.shape[0]
    if trace1_n_samples < n_samples:
        print(
            f"Warning: {path1} has less than {n_samples} samples: {trace1_n_samples}."
        )

    trace2_n_samples = trace2.shape[0]
    if trace2_n_samples < n_samples:
        print(
            f"Warning: {path2} has less than {n_samples} samples: {trace2_n_samples}."
        )

    model_name1 = "/".join(path1.split("/")[2:4])
    model_name2 = "/".join(path2.split("/")[2:4])

    n_samples = min(trace1_n_samples, trace2_n_samples, n_samples)

    sample1 = trace1[MIGRATION_FIELDS[model_name1]].sample(n_samples).values
    sample2 = trace2[MIGRATION_FIELDS[model_name2]].sample(n_samples).values

    return sample1, sample2


def compute_overestimation_selection_stats(group):
    neutral = group[group["selection"] == "neutral"]["path"].values
    selection = group[group["selection"] == "lognormal"]["path"].values

    overestimations = []
    for n, s in zip(neutral, selection):
        neutral_sample, selection_sample = load_and_sample_pair(n, s)
        overestimation = (selection_sample > neutral_sample).mean()
        overestimations.append(overestimation)

    return sum(overestimations) / len(overestimations)


overestimation_selection_stats = (
    data.groupby(list(set(experiment_group) - {"selection"}))[["path", "selection"]]
    .progress_apply(compute_overestimation_selection_stats)
    .rename("overestimation")
    .reset_index()
)

# %%
os.makedirs(f"{prefix}/figures/overestimation_selection_stats", exist_ok=True)
overestimation_selection_stats.to_csv(
    f"{prefix}/figures/overestimation_selection_stats/{selector_label}.csv", index=False
)

# %%
overestimation_selection_stats

average_overestimation_selection_stats = overestimation_selection_stats.groupby(
    list(set(experiment_group) - {"selection", "migration_rate"})
).mean()

average_overestimation_selection_stats.to_csv(
    f"{prefix}/figures/overestimation_selection_stats/avg_{selector_label}.csv",
    index=False,
)

# %%
average_overestimation_stats = overestimation_stats.groupby(
    list(set(experiment_group) - {"migration_rate"})
).mean()
average_overestimation_stats.to_csv(
    f"{prefix}/figures/overestimation_stats/avg_{selector_label}.csv", index=False
)

# %%
average_overestimation_stats.reset_index().set_index(
    ["method", "selection", "tree"]
).sort_index()

# %%
average_overestimation_selection_stats.reset_index().set_index(["method", "tree"])
