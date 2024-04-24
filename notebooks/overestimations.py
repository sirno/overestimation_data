# %%
import random

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

try:
    from IPython.display import Markdown
except ImportError:
    Markdown = lambda x: print(x)
from tqdm import tqdm

from load_data import (
    fast_scandir,
    is_notebook,
    load_data,
    generate_hash,
    get_function_hash,
    get_migration_rates,
)

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
    fig.legend(
        handles[:-2],
        labels[:-2],
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
    labels[0] = "Regime"
    fig.legend(
        handles[1:-2],
        labels[1:-2],
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

plt.savefig(
    f"{prefix}/figures/migration_rates_{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## Lemey

# %%
plot_migration_rates("beast/lemey", data, ylim=(0.0, 0.02), xlim=(0.0, 0.01))

plt.savefig(
    f"{prefix}/figures/lemey_migration_rates_{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## Mascot

# %%
plot_migration_rates("beast/mascot", data, ylim=(0.0, 0.02), xlim=(0.0, 0.01))

plt.savefig(
    f"{prefix}/figures/mascot_migration_rates_{selector_label}.pdf", bbox_inches="tight"
)

# %% [markdown]
## BDMM

# %%
# plot_migration_rates("beast/bdmm", data, ylim=(0.0, 0.5))


# %%
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

plt.savefig(f"{prefix}/figures/posterior_{selector_label}.pdf", bbox_inches="tight")

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
    f"{prefix}/figures/lemey_posterior_{selector_label}.pdf", bbox_inches="tight"
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
    f"{prefix}/figures/mascot_posterior_{selector_label}.pdf", bbox_inches="tight"
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
sns.boxplot(data=data, x="template", y="ess", hue="selection")

# %%
sns.boxplot(data=data, x="template", y="n_samples", hue="tree")

# %%
data.loc[data["sample"].str.split("/").apply(lambda x: int(x[-1])) > 10]

# %%
data.groupby(["template", "selection", "tree", "migration_rate"])[
    ["migration_rate", "mean"]
].apply(lambda group: int(group["mean"].mean() >= group["migration_rate"].iloc[0]))

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
    pvals.to_csv(f"{prefix}/figures/{selector_label}_pvals.csv", index=False)

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
    for ax in axes:
        ax.set_xlim(-0.0005, 0.0105)
        ax.fill_between(
            [-1, 1],
            [0.05, 0.05],
            [1e-6, 1e-6],
            color="red",
            alpha=0.1,
            label="p < 0.05",
        )
        ax.set_xlabel("migration rate")
        ax.set_yscale("log")


# %%
pvals

# %%
plot_pvals(pvals, "lemey")
plt.savefig(f"{prefix}/figures/lemey_pvals_{selector_label}.pdf", bbox_inches="tight")

# %%
plot_pvals(pvals, "mascot")
plt.savefig(f"{prefix}/figures/mascot_pvals_{selector_label}.pdf", bbox_inches="tight")

# %%
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
mse = (
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
        lambda group: ((group["mean"] - group["migration_rate"]) ** 2).mean(),
    )
    .rename("mse")
    .reset_index()
)
mse = (
    mse.set_index(
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
    mse.to_csv(f"{prefix}/figures/{selector_label}_mse.csv", index=False)

g = sns.FacetGrid(mse, col="tree", row="selection", hue="method")
g.map(sns.scatterplot, "migration_rate", "mse", alpha=0.7)
g.add_legend()
# plt.xlim(0, 0.01)
# plt.ylim(0, 0.0001)


# %%
def selection_pairwise_difference(group):
    neutral = group[group["selection"] == "neutral"]["mean"].values
    selection = group[group["selection"] == "lognormal"]["mean"].values
    return np.subtract.outer(selection, neutral)


pairwise_differences = (
    (
        data.groupby(
            [
                "method",
                "mutation_rate",
                "population_size",
                "tree",
                "migration_rate",
            ]
        )
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


# %%
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
    l = min(len(neutral), len(selection))
    distance = selection[:l] - neutral[:l]
    pval = sp.stats.wilcoxon(distance, alternative="greater").pvalue
    if np.isnan(pval):
        print(
            "NaN:",
            group[
                [
                    "method",
                    "mutation_rate",
                    "population_size",
                    "migration_rate",
                    "selection",
                    "tree",
                ]
            ],
        )
    return pval


selection_pvals = (
    (
        data.groupby(
            [
                "method",
                "mutation_rate",
                "population_size",
                "tree",
                "migration_rate",
            ]
        )[
            [
                "method",
                "mutation_rate",
                "population_size",
                "tree",
                "migration_rate",
                "selection",
                "mean",
            ]
        ]
        .apply(selection_wilcoxon)
        .rename("pval")
    )
    .to_frame()
    .reset_index()
)

selection_pvals.to_csv(
    f"{prefix}/figures/{selector_label}_selection_pvals_gt.csv", index=False
)
