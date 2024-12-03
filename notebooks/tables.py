# %%
import argparse
import os
from glob import glob

import dataframe_image as dfi
import numpy as np
import pandas as pd
from load_data import is_notebook


# %%
def table_styles(width):
    return [
        {
            "selector": "",
            "props": [("border", "1px solid grey !important")],
        },
        {
            "selector": "tbody td",
            "props": [
                ("border", "1px solid grey !important"),
                ("text-align", "center"),
            ],
        },
        {
            "selector": "th",
            "props": [
                ("border", "1px solid grey !important"),
                ("min-width", f"{width}px"),
            ],
        },
    ]


# %%
parser = argparse.ArgumentParser()

parser.add_argument(
    "--glob-type", type=str, help="Type of result file to glob, e.g. relaxed_clock"
)

if is_notebook():
    args = parser.parse_args(
        args=[
            "--glob-type",
            "relaxed_clock",
        ]
    )
    prefix = ".."
else:
    args = parser.parse_args()
    prefix = "."


glob_type = args.glob_type
ancestry_name = "true"
sampled_name = "recon"


# %%
def load_data(glob_type, data_type):
    files = glob(f"{prefix}/figures/{data_type}/{glob_type}_*.csv")
    data = pd.concat([pd.read_csv(f, dtype={"mutation_rate": str}) for f in files])
    data = data[data.migration_rate <= 1e-2]
    data["migration_rate"] = data.migration_rate.apply(lambda x: f"{x:.0e}")
    data["tree"] = data.tree.apply(
        lambda x: ancestry_name if x == "ancestry" else sampled_name
    )
    return data


# %%
pvals = load_data(glob_type, "pvals_mean")
pvals_selection = load_data(glob_type, "pvals_selection")
mse = load_data(glob_type, "mse_mean")
overestimation = load_data(glob_type, "overestimation_stats")
overestimation_selection = load_data(glob_type, "overestimation_selection_stats")


# %%
# Styling function to color cells based on p-value
def color_pvalues(val):
    """
    Colors elements in a DataFrame
    red if p-value < 0.001,
    orange if 0.001 <= p-value < 0.005,
    green if p-value >= 0.005.
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


def style_pvals(pvals, index):
    # Reshape the Dataframe to have p-values in grid form
    pvals_grid = pvals.pivot(
        index=index,
        columns=["selection", "migration_rate"],
        values="pval",
    )
    pvals_grid.columns.names = ["selection", "migration rate"]
    pvals_grid.index.rename([x.replace("_", " ") for x in index], inplace=True)
    pvals_grid.sort_index(inplace=True, ascending=[True, False, True])

    # Apply the styling
    pvals_styled = (
        pvals_grid.style.map(color_pvalues)
        .format(lambda x: "-" if x < 0 else "+")
        .set_table_styles(table_styles(width=40))
    )
    return pvals_styled


# %%
pvals.reset_index(inplace=True, drop=True)
pvals_styled = style_pvals(
    pvals, index=["population_size", "mutation_rate", "method", "tree"]
)
dfi.export(pvals_styled, f"{prefix}/figures/pvals_mean/{glob_type}.png", max_cols=-1)
pvals_styled

# %%
pvals_tiny = pvals[pvals.population_size == 100].copy()
pvals_tiny.drop(columns=["population_size"], inplace=True)
pvals_tiny.reset_index(inplace=True, drop=True)

pvals_styled = style_pvals(pvals_tiny, index=["mutation_rate", "method", "tree"])
dfi.export(
    pvals_styled, f"{prefix}/figures/pvals_mean/{glob_type}_tiny.png", max_cols=-1
)
pvals_styled

# %%
pvals_small = pvals[pvals.population_size == 1000].copy()
pvals_small.drop(columns=["population_size"], inplace=True)
pvals_small.reset_index(inplace=True, drop=True)

pvals_styled = style_pvals(pvals_small, index=["mutation_rate", "method", "tree"])
dfi.export(
    pvals_styled, f"{prefix}/figures/pvals_mean/{glob_type}_small.png", max_cols=-1
)
pvals_styled

# %%
mse_width = 55

mse.reset_index(inplace=True, drop=True)
mse["migration_rate_value"] = mse.migration_rate.map(lambda x: float(x))
mse.sort_values(
    ["selection", "migration_rate_value"],
    inplace=True,
    ascending=[False, True],
)
mse.drop("migration_rate_value", axis=1, inplace=True)
mse_grid = mse.pivot(
    index=["population_size", "mutation_rate", "method", "tree"],
    columns=["selection", "migration_rate"],
    values="mse",
)
mse_grid.index.rename(
    ["population size", "mutation rate", "method", "tree"], inplace=True
)
mse_grid.sort_index(inplace=True, ascending=[True, False, True])
mse_grid_styled = (
    mse_grid.style.format("{:.2e}")
    .set_table_styles(table_styles(width=mse_width))
    .background_gradient(
        cmap="YlGnBu",
        axis=None,
        vmax=-4,
        gmap=mse_grid.map(np.log10),
    )
)
dfi.export(mse_grid_styled, f"{prefix}/figures/mse_mean/{glob_type}.png", max_cols=-1)
mse_grid_styled

# %%
mse_tiny = mse[mse.population_size == 100].copy()
mse_tiny.drop(columns=["population_size"], inplace=True)
mse_tiny.reset_index(inplace=True, drop=True)

mse_grid = mse_tiny.pivot(
    index=["mutation_rate", "method", "tree"],
    columns=["selection", "migration_rate"],
    values="mse",
)
mse_grid.sort_index(inplace=True, ascending=[True, False, True])
mse_grid_styled = (
    mse_grid.style.format("{:.2e}")
    .set_table_styles(table_styles(width=mse_width))
    .background_gradient(
        cmap="YlGnBu",
        axis=None,
        vmax=-4,
        gmap=mse_grid.map(np.log10),
    )
)
os.makedirs(f"{prefix}/figures/mse", exist_ok=True)
dfi.export(
    mse_grid_styled, f"{prefix}/figures/mse_mean/{glob_type}_tiny.png", max_cols=-1
)
mse_grid_styled

# %%
mse_small = mse[mse.population_size == 1000].copy()
mse_small.drop(columns=["population_size"], inplace=True)
mse_small.reset_index(inplace=True, drop=True)

mse_grid = mse_small.pivot(
    index=["mutation_rate", "method", "tree"],
    columns=["selection", "migration_rate"],
    values="mse",
)
mse_grid.sort_index(inplace=True, ascending=[True, False, True])
mse_grid_styled = (
    mse_grid.style.format("{:.2e}")
    .set_table_styles(table_styles(width=mse_width))
    .background_gradient(
        cmap="YlGnBu",
        axis=None,
        vmax=-4,
        gmap=mse_grid.map(np.log10),
    )
)
dfi.export(
    mse_grid_styled, f"{prefix}/figures/mse_mean/{glob_type}_small.png", max_cols=-1
)
mse_grid_styled

# %%


# Styling function to color cells based on p-value
def color_pvalues_transparent(val):
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

    color = "transparent"
    return f"background-color: {bgcolor}; color: {color};"


pvals_selection["migration_rate"] = pvals_selection.migration_rate.map(
    lambda x: float(x)
)
pvals_selection_table = pvals_selection.pivot(
    index=["population_size", "mutation_rate", "method", "tree"],
    columns=["migration_rate"],
    values="pval",
)

pvals_selection_table.columns = pvals_selection_table.columns.map(lambda x: f"{x:.0e}")
pvals_selection_table.columns.name = "migration rate"
pvals_selection_table.index.rename(
    ["population size", "mutation rate", "method", "tree"], inplace=True
)

pvals_selection_table_styled = (
    pvals_selection_table.style.format("{:.2e}")
    .map(color_pvalues_transparent)
    .set_table_styles(table_styles(width=40))
)

dfi.export(
    pvals_selection_table_styled,
    f"{prefix}/figures/pvals_selection/{glob_type}.png",
    max_cols=-1,
    dpi=200,
)
pvals_selection_table_styled

# %%
overestimation

# %%
overestimation_width = 38
overestimation.reset_index(inplace=True, drop=True)
overestimation["migration_rate_value"] = overestimation.migration_rate.map(
    lambda x: float(x)
)
overestimation.sort_values(
    ["selection", "migration_rate_value"],
    inplace=True,
    ascending=[False, True],
)
overestimation.drop("migration_rate_value", axis=1, inplace=True)
overestimation_grid = overestimation.pivot(
    index=["population_size", "mutation_rate", "method", "tree"],
    columns=["selection", "migration_rate"],
    values="overestimation",
)
overestimation_grid.index.rename(
    ["population size", "mutation rate", "method", "tree"], inplace=True
)
overestimation_grid.sort_index(inplace=True, ascending=[True, False, True])
overestimation_grid_styled = (
    overestimation_grid.style.format("{:.2f}")
    .set_table_styles(table_styles(width=overestimation_width))
    .background_gradient(
        cmap="coolwarm",
        axis=None,
        vmin=0,
        vmax=1,
    )
)
dfi.export(
    overestimation_grid_styled,
    f"{prefix}/figures/overestimation_stats/{glob_type}.png",
    max_cols=-1,
    dpi=200,
)
overestimation_grid_styled

# %%
overestimation_width = 38
overestimation_selection.reset_index(inplace=True, drop=True)
overestimation_selection["migration_rate"] = (
    overestimation_selection.migration_rate.map(lambda x: float(x))
)
overestimation_selection_grid = overestimation_selection.pivot(
    index=["population_size", "mutation_rate", "method", "tree"],
    columns=["migration_rate"],
    values="overestimation",
)
overestimation_selection_grid.index.rename(
    ["population size", "mutation rate", "method", "tree"], inplace=True
)
overestimation_selection_grid.columns = overestimation_selection_grid.columns.map(
    lambda x: f"{x:.0e}"
)
overestimation_selection_grid.columns.name = "migration rate"
overestimation_selection_grid.sort_index(inplace=True, ascending=[True, False, True])
overestimation_selection_grid_styled = (
    overestimation_selection_grid.style.format("{:.2f}")
    .set_table_styles(table_styles(width=overestimation_width))
    .background_gradient(
        cmap="coolwarm",
        axis=None,
        vmin=0,
        vmax=1,
    )
)
dfi.export(
    overestimation_selection_grid_styled,
    f"{prefix}/figures/overestimation_selection_stats/{glob_type}.png",
    max_cols=-1,
    dpi=200,
)
o
