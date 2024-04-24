# %%
from glob import glob

import argparse
import dataframe_image as dfi
import numpy as np
import pandas as pd

from load_data import is_notebook


# %%
parser = argparse.ArgumentParser()
parser.add_argument(
    "--selector",
    type=str,
    help="The selector to use",
)

if is_notebook():
    args = parser.parse_args(
        args=[
            "--selector",
            "relaxed_clock",
        ]
    )
    prefix = ".."
else:
    args = parser.parse_args()
    prefix = "."

# %%
glob_type = args.selector
ancestry_name = "true"
sampled_name = "recon"

# %%
# load two-sided p-values
pval_files = glob(f"{prefix}/figures/{glob_type}_*_pvals.csv")
pval_files = filter(lambda x: "selection" not in x, pval_files)
pvals = pd.concat([pd.read_csv(f, dtype={"mutation_rate": str}) for f in pval_files])

pvals = pvals[pvals.migration_rate <= 1e-2]
pvals["migration_rate"] = pvals.migration_rate.apply(lambda x: f"{x:.0e}")
pvals["tree"] = pvals.tree.apply(
    lambda x: ancestry_name if x == "ancestry" else sampled_name
)

# %%
# load mean squared errors
mse_files = glob(f"{prefix}/figures/{glob_type}_*_mse.csv")
mse = pd.concat([pd.read_csv(f, dtype={"mutation_rate": str}) for f in mse_files])

mse = mse[mse.migration_rate <= 1e-2]
mse["migration_rate"] = mse.migration_rate.apply(lambda x: f"{x:.0e}")
mse["tree"] = mse.tree.apply(
    lambda x: ancestry_name if x == "ancestry" else sampled_name
)

# %%
# load p-values of one-sided tests (with vs without selection)
pvals_selection_files = glob(f"{prefix}/figures/{glob_type}_*_selection_pvals_gt.csv")
pvals_selection = pd.concat(
    [pd.read_csv(f, dtype={"mutation_rate": str}) for f in pvals_selection_files]
)
pvals_selection = pvals_selection[pvals_selection.migration_rate <= 1e-2]
pvals_selection["tree"] = pvals_selection.tree.apply(
    lambda x: ancestry_name if x == "ancestry" else sampled_name
)


# %%
# styling table styles specifically needed for export
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


# styling function to color cells based on p-value
def color_pvalues_transparent(val):
    """
    Colors elements in a DataFrame
    red if p-value < 0.001,
    orange if 0.01 <= p-value < 0.005,
    green if p-value >= 0.005.
    color: transparent.
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


# styling function to color cells based on p-value
def color_pvalues(val):
    """
    Colors elements in a DataFrame
    red if p-value < 0.001,
    orange if 0.001 <= p-value < 0.005,
    green if p-value >= 0.005.
    color: black if p-value < 0.005 else transparent.
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


# %%
def style_pvals(pvals, index):
    # Reshape the Dataframe to have p-values in grid form
    pvals_grid = pvals.pivot(
        index=index,
        columns=["selection", "migration_rate"],
        values="pval",
    )
    pvals_grid.columns.names = ["selection", "migration rate"]
    pvals_grid.index.rename([x.replace("_", " ") for x in index], inplace=True)

    # Apply the styling
    pvals_styled = (
        pvals_grid.style.map(color_pvalues)
        .format(lambda x: "-" if x < 0 else "+")
        .set_table_styles(table_styles(width=40))
    )
    return pvals_styled


# %%
# Export plots for p-values of two-sided tests
pvals.reset_index(inplace=True, drop=True)
pvals_styled = style_pvals(
    pvals, index=["population_size", "mutation_rate", "method", "tree"]
)
dfi.export(pvals_styled, f"{prefix}/figures/pvals_{glob_type}.png", max_cols=-1)
pvals_styled

# %%
pvals_tiny = pvals[pvals.population_size == 100].copy()
pvals_tiny.drop(columns=["population_size"], inplace=True)
pvals_tiny.reset_index(inplace=True, drop=True)

pvals_styled = style_pvals(pvals_tiny, index=["mutation_rate", "method", "tree"])
dfi.export(pvals_styled, f"{prefix}/figures/pvals_{glob_type}_tiny.png", max_cols=-1)
pvals_styled

# %%
pvals_small = pvals[pvals.population_size == 1000].copy()
pvals_small.drop(columns=["population_size"], inplace=True)
pvals_small.reset_index(inplace=True, drop=True)

pvals_styled = style_pvals(pvals_small, index=["mutation_rate", "method", "tree"])
dfi.export(pvals_styled, f"{prefix}/figures/pvals_{glob_type}_small.png", max_cols=-1)
pvals_styled

# %%
# Export plots for mean squared errors
mse_width = 55

mse.reset_index(inplace=True, drop=True)
mse_grid = mse.pivot(
    index=["population_size", "mutation_rate", "method", "tree"],
    columns=["selection", "migration_rate"],
    values="mse",
)
mse_grid.index.rename(
    ["population size", "mutation rate", "method", "tree"], inplace=True
)
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
dfi.export(mse_grid_styled, f"{prefix}/figures/mse_{glob_type}.png", max_cols=-1)
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
dfi.export(mse_grid_styled, f"{prefix}/figures/mse_{glob_type}_tiny.png", max_cols=-1)
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
dfi.export(mse_grid_styled, f"{prefix}/figures/mse_{glob_type}_small.png", max_cols=-1)
mse_grid_styled


# %%
# Export plots for p-values of one-sided tests (with vs without selection)
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
    f"{prefix}/figures/pvals_selection_{glob_type}.png",
    max_cols=-1,
    dpi=200,
)
pvals_selection_table_styled
