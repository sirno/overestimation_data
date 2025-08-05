# HIV-1 Simulation Reveals Overestimation Bias in Within-Host Phylodynamic Migration Rate Estimates Under Selection

This repository contains the code and when possible data used in the manuscript
"HIV-1 Simulation Reveals Overestimation Bias in Within-Host Phylodynamic
Migration Rate Estimates Under Selection".

## Abstract

Phylodynamic methods are widely used to infer the population dynamics of
viruses between and within hosts. For HIV-1, these methods have been used to
estimate migration rates between different anatomical compartments within a
host. These methods typically assume that the genomic regions used for
reconstruction are evolving without selective pressure, even though other parts
of the viral genome are known to experience strong selection. In this study, we
investigate how selection affects phylodynamic migration rate estimates. To
this end, we developed a novel agent-based simulation tool, `virolution`, to
simulate the evolution of virus within two anatomical compartments of a host.
Using this tool, we generated viral sequences and genealogies assuming both,
neutral evolution and purifying selection that is concordant in both
compartments. We found that, under the selection regime, migration rates are
significantly overestimated with a stochastic mixture model and a structured
coalescent model in the Bayesian inference framework BEAST2. Our results reveal
that commonly used phylogeographic methods, which assume neutral evolution, can
significantly bias migration rate estimates in selective regimes. This study
underscores the need for assessing the robustness of phylodynamic analysis with
respect to more realistic selection regimes.

## Data

The raw data used in the manuscript constitutes to about 600 GB of log files
from the BEAST2 analysis. We therefore provide access to the template XML
files, the scripts used to generate the XML files, examples of BEAST2 XML and
log files, and the scripts used to analyze the log files.

Additionally, raw and processed data is available in the research collection at
[doi:10.3929/ethz-b-000749260](https://doi.org/10.3929/ethz-b-000749260).

The simulations were conducted with `virolution` as archived at
[doi:10.5281/zenodo.15827569](https://doi.org/10.5281/zenodo.15827569).

## Overview

The repository is organized as follows:

- `data/`: Configurations and data used in the analysis.
- `figures/`: Figures generated from the analysis.
- `notebooks/`: Python notebooks to inspect the data and results.
- `out/`: Any output files generated during the analysis.
- `remote/`: Any remote scripts and wrappers used for the simulations and
  BEAST2 analysis.
- `scripts/`: Scripts used to generate the XML files and run the BEAST2 files.

## Usage

This repository can be used with [`uv`](https://github.com/astral-sh/uv).

```bash
pip install uv
```

### Data Analysis

To run the data analysis and generate all figures based on the processed data:

```bash
uv run scripts/create_figures.sh --load-store
```

This command will generate the figures in the `figures/` folder, based on the
aggregated analysis results in `out/cache/`.

### Inspect the data

Open the notebooks in the notebooks folder to inspect the data and the results.
The notebooks have been edited and executed in VS Code. They can be converted
to jupyter notebooks with `jupytext`.

### Template XML files

The template XML files are located in `data/beast/`. The XML files are
organized in subdirectories according to the method and scenario used in the
analysis. Get an overview of the XML files with:

```bash
find data -name "template.xml" | sort
```

The template XML files contain placeholders for the sequences, times, and types
data. The placeholders are formatted as `{placeholder}`. Further, any tree
operators are wrapped between a `<!--tree-operators-->` and
`<!--/tree-operators/-->` tag, such that they can be disabled when needed.

### BEAST2 Analysis

The BEAST2 analysis was performed on a high-performance computing cluster. The
analysis was performed using the BEAST2 version 2.7.4 (MASCOT version 3.0.0,
BEAST_CLASSIC version 1.6.2).

The code used to create the XML files and run the BEAST2 analysis is located in
`scripts/` and any server-side scripts are located in `remote/`. The scripts
run with the SLURM scheduler and may need to be adjusted for successful execution
on other systems.

