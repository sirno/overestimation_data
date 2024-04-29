# HIV-1 Simulation Reveals Overestimation Bias in Within-Host Phylodynamic Migration Rate Estimates Under Selection

This repository contains the code and when possible data used in the manuscript
"HIV-1 Simulation Reveals Overestimation Bias in Within-Host Phylodynamic
Migration Rate Estimates Under Selection".

## Abstract

Phylodynamic methods are widely used to analyze the spread of viruses between
and within hosts. For HIV-1, these methods have been used to estimate migration
rates between different anatomical compartments within a host. However, these
methods often assume that there is no selective pressure acting on the virus.
Even though it is known that selection will affect the spread of viruses within
a host, the impact of selection on the estimation of migration rates has not
been tested. In this study, we used a novel agent-based simulation tool,
`virolution`, to simulate the evolution of HIV-1 within different anatomical
compartments of a host. We generated viral sequences and genealogies under two
evolutionary scenarios: neutral evolution and selection-driven evolution.
Through the application of Discrete Trait Analysis (DTA) and the MASCOT model
in the Bayesian inference framework BEAST2, we identified significant biases in
migration rate estimates induced by selection. Our results reveal that commonly
used phylogeographic methods, which assume neutral selection, can significantly
overestimate migration rates when selection is present. This study underscores
the need for rigorous testing of phylodynamic models on data sets with
realistic assumptions before their application.

## Data

The raw data used in the manuscript constitutes to about 600 GB of log files
from the BEAST2 analysis. We therefore provide access to the template XML
files, the scripts used to generate the XML files, examples of BEAST2 XML and
log files, and the scripts used to analyze the log files.

Additionally, we provide the simulated sequences and genealogies used in the
manuscript.

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

Install dependencies:

```bash
pip install -r requirements.txt
```

### Data Analysis

If you have [just](https://github.com/caseyjust/just) installed, you can run
the analysis with:

```bash
just create_figures
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
just list_xml
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

