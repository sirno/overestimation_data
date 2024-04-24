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

The raw data used in the manuscript constitutes to about 600GB of log files
that the beast analysis generated. We therefore provide access to the template
XML files, the scripts used to generate the XML files, examples of BEAST2 XML
and log files, and the scripts used to analyze the log files.

Additionally, we provide the simulated sequences and genealogies used in the
manuscript.

## Usage

### Data Analysis

Install dependencies:

```bash
pip install -r requirements.txt
```

If you have [just](https://github.com/caseyjust/just) installed, you can run
the analysis with:

```bash
just create_figures
```



