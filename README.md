# Estimation of the transmission dynamics of H5N1 HPAI outbreak in a dairy herd using a modeling approach

## ABSTRACT

The emergence of Highly Pathogenic Avian Influenza (HPAI 2.3.4.4b) in dairy herds in 2024 across 19 states in the United States of America has raised concerns regarding the potential national and global zoonotic impact. All recent modeling efforts implemented homogenous cattle-to-cattle (both intra and inter-herd) transmission models which did not capture the real-world heterogeneity in mixing of animals, individual variations in susceptibility and infectiousness and clinical incidences across pens and lactation groups. The aim of this study was to develop a heterogenous transmission model to estimate the epidemiological parameters for intra-herd HPAI transmission on Californian dairies. We developed a validated stochastic agent-based model to estimate the epidemiological parameters for intra-herd HPAI transmission on California dairies. The hierarchical agent-based model also parameterized stochastic cattle movements within-herd to simulate real dairy management practices. A cow-level SEIR transmission approach was assumed during the outbreak. A novel Bayesian Optimizer with Gaussian Process (BO-GP) was fitted to the agent-based model for validation which converged within 25-40 iterations (out of 100 per farm) with minimal loss over two distinct error metrics, namely, Poisson loss function and temporal distance metric. Our optimized simulations estimated an average R0 was around 10.7-10.8 across all farms within the first 15 days of observed outbreak on four dairy farms with a mean effective transmission rate of 4.5% per contact between susceptible and infectious cows within each pen. Our model demonstrated that the movement of cows between pens ensured localized clusters of outbreaks within the sub-herds (pen population) that prolonged the overall outbreak within farms. We estimated the total duration of infection between 14.5 and 28 days, which is higher than the estimates from the homogenous models. With an integrated hierarchical agent-based model combined with Bayesian approximation, we produced actionable insights on the epidemiology of intra-farm spread of HPAI within cow herds, thereby guiding both future model development and applied disease control strategy.

## Keywords

Avian Influenza, dairy, agent-based modeling, transmission dynamics


## Project Organization
------------

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── Inputs         <- Data from third-party sources that is shareable.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── Outputs        <- Outputs of the modeling.
    ├── notebooks          <- Jupyter notebooks. The naming convention is a date (for ordering),
    │                         the creator's initials, and a short `_` delimited description, e.g.
    │                         `20260101_psk_data-exploration.ipynb`.
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make predictions
    │   │
    │   └── visualization  <- Scripts to create exploratory and results-oriented visualizations
--------------
