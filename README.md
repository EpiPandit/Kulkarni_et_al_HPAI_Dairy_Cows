# Estimation of the transmission dynamics of H5N1 HPAI outbreak in a dairy herd using a modeling approach
[![DOI](https://zenodo.org/badge/1380203248.svg)](https://doi.org/10.5281/zenodo.24292683)
## ABSTRACT

This repository contains a pen-structured, stochastic agent-based model of HPAI H5N1 (clade 2.3.4.4b) transmission within California dairy herds. Individual cows move between functional pens according to dairy management practices, and transmission follows a cow-level SEIR process within each pen. The model was calibrated to clinical cases recorded on four dairy farms (September–November 2024), using approximate Bayesian computation followed by Bayesian optimization with a Gaussian process. Two error metrics were used: a distance-based metric and a Poisson loss function.

Main estimates across the four farms:

- an average R0 of 4.3 (distance-based metric) and 3.2 (Poisson loss function), with a range of 1.2–6.3;
- a daily transmission probability of 0.0006–0.0038 between an infectious and a susceptible cow in the same pen;
- a duration of infection of 13.6–20.0 days per cow.
- Outbreak size was most sensitive to herd size and transmission rate, and on the larger farms infection persisted at low levels for months after the observed outbreak.

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
    ├── notebooks          <- Jupyter notebooks
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    ├── src                <- Source code for use in this project.
