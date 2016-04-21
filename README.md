# BBMV
This repository provides a set of R functions to fit general macroevolutionary models for continuous traits evolving in adaptive landscapes of any shape.
This new model is based on [bounded Brownian motion](https://github.com/fcboucher/BBM) (BBM), in which a continuous trait undergoes constant-rate diffusion between two reflective bounds. In addition to this random component, the trait evolves in a potential and is thus subject to a force that pulls it towards specific value - this force can be of any shape. We label this model BBM+V, for BBM + potential.

Functions of the BBMV package only depend on the {ape} package in R and likelihoods are compatible with those of other models fitted by the 'fitContinuous' function in package {geiger}. The package implements both maximum likelihood and MCMC estimation of model parameters.

Functions were written by [Florian Boucher](https://sites.google.com/site/floriaboucher/) based on equations from [Vincent Démery](https://www.pct.espci.fr/~vdemery/).

Final versions of R functions as well as the paper describing the model are still work in progress...
