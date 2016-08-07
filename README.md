# BBMV
This repository provides a set of R functions to fit general macroevolutionary models for continuous traits evolving in adaptive landscapes of any shape.
This new model is based on [bounded Brownian motion](https://github.com/fcboucher/BBM) (BBM), in which a continuous trait undergoes constant-rate diffusion between two reflective bounds. In addition to this random component, the trait evolves in a potential and is thus subject to a force that pulls it towards specific values - this force can be of any shape. You can see an example of the behaviour of the model [here](https://github.com/fcboucher/BBMV/blob/master/BBM%2BV%20basics%20Figure.png). We label this model *BBM+V*, for BBM + potential.

Functions of the **BBMV** package only depend on the **ape** package in R and likelihoods are compatible with those of other models fitted by the *fitContinuous* function in package **geiger**. The package implements both maximum likelihood and MCMC estimation of model parameters.

Functions were written by [Florian Boucher](https://sites.google.com/site/floriaboucher/) based on equations from [Vincent Démery](https://www.pct.espci.fr/~vdemery/).

The [R folder](https://github.com/fcboucher/BBMV/tree/master/R) contains functions to simulate traits evolving under BBM+V, plot adaptive landscapes, and fit the model to empirical data using either maximum-likelihood or MCMC estimation. 
The [tutorial](https://github.com/fcboucher/BBMV/blob/master/Tutorial-BBMV.md) shows basic examples of use of the functions to simulate and infer under BBM+V and in case you need it there is also a full [R script of the tutorial](https://github.com/fcboucher/BBMV/blob/master/Example_ML_MCMC.R). 

Help files for each function  in the package can be found in the [manual](https://github.com/fcboucher/BBMV/blob/master/BBMV-manual.pdf). In addition, function names should be rather self-explanatory and R scripts are heavily commented: they should give you all the information needed on parameters, outputs, etc.

The paper describing the model is still work in progress...
