# BBMV
This repository provides a set of R functions to fit general macroevolutionary models for continuous traits evolving in macroevolutionary landscapes of any shape. As such, the model can accomodate a large variety of evolutionary scenarios, including hard constraints on trait evolution and selection of any possible kind.

The latest version of the **BBMV** package can also be found on [CRAN](https://CRAN.R-project.org/package=BBMV). 

This new model is based on [bounded Brownian motion](https://github.com/fcboucher/BBM) (BBM), in which a continuous trait undergoes constant-rate diffusion between two reflective bounds. In addition to this random component, the trait evolves in a potential and is thus subject to a force that pulls it towards specific values - this force can be of any shape. You can see an example of the behaviour of the model [here](https://github.com/fcboucher/BBMV/blob/master/FPK_basics_figure.pdf). We label this model *FPK* since it is based on the Fokker-Planck equation, also known as the Kolmogorov forward equation in population genetics. The *FPK* model is a generalization of other classic models for continuous trait evolution, namely Brownian Motion and the Ornstein-Uhlenbeck model.  The *FPK* model has a special case in which hard, reflective, bounds exist on the trait interval, which we label *BBMV* for BBM + potential.

Functions of the **BBMV** package only depend on the **ape** package in R and likelihoods are compatible with those of other models fitted by the *fitContinuous* function in package **geiger**. The package implements both maximum likelihood and MCMC estimation of model parameters.

Functions were written by [Florian Boucher](https://sites.google.com/site/floriaboucher/) based on equations from [Vincent Démery](https://www.pct.espci.fr/~vdemery/).

The [R folder](https://github.com/fcboucher/BBMV/tree/master/R) contains functions to simulate traits evolving under the *FPK* (or *BBMV*) model, plot macroevolutionary landscapes, and fit the model to empirical data using either maximum-likelihood or MCMC estimation. 
The [tutorial](https://github.com/fcboucher/BBMV/blob/master/Tutorial-BBMV.md) shows basic examples of use of the functions to simulate and infer under *FPK*. In case you need it there is also a full [R script of the tutorial](https://github.com/fcboucher/BBMV/blob/master/Example_script_BBMV_package.r). 

The paper describing the model has been accepted in *Systematic Biology* and is currently [in press](https://academic.oup.com/sysbio/article-abstract/doi/10.1093/sysbio/syx075/4210009/A-General-Model-for-Estimating-Macroevolutionary?redirectedFrom=fulltext), but you can find a preprint of the latest version [here](https://github.com/fcboucher/BBMV/blob/master/Manuscript/). 

Help files for all functions can be found in this [pdf](https://github.com/fcboucher/BBMV/blob/master/BBMV-manual.pdf) and there is also a troubleshooting section at the end of the [tutorial](https://github.com/fcboucher/BBMV/blob/master/Tutorial-BBMV.md).

### Latest additions

- 2017.10.23: it is now possible to fit the *FPK* model on multiple clades together in order to statistically test whether they share a similar macroevolutionary landscape. Examples of use of these new functions are given in the [tutorial](https://github.com/fcboucher/BBMV/blob/master/Tutorial-BBMV.md) and the [example R script](https://github.com/fcboucher/BBMV/blob/master/Example_script_BBMV_package.r). Be careful when using these functions since optimizing the likelihood is often difficult.

- 2017.09.29: all functions (ML and MCMC) can now take measurement error in trait data at the tips of the tree into account. The format required for the trait data in this case is explained in the [tutorial](https://github.com/fcboucher/BBMV/blob/master/Tutorial-BBMV.md) and the [example R script](https://github.com/fcboucher/BBMV/blob/master/Example_script_BBMV_package.r). Be cautious when fitting the model with measurement error since the behaviour of the model has not been extensively tested yet.
