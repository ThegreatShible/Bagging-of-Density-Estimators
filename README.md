# Bagging-of-Density-Estimators

Implementation and optimisation of the model introduced in [Bagging of Density Estimators](https://hal.archives-ouvertes.fr/hal-01856183v2/document).

We obtained an optimisation of x32 in some functions by : 
* Implement critical functions in Rcpp
* Vectorize some computations
* Code level optimization 
* Parallelisation via `parallel` library
