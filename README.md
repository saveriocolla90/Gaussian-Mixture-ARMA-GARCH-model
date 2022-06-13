# Gaussian Mixture ARMA-GARCH model

The intent of this project is to produce a stable and reusable module that implements a Gaussian Mixture ARMA-GARCH model 
to fit a generic time series with zero mean and gaussian distributed standard deviations.

An initial implementation using Python3 language has been performed in a semi-structured way, by adopting the Expectation
Maximization algorithm as core strategy for computing the Maximum Likelihood Estimator in an iteratively way, leading to
an initial satisfactory data fit for the autoregressive ARMA-GARCH model. 

The goal is to implement the Python core code in a much more robust and optimized way and eventually translate the main
arithmetical part in Fortran language, which is much more comfortable in doing expensive computational operations.

The idea follows the paper of reference "Finite Mixture of ARMA-GARCH Model for Stock Price Prediction" by H. Tang, K. Chiu,
and L. Xu, University of Hong Kong (2003) that can be find in the /Docs folder.
