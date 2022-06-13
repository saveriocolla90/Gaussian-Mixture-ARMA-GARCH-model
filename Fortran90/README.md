# State of the art

The Fortran90 folder contains a partial initial implemantaion of the Gaussian Mixture ARMA-GARCH model fit algorithm.

The aim of this sub-project is to reach a full stable, reusable and optimized module that will provides a general 
model fit implementation using the Expectation Maximization MLE method, in order to substitue the core algorithm
of the main Python3 project with a much faster one. 
Eventually there will be a sort of inter-relationship between the two codes, letting the Fortran core alghoritm invoked by
the more general Python interface; the latter can handle the dataframe manipulation much better and retrive/store data 
in different type of formats, adopting the large malleable pool of different python's packages.

The state of the art is outlined by the following:
* General model parameters initialization has been achived with 3 derived types, expanding each one from the base ARMA(R,S)
    model, through the GARCH(Q,P) to reach a general Gaussian Mixture(K) of different components.
* Parameters values initialization has been set to a random state, which has been performed by 3 different class subroutines. 
* A first function to check the initial stability of the model has been implemented.
* The two main function to compute the ARMA and GARCH time series of residuals and variances (respectively) has been implemented 
    with a correspondence in the class of the same derived type. (mod_MODELS.f90)
* A TIME_SERIES object has been built with a class of two I/O subroutines. (mod_IO.f90)
* All the relevant fixed parameters values are stored in a different module (time series length, models ranks, ecc.). (mod_PARAMS.f90) 
* The main program provides a simple introduction to the general model initialization, loading the time series data from file, 
    and compute the first residuals and variances time series using the random initial parameters (a,b,delta_0,delta,beta). (EM_core.f90)

All the time series indices are managed with a datetime derived type, provided by a side fortran module located in the "datetime-fortran" folder.
(credits: https://github.com/wavebitscientific/datetime-fortran)
