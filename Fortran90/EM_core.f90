!!-----######################################################################################################################-----!!
!!---- The Mixture Arma-Garch Expectation-Maximization (MAGEM) ALTernative Unified Model (ALTUM) code implements a fitting -------!!
!!---- procedure for the ARMA(R,S)-GARCH(Q,P) model on a zero mean time series with K gaussian distributed unobserved compon -----!!
!!---- ents. It computes the maximum likelihood estimator in an iterative way, adopting the Expectation Maximization algorithm ---!!
!!---- to find the maximum Q value assuming known model's parameters values (to compute p(Z|y)) and known likelihood probabili ---!!
!!---- ties for the unobserved components Z (to compute parameters values), alternatively at each iteration. ---------------------!!
!!---- The final outputs provides the parameters values estimated for the model fit, as well as the time series of the variances -!!
!!---- representing the volatility of the data, which are somehow usefull to predict the behavior of volatility clusters in the --!!
!!---- near future.													----------!!
!!---- S.F. Colla													----------!!   
!!-----######################################################################################################################-----!!


!!! START PROGRAM !!!
PROGRAM MAGEM_ALTUM

	!--Modules used--!
	USE IO
	USE TSERIES
	USE MODELS
	USE FITTING
	USE STAT
	USE PARAMS
	IMPLICIT NONE
	
	!-Input file name for reading the series to be analyzed
	CHARACTER(LEN=:), ALLOCATABLE	:: file_name

	!-Time series declaration
	type(TIME_SERIES)		:: price_returns
	type(TIME_SERIES),ALLOCATABLE	:: rsq(:)

	!-Model variable declaration
	type(GMIX_MODEL) 	:: mix_model

	!-Aux for parameters fitting
	REAL(r18)		:: Q_value
	INTEGER		:: i,M,j
	
	
	!-Read data from file
	file_name = "./RUNEUSDT_4H_conv.txt"
	CALL price_returns % load(file_name)

	!-Initialize the Mixture model name and number of components
	print*,NL,"Initialize the Gaussian Mixture ARMA-GARCH model:"
	mix_model = GMIX_MODEL(K,R,S,Q,P) 
	mix_model % name = "Model 1"
	
	!-Display initial parameters values
	print*,NL,"INITIALIZED PARAMETER VALUES:"
	CALL display(mix_model)
	
	!-First compute the residuals from the ARMA model (given the time series to be analized) ...
	mix_model % resid = mix_model % residuals(price_returns)
	print*,NL,"Initial random ARMA modelled Y_HAT:"
	CALL mix_model % y_hat(1) % display()

	print*,NL,"Initial random ARMA model RESIDUALS:"
	CALL mix_model % resid(1) % display()
	
	!...then compute the variance from the GARCH module (given the residuals) 
	print*,NL,"Initial random GARCH model VARIANCIES:"
	rsq = mix_model % resid**2
	mix_model % var = mix_model % variances(rsq)
	CALL mix_model % var(1) % display()

	

	!!---------------- This part needs to be iterated untill the convergence of Q value --------------------!!

		!-Compute likelihood for unobserved components Z
		print*,NL,"Initial likelihoods P(Z|y):"	
		mix_model % pZ = prob_pZ(mix_model % resid, mix_model % var, mix_model % alpha) 
		CALL mix_model % pZ(1) % display()
		
		!-Fit the Gaussian Mixture ARMA-GARCH model on the price return time series
		CALL fitting_procedure(mix_model, price_returns)

		!-Compute the total log-likelihood (Q func)
		Q_value = Q_func(mix_model % resid, mix_model % var, mix_model % alpha, mix_model % pZ)		
		print*,NL,"Q value:",Q_value,NL
		
		!-Compute new Z probabilities (alpha)
		DO i=1, mix_model % K
			mix_model % alpha(i) = SUM(mix_model % pZ(i) % value) / SIZE(mix_model % var(i) % idx)
		END DO
	!!------------------------------------------------------------------------------------------------------!!
		
	!-Display final parameters estimated
	print*,NL,"FINAL PARAMETER VALUES:"
	CALL display(mix_model)
	
	print*,NL,"Everything OK!"


END PROGRAM MAGEM_ALTUM


