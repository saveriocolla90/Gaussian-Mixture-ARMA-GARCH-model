MODULE FITTING
	USE PARAMS
	USE TSERIES
	USE MODELS
	USE MATH
	USE DERIVATIVES
	IMPLICIT NONE

	TYPE :: PARAM_ARMA
		REAL(r9),ALLOCATABLE	:: p1(:,:)	!-Equivalent parameter (b or delta)
		REAL(r9),ALLOCATABLE	:: p2(:,:)	!-Equivalent parameter (a or beta)
	END TYPE PARAM_ARMA

	TYPE, EXTENDS(PARAM_ARMA) :: PARAM_GARCH
		REAL(r9),ALLOCATABLE	:: p0(:,:)	!-Equivalent parameter (delta_0)
	END TYPE PARAM_GARCH

	
CONTAINS

	!-Here define a subroutine to put all the derivatives together and compute the parameters increment
	!-First compute the partial derivatives for residuals (dr) and variancies (dv) , then compute the hessian matrix (H)
	!-Finally compute the matrix products: H^-1 * dr , H^-1 * dv ....(for ARMA and GARCH parameters respectively)  
	SUBROUTINE fitting_procedure(model,y)
		type(GMIX_MODEL),INTENT(INOUT)	:: model	!-Model to be considered
		type(TIME_SERIES)			:: y		!-Data to be modelled
		type(PARAM_ARMA)			:: theta	!-Parameters set of the ARMA model 
		type(PARAM_GARCH)			:: omega	!-Transformed parameters set of the GARCH model
		type(TIME_SERIES),ALLOCATABLE		:: rsq(:)	!-Squared values of residuals 
		REAL,ALLOCATABLE			:: ARMA_eps(:,:,:), GARCH_eps(:,:,:)	!-Tollerance variables to check the covergence
		INTEGER				:: i,c
		
		
		!-Set tollerance to initial unity value
		ALLOCATE(ARMA_eps(model % K, 2, max(model % R, model % S)))
		ALLOCATE(GARCH_eps(model % K, 2, max(model % Q, model % P)))
		DO c=1, model % K
			ARMA_eps(c,:,:) = 1.
			GARCH_eps(c,:,:) = 1.
		END DO
		
		!-Compute parameter's increment untill convergence (using current pZ)
		!DO WHILE(ALL(abs(ARMA_eps) > PARAM_TOLL) .OR. ALL(abs(GARCH_eps) > PARAM_TOLL))
		!DO i=1,4
			!-Call the function to increment ARMA parameters
			theta = param_incr_ARMA(model)

			!-Compute the precent change for ARMA parameters
			DO c=1, model % K
				!-Only two parameters (b,a)
				ARMA_eps(c,1,:) = ((theta % p1(c,:) - model % b(c,:)) / model % b(c,:))*100. 
				ARMA_eps(c,2,:) = ((theta % p2(c,:) - model % a(c,:)) / model % a(c,:))*100.				 
			END DO		
			print*,ARMA_eps(1,:,:)
							
			!-Set the new value 
			CALL model % set_params_ARMA(theta % p1, theta % p2)
			
			!-Update the current model series (y_hat, residuals and variances)
			model % resid = model % residuals(y)
			rsq = model % resid**2
			model % var = model % variances(rsq)
			
			!-Call the function to increment GARCH parameters
			omega = param_incr_GARCH(model)

			IF (stability(omega % p1, omega % p2)) THEN
				!-Compute the precent change for GARCH parameters
				DO c=1, model % K
					!-Only 2 parameters (gam, rho)
					GARCH_eps(c,1,:) = ((omega % p1(c,:) - model % gam(c,:)) / model % gam(c,:))*100. 
					GARCH_eps(c,2,:) = ((omega % p2(c,:) - model % rho(c,:)) / model % rho(c,:))*100.				 
				END DO		
				print*,GARCH_eps(1,:,:)
							
				!-Set the new value 
				CALL model % set_params_GARCH(omega % p0, omega % p1, omega % p2)
				CALL model % detransform()
			END IF
				
			!-Update the current model series (variances)
			model % var = model % variances(rsq)	
		!END DO
		!END DO
	END SUBROUTINE fitting_procedure


	FUNCTION param_incr_ARMA(model) RESULT(theta)
		type(GMIX_MODEL),INTENT(INOUT)		:: model
		type(PARAM_ARMA)				:: theta
		REAL,ALLOCATABLE				:: hessian(:,:,:)	!-Hessian matrix (K x M x M) with M=[R or S]	
		REAL,ALLOCATABLE				:: d_loglike(:,:)	!-Log-likelihood gradient (K x M) with M=[R or S]
		INTEGER					:: c

		!-Assign current parameter values to output parameter to be incremented 
		theta % p1 = model % b
		theta % p2 = model % a 
		
		!-Check whether the rank of the submodel's parameter B is grather then zero
		IF (model % R > 0) THEN		
			!-Compute the Hessian matrix of 2nd order log-likelihood derivatives respect parameter B 
			hessian = mix_derivate_ARMA(model, model % y_hat, model % a, model % d_b % d_resid, model % d_b % d_var)
			
			!-Compute the loglikelihood gradient respect parameter B
			d_loglike = ll_gradient(model, model % d_b % d_resid, model % d_b % d_var)
			
			!-Increment parameter B for every K component
			DO c=1, model % K
				IF (model % R > 1) THEN
					theta % p1(c,:) = theta % p1(c,:) + MATMUL(inverse(hessian(c,:,:)), d_loglike(c,:))										
				ELSE
					theta % p1(c,1) = theta % p1(c,1) + (hessian(c,1,1)**(-1.) * d_loglike(c,1))				
				END IF
			END DO
		END IF

		!-Set the new value 
		!CALL model % set_params_ARMA(theta % p1, theta % p2)

		!-Check whether the rank of the submodel's parameter A is grather then zero			
		IF (model % S > 0) THEN
			!-Compute the Hessian matrix of 2nd order log-likelihood derivatives respect parameter A 
			hessian = mix_derivate_ARMA(model, model % resid, model % a, model % d_a % d_resid, model % d_a % d_var)
			
			!-Compute the loglikelihood gradient respect parameter A
			d_loglike = ll_gradient(model, model % d_a % d_resid, model % d_a % d_var)

			!-Increment parameter A for every K component
			DO c=1, model % K
				IF (model % S > 1) THEN
					theta % p2(c,:) = theta % p2(c,:) + MATMUL(inverse(hessian(c,:,:)), d_loglike(c,:))					
				ELSE
					theta % p2(c,1) = theta % p2(c,1) + (hessian(c,1,1)**(-1.) * d_loglike(c,1))				
				END IF
			END DO
		END IF				
	END FUNCTION param_incr_ARMA	

		
		
	FUNCTION param_incr_GARCH(model) RESULT(omega)
		type(GMIX_MODEL),INTENT(INOUT)		:: model
		type(PARAM_GARCH)				:: omega
		type(TIME_SERIES),ALLOCATABLE			:: ones_series(:), resid_squared(:), zeros_d_resid(:,:)
		REAL,ALLOCATABLE				:: hessian(:,:,:)	!-Hessian matrix (K x M x M) with M=[1, Q or P]
		REAL,ALLOCATABLE				:: d_loglike(:,:)	!-Log-likelihood gradient (K x M) with M=[1, Q or P]
		INTEGER					:: c,i

		!-Assign current values for GARCH parameters to omega variable 
		omega % p0 = model % gam_0
		omega % p1 = model % gam
		omega % p2 = model % rho 
		
		!-Check whether the rank of the submodel's parameter DELTA_0 is grather then zero
		IF (ALL(model % delta_0 > 0)) THEN

			!-Set the input time series for the hessian calculation
			ALLOCATE(ones_series(model % K))
			DO c=1, model % K
				!-Set the index equivalent to index of variances
				ones_series(c) % idx = model % resid(c) % idx 
				CALL ones_series(c) % fill_ones() 
			END DO
			
			!-Compute the Hessian matrix of 2nd order log-likelihood derivatives respect parameter DELTA_0 
			hessian = mix_derivate_GARCH(model=model, ts_in=ones_series, d_var=model % d_d0 % d_var, p0=omega % p0)
			
			!-Initialize the derivatives of residuals with zeros to eliminate the residuals part in gradient computation
			ALLOCATE(zeros_d_resid(model % K,1))
			DO c=1, model % K
				!-Set the index equivalent to index of variances
				zeros_d_resid(c,1) % idx = model % resid(c) % idx 
				CALL zeros_d_resid(c,1) % fill_zeros() 
			END DO
			
			!-Compute the loglikelihood gradient respect parameter DELTA_0
			d_loglike = ll_gradient(model, zeros_d_resid, model % d_d0 % d_var)
			DEALLOCATE(zeros_d_resid)
			
			!-Increment parameter DELTA_0 for every K component
			DO c=1, model % K
				omega % p0(c,1) = omega % p0(c,1) - (hessian(c,1,1)**(-1.) * d_loglike(c,1))
			END DO			
		END IF

		!-Set the new value 
		!CALL model % set_params_GARCH(omega % p0, omega % p1, omega % p2)
		!CALL model % detransform()


		!-Check whether the rank of the submodel's parameter DELTA is grather then zero
		IF (model % Q > 0) THEN
			!-Set the input time series for the hessian calculation
			resid_squared = model % resid**2
			
			!-Compute the Hessian matrix of 2nd order log-likelihood derivatives respect parameter DELTA
			hessian = mix_derivate_GARCH(model=model, ts_in=resid_squared, d_var=model % d_delta % d_var, p1=omega % p1)
			
			!-Initialize the derivatives of residuals with zeros to eliminate the residuals part in gradient computation
			ALLOCATE(zeros_d_resid(model % K, model % Q))
			DO c=1, model % K
				DO i=1, model % Q
					!-Set the index equivalent to index of variances
					zeros_d_resid(c,i) % idx = model % resid(c) % idx 
					CALL zeros_d_resid(c,i) % fill_zeros() 
				END DO
			END DO

			!-Compute the loglikelihood gradient respect parameter DELTA
			d_loglike = ll_gradient(model, zeros_d_resid, model % d_delta % d_var)
			DEALLOCATE(zeros_d_resid)

			!-Increment parameter DELTA for every K component
			DO c=1, model % K
				IF (model % Q > 1) THEN
					omega % p1(c,:) = omega % p1(c,:) - MATMUL(inverse(hessian(c,:,:)), d_loglike(c,:))
				ELSE
					omega % p1(c,1) = omega % p1(c,1) - (hessian(c,1,1)**(-1.) * d_loglike(c,1))				
				END IF
			END DO			
		END IF


		!-Set the new value 
		!CALL model % set_params_GARCH(omega % p0, omega % p1, omega % p2)
		!CALL model % detransform()


		!-Check whether the rank of the submodel's parameter BETA is grather then zero
		IF (model % P > 0) THEN
			!-Compute the Hessian matrix of 2nd order log-likelihood derivatives respect parameter BETA 
			hessian = mix_derivate_GARCH(model=model, ts_in=model % var, d_var=model % d_beta % d_var)
			
			!-Initialize the derivatives of residuals with zeros to eliminate the residuals part in gradient computation
			ALLOCATE(zeros_d_resid(model % K, model % P))
			DO c=1, model % K
				DO i=1, model % P
					!-Set the index equivalent to index of variances
					zeros_d_resid(c,i) % idx = model % resid(c) % idx 
					CALL zeros_d_resid(c,i) % fill_zeros() 
				END DO
			END DO

			!-Compute the loglikelihood gradient respect parameter BETA
			d_loglike = ll_gradient(model, zeros_d_resid, model % d_beta % d_var)
			DEALLOCATE(zeros_d_resid)

			!-Increment parameter BETA for every K component
			DO c=1, model % K
				IF (model % P > 1) THEN
					omega % p2(c,:) = omega % p2(c,:) - MATMUL(inverse(hessian(c,:,:)), d_loglike(c,:))
				ELSE
					omega % p2(c,1) = omega % p2(c,1) - (hessian(c,1,1)**(-1.) * d_loglike(c,1))				
				END IF
			END DO
		END IF		
	END FUNCTION param_incr_GARCH
	
	
	
	FUNCTION ll_gradient(model, d_resid, d_var)
		type(GMIX_MODEL),INTENT(IN)	:: model
		type(TIME_SERIES),INTENT(IN)	:: d_resid(:,:), d_var(:,:)	!-Partial derivatives for residuals and variances respect an implicit given parameter
		REAL,ALLOCATABLE		:: ll_gradient(:,:)		!-Log-likelihood gradient (K x M) with M=[R or S]
		type(TIME_SERIES)		:: part_resid, part_var
		REAL(selected_real_kind(18))	:: total_resid, total_var
		INTEGER			:: c,i,K,M,N
		
		!-Set dimensions
		K = model % K
		M = SIZE(d_resid,DIM=2)			!-Derived parameter dimension
		N = SIZE(model % resid(1) % idx) - model % Q	!-Total observations(n) - (AR + GA) dimensions ==> n - R - Q 

		!-Allocate space for output
		ALLOCATE(ll_gradient(K,M))
		
		!-Compute the gradient for every dimension
		DO c=1, K
			DO i=1, M
				!-Cpmpute the residuals part, multiply by likelihood pZ...
				part_resid = (model % pZ(c) * model % resid(c)) &
					& * (model % var(c)**(-1.)) &
					& * d_resid(c,i)  

				!... and sum all over the time dimension
				total_resid = part_resid % sum_all()
								
				!-Compute the variance part, multiply by likelihood pZ...
				part_var = (model % pZ(c)*0.5) * (model % var(c)**(-1.)) &
					& * (((model % resid(c)**2) * (model % var(c)**(-1.))) - 1.) &
					& * d_var(c,i) 
				
				!... and sum all over the time dimension
				total_var = part_var % sum_all()
				
				!-Add together and normalize over the reduced number of observations
				ll_gradient(c,i) = (1./N)*(total_var - total_resid)
			END DO
		END DO		
	END FUNCTION ll_gradient
	
END MODULE FITTING









