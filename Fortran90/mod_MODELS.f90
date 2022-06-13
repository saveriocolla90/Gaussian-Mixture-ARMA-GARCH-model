MODULE MODELS
	USE TSERIES
	USE PARAMS
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: stability
	
	! RESIDUALS and VARIANCES (derivatives)	
	TYPE :: DERIVATE
		type(TIME_SERIES),ALLOCATABLE	:: d_resid(:,:) !-Time series of residuals
		type(TIME_SERIES),ALLOCATABLE	:: d_var(:,:)	!-Time series of variances				
	END TYPE DERIVATE
		
	
	! ARMA model ranks data type
	TYPE, PUBLIC :: ARMA
		INTEGER			:: R, S
		REAL(r9),ALLOCATABLE		:: a(:,:) 
		REAL(r9),ALLOCATABLE		:: b(:,:)
		type(TIME_SERIES),ALLOCATABLE	:: y_hat(:)
		type(TIME_SERIES),ALLOCATABLE	:: resid(:)
		type(DERIVATE)			:: d_a, d_b
	CONTAINS
		PROCEDURE :: residuals => ARMA_residuals
		PROCEDURE :: set_params_ARMA
	END TYPE ARMA


	! GARCH model ranks data type
	TYPE, PUBLIC, EXTENDS(ARMA) :: ARMA_GARCH
		INTEGER			:: Q, P
		REAL(r9),ALLOCATABLE		:: delta_0(:,:), gam_0(:,:)
		REAL(r9),ALLOCATABLE		:: delta(:,:), gam(:,:)
		REAL(r9),ALLOCATABLE		:: beta(:,:), rho(:,:)
		type(TIME_SERIES),ALLOCATABLE	:: var(:)
		type(DERIVATE)			:: d_d0, d_delta, d_beta
	CONTAINS
		PROCEDURE :: init_ARMA_GARCH => random_init_ARMA_GARCH
		PROCEDURE :: variances => GARCH_variances
		PROCEDURE :: set_params_GARCH
		PROCEDURE :: transform, detransform
	END TYPE ARMA_GARCH


	! Gaussian Mixture ARMA-GARCH model with K components 
	TYPE, PUBLIC, EXTENDS(ARMA_GARCH) :: GMIX_MODEL
		CHARACTER(LEN=:), ALLOCATABLE	:: name
		type(TIME_SERIES),ALLOCATABLE	:: pZ(:)	!-Likelihood for unobserved component Z	
		REAL, ALLOCATABLE		:: alpha(:)	!-Probabilities for unobserved component Z
		INTEGER			:: K 
	END TYPE GMIX_MODEL
	INTERFACE GMIX_MODEL
		MODULE PROCEDURE random_init_GMIX
	END INTERFACE GMIX_MODEL


CONTAINS
!!--- SUBROUTINES & FUNCTIONS ---!!	
	
	!!-Set the parameters for ARMA model to the input ones
	SUBROUTINE set_params_ARMA(self,b,a)
		class(ARMA)		:: self
		REAL(r9),ALLOCATABLE	:: a(:,:) 
		REAL(r9),ALLOCATABLE	:: b(:,:)
		
		self % a = a
		self % b = b
	END SUBROUTINE set_params_ARMA
	
	!!-Set the parameters for GARCH model to the input ones
	SUBROUTINE set_params_GARCH(self,gam_0,gam,rho)
		class(ARMA_GARCH)	:: self
		REAL(r9),ALLOCATABLE	:: gam_0(:,:)
		REAL(r9),ALLOCATABLE	:: gam(:,:)
		REAL(r9),ALLOCATABLE	:: rho(:,:)
		
		self % gam_0 = gam_0
		self % gam = gam
		self % rho = rho
	END SUBROUTINE set_params_GARCH


	!!-Transform to non-negative parameters	
	SUBROUTINE transform(self)
		class(ARMA_GARCH), INTENT(INOUT)	:: self

		!-Transform to exponential indices
		self % gam_0 = LOG(self % delta_0)
		self % gam = LOG(self % delta)
		self % rho = LOG(self % beta)		
	END SUBROUTINE transform


	!!-Transform to non-negative parameters	
	SUBROUTINE detransform(self)
		class(ARMA_GARCH), INTENT(INOUT)	:: self

		!-Transform to logarithmic base
		self % delta_0 = EXP(self % gam_0)
		self % delta = EXP(self % gam)
		self % beta = EXP(self % rho)	
	END SUBROUTINE detransform


	!-ARMA-GARCH model parameters random initialization
	SUBROUTINE random_init_ARMA_GARCH(self,R,S,Q,P,K)
		class(ARMA_GARCH), INTENT(INOUT)	:: self
		INTEGER,INTENT(IN)			:: R,S,Q,P,K
		LOGICAL				:: check_stability = .False.	!-Will be used to get result for stability check
		INTEGER 				:: i

		!-Fix parameter dimensions
		self % R = R
		self % S = S
		self % Q = Q
		self % P = P

		!-Allocate memory for parameters...
		ALLOCATE(self % a(K,S))
		ALLOCATE(self % b(K,R))		
		ALLOCATE(self % delta_0(K,1))
		ALLOCATE(self % delta(K,Q))
		ALLOCATE(self % beta(K,P))

		!... for readisula and variances...
		ALLOCATE(self % y_hat(K))
		ALLOCATE(self % resid(K))
		ALLOCATE(self % var(K))

		!... and their derivatives (ARMA)
		ALLOCATE(self % d_a % d_resid(K,S))
		ALLOCATE(self % d_b % d_resid(K,R))
		ALLOCATE(self % d_a % d_var(K,S))
		ALLOCATE(self % d_b % d_var(K,R))

		!... and their derivatives (GARCH)
		ALLOCATE(self % d_d0 % d_var(K,1))
		ALLOCATE(self % d_delta % d_var(K,Q))
		ALLOCATE(self % d_beta % d_var(K,P))
		
		!-Finally fill them with random numbers (0,1)
		CALL RANDOM_NUMBER(self % a)
		CALL RANDOM_NUMBER(self % b)
		CALL RANDOM_NUMBER(self % delta_0)
		
		self % a = abs(self % a - 0.5)
		self % b = abs(self % b - 0.5)
		
		!-Call RANDOM_NUMBER until match for stability condition (delta + beta < 1) for GARCH parameters (over K components)
		DO WHILE(check_stability .eqv. .False.)
			CALL RANDOM_NUMBER(self % delta)
			CALL RANDOM_NUMBER(self % beta)

			!-Rescale over the K components
			self % delta = self % delta / Q
			self % beta = self % beta / P
			
			!-Check stability	
			check_stability = stability(self % delta, self % beta)
		END DO
				
		!-Transform to non-negative parameters with exponential indices		
		CALL self % transform()
	END SUBROUTINE random_init_ARMA_GARCH



	!-GAUSSIAN MIXTURE model parameters initialization
	!-Call the subroutine for ARMA-GARCH random initialization 
	FUNCTION random_init_GMIX(K,R,S,Q,P) RESULT(self)
		TYPE(GMIX_MODEL)		:: self
		INTEGER,INTENT(IN)		:: K,R,S,Q,P

		!-Allocate memory for each Gaussian component
		ALLOCATE(self % pZ(K)) 
		ALLOCATE(self % alpha(K))

		!-Initialize the number of the Gaussian mixture components 
		self % K = K
				
		!-Initialize the components probabilities with equal values
		self % alpha = 1. / K
		
		!-Finally initialize the sub-models parameters
		CALL self % init_ARMA_GARCH(R,S,Q,P,K)
	END FUNCTION random_init_GMIX



	!-Compute residuals time series from ARMA model with current parameters values
	FUNCTION ARMA_residuals(model,series) RESULT(residuals)
		class(ARMA),INTENT(INOUT)	:: model				!-Model ARMA: containing parameters (a,b,R,S) and y_hat series
		type(TIME_SERIES),INTENT(IN)	:: series				!-Time series to be modelled
		type(TIME_SERIES)		:: residuals(SIZE(model % a,DIM=1))	!-Residuals terms: difference between real data and the model
		INTEGER			:: i,t,M,K
				
		!-Start two nasted loops: over K components and along M size
		M = SIZE(series % idx)
		K = SIZE(model % a,DIM=1)
		DO i=1,K
			!-Allocate space for the output series of rank K
			ALLOCATE(residuals(i) % idx(M - model % R))
			ALLOCATE(residuals(i) % value(M - model % R))
			
			!-Allocate space for the modelled data predictions time series (Y_HAT)
			IF (.not. ALLOCATED(model % y_hat(i) % idx)) ALLOCATE(model % y_hat(i) % idx(M - model % R))
			IF (.not. ALLOCATED(model % y_hat(i) % value)) ALLOCATE(model % y_hat(i) % value(M - model % R))
						
			!-Iterate until the end of the input series
			DO t = model % R, M-1
			
				!-AR part of the model
				residuals(i) % idx(t-model % R+1) = series % idx(t+1)
				residuals(i) % value(t-model % R+1) = series % value(t+1) &
									&- DOT_PRODUCT(model % b(i,:), series % value(t:t-model % R+1:-1))
				
				!-MA part of the model
				IF (t .ge. (model % R + model % S)) THEN		
					residuals(i) % value(t-model % R+1) =  residuals(i) % value(t-model % R+1) & 
										&- DOT_PRODUCT(model % a(i,:), &
										& residuals(i) % value(t-model % R:t-model % R-model % S+1:-1))
				END IF
				
				!-Fill the Y_HAT series (i.e.: modelled data predictions)
				model % y_hat(i) % idx(t-model % R+1) = series % idx(t+1)
				model % y_hat(i) % value(t-model % R+1) = series % value(t+1) - residuals(i) % value(t-model % R+1)
			END DO 
		END DO
	END FUNCTION ARMA_residuals



	!-Compute variances time series for GARCH model with current parameter values 
	FUNCTION GARCH_variances(model,series) RESULT(var)
		class(ARMA_GARCH),INTENT(IN)	:: model				!-Model GARCH: containing parameters (delta_0,delta,beta,Q,P)
		type(TIME_SERIES),INTENT(IN)	:: series(SIZE(model % delta,DIM=1))	!-Time series to be modelled (...residuals from ARMA model)
		type(TIME_SERIES)		:: var(SIZE(model % delta,DIM=1))	!-Modelled variances 
		INTEGER			:: i,t,M,K
				
		!-Start two nasted loops: over K components and along M size
		M = SIZE(series(1) % idx)
		K = SIZE(model % delta,DIM=1)
		DO i=1,K
			 
			!-Allocate space for the output series of rank K
			ALLOCATE(var(i) % idx(M - model % Q))
			ALLOCATE(var(i) % value(M - model % Q))

			!-Iterate until the end of the input series
			DO t = model % Q, M-1
			
				!-GAR part of the model
				var(i) % idx(t-model % Q+1) = series(i) % idx(t+1)
				var(i) % value(t-model % Q+1) = model % delta_0(i,1) &
								   &+ DOT_PRODUCT(model % delta(i,:), series(i) % value(t:t-model % Q+1:-1))
				
				!-CH part of the model
				IF (t .ge. (model % Q + model % P)) THEN
					var(i) % value(t-model % Q+1) =  var(i) % value(t-model % Q+1) &
									    &+ DOT_PRODUCT(model % beta(i,:), &
									    & var(i) % value(t-model % Q:t-model % Q-model % P+1:-1))
				END IF
			END DO
		END DO
	END FUNCTION GARCH_variances		
	
	
		
	!-Return the stability logical value for the GARCH model parameters (delta + beta < 1)
	FUNCTION stability(d,b)
		REAL(r9),INTENT(IN)	:: d(:,:), b(:,:)		!-Delta and Beta alias
		REAL(r9)		:: db_sum(SIZE(d,DIM=1))	!-To compute their total sum
		LOGICAL		:: stability			!-Returned logical value
		INTEGER		:: i,j
		
		!-Compute the sum of all delta and beta values 
		db_sum = 0		
		db_sum = SUM(d, DIM=2)
		db_sum = db_sum + SUM(b, DIM=2)
		
		!-Check if all sum elements are < 1
		stability = ALL(db_sum < 1)
	END FUNCTION stability


END MODULE MODELS








