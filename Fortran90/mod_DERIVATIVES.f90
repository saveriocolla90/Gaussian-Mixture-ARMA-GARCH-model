MODULE DERIVATIVES
	USE TSERIES
	USE MODELS
	USE PARAMS
	IMPLICIT NONE
	PRIVATE
	
	PRIVATE :: d_resid_ARMA, d_var_ARMA, d_var_GARCH
	PUBLIC :: mix_derivate_ARMA, mix_derivate_GARCH
	
	INTERFACE partial_derivate
		MODULE PROCEDURE d_resid_ARMA, d_var_ARMA, d_var_GARCH
	END INTERFACE partial_derivate
			

CONTAINS

	!!-Derivatives of residuals respect to an ARMA parameter
	FUNCTION d_resid_ARMA(z, a, M)
		type(TIME_SERIES),INTENT(IN)	:: z(:)	!-Reference time series (original data or residuals)
		REAL(r9), INTENT(IN)		:: a(:,:)	!-Parameters matrix (K x parameter rank S (ARMA parameter "a"))
		INTEGER,INTENT(IN)		:: M		!-Reference parameter rank (R or S)
		type(TIME_SERIES)		:: d_resid_ARMA(SIZE(a,DIM=1),SIZE(a,DIM=2))	!-Output derivatives matrix 
		INTEGER			:: i,j,t,K,S,N
		
		!-Set the dimensions
		K = SIZE(z,DIM=1)				!-Unobserved varibales
		S = SIZE(a,DIM=2)				!-Parameters rank for moveing avarage part of ARMA model
		N = SIZE(z(1) % idx)				!-Time index (total number of reference data)
		
		!-Iterate over K components
		DO i=1,K
			!-Iterate over M parameters
			DO j=1,M				
				!-Allocate time series memory based on parameters rank M
				ALLOCATE(d_resid_ARMA(i,j) % idx(N-M))			
				ALLOCATE(d_resid_ARMA(i,j) % value(N-M))

				!-Compute the derivate and fill the series values
				DO t=M, N-1
					d_resid_ARMA(i,j) % value(t-M+1) = -z(i) % value(t-j+1) 
					
					!-Check to include the autoregressive part
					IF (t .ge. (M+S)) THEN
						d_resid_ARMA(i,j) % value(t-M+1) = d_resid_ARMA(i,j) % value(t-M+1) &
										& - DOT_PRODUCT(a(i,:), d_resid_ARMA(i,j) % value(t-M:t-M-S+1:-1))
					END IF
				END DO

				!-Finaly set the datetime index
				d_resid_ARMA(i,j) % idx = z(i) % idx(M+1:)
			END DO
		END DO
	END FUNCTION d_resid_ARMA
	
	
	!!-Derivatives of variances respect to an ARMA parameter
	FUNCTION d_var_ARMA(z, dz, gam, rho)
		type(TIME_SERIES),INTENT(IN)	:: z(:), dz(:,:)	!-Reference time series and its derivate respect of ARMA parameters (residuals terms)
		REAL(r9),INTENT(IN)		:: gam(:,:), rho(:,:)	!-Parameters matrices (K x parameter rank (Q or P))
		type(TIME_SERIES)		:: d_var_ARMA(SIZE(dz,DIM=1),SIZE(dz,DIM=2))	!-Output derivatives matrix 
		INTEGER			:: i,j,t,K,Q,P,M,N
		
		!-Set the dimensions
		K = SIZE(z,DIM=1)					!-Unobserved varibales
		M = SIZE(dz,DIM=2)					!-Reference parameter rank from ARMA model (R or S)
		Q = SIZE(gam,DIM=2)					!-Parameters rank for autoregressive part of GARCH model
		P = SIZE(rho,DIM=2)					!-Parameters rank for conditional heteroskedasticity part of GARCH model		
		N = SIZE(dz(1,1) % idx)					!-Time index (total number of reference data)
		
		!-Iterate over K components
		DO i=1,K
			!-Iterate over M parameters
			DO j=1,M				
				!-Allocate time series memory based on parameters rank Q
				ALLOCATE(d_var_ARMA(i,j) % idx(N-Q))			
				ALLOCATE(d_var_ARMA(i,j) % value(N-Q))

				!-Compute the derivate and fill the series values
				DO t=Q, N-1
					d_var_ARMA(i,j) % value(t-Q+1) = 2.*DOT_PRODUCT(z(i) % value(t:t-Q+1:-1) * exp(gam(i,:)), dz(i,j) % value(t:t-Q+1:-1)) 

					!-Check to include the autoregressive part
					IF (t .ge. (Q+P)) THEN
						d_var_ARMA(i,j) % value(t-Q+1) = d_var_ARMA(i,j) % value(t-Q+1) &
										& + DOT_PRODUCT(exp(rho(i,:)), d_var_ARMA(i,j) % value(t-Q:t-Q-P+1:-1)) 
					END IF
				END DO

				!-Finaly set the datetime index
				d_var_ARMA(i,j) % idx = dz(i,j) % idx(Q+1:)
			END DO
		END DO		
	END FUNCTION d_var_ARMA
	
	
	!!-Derivatives of variances respect to a GARCH parameter
	FUNCTION d_var_GARCH(z, rho, p0, p1)
		type(TIME_SERIES),INTENT(IN)	:: z(:)
		REAL(r9),INTENT(IN)		:: rho(:,:)		!-Parameter for partial computation respect RHO
		REAL(r9),INTENT(IN),OPTIONAL	:: p0(:,:), p1(:,:)	!-Parameter for partial computation respect GAM_0 and GAM   
		type(TIME_SERIES),ALLOCATABLE	:: d_var_GARCH(:,:)
		INTEGER			:: i,j,t,K,P,M,N
		
		!-Allocate memory for output series based on optional parameters dimension
		IF (PRESENT(p0))  ALLOCATE(d_var_GARCH(SIZE(p0,DIM=1),1))							!-Dim: K x 1
		IF (PRESENT(p1))  ALLOCATE(d_var_GARCH(SIZE(p1,DIM=1),SIZE(p1,DIM=2)))					!-Dim: K x Q
		IF ((.not. PRESENT(p0)) .and. (.not. PRESENT(p1))) ALLOCATE(d_var_GARCH(SIZE(rho,DIM=1),SIZE(rho,DIM=2)))	!-dim: K x P

		!-Set the dimensions
		K = SIZE(z,DIM=1)					!-Unobserved varibales
		M = SIZE(d_var_GARCH,DIM=2)				!-Reference parameter rank from GARCH model (1, Q or P)
		P = SIZE(rho,DIM=2)					!-Parameters rank for conditional heteroskedasticity part of GARCH model		
		N = SIZE(z(1) % idx)					!-Time index (total number of reference data)

		
		!-Iterate over K components
		DO i=1,K
			!-Iterate over M parameters
			DO j=1,M				
				!-Allocate time series memory based on parameters rank P
				ALLOCATE(d_var_GARCH(i,j) % idx(N-M))			
				ALLOCATE(d_var_GARCH(i,j) % value(N-M))

				!-Compute the derivate and fill the series values
				DO t=M, N-1
					!-First series elements computation depends on the input parameter
					IF (PRESENT(p0)) d_var_GARCH(i,j) % value(t-M+1) = exp(p0(i,1))
					IF (PRESENT(p1)) d_var_GARCH(i,j) % value(t-M+1) = z(i) % value(t-j+1) * exp(p1(i,j))
					IF ((.not. PRESENT(p0)) .and. (.not. PRESENT(p1))) d_var_GARCH(i,j) % value(t-M+1) = z(i) % value(t-j+1) * exp(rho(i,j))
					
					!-Check to include the autoregressive part
					IF (t .ge. (M+P)) THEN
						d_var_GARCH(i,j) % value(t-M+1) = d_var_GARCH(i,j) % value(t-M+1) &
										& + DOT_PRODUCT(exp(rho(i,:)), d_var_GARCH(i,j) % value(t-M:t-M-P+1:-1))
					END IF
				END DO
				!-Finaly set the datetime index
				d_var_GARCH(i,j) % idx = z(i) % idx(M+1:)
			END DO
		END DO				
	END FUNCTION d_var_GARCH
	
	
	!!-Hessian matrix for ARMA parameters 
	FUNCTION mix_derivate_ARMA(model, ts_in, a, d_resid, d_var)
		type(GMIX_MODEL),INTENT(IN)	:: model			!-The Mixture model to be considered
		type(TIME_SERIES),INTENT(IN)	:: ts_in(:)			!-Input time series of reference (residuals or original data)
		REAL(r9),INTENT(IN)		:: a(:,:)			!-Parameter from MA part of the model
		type(TIME_SERIES),INTENT(OUT)	:: d_resid(:,:), d_var(:,:)	!-Partial derivatives to be filled	
		REAL(r9),ALLOCATABLE		:: mix_derivate_ARMA(:,:,:)	!-Hessian matrix 
		type(TIME_SERIES)		:: part_resid, part_var, total !-Partial computation variables for Hessian matrix values
		REAL(r9)			:: total_sum			!-Sum over the time dimension of the series
		INTEGER			:: c,i,j,K,M,N
		
		!-Set the dimensions
		K = model % K			
		M = SIZE(d_resid,DIM=2)			!-Derived parameter dimension [R or S]
		N = SIZE(model % resid(1) % idx) - model % Q	!-Total observations(n) - (AR + GA) dimensions ==> n - R - Q 
		
		!-Allocate matrix dimensions for the output result... 
		!... (K x M x M where M=[R or S] depending the parameter under consideration) 
		ALLOCATE(mix_derivate_ARMA(K,M,M))

		!-Compute partial derivatives of residuals and variances respect the input parameter
		d_resid = partial_derivate(ts_in, a, SIZE(d_resid,DIM=2))		
		d_var = partial_derivate(model % resid, d_resid, model % gam, model % rho)
		
		!-Compute every element of the Hessian and fill the matrix for each of the K components
		DO c=1,K
			DO i=1,M
				DO j=1,M
					!-Residuals derived components (i,j) combiantion
					part_resid = d_resid(c,i) * d_resid(c,j) * (model % var(c)**(-1))
					
					!-Variances derived components (i,j) combiantion
					part_var = d_var(c,i) * d_var(c,j) * ((model % var(c)**(-2)) * 0.5)
					
					!-Multiply by component likelihood pZ  
					total = (part_resid + part_var) * model % pZ(c)
					
					!... and Sum over all observations
					total_sum = total % sum_all()
					
					!-Rescale the total sum by the effective number of data (N)
					mix_derivate_ARMA(c,i,j) = (1./N)*total_sum
				END DO
			END DO
		END DO
	END FUNCTION mix_derivate_ARMA


	!!-Hessian matrix for GARCH parameters 
	FUNCTION mix_derivate_GARCH(model, ts_in, d_var, p0, p1)
		type(GMIX_MODEL),INTENT(IN)	:: model			!-Model to be considered
		type(TIME_SERIES),INTENT(IN)	:: ts_in(:)			!-Input time series of reference ({1}, residuals**2 or variance)
		REAL(r9),INTENT(IN),OPTIONAL	:: p0(:,:), p1(:,:)		!-Parameters determining the behavior of the function (GAM_0 or GAM)
		type(TIME_SERIES),INTENT(OUT)	:: d_var(:,:)			!-Partial derivative to be filled
		REAL(r9),ALLOCATABLE		:: mix_derivate_GARCH(:,:,:)	!-Hessian matrix 
		type(TIME_SERIES)		:: part_var, total 		!-Partial computation variables for Hessian matrix values
		REAL(r9)			:: total_sum			!-Sum over the time dimension of the series
		INTEGER			:: c,i,j,K,M,N
		
		!-Set the dimensions
		K = model % K			
		M = SIZE(d_var,DIM=2)				!-Derived parameter dimension [1, Q or P]
		N = SIZE(model % resid(1) % idx) - model % Q	!-Total observations(n) - (AR + GA) dimensions ==> n - R - Q 
		
		!-Allocate matrix dimensions for the output result... 
		!... (K x M x M where M=[1, Q or P] depending the parameter under consideration) 
		ALLOCATE(mix_derivate_GARCH(K,M,M))
				
		!-Compute partial derivatives of residuals and variances respect the input parameter
		IF (PRESENT(p0)) d_var = partial_derivate(ts_in, model % rho, p0=p0)
		IF (PRESENT(p1)) d_var = partial_derivate(ts_in, model % rho, p1=p1)		
		IF ((.not. PRESENT(p0)) .and. (.not. PRESENT(p1))) d_var = partial_derivate(ts_in, model % rho)
		
		!-Compute every element of the Hessian and fill the matrix for each of the K components
		DO c=1,K
			DO i=1,M
				DO j=1,M
					!-Variances derived components (i,j) combiantion
					part_var = d_var(c,i) * d_var(c,j) * ((model % var(c)**(-2)) * 0.5)
					
					!-Multiply by component likelihood pZ  
					total = part_var * model % pZ(c)
					
					!... and Sum over all observations
					total_sum = total % sum_all()
					
					!-Rescale the total sum by the effective number of data (N)
					mix_derivate_GARCH(c,i,j) = (1./N)*total_sum				
				END DO
			END DO
		END DO		
	END FUNCTION mix_derivate_GARCH

	
END MODULE DERIVATIVES
