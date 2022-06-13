MODULE STAT
	USE TSERIES
	USE PARAMS
	IMPLICIT NONE
	
	TYPE, PRIVATE :: NORMAL
		REAL	:: mean, std
	CONTAINS
		PROCEDURE	:: pdf => normal_pdf	
		!PROCEDURE	:: cdf	=> normal_cdf
	END TYPE

	TYPE, PRIVATE :: T_STUDENT
		REAL	:: dof
	CONTAINS
		PROCEDURE	:: pdf => t_student_pdf	
		!PROCEDURE	:: cdf	=> t_student_cdf	
	END TYPE


CONTAINS
	FUNCTION prob_pZ(resid, var, alpha) RESULT(pZ)
		type(TIME_SERIES),INTENT(IN)	:: resid(:), var(:)	
		REAL,INTENT(IN)		:: alpha(:)
		type(TIME_SERIES),ALLOCATABLE	:: pZ(:)
		type(NORMAL)			:: norm
		REAL(r27),ALLOCATABLE		:: normalization(:)
		INTEGER			:: i,K
		
		!-Set standard Gaussian
		norm % mean = 0.
		norm % std = 1.
		
		!-Allocate memory for output results
		K = SIZE(alpha)
		ALLOCATE(pZ(K))
		ALLOCATE(normalization(SIZE(var(1) % idx)))
		
		!-Compute normalization factors
		DO i=1,K
			normalization = normalization &
					+ (alpha(i) * norm % pdf(resid(i) % value(R+1:) / sqrt(var(i) % value)) &
					/ sqrt(var(i) % value))
		END DO		

		!-Iterate over all K components to get every single likelihood for unobserved variable Z
		DO i=1,K
			!-Compute the Z likelihood values 
			pZ(i) % value = ((alpha(i) / sqrt(var(i) % value)) &
					* norm % pdf(resid(i) % value(R+1:) / sqrt(var(i) % value))) &
					/ normalization

			!-Set the datetime index
			pZ(i) % idx = var(i) % idx 	
		END DO	
		DEALLOCATE(normalization)
	END FUNCTION prob_pZ
	
	
	FUNCTION Q_func(resid, var, alpha, pZ)
		type(TIME_SERIES),INTENT(IN)	:: resid(:), var(:), pZ(:)	
		REAL,INTENT(IN)		:: alpha(:)
		REAL(r18)			:: Q_func
		REAL(r18),ALLOCATABLE 		:: Q_star(:)
		type(NORMAL)			:: norm
		INTEGER			:: c

		!-Set standard Gaussian
		norm % mean = 0.
		norm % std = 1.

		!-Allocate the Q* series of N - Q - P observations
		ALLOCATE(Q_star(SIZE(var(1) % idx)))
		Q_star = 0.
		
		!-Add together the Q* for every K component...
		DO c=1, SIZE(alpha)
			Q_star = Q_star + (pZ(c) % value * log((alpha(c) * norm % pdf(resid(c) % value(R+1:) / sqrt(var(c) % value)) &
					/ sqrt(var(c) % value))))
		END DO
		
		!... and sum all over the total observations
		Q_func = SUM(Q_star)
		DEALLOCATE(Q_star)
	END FUNCTION Q_func
	
	
	FUNCTION normal_pdf(self,x)
		class(NORMAL)		:: self
		REAL(r27),INTENT(IN)	:: x(:)
		REAL(r27)		:: normal_pdf(SIZE(x))
		
		!-Probability density function for a normal distributed random variable
		normal_pdf = exp(-((x - self % mean) / self % std)**2 / 2.) / (self % std * sqrt(2.*PI))
	END FUNCTION normal_pdf
	

	FUNCTION t_student_pdf(self,x)
		class(T_STUDENT)	:: self
		REAL,INTENT(IN)	:: x(:)
		REAL			:: t_student_pdf(SIZE(x))
	END FUNCTION t_student_pdf
	
END MODULE STAT
