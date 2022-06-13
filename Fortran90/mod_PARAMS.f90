MODULE PARAMS
	IMPLICIT NONE
	PUBLIC 	
	
	!--MODEL PARAMETERS--!
	INTEGER, PARAMETER		:: N = 300	!-Max number of observations
	INTEGER, PARAMETER		:: K = 2	!-Gaussian mixture number of components
	INTEGER, PARAMETER		:: R = 2
	INTEGER, PARAMETER		:: S = 2
	INTEGER, PARAMETER		:: Q = 2
	INTEGER, PARAMETER		:: P = 2

	!--TOLLERANCE--!
	REAL, PARAMETER		:: PARAM_TOLL = 0.1	!-Tollerance (%) for parameters convergence 

	!--MATH CONSTANTS--!
	REAL, PARAMETER		:: PI = ACOS(-1.)	!-Pi greco

	!--STRINGS FORMAT--!
	CHARACTER(LEN=80)		:: LINE				!-Length of a readble line		
	CHARACTER(LEN=13), PARAMETER	:: DATETIME_FMT = '%Y%m%d %H%M%S'	!-Format for readble date-time 

	!!--- SPECIAL CHARACTERS --!!
	CHARACTER (LEN=1) :: TAB = char(9) 		!-Tab space
	CHARACTER (LEN=1) :: NL = NEW_LINE(NL)	!-New line
	
	!--KIND CONVERSIONS--!
	INTEGER, PARAMETER		:: i9 = selected_int_kind(9)
	INTEGER, PARAMETER		:: i18 = selected_int_kind(18)
	INTEGER, PARAMETER		:: i27 = selected_int_kind(27)
	
	INTEGER, PARAMETER		:: r9 = selected_real_kind(9)
	INTEGER, PARAMETER		:: r18 = selected_real_kind(18)
	INTEGER, PARAMETER		:: r27 = selected_real_kind(27)
	
	!--DATETIME CONSTANTS--!
	REAL(r9),PARAMETER		:: DAY_SEC = 86400.
	REAL(r9),PARAMETER		:: HOUR_SEC = 3600.
	REAL(r9),PARAMETER		:: MINUTE_SEC = 60. 

	REAL(r9),PARAMETER		:: REF_SEC = 621356832d2	 	!-Reference seconds:  1970-01-01 00:00:00   minus   0001-01-01 00:00:00 UTC
	REAL(r9),PARAMETER		:: STANDARD_DAYS = REF_SEC / DAY_SEC 	!-Total days in reference seconds
	
 
END MODULE PARAMS
