MODULE TSERIES
	USE PARAMS
	USE datetime_module
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: operator(+), operator(-), operator(*), operator(/), operator(**)
	
	
	!!--- DERIVIED TYPES ---!!
	TYPE, PUBLIC :: TIME_SERIES
		type(datetime), ALLOCATABLE	:: idx(:)		!-Date-time index
		INTEGER,ALLOCATABLE		:: timestamp(:)	!-Timestamp index
		REAL(r27), ALLOCATABLE		:: value(:)	
	CONTAINS
		PROCEDURE :: load => load_data
		PROCEDURE :: display => print_data
		PROCEDURE :: to_timestamp
		PROCEDURE :: to_datetime
		PROCEDURE :: sum_all
		PROCEDURE :: fill_ones
		PROCEDURE :: fill_zeros		
	END TYPE TIME_SERIES



	!!--- INTERFACE FOR BINARY OPERATORS ---!!
	INTERFACE OPERATOR(+)
		MODULE PROCEDURE add
		!MODULE PROCEDURE MD_add
		MODULE PROCEDURE int_add
		MODULE PROCEDURE real_add
	END INTERFACE OPERATOR(+)

	INTERFACE OPERATOR(-)
		MODULE PROCEDURE sub
		!MODULE PROCEDURE MD_sub
		MODULE PROCEDURE int_sub
		MODULE PROCEDURE real_sub
	END INTERFACE OPERATOR(-)

	INTERFACE OPERATOR(*)
		MODULE PROCEDURE mult
		!MODULE PROCEDURE MD_mult
		MODULE PROCEDURE int_mult
		MODULE PROCEDURE real_mult
	END INTERFACE OPERATOR(*)

	INTERFACE OPERATOR(/)
		MODULE PROCEDURE div
		!MODULE PROCEDURE MD_div
		MODULE PROCEDURE int_div
		MODULE PROCEDURE real_div
	END INTERFACE OPERATOR(/)

	INTERFACE OPERATOR(**)
		MODULE PROCEDURE int_power
		MODULE PROCEDURE real_power
		MODULE PROCEDURE int_MD_power
		MODULE PROCEDURE real_MD_power
	END INTERFACE OPERATOR(**)



CONTAINS
!!=============================================================!!
!!---------------- SUBROUTINES & FUNCTIONS --------------------!!
!!=============================================================!!
	



	!!------------------------------------------!!
	!!------- INPUT / OUTPUT SUBROUTINES -------!!
	!!------------------------------------------!!
		
	!!---- Load time series data from a formatted .txt file composed by 3 columns (date, time, value)  ---------!!
	!!---- and allocate space for the time series derived data type based on the file's length (limited by N) --!!
	SUBROUTINE load_data(self,file_name)
		class(TIME_SERIES), INTENT(OUT)	:: self
		CHARACTER(LEN=:), ALLOCATABLE		:: file_name
		CHARACTER(LEN=8)			:: date
		CHARACTER(LEN=6)			:: time
		INTEGER				:: i,res,M
		
		!-Open the formatted file to read from 
		OPEN(UNIT=1, FILE=file_name, FORM="FORMATTED", IOSTAT=res)
			IF (res /= 0) THEN
				print*,"Error in opening file:",file_name
				print*,"Err status:",res
				STOP
			END IF
			
			!-Count the total lines in the .txt file
			i = 0
			DO WHILE(.true.)
				READ(1,'(A)',END=99) LINE
				i = i + 1
			END DO
			99 CONTINUE
			CLOSE(1)
		print*,"Reading file:",file_name
		
		!-Set the lenght of the time series (truncate to N if longer)
		IF (i > N) THEN
			ALLOCATE(self % idx(N))
			ALLOCATE(self % value(N))
			M = N
		ELSE
			ALLOCATE(self % idx(i))
			ALLOCATE(self % value(i))
			M = i
		END IF
		
		!-Read datetime index and the corresponding real value from file		
		OPEN(UNIT=1, FILE=file_name, FORM="FORMATTED", IOSTAT=res)
			DO i=1,M
				READ(UNIT=1, FMT=*, IOSTAT=res, END=100) date, time, self % value(i)
				IF (res /= 0) THEN
					print*,"Error in reading file:",file_name
					print*,"Err status:",res
					CLOSE(1)
					STOP
				END IF

				!-Fill the TIME_SERIES datetime index by converting the date_time string in input
				self % idx(i) = strptime(date//" "//time,DATETIME_FMT)
			END DO		
			100 CONTINUE
			CLOSE(1)
		print*,"Allocated space for time series of length:",M
	END SUBROUTINE load_data



	!!--- Print to screen the INPUT time series ---!!
	SUBROUTINE print_data(self)
		class(TIME_SERIES),INTENT(IN)	:: self
		LOGICAL			:: compact	!-Used to give a compact displaied series 
		INTEGER			:: i

		print*,"------------------------------------"
		print*,"  DATE",TAB,TAB,"TIME",TAB,"    VALUE"

		compact = .True.
		DO i=1,SIZE(self % idx)
					
			!-Try to print a compact view of the time series
			IF (i .gt. 5 .and. i .lt. SIZE(self % idx) - 5) THEN
				IF (compact .eqv. .True.) THEN
					print*,TAB,TAB,TAB,"    ..."
					IF (i .eq. 6) THEN 
						compact = .False.
					END IF
				END IF
			ELSE
				write(*,'(A21,2x,ES12.4)'),self % idx(i) % strftime('%Y-%m-%d %H:%M:%S'), self % value(i)
			END IF			
		END DO
		print*,"Length:",SIZE(self % idx)
		print*,"------------------------------------"
		print*,NL
	END SUBROUTINE print_data




	!!-------------------------------------------------!!
	!!----- CONVSERSION OPERATIONS ON TIME SERIES -----!!
	!!-------------------------------------------------!!

	SUBROUTINE to_timestamp(self)
		class(TIME_SERIES)	:: self
		INTEGER		:: i
		
		!-Check if the time series has already allocated a timestamp index
		IF (.not. ALLOCATED(self % timestamp)) ALLOCATE(self % timestamp(SIZE(self % idx)))
		DO i=1,SIZE(self % idx)
			self % timestamp(i) = self % idx(i) % secondsSinceEpoch()
		END DO		
	END SUBROUTINE to_timestamp


	SUBROUTINE to_datetime(self)
		class(TIME_SERIES)	:: self
		type(timedelta)	:: incr
		REAL(r9)		:: num_days,left_sec
		REAL(r9)		:: hours, minutes, seconds, milliseconds
		INTEGER		:: i
		
		!-Check if the time series has already allocated a datetime index
		IF (.not. ALLOCATED(self % idx))  ALLOCATE(self % idx(SIZE(self % timestamp)))	
		DO i=1,SIZE(self % timestamp)		
			!-Subdivide seconds in days since 1970-01-01 00:00:00 (Unix epoch)
			num_days = self % timestamp(i) / DAY_SEC
			
			!-Get the remainder and convert it to hours, minutes, seconds and milliseconds
			hours = MOD(num_days,1.)*24
			minutes = MOD(hours,1.)*60
			seconds = MOD(minutes,1.)*60
			milliseconds = 0

			!-Correct for errors
			IF ((seconds + 0.5) > 60.) THEN
				seconds = 0.
				minutes = minutes + 1.
				IF ((minutes + 0.5) > 60.) THEN
					minutes = 0.
					hours = hours + 1.
				END IF
			END IF
						
			!-Set the time delta to increment the datetime value
			incr = timedelta(0,int(hours),int(minutes),int(seconds),int(milliseconds))
			
			!-Put all together
			self % idx(i) = num2date(STANDARD_DAYS + int(num_days)) + incr
		END DO		
	END SUBROUTINE to_datetime
	



	!!------------------------------------------------!!
	!!----- ARITHMETIC OPERATIONS ON TIME SERIES -----!!
	!!------------------------------------------------!!

	!!--- Fill every element value in the series with ONES ---!!	
	SUBROUTINE fill_ones(self)
		class(TIME_SERIES)	:: self
		INTEGER		:: i,N
		
		!-Length 
		N = SIZE(self % idx)
		IF (.not. ALLOCATED(self % value)) ALLOCATE(self % value(N)) 
		
		!-Iterate over the index dimension
		DO i=1,N
			self % value(i) = 1.
		END DO
	END SUBROUTINE fill_ones



	!!--- Fill every element value in the series with ZEROS ---!!	
	SUBROUTINE fill_zeros(self)
		class(TIME_SERIES)	:: self
		INTEGER		:: i,N

		!-Length 
		N = SIZE(self % idx)
		IF (.not. ALLOCATED(self % value)) ALLOCATE(self % value(N)) 
		
		!-Iterate over the index dimension
		DO i=1,N
			self % value(i) = 0.
		END DO
	END SUBROUTINE fill_zeros	



	!!--- SUM all time series values over the time series length 
	FUNCTION sum_all(self)
		class(TIME_SERIES)	:: self
		REAL			:: sum_all
		
		!-Return the sum of all the values in the series
		sum_all = SUM(self % value)
	END FUNCTION
	
	
	!!--- FIND the INTERSECT array of datetime indices between two time series ---!!
	FUNCTION intersect(idx1, idx2)
		INTEGER,INTENT(IN)	:: idx1(:), idx2(:)
		INTEGER,ALLOCATABLE	:: ts1(:),ts2(:),intersect(:)
		INTEGER		:: i,j,m,k,len1,len2
		
		!-Minimum length between the two indices array
		len1 = SIZE(idx1)
		len2 = SIZE(idx2)
		m = MIN(len1,len2)

		!-Swap the arrays if necessary (ts1 has to be the smallest)
		IF (len1 > len2) THEN
			ts2 = idx1
			ts1 = idx2
		ELSE
			ts1 = idx1
			ts2 = idx2			
		END IF
		
		!-Allocate memory for the temporaneus result
		ALLOCATE(intersect(m))
		
		!-Find the common elements comparing every element in the smallest array with the larger one
		k=0
		DO i=1,SIZE(ts1)
			DO j=1,SIZE(ts2)
				IF (ts1(i) == ts2(j)) THEN
					k = k+1
					intersect(k) = ts1(i) 
				END IF					
			END DO
		END DO
		IF (k < m) intersect = intersect(:k) 
	END FUNCTION intersect






	!!--------------------------------------------!!
	!!----- BINARY OPERATIONS ON TIME SERIES -----!!
	!!--------------------------------------------!!
	
	!!--- SUM operator between two time series ---!!
	FUNCTION add(ts1,ts2)
		type(TIME_SERIES),INTENT(IN)	:: ts1, ts2
		type(TIME_SERIES)		:: add	
		type(TIME_SERIES)		:: inter
		INTEGER			:: i, idx1, idx2
		
		!-Convert datetime to seconds (timestamp Unix format)
		CALL ts1 % to_timestamp()
		CALL ts2 % to_timestamp()

		!-Find the common indices and set the output time series dimension
		inter % timestamp = intersect(ts1 % timestamp, ts2 % timestamp)

		!-Allocate memory for the result
		ALLOCATE(add % idx(SIZE(inter % timestamp)))
		ALLOCATE(add % value(SIZE(inter % timestamp)))
				
		!-Sum element by element (with common indices)
		DO i=1,SIZE(inter % timestamp)			
			
			!-Find the corresponding indices in the series
			idx1 = FINDLOC(ts1 % timestamp, inter % timestamp(i), DIM=1)
			idx2 = FINDLOC(ts2 % timestamp, inter % timestamp(i), DIM=1)
			
			!-ADD the corresponding values
			add % value(i) = ts1 % value(idx1) + ts2 % value(idx2)
		END DO

		!-Finally set the index to datetime values
		CALL inter % to_datetime()
		add % idx = inter % idx
	END FUNCTION add

	!!--- SUM operator between a time series and a value ---!!
	FUNCTION int_add(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		INTEGER,INTENT(IN)		:: val
		type(TIME_SERIES)		:: int_add
		INTEGER			:: i
		
		!-Allocate memory for result
		ALLOCATE(int_add % value(SIZE(ts % value)))
		
		!-Add the integer value to every element of the time series
		DO i=1,SIZE(ts % idx)
			int_add % value(i) = ts % value(i) + val
		END DO

		!-Finally set the index
		int_add % idx = ts % idx
	END FUNCTION int_add

	FUNCTION real_add(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		REAL,INTENT(IN)		:: val
		type(TIME_SERIES)		:: real_add
		INTEGER			:: i

		!-Allocate memory for result
		ALLOCATE(real_add % value(SIZE(ts % value)))
		
		!-Add the integer value to every element of the time series
		DO i=1,SIZE(ts % idx)
			real_add % value(i) = ts % value(i) + val
		END DO

		!-Finally set the index
		real_add % idx = ts % idx
	END FUNCTION real_add


	!!======================================================!!


	!!--- SUBTRACTION operator between two time series ---!!
	FUNCTION sub(ts1,ts2)
		type(TIME_SERIES),INTENT(IN)	:: ts1, ts2
		type(TIME_SERIES)		:: sub		
		type(TIME_SERIES)		:: inter
		INTEGER			:: i, idx1, idx2
		
		!-Convert datetime to seconds (timestamp Unix format)
		CALL ts1 % to_timestamp()
		CALL ts2 % to_timestamp()

		!-Find the common indices and set the output time series dimension
		inter % timestamp = intersect(ts1 % timestamp, ts2 % timestamp)

		!-Allocate memory for the result
		ALLOCATE(sub % idx(SIZE(inter % timestamp)))
		ALLOCATE(sub % value(SIZE(inter % timestamp)))
				
		!-Sum element by element (with common indices)
		DO i=1,SIZE(inter % timestamp)			
			
			!-Find the corresponding indices in the series
			idx1 = FINDLOC(ts1 % timestamp, inter % timestamp(i), DIM=1)
			idx2 = FINDLOC(ts2 % timestamp, inter % timestamp(i), DIM=1)
			
			!-ADD the corresponding values
			sub % value(i) = ts1 % value(idx1) - ts2 % value(idx2)
		END DO

		!-Finally set the index to datetime values
		CALL inter % to_datetime()
		sub % idx = inter % idx
	END FUNCTION sub

	!!--- SUBTRACTION operator between a time series and a value ---!!
	FUNCTION int_sub(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		INTEGER,INTENT(IN)		:: val
		type(TIME_SERIES)		:: int_sub
		INTEGER			:: i
		
		!-Allocate memory for result
		ALLOCATE(int_sub % value(SIZE(ts % value)))
		
		!-Add the integer value to every element of the time series
		DO i=1,SIZE(ts % idx)
			int_sub % value(i) = ts % value(i) - val
		END DO

		!-Finally set the index
		int_sub % idx = ts % idx
	END FUNCTION int_sub

	FUNCTION real_sub(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		REAL,INTENT(IN)		:: val
		type(TIME_SERIES)		:: real_sub
		INTEGER			:: i

		!-Allocate memory for result
		ALLOCATE(real_sub % value(SIZE(ts % value)))
		
		!-Add the integer value to every element of the time series
		DO i=1,SIZE(ts % idx)
			real_sub % value(i) = ts % value(i) - val
		END DO

		!-Finally set the index
		real_sub % idx = ts % idx
	END FUNCTION real_sub


	!!======================================================!!


	!!--- MOLTIPLICATION operator between two time series ---!!
	FUNCTION mult(ts1,ts2)
		type(TIME_SERIES),INTENT(IN)	:: ts1, ts2
		type(TIME_SERIES)		:: mult		
		type(TIME_SERIES)		:: inter
		INTEGER			:: i, idx1, idx2
		
		!-Convert datetime to seconds (timestamp Unix format)
		CALL ts1 % to_timestamp()
		CALL ts2 % to_timestamp()
		
		!-Find the common indices and set the output time series dimension
		inter % timestamp = intersect(ts1 % timestamp, ts2 % timestamp)

		!-Allocate memory for the result
		ALLOCATE(mult % idx(SIZE(inter % timestamp)))
		ALLOCATE(mult % timestamp(SIZE(inter % timestamp)))
		ALLOCATE(mult % value(SIZE(inter % timestamp)))
				
		!-Sum element by element (with common indices)
		DO i=1,SIZE(inter % timestamp)									
			!-Find the corresponding indices in the series
			idx1 = FINDLOC(ts1 % timestamp, inter % timestamp(i), DIM=1)
			idx2 = FINDLOC(ts2 % timestamp, inter % timestamp(i), DIM=1)
			
			!-ADD the corresponding values
			mult % timestamp(i) = inter % timestamp(i)
			mult % value(i) = ts1 % value(idx1) * ts2 % value(idx2)
		END DO

		!-Finally set the index to datetime values
		CALL mult % to_datetime()
	END FUNCTION mult

	!!--- MOLTIPLICATION operator between a time series and a value ---!!
	FUNCTION int_mult(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		INTEGER,INTENT(IN)		:: val
		type(TIME_SERIES)		:: int_mult
		INTEGER			:: i

		!-Allocate memory for result
		ALLOCATE(int_mult % value(SIZE(ts % value)))
		
		!-Add the integer value to every element of the time series
		DO i=1,SIZE(ts % idx)
			int_mult % value(i) = ts % value(i) * val
		END DO

		!-Finally set the index
		int_mult % idx = ts % idx
	END FUNCTION int_mult

	FUNCTION real_mult(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		REAL,INTENT(IN)		:: val
		type(TIME_SERIES)		:: real_mult
		INTEGER			:: i

		!-Allocate memory for result
		ALLOCATE(real_mult % value(SIZE(ts % value)))
		
		!-Add the integer value to every element of the time series
		DO i=1,SIZE(ts % idx)
			real_mult % value(i) = ts % value(i) * val
		END DO

		!-Finally set the index
		real_mult % idx = ts % idx
	END FUNCTION real_mult


	!!======================================================!!


	!!--- DIVISION operator between two time series ---!!
	FUNCTION div(ts1,ts2)
		type(TIME_SERIES),INTENT(IN)	:: ts1, ts2
		type(TIME_SERIES)		:: div		
	END FUNCTION div

	!!--- DIVISION operator between a time series and a value ---!!
	FUNCTION int_div(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		INTEGER,INTENT(IN)		:: val
		type(TIME_SERIES)		:: int_div
	END FUNCTION int_div

	FUNCTION real_div(ts,val)
		type(TIME_SERIES),INTENT(IN)	:: ts
		REAL,INTENT(IN)		:: val
		type(TIME_SERIES)		:: real_div
	END FUNCTION real_div


	!!======================================================!!

	
	!!--- EXPONENTIATION operator on a single time series ---!!
	FUNCTION int_power(ts,pot)
		type(TIME_SERIES),INTENT(IN)	:: ts
		INTEGER,INTENT(IN)		:: pot
		type(TIME_SERIES)		:: int_power
				
		!-Return the time series value to the power of "pot"
		int_power % idx = ts % idx
		int_power % value = ts % value**pot
	END FUNCTION int_power

	FUNCTION real_power(ts,pot)
		type(TIME_SERIES),INTENT(IN)	:: ts
		REAL,INTENT(IN)		:: pot
		type(TIME_SERIES)		:: real_power
				
		!-Return the time series value to the power of "pot"
		real_power % idx = ts % idx
		real_power % value = ts % value**pot
	END FUNCTION real_power


	!!--- EXPONENTIATION operator on Multi-Dimensional time series ---!!
	FUNCTION int_MD_power(ts,pot)
		type(TIME_SERIES),INTENT(IN)	:: ts(:)
		INTEGER,INTENT(IN)		:: pot
		type(TIME_SERIES)		:: int_MD_power(SIZE(ts))
		INTEGER			:: i
				
		!-Return the time series value to the power of "pot" (for every dimension)
		DO i=1,size(ts)
			int_MD_power(i) % idx = ts(i) % idx
			int_MD_power(i) % value = ts(i) % value**pot
		END DO 
	END FUNCTION int_MD_power

	FUNCTION real_MD_power(ts,pot)
		type(TIME_SERIES),INTENT(IN)	:: ts(:)
		REAL,INTENT(IN)		:: pot
		type(TIME_SERIES)		:: real_MD_power(SIZE(ts))
		INTEGER			:: i
				
		!-Return the time series value to the power of "pot" (for every dimension)
		DO i=1,size(ts)
			real_MD_power(i) % idx = ts(i) % idx
			real_MD_power(i) % value = ts(i) % value**pot
		END DO 
	END FUNCTION real_MD_power
	
END MODULE TSERIES

