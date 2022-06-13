MODULE IO
	USE MODELS
	USE PARAMS
	IMPLICIT NONE

CONTAINS

	SUBROUTINE display(model)
		type(GMIX_MODEL),INTENT(IN)	:: model
		INTEGER			:: c,i
		
		print*,"================================================================================="
		print*,"Model NAME |",TAB,TAB,"K",TAB,"    R",TAB,TAB,"S",TAB,"    Q",TAB,TAB,"P"
		print*,"---------------------------------------------------------------------------------"
		print*,"   ",model % name," |",model % K, model % R, model % S, model % Q, model % P
		print*,"---------------------------------------------------------------------------------"

		DO c=1, model % K
			write(*,'(14x)',advance='no')		
			write(*,'(A3,1x,I1)',advance='no') "k =",c
		END DO
		print*,""
		
		IF (model % S > 0) THEN
			DO i=1, model % S
				write(*,100,advance='no') "a",i,":"
				DO c=1, model % K
					write(*,'(3x,E12.4,5x,A1)',advance='no') model % a(c,i),"|"
				END DO
				print*,""
			END DO
		END IF
		
		print*,NL
		
		IF (model % R > 0) THEN
			DO i=1, model % R
				write(*,100,advance='no') "b",i,":"
				DO c=1, model % K
					write(*,'(3x,E12.4,5x,A1)',advance='no') model % b(c,i),"|"
				END DO
				print*,""
			END DO
		END IF

		print*,NL

		IF (model % Q > 0) THEN
			DO i=1, model % Q
				write(*,100,advance='no') "D",i,":"
				DO c=1, model % K
					write(*,'(3x,E12.4,5x,A1)',advance='no') model % delta(c,i),"|"
				END DO
				print*,""
			END DO
		END IF

		print*,NL

		IF (model % P > 0) THEN
			DO i=1, model % P
				write(*,100,advance='no') "B",i,":"
				DO c=1, model % K
					write(*,'(3x,E12.4,5x,A1)',advance='no') model % beta(c,i),"|"
				END DO
				print*,""
			END DO
		END IF

		print*,NL

		write(*,100,advance='no') "D0 :"
		DO c=1, model % K
			write(*,'(3x,E12.4,5x,A1)',advance='no') model % delta_0(c,1),"|"
		END DO
		print*,""

		write(*,100,advance='no') "alpha :"
		DO c=1, model % K
			write(*,'(3x,F10.4,6x,A1)',advance='no') model % alpha(c),"|"
		END DO
		
		100 FORMAT(A7,I1,1x,A1)
		print*,""
		print*,"================================================================================="
			
	END SUBROUTINE	display
END MODULE IO
