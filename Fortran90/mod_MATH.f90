MODULE MATH
	USE, INTRINSIC :: IEEE_ARITHMETIC
	IMPLICIT NONE
	PRIVATE
	
	PUBLIC :: inverse
	
CONTAINS

	FUNCTION inverse(A) RESULT(I)
		REAL, INTENT(IN)	 	:: A(:,:)				!-INPUT matrix
		REAL			:: I(SIZE(A,DIM=1),SIZE(A,DIM=1))	!-Inverse matrix
		REAL			:: detInv

		!-Set Inverse matrix to zero value in case of return error
		I = 0
		IF (SIZE(A,DIM=1) == SIZE(A,DIM=2)) THEN 

			!-Direct calculation of the inverse 2 x 2 matrix
			IF (SIZE(A,DIM=1) == 2) THEN
				!-Calculate the inverse determinant 
				detInv = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

				!-Check whether the matrix has finite != 0 determinant value
				IF ((detInv == 0.) .or. (.not. IEEE_IS_FINITE(detInv))) THEN
					print*,"Matrix NOT invertible !"
					RETURN
				END IF
		
				!-Calculate the inverse matrix
				I(1,1) = +detInv * A(2,2)
				I(2,1) = -detInv * A(2,1)
				I(1,2) = -detInv * A(1,2)
				I(2,2) = +detInv * A(1,1)
			END IF

			
			!-Direct calculation of the inverse 3 x 3 matrix
			IF (SIZE(A,DIM=1) == 3) THEN
				!-Calculate the inverse determinant of the matrix
				detInv = 1./(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
				      	 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
				         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

				!-Check whether the matrix has finite != 0 determinant value
				IF ((detInv == 0.) .or. (.not. IEEE_IS_FINITE(detInv))) THEN
					print*,"Matrix NOT invertible !"
					RETURN
				END IF
							
				!-Calculate the inverse of the matrix
				I(1,1) = +detInv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
				I(2,1) = -detInv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
				I(3,1) = +detInv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
				I(1,2) = -detInv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
				I(2,2) = +detInv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
				I(3,2) = -detInv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
				I(1,3) = +detInv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
				I(2,3) = -detInv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
				I(3,3) = +detInv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

			END IF


			!-Direct calculation of the inverse 4 x 4 matrix
			IF (SIZE(A,DIM=1) == 4) THEN

				!-Calculate the inverse determinant of the matrix
				detInv = &
				1./(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
				- A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
				+ A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
				- A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

				!-Check whether the matrix has finite != 0 determinant value
				IF ((detInv == 0.) .or. (.not. IEEE_IS_FINITE(detInv))) THEN
					print*,"Matrix NOT invertible !"
					RETURN
				END IF
				 	
				! Calculate the inverse of the matrix
				I(1,1) = detInv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
				I(2,1) = detInv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
				I(3,1) = detInv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
				I(4,1) = detInv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
				I(1,2) = detInv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
				I(2,2) = detInv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
				I(3,2) = detInv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
				I(4,2) = detInv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
				I(1,3) = detInv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
				I(2,3) = detInv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
				I(3,3) = detInv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
				I(4,3) = detInv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
				I(1,4) = detInv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
				I(2,4) = detInv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
				I(3,4) = detInv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
				I(4,4) = detInv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
			END IF
		ELSE
			print*,"Input matrix is not squared:",SIZE(A,DIM=1),"x",SIZE(A,DIM=2)
			RETURN
		END IF			
	END FUNCTION inverse
	
END MODULE MATH

