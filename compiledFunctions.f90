!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! HELIUM !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine helium_test_function(N, r, alpha, psi)
	implicit none
	INTEGER(8), INTENT(IN)				 :: N
	REAL(8), DIMENSION(N,6), INTENT(IN)  :: r
	REAL(8),				 INTENT(IN)  :: alpha

	REAL(8), DIMENSION(N)			   	 :: r1
	REAL(8), DIMENSION(N)			     :: r2
	REAL(8), DIMENSION(N) 			     :: r12

	REAL(8), DIMENSION(N), 	 INTENT(OUT) :: psi

	r1  = sqrt(sum(r(:,1:3)**2,2))
	r2  = sqrt(sum(r(:,4:6)**2,2))
	r12 = sqrt(sum((r(:,1:3)-r(:,4:6))*(r(:,1:3)-r(:,4:6)),2))

	psi = exp(-2*(r1+r2) + r12/(2*(1 + alpha*r12)))
end subroutine helium_test_function




subroutine helium_local_energy(N, r, alpha, E)
	implicit none
	INTEGER(8), INTENT(IN)				 :: N
	REAL(8), DIMENSION(N,6), INTENT(IN)  :: r
	REAL(8),				 INTENT(IN)  :: alpha

	REAL(8), DIMENSION(N)			   	 :: r1
	REAL(8), DIMENSION(N)			     :: r2
	REAL(8), DIMENSION(N) 			     :: r12
	REAL(8), DIMENSION(N) 			     :: r12special

	REAL(8), DIMENSION(N), 	 INTENT(OUT) :: E

	r1  = sqrt(sum(r(:,1:3)**2,2))
	r2  = sqrt(sum(r(:,4:6)**2,2))
	r12 = sqrt(sum((r(:,1:3)-r(:,4:6))*(r(:,1:3)-r(:,4:6)),2))
	r12special = sum((r(:,1:3)/spread(r1,2,3) - r(:,4:6)/spread(r2,2,3))* &
		(r(:,1:3)-r(:,4:6)),2)/r12
	E = (-4 + r12special/(1 + alpha*r12)**2 - 1/(r12*(1 + alpha*r12)**3) -&
	 1/(4*(1 + alpha*r12)**4) + 1/r12)
end subroutine helium_local_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! HYDROGEN !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hydrogen_local_energy(N, v_r, beta, s, a, E, steep_descent_beta)
	implicit none
	!INPUT
	INTEGER(8), INTENT(IN)				 :: N
	REAL(8), DIMENSION(N,6), INTENT(IN)  :: v_r
	REAL(8),				 INTENT(IN)  :: beta
	REAL(8),				 INTENT(IN)  :: s
	REAL(8),				 INTENT(IN)  :: a
	!OUTPUT
	REAL(8), DIMENSION(N), 	 INTENT(OUT) :: E
	REAL(8), DIMENSION(N), 	 INTENT(OUT) :: steep_descent_beta	
	!DUMMYS

	!POSITIONS
	REAL(8), DIMENSION(N,3)			   	 :: v_r1
	REAL(8), DIMENSION(N,3)			   	 :: v_r2
	REAL(8), DIMENSION(N,3)			   	 :: v_h

	REAL(8), DIMENSION(N,3)			   	 :: v_r1L
	REAL(8), DIMENSION(N,3)			   	 :: v_r1R
	REAL(8), DIMENSION(N,3)			     :: v_r2L
	REAL(8), DIMENSION(N,3)			   	 :: v_r2R
	REAL(8), DIMENSION(N,3)			     :: v_r12
	REAL(8), DIMENSION(N,3)			     :: v_special

	!DISTANCES
	REAL(8), DIMENSION(N)			   	 :: r1
	REAL(8), DIMENSION(N)	 		   	 :: r2

	REAL(8), DIMENSION(N)			   	 :: r1L
	REAL(8), DIMENSION(N)			   	 :: r1R
	REAL(8), DIMENSION(N)			     :: r2L
	REAL(8), DIMENSION(N)			   	 :: r2R
	REAL(8), DIMENSION(N) 			     :: r12
	REAL(8), DIMENSION(N) 			     :: special

	!INTERMEDIATE EXPRESSIONS
	REAL(8), DIMENSION(N)			   	 :: phi1
	REAL(8), DIMENSION(N)			   	 :: phi2
	REAL(8), DIMENSION(N)			   	 :: phi1L
	REAL(8), DIMENSION(N)			   	 :: phi1R
	REAL(8), DIMENSION(N)			   	 :: phi2L
	REAL(8), DIMENSION(N)			   	 :: phi2R


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	v_r1 = v_r(:, 1:3)
	v_r2 = v_r(:, 4:6)
	
	v_h  = spread((/ s/2.0_8, 0.0_8, 0.0_8 /),1,N)
	v_r1L = v_r1 + v_h
	v_r1R = v_r1 - v_h
	v_r2L = v_r2 + v_h
	v_r2R = v_r2 - v_h
	v_r12 = v_r1 - v_r2

	r1  = sqrt(sum(v_r1**2,2))
	r2  = sqrt(sum(v_r2**2,2))

	r1L = sqrt(sum(v_r1L**2,2))
	r1R = sqrt(sum(v_r1R**2,2))
	r2L = sqrt(sum(v_r2L**2,2))
	r2R = sqrt(sum(v_r2R**2,2))
	r12 = sqrt(sum(v_r12**2,2))

	phi1L = exp(-r1L/a)
	phi1R = exp(-r1R/a)
	phi2L = exp(-r2L/a)
	phi2R = exp(-r2R/a)
	phi1  = phi1L + phi1R
	phi2  = phi2L + phi2R

	v_special = spread(phi1L/r1L/phi1,2,3) *v_r1L + spread(phi1R/r1R/phi1,2,3) *v_r1R &
		- spread(phi2L/r2L/phi2,2,3) *v_r2L - spread(phi2R/r2R/phi2,2,3) *v_r2R
	special = sum(v_special*v_r12,2)/(r12*2*a*(1+beta*r12)**2)

	E = -1/a**2 + &
		 1/(a*phi1) * (phi1L/r1L + phi1R/r1R) + &
		 1/(a*phi2) * (phi2L/r2L + phi2R/r2R) + &
		-(1/r1L + 1/r1R + 1/r2L + 1/r2R) + &
		 1/r12 + &
		 special + &
		-((4*beta + 1)*r12 + 4)/(4*(1+beta*r12)**4 * r12)

	steep_descent_beta =  - r12**2 /(2 * (1 + beta*r12)**2 ) &
	* phi1 * phi2 * exp(r12/(2* (1 + beta*r12)))
	end subroutine hydrogen_local_energy


subroutine hydrogen_test_function(N, v_r, beta, s, a, psi)
	implicit none
	!INPUT
	INTEGER(8), INTENT(IN)				 :: N
	REAL(8), DIMENSION(N,6), INTENT(IN)  :: v_r
	REAL(8),				 INTENT(IN)  :: beta
	REAL(8),				 INTENT(IN)  :: s
	REAL(8),				 INTENT(IN)  :: a
	!OUTPUT
	REAL(8), DIMENSION(N), 	 INTENT(OUT) :: psi
	
	!DUMMYS

	!POSITIONS
	REAL(8), DIMENSION(N,3)			   	 :: v_r1
	REAL(8), DIMENSION(N,3)			   	 :: v_r2
	REAL(8), DIMENSION(N,3)			   	 :: v_h

	REAL(8), DIMENSION(N,3)			   	 :: v_r1L
	REAL(8), DIMENSION(N,3)			   	 :: v_r1R
	REAL(8), DIMENSION(N,3)			     :: v_r2L
	REAL(8), DIMENSION(N,3)			   	 :: v_r2R
	REAL(8), DIMENSION(N,3)			     :: v_r12

	!DISTANCES
	REAL(8), DIMENSION(N)			   	 :: r1
	REAL(8), DIMENSION(N)	 		   	 :: r2

	REAL(8), DIMENSION(N)			   	 :: r1L
	REAL(8), DIMENSION(N)			   	 :: r1R
	REAL(8), DIMENSION(N)			     :: r2L
	REAL(8), DIMENSION(N)			   	 :: r2R
	REAL(8), DIMENSION(N) 			     :: r12

	!INTERMEDIATE EXPRESSIONS
	REAL(8), DIMENSION(N)			   	 :: phi1L
	REAL(8), DIMENSION(N)			   	 :: phi1R
	REAL(8), DIMENSION(N)			   	 :: phi2L
	REAL(8), DIMENSION(N)			   	 :: phi2R


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	v_r1 = v_r(:, 1:3)
	v_r2 = v_r(:, 4:6)
	
	v_h  = spread((/ s/2.0_8, 0.0_8, 0.0_8 /),1,N)
	v_r1L = v_r1 + v_h
	v_r1R = v_r1 - v_h
	v_r2L = v_r2 + v_h
	v_r2R = v_r2 - v_h
	v_r12 = v_r1 - v_r2

	r1  = sqrt(sum(v_r1**2,2))
	r2  = sqrt(sum(v_r2**2,2))

	r1L = sqrt(sum(v_r1L**2,2))
	r1R = sqrt(sum(v_r1R**2,2))
	r2L = sqrt(sum(v_r2L**2,2))
	r2R = sqrt(sum(v_r2R**2,2))
	r12 = sqrt(sum(v_r12**2,2))

	phi1L = exp(-r1L/a)
	phi1R = exp(-r1R/a)
	phi2L = exp(-r2L/a)
	phi2R = exp(-r2R/a)


	psi = (phi1L + phi1R)* (phi2L + phi2R)* &
		exp(r12/(2* (1 + beta*r12)))
	end subroutine hydrogen_test_function
