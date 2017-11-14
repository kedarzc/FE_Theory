program jac

	implicit none

	integer :: i
	real :: eta, xi, detJac
	real, dimension(4) :: x,y
	real :: coordinates(4,2), A(2,4), jacobian(2,2)

	x = (/0,2,2,0/)
	y = (/0,0,4,4/)

	do i=1,4
		coordinates(i,1) = x(i)
		coordinates(i,2) = y(i)
	enddo

	A(1,1) = -(1-eta)
	A(1,2) = +(1-eta)
	A(1,3) = +(1+eta)
	A(1,4) = -(1+eta)

	A(2,1) = -(1-xi)
	A(2,2) = -(1+xi)
	A(2,3) = +(1+xi)
	A(2,4) = +(1-xi)

	jacobian = (1.0/4.0)*matmul(A,coordinates)

	detJac = jacobian(1,1)*jacobian(1,1) - jacobian(1,2)*jacobian(1,2)


	write(*,*) jacobian(:,1)
	write(*,*) jacobian(:,2)
	write(*,*) ''
	write(*,*) detJac

end program jac