 program main_program
implicit none
integer :: i, j, k, l, i1, i2, j1, j2, n, inverse_K_global
real :: q(29), x(29), y(29), member(52), joint_1(52), joint_2(52)
real :: dx, dy, dl, s, c
real, dimension(29,3) :: matrix
   character(len=20) :: filename 
real, dimension(52,3) :: matrix_2
   character(len=20) :: filename_2
real, dimension(58,58) :: K_global
real, dimension(52,4,4) :: K_member 
real,dimension(:,:),allocatable::A
real :: load_group_pos, load_group_magnitude 
real :: joint_forces(29) 
  ! Open the file for reading
   open(unit=10,file="coordinates.txt",status='old')

   ! Read in the matrix values from the file
   do i = 1, 29
         read(10,*) (matrix(i,j), j=1,3)
   end do

   close(10)
   
   ! Print the matrix
   do i = 1, 29
         write(*,*) (matrix(i,j), j=1,3)
   end do
   
   
    !saving values from matrices where the first column is q(i)=joint number
	!x(i) = coordinate in x direction
	!y(i) = coordinate in y direction
    do i = 1, 29
    q(i) = matrix(i, 1)
    x(i) = matrix(i, 2)
    y(i) = matrix(i, 3)
   end do
   
    ! Open the file for reading
   open(unit=10,file="truss_members.txt",status='old')

   ! Read in the matrix values from the file
   do k = 1, 52
         read(10,*) (matrix_2(k,l), l=1,3)
   end do
 
  
   close(10)
   
   ! Print the matrix
   do k = 1, 52
         write(*,*) (matrix_2(k,l), l=1,3) 
   end do
   
  
    !saving values from matrices where the first column is member(k)=member number
	!joint_1(k) = first joint  number from second column 
	!joint_2(k) = second joint number from third column
   do k = 1, 52
    member(k) = matrix_2(k, 1)
    joint_1(k) = matrix_2(k, 2)
    joint_2(k) = matrix_2(k, 3)
   end do
   
   ! calculating the length of each member
   do k = 1, 52
    dx = x(int(joint_2(k))) - x(int(joint_1(k)))
    dy = y(int(joint_2(k))) - y(int(joint_1(k)))
    dl = sqrt(dx**2 + dy**2)
    member(k) = dl
	s=dy/dl
	c=dx/dl
	K_member(k,1,1) = c**2
    K_member(k,1,2) = c*s
    K_member(k,1,3) = -c**2
    K_member(k,1,4) = -c*s
    K_member(k,2,1) = c*s
    K_member(k,2,2) = s**2
    K_member(k,2,3) = -c*s
    K_member(k,2,4) = -s**2
    K_member(k,3,1) = -c**2
    K_member(k,3,2) = -c*s
    K_member(k,3,3) = c**2
    K_member(k,3,4) = c*s
    K_member(k,4,1) = -c*s
    K_member(k,4,2) = -s**2
    K_member(k,4,3) = c*s
    K_member(k,4,4) = s**2
	    write(*,*) "Length of member ", k, " is ", dl
        write(*,*) "sin of member    ", k, " is ", s
        write(*,*) "cos of member    ", k, " is ", c		
	do i=1,4
		write(*,*) "stiffness matrix    ", k, " is ", (K_member(k,i,j),j=1,4)
	  
    		
	end do
	
   end do



    K_global=0
    ! Assemble the global stiffness matrix
    do k = 1, 52
        i1 = 2*(int(joint_1(k))-1)+1
        i2 = 2*(int(joint_2(k))-1)+1
        j1 = i1+1
        j2 = i2+1
        K_global(i1:i2,j1:j2) = K_global(i1:i2,j1:j2) + K_member(k,:,:)
        K_global(j1:j2,i1:i2) = K_global(j1:j2,i1:i2) + K_member(k,:,:)
    end do

    ! Print the global stiffness matrix
    do i = 1, 58
        write(*,*) (K_global(i,j), j=1, 58)
    end do
	
	
	
	
	
	! ! Open the file
    ! open(unit=10,file="load_group.txt",status='old')

    ! ! ! Get the number of rows in the matrix
     ! read(10,*) n

    ! ! ! Allocate memory for the matrix
    ! allocate(A(n,2))

    ! ! ! Read in the matrix values from the file
    ! do i = 1, n
         ! read(10,*) (A(i,j), j=1,2)
    ! end do

    ! ! Close the file
    ! close(10)
	
    ! Print the matrix
	
    ! do i = 1, n
        ! write(*,*) (A(i,j), j=1,2)
    ! end do
	
	open(unit=10,file="load_group.txt",status='old')

    ! Read in the matrix values from the file
    do i = 1, 5
         read(10,*) (matrix(i,j), j=1,2)
    end do

    close(10)
	
    joint_forces=0

    ! Print the matrix
	
    do i = 1, 5
        write(*,*) (matrix(i,j), j=1,2)
    end do
	do i = 1, 5
    load_group_pos = matrix(i, 1)
    load_group_magnitude = matrix(i, 2)
   
	
	    ! Apply the load group to the appropriate member
        do k = 1, 52
	
            ! Check if the member is a horizontal member
	
            if (y(int(joint_2(k)))==0 .and. y(int(joint_1(k)))==0) then
		
            ! Check if the load group position is within the x-coordinates of the member's joints
			
                if (x(int(joint_1(k)))== load_group_pos) then 	
			
                    joint_forces(int(joint_1(k))) = joint_forces(int(joint_1(k)))+ load_group_magnitude
				
                else if (x(int(joint_2(k)))== load_group_pos) then 	
			
                    joint_forces(int(joint_2(k))) = joint_forces(int(joint_2(k)))+ load_group_magnitude
					
					! Calculate the joint force for each joint
			    else if (x(int(joint_2(k)))>load_group_pos .and. load_group_pos>  x(int(joint_1(k)))) then
                    joint_forces(int(joint_2(k))) = joint_forces(int(joint_2(k)))+&
					load_group_magnitude*(load_group_pos-x(int(joint_1(k))))/member(k)
                    joint_forces(int(joint_1(k))) = joint_forces(int(joint_1(k)))+ &
					load_group_magnitude*(x(int(joint_2(k)))-load_group_pos)/member(k)
                end if
				
            end if
        end do
	end do
	
	do i = 1, 29
        write(*,*) "joint force for joint ",i," is : ",joint_forces(i)
    end do
    
   
	
end program main_program