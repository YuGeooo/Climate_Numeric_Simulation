! Vorticity_equation.f95 -- by GeYq 2021.05.26

module global

    implicit none

    real, dimension(6,5):: f
    real, dimension(6,5,960):: vorticity, u, v
    real:: ds, dt
    integer:: m, n, mn

end module global


program Vorticity_equation

    !-----------------------------------

    use global
    implicit none

    real, external:: start_vorticity, time_difference_equation
    integer:: i, j, t
    !--- i for loop(n); j for loop(m)
    !--- So that loop 'j' outer ,and 'i' iner.

    !----------------------------------

    ds = 30  ! km
    dt = 0.5 ! min
    m = 5
    n = 6
    1001 format(6f8.2)

    mn = m * n 
    
    call read_array_from_txt('U.txt', u(:,:,1)) 
    call read_array_from_txt('V.txt', v(:,:,1))
    call read_array_from_txt('F.txt', f)
    
    !-----------------------------------
    write(*,*)'Start vorticity are:'
    t = 1
    do j = 1, m 
        do i = 1, n 
            vorticity(i,j,t) = start_vorticity(i,j,t)    
        end do

        !write(*,1001)u(:,j,t)
        write(*,1001)vorticity(:,j,t)
    end do
    
    call solution_of_u_v(t)
    !------------------------------------

    write(*,*)
    write(*,*)'t = 2'
    t = 2
    !u(:,:,t) = u(:,:,t)
    !v(:,:,t-1) = v(:,:,t-1)

    do j = 1, m 
        do i = 1, n 
            vorticity(i,j,t) = vorticity(i,j,t-1) - &
                 time_difference_equation(i,j,t-1)    
        end do

        write(*,1001)vorticity(:,j,t)
    end do

    call solution_of_u_v(t)
    !------------------------------------
    
    do t = 3, 960
    
        !u(:,:,t) = u(:,:,t)
        !v(:,:,t-1) = v(:,:,t-1)
        call solution_of_u_v(t)

        write(*,*)
        do j = 1, m 
            do i = 1, n 
            vorticity(i,j,t) = vorticity(i,j,t-2) - &
                2 * time_difference_equation(i,j,t-1)    
            end do

            !write(*,1001)vorticity(:,j,t)
        end do
        
    call solution_of_u_v(t)
    end do

    !------------------------------------
    

    read(*,*)

end program Vorticity_equation


subroutine solution_of_u_v(t)

    use global
    implicit none
    integer, intent(in):: t
    real, dimension(2 * mn + 1, 2 * mn):: matrix
    real, dimension(2 * mn):: x
    integer:: i, j

    call Augmented_Matrix_of_u_v(t, matrix)
    call solve_equation(matrix, 2 * mn, x)

    do j = 1, m
        do i = 1, n 
            u(i,j,t+1) = x((j-1) * n + i)
            v(i,j,t+1) = x((j-1) * n + i + mn)
        end do
    end do

end subroutine solution_of_u_v


subroutine Augmented_Matrix_of_u_v(t, matrix)

    use global
    implicit none
    integer, intent(in):: t
    real, dimension(2 * mn + 1, 2 * mn), intent(out):: matrix
    integer:: i, j, k, h

    write(*,*)'Building augmented matrix of u, v...'
    
    do j = 1, 2 * mn
        do i = 1, 2 * mn + 1
        matrix(i,j) = 0
        end do 
    end do 

    h = 1
    do j = 1, m
        do i = 1, n
            matrix(2 * mn + 1, h)    = vorticity(i,j,t) * ds
            matrix(2 * mn + 1, mn+h) =  0

            k = (j-1) * n + i 
            ! k: Num code of each point's u, v: u in 1~30, v in 31~60
                
            if(i > 1 .and. i < n .and. j > 1 .and. j < m)then   ! Centre
                
                matrix(mn + k + 1   , h) =  0.5
                matrix(mn + k - 1   , h) = -0.5
                matrix(k + n        , h) = -0.5
                matrix(k - n        , h) =  0.5 

                matrix(k + 1        , mn+h) =  0.5
                matrix(k - 1        , mn+h) = -0.5
                matrix(mn + k + n   , mn+h) =  0.5
                matrix(mn + k - n   , mn+h) = -0.5

            else if(i > 1 .and. i < n .and. j == 1)then         ! Top centre
                    
                matrix(mn + k + 1   , h) =  0.5
                matrix(mn + k - 1   , h) = -0.5
                matrix(k + n        , h) = -1 
                matrix(k            , h) =  1  

                matrix(k + 1        , mn+h) =  0.5
                matrix(k - 1        , mn+h) = -0.5
                matrix(mn + k + n   , mn+h) =  1
                matrix(mn + k       , mn+h) = -1

            else if(i > 1 .and. i < n .and. j == m)then         ! Bottom centre

                matrix(mn + k + 1   , h) =  0.5
                matrix(mn + k - 1   , h) = -0.5
                matrix(k            , h) = -1 
                matrix(k - n        , h) =  1  

                matrix(k + 1        , mn+h) =  0.5
                matrix(k - 1        , mn+h) = -0.5
                matrix(mn + k       , mn+h) =  1 
                matrix(mn + k - n   , mn+h) = -1 

            else if(i == 1 .and. j > 1 .and. j < m)then         ! Left centre
            
                matrix(mn + k + 1   , h) =  1 
                matrix(mn + k       , h) = -1 
                matrix(k + n        , h) = -0.5
                matrix(k - n        , h) =  0.5

                matrix(k + 1        , mn+h) =  1 
                matrix(k            , mn+h) = -1 
                matrix(mn + k + n   , mn+h) =  0.5
                matrix(mn + k - n   , mn+h) = -0.5

            else if(i == n .and. j > 1 .and. j < m)then         ! Right centre   

                matrix(mn + k       , h) =  1
                matrix(mn + k - 1   , h) = -1
                matrix(k + n        , h) = -0.5
                matrix(k - n        , h) =  0.5 

                matrix(k            , mn+h) =  1
                matrix(k - 1        , mn+h) = -1
                matrix(mn + k + n   , mn+h) =  0.5
                matrix(mn + k - n   , mn+h) = -0.5

            else if(i == 1 .and. j == 1)then                    ! Top left

                matrix(mn + k + 1   , h) =  1 
                matrix(mn + k       , h) = -1 
                matrix(k + n        , h) = -1 
                matrix(k            , h) =  1 

                matrix(k + 1        , mn+h) =  1
                matrix(k            , mn+h) = -1
                matrix(mn + k + n   , mn+h) =  1
                matrix(mn + k       , mn+h) = -1

            else if(i == 1 .and. j == m)then                    ! Bottom left
            
                matrix(mn + k + 1   , h) =  1
                matrix(mn + k       , h) = -1
                matrix(k            , h) = -1
                matrix(k - n        , h) =  1

                matrix(k + 1        , mn+h) =  1
                matrix(k            , mn+h) = -1
                matrix(mn + k       , mn+h) =  1
                matrix(mn + k - n   , mn+h) = -1

            else if(i == n .and. j == 1)then                    ! Top right
                        
                matrix(mn + k       , h) =  1
                matrix(mn + k - 1   , h) = -1
                matrix(k + n        , h) = -1
                matrix(k            , h) =  1

                matrix(k            , mn+h) =  1
                matrix(k - 1        , mn+h) = -1
                matrix(mn + k + n   , mn+h) =  1
                matrix(mn + k       , mn+h) = -1

            else if(i == n .and. j == m)then                    ! Bottom right
                
                matrix(mn + k       , h) =  1
                matrix(mn + k - 1   , h) = -1
                matrix(k            , h) = -1
                matrix(k - n        , h) =  1
                
                matrix(k            , mn+h) =  1
                matrix(k - 1        , mn+h) = -1
                matrix(mn + k       , mn+h) =  1
                matrix(mn + k - n   , mn+h) = -1

            end if
            h = h + 1
        end do
    end do

end subroutine Augmented_Matrix_of_u_v


subroutine solve_equation(matrix, mn, x)

    implicit none
    
    integer, intent(in):: mn
    real, dimension(mn + 1, mn), intent(in):: matrix
    real, dimension(mn + 1, mn):: a
    real, dimension(mn), intent(out):: x
    real, dimension(mn + 1):: b
    integer:: i, j, k
    real:: sum = 0
    
    a = matrix

    write(*,*)'Solving equation...'

    do k = 1, mn-1
        write(*,*)k, ' of', mn-1
        !write(*,*)a(30,:)

        ! A ranking
        do while(a(k,k) == 0)
            b(:) = a(:,k)
            do i = k, mn-1
                a(:,i) = a(:,i+1)
            end do
            a(:,mn) = b(:)
        end do

        ! Simplify 
        do i = k+1, mn
            a(:,i) = a(:,i) - a(k,i) / a(k,k) * a(:,k)
        end do
        write(*,*)a(61:,k)/30

    end do

    ! Get a solution
    x(mn) = a(mn+1,mn) / a(mn,mn)

    ! Get other solutions
    do i = mn-1, 1, -1
        do j = i+1, mn
            sum = sum + a(j,i) * x(j)
        end do
        x(i) = (a(mn+1,i) - sum) / a(i,i)
        sum = 0
    end do  

    write(*,*)x
        
end subroutine solve_equation


subroutine read_array_from_txt(filename, array)

    use global
    implicit none
    
    character(len = *), intent(in):: filename
    real, dimension(n,m), intent(out):: array
    
    open(101, file = filename, status = 'old', action = 'read')
        read(101,*)array(1:n,1:m)
    close(101)

end subroutine read_array_from_txt


function start_vorticity(i,j,t)
    
    use global
    implicit none

    integer, intent(in):: i, j, t
    real:: start_vorticity

    if(i > 1 .and. i < n .and. j > 1 .and. j < m)then   ! Centre
        start_vorticity = (( v(i+1, j, t) - v(i-1, j, t) ) / 2 &
                         - ( u(i, j+1, t) - u(i, j-1, t) ) / 2) / ds
        
    else if(i > 1 .and. i < n .and. j == 1)then         ! Top centre
        start_vorticity = (( v(i+1, j, t) - v(i-1, j, t) ) / 2 &
                         - ( u(i, j+1, t) - u(i, j, t) ) ) / ds

    else if(i > 1 .and. i < n .and. j == m)then         ! Bottom centre
        start_vorticity = (( v(i+1, j, t) - v(i-1, j, t) ) / 2 &
                           - ( u(i, j, t) - u(i, j-1, t) )) / ds

    else if(i == 1 .and. j > 1 .and. j < m)then         ! Left centre
        start_vorticity = (( v(i+1, j, t) - v(i, j, t) ) &
                         - ( u(i, j+1, t) - u(i, j-1, t) ) / 2) / ds
            
    else if(i == n .and. j > 1 .and. j < m)then         ! Right centre   
        start_vorticity = (( v(i, j, t) - v(i-1, j, t) )  &
                         - ( u(i, j+1, t) - u(i, j-1, t) ) / 2) / ds

    else if(i == 1 .and. j == 1)then                    ! Top left
        start_vorticity = (( v(i+1, j, t) - v(i, j, t) ) &
                         - ( u(i, j+1, t) - u(i, j, t) )) / ds 

    else if(i == 1 .and. j == m)then                    ! Bottom left
        start_vorticity = (( v(i+1, j, t) - v(i, j, t) ) &
                         - ( u(i, j, t) - u(i, j-1, t) )) / ds 

    else if(i == n .and. j == 1)then                    ! Top right
        start_vorticity = (( v(i, j, t) - v(i-1, j, t) ) &
                         - ( u(i, j+1, t) - u(i, j, t) )) / ds 
                        
    else if(i == n .and. j == m)then                    ! Bottom right
        start_vorticity = (( v(i, j, t) - v(i-1, j, t) ) &
                         - ( u(i, j, t) - u(i, j-1, t) )) / ds    
    end if

end function start_vorticity


function time_difference_equation(i,j,t)

    use global
    implicit none
    integer, intent(in):: i, j, t
    real:: time_difference_equation

    if(i > 1 .and. i < n .and. j > 1 .and. j < m)then   ! Centre

        time_difference_equation = dt / ds * &
        ( u(i,j,t) /2 * ( vorticity(i+1,j,t) + f(i+1,j) - vorticity(i-1,j,t) - f(i-1,j) ) &
        + v(i,j,t) /2 * ( vorticity(i,j+1,t) + f(i,j+1) - vorticity(i,j-1,t) - f(i,j-1) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i-1,j) ) /2 + ( u(i,j+1) - u(i,j-1) ) /2 ) )

    else if(i > 1 .and. i < n .and. j == 1)then         ! Top centre
    
        time_difference_equation = dt / ds * &
        ( u(i,j,t) /2 * ( vorticity(i+1,j,t) + f(i+1,j) - vorticity(i-1,j,t) - f(i-1,j) ) &
        + v(i,j,t) * ( vorticity(i,j+1,t) + f(i,j+1) - vorticity(i,j,t) - f(i,j) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i-1,j) ) /2 + ( u(i,j+1) - u(i,j) ) ) )

    else if(i > 1 .and. i < n .and. j == m)then         ! Bottom centre

        time_difference_equation = dt / ds * &
        ( u(i,j,t) /2 * ( vorticity(i+1,j,t) + f(i+1,j) - vorticity(i-1,j,t) - f(i-1,j) ) &
        + v(i,j,t) * ( vorticity(i,j,t) + f(i,j) - vorticity(i,j-1,t) - f(i,j-1) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i-1,j) ) /2 + ( u(i,j) - u(i,j-1) ) ) )

    else if(i == 1 .and. j > 1 .and. j < m)then         ! Left centre

        time_difference_equation = dt / ds * &
        ( u(i,j,t) * ( vorticity(i+1,j,t) + f(i+1,j) - vorticity(i,j,t) - f(i,j) ) &
        + v(i,j,t) /2 * ( vorticity(i,j+1,t) + f(i,j+1) - vorticity(i,j-1,t) - f(i,j-1) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i,j) ) + ( u(i,j+1) - u(i,j-1) ) /2 ) )

    else if(i == n .and. j > 1 .and. j < m)then         ! Right centre

        time_difference_equation = dt / ds * &
        ( u(i,j,t) * ( vorticity(i,j,t) + f(i,j) - vorticity(i-1,j,t) - f(i-1,j) ) &
        + v(i,j,t) /2 * ( vorticity(i,j+1,t) + f(i,j+1) - vorticity(i,j-1,t) - f(i,j-1) ) )
        !+ f(i,j) * ( ( v(i,j) - v(i-1,j) ) + ( u(i,j+1) - u(i,j-1) ) /2 ) )
        
    else if(i == 1 .and. j == 1)then                    ! Top left

        time_difference_equation = dt / ds * &
        ( u(i,j,t) * ( vorticity(i+1,j,t) + f(i+1,j) - vorticity(i,j,t) - f(i,j) ) &
        + v(i,j,t) * ( vorticity(i,j+1,t) + f(i,j+1) - vorticity(i,j,t) - f(i,j) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i,j) ) + ( u(i,j+1) - u(i,j) ) ) )
        
    else if(i == 1 .and. j == m)then                    ! Bottom left

        time_difference_equation = dt / ds * &
        ( u(i,j,t) * ( vorticity(i+1,j,t) + f(i+1,j) - vorticity(i,j,t) - f(i,j) ) &
        + v(i,j,t) * ( vorticity(i,j,t) + f(i,j) - vorticity(i,j-1,t) - f(i,j-1) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i,j) ) + ( u(i,j) - u(i,j-1) ) ) )

    else if(i == n .and. j == 1)then                    ! Top right
                        
        time_difference_equation = dt / ds * &
        ( u(i,j,t) * ( vorticity(i,j,t) + f(i,j) - vorticity(i-1,j,t) - f(i-1,j) ) &
        + v(i,j,t) * ( vorticity(i,j+1,t) + f(i,j+1) - vorticity(i,j,t) - f(i,j) ) )
        !+ f(i,j) * ( ( v(i+1,j) - v(i-1,j) ) + ( u(i,j+1) - u(i,j) ) ) )

    else if(i == n .and. j == m)then                    ! Bottom right

        time_difference_equation = dt / ds * &
        ( u(i,j,t) * ( vorticity(i,j,t) + f(i,j) - vorticity(i-1,j,t) - f(i-1,j) ) &
        + v(i,j,t) * ( vorticity(i,j,t) + f(i,j) - vorticity(i,j-1,t) - f(i,j-1) ) )
        !+ f(i,j) * ( ( v(i,j) - v(i-1,j) ) + ( u(i,j) - u(i,j-1) ) ) )

    end if

end function time_difference_equation
