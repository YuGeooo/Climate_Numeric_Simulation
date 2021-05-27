! Map_projection.f95 -- by GeYq 2021.05.24

module global

    implicit none 

    integer:: m, n    

end module global

program Map_projection

    use global
    implicit none

    real:: lat1, lat2, lon1, lon2, degree
    real, dimension(99,99):: lanti, x
    integer:: i, j

    !-------------- The extent of lattice point area

    lat1   = 30
    lat2   = 75
    lon1   = 85
    lon2   = 155
    degree = 2.5
    
    !--- long / lat dimension of latttice point area

    m = (lon2 - lon1) / degree + 1
    n = (lat2 - lat1) / degree + 1

    !-------------- Lantitudes of each lattice point

    write(*,*)'Calculate lantitudes...'
    
    do i = 1, n
        lanti(1:m, i) = lat2 - (i-1) * degree
    end do

    !--------------------------- Calculate and output

    call map_projection_coefficient(lanti, x)
    call save_array_as_csv(x,'Map_projiction.csv')

    write(*,*)'Program finished!'
    read(*,*)

end program Map_projection

!----------------------- A subroutine that calculate map projection coefficient

subroutine map_projection_coefficient(lantitudes, mpc)

    ! lantitudes : lantitudes of each lattice point
    ! mpc        : result(output), map projection coefficient

    use global
    implicit none
    
    real, dimension(99,99), intent(in):: lantitudes
    real, dimension(99,99), intent(out):: mpc
    real:: pi = 3.1415926
    integer:: i, j

    write(*,*)'Calculating map projection coefficient of each lattice point...'

    do i = 1, n
        do j = 1, m

            mpc(j,i) = 1 / (2 * sin( lantitudes(j,i) * pi/180 ) ) &
             * ( tan( lantitudes(j,i) / 2 * pi/180 ) / tan(15 * pi/180) )**0.7156
        
        end do
    end do

    return 
end subroutine map_projection_coefficient

!--------------------- A subroutine that save array as CSV file

subroutine save_array_as_csv(array, filename)

    ! array    : The array we need to output
    ! filename : The name of outputing CSV file
    
    use global
    implicit none
    real, dimension(99,99), intent(in):: array
    character(len = *), intent(in):: filename
    integer:: i, j

    write(*,*)'Writing results as CSV file...'
    
    open(101, file = filename, status = 'replace', action = 'write')  
        do i = 1, n
            do j = 1, m
                write(101, 1001, advance = 'no')array(j,i), ","
            end do
            write(101,*)
        end do
    close(101)

    1001 format(f8.3,a)

end subroutine save_array_as_csv