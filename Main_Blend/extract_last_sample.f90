program extract_last_sample
    implicit none
    integer, parameter :: Nx = 16, Ny = 16, Nz = 16
    integer, parameter :: N = 1000
    

  ! Total number of samples
    real :: col1, col2,col3,col4
    integer :: i, j, k, n1
    integer :: io_status
    character(len=100) :: input_filename, output_filename

    ! Read command-line arguments for input and output file names
    call get_command_argument(1, input_filename)
    call get_command_argument(2, output_filename)

    ! Open the input file (assumed formatted file with column data)
    open(unit=10, file=trim(input_filename), status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
        print *, 'Error opening input file: ', trim(input_filename)
        stop
    end if

    ! Skip all samples except the last one
    do n1 = 1, N - 1
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx
                    ! Skip reading the data for previous samples
                    read(10, *, iostat=io_status) col1, col2!,col3,col4
                    if (io_status /= 0) then
                        print *, 'Error reading data from file!'
                        stop
                    end if
                end do
            end do
        end do
    end do

    ! Open the output file to write the last sample (Nth sample)
    open(unit=20, file=trim(output_filename), status='replace', action='write', iostat=io_status)
    if (io_status /= 0) then
        print *, 'Error opening output file: ', trim(output_filename)
        stop
    end if

    ! Now read and write the last sample (Nth sample)
    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                ! Read data for the last sample
                read(10, *, iostat=io_status) col1, col2!,col3,col4
                if (io_status /= 0) then
                    print *, 'Error reading data from file!'
                    stop
                end if
                ! Write the last sample to output file
                write(20, *) col1, col2
            end do
        end do
    end do

    ! Close the input and output files
    close(10)
    close(20)

    print *, 'Last sample from ', trim(input_filename), ' has been successfully written to ', trim(output_filename)

end program extract_last_sample

