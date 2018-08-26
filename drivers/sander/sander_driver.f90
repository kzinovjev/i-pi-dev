! I-PI driver for sander (AmberTools)
! Based on drivers/driver.f90 driver for test potentials
! Kirill Zinovjev, 08.2018
! Usage: ./sander_driver input_file prmtop inpcrd
program sander_driver

    use f90sockets, only : open_socket, writebuffer, readbuffer

    use sander_api, only: sander_input, pme_sander_input, &
                          sander_setup, potential_energy_rec, energy_forces, &
                          qmmm_input_options, qm_sander_input, &
                          read_inpcrd_file, get_inpcrd_natom, &
                          set_positions, set_box

    implicit none

    integer, parameter :: MSGLEN = 12
    double precision, parameter :: BOHR_TO_ANGSTROM = 0.529177249
    double precision, parameter :: KCALMOL_TO_HARTREE = 1.593601D-3

    integer socket, inet, port, ierr, i
    character(len=1024) :: host

    ! Socket communication
    character(len=MSGLEN) :: header
    logical :: hasdata = .false.
    integer :: ipi_natom

    ! Sander structures
    type(potential_energy_rec) :: energies
    type(sander_input)         :: options
    type(qmmm_input_options)   :: qmmm_options

    ! Input files
    character(len=1024) :: input_file, inpcrd, prmtop

    ! System
    integer :: natom
    double precision :: energy
    double precision, dimension(6) :: box
    double precision, dimension(:), allocatable :: coordinates, forces

    double precision, dimension(9) :: virial = 0 ! Just a hardcoded 0 vector

    ! Cell
    double precision :: a, b, c
    double precision :: alpha = 90., beta = 90., gamma = 90. ! Assumed to be orthorhombic
    double precision, dimension(9) :: cell

    ! TODO: hardcoded for now
    inet = 0
    port = 21142
    host = "localhost"
    input_file = "in"
    inpcrd = "../../ch3cl/react.rst"
    prmtop = "../../ch3cl/ch3cl.prmtop"

    ! Prepare sander options
    call pme_sander_input(options)
    options%cut = 8.d0
    options%ifqnt=1

    ! Prepare qmmm options
    call qm_sander_input(qmmm_options)
    qmmm_options%qmmask = ':1-2'
    qmmm_options%qmcharge = -1
    qmmm_options%qm_theory = 'AM1'
    qmmm_options%qmcut = 8.d0

    ! Read coordinates and parameters
    call get_inpcrd_natom(INPCRD, natom)
    allocate(coordinates(3*natom), forces(3*natom))
    call read_inpcrd_file(INPCRD, coordinates, box, ierr)
    call sander_setup(PRMTOP, coordinates, box, options, qmmm_options, ierr)

    call open_socket(socket, inet, port, trim(host)//achar(0))

    do while (.true.)

        call readbuffer(socket, header, MSGLEN)

        select case (header)

            ! The wrapper is inquiring on what we are doing
            case ("STATUS")
                if (hasdata) then
                    ! Signals that we are done computing and can return forces
                    call writebuffer(socket, "HAVEDATA    ", MSGLEN)
                else
                    ! We are idling and eager to compute something
                    CALL writebuffer(socket, "READY       ", MSGLEN)
                end if

            ! The driver is sending the positions of the atoms. Here is where we do the calculation!
            case ("POSDATA")
                call readbuffer(socket, cell, 9)
                !Check that all off-diagonal elements are 0
                if (maxval(abs(cell), mask=(/ (mod(i, 4) /= 1, i = 1, 9) /)) > 1D-10 ) then
                    write(*,*) " Box not orthorhombic ", cell
                    STOP "ENDED"
                end if
                ! Take a, b, c from the diagonal
                a = cell(1) * BOHR_TO_ANGSTROM
                b = cell(5) * BOHR_TO_ANGSTROM
                c = cell(9) * BOHR_TO_ANGSTROM
                call readbuffer(socket, cell, 9) ! Inverse is not needed, read and forget

                ! Read number of atoms from the wrapper and compare
                call readbuffer(socket, ipi_natom)
                if (ipi_natom /= natom) then
                    write(*,*) " Inconsistent number of atoms ", natom, ipi_natom
                    STOP "ENDED"
                end if

                call readbuffer(socket, coordinates, 3*natom)
                call set_positions(coordinates * BOHR_TO_ANGSTROM)
                call set_box(a, b, c, alpha, beta, gamma)
                call energy_forces(energies, forces)
                hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper

            ! The driver calculation is finished, it's time to send the results back to the wrapper
            case ("GETFORCE")
                call writebuffer(socket, "FORCEREADY  ", MSGLEN)
                call writebuffer(socket, energies%tot * KCALMOL_TO_HARTREE)
                call writebuffer(socket, natom)
                call writebuffer(socket, forces * KCALMOL_TO_HARTREE * BOHR_TO_ANGSTROM, 3*natom)
                call writebuffer(socket, virial, 9)
                call writebuffer(socket, 0)
                hasdata = .false.

            case default
                write(*,*) " Unexpected header ", header
                stop "ENDED"

        end select

    end do

end program sander_driver