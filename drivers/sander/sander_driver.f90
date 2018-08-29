! I-PI driver for sander (AmberTools)
! Based on drivers/driver.f90 driver for test potentials
! Kirill Zinovjev, 08.2018
! Usage: sander_driver input_file prmtop inpcrd host port inet
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
    character(len=1024) :: tmp
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

    if (command_argument_count() < 6) then
        write(*,*) "Usage: sander_driver input_file prmtop inpcrd host port inet"
        stop
    end if

    call get_command_argument(1, input_file)
    call get_command_argument(2, prmtop)
    call get_command_argument(3, inpcrd)
    call get_command_argument(4, host)
    call get_command_argument(5, tmp)
    read(tmp,*) port
    call get_command_argument(6, tmp)
    read(tmp,*) inet

    call read_sander_input(input_file, options, qmmm_options)

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
                    stop "ENDED"
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
                    stop "ENDED"
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

contains

    subroutine read_sander_input (filename, options, qmmm_options)

        character(len=*), intent(in)         :: filename
        type(sander_input), intent(out)       :: options
        type(qmmm_input_options), intent(out) :: qmmm_options

        integer, parameter :: MAX_QUANTUM_ATOMS = size(qmmm_options%iqmatoms)

        ! sander_input variables
        double precision :: extdiel, intdiel, rgbmax, saltcon, cut, &
                dielc, rdt, fswitch, restraint_wt
        integer :: igb, alpb, gbsa, lj1264, ipb, inp, vdwmeth, ew_type, ntb, ifqnt, &
                jfastw, ntf, ntc, ntr, ibelly, mask_from_ref
        character(len=256) :: restraintmask, bellymask, refc

        ! qmmm_input_options variables
        double precision :: qmcut, lnk_dis, scfconv, errconv, dftb_telec, dftb_telec_step, &
                fockp_d1, fockp_d2, fockp_d3, fockp_d4, damp, vshift, kappa, &
                pseudo_diag_criteria, min_heavy_mass, r_switch_hi, r_switch_lo
        integer :: iqmatoms(MAX_QUANTUM_ATOMS), qmgb, lnk_atomic_no, &
                ndiis_matrices, ndiis_attempts, lnk_method, qmcharge, &
                corecharge, buffercharge, spin, qmqmdx, verbosity, &
                printcharges, printdipole, print_eigenvalues, peptide_corr, &
                itrmax, printbondorders, qmshake, qmmmrij_incore, &
                qmqm_erep_incore, pseudo_diag, qm_ewald, qm_pme, kmaxqx, &
                kmaxqy, kmaxqz, ksqmaxq, qmmm_int, adjust_q, tight_p_conv, &
                diag_routine, density_predict, fock_predict, vsolv, &
                dftb_maxiter, dftb_disper, dftb_chg, abfqmmm, hot_spot, &
                qmmm_switch, core_iqmatoms(MAX_QUANTUM_ATOMS), &
                buffer_iqmatoms(MAX_QUANTUM_ATOMS)
        character(len=8192) :: qmmask, coremask, buffermask, centermask
        character(len=256) :: dftb_3rd_order, dftb_slko_path
        character(len=12) :: qm_theory

        namelist /cntrl/ extdiel, intdiel, rgbmax, saltcon, cut, dielc, rdt, fswitch, &
                         restraint_wt,igb, alpb, gbsa, lj1264, ipb, inp, vdwmeth, ew_type, &
                         ntb, ifqnt, jfastw, ntf, ntc, ntr, ibelly, mask_from_ref

        namelist /qmmm/ qmcut, lnk_dis, scfconv, errconv, dftb_telec, dftb_telec_step, &
                fockp_d1, fockp_d2, fockp_d3, fockp_d4, damp, vshift, kappa, &
                pseudo_diag_criteria, min_heavy_mass, r_switch_hi, r_switch_lo, &
                iqmatoms, qmgb, lnk_atomic_no, ndiis_matrices, ndiis_attempts, lnk_method, &
                qmcharge, corecharge, buffercharge, spin, qmqmdx, verbosity, &
                printcharges, printdipole, print_eigenvalues, peptide_corr, &
                itrmax, printbondorders, qmshake, qmmmrij_incore, &
                qmqm_erep_incore, pseudo_diag, qm_ewald, qm_pme, kmaxqx, &
                kmaxqy, kmaxqz, ksqmaxq, qmmm_int, adjust_q, tight_p_conv, &
                diag_routine, density_predict, fock_predict, vsolv, &
                dftb_maxiter, dftb_disper, dftb_chg, abfqmmm, hot_spot, &
                qmmm_switch, core_iqmatoms, buffer_iqmatoms, qmmask, coremask, &
                buffermask, centermask, dftb_3rd_order, dftb_slko_path, qm_theory


        call pme_sander_input(options)

        igb = options%igb
        alpb = options%alpb
        gbsa = options%gbsa
        lj1264 = options%lj1264
        ipb = options%ipb
        inp = options%inp
        vdwmeth = options%vdwmeth
        ew_type = options%ew_type
        extdiel = options%extdiel
        intdiel = options%intdiel
        rgbmax = options%rgbmax
        saltcon = options%saltcon
        cut = options%cut
        dielc = options%dielc
        ifqnt = options%ifqnt
        jfastw = options%jfastw
        ntf = options%ntf
        ntc = options%ntc
        fswitch = options%fswitch
        ntr = options%ntr
        ibelly = options%ibelly
        restraint_wt = options%restraint_wt
        restraintmask = options%restraintmask
        bellymask = options%bellymask
        refc = options%refc
        mask_from_ref = options%mask_from_ref


        call qm_sander_input(qmmm_options)

        qmcut = qmmm_options%qmcut
        lnk_dis = qmmm_options%lnk_dis
        lnk_atomic_no = qmmm_options%lnk_atomic_no
        lnk_method = qmmm_options%lnk_method
        qmgb = qmmm_options%qmgb
        qm_theory = qmmm_options%qm_theory
        qmcharge = qmmm_options%qmcharge
        corecharge = qmmm_options%corecharge
        buffercharge = qmmm_options%buffercharge
        spin = qmmm_options%spin
        qmqmdx = qmmm_options%qmqmdx
        verbosity = qmmm_options%verbosity
        scfconv = qmmm_options%scfconv
        errconv = qmmm_options%errconv
        ndiis_matrices = qmmm_options%ndiis_matrices
        ndiis_attempts = qmmm_options%ndiis_attempts
        printcharges = qmmm_options%printcharges
        printbondorders = qmmm_options%printbondorders
        printdipole = qmmm_options%printdipole
        print_eigenvalues = qmmm_options%print_eigenvalues
        peptide_corr = qmmm_options%peptide_corr
        itrmax = qmmm_options%itrmax
        qmshake = qmmm_options%qmshake
        qmmask = qmmm_options%qmmask
        coremask = qmmm_options%coremask
        buffermask = qmmm_options%buffermask
        iqmatoms(1:MAX_QUANTUM_ATOMS) = qmmm_options%iqmatoms(1:MAX_QUANTUM_ATOMS)
        core_iqmatoms(1:MAX_QUANTUM_ATOMS) = qmmm_options%core_iqmatoms(1:MAX_QUANTUM_ATOMS)
        buffer_iqmatoms(1:MAX_QUANTUM_ATOMS) = qmmm_options%buffer_iqmatoms(1:MAX_QUANTUM_ATOMS)
        centermask = qmmm_options%centermask
        qmmmrij_incore = qmmm_options%qmmmrij_incore
        qmqm_erep_incore = qmmm_options%qmqm_erep_incore
        pseudo_diag = qmmm_options%pseudo_diag
        pseudo_diag_criteria = qmmm_options%pseudo_diag_criteria
        qm_ewald = qmmm_options%qm_ewald
        qm_pme = qmmm_options%qm_pme
        kmaxqx = qmmm_options%kmaxqx
        kmaxqy = qmmm_options%kmaxqy
        kmaxqz = qmmm_options%kmaxqz
        kappa = qmmm_options%kappa
        ksqmaxq = qmmm_options%ksqmaxq
        qmmm_int = qmmm_options%qmmm_int
        adjust_q = qmmm_options%adjust_q
        diag_routine = qmmm_options%diag_routine
        density_predict = qmmm_options%density_predict
        fock_predict = qmmm_options%fock_predict
        fockp_d1 = qmmm_options%fockp_d1
        fockp_d2 = qmmm_options%fockp_d2
        fockp_d3 = qmmm_options%fockp_d3
        fockp_d4 = qmmm_options%fockp_d4
        damp = qmmm_options%damp
        vshift = qmmm_options%vshift
        vsolv = qmmm_options%vsolv
        qmmm_switch = qmmm_options%qmmm_switch
        r_switch_hi = qmmm_options%r_switch_hi
        r_switch_lo = qmmm_options%r_switch_lo
        dftb_maxiter = qmmm_options%dftb_maxiter
        dftb_disper = qmmm_options%dftb_disper
        dftb_chg = qmmm_options%dftb_chg
        dftb_telec = qmmm_options%dftb_telec
        dftb_telec_step = qmmm_options%dftb_telec_step
        dftb_3rd_order = qmmm_options%dftb_3rd_order
        dftb_slko_path = qmmm_options%dftb_slko_path
        abfqmmm = qmmm_options%abfqmmm
        hot_spot = qmmm_options%hot_spot
        min_heavy_mass = qmmm_options%min_heavy_mass
        tight_p_conv = qmmm_options%tight_p_conv

        open(unit=10, file=filename, status="old")
        read(10, nml=cntrl)
        options%igb = igb
        options%alpb = alpb
        options%gbsa = gbsa
        options%lj1264 = lj1264
        options%ipb = ipb
        options%inp = inp
        options%vdwmeth = vdwmeth
        options%ew_type = ew_type
        options%extdiel = extdiel
        options%intdiel = intdiel
        options%rgbmax = rgbmax
        options%saltcon = saltcon
        options%cut = cut
        options%dielc = dielc
        options%ifqnt = ifqnt
        options%jfastw = jfastw
        options%ntf = ntf
        options%ntc = ntc
        options%fswitch = fswitch
        options%ntr = ntr
        options%ibelly = ibelly
        options%restraint_wt = restraint_wt
        options%restraintmask = restraintmask
        options%bellymask = bellymask
        options%refc = refc
        options%mask_from_ref = mask_from_ref

        qmcut = cut !in normal sander qmcut == cut by default, but not when API is used
        if (ifqnt > 0) then
            read(10, nml=qmmm)
            qmmm_options%qmcut = qmcut
            qmmm_options%lnk_dis = lnk_dis
            qmmm_options%lnk_atomic_no = lnk_atomic_no
            qmmm_options%lnk_method = lnk_method
            qmmm_options%qmgb = qmgb
            qmmm_options%qm_theory = qm_theory
            qmmm_options%qmcharge = qmcharge
            qmmm_options%corecharge = corecharge
            qmmm_options%buffercharge = buffercharge
            qmmm_options%spin = spin
            qmmm_options%qmqmdx = qmqmdx
            qmmm_options%verbosity = verbosity
            qmmm_options%scfconv = scfconv
            qmmm_options%errconv = errconv
            qmmm_options%ndiis_matrices = ndiis_matrices
            qmmm_options%ndiis_attempts = ndiis_attempts
            qmmm_options%printcharges = printcharges
            qmmm_options%printbondorders = printbondorders
            qmmm_options%printdipole = printdipole
            qmmm_options%print_eigenvalues = print_eigenvalues
            qmmm_options%peptide_corr = peptide_corr
            qmmm_options%itrmax = itrmax
            qmmm_options%qmshake = qmshake
            qmmm_options%qmmask = qmmask
            qmmm_options%coremask = coremask
            qmmm_options%buffermask = buffermask
            qmmm_options%iqmatoms(1:MAX_QUANTUM_ATOMS) = iqmatoms(1:MAX_QUANTUM_ATOMS)
            qmmm_options%core_iqmatoms(1:MAX_QUANTUM_ATOMS) = core_iqmatoms(1:MAX_QUANTUM_ATOMS)
            qmmm_options%buffer_iqmatoms(1:MAX_QUANTUM_ATOMS) = buffer_iqmatoms(1:MAX_QUANTUM_ATOMS)
            qmmm_options%centermask = centermask
            qmmm_options%qmmmrij_incore = qmmmrij_incore
            qmmm_options%qmqm_erep_incore = qmqm_erep_incore
            qmmm_options%pseudo_diag = pseudo_diag
            qmmm_options%pseudo_diag_criteria = pseudo_diag_criteria
            qmmm_options%qm_ewald = qm_ewald
            qmmm_options%qm_pme = qm_pme
            qmmm_options%kmaxqx = kmaxqx
            qmmm_options%kmaxqy = kmaxqy
            qmmm_options%kmaxqz = kmaxqz
            qmmm_options%kappa = kappa
            qmmm_options%ksqmaxq = ksqmaxq
            qmmm_options%qmmm_int = qmmm_int
            qmmm_options%adjust_q = adjust_q
            qmmm_options%diag_routine = diag_routine
            qmmm_options%density_predict = density_predict
            qmmm_options%fock_predict = fock_predict
            qmmm_options%fockp_d1 = fockp_d1
            qmmm_options%fockp_d2 = fockp_d2
            qmmm_options%fockp_d3 = fockp_d3
            qmmm_options%fockp_d4 = fockp_d4
            qmmm_options%damp = damp
            qmmm_options%vshift = vshift
            qmmm_options%vsolv = vsolv
            qmmm_options%qmmm_switch = qmmm_switch
            qmmm_options%r_switch_hi = r_switch_hi
            qmmm_options%r_switch_lo = r_switch_lo
            qmmm_options%dftb_maxiter = dftb_maxiter
            qmmm_options%dftb_disper = dftb_disper
            qmmm_options%dftb_chg = dftb_chg
            qmmm_options%dftb_telec = dftb_telec
            qmmm_options%dftb_telec_step = dftb_telec_step
            qmmm_options%dftb_3rd_order = dftb_3rd_order
            qmmm_options%dftb_slko_path = dftb_slko_path
            qmmm_options%abfqmmm = abfqmmm
            qmmm_options%hot_spot = hot_spot
            qmmm_options%min_heavy_mass = min_heavy_mass
            qmmm_options%tight_p_conv = tight_p_conv
        end if

        close(10)

    end subroutine read_sander_input

end program sander_driver