program circular_inductors
    
    use hdf5
    
    use utils_library

    use materials_library
    use sources_library
    use components_library
    use circular_inductors_library
    use inductor_solvers 

    implicit none

    type(Cu_Material), target :: Cu

    type(port_t) :: port1, port2
    type(solver) :: solver1
    type(sinusoidal_source) :: source1
    type(sresistor) :: resistive_load
    type(circularcs_rldowell_inductor) :: inductor1, inductor2

    integer, parameter :: npoints = 200
    integer :: unit, i, j
    real(dp) :: inductances(2,2), xyz(3,npoints**2) 
    complex(dp) :: currents(2,1), bfield(3,npoints**2,1)

    ! source configuration

    call source1%set_resistance(50.0_dp) ! 50 ohm source Rs
    call source1%configure_source(cmplx(1.0_dp,0.0_dp,kind=dp), 200E3_dp) ! 1V, 200kHz sinusoidal source
    
    ! load configuration

    call resistive_load%set_resistance(50.0_dp) ! 50 ohm load impedance

    ! TX Inductor

    inductor1%material => Cu
    inductor1%center = (/0.0_dp, 0.0_dp, 0.0_dp/)
    inductor1%normal = (/0.0_dp, 0.0_dp, 1.0_dp/)
    inductor1%radius = 5E-2_dp
    inductor1%radius_section = 5E-4_dp

    ! RX Inductor

    inductor2%material => Cu
    inductor2%center = (/0.0_dp, 0.0_dp, 5E-2_dp/)
    inductor2%normal = (/0.0_dp, 0.0_dp, 1.0_dp/)
    inductor2%radius = 5E-2_dp
    inductor2%radius_section = 5E-4_dp

    ! mesh setup

    call port1%add_inductor(inductor1)
    call port1%add_load(source1)
    
    call port2%add_inductor(inductor2)
    call port2%add_load(resistive_load)

    ! solver setup

    call solver1%add_port(port1)
    call solver1%add_port(port2) 
    
    ! simulation

    call solver1%set_samples(10000)
    call solver1%solve_inductances()

    call solver1%get_inductances(inductances)

    open(newunit=unit, file="sim.csv", status='replace', action='write', iostat=i)
    if (i /= 0) then
        print *, "Error abriendo archivo"
        stop
    end if

    do i = 1, size(inductances,1)
        do j = 1, size(inductances,2)
            write(unit, '(I3,",",I3,",",E12.5)') i,j,inductances(i,j)
        end do
    end do

    close(unit)

    call solver1%run()
    call solver1%get_currents(currents)

    open(newunit=unit, file="current.csv", status='replace', action='write', iostat=i)
    if (i /= 0) then
        print *, "Error abriendo archivo"
        stop
    end if

    do i=1,2
        write(unit,'(I3,",",ES12.5,",",ES12.5)') i, real(currents(i,1), kind=dp), aimag(currents(i,1))
    end do

    close(unit)
    
    call generate_xzmesh(npoints, -20E-2_dp, 20E-2_dp, -10E-2_dp, 10E-2_dp, 0.0_dp, xyz)
    call solver1%get_bfield(npoints*npoints, xyz, bfield)

    open(newunit=unit, file="bfield.csv", status='replace', action='write', iostat=i)
    if (i /= 0) then
        print *, "Error abriendo archivo"
        stop
    end if

    do i=1,npoints*npoints
        write(unit,'(ES12.5,",",ES12.5,",",ES12.5,",",ES12.5,",",ES12.5,",",ES12.5)') xyz(1,i), xyz(2,i), &
            & xyz(3,i), abs(bfield(1,i,1)), abs(bfield(2,i,1)), abs(bfield(3,i,1))
    end do

    close(unit)
    
    ! data saving & memory deallocation
    
    !call solver1%kill()
end program
