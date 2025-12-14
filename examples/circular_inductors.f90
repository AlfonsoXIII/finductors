program circular_inductors
    
    use hdf5

    use materials_library
    use sources_library
    use components_library
    use circular_inductors_library
    use inductor_solvers 
    implicit none
    
    type(port_t) :: port1, port2
    type(solver_t) :: solver
    type(sinusoidal_source) :: source1
    type(presistor) :: resistive_load
    type(circularcs_rldowell_inductor) :: inductor1, inductor2

    real(dp) :: inductances(2,2)

    integer(HID_T) :: file_id, dataspace_id, dataset_id
    integer(HSIZE_T), dimensions(2) :: dims
    integer :: error

    ! source configuration

    call source1%set_resistance(50) ! 50 ohm source Rs
    call source1%configure_source(1, 200E3_dp) ! 1V, 200kHz sinusoidal source
    
    ! load configuration

    call resistive_load%set_resistance(50) ! 50 ohm load impedance

    ! TX Inductor

    inductor1%material = Cu 
    inductor1%center = (/0.0_dp, 0.0_dp, 0.0_dp/)
    inductor1%normal = (/0.0_dp, 0.0_dp, 1.0_dp/)
    inductor1%radius = 5E-2_dp
    inductor1%radius_section = 5E-4_dp

    ! RX Inductor

    inductor2%material = Cu
    inductor2%center = (/0.0_dp, 0.0_dp, 5E-2_dp/)
    inductor2%normal = (/0.0_dp, 0.0_dp, 1.0_dp/)
    inductor2%radius = 5E-2_dp
    inductor2%radius_section = 5E-4_dp

    ! mesh setup

    call port1%init()
    call port1%add_inductor(inductor1)
    call port1%add_load(source1)

    call port2%init()
    call port2%add_inductor(inductor2)
    call port2%add_load(resistive_load)

    ! solver setup

    call solver%init()

    call solver%add_port(port1)
    call solver%add_port(port2) 
    
    ! simulation

    call solver%set_samples(10E3)
    call solver%solve_inductances()

    call solver%get_inductances(inductances)
    
    ! hdf5 data dumping

    dims = (/2,2/)

    file_id = H5Fcreate("simulation.h5", H5F_ACC_TRUNC_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
    dataspace_id = H5Screate_simple(2, dims, dims)
    dataset_id = H5Dcreate(file_id, "inductances", H5T_NATIVE_DOUBLE, dataspace_id, &
        H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)

    call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, inductances, dims, error)
    if (error /= 0) print *, "Error opening the dataset!"

    call H5Dclose_f(dataset_id, error)
    call H5Dclose_f(dataspace_id, error)
    call H5Dclose_f(file_id, error)

    ! data saving & memory deallocation
    
    call solver%kill()
end program
