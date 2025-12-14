module inductor_solvers
    use types
    use math_library
    use components_library
    use sources_library
    use base_inductors_library
    implicit none

    private
    public :: port_t, solver

    type :: inductor_ptr
        class(inductor_t), pointer :: ptr => null()
    end type inductor_ptr

    type :: component_ptr
        class(component_t), pointer :: ptr => null()
    end type component_ptr

    type :: port_t
        integer :: n_components = 0
        type(inductor_ptr) :: inductor
        type(component_ptr) :: load
        type(component_ptr), allocatable :: components(:)
    contains
        procedure :: kill => port_kill
        procedure :: add_component => port_add_component
        procedure :: add_load => port_add_load
        procedure :: add_inductor => port_add_inductor
    end type port_t

    type :: port_ptr
        type(port_t), pointer :: ptr => null()
    end type port_ptr

    type :: solver
    private
        real(dp) :: temperature = 20
        integer :: n_ports = 0, n_samples = 0, n_harmonics = 0
        type(port_ptr), allocatable :: ports(:)
        complex(dp), allocatable :: v_vector(:,:)
        real(dp), allocatable :: harmonics(:), samples(:,:), inductances(:,:)
    contains
        procedure :: kill => solver_kill
        procedure :: set_temperature => solver_set_temperature
        procedure :: set_samples => solver_set_samples
        procedure :: add_port => solver_add_port
        procedure :: run => solver_run
        procedure :: solve_inductances => solver_compute_inductances
        procedure :: get_inductances => solver_get_inductances
        procedure :: get_currents => solver_get_currents
        !procedure :: get_harmonics => solver_get_harmonics
        !procedure :: get_efield => solver_get_efield
        !procedure :: get_bfield => solver_get_bfield
        !procedure :: get_power_distribution => solver_get_power_distribution

        procedure, private :: mutual_inductance => solve_mutual_inductance
        procedure, private :: inductance => solve_inductance
 
    end type solver

contains

subroutine port_kill(self)
    class(port_t), intent(inout) :: self
    deallocate(self%components)
end subroutine

subroutine port_add_component(self, new_component)
    class(port_t), intent(inout) :: self
    class(component_t), target, intent(inout) :: new_component
    
    integer :: idx  
    type(component_ptr), allocatable :: tmp(:)
    
    if (allocated(self%components)) then
        allocate(tmp(self%n_components))
        do idx=1,self%n_components
            tmp(idx)%ptr => self%components(idx)%ptr
        end do
        deallocate(self%components)
    end if

    allocate(self%components(self%n_components+1))

    do idx=1,self%n_components
        self%components(idx)%ptr => tmp(idx)%ptr
    end do

    if (allocated(tmp)) deallocate(tmp)

    self%n_components = self%n_components + 1
    self%components(self%n_components)%ptr => new_component
end subroutine

subroutine port_add_load(self, new_load)
    class(port_t), intent(inout) :: self
    class(component_t), intent(in), target :: new_load
    self%load%ptr => new_load
end subroutine

subroutine port_add_inductor(self, new_inductor)
    class(port_t), intent(inout) :: self
    class(inductor_t), intent(in), target :: new_inductor
    self%inductor%ptr => new_inductor
end subroutine

subroutine solver_kill(self)
    class(solver), intent(inout) :: self
    integer :: idx
    !do idx=1,self%n_ports
    !    call self%ports(idx)%ptr%kill()
    !end do
    deallocate(self%ports, self%v_vector)
    deallocate(self%harmonics, self%samples, self%inductances)
end subroutine

subroutine solver_set_temperature(self, new_temperature)
    class(solver), intent(inout) :: self
    real(dp), intent(in) :: new_temperature
    self%temperature = new_temperature
end subroutine

subroutine solver_set_samples(self, new_samples)
    class(solver), intent(inout) :: self
    integer, intent(in) :: new_samples
    
    self%n_samples = new_samples
    
    if (allocated(self%samples)) deallocate(self%samples)
    allocate(self%samples(new_samples,2))

    call hammersley2d_sequence(new_samples, self%samples)
end subroutine

subroutine solver_add_port(self, new_port)
    class(solver), intent(inout) :: self
    type(port_t), intent(inout), target :: new_port
    
    integer :: idx
    type(port_ptr), allocatable :: tmp(:)

    if (allocated(self%ports)) then
        allocate(tmp(self%n_ports))
        do idx=1,self%n_ports
            tmp(idx)%ptr => self%ports(idx)%ptr
        end do
        deallocate(self%ports)
    end if

    allocate(self%ports(self%n_ports+1))

    do idx=1,self%n_ports
        self%ports(idx)%ptr => tmp(idx)%ptr
    end do

    if (allocated(tmp)) deallocate(tmp)

    self%n_ports = self%n_ports + 1
    self%ports(self%n_ports)%ptr => new_port
end subroutine

subroutine solver_run(self)
    class(solver), intent(inout) :: self
    
    integer :: idx, jdx, ref, counter
    real(dp) :: harmonic_buff
    integer, allocatable :: sources_buff(:), ipiv(:)
    real(dp), allocatable :: harmonics_buff(:)
    complex(dp) :: z_buff(4), z_buff_aux(4)
    complex(dp), allocatable :: amplitudes_buff(:), z_matrix(:,:)

    ! INDUCTANCES INITIAL CHECK

    if (size(self%inductances) /= self%n_ports*self%n_ports) then
        print *, "Recalculating Inductances..."
        call self%solve_inductances()
    end if

    ! HARMONICS FETCHING
    
    counter = 0
    do idx=1,self%n_ports
        select type(src => self%ports(idx)%ptr%load%ptr)
            class is (source_t)
                counter = counter + src%n_harmonics
        end select
    end do

    if (counter == 0) then
        print *, "No Sources Configured! Aborting..."
        return
    end if
    
    allocate(sources_buff(counter), harmonics_buff(counter), amplitudes_buff(counter))

    ref = 0
    do idx=1,self%n_ports
        select type(src => self%ports(idx)%ptr%load%ptr)
            class is (source_t)
                !$omp simd aligned(sources_buff, harmonics_buff, amplitudes_buff:64)
                do jdx=1,src%n_harmonics
                    sources_buff(jdx+ref) = idx
                    harmonics_buff(jdx+ref) = src%harmonics(jdx)%frequency
                    amplitudes_buff(jdx+ref) = src%harmonics(jdx)%amplitude
                end do
                ref = ref + src%n_harmonics
        end select
    end do

    call hsa_quicksort(counter, harmonics_buff, sources_buff, amplitudes_buff, 1, counter)

    ! HARMONICS GROUPING & DENORMALIZED EXCITATION VECTOR ASSEMBLY
    
    self%n_harmonics = 1
    harmonic_buff = harmonics_buff(1)

    do idx=1,counter
        if (harmonics_buff(idx) /= harmonic_buff) then
            self%n_harmonics = self%n_harmonics + 1
            harmonic_buff = harmonics_buff(idx)
        end if
    end do

    if (allocated(self%v_vector)) deallocate(self%v_vector)
    allocate(self%v_vector(self%n_ports,self%n_harmonics), self%harmonics(self%n_harmonics))
    
    self%v_vector = CZERO

    ref = 0
    harmonic_buff = -1
    do idx=1,counter
        
        if (harmonics_buff(idx) /= harmonic_buff) then
            ref = ref + 1
            harmonic_buff = harmonics_buff(idx)
            self%harmonics(ref) = harmonic_buff
        end if
        
        self%v_vector(sources_buff(idx),ref) = amplitudes_buff(idx)

    end do

    ! SIMULATION SWEEP
    
    !$omp parallel do schedule(static) default(shared) private(idx, jdx, z_matrix, z_buff, z_buff_aux, ipiv)
    do idx=1,self%n_harmonics
        
        ! impedance matrix assembly
        
        allocate(ipiv(self%n_ports))
        allocate(z_matrix(self%n_ports,self%n_ports)) 

        z_matrix = self%inductances * CJ * 2.0_dp * PI * self%harmonics(idx)

        do jdx=1,self%n_ports
            
            z_buff = (/CONE, CZERO, CONE, CZERO/)
            
            ! source / load

            call self%ports(jdx)%ptr%load%ptr%get_abcd(self%temperature, self%harmonics(idx), z_buff_aux)
            call cmmp(z_buff, z_buff_aux)

            ! components

            do ref=1,self%ports(jdx)%ptr%n_components
                call self%ports(jdx)%ptr%components(ref)%ptr%get_abcd(self%temperature, self%harmonics(idx), z_buff_aux)
                call cmmp(z_buff, z_buff_aux)
            end do

            ! inductor parasitic model

            call self%ports(jdx)%ptr%inductor%ptr%get_abcd(self%temperature, self%harmonics(idx), z_buff_aux)
            call cmmp(z_buff, z_buff_aux)

            ! voltage normalization

            self%v_vector(jdx,idx) = self%v_vector(jdx,idx) / z_buff(1)

            ! series impedance

            z_matrix(jdx,jdx) = z_matrix(jdx,jdx) + z_buff(2) / z_buff(1)

        end do

        ! linear system solving ( lapack call )

        call zgesv(self%n_ports, 1, z_matrix, self%n_ports, ipiv, self%v_vector(:,idx), self%n_ports, ref)
        
        if (ref /= 0) then
            print *, "Matrix Inversion Failed: Singular Matrix!"
        end if 

        deallocate(ipiv, z_matrix)

    end do

    deallocate(harmonics_buff, sources_buff, amplitudes_buff)
end subroutine

subroutine solver_get_currents(self, currents)
    class(solver), intent(inout) :: self
    complex(dp), intent(out) :: currents(self%n_ports,self%n_harmonics)
    
    integer :: idx, jdx
    
    do jdx=1,self%n_harmonics
        !$omp simd
        do idx=1,self%n_ports
            currents(idx,jdx) = self%v_vector(idx,jdx)
        end do
    end do
end subroutine

subroutine solver_compute_inductances(self)
    class(solver), intent(inout) :: self

    integer :: idx, jdx

    if (.not. allocated(self%inductances)) then
        allocate(self%inductances(self%n_ports,self%n_ports))
    else if (size(self%inductances) /= self%n_ports*self%n_ports) then
        deallocate(self%inductances)
        allocate(self%inductances(self%n_ports,self%n_ports))
    end if

    do idx=1,self%n_ports
        do jdx=1,idx
            if (jdx==idx) then
                call self%inductance(idx)
            else
                call self%mutual_inductance(idx,jdx)
                self%inductances(jdx,idx) = self%inductances(idx,jdx)
            end if
        end do
    end do
end subroutine

subroutine solver_get_inductances(self, inductances)
    class(solver), intent(inout) :: self
    real(dp), intent(out) :: inductances(self%n_ports,self%n_ports)
    
    integer :: idx, jdx

    do idx=1,self%n_ports
        !$omp simd
        do jdx=1,idx
            inductances(jdx,idx) = self%inductances(jdx,idx)
            inductances(idx,jdx) = self%inductances(jdx,idx)
        end do
    end do
end subroutine

subroutine solve_mutual_inductance(self, ind1, ind2)
    class(solver), intent(inout) :: self
    integer, intent(in) :: ind1, ind2
    
    integer :: i,j
    real(dp) :: limits1(2), limits2(2), dist, buffer
    real(dp), allocatable :: theta1(:), theta2(:)
    real(dp), allocatable :: r1(:,:), r2(:,:)
    real(dp), allocatable :: dl1(:,:), dl2(:,:)
    
    allocate(theta1(self%n_samples), theta2(self%n_samples))
    allocate(r1(3,self%n_samples), r2(3,self%n_samples))
    allocate(dl1(3,self%n_samples), dl2(3,self%n_samples))

    call self%ports(ind1)%ptr%inductor%ptr%get_limits(limits1)
    call self%ports(ind2)%ptr%inductor%ptr%get_limits(limits2)
    
    !$omp simd
    do i=1,self%n_samples
        theta1(i) = limits1(1) + self%samples(i,1)*(limits1(2)-limits1(1))
        theta2(i) = limits2(1) + self%samples(i,2)*(limits2(2)-limits2(1))
    end do

    call self%ports(ind1)%ptr%inductor%ptr%get_r(self%n_samples,theta1,r1)
    call self%ports(ind2)%ptr%inductor%ptr%get_r(self%n_samples,theta2,r2)
    call self%ports(ind1)%ptr%inductor%ptr%get_dl(self%n_samples,theta1,dl1)
    call self%ports(ind2)%ptr%inductor%ptr%get_dl(self%n_samples,theta2,dl2)

    buffer = 0.0_dp
    !$omp parallel do private(i,dist) reduction(+:buffer) schedule(static)
    do i=1,self%n_samples
        dist = max(sqrt(sum((r1(:,i)-r2(:,i))**2)), 1.0E-12_dp)
        buffer = buffer + dot_product(dl1(:,i),dl2(:,i)) / dist
    end do
    !$omp end parallel do
   
    buffer = ( buffer * 1.0E-7_dp * limits1(2) * limits2(2) ) / real(self%n_samples, kind=dp)
    self%inductances(ind1,ind2) = buffer

    deallocate(theta1,theta2,r1,r2,dl1,dl2)
end subroutine solve_mutual_inductance

subroutine solve_inductance(self, ind)
    class(solver), intent(inout) :: self
    integer, intent(in) :: ind
    
    integer :: i
    real(dp) :: limits(2), dist, buffer
    real(dp), allocatable :: r1(:,:), r2(:,:)
    real(dp), allocatable :: theta1(:), theta2(:)
    real(dp), allocatable :: dl1(:,:), dl2(:,:)
    
    allocate(theta1(self%n_samples), theta2(self%n_samples))
    allocate(dl1(3,self%n_samples), dl2(3,self%n_samples))
    allocate(r1(3,self%n_samples), r2(3,self%n_samples))

    call self%ports(ind)%ptr%inductor%ptr%get_limits(limits)

    !$omp simd
    do i=1,self%n_samples
        theta1(i) = limits(1) + self%samples(i,1) * (limits(2)-limits(1))
        theta2(i) = limits(1) + self%samples(i,2) * (limits(2)-limits(1))
    end do

    call self%ports(ind)%ptr%inductor%ptr%get_r(self%n_samples,theta1,r1) 
    call self%ports(ind)%ptr%inductor%ptr%get_r(self%n_samples,theta2,r2)
    call self%ports(ind)%ptr%inductor%ptr%get_dl(self%n_samples,theta1,dl1)
    call self%ports(ind)%ptr%inductor%ptr%get_dl(self%n_samples,theta2,dl2)
    call self%ports(ind)%ptr%inductor%ptr%get_gmd(dist)
    
    !$omp simd
    do i=1,self%n_samples
        r2(:,i) = r2(:,i) + dist * self%ports(ind)%ptr%inductor%ptr%normal
    end do

    buffer = 0.0_dp
    !$omp parallel do private(i) reduction(+:buffer) schedule(static)
    do i=1,self%n_samples
        dist = max(sqrt(sum((r1(:,i)-r2(:,i))**2)),1.0E-12_dp)
        buffer = buffer + dot_product(dl1(:,i),dl2(:,i)) / dist
    end do
    !$omp end parallel do

    buffer = ( buffer * 1.0E-7_dp * limits(2)**2 ) / real(self%n_samples, kind=dp)
    self%inductances(ind,ind) = buffer

    deallocate(r1,r2,theta1,theta2,dl1,dl2)
end subroutine solve_inductance

end module inductor_solvers
