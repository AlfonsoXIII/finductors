module inductor_solvers
    use types
    use math_library
    use components_library
    use sources_library
    use base_inductors_library
    implicit none

    private
    public :: port_t, solver

    type :: port_t
        integer :: n_components
        class(inductor_t), allocatable :: inductor
        class(component_t), allocatable :: load
        class(component_t), allocatable :: components(:)
    contains
        procedure :: add_component => port_add_component
        procedure :: add_load => port_add_load
        procedure :: add_inductor => port_add_inductor
    end type port_t

    type :: solver
    private
        real(dp) :: temperature
        integer :: n_ports, n_samples, n_harmonics
        type(port_t), allocatable :: ports(:)
        complex(dp), allocatable :: v_vector(:,:)
        real(dp), allocatable :: harmonics(:), samples(:,:), inductances(:,:)

    contains
        procedure :: init => solver_init
        procedure :: kill => solver_kill
        procedure :: set_temperature => solver_set_temperature
        procedure :: set_samples => solver_set_samples
        procedure :: add_port => solver_add_port
        procedure :: run => solver_run
        procedure :: solve_inductances => solver_compute_inductances
        procedure :: get_inductances => solver_get_inductances
        !procedure :: get_harmonics => solver_get_harmonics
        !procedure :: get_currents => solver_get_currents
        !procedure :: get_efield => solver_get_efield
        !procedure :: get_bfield => solver_get_bfield
        !procedure :: get_power_distribution => solver_get_power_distribution

        procedure, private :: mutual_inductance => solve_mutual_inductance
        procedure, private :: inductance => solve_inductance
 
    end type solver

contains

subroutine port_add_component(self, new_component)
    class(port_t), intent(inout) :: self
    class(component_t), intent(in) :: new_component


end subroutine

subroutine port_add_load(self, new_load)
    class(port_t), intent(inout) :: self
    class(component_t), intent(in) :: new_load


end subroutine

subroutine port_add_inductor(self, new_inductor)
    class(port_t), intent(inout) :: self
    class(inductor_t), intent(in) :: new_inductor


end subroutine

subroutine solver_init(self)
    class(solver), intent(inout) :: self
    
    self%temperature = 20
    
    self%n_harmonics = 0
    allocate(self%harmonics(0))

    self%n_ports = 0
    allocate(self%ports(0))

    allocate(self%v_vector(0,0))

    self%n_samples = 0
    allocate(self%samples(0,0))

    allocate(self%inductances(0,0))
end subroutine

subroutine solver_kill(self)
    class(solver), intent(inout) :: self
    
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
    integer :: new_samples

    self%n_samples = new_samples

    allocate(self%samples(new_samples,2))
    call hammersley2d_sequence(new_samples, self%samples)
end subroutine

subroutine solver_add_port(self, new_port)
    class(solver), intent(inout) :: self
    type(port_t), intent(inout) :: new_port
    
    self%n_ports = self%n_ports + 1
    self%ports = reshape(self%ports, [self%n_ports])

    self%ports(self%n_ports) = new_port
end subroutine

subroutine solver_run(self)
    class(solver), intent(inout) :: self
    
    integer :: idx, jdx, ref, counter
    integer, allocatable :: sources_buff(:), ipiv(:)
    real(dp) :: harmonic_buff
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
        select type(src => self%ports(idx)%load)
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
        select type(src => self%ports(idx)%load)
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

    deallocate(self%v_vector)
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
    
    !$omp parallel do schedule(static) default(shared) private(idx, jdx, z_matrix, z_buff, z_buff_aux)
    do idx=1,self%n_harmonics
        
        ! impedance matrix assembly
        
        allocate(ipiv(self%n_ports))
        allocate(z_matrix(self%n_ports,self%n_ports)) 

        z_matrix = self%inductances * CJ * 2.0_dp * PI * self%harmonics(idx)
        
        do jdx=1,self%n_ports
            
            z_buff = (/CONE, CZERO, CONE, CZERO/)
            
            ! source / load

            call self%ports(jdx)%load%get_abcd(self%temperature, self%harmonics(idx), z_buff_aux)
            call cmmp(z_buff, z_buff_aux)
            
            ! components

            do ref=1,self%ports(jdx)%n_components
                call self%ports(jdx)%components(ref)%get_abcd(self%temperature, self%harmonics(idx), z_buff_aux)
                call cmmp(z_buff, z_buff_aux)
            end do

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

!subroutine solver_get_currents(self, i_vector)
!    class(solver), intent(inout) :: self
!    complex(dp), intent(out) :: i_vector(self%n_ports)
!    
!    integer :: idx
!    
!    !$omp simd
!    do idx=1,self%n_ports
!        i_vector(idx) = self%v_vector(idx)
!    end do
!end subroutine

subroutine solver_compute_inductances(self)
    class(solver), intent(inout) :: self

    integer :: i, j

    if (size(self%inductances) /= self%n_ports**2) then
       self%inductances = reshape(self%inductances,[self%n_ports,self%n_ports]) 
    end if

    do i=1,self%n_ports
        do j=1,i
            if (j==i) then
                call self%inductance(i)
            else
                call self%mutual_inductance(i,j)
            end if
        end do
    end do
end subroutine

subroutine solver_get_inductances(self, inductances)
    class(solver), intent(inout) :: self
    real(dp), intent(out) :: inductances(self%n_ports,self%n_ports)
    
    integer :: i, j

    do i=1,self%n_ports
        !$omp simd
        do j=1,i
            inductances(j,i) = self%inductances(j,i)
            inductances(i,j) = self%inductances(j,i)
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

    call self%ports(ind1)%inductor%get_limits(limits1)
    call self%ports(ind2)%inductor%get_limits(limits2)
    
    !$omp simd
    do i=1,self%n_samples
        theta1(i) = limits1(1) + self%samples(i,1)*(limits1(2)-limits1(1))
        theta2(i) = limits2(1) + self%samples(i,2)*(limits2(2)-limits2(1))
    end do

    call self%ports(ind1)%inductor%get_r(self%n_samples,theta1,r1)
    call self%ports(ind2)%inductor%get_r(self%n_samples,theta2,r2)
    call self%ports(ind1)%inductor%get_dl(self%n_samples,theta1,dl1)
    call self%ports(ind2)%inductor%get_dl(self%n_samples,theta2,dl2)

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

    call self%ports(ind)%inductor%get_limits(limits)

    !$omp simd
    do i=1,self%n_samples
        theta1(i) = limits(1) + self%samples(i,1) * (limits(2)-limits(1))
        theta2(i) = limits(1) + self%samples(i,2) * (limits(2)-limits(1))
    end do

    call self%ports(ind)%inductor%get_r(self%n_samples,theta1,r1) 
    call self%ports(ind)%inductor%get_r(self%n_samples,theta2,r2)
    call self%ports(ind)%inductor%get_dl(self%n_samples,theta1,dl1)
    call self%ports(ind)%inductor%get_dl(self%n_samples,theta2,dl2)
    call self%ports(ind)%inductor%get_gmd(dist)
    
    !$omp simd
    do i=1,self%n_samples
        r2(:,i) = r2(:,i) + dist * self%ports(ind)%inductor%normal
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
