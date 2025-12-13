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
        class(inductor_t), allocatable :: inductor
        class(component_t), allocatable :: load
        class(component_t), allocatable :: components(:)
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
        procedure :: compute_inductances => solver_compute_inductances
        procedure :: get_inductances => solver_get_inductances
        !procedure :: get_currents => solver_get_currents
        !procedure(solver_get_efield) :: get_efield
        !procedure(solver_get_bfield) :: get_bfield

        procedure, private :: mutual_inductance => solve_mutual_inductance
        procedure, private :: inductance => solve_inductance
 
    end type solver

contains

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
    
    deallocate(self%ports, self%v_vector, self%z_matrix)
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
    
    real(dp) :: harmonic_buff
    real(dp), allocatable :: harmonics_buff(:)
    complex(dp), allocatable :: complex_buff(:)

    integer :: idx, counter
    integer, allocatable :: sources_buff(:)

    ! harmonics fetching
    
    counter = 0
    do idx=1,self%n_ports
        select type(src => self%ports(idx)%load)
            class is (source_t)
                counter = counter + src%n_harmonics
        end select
    end do

    if (counter == 0) then
        print *, "No sources configured!"
        return
    end if
    
    allocate(sources_buff(counter), harmonics_buff(counter), complex_buff(counter))

    do idx=1,self%n_ports
        select type(src => self%ports(idx)%load)
            class is (source_t)

        end select
    end do

    call hsa_quicksort(counter, harmonics_buff, sources_buff, complex_buff, 1, counter)

    ! grouping
    
    self%n_harmonics = 1
    harmonic_buff = harmonics_buff(1)

    do idx=1,counter
        if (harmonics_buff(idx) /= harmonic_buff)
            self%n_harmonics = counter_group + 1
            harmonic_buff = harmonics_buff(idx)
        end if
    end do

    allocate(self%harmonics(self%n_harmonics))

    harmonic_buff = harmonics_buff(1)
    do idx=1,counter
        if (harmonics_buff(idx) /= harmonic_buff)
            harmonic_buff = harmonics_buff(idx)
        end if

    end do

    ! simulation sweep

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
