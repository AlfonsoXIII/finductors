module polygonal_inductors_library
    use types
    use base_inductors_library
    implicit none

    ! polygonal inductor with circular section

    type, extends(inductor_t) :: polygon_inductor
        integer :: N
        real(dp) :: side
    contains
        procedure :: get_r => polygon_get_r
        procedure :: get_dl => polygon_get_dl
        procedure :: get_limits => polygon_get_limits
    end type polygon_inductor

contains

    ! POLYGONAL INDUCTOR - CIRCULAR SECTION

    subroutine polygon_get_r(self, N, theta, res)
        class(polygon_inductor), intent(in) :: self
        integer, intent(in) :: N
        real(dp), intent(in) :: theta(N)
        real(dp), intent(out) :: res(3,N)
        
        real(dp) :: theta_center, alpha
        real(dp), allocatable :: r(:), phi(:)
        
        allocate(phi(N), r(N))
        if (mod(self%N,2) == 0) then
            theta_center = PI / real(self%N, kind=dp)
        else
            theta_center = 0.0_dp
        end if
    
        alpha = 2.0_dp*PI/real(self%N,dp)
        phi(:) = theta(:) - theta_center
        phi(:) = phi(:) - alpha * floor(phi(:)/alpha + 0.5_dp)
        r(:) = self%side * cos(alpha/2) / (2*sin(alpha/2)*cos(phi(:)))
        
        res(1,:) = r * cos(theta(:)) + self%center(1)
        res(2,:) = r * sin(theta(:)) + self%center(2)
        res(3,:) = self%center(3)
        call rotate_vector(N, self%normal, res)
        deallocate(phi, r)
    end subroutine polygon_get_r

    subroutine polygon_get_dl(self, N, theta, res)
        class(polygon_inductor), intent(in) :: self
        integer, intent(in) :: N
        real(dp), intent(in) :: theta(N)
        real(dp), intent(out) :: res(3,N)

        real(dp) :: theta_center, alpha
        real(dp), allocatable :: phi(:), r(:), rp(:)
        
        allocate(phi(N), r(N), rp(N))
        if (mod(self%N,2) == 0) then
            theta_center = PI / real(self%N, kind=dp)
        else
            theta_center = 0.0_dp
        end if
        
        alpha = 2.0_dp * PI / real(self%N, dp)
        phi(:) = theta(:) - theta_center
        phi(:) = phi(:) - alpha * floor(phi(:)/alpha + 0.5_dp)
        rp(:) = self%side * sin(phi(:)) * cos(alpha/2) / (2*sin(alpha/2)*cos(phi(:))**2)
        r(:) = self%side * cos(alpha/2) / (2*sin(alpha/2)*cos(phi(:)))
        
        res(1,:) = rp(:) * cos(theta(:)) - r(:) * sin(theta(:))
        res(2,:) = rp(:) * sin(theta(:)) + r(:) * cos(theta(:))
        res(3,:) = 0.0_dp
        call rotate_vector(N, self%normal, res)
        deallocate(phi, r, rp)
    end subroutine polygon_get_dl

    subroutine polygon_get_limits(self, res)
        class(polygon_inductor), intent(in) :: self
        real(dp), intent(out) :: res(2)
        res(:) = (/0.0_dp, 2.0_dp*PI/)
    end subroutine polygon_get_limits
end module 
