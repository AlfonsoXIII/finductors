module circular_inductors_library
    use types
    use base_inductors_library
    implicit none
   
    private
    public :: circularcs_rldrude_inductor, circularfs_rldrude_inductor, &
        & circularcs_rldowell_inductor, circularfs_rldowell_inductor, &
        & lcspiralcs_rldrude_inductor, lcspiralfs_rldrude_inductor, &
        & lcspiralcs_rldowell_inductor, lcspiralfs_rldowell_inductor
    
    ! SIMPLE CIRCULAR INDUCTOR

    type, extends(inductor_t), abstract :: circular_inductor
        real(dp) :: radius
    contains
        procedure :: get_r => circular_get_r
        procedure :: get_dl => circular_get_dl
        procedure :: get_limits => circular_get_limits
    end type circular_inductor

    type, extends(circular_inductor), abstract :: circularcs_inductor
        real(dp) :: radius_section
    contains
        procedure :: get_gmd => circularcs_get_gmd
    end type circularcs_inductor

    type, extends(circular_inductor), abstract :: circularfs_inductor
        real(dp) :: width_section, thickness_section
    contains
        procedure :: get_gmd => circularfs_get_gmd
    end type circularfs_inductor

    ! Circular Section

    type, extends(circularcs_inductor) :: circularcs_rldrude_inductor
    contains
        procedure :: get_abcd => circularcs_rldrude_get_abcd
    end type circularcs_rldrude_inductor

    type, extends(circularcs_inductor) :: circularcs_rldowell_inductor
    contains
        procedure :: get_abcd => circularcs_rldowell_get_abcd
    end type circularcs_rldowell_inductor

    ! Flat Section

    type, extends(circularfs_inductor) :: circularfs_rldrude_inductor
    contains
        procedure :: get_abcd => circularfs_rldrude_get_abcd
    end type circularfs_rldrude_inductor

    type, extends(circularfs_inductor) :: circularfs_rldowell_inductor
    contains
        procedure :: get_abcd => circularfs_rldowell_get_abcd
    end type circularfs_rldowell_inductor

    ! LINEAR CIRCULAR SPIRAL

    type, extends(inductor_t), abstract :: lcspiral_inductor
        real(dp) :: radius, pitch, turns
    contains
        procedure :: get_r => lcspiral_get_r
        procedure :: get_dl => lcspiral_get_dl
        procedure :: get_limits => lcspiral_get_limits
    end type lcspiral_inductor

    type, extends(lcspiral_inductor), abstract :: lcspiralcs_inductor
        real(dp) :: radius_section
    contains
        procedure :: get_gmd => lcspiralcs_get_gmd
    end type lcspiralcs_inductor

    type, extends(lcspiral_inductor), abstract :: lcspiralfs_inductor
        real(dp) :: width_section, thickness_section
    contains
        procedure :: get_gmd => lcspiralfs_get_gmd
    end type lcspiralfs_inductor
    
    ! Circular Section

    type, extends(lcspiralcs_inductor) :: lcspiralcs_rldrude_inductor
    contains
        procedure :: get_abcd => lcspiralcs_rldrude_get_abcd
    end type lcspiralcs_rldrude_inductor

    type, extends(lcspiralcs_inductor) :: lcspiralcs_rldowell_inductor
    contains
        procedure :: get_abcd => lcspiralcs_rldowell_get_abcd
    end type lcspiralcs_rldowell_inductor

    ! Flat Section

    type, extends(lcspiralfs_inductor) :: lcspiralfs_rldrude_inductor
    contains
        procedure :: get_abcd => lcspiralfs_rldrude_get_abcd
    end type lcspiralfs_rldrude_inductor

    type, extends(lcspiralfs_inductor) :: lcspiralfs_rldowell_inductor
    contains
        procedure :: get_abcd => lcspiralfs_rldowell_get_abcd
    end type lcspiralfs_rldowell_inductor

contains
    
    ! CIRCULAR INDUCTOR

    subroutine circular_get_r(self, N, theta, res)
        class(circular_inductor), intent(in) :: self
        integer, intent(in) :: N
        real(dp), intent(in) :: theta(N)
        real(dp), intent(out) :: res(3,N)
        
        res(1,:) = self%radius * cos(theta(:)) + self%center(1)
        res(2,:) = self%radius * sin(theta(:)) + self%center(2)
        res(3,:) = self%center(3)
        call rotate_vector(N, self%normal, res)
    end subroutine circular_get_r

    subroutine circular_get_dl(self, N, theta, res)
        class(circular_inductor), intent(in) :: self
        integer, intent(in) :: N
        real(dp), intent(in) :: theta(N)
        real(dp), intent(out) :: res(3,N)

        res(1,:) = -self%radius * sin(theta(:))
        res(2,:) =  self%radius * cos(theta(:))
        res(3,:) = 0.0_dp
        call rotate_vector(N, self%normal, res)
    end subroutine circular_get_dl

    subroutine circular_get_limits(self, res)
        class(circular_inductor), intent(in) :: self
        real(dp), intent(out) :: res(2)
        res(:) = (/0.0_dp, 2.0_dp*PI/)
    end subroutine circular_get_limits

    ! CIRCULAR GEOMETRICAL MEAN DISTANCES
    
    subroutine circularcs_get_gmd(self, gmd)
        class(circularcs_inductor), intent(in) :: self
        real(dp), intent(out) :: gmd
        gmd = 0.7788_dp * self%radius_section
    end subroutine circularcs_get_gmd

    subroutine circularfs_get_gmd(self, gmd)
        class(circularfs_inductor), intent(in) :: self
        real(dp), intent(out) :: gmd
        gmd = 0.22313 * self%width_section
    end subroutine circularfs_get_gmd

    ! CIRCULAR PARASITIC MODELS
    
    subroutine circularcs_rldrude_get_abcd(self, t, f, abcd)
        class(circularcs_rldrude_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)

        real(dp) :: rho

        call self%material%get_rho(f,t,rho)

        abcd(1) = CONE
        abcd(2) = (2 * self%radius * rho) / (self%radius_section * self%radius_section)
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine circularcs_rldrude_get_abcd

    subroutine circularfs_rldrude_get_abcd(self, t, f, abcd)
        class(circularfs_rldrude_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        real(dp) :: rho

        call self%material%get_rho(f,t,rho)

        abcd(1) = CONE
        abcd(2) = (2 * PI * self%radius * rho) / (self%width_section * self%thickness_section)
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine circularfs_rldrude_get_abcd

    subroutine circularcs_rldowell_get_abcd(self, t, f, abcd)
        class(circularcs_rldowell_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)

        real(dp) :: rho, x

        call self%material%get_rho(f,t,rho)
        
        abcd(1) = CONE
        abcd(3) = CZERO
        abcd(4) = CONE

        x = 1.77_dp * self%radius_section * sqrt(PI * f * MU0 / rho) 
        abcd(2) = (2 * self%radius * rho * x * (sinh(2*x) + sin(2*x))) &
            & / ((cosh(2*x) - cos(2*x)) * self%radius_section * self%radius_section)
    end subroutine circularcs_rldowell_get_abcd

    subroutine circularfs_rldowell_get_abcd(self, t, f, abcd)
        class(circularfs_rldowell_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)

        real(dp) :: rho, x

        call self%material%get_rho(f,t,rho)
        
        abcd(1) = CONE
        abcd(3) = CZERO
        abcd(4) = CONE

        x = self%thickness_section * sqrt(PI * f * MU0 / rho) 
        abcd(2) = (2 * self%radius * rho * x * (sinh(2*x) + sin(2*x))) &
            & / ((cosh(2*x) - cos(2*x)) * self%width_section * self%thickness_section)
    end subroutine circularfs_rldowell_get_abcd

    ! LINEAR SPIRAL CIRCULAR INDUCTOR

    subroutine lcspiral_get_r(self, N, theta, res)
        class(lcspiral_inductor), intent(in) :: self
        integer, intent(in) :: N
        real(dp), intent(in) :: theta(N)
        real(dp), intent(out) :: res(3,N)
        
        res(1,:) = (self%radius - ( self%pitch*theta(:) / (2*PI) )) * cos(theta(:)) + self%center(1)
        res(2,:) = (self%radius - ( self%pitch*theta(:) / (2*PI) )) * sin(theta(:)) + self%center(2)
        res(3,:) = self%center(3)
        call rotate_vector(N, self%normal, res)
    end subroutine lcspiral_get_r

    subroutine lcspiral_get_dl(self, N, theta, res)
        class(lcspiral_inductor), intent(in) :: self
        integer, intent(in) :: N
        real(dp), intent(in) :: theta(N)
        real(dp), intent(out) :: res(3,N)

        res(1,:) = -(self%pitch / (2*PI)) * cos(theta(:)) - (self%radius - ( self%pitch*theta(:) / (2*PI) )) * sin(theta(:))
        res(2,:) =  -(self%pitch / (2*PI)) * sin(theta(:)) + (self%radius - ( self%pitch*theta(:) / (2*PI) )) * cos(theta(:)) 
        res(3,:) = 0.0_dp
        call rotate_vector(N, self%normal, res)
    end subroutine lcspiral_get_dl

    subroutine lcspiral_get_limits(self, res)
        class(lcspiral_inductor), intent(in) :: self
        real(dp), intent(out) :: res(2)
        res(:) = (/0.0_dp, 2*PI*self%turns/)
    end subroutine lcspiral_get_limits

    ! LINEAR SPIRAL GEOMETRICAL MEAN DISTANCES

    subroutine lcspiralcs_get_gmd(self, gmd)
        class(lcspiralcs_inductor), intent(in) :: self
        real(dp), intent(out) :: gmd
        gmd = 0.7788_dp * self%radius_section
    end subroutine lcspiralcs_get_gmd

    subroutine lcspiralfs_get_gmd(self, gmd)
        class(lcspiralfs_inductor), intent(in) :: self
        real(dp), intent(out) :: gmd
        gmd = 0.22313 * self%width_section
    end subroutine lcspiralfs_get_gmd
    
    ! LINEAR SPIRAL PARASITIC MODELS

    subroutine lcspiralcs_rldrude_get_abcd(self, t, f, abcd)
        class(lcspiralcs_rldrude_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)

        real(dp) :: rho

        call self%material%get_rho(f,t,rho)

        abcd(1) = CONE
        abcd(2) = (self%turns * (2*self%radius - self%turns * self%pitch) * rho) &
            & / (self%radius_section * self%radius_section)
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine lcspiralcs_rldrude_get_abcd

    subroutine lcspiralfs_rldrude_get_abcd(self, t, f, abcd)
        class(lcspiralfs_rldrude_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        real(dp) :: rho

        call self%material%get_rho(f,t,rho)

        abcd(1) = CONE
        abcd(2) = (self%turns * (2*self%radius - self%turns * self%pitch) * rho) &
            & / (self%width_section * self%thickness_section)
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine lcspiralfs_rldrude_get_abcd

    subroutine lcspiralcs_rldowell_get_abcd(self, t, f, abcd)
        class(lcspiralcs_rldowell_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)

        real(dp) :: rho, x

        call self%material%get_rho(f,t,rho)
        
        abcd(1) = CONE
        abcd(3) = CZERO
        abcd(4) = CONE

        x = 1.77_dp * self%radius_section * sqrt(self%radius_section * PI * f * MU0 / (rho * self%pitch)) 
        abcd(2) = (self%turns * (2*self%radius - self%turns * self%pitch) * rho * x * (sinh(2*x) + sin(2*x))) &
            & / ((cosh(2*x) - cos(2*x)) * self%radius_section * self%radius_section)
    end subroutine lcspiralcs_rldowell_get_abcd

    subroutine lcspiralfs_rldowell_get_abcd(self, t, f, abcd)
        class(lcspiralfs_rldowell_inductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)

        real(dp) :: rho, x

        call self%material%get_rho(f,t,rho)
        
        abcd(1) = CONE
        abcd(3) = CZERO
        abcd(4) = CONE

        x = self%thickness_section * sqrt(self%width_section * PI * f * MU0 / (rho * self%pitch))
        abcd(2) = (self%turns * (2*self%radius - self%turns * self%pitch) * rho * x * (sinh(2*x) + sin(2*x))) &
            & / ((cosh(2*x) - cos(2*x)) * self%width_section * self%thickness_section)
    end subroutine lcspiralfs_rldowell_get_abcd

end module circular_inductors_library
