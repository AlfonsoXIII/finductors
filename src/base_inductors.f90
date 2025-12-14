module base_inductors_library
    use types
    use materials_library
    implicit none

    type, abstract :: inductor_t
        real(dp) :: center(3), normal(3)
        class(material_t), pointer :: material
    contains
        procedure(r_sub), deferred :: get_r
        procedure(dl_sub), deferred :: get_dl
        procedure(limits_sub), deferred :: get_limits
        procedure(gmd_sub), deferred :: get_gmd
        procedure(abcd_sub), deferred :: get_abcd
    end type inductor_t

    abstract interface
        subroutine r_sub(self, N, theta, res)
            import :: dp, inductor_t
            class(inductor_t), intent(in) :: self
            integer, intent(in) :: N
            real(dp), intent(in) :: theta(N)
            real(dp), intent(out) :: res(3,N)
        end subroutine r_sub

        subroutine dl_sub(self, N, theta, res)
            import :: dp, inductor_t
            class(inductor_t), intent(in) :: self
            integer, intent(in) :: N
            real(dp), intent(in) :: theta(N)
            real(dp), intent(out) :: res(3,N)
        end subroutine dl_sub

        subroutine limits_sub(self, res)
            import :: dp, inductor_t
            class(inductor_t), intent(in) :: self
            real(dp), intent(out) :: res(2)
        end subroutine limits_sub

        subroutine gmd_sub(self, gmd)
            import :: dp, inductor_t
            class(inductor_t), intent(in) :: self
            real(dp), intent(out) :: gmd
        end subroutine

        subroutine abcd_sub(self, t, f, abcd)
            import :: dp, inductor_t
            class(inductor_t), intent(in) :: self
            real(dp), intent(in) :: t, f
            complex(dp), intent(out) :: abcd(4)
        end subroutine
    end interface

contains

    subroutine rotate_vector(N, normal, vector)
        integer, intent(in) :: N
        real(dp), intent(in) :: normal(3)
        real(dp), intent(inout) :: vector(3,N)

        real(dp), allocatable :: buffer(:,:)
        real(dp) :: axis(2), theta, costheta, sintheta
        
        if (abs(normal(3) - 1.0_dp) < 1e-12_dp) return       
        
        allocate(buffer(3,N))
        buffer = vector
        axis = (/-normal(2), normal(1)/) 
        axis = axis / sqrt(axis(1)**2 + axis(2)**2) 

        theta = acos(normal(3))
        costheta = cos(theta)
        sintheta = sin(theta)

        vector(1,:) = (costheta + axis(1)**2 * (1 - costheta)) * buffer(1,:) + & 
                      axis(1)*axis(2)*(1-costheta) * buffer(2,:) + & 
                      axis(2)*sintheta * buffer(3,:)
        vector(2,:) = axis(2)*axis(1)*(1-costheta) * buffer(1,:) + &
                      (costheta + axis(2)**2 * (1-costheta)) * buffer(2,:) - &
                      axis(1)*sintheta * buffer(3,:)
        vector(3,:) = -axis(2)*sintheta * buffer(1,:) + &
                      axis(1)*sintheta * buffer(2,:) + &
                      costheta * buffer(3,:)
        deallocate(buffer)
    end subroutine rotate_vector
end module base_inductors_library
