module materials_library
    use types
    implicit none
    
    ! BASE MATERIAL

    type, abstract :: material_t
    contains         
        procedure(material_get_rho), deferred :: get_rho
        procedure(material_get_mu), deferred :: get_mu
    end type material_t

    abstract interface
        subroutine material_get_rho(self, f, t, rho)
            import :: dp, material_t
            class(material_t), intent(in) :: self
            real(dp), intent(in) :: f, t
            real(dp), intent(out) :: rho
        end subroutine material_get_rho

        subroutine material_get_mu(self, f, t, mu)
            import :: dp, material_t
            class(material_t), intent(in) :: self
            real(dp), intent(in) :: f, t
            real(dp), intent(out) :: mu
        end subroutine material_get_mu

    end interface

    ! MATERIALS

    type, extends(material_t) :: Cu
    contains
        procedure :: get_rho => cu_get_rho
        procedure :: get_mu => cu_get_mu
    end type Cu

contains

    subroutine cu_get_rho(self, f, t, rho)
        class(Cu), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        
        rho = 1.72E-8_dp * (1.0_dp + 0.00393_dp * (t - 20.0_dp))
    end subroutine cu_get_rho

    subroutine cu_get_mu(self, f, t, mu)
        class(Cu), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        
        mu = 1.0_dp
    end subroutine cu_get_mu

end module materials_library
