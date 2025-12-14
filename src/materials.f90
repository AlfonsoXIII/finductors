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

    type, extends(material_t) :: Cu_Material
    contains
        procedure :: get_rho => cu_get_rho
        procedure :: get_mu => cu_get_mu
    end type

    type, extends(material_t) :: Cu_ETP_Material
    contains
        procedure :: get_rho => cu_etp_get_rho
        procedure :: get_mu => cu_etp_get_mu
    end type

    type, extends(material_t) :: Cu_OF_Material
    contains
        procedure :: get_rho => cu_of_get_rho
        procedure :: get_mu => cu_of_get_mu
    end type

    type, extends(material_t) :: Ag_Material
    contains
        procedure :: get_rho => ag_get_rho
        procedure :: get_mu => ag_get_mu
    end type

    type, extends(material_t) :: Ag_6061_Material
    contains
        procedure :: get_rho => ag_6061_get_rho
        procedure :: get_mu => ag_6061_get_mu
    end type

    type, extends(material_t) :: Ag_1350_Material
    contains
        procedure :: get_rho => ag_1350_get_rho
        procedure :: get_mu => ag_1350_get_mu
    end type

    type, extends(material_t) :: Au_Material
    contains
        procedure :: get_rho => au_get_rho
        procedure :: get_mu => au_get_mu
    end type

contains

    ! PURE Cu

    subroutine cu_get_rho(self, f, t, rho)
        class(Cu_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 1.68E-8_dp * (1.0_dp + 0.00404_dp * (t - 20.0_dp))
    end subroutine cu_get_rho

    subroutine cu_get_mu(self, f, t, mu)
        class(Cu_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp
    end subroutine cu_get_mu

    ! Annealed Cu

    subroutine cu_etp_get_rho(self, f, t, rho)
        class(Cu_ETP_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 1.72E-8_dp * (1.0_dp + 0.00393_dp * (t - 20.0_dp))
    end subroutine

    subroutine cu_etp_get_mu(self, f, t, mu)
        class(Cu_ETP_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp 
    end subroutine
    
    ! Oxygen Free Cu

    subroutine cu_of_get_rho(self, f, t, rho)
        class(Cu_OF_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 1.707E-8_dp * (1.0_dp + 0.00400_dp * (t - 20.0_dp))
    end subroutine

    subroutine cu_of_get_mu(self, f, t, mu)
        class(Cu_OF_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp 
    end subroutine

    ! PURE Ag

    subroutine ag_get_rho(self, f, t, rho)
        class(Ag_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 1.59E-8_dp * (1.0_dp + 0.0038_dp * (t - 20.0_dp))
    end subroutine

    subroutine ag_get_mu(self, f, t, mu)
        class(Ag_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp 
    end subroutine

    ! Ag 1350

    subroutine ag_1350_get_rho(self, f, t, rho)
        class(Ag_1350_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 2.82E-8_dp * (1.0_dp + 0.00403_dp * (t - 20.0_dp))
    end subroutine

    subroutine ag_1350_get_mu(self, f, t, mu)
        class(Ag_1350_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp 
    end subroutine

    ! Ag 6061

    subroutine ag_6061_get_rho(self, f, t, rho)
        class(Ag_6061_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 4.00E-8_dp * (1.0_dp + 0.0038_dp * (t - 20.0_dp))
    end subroutine

    subroutine ag_6061_get_mu(self, f, t, mu)
        class(Ag_6061_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp 
    end subroutine

    ! PURE Au

    subroutine au_get_rho(self, f, t, rho)
        class(Au_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: rho
        rho = 2.44E-8_dp * (1.0_dp + 0.0034_dp * (t - 20.0_dp))
    end subroutine

    subroutine au_get_mu(self, f, t, mu)
        class(Au_Material), intent(in) :: self
        real(dp), intent(in) :: f, t
        real(dp), intent(out) :: mu
        mu = 1.0_dp 
    end subroutine

end module materials_library
