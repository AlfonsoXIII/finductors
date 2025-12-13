module components_library
    use types
    use materials_library
    implicit none

    private
    public :: component_t, sresistor, presistor, scapacitor, pcapacitor, &
        & sinductor, pinductor
   
    type, abstract :: component_t
    contains         
        procedure(comp_get_abcd), deferred :: get_abcd
    end type component_t
    
    abstract interface
        subroutine comp_get_abcd(self, t, f, abcd)
            import :: dp, component_t
            class(component_t), intent(in) :: self
            real(dp), intent(in) :: t, f
            complex(dp), intent(out) :: abcd(4)
        end subroutine comp_get_abcd
    end interface

    ! PASSIVE COMPONENTS

    type, extends(component_t) :: sresistor
    private
        real(dp) :: resistance
    contains
        procedure :: set_resistance => sresistor_set_resistance
        procedure :: get_abcd => sresistor_get_abcd
    end type sresistor

    type, extends(component_t) :: presistor
    private
        real(dp) :: resistance
    contains
        procedure :: set_resistance => presistor_set_resistance
        procedure :: get_abcd => presistor_get_abcd
    end type presistor

    type, extends(component_t) :: scapacitor
    private
        real(dp) :: capacitance 
    contains
        procedure :: set_capacitance => scapacitor_set_capacitance
        procedure :: get_abcd => scapacitor_get_abcd
    end type scapacitor

    type, extends(component_t) :: pcapacitor
    private
        real(dp) :: capacitance
    contains
        procedure :: set_capacitance => pcapacitor_set_capacitance
        procedure :: get_abcd => pcapacitor_get_abcd
    end type pcapacitor

    type, extends(component_t) :: sinductor
    private
        real(dp) :: inductance
    contains
        procedure :: set_inductance => sinductor_set_inductance
        procedure :: get_abcd => sinductor_get_abcd
    end type sinductor

    type, extends(component_t) :: pinductor
    private
        real(dp) :: inductance
    contains
        procedure :: set_inductance => pinductor_set_inductance
        procedure :: get_abcd => pinductor_get_abcd
    end type pinductor

contains

    ! COMPONENTS

    ! series resistor

    subroutine sresistor_set_resistance(self, resistance)
        class(sresistor), intent(inout) :: self
        real(dp), intent(in) :: resistance

        self%resistance = resistance
    end subroutine

    subroutine sresistor_get_abcd(self, t, f, abcd)
        class(sresistor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = self%resistance
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine
 
    ! parallel resistor

    subroutine presistor_set_resistance(self, resistance)
        class(presistor), intent(inout) :: self
        real(dp), intent(in) :: resistance

        self%resistance = resistance
    end subroutine

    subroutine presistor_get_abcd(self, t, f, abcd)
        class(presistor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = CZERO
        abcd(3) = self%resistance
        abcd(4) = CONE
    end subroutine

    ! series capacitor
    
    subroutine scapacitor_set_capacitance(self, capacitance)
        class(scapacitor), intent(inout) :: self
        real(dp), intent(in) :: capacitance

        self%capacitance = capacitance
    end subroutine

    subroutine scapacitor_get_abcd(self, t, f, abcd)
        class(scapacitor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = 1.0_dp / cmplx(0.0_dp,2*PI*f*self%capacitance)
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine
    
    ! parallel capacitor

    subroutine pcapacitor_set_capacitance(self, capacitance)
        class(pcapacitor), intent(inout) :: self
        real(dp), intent(in) :: capacitance

        self%capacitance = capacitance
    end subroutine
 
    subroutine pcapacitor_get_abcd(self, t, f, abcd)
        class(pcapacitor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = CZERO
        abcd(3) = cmplx(0.0_dp,2*PI*f*self%capacitance)
        abcd(4) = CONE
    end subroutine

    ! series inductor
    
    subroutine sinductor_set_inductance(self, inductance)
        class(sinductor), intent(inout) :: self
        real(dp), intent(in) :: inductance

        self%inductance = inductance
    end subroutine

    subroutine sinductor_get_abcd(self, t, f, abcd)
        class(sinductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = cmplx(0.0_dp,2*PI*f*self%inductance)
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine
 
    ! parallel inductor

    subroutine pinductor_set_inductance(self, inductance)
        class(pinductor), intent(inout) :: self
        real(dp), intent(in) :: inductance

        self%inductance = inductance
    end subroutine

    subroutine pinductor_get_abcd(self, t, f, abcd)
        class(pinductor), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = CZERO
        abcd(3) = 1.0_dp / cmplx(0.0_dp,2*PI*f*self%inductance)
        abcd(4) = CONE
    end subroutine

end module components_library
