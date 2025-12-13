module sources_library
    use types
    use components_library
    implicit none

    type :: harmonic_t
        complex(dp) :: amplitude
        real(dp) :: frequency
    end type

    type, extends(component_t), abstract :: source_t
        integer :: n_harmonics = 1
        type(harmonic_t), allocatable :: harmonics(:)
        real(dp) :: resistance
    contains
        procedure :: set_resistance => source_set_resistance
        procedure :: get_abcd => source_get_abcd
    end type source_t

    type, extends(source_t) :: sinusoidal_source
    contains
        procedure :: configure_source => sinusoidal_configure_source
    end type sinusoidal_source

    type, extends(source_t) :: pulsed_source
    contains
        procedure :: set_harmonics => pulsed_set_harmonics
        procedure :: configure_source => pulsed_configure_source
    end type pulsed_source

    type, extends(source_t) :: sawtooth_source
    contains
        procedure :: set_harmonics => sawtooth_set_harmonics
        procedure :: configure_source => sawtooth_configure_source
    end type sawtooth_source

contains

    subroutine source_set_resistance(self, new_resistance)
        class(source_t), intent(inout) :: self
        real(dp), intent(in) :: new_resistance
        self%resistance = new_resistance
    end subroutine

    subroutine source_get_abcd(self, t, f, abcd)
        class(source_t), intent(in) :: self
        real(dp), intent(in) :: t, f
        complex(dp), intent(out) :: abcd(4)
        
        abcd(1) = CONE
        abcd(2) = self%resistance
        abcd(3) = CZERO
        abcd(4) = CONE
    end subroutine

    ! SINUSOIDAL SOURCE

    subroutine sinusoidal_configure_source(self, amplitude, frequency)
        class(sinusoidal_source), intent(inout) :: self
        complex(dp), intent(in) :: amplitude
        real(dp), intent(in) :: frequency

        if (.not. allocated(self%harmonics)) allocate(self%harmonics(1))

        self%harmonics(1)%amplitude = amplitude
        self%harmonics(1)%frequency = frequency
 
    end subroutine sinusoidal_configure_source

    ! PULSED SOURCE

    subroutine pulsed_configure_source(self, amplitude, frequency, duty_cicle)
        class(pulsed_source), intent(inout) :: self
        complex(dp), intent(in) :: amplitude
        real(dp), intent(in) :: frequency, duty_cicle

        integer :: idx

        if (.not. allocated(self%harmonics)) then
            allocate(self%harmonics(self%n_harmonics))
        else if (size(self%harmonics) /= self%n_harmonics) then
            deallocate(self%harmonics)
            allocate(self%harmonics(self%n_harmonics))
        end if
        
        !$omp simd
        do idx=1,self%n_harmonics
            self%harmonics(idx)%amplitude = (amplitude * sin(idx*PI*duty_cicle) * &
                & exp(-CJ*(idx*PI*duty_cicle))) / (idx*PI)
            self%harmonics(idx)%frequency = idx * frequency
        end do
    end subroutine pulsed_configure_source

    subroutine pulsed_set_harmonics(self, n_harmonics)
        class(pulsed_source), intent(inout) :: self
        integer, intent(in) :: n_harmonics
        
        self%n_harmonics = n_harmonics
    end subroutine pulsed_set_harmonics

    ! SAWTOOTH SOURCE

    subroutine sawtooth_configure_source(self, amplitude, frequency, duty_cicle)
        class(sawtooth_source), intent(inout) :: self
        complex(dp), intent(in) :: amplitude
        real(dp), intent(in) :: frequency, duty_cicle
        
        integer :: idx

        if (.not. allocated(self%harmonics)) then
            allocate(self%harmonics(self%n_harmonics))
        else if (size(self%harmonics) /= self%n_harmonics) then
            deallocate(self%harmonics)
            allocate(self%harmonics(self%n_harmonics))
        end if
 
        !$omp simd
        do idx=1,self%n_harmonics
            self%harmonics(idx)%amplitude = amplitude * (exp(-CJ*2*PI*idx*duty_cicle) - 1) &
                & / (2*PI*PI*idx*idx*duty_cicle*(1-duty_cicle))
            self%harmonics(idx)%frequency = idx * frequency
        end do
    end subroutine sawtooth_configure_source

    subroutine sawtooth_set_harmonics(self, n_harmonics)
        class(sawtooth_source), intent(inout) :: self
        integer, intent(in) :: n_harmonics

        self%n_harmonics = n_harmonics
    end subroutine sawtooth_set_harmonics

end module sources_library
