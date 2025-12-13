module math_library
    use types
    implicit none
    
    private
    public hammersley2d_sequence, hsa_quicksort, cmmp

contains

    ! MATRIX ARITHMETIC

    pure subroutine cmmp(a, b)
        complex(dp), intent(inout) :: a(4)
        complex(dp), intent(in) :: b(4)
    
        complex(dp) :: c11, c12, c21, c22

        c11 = a(1)*b(1) + a(2)*b(3)
        c12 = a(1)*b(2) + a(2)*b(4)
        c21 = a(3)*b(1) + a(4)*b(3)
        c22 = a(3)*b(2) + a(4)*b(4)

        a(1) = c11
        a(2) = c12
        a(3) = c21
        a(4) = c22
    end subroutine
    
    ! HAMMERSLEY SEQUENCE

    pure subroutine radical_inverse(N, base, samples)
        integer, intent(in) :: N, base
        real(dp), intent(out) :: samples(N)

        real(dp) :: phi, f
        integer :: idx, i

        do i=0,N-1
            phi = 0.0_dp
            f = 1.0_dp / real(base, kind=dp)
            idx = i

            do while (idx > 0)
                phi = phi + mod(idx, base)* f
                idx = idx / base
                f = f / real(base, kind=dp)
            end do
            samples(i+1) = phi
        end do
    end subroutine

    pure subroutine hammersley2d_sequence(N, samples)
        integer, intent(in) :: N
        real(dp), intent(out) :: samples(N,2)

        call radical_inverse(N, 2, samples(:,1))
        call radical_inverse(N, 3, samples(:,2))
    end subroutine
    
    ! QUICKSORT

    pure subroutine hsa_quicksort(n_elements, harmonics, sources, amplitudes, low, high)
        integer, intent(in) :: n_elements, low, high
        real(dp), intent(inout) :: harmonics(n_elements)
        integer, intent(inout) :: sources(n_elements)
        complex(dp), intent(inout) :: amplitudes(n_elements)

        integer :: pivot, idx, jdx, source_buffer
        real(dp) :: pivot_eval, harmonic_buffer
        complex(dp) :: amplitude_buffer

        if (low > high) then
            return
        end if

        ! pivot
        
        pivot = low + (high - low) / 2 

        ! partition
        
        idx = low
        jdx = high
        pivot_eval = harmonics(pivot)
        do while(idx <= jdx)
            
            do while(harmonics(idx) < pivot_eval)
                idx = idx + 1
            end do
            
            do while(harmonics(jdx) > pivot_eval)
                jdx = jdx - 1
            end do

            if (idx <= jdx) then
                harmonic_buffer = harmonics(idx)
                harmonics(idx) = harmonics(jdx)
                harmonics(jdx) = harmonic_buffer

                source_buffer = sources(idx)
                sources(idx) = sources(jdx)
                sources(jdx) = source_buffer

                amplitude_buffer = amplitudes(idx)
                amplitudes(idx) = amplitudes(jdx)
                amplitudes(jdx) = amplitude_buffer
                
                idx = idx + 1
                jdx = jdx - 1

            end if
        end do

        if (low < jdx) call hsa_quicksort(n_elements, harmonics, sources, amplitudes, low, jdx)
        if (idx < high) call hsa_quicksort(n_elements, harmonics, sources, amplitudes, idx, high)
    end subroutine

end module math_library
