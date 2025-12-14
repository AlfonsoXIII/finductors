module utils_library
    use types
    implicit none
    
    private
    public :: generate_xzmesh

contains
    
    subroutine generate_xzmesh(npoints, xmin, xmax, zmin, zmax, yval, xyz)
        integer, intent(in) :: npoints
        real(dp), intent(in) :: xmin, xmax, zmin, zmax, yval
        real(dp), intent(out) :: xyz(3,npoints**2)

        integer :: idx, jdx, k
        real(dp) :: dx, dz

        dx = (xmax - xmin) / (npoints - 1)
        dz = (zmax - zmin) / (npoints - 1)
        
        do idx=1,npoints
            !$omp simd
            do jdx=1,npoints
                k = (idx-1)*npoints + jdx
                xyz(1,k) = xmin + (idx-1) * dx
                xyz(2,k) = yval
                xyz(3,k) = zmin + (jdx-1) * dz
            end do
        end do
    end subroutine
end module
