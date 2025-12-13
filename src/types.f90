module types
    implicit none
    
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: MU0 = 4*PI*1.0E-7_dp

    complex(dp), parameter :: CONE = cmplx(1.0_dp, 0.0_dp, kind=dp)
    complex(dp), parameter :: CZERO = cmplx(0.0_dp, 0.0_dp, kind=dp)
    complex(dp), parameter :: CJ = cmplx(0.0_dp,1.0_dp, kind=dp)

end module types
