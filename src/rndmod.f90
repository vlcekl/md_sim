!===============================================================
!   MODULES & FUNCTIONS INVOLVED IN RANDOM NUMBER GENERATION  
!===============================================================

!
!====== Function prototypes of random generators
!
module rndmod

    real*8 :: seed

    contains

!
!====== Initializes random seed
!
subroutine rnd_init(start)
    implicit none
    real*8 :: start
    seed = start
    return
end subroutine rnd_init


!
!===== Congruential random nuber generator <0; 1) =====!
!
function rnd(dummy)
    implicit none
    real*8 :: rnd, dummy

    seed = mod(0.8192d4*seed, 0.67099547d8)
    rnd = seed*1.4903230270690d-8
    return
end function rnd


!
!===== Generator of random numbers with the Gaussian distribution =====!
!
function rndgauss(dummy)

    implicit none
    real*8, parameter :: M_PI = 3.1415926535898  !/* pie */
    real*8 :: rndgauss, dummy
    real*8, SAVE :: ra, rb
    integer*4, SAVE :: io = 0

    if (io == 1) then
        io = 0
        rndgauss = ra*sin(rb)
        return
    end if
    io = 1

    ra = sqrt(-2.0*log(rnd(dummy)));
    rb = 2.*M_PI*rnd(dummy);

    rndgauss = ra*cos(rb)
    return

end function rndgauss

end module rndmod
