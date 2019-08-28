module global          ! global parameters of simulation
    implicit none

    ! Limits
    integer*4, parameter :: TMX = 12             ! max. # of atom types
    integer*4, parameter :: RNMX = 2048          ! max. # of rigid molecules
    integer*4, parameter :: SMX = 3              ! max. # of sites per rigid mol
    integer*4, parameter :: TNMI = (RNMX-200)*SMX ! min. # of rigid atoms
    integer*4, parameter :: FNMX = 650           ! max. # of moving flexible atoms
    integer*4, parameter :: MNMX = RNMX*SMX+FNMX ! max. # of atoms
    integer*4, parameter :: WNMX = 3456          ! max. # of fixed atoms
    integer*4, parameter :: NMX = MNMX+WNMX      ! max. # of atoms
    integer*4, parameter :: BMX = FNMX           ! max. # of intramol. bonds
    integer*4, parameter :: KMX = 5              ! max. # of k-vectors
    integer*4, parameter :: LMX = 5000           ! max. # of histogram windows

    ! Mathematical and physical constants (copied from unitsmod.f90)
    real*8, parameter ::  M_PI      = 3.1415926535898    ! pie
    real*8, parameter ::  PH_e      = 1.6021892E-19      ! elementary charge
    real*8, parameter ::  PH_NA     = 6.022045E23        ! Avogadro's number
    real*8, parameter ::  PH_R      = 8.31441                  !/* universal gas constant */
    real*8, parameter ::  PH_eps0   = 8.85418782E-12           !/* vacuum permitivity */
    real*8, parameter ::  PH_4epi   = 8.987551784952e+9  ! 1./(4.*PI*eps0)
    real*8, parameter ::  PH_A2m    = 1.0E-10            ! Angstrom to SI: 1 A = 1.0E-10 m

    ! MPI variables
    integer*4 :: nproc, me, nodeproc                           ! number of processors
end module global

