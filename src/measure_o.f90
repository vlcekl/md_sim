module measure
    use config
    use rndmod
    implicit none

    ! files
    character*40 :: filmeas
    ! time
    real*8 :: timprn
    integer*4 :: itprn
    ! length parameters 
    !real*8 :: rdel, rdeli
    real*8 :: rdel, rdeli, rfrst, rscnd
    integer*8 :: lmax

    real*8, dimension(:), allocatable :: sumx, sumy
    ! correaltions functions
    integer*8, dimension(:, :), allocatable :: ngc, ngang, ngz

    ! diffusion
    logical :: l2d, ldiff, ldiff_r, lvisc, lhb
    real*8, dimension(:, :), allocatable :: diffx, diffy, diffz
    real*8, dimension(:), allocatable :: xo, yo, zo
    integer*4 :: igmin, igmax

    integer*8 :: rbin
    real*8 :: diftime, diftimax, timdif
    real*8 :: rishell 

    real*8, dimension(-1:1) :: fwxmeas, fwymeas, fwzmeas, ftotx, ftoty, ftotz

    ! chemical potential
    real*8, dimension(3, 3) :: rins
    real*8 :: lwmin, lwmax, facbin, chtemper
    integer*8 :: ichmax, lchmax, ilmid
    real*8, dimension(:), allocatable :: pave, pavesum, xarr, yarr, zarr
    integer*8, dimension(:), allocatable :: inave, inavesum, iarr
    contains

    subroutine measure_zero
        implicit none

        ngc = 0
        ngang = 0

        ngz = 0

        sumx = 0.0
        sumy = 0.0

        if (.true.) then
            xo(igmin:igmax) = x0(igmin:igmax)
            yo(igmin:igmax) = y0(igmin:igmax)
            zo(igmin:igmax) = z0(igmin:igmax)
            diffx = 0.0 ; diffy = 0.0 ; diffz = 0.0
        end if

        return
    end subroutine measure_zero

    subroutine measure_init
        implicit none
        integer*4 :: ig

        itprn = nint(timprn/dt)  ! interval for measurement

        ig = nmove/nproc
        igmin = ig*me + 1
        igmax = ig*(me + 1)
        if (me == nproc-1) igmax = nmove

        rdeli = 1.0/rdel
        lmax = int((rcut+2.0)*rdeli + 0.5) + 1
        allocate(ngc(nijt, 0:lmax), ngang(nit, 0:201))
        allocate(ngz(nit, 0:lmax))
        allocate(sumx(nit), sumy(nit))

        if (.true.) then
            allocate(xo(igmin:igmax), yo(igmin:igmax), zo(igmin:igmax))
            allocate(diffx(itrun, nit), diffy(itrun, nit), diffz(itrun, nit))
            diffx = 0.0 ; diffy = 0.0 ; diffz = 0.0
        end if

        call measure_zero

        return
    end subroutine measure_init

    subroutine measure_do
        implicit none
        integer*4 :: k, l, i, it
        real*8 :: xx, yy
        integer*8, save :: ii = 0

        ftotx = ftotx + fwxmeas
        ftoty = ftoty + fwymeas
        ftotz = ftotz + fwzmeas

        ii = ii + 1

        do i = 1, nmove
            it = itype(i)
            xx = x0(i)
            yy = y0(i)
            !z0(i) = z0(i) - lz*floor(z0(i)*lzi+0.5)
            l = int(sqrt(xx**2 + yy**2)*rdeli) + 1
            if (l <= lmax) ngz(it, l) = ngz(it, l) + 1
            sumx(it) = sumx(it) + x0(i)
            sumy(it) = sumy(it) + y0(i)
        end do
!        do i = 1, 8
!            print *, i, sumx(i)/real(ii), sumy(i)/real(ii)
!        end do

        if (.false.) then
            do i = igmin, igmax
                it = itype(i)
                diffx(ii, it) = diffx(ii, it) + (x0(i) - xo(i))**2
                diffy(ii, it) = diffy(ii, it) + (y0(i) - yo(i))**2
                diffz(ii, it) = diffz(ii, it) + (z0(i) - zo(i))**2
            end do
        end if

        return
    end subroutine measure_do

    subroutine measure_output
        implicit none
        integer*4 :: l, it, ii
        real*8 :: r, fac, facv, rh2o

        if (.false.) then
            if (nproc > 1) then
                call measure_gather_dble2d(diffx)
                call measure_gather_dble2d(diffy)
                call measure_gather_dble2d(diffz)
            end if
            do ii = 1, itrun
                do it = 1, nit
                    diffx(ii, it) = diffx(ii, it)/(2.0*dt*real(ii)*nt(it))
                    diffy(ii, it) = diffy(ii, it)/(2.0*dt*real(ii)*nt(it))
                    diffz(ii, it) = diffz(ii, it)/(2.0*dt*real(ii)*nt(it))
                end do
            end do
        end if

        if (nproc > 1) then
            call measure_gather_int2d(ngc)
            call measure_gather_int2d(ngang)
            call measure_gather_int2d(ngz)
        end if

        if (me == 0) then
            open(4, file = filmeas, status='unknown')

            fac = 3*lx*ly*lz/(2*M_PI*real(itrun)) !*real(nmol*nmol))

            do l = 1, lmax
                r = (real(l)+0.5)*rdel
                facv = fac/((r + rdel)**3 - r**3)
                write(4, 10) r, (facv*real(ngc(it, l)/(ntota*ntota/2.)), it = 1, nijt)
            end do

            write(4, *) '# angle'
            fac = 1.0/real(itrun)
            do l = 1, 200
                r = real(l-1)/100 - 1.0
                write(4, *) r, (fac*ngang(it, l), it = 1, ntype)
            end do

            write(4, *) '# crossection'
            fac = 3*lx*ly*lz/(2*M_PI*real(itrun)) !*real(nmol*nmol))
            do l = 1, lmax
                r = (real(l)+0.5)*rdel
                facv = fac/((r + rdel)**2 - r**2)
                write(4, 10) r, (facv*real(ngz(it, l)/(ntota*ntota/2.)), it = 1, nit)
            end do

            ! Diffusion
            if (.false.) then
                write(4, *) '# Diffusion'
                fac = (1e-8)**2/(1e-15)*1e5  
                diffx = fac*(diffx + diffy)
                diffz = fac*(diffz)
                do ii = 1, itrun
                    write(4, *) (diffx(ii, it), diffz(ii, it), it=1, nit)
                end do 
            end if

            close(4, status='keep')
        end if

        call measure_deallocate

        return
10      format(f7.3,5(1x,f11.6))
!60      format(f5.2,6(1x,f11.6))
    end subroutine measure_output

    subroutine measure_gather_int2d(iarr)
        implicit none
        integer*8 :: ierr, siz
        integer*8, dimension(:, :) :: iarr
        integer*8, dimension(:), allocatable :: ibufi, ibufo
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz), ibufo(siz))
        ibufi = 0 ; ibufo = 0
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_allreduce(ibufi, ibufo, siz, mpi_integer8, mpi_sum, mpi_comm_world, ierr)
        iarr = reshape(ibufo, shape(iarr))
        deallocate(ibufi, ibufo)

        return
    end subroutine measure_gather_int2d

    subroutine measure_gather_dble2d(iarr)
        implicit none
        integer*8 :: ierr, siz
        real*8, dimension(:, :) :: iarr
        real*8, dimension(:), allocatable :: ibufi, ibufo
        include 'mpif.h'
        
        siz = size(iarr)
        allocate(ibufi(siz), ibufo(siz))
        ibufi = 0 ; ibufo = 0
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_allreduce(ibufi, ibufo, siz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        iarr = reshape(ibufo, shape(iarr))
        deallocate(ibufi, ibufo)

        return
    end subroutine measure_gather_dble2d


    subroutine measure_deallocate
        implicit none

        deallocate(ngc, ngang)
        deallocate(xo, yo, zo, diffx, diffy, diffz)

        return
    end subroutine measure_deallocate

end module measure
