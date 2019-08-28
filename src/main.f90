!===========================================================================+
!

!        MD simulation of water - metal oxide surface
!        EW3DC summation used for long range interactions
!        Simulation in NVT ensemble with Berendsen thermostat
!       
!        Written in Fortran 90/95 (free format)
!
!        L. Vlcek (November 2004 - June 2005)
!
!===========================================================================+


!===========================================================================+
!                               Main program                                |
!===========================================================================+

program main
    use move
    implicit none
    character*40 :: fil, fil2, fil3
    integer :: ierr

    include 'mpif.h'
    
    nproc = 1
    me = 0
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)
    call mpi_comm_rank(mpi_comm_world, me, ierr)

    call input
    call init

    ! equilibration
    uavextot = 0.0
    do itime = 1, iteq                                   ! time steps
!        print *, 'a'
        call move_predict                                  ! predictor
!        print *, 'b'
        call force_get                                    ! compute forces
!        print *, 'c'
        call move_correct                                  ! move particles
!        print *, 'd'
        if (mod(itime, itprn) == 0) call prntout(itime)     ! runtime print
!        print *, 'e'
    end do

!    fil = 'eqil.cfg'
!    call config_dump(fil)

    call measure_zero
    call force_zero

    ! production
    fil3 = filcfo
    do itime = 1, itrun                                  ! time steps
!        print *, 'a'
        call move_predict                                  ! predictor
!        print *, 'b'
        call force_get                                    ! compute forces
!        print *, 'c'
        call move_correct                                  ! move particles
!        print *, 'd'
        call measure_do                                     ! measurement
!        print *, 'e'
        if (mod(itime, itprn) == 0) call prntout(itime)     ! runtime print
!        print *, 'f'
        if (mod(itime, 1000000) == 0) then !save configurations
            write(fil, '(i4)') itime/1000000
            fil2 = adjustl(fil)
            ierr = len_trim(fil2)
            fil2(ierr+1:ierr+6) = 'xx.cfg'
            filcfo = fil2
            if (me == 0) then
                call config_output
            end if
        end if
!
        if (mod(itime, 2000) == 0) then !save configurations
            if (me == 0) then
                write(fil, '(i5)') itime/2000
                fil2 = adjustl(fil)
                ierr = len_trim(fil2)
                fil2(ierr+1:ierr+4) = '.cfg'
                call config_wdump(fil2)
            end if
        end if

!        if (mod(itime, 1000) == 0) then !save configurations
!            if (me == 0) then
!                write(fil, '(i5)') itime/10000
!                fil2 = adjustl(fil)
!                ierr = len_trim(fil2)
!                fil2(ierr+1:ierr+4) = '.cfg'
!                call measure_idump(fil2)
!            end if
!        end if
    end do
    filcfo = fil3

    call output

    call mpi_finalize(ierr)

    stop
end program main

!
!======================= Reads input data ====================================
!
subroutine input
    use move                ! rcut and similar
    use comm
    implicit none
    integer*4 :: it
    character*3 :: enstype
    character*40 :: str, filname

    if (me == 0) then
        open(1, file='tot.inp', status='old')
        read(1, *) str, enstype                         ! temperature
        if (enstype .eq. 'nvt') then
            read(1, *) str, temper                         ! temperature
        else
            temper = -1.0
        end if
        read(1, *) str, timeq, timrun, dt   ! number of cycles
        read(1, *) str, timprn              ! print interval

        read(1, *) str, rcut          ! ewald params
        read(1, *) str, kmax
        read(1, *) str, skin          ! width of the skin of the nnlist

        read(1, *) str, rdel
        read(1, *) str, l2d    ! 2d distribution
        read(1, *) str, lhb
        if (l2d .or. lhb) read(1, *) str, rfrst, rscnd
        if (lhb) read(1, *) str, rishell
            read(1, *) str, ldiff  ! diffusion
        if (ldiff) read(1, *) str, rbin, timdif, diftimax
        !print *, rbin, timdif
            read(1, *) str, lcouett ! Couette flow - move one wall
        if (lcouett) read(1, *) str, dshx, dshy, dshz
            read(1, *) str, lfext   ! Pouisseuille flow - push particles
        if (lfext) read(1, *) str, fxext, fyext, fzext
            read(1, *) str, lchempot
        if (lchempot) then
            read(1, *) str, ichmax
            read(1, *) str, lwmin, lwmax, lchmax, ilmid ! minimal z-coord, maximal z-coord, number of bins
        end if

        read(1, *) str, filname              ! output configuration
        close(1, status='keep')
    end if

    print *, 'l2d, ldiff, lcouett, lfext, lhb'
    print *, l2d, ldiff, lcouett, lfext, lhb

!    read *, str, enstype                         ! temperature
!    if (enstype .eq. 'nvt') then
!        read *, str, temper                         ! temperature
!    else
!        temper = -1.0
!    end if
!    read *, str, timeq, timrun, dt   ! number of cycles
!    read *, str, timprn              ! print interval
!
!    read *, str, rcut          ! ewald params
!    read *, str, kmax
!    read *, str, skin                ! width of the skin of the nnlist
!    read *, str, rdel
!    read *, str, rfrst, rscnd
!
!    read *, str, filname              ! output configuration

    ! name files
    it = len_trim(filname)
    filfld(1:it) = filname
    filtop(1:it) = filname
    filcfi(1:it) = filname
    filcfo(1:it) = filname
    filmeas(1:it) = filname
    filfld(it+1:it+4) = '.fld'
    filtop(it+1:it+4) = '.top'
    filcfi(it+1:it+4) = '.cfg'
    filcfo(it+1:it+6) = '_n.cfg'
    filmeas(it+1:it+4) = '.crl'


    ! communicate input variables
    call bcast_dbl(temper)
    call bcast_dbl(timeq)
    call bcast_dbl(timrun)
    call bcast_dbl(dt)
    call bcast_dbl(timprn)
    call bcast_dbl(rcut)
    call bcast_int(kmax)
    call bcast_dbl(skin)
    call bcast_dbl(rdel)

    call bcast_log(l2d)
    call bcast_log(lhb)
    if (l2d .or. lhb) then
        call bcast_dbl(rfrst)
        call bcast_dbl(rscnd)
    end if
    if (lhb) then
        call bcast_dbl(rishell)
    end if

    call bcast_log(ldiff)
    if (ldiff) then
        call bcast_int(rbin)
        call bcast_dbl(timdif)
        call bcast_dbl(diftimax)
    endif

    call bcast_log(lcouett)
    if (lcouett) then
        call bcast_dbl(dshx)
        call bcast_dbl(dshy)
        call bcast_dbl(dshz)
    end if
    call bcast_log(lfext)
    if (lfext) then
        call bcast_dbl(fxext)
        call bcast_dbl(fyext)
        call bcast_dbl(fzext)
    end if

    call bcast_log(lchempot)
    if (lchempot) then
        call bcast_int(ichmax)
        call bcast_dbl(lwmin)
        call bcast_dbl(lwmax)
        call bcast_int(lchmax)
        call bcast_int(ilmid)
    end if

    ! force field
    call sysdef_zero
    call sysdef_input
    call sysdef_alloc

    ! configuration
    call config_zero
    if (me == 0) call config_input

    call bcast_dbl(lx)
    call bcast_dbl(ly)
    call bcast_dbl(lz)
    call bcast_int1d(itype)
    call bcast_dbl1d(x0) ; call bcast_dbl1d(y0) ; call bcast_dbl1d(z0) 
    call bcast_dbl1d(x1) ; call bcast_dbl1d(y1) ; call bcast_dbl1d(z1) 
    call bcast_dbl1d(x2) ; call bcast_dbl1d(y2) ; call bcast_dbl1d(z2) 
    call bcast_dbl1d(x3) ; call bcast_dbl1d(y3) ; call bcast_dbl1d(z3) 
    call bcast_dbl1d(x4) ; call bcast_dbl1d(y4) ; call bcast_dbl1d(z4) 

    call bcast_dbl1d(cx0) ; call bcast_dbl1d(cy0) ; call bcast_dbl1d(cz0) 
    call bcast_dbl1d(cx1) ; call bcast_dbl1d(cy1) ; call bcast_dbl1d(cz1) 
    call bcast_dbl1d(cx2) ; call bcast_dbl1d(cy2) ; call bcast_dbl1d(cz2) 
    call bcast_dbl1d(cx3) ; call bcast_dbl1d(cy3) ; call bcast_dbl1d(cz3) 
    call bcast_dbl1d(cx4) ; call bcast_dbl1d(cy4) ; call bcast_dbl1d(cz4) 

    call bcast_dbl1d(wx0) ; call bcast_dbl1d(wy0) ; call bcast_dbl1d(wz0) 
    call bcast_dbl1d(wx1) ; call bcast_dbl1d(wy1) ; call bcast_dbl1d(wz1) 
    call bcast_dbl1d(wx2) ; call bcast_dbl1d(wy2) ; call bcast_dbl1d(wz2) 
    call bcast_dbl1d(wx3) ; call bcast_dbl1d(wy3) ; call bcast_dbl1d(wz3) 
    call bcast_dbl1d(wx4) ; call bcast_dbl1d(wy4) ; call bcast_dbl1d(wz4) 

    call bcast_dbl1d(q10) ; call bcast_dbl1d(q20) ; call bcast_dbl1d(q30)  ; call bcast_dbl1d(q40) 
    call bcast_dbl1d(q11) ; call bcast_dbl1d(q21) ; call bcast_dbl1d(q31)  ; call bcast_dbl1d(q41) 
    call bcast_dbl1d(q12) ; call bcast_dbl1d(q22) ; call bcast_dbl1d(q32)  ; call bcast_dbl1d(q42) 
    call bcast_dbl1d(q13) ; call bcast_dbl1d(q23) ; call bcast_dbl1d(q33)  ; call bcast_dbl1d(q43) 
    call bcast_dbl1d(q14) ; call bcast_dbl1d(q24) ; call bcast_dbl1d(q34)  ; call bcast_dbl1d(q44) 

    call bcast_dbl1d(a11) ; call bcast_dbl1d(a12) ; call bcast_dbl1d(a13) 
    call bcast_dbl1d(a21) ; call bcast_dbl1d(a22) ; call bcast_dbl1d(a23) 
    call bcast_dbl1d(a31) ; call bcast_dbl1d(a32) ; call bcast_dbl1d(a33) 

    call bcast_dbl1d(rsx) ; call bcast_dbl1d(rsy) ; call bcast_dbl1d(rsz) 
    call bcast_dbl1d(xold) ; call bcast_dbl1d(zold) ; call bcast_dbl1d(zold) 

    call bcast_dbl(s0)
    call bcast_dbl(s1)
    call bcast_dbl(s2)
    call bcast_dbl(s3)
    call bcast_dbl(s4)
    call bcast_dbl(rs0)
    call bcast_dbl(rs1)
    call bcast_dbl(rs2)
    call bcast_dbl(rs3)
    call bcast_dbl(rs4)


    return
end subroutine input

!
!====== Initialization of global variables + precomputation ====================
!
subroutine init
    use move
    implicit none
 
    call sysdef_init
    call config_init
    call force_init
    call move_init
    call measure_init

    return
end subroutine init

!
!======================= Prints runtime data ==================================
!
subroutine prntout(it)
    use move
    implicit none
    integer*4 :: it

    if (me == 0) then
        write(*, 10) dt*real(it), tempnow, tkps, press, uavex/1000, htot/1000
        write(*, *) 'fz', ftotz(1)/real(it), ftotz(-1)/real(it)
        write(*, *) 'uavex', uavextot/real(it), uavex, uconf
        !write(*, *) 'chp', -log(pavesum/real(ichmax*it*nproc))*PH_R*temper/1000.0, nproc
        !write(*, *) 'u,ch', uavextot/real(it), -log(sum(pavesum(:))/real(sum(inavesum(:))))*PH_R*chtemper/1000.0
    end if

    return
10  format(f9.1,5(1x,f10.3))
end subroutine prntout

!
!======================= Writes output data ==================================
!
subroutine output
    use move
    implicit none

    if (me == 0) call config_output
    call config_dealloc
    call measure_output

    return
end subroutine output
