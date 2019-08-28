module config  ! state of the system : configuration, energy, forces, ...
    use sysdef
    implicit none

    ! files
    character*40 :: filcfi, filcfo

    ! time
    real*8 :: timeq, timrun                        ! time & time step
    integer*8 :: itime, iteq, itrun                       ! counting timesteps

    ! volume
    real*8 :: lx, ly, lz, lxh, lyh, lzh, lxi, lyi, lzi, lxhi, lyhi, lzhi

    ! temperature and thermostat
    real*8 :: temper, htot, tempnow, qmass, rmass
    real*8 :: s0, s1, s2, s3, s4, rs0, rs1, rs2, rs3, rs4, ke, sft, sfr

    ! all atoms
    real*8, dimension(:), allocatable :: x0, y0, z0          ! positions of all atoms

    ! atmol movement
    real*8, dimension(:), allocatable :: x1, y1, z1          ! velocities of com
    real*8, dimension(:), allocatable :: x2, y2, z2          ! accelerations of com
    real*8, dimension(:), allocatable :: x3, y3, z3          ! da/dt of com
    real*8, dimension(:), allocatable :: x4, y4, z4          ! da/dt^2 of com

    ! rigmol movement
    real*8, dimension(:), allocatable :: cx0, cy0, cz0          ! positions of com
    real*8, dimension(:), allocatable :: cx1, cy1, cz1          ! velocities of com
    real*8, dimension(:), allocatable :: cx2, cy2, cz2          ! accelerations of com
    real*8, dimension(:), allocatable :: cx3, cy3, cz3          ! da/dt of com
    real*8, dimension(:), allocatable :: cx4, cy4, cz4          ! da/dt^2 of com

    real*8, dimension(:), allocatable :: q10, q20, q30, q40      ! orientation
    real*8, dimension(:), allocatable :: q11, q21, q31, q41      ! 1.st deriv
    real*8, dimension(:), allocatable :: q12, q22, q32, q42      ! 2.st deriv
    real*8, dimension(:), allocatable :: q13, q23, q33, q43      ! 3.st deriv
    real*8, dimension(:), allocatable :: q14, q24, q34, q44      ! 4.st deriv

    real*8, dimension(:), allocatable :: wx0, wy0, wz0          ! angular momentum
    real*8, dimension(:), allocatable :: wx1, wy1, wz1          ! dw/dt
    real*8, dimension(:), allocatable :: wx2, wy2, wz2          ! dw/dt^2
    real*8, dimension(:), allocatable :: wx3, wy3, wz3          ! dw/dt^3
    real*8, dimension(:), allocatable :: wx4, wy4, wz4          ! dw/dt^4

    real*8, dimension(:), allocatable :: xold, yold, zold

    ! auxiliary site positions
    real*8, dimension(:), allocatable :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    real*8, dimension(:), allocatable :: rsx, rsy, rsz  ! site position

    contains

    subroutine config_init
        implicit none

        ! time
        iteq = nint(timeq/dt)
        itrun = nint(timrun/dt)

        ! volume
        lxh = 0.5*lx
        lyh = 0.5*ly
        lzh = 0.5*lz
 
        lxi = 1.0/lx
        lyi = 1.0/ly
        lzi = 1.0/lz

        lxhi = 1.0/lxh
        lyhi = 1.0/lyh
        lzhi = 1.0/lzh

    ! from 1fs to required time units
        x1 = x1*dt    ; y1 = y1*dt    ; z1 = z1*dt    
        x2 = x2*dt**2 ; y2 = y2*dt**2 ; z2 = z2*dt**2 
        x3 = x3*dt**3 ; y3 = y3*dt**3 ; z3 = z3*dt**3 
        x4 = x4*dt**4 ; y4 = y4*dt**4 ; z4 = z4*dt**4 

        cx1 = cx1*dt    ; cy1 = cy1*dt    ; cz1 = cz1*dt    
        cx2 = cx2*dt**2 ; cy2 = cy2*dt**2 ; cz2 = cz2*dt**2 
        cx3 = cx3*dt**3 ; cy3 = cy3*dt**3 ; cz3 = cz3*dt**3 
        cx4 = cx4*dt**4 ; cy4 = cy4*dt**4 ; cz4 = cz4*dt**4 

        wx0 = wx0*dt    ; wy0 = wy0*dt    ; wz0 = wz0*dt    
        wx1 = wx1*dt**2 ; wy1 = wy1*dt**2 ; wz1 = wz1*dt**2 
        wx2 = wx2*dt**3 ; wy2 = wy2*dt**3 ; wz2 = wz2*dt**3 
        wx3 = wx3*dt**4 ; wy3 = wy3*dt**4 ; wz3 = wz3*dt**4 
        wx4 = wx4*dt**5 ; wy4 = wy4*dt**5 ; wz4 = wz4*dt**5 

        call config_quatsite(1_8, nrig)

        return
    end subroutine config_init

    subroutine config_zero
        implicit none

        ! allocation 

        allocate(x0(ntota), y0(ntota), z0(ntota))

        allocate(x1(1:ntota), y1(1:ntota), z1(1:ntota))
        allocate(x2(1:ntota), y2(1:ntota), z2(1:ntota))
        allocate(x3(1:ntota), y3(1:ntota), z3(1:ntota))
        allocate(x4(1:ntota), y4(1:ntota), z4(1:ntota))

        allocate(cx0(nrig), cy0(nrig), cz0(nrig))
        allocate(cx1(nrig), cy1(nrig), cz1(nrig))
        allocate(cx2(nrig), cy2(nrig), cz2(nrig))
        allocate(cx3(nrig), cy3(nrig), cz3(nrig))
        allocate(cx4(nrig), cy4(nrig), cz4(nrig))

        allocate(q10(nrig), q20(nrig), q30(nrig), q40(nrig))
        allocate(q11(nrig), q21(nrig), q31(nrig), q41(nrig))
        allocate(q12(nrig), q22(nrig), q32(nrig), q42(nrig))
        allocate(q13(nrig), q23(nrig), q33(nrig), q43(nrig))
        allocate(q14(nrig), q24(nrig), q34(nrig), q44(nrig))

        allocate(wx0(nrig), wy0(nrig), wz0(nrig))
        allocate(wx1(nrig), wy1(nrig), wz1(nrig))
        allocate(wx2(nrig), wy2(nrig), wz2(nrig))
        allocate(wx3(nrig), wy3(nrig), wz3(nrig))
        allocate(wx4(nrig), wy4(nrig), wz4(nrig))

        allocate(a11(nrig), a12(nrig), a13(nrig))
        allocate(a21(nrig), a22(nrig), a23(nrig))
        allocate(a31(nrig), a32(nrig), a33(nrig))

        allocate(rsx(nmove), rsy(nmove), rsz(nmove))

        allocate(xold(nliqa+1:ntota), yold(nliqa+1:ntota), zold(nliqa+1:ntota))

        ! zeroing

        x0 = 0 ; y0 = 0 ; z0 = 0          
        x1 = 0 ; y1 = 0 ; z1 = 0          
        x2 = 0 ; y2 = 0 ; z2 = 0          
        x3 = 0 ; y3 = 0 ; z3 = 0
        x4 = 0 ; y4 = 0 ; z4 = 0          

        cx0 = 0 ; cy0 = 0 ; cz0 = 0          
        cx1 = 0 ; cy1 = 0 ; cz1 = 0          
        cx2 = 0 ; cy2 = 0 ; cz2 = 0          
        cx3 = 0 ; cy3 = 0 ; cz3 = 0
        cx4 = 0 ; cy4 = 0 ; cz4 = 0          
           
        q10 = 0 ; q20 = 0 ; q30 = 0 ; q40 = 0   
        q11 = 0 ; q21 = 0 ; q31 = 0 ; q41 = 0   
        q12 = 0 ; q22 = 0 ; q32 = 0 ; q42 = 0   
        q13 = 0 ; q23 = 0 ; q33 = 0 ; q43 = 0   
        q14 = 0 ; q24 = 0 ; q34 = 0 ; q44 = 0   
                           
        wx0 = 0 ; wy0 = 0 ; wz0 = 0        
        wx1 = 0 ; wy1 = 0 ; wz1 = 0        
        wx2 = 0 ; wy2 = 0 ; wz2 = 0        
        wx3 = 0 ; wy3 = 0 ; wz3 = 0        
        wx4 = 0 ; wy4 = 0 ; wz4 = 0        

        xold = 0 ; yold = 0 ; zold = 0

        s0 = 0 ; s1 = 0 ; s2 = 0 ; s3 = 0 ; s4 = 0
        rs0 = 0 ; rs1 = 0 ; rs2 = 0 ; rs3 = 0 ; rs4 = 0

        rsx = 0.0 ; rsy = 0.0; rsz = 0.0

        return
    end subroutine config_zero

    !
    ! ====== READS CONFIGURATION FROM A TEXT FILE ======
    !
    subroutine config_input
        implicit none
        character*80 :: aux
        character*40 :: filcfg
        character*20 :: mname, mtype
        integer*8 :: iu, it, ntx, inow

        iu = 4
        open(iu, file = filcfi, status = 'unknown')

            read(iu, *) aux
            read(iu, *) lx, ly, lz
            read(iu, *)

            read(iu, *) aux
            read(iu, *) ntx
            if (ntx /= ntype) stop 'config_input(): wrong number of types!'
            read(iu, *)

            read(iu, *) aux
            nt(0) = 0
            inow = 0
!                end if
            do it = 1, ntype
                read(iu, *) mname, mtype, ntx
                if ((mname .ne. molname(it)) .or. (mtype .ne. moltype(it)) .or. (ntx /= nt(it))) then
                    print *, mname, molname(it)
                    print *, mtype, moltype(it)
                    print *, ntx, nt(it)
                    stop 'config_input(): Error in configfile!'
                end if
                select case (moltype(it)) 
                    case ('rigmol')
                        call read_rigmol(iu, inow + 1, inow + nt(it), ns(it))
                    case ('atmol')
                        call read_atmol(iu, inow + 1, inow + nt(it)*ns(it))
                    case ('surfmol')
                        call read_surfmol(iu, inow + 1, inow + nt(it)*ns(it))
                    case ('crystal')
                        call read_crystal(iu, inow + 1, inow + nt(it)*ns(it))
                    case default
                        stop 'cfg_out(): unknown moltype'
                end select
                inow = inow + nt(it)
                read(iu, *)
            end do

            read(iu, *) aux    ! thermostat
            call read_thermostat(iu)

        close(iu, status = 'keep')

!        x1 = x1*s0
!        y1 = y1*s0
!        z1 = z1*s0

!        cx1 = cx1*s0
!        cy1 = cy1*s0
!        cz1 = cz1*s0

!        wx0 = wx0*rs0
!        wy0 = wy0*rs0
!        wz0 = wz0*rs0

!        s0 = 1.0 ; s1 = 0.0 ; s2 = 0.0 ; s3 = 0.0 ; s4 = 0.0
!        rs0 = 1.0 ; rs1 = 0.0 ; rs2 = 0.0 ; rs3 = 0.0 ; rs4 = 0.0

        return
    end subroutine config_input 

    subroutine config_output
        implicit none
        character*20 :: filcfg
        integer*8 :: i, iu, it, inow

    ! to 1fs units
        x1 = x1/dt    ; y1 = y1/dt    ; z1 = z1/dt    
        x2 = x2/dt**2 ; y2 = y2/dt**2 ; z2 = z2/dt**2 
        x3 = x3/dt**3 ; y3 = y3/dt**3 ; z3 = z3/dt**3 
        x4 = x4/dt**4 ; y4 = y4/dt**4 ; z4 = z4/dt**4 

        cx1 = cx1/dt    ; cy1 = cy1/dt    ; cz1 = cz1/dt    
        cx2 = cx2/dt**2 ; cy2 = cy2/dt**2 ; cz2 = cz2/dt**2 
        cx3 = cx3/dt**3 ; cy3 = cy3/dt**3 ; cz3 = cz3/dt**3 
        cx4 = cx4/dt**4 ; cy4 = cy4/dt**4 ; cz4 = cz4/dt**4 

        wx0 = wx0/dt    ; wy0 = wy0/dt    ; wz0 = wz0/dt    
        wx1 = wx1/dt**2 ; wy1 = wy1/dt**2 ; wz1 = wz1/dt**2 
        wx2 = wx2/dt**3 ; wy2 = wy2/dt**3 ; wz2 = wz2/dt**3 
        wx3 = wx3/dt**4 ; wy3 = wy3/dt**4 ; wz3 = wz3/dt**4 
        wx4 = wx4/dt**5 ; wy4 = wy4/dt**5 ; wz4 = wz4/dt**5 

        iu = 4
        open(iu, file = filcfo, status = 'unknown')
            write(iu, *) '# SIMULATION BOX SIZE'
            write(iu, 10) lx, ly, lz
            write(iu, *)
            write(iu, *) '# NUMBER OF PARTICLE TYPES'
            write(iu, *) ntype
            write(iu, *)
            write(iu, *) '# PARTICLE COORDINATES'
            inow = 0
            do it = 1, ntype           ! molecular types
                write(iu, *) molname(it), moltype(it), nt(it) ! nmt - length of the list
                select case (moltype(it)) 
                    case ('rigmol')
                        call write_rigmol(iu, inow + 1, inow + nt(it), ns(it))
                    case ('atmol')
                        call write_atmol(iu, inow + 1, inow + nt(it)*ns(it))
                    case ('surfmol')
                        call write_surfmol(iu, inow + 1, inow + nt(it)*ns(it))
                    case ('crystal')
                        call write_crystal(iu, inow + 1, inow + nt(it)*ns(it))
                    case default
                        stop 'cfg_out(): unknown moltype'
                end select
                inow = inow + nt(it)*ns(it)
                write(iu, *)
            end do

            write(iu, *) '# THERMOSTAT'
            call write_thermostat(iu)

        close(iu, status = 'keep')

        return
10      format(3(1x,f14.10))
    end subroutine config_output  

    subroutine config_quatsite(imin, imax)
        implicit none
        integer*8 :: ix, i, j, imin, imax
 
        forall (i = imin:imax)
            a11(i) = -q10(i)*q10(i) + q20(i)*q20(i) - q30(i)*q30(i) + q40(i)*q40(i)
            a12(i) =  2.0*(q30(i)*q40(i) - q10(i)*q20(i))
            a13(i) =  2.0*(q20(i)*q30(i) + q10(i)*q40(i))
            a21(i) = -2.0*(q10(i)*q20(i) + q30(i)*q40(i))
            a22(i) =  q10(i)*q10(i) - q20(i)*q20(i) - q30(i)*q30(i) + q40(i)*q40(i)
            a23(i) =  2.0*(q20(i)*q40(i) - q10(i)*q30(i))
            a31(i) =  2.0*(q20(i)*q30(i) - q10(i)*q40(i))
            a32(i) = -2.0*(q10(i)*q30(i) + q40(i)*q20(i))
            a33(i) = -q10(i)*q10(i) - q20(i)*q20(i) + q30(i)*q30(i) + q40(i)*q40(i)
        end forall
 
        do i = imin, imax
            ix = 3*(i-1)            ! just for molecules with 3 sites (spce water)
            forall (j = 1:3)
                x0(ix+j) = cx0(i) + a11(i)*xbody(j) + a21(i)*ybody(j) + a31(i)*zbody(j)
                y0(ix+j) = cy0(i) + a12(i)*xbody(j) + a22(i)*ybody(j) + a32(i)*zbody(j)
                z0(ix+j) = cz0(i) + a13(i)*xbody(j) + a23(i)*ybody(j) + a33(i)*zbody(j)
            end forall
        end do
  
        return
    end subroutine config_quatsite

    subroutine config_wdump(fil) 
        implicit none
        character*40 :: fil

        call config_quatsite(1_8,nrig)

        open(1, file=fil, form='unformatted')
            write(1) x0(1:nmove), y0(1:nmove), z0(1:nmove)
        close(1, status='keep')

        return
    end subroutine config_wdump 

    subroutine config_wdedump(fil) 
        implicit none
        character*40 :: fil

        open(1, file=fil, form='unformatted', status='old')
            read(1) x0(1:nmove), y0(1:nmove), z0(1:nmove)
        close(1, status='keep')

        return
    end subroutine config_wdedump 

    subroutine config_dump(fil) 
        implicit none
        character*40 :: fil

        open(1, file=fil, form='unformatted')
            write(1) cx0,cy0,cz0,cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3,cx4,cy4,cz4
            write(1) q10,q20,q30,q40,q11,q21,q31,q41,q12,q22,q32,q42,q13,q23,q33,q43,q14,q24,q34,q44
            write(1) wx0,wy0,wz0,wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3,wx4,wy4,wz4
            write(1) x0(nriga+1:nmove), y0(nriga+1:nmove), z0(nriga+1:nmove)
            write(1) x1(nriga+1:nliqa), y1(nriga+1:nliqa), z1(nriga+1:nliqa)
            write(1) x2(nriga+1:nliqa), y2(nriga+1:nliqa), z2(nriga+1:nliqa)
            write(1) x3(nriga+1:nliqa), y3(nriga+1:nliqa), z3(nriga+1:nliqa)
            write(1) x4(nriga+1:nliqa), y4(nriga+1:nliqa), z4(nriga+1:nliqa)
            write(1) x1(nliqa+1:nmove), y1(nliqa+1:nmove), z1(nliqa+1:nmove)
            write(1) xold(nliqa+1:nmove), yold(nliqa+1:nmove), zold(nliqa+1:nmove)
        close(1, status='keep')

        return
    end subroutine config_dump 

    subroutine config_dedump(fil) 
        implicit none
        character*40 :: fil

        open(1, file=fil, form='unformatted', status='old')
            read(1) cx0,cy0,cz0,cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3,cx4,cy4,cz4
            read(1) q10,q20,q30,q40,q11,q21,q31,q41,q12,q22,q32,q42,q13,q23,q33,q43,q14,q24,q34,q44
            read(1) wx0,wy0,wz0,wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3,wx4,wy4,wz4
            read(1) x0(nriga+1:nmove), y0(nriga+1:nmove), z0(nriga+1:nmove)
            read(1) x1(nriga+1:nliqa), y1(nriga+1:nliqa), z1(nriga+1:nliqa)
            read(1) x2(nriga+1:nliqa), y2(nriga+1:nliqa), z2(nriga+1:nliqa)
            read(1) x3(nriga+1:nliqa), y3(nriga+1:nliqa), z3(nriga+1:nliqa)
            read(1) x4(nriga+1:nliqa), y4(nriga+1:nliqa), z4(nriga+1:nliqa)
            read(1) x1(nliqa+1:nmove), y1(nliqa+1:nmove), z1(nliqa+1:nmove)
            write(1) xold(nliqa+1:nmove), yold(nliqa+1:nmove), zold(nliqa+1:nmove)
        close(1, status='keep')

        return
    end subroutine config_dedump 

    subroutine read_rigmol(un, istart, iend, nsit)
        implicit none
        integer*8 :: un, istart, iend, i, is, nsit, j

        do j = istart, iend
            read(un, *) i, (itype(nsit*(i - 1) + is), is = 1, nsit)

            read(un, *) cx0(i), cy0(i), cz0(i)         ! positions and derivatives
            read(un, *) cx1(i), cy1(i), cz1(i)
            read(un, *) cx2(i), cy2(i), cz2(i)
            read(un, *) cx3(i), cy3(i), cz3(i)
            read(un, *) cx4(i), cy4(i), cz4(i)
           
            read(un, *) q10(i), q20(i), q30(i), q40(i)      ! quaternions and deriv.
            read(un, *) q11(i), q21(i), q31(i), q41(i)
            read(un, *) q12(i), q22(i), q32(i), q42(i)
            read(un, *) q13(i), q23(i), q33(i), q43(i)
            read(un, *) q14(i), q24(i), q34(i), q44(i)
           
            read(un, *) wx0(i), wy0(i), wz0(i)      ! angular momenta and derivs.
            read(un, *) wx1(i), wy1(i), wz1(i)
            read(un, *) wx2(i), wy2(i), wz2(i)
            read(un, *) wx3(i), wy3(i), wz3(i)
            read(un, *) wx4(i), wy4(i), wz4(i)
        enddo

        return
    end subroutine read_rigmol


    subroutine read_atmol(un, istart, iend)
        implicit none
        integer*8 :: un, istart, iend, i, j

        do j = istart, iend
            read(un, *) i, itype(i)
            read(un, *) x0(i), y0(i), z0(i)         ! r = positions (suitable for wall)
            read(un, *) x1(i), y1(i), z1(i)         ! r' = velocities
            read(un, *) x2(i), y2(i), z2(i)         ! r'' = acceleration
            read(un, *) x3(i), y3(i), z3(i)         ! r'''
            read(un, *) x4(i), y4(i), z4(i)         ! r''''
        enddo

        return
    end subroutine read_atmol

    subroutine read_surfmol(un, istart, iend)
        implicit none
        integer*8 :: un, istart, iend, i, j

        do j = istart, iend
            read(un, *) i, itype(i)
            read(un, *) x0(i), y0(i), z0(i)         ! r = positions (suitable for wall)
            read(un, *) xold(i), yold(i), zold(i)         ! r = positions (suitable for wall)
            xold(i) = x0(i) + (xold(i) - x0(i))*dt
            yold(i) = y0(i) + (yold(i) - y0(i))*dt
            zold(i) = z0(i) + (zold(i) - z0(i))*dt
!            x1(i) = x0(i)
!            y1(i) = y0(i)
!            z1(i) = z0(i)
        enddo

        return
    end subroutine read_surfmol

    subroutine read_crystal(un, istart, iend)
        implicit none
        integer*8 :: un, istart, iend, i, j

        do j = istart, iend
            read(un, *) i, itype(i)                  ! number and charge
            read(un, *) x0(i), y0(i), z0(i)         ! r = positions (suitable for wall)
            walltype(i) = itype(i)
        enddo

        return
    end subroutine read_crystal

    subroutine read_thermostat(iu) 
        implicit none
        integer*8 :: iu, der

        read(iu, *) s0, s1, s2, s3, s4
        read(iu, *) rs0, rs1, rs2, rs3, rs4

        return
    end subroutine read_thermostat


    subroutine write_rigmol(un, istart, iend, nsit)
        implicit none
        integer*8 :: un, istart, iend, i, nsit, is

        do i = istart, iend
            write(un, 30) i, (itype(nsit*(i - 1) + is), is = 1, nsit)
           
            write(un, 10) cx0(i), cy0(i), cz0(i)         ! positions and derivatives
            write(un, 10) cx1(i), cy1(i), cz1(i)
            write(un, 10) cx2(i), cy2(i), cz2(i)
            write(un, 10) cx3(i), cy3(i), cz3(i)
            write(un, 10) cx4(i), cy4(i), cz4(i)
           
            write(un, 20) q10(i), q20(i), q30(i), q40(i)      ! quaternions and deriv.
            write(un, 20) q11(i), q21(i), q31(i), q41(i)
            write(un, 20) q12(i), q22(i), q32(i), q42(i)
            write(un, 20) q13(i), q23(i), q33(i), q43(i)
            write(un, 20) q14(i), q24(i), q34(i), q44(i)
           
            write(un, 10) wx0(i), wy0(i), wz0(i)      ! angular momenta and derivs.
            write(un, 10) wx1(i), wy1(i), wz1(i)
            write(un, 10) wx2(i), wy2(i), wz2(i)
            write(un, 10) wx3(i), wy3(i), wz3(i)
            write(un, 10) wx4(i), wy4(i), wz4(i)
        enddo

        return
10     format(3(1x,e15.8))
20     format(4(1x,e15.8))
30     format(4(1x,i5))
    end subroutine write_rigmol

    subroutine write_atmol(un, istart, iend)
        implicit none
        integer*8 :: un, istart, iend, i, it

        it = itype(istart) - 1
        do i = istart, iend
            write(un, 30) i, itype(i) - it             ! numbe and type
            write(un, 10) x0(i), y0(i), z0(i)         ! r = positions (suitable for wall)
            write(un, 10) x1(i), y1(i), z1(i)         ! r' = velocities
            write(un, 10) x2(i), y2(i), z2(i)         ! r'' = acceleration
            write(un, 10) x3(i), y3(i), z3(i)         ! r'''
            write(un, 10) x4(i), y4(i), z4(i)         ! r''''
        end do

        return
10     format(3(1x,e15.8))
30     format(2(1x,i5))
    end subroutine write_atmol

    subroutine write_surfmol(un, istart, iend)
        implicit none
        integer*8 :: un, istart, iend, i, it

        it = itype(istart) - 1
        do i = istart, iend
            write(un, 30) i, itype(i) - it             ! numbe and type
            write(un, 10) x0(i), y0(i), z0(i)         ! r = positions (suitable for wall)
            xold(i) = x0(i) + (xold(i) - x0(i))/dt
            yold(i) = y0(i) + (yold(i) - y0(i))/dt
            zold(i) = z0(i) + (zold(i) - z0(i))/dt
            write(un, 10) xold(i), yold(i), zold(i)         ! r = positions (suitable for wall)
        end do

        return
10     format(3(1x,e15.8))
30     format(2(1x,i5))
    end subroutine write_surfmol

    subroutine write_crystal(un, istart, iend)
        implicit none
        integer*8 :: un, istart, iend, i, it

        it = itype(istart) - 1
        do i = istart, iend
!            walltype(i) = itype(i)
!            write(un, *) i, itype(i) - it  ! number, type, and charge
            write(un, 30) i, walltype(i)  ! number, type, and charge
            write(un, 10) x0(i), y0(i), z0(i)              ! r = positions (suitable for wall)
        enddo

30     format(2(1x,i5))
10     format(3(1x,e15.8))
        return
    end subroutine write_crystal

    subroutine write_thermostat(iu) 
        implicit none
        integer*8 :: iu, der

        write(iu, 10) s0, s1, s2, s3, s4
        write(iu, 10) rs0, rs1, rs2, rs3, rs4

        return
10      format(5(1x,e15.8))
    end subroutine write_thermostat

    subroutine config_dealloc
        implicit none

        deallocate(x0, y0, z0)
        deallocate(x1, y1, z1)
        deallocate(x2, y2, z2)
        deallocate(x3, y3, z3)
        deallocate(x4, y4, z4)
        deallocate(cx0, cy0, cz0)
        deallocate(cx1, cy1, cz1)
        deallocate(cx2, cy2, cz2)
        deallocate(cx3, cy3, cz3)
        deallocate(cx4, cy4, cz4)
        deallocate(q10, q20, q30, q40)
        deallocate(q11, q21, q31, q41)
        deallocate(q12, q22, q32, q42)
        deallocate(q13, q23, q33, q43)
        deallocate(q14, q24, q34, q44)
        deallocate(wx0, wy0, wz0)
        deallocate(wx1, wy1, wz1)
        deallocate(wx2, wy2, wz2)
        deallocate(wx3, wy3, wz3)
        deallocate(wx4, wy4, wz4)
        deallocate(a11, a12, a13)
        deallocate(a21, a22, a23)
        deallocate(a31, a32, a33)
        deallocate(rsx, rsy, rsz)
        deallocate(xold)
        deallocate(yold)
        deallocate(zold)

        return
    end subroutine config_dealloc

end module config

