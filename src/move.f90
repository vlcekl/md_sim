module move
    use force

    implicit none
    real*8, dimension(:), allocatable :: cort, cor1, cor2, cor3, cor4
    real*8 :: kp, kw, ket, ker, inxxi, inyyi, inzzi
    integer*4 :: tmin, tmax, wmin, wmax, ifi, iti, nfree

    ! shake variables
    real*8 :: tol2, rptol, tkps, kps(TMX)
    integer*4, dimension(:), allocatable :: moving, moved
    integer*4 :: ismin, ismax, maxit, nfreed

    contains

    subroutine move_init
        implicit none

        inxxi = 1.0/inxx
        inyyi = 1.0/inyy
        inzzi = 1.0/inzz

        ! can be easily adjusted for MPI (sharing within a node)
        tmin = nriga + 1
        tmax = nliqa
!        if ((tmin < TNMI) .or. (tmax > MNMX)) then
!            write(*, *) tmin, tmax, TNMI, MNMX
!            stop 'move_init(): Wrong range for moving atoms!'
!        end if
        wmin = 1
        wmax = nrig       ! suppose nrig == nwat
        nfree = nliqa - nriga

        ! transformation from atomic fx array to rigmol form (just for water)
        ifi = (3-1)*nrig  ! fx(at) -> fx(mol)
        iti = (3-2)*nrig  ! fx(at) -> tx(mol)


        allocate(cort(nriga+1:nmove))
        allocate(cor1(nrig), cor2(nrig), cor3(nrig), cor4(nrig))


        call move_shake_init

        return
    end subroutine move_init

    subroutine move_predict
        implicit none
        integer*4 ::  i, j, ix
        real*8 :: tstar


        ! translation
        forall (i = tmin:tmax)
            x0(i) = x0(i) +   x1(i) +   x2(i) +   x3(i) + x4(i)
            x1(i) = x1(i) + 2*x2(i) + 3*x3(i) + 4*x4(i)
            x2(i) = x2(i) + 3*x3(i) + 6*x4(i)
            x3(i) = x3(i) + 4*x4(i)
            y0(i) = y0(i) +   y1(i) +   y2(i) +   y3(i) + y4(i)
            y1(i) = y1(i) + 2*y2(i) + 3*y3(i) + 4*y4(i)
            y2(i) = y2(i) + 3*y3(i) + 6*y4(i)
            y3(i) = y3(i) + 4*y4(i)
            z0(i) = z0(i) +   z1(i) +   z2(i) +   z3(i) + z4(i)
            z1(i) = z1(i) + 2*z2(i) + 3*z3(i) + 4*z4(i)
            z2(i) = z2(i) + 3*z3(i) + 6*z4(i)
            z3(i) = z3(i) + 4*z4(i)
        end forall
        kp = sum(mi(tmin:tmax)*(x1(tmin:tmax)**2 + y1(tmin:tmax)**2 + z1(tmin:tmax)**2))

        ! rigmol translation and rotation (just for water)
        forall (i = wmin:wmax)
            cx0(i) = cx0(i) +   cx1(i) +   cx2(i) +   cx3(i) + cx4(i)
            cx1(i) = cx1(i) + 2*cx2(i) + 3*cx3(i) + 4*cx4(i)
            cx2(i) = cx2(i) + 3*cx3(i) + 6*cx4(i)
            cx3(i) = cx3(i) + 4*cx4(i)
            cy0(i) = cy0(i) +   cy1(i) +   cy2(i) +   cy3(i) + cy4(i)
            cy1(i) = cy1(i) + 2*cy2(i) + 3*cy3(i) + 4*cy4(i)
            cy2(i) = cy2(i) + 3*cy3(i) + 6*cy4(i)
            cy3(i) = cy3(i) + 4*cy4(i)
            cz0(i) = cz0(i) +   cz1(i) +   cz2(i) +   cz3(i) + cz4(i)
            cz1(i) = cz1(i) + 2*cz2(i) + 3*cz3(i) + 4*cz4(i)
            cz2(i) = cz2(i) + 3*cz3(i) + 6*cz4(i)
            cz3(i) = cz3(i) + 4*cz4(i)
        end forall
        kp = kp + mw*(sum(cx1(wmin:wmax)**2) + sum(cy1(wmin:wmax)**2) + sum(cz1(wmin:wmax)**2))

        forall (i = wmin:wmax)
            wx0(i) = wx0(i) +   wx1(i) +   wx2(i) +   wx3(i) + wx4(i)
            wx1(i) = wx1(i) + 2*wx2(i) + 3*wx3(i) + 4*wx4(i)
            wx2(i) = wx2(i) + 3*wx3(i) + 6*wx4(i)
            wx3(i) = wx3(i) + 4*wx4(i)
            wy0(i) = wy0(i) +   wy1(i) +   wy2(i) +   wy3(i) + wy4(i)
            wy1(i) = wy1(i) + 2*wy2(i) + 3*wy3(i) + 4*wy4(i)
            wy2(i) = wy2(i) + 3*wy3(i) + 6*wy4(i)
            wy3(i) = wy3(i) + 4*wy4(i)
            wz0(i) = wz0(i) +   wz1(i) +   wz2(i) +   wz3(i) + wz4(i)
            wz1(i) = wz1(i) + 2*wz2(i) + 3*wz3(i) + 4*wz4(i)
            wz2(i) = wz2(i) + 3*wz3(i) + 6*wz4(i)
            wz3(i) = wz3(i) + 4*wz4(i)
        end forall
        kw = inxx*sum(wx0(wmin:wmax)**2) + inyy*sum(wy0(wmin:wmax)**2) + inzz*sum(wz0(wmin:wmax)**2)
            
        forall (i = wmin:wmax)
            q10(i) = q10(i) +   q11(i) +   q12(i) +   q13(i) + q14(i)
            q11(i) = q11(i) + 2*q12(i) + 3*q13(i) + 4*q14(i)
            q12(i) = q12(i) + 3*q13(i) + 6*q14(i)
            q13(i) = q13(i) + 4*q14(i)
            q20(i) = q20(i) +   q21(i) +   q22(i) +   q23(i) + q24(i)
            q21(i) = q21(i) + 2*q22(i) + 3*q23(i) + 4*q24(i)
            q22(i) = q22(i) + 3*q23(i) + 6*q24(i)
            q23(i) = q23(i) + 4*q24(i)
            q30(i) = q30(i) +   q31(i) +   q32(i) +   q33(i) + q34(i)
            q31(i) = q31(i) + 2*q32(i) + 3*q33(i) + 4*q34(i)
            q32(i) = q32(i) + 3*q33(i) + 6*q34(i)
            q33(i) = q33(i) + 4*q34(i)
            q40(i) = q40(i) +   q41(i) +   q42(i) +   q43(i) + q44(i)
            q41(i) = q41(i) + 2*q42(i) + 3*q43(i) + 4*q44(i)
            q42(i) = q42(i) + 3*q43(i) + 6*q44(i)
            q43(i) = q43(i) + 4*q44(i)
         
            ! quaternion normalization (mere rescaling is not optimal!)
            cor1(i) = 1.0/sqrt(q10(i)**2 + q20(i)**2 + q30(i)**2 + q40(i)**2)
            q10(i) = q10(i)*cor1(i)
            q20(i) = q20(i)*cor1(i)
            q30(i) = q30(i)*cor1(i)
            q40(i) = q40(i)*cor1(i)
            
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

!        if (((.not.(all((abs(fx)+abs(fy)+abs(fz)) < 1e6))) .or. (.not.(all(abs(z0) < 50.0)))) .and. itime > 1) then
!            do i = 1, nmove
!                if(.not. ((abs(fx(i))+abs(fy(i))+abs(fz(i))) < 1e6)) print *, 'yf',itime, i, fx(i), fy(i), fz(i)
!                if(.not. (abs(z0(i)) < 50.0)) print *, 'yz', itime, i, z0(i)
!            end do
!            stop
!        end if
!        print *,'p'
        do i = wmin, wmax
!            cx0(i) = cx0(i) - lx*aint(cx0(i)*lxhi)
!            cy0(i) = cy0(i) - ly*aint(cy0(i)*lyhi)
            ix = 3*(i-1)
            forall (j = 1:3)
                x0(ix+j) = cx0(i) + a11(i)*xbody(j) + a21(i)*ybody(j) + a31(i)*zbody(j)
                y0(ix+j) = cy0(i) + a12(i)*xbody(j) + a22(i)*ybody(j) + a32(i)*zbody(j)
                z0(ix+j) = cz0(i) + a13(i)*xbody(j) + a23(i)*ybody(j) + a33(i)*zbody(j)
            end forall
            forall (j = 2:3)
                rsx(ix+j) = x0(ix+j) - x0(ix+1)
                rsy(ix+j) = y0(ix+j) - y0(ix+1)
                rsz(ix+j) = z0(ix+j) - z0(ix+1)
            end forall
        end do
!        print *, 'q'
!        if (((.not.(all((abs(fx)+abs(fy)+abs(fz)) < 1e6))) .or. (.not.(all((abs(tx)+abs(ty)+abs(tz)) < 1e6))) .or. (.not.(all(abs(z0) < 50.0)))) .and. itime > 1) then
!            do i = 1, nmove
!                if(.not. ((abs(fx(i))+abs(fy(i))+abs(fz(i))) < 1e6)) print *, 'zf',itime, i, fx(i), fy(i), fz(i)
!                if(.not. (abs(z0(i)) < 50.0)) print *, 'zz', itime, i, z0(i)
!                if(.not. (abs(z0(i)) < 50.0)) print *, 'zzx', sum(q10+q20+q30+q40),sum(abs(tx)+abs(ty)+abs(tz))
!            end do
!            stop
!        end if

        if (temper == -1.0) then        ! nve
            s0 = 0
            s1 = 0
            rs0 = 0
            rs1 = 0
        else                             ! nvt
            s0 = s0 +   s1 +   s2 +   s3 + s4
            s1 = s1 + 2*s2 + 3*s3 + 4*s4
            s2 = s2 + 3*s3 + 6*s4
            s3 = s3 + 4*s4

            rs0 = rs0 +   rs1 +   rs2 +   rs3 + rs4
            rs1 = rs1 + 2*rs2 + 3*rs3 + 4*rs4
            rs2 = rs2 + 3*rs3 + 6*rs4
            rs3 = rs3 + 4*rs4
        end if

        ke = kp + kw

        tstar = PH_R*temper/1000.0      ! kJ/mol
        qmass = real(3*(nfree+nrig))*tstar*2e3
        rmass = real(3*nrig)*tstar*1e3
        qmass = qmass/10.0
        rmass = rmass/10.0
!        print *, 'qmass', qmass/nrig
!        print *, 'rmass', rmass/nrig

        !  total hamiltonian
!        print *, ke, s0, rs0, uconf
        htot = 0.5*(ke + qmass*s0*s0 + rmass*rs0*rs0)
!        htot = htot + tstar*(real(3*(nfree + nrig) - 2)*log(s0) + real(3*nrig + 1)*log(rs0))
        htot = htot + uconf

        sft = (kp - real(3*(nfree+nrig) - 2)*tstar)/qmass
        sfr = (kw - real(3*nrig + 1)*tstar)/rmass
 
!        tempnow = kp*1000.0/(3.0*PH_R*real(nfree + nrig))
!        print *, 'trans', tempnow, sft
!        tempnow = kw*1000.0/(3.0*PH_R*real(nrig))
!        print *, 'rot', tempnow, sfr

!        print *, "mass: ", mi(tmin)

        return
    end subroutine move_predict


    subroutine move_correct
        implicit none
        integer*4 :: i
        real*8 :: rcorr, s0i, rs0i
        real*8, parameter :: f01 = 0.348611111, f11 = 1.0, f21 = 0.916666667, f31 = 0.333333333, f41 = 0.041666667
        real*8, parameter :: f02 = 0.158300000, f12 = 0.757000000, f22 = 1.0, f32 = 0.5, f42 = 0.083333333
        integer*4, save :: icnt = 0
        real*8 :: co

    ! translation 
!        s0i = 1.0/s0                        ! just auxiliaries
!        rcorr = 0.5*s0i**2
!        s0i = s1*s0i

        rcorr = 0.5
        s0i = s0
        forall (i = tmin:tmax)
            cort(i) = rcorr*ftx(i)*mii(i) - x1(i)*s0i - x2(i)
            x0(i) = x0(i) + f02*cort(i)
            x1(i) = x1(i) + f12*cort(i)
            x2(i) = x2(i) + f22*cort(i)
            x3(i) = x3(i) + f32*cort(i)
            x4(i) = x4(i) + f42*cort(i)
 
            cort(i) = rcorr*fty(i)*mii(i) - y1(i)*s0i - y2(i)
            y0(i) = y0(i) + f02*cort(i)
            y1(i) = y1(i) + f12*cort(i)
            y2(i) = y2(i) + f22*cort(i)
            y3(i) = y3(i) + f32*cort(i)
            y4(i) = y4(i) + f42*cort(i)
 
            cort(i) = rcorr*ftz(i)*mii(i) - z1(i)*s0i - z2(i)
            z0(i) = z0(i) + f02*cort(i)
            z1(i) = z1(i) + f12*cort(i)
            z2(i) = z2(i) + f22*cort(i)
            z3(i) = z3(i) + f32*cort(i)
            z4(i) = z4(i) + f42*cort(i)
        end forall

        ! rigmol (just for water)
        rcorr = rcorr/mw
        forall (i = wmin:wmax)
            cor1(i) = rcorr*fcx(i) - cx1(i)*s0i - cx2(i)
            cx0(i) = cx0(i) + f02*cor1(i)
            cx1(i) = cx1(i) + f12*cor1(i)
            cx2(i) = cx2(i) + f22*cor1(i)
            cx3(i) = cx3(i) + f32*cor1(i)
            cx4(i) = cx4(i) + f42*cor1(i)
 
            cor1(i) = rcorr*fcy(i) - cy1(i)*s0i - cy2(i)
            cy0(i) = cy0(i) + f02*cor1(i)
            cy1(i) = cy1(i) + f12*cor1(i)
            cy2(i) = cy2(i) + f22*cor1(i)
            cy3(i) = cy3(i) + f32*cor1(i)
            cy4(i) = cy4(i) + f42*cor1(i)
 
            cor1(i) = rcorr*fcz(i) - cz1(i)*s0i - cz2(i)
            cz0(i) = cz0(i) + f02*cor1(i)
            cz1(i) = cz1(i) + f12*cor1(i)
            cz2(i) = cz2(i) + f22*cor1(i)
            cz3(i) = cz3(i) + f32*cor1(i)
            cz4(i) = cz4(i) + f42*cor1(i)
        end forall

        forall (i = wmin:wmax)
            cor1(i) = 0.5*(-q30(i)*wx0(i) - q40(i)*wy0(i) + q20(i)*wz0(i)) !cq1
            cor2(i) = 0.5*( q40(i)*wx0(i) - q30(i)*wy0(i) - q10(i)*wz0(i)) !cq2
            cor3(i) = 0.5*( q10(i)*wx0(i) + q20(i)*wy0(i) + q40(i)*wz0(i)) !cq3
            cor4(i) = 0.5*(-q20(i)*wx0(i) + q10(i)*wy0(i) - q30(i)*wz0(i)) !cq4
 
            cor1(i) = cor1(i) - q11(i)
            q10(i) = q10(i)+f01*cor1(i)
            q11(i) = q11(i)+f11*cor1(i)
            q12(i) = q12(i)+f21*cor1(i)
            q13(i) = q13(i)+f31*cor1(i)
            q14(i) = q14(i)+f41*cor1(i)
 
            cor2(i) = cor2(i) - q21(i)
            q20(i) = q20(i)+f01*cor2(i)
            q21(i) = q21(i)+f11*cor2(i)
            q22(i) = q22(i)+f21*cor2(i)
            q23(i) = q23(i)+f31*cor2(i)
            q24(i) = q24(i)+f41*cor2(i)
 
            cor3(i) = cor3(i) - q31(i)
            q30(i) = q30(i)+f01*cor3(i)
            q31(i) = q31(i)+f11*cor3(i)
            q32(i) = q32(i)+f21*cor3(i)
            q33(i) = q33(i)+f31*cor3(i)
            q34(i) = q34(i)+f41*cor3(i)

            cor4(i) = cor4(i) - q41(i)
            q40(i) = q40(i)+f01*cor4(i)
            q41(i) = q41(i)+f11*cor4(i)
            q42(i) = q42(i)+f21*cor4(i)
            q43(i) = q43(i)+f31*cor4(i)
            q44(i) = q44(i)+f41*cor4(i)
        end forall

!        rs0i = 1.0/rs0                  ! auxiliaries
!        rcorr = rs0i**2
!        rs0i = 2.0*rs1*rs0i

        rcorr = 1.0
        rs0i = rs0
        forall (i = wmin:wmax)
            cor1(i) = a11(i)*tx(i) + a12(i)*ty(i) + a13(i)*tz(i)        !wxcorr
            cor1(i) = ((cor1(i)*rcorr + wy0(i)*wz0(i)*(inyy - inzz))*inxxi) - wx1(i) - rs0i*wx0(i)
            cor2(i) = a21(i)*tx(i) + a22(i)*ty(i) + a23(i)*tz(i)        !wycorr
            cor2(i) = ((cor2(i)*rcorr + wz0(i)*wx0(i)*(inzz - inxx))*inyyi) - wy1(i) - rs0i*wy0(i)
            cor3(i) = a31(i)*tx(i) + a32(i)*ty(i) + a33(i)*tz(i)        !wzcorr
            cor3(i) = ((cor3(i)*rcorr + wx0(i)*wy0(i)*(inxx - inyy))*inzzi) - wz1(i) - rs0i*wz0(i)
 
            wx0(i) = wx0(i) + f01*cor1(i)
            wx1(i) = wx1(i) + f11*cor1(i)
            wx2(i) = wx2(i) + f21*cor1(i)
            wx3(i) = wx3(i) + f31*cor1(i)
            wx4(i) = wx4(i) + f41*cor1(i)

            wy0(i) = wy0(i) + f01*cor2(i)
            wy1(i) = wy1(i) + f11*cor2(i)
            wy2(i) = wy2(i) + f21*cor2(i)
            wy3(i) = wy3(i) + f31*cor2(i)
            wy4(i) = wy4(i) + f41*cor2(i)

            wz0(i) = wz0(i) + f01*cor3(i)
            wz1(i) = wz1(i) + f11*cor3(i)
            wz2(i) = wz2(i) + f21*cor3(i)
            wz3(i) = wz3(i) + f31*cor3(i)
            wz4(i) = wz4(i) + f41*cor3(i)
        end forall
 
        if (temper /= -1.0) then
            rcorr = sft - s1
            s0 = s0 + f01*rcorr
            s1 = sft
            s2 = s2 + f21*rcorr
            s3 = s3 + f31*rcorr
            s4 = s4 + f41*rcorr
  
            rcorr = sfr - rs1
            rs0 = rs0 + f01*rcorr
            rs1 = sfr
            rs2 = rs2 + f21*rcorr
            rs3 = rs3 + f31*rcorr
            rs4 = rs4 + f41*rcorr
        end if

!        if ((.not.(all((abs(fx)+abs(fy)+abs(fz)) < 1e6))) .or. (.not.(all(abs(z0) < 50.0)))) then
!            do i = 1, nmove
!                if(.not. ((abs(fx(i))+abs(fy(i))+abs(fz(i))) < 1e6)) print *, 'mf', i, fx(i), fy(i), fz(i)
!                if(.not. (abs(z0(i)) < 50.0)) print *, 'mz', i,z0(i)
!            end do
!            stop
!        end if
!        print *, 'r'
        if (nsurf > 0) call move_shake

        ! instantaneous temperature
        tempnow = ke*1000.0/(3.0*PH_R*real(nfree + 2*nrig))
!        print *, 'srs', s0, rs0
 
        return
    end subroutine move_correct

    subroutine move_shake_init
        implicit none

        ismin = nliqa + 1
        ismax = nmove

        allocate(moving(ismin:ntota), moved(ismin:ntota))
        moved(ismax+1:ntota) = 1

        xold(ismin:ntota) = x0(ismin:ntota)
        yold(ismin:ntota) = y0(ismin:ntota)
        zold(ismin:ntota) = z0(ismin:ntota)

        maxit = 10
        tol2 = 1.e-8
        rptol = 0.85

        nfreed = 3*(nmove - nliqa) - sum(nt(nltype+1:nstype)*nconst(nltype+1:nstype))
        print *, 'nrfreed', nfreed, nltype, nstype, nconst(nltype+1:nstype)

        return
    end subroutine move_shake_init
    !
    subroutine move_shake
        implicit none
        logical :: done
        integer*4 :: i, it, ia, ib, ite, ic, is
        real*8 :: aux, kaux, pxab, pyab, pzab, rxab, ryab, rzab, rabsq, diffsq, rma, rmb, rpab

        ! basic verlet + kinetic energy
        kps = 0.0
        do i = ismin, ismax
            it = itype(i)
            aux = xold(i)
            xold(i) = x0(i)
            aux = x0(i) - aux
            aux = aux - lx*anint(aux*lxi)
            x0(i) = x0(i) + aux + ftx(i)*mii(i)
            kaux = aux*aux

            aux = yold(i)
            yold(i) = y0(i)
            aux = y0(i) - aux
            aux = aux - ly*anint(aux*lyi)
            y0(i) = y0(i) + aux + fty(i)*mii(i)
            kaux = kaux + aux*aux

            aux = zold(i)
            zold(i) = z0(i)
            aux = z0(i) - aux
            z0(i) = z0(i) + aux + ftz(i)*mii(i)
            kaux = kaux + aux*aux

            kps(it) = kps(it) + mi(i)*kaux
        end do
!        print *, 'ft', ftx(ismin+1)*mii(ismin+1), fty(ismin+1)*mii(ismin+1), ftz(ismin+1)*mii(ismin+1)
!        print *, 'fb', ftx(ismax)*mii(ismin+1), fty(ismax)*mii(ismin+1), ftz(ismax)*mii(ismin+1)

!        do i = 1, nit
!            if (kps(i) /= 0.0) print *,'kp', i, kps(i)*1000.0/(PH_R*3.0*real((ismax-ismin+1)/4))
!        end do
    
        tkps = sum(kps)*1000.0/(PH_R*real(nfreed))

        moving(ismin:ismax) = 1
        moved(ismin:ismax) = 0

        ! apply all constraints at once
        done = .false.
        ite = 0
        do 
            if (done .or. (ite > maxit)) exit
            do ic = 1, ncons
                done = .true.
                it = ictype(ic)
                do is = 1, nconst(it)
                    ia = constra(is, ic)
                    ib = constrb(is, ic)
 
                    if (moved(ia)*moved(ib) == 0) then
                 
                        pxab = x0(ia) - x0(ib)
                        pxab = pxab - lx*anint(pxab*lxi)
                        pyab = y0(ia) - y0(ib)
                        pyab = pyab - ly*anint(pyab*lyi)
                        pzab = z0(ia) - z0(ib)
                        pzab = pzab - lz*anint(pzab*lzi)
                 
                        rabsq = length(it, is)**2
                        diffsq = rabsq - (pxab*pxab + pyab*pyab + pzab*pzab)
                 
                        if (abs(diffsq) > rabsq*tol2) then
                 
                            rxab = xold(ia) - xold(ib)
                            rxab = rxab - lx*anint(rxab*lxi)
                            ryab = yold(ia) - yold(ib)
                            ryab = ryab - ly*anint(ryab*lyi)
                            rzab = zold(ia) - zold(ib)
                            rzab = rzab - lz*anint(rzab*lzi)

                            rpab = rxab*pxab + ryab*pyab + rzab*pzab

                            if (rpab < rabsq*rptol) then
                                write(*, *) ia, ib
                                write(*, *) x0(ia), x0(ib)
                                write(*, *) x0(ia+1), x0(ib+1)
                                write(*, *) y0(ia), y0(ib)
                                write(*, *) z0(ia), z0(ib)
                                write(*, *) sqrt(pxab*pxab + pyab*pyab + pzab*pzab)
                                write(*, *) rpab, rabsq, rptol
                                stop 'Constraint failure'
                            end if

                            rma = mii(ia)
                            rmb = mii(ib)
                            aux = diffsq/(2.0*(rma + rmb)*rpab)
                            rma = rma*aux
                            rmb = rmb*aux

                            x0(ia) = x0(ia) + rxab*rma
                            y0(ia) = y0(ia) + ryab*rma
                            z0(ia) = z0(ia) + rzab*rma

                            x0(ib) = x0(ib) - rxab*rmb
                            y0(ib) = y0(ib) - ryab*rmb
                            z0(ib) = z0(ib) - rzab*rmb

                            moving(ia) = 0
                            moving(ib) = 0
                            done = .false.
                        end if
                    end if
                end do
            end do

            ite = ite + 1
            moved(ismin:ismax) = moving(ismin:ismax)
            moving(ismin:ismax) = 1
        end do

!        print *, 'ite', ite
        if (.not. done) stop 'shake(): Too many constraint iterations!'

        return
    end subroutine move_shake

end module move

