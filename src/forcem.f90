!----------------------------------------------------------------------------
!       Module CONFIG contains molecular configuration
!       and functions used for its input and output
!--------------------------------------------------------------------------

module force  ! state of the system : configuration, energy, forces, ...
    use nnlist
    use measure
    implicit none

    ! rational approximation for the erf(x), abramowitz & stegun, eq 7.1.26
    real*8, parameter, private :: a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741
    real*8, parameter, private :: a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911

    ! basic variables
    real*8, dimension(:), allocatable :: fx, fy, fz
    real*8, dimension(:), allocatable :: ftx, fty, ftz, fcx, fcy, fcz, tx, ty, tz, btx, bty, btz
    real*8, dimension(-1:1) :: fwsx, fwsy, fwsz
    real*8, dimension(:), allocatable :: stat, bufs  ! mpi buffer for statistics (ene, pressure, ...)
    real*8, dimension(:), allocatable :: elj, eqq
    real*8 :: alfa, rcutsq, spii, facq, press, uconf, uewa, urf, uself, ewqWW, ewqLL, ewqWL, uselfw, ulr, vir
    real*8 :: totpres
    integer*4 :: iamin, iamax
    ! fouriere part
    integer*8 :: kmax, klmax, kmmax, knmax, iqmax ! number of k-vectors and charges
    integer*4 :: n1min, n1max, n2min, n2max
    complex*16, dimension(:, :), allocatable :: el, em, en           ! exp(i2Pkr/L)
    complex*16, dimension(:), allocatable :: elm, elmn           ! el*em, elm*en
    complex*16, dimension(:), allocatable :: Qsum
    real*8, dimension(:), allocatable :: ft   ! fourier sum
    complex*16, dimension(:), allocatable :: Qsumw1, Qsumw2   ! wall Qsum, Ak
    real*8, dimension(:), allocatable :: Ak
    integer*4, dimension(:), allocatable :: myk
    real*8 :: pipix, pipiy, pipiz 
    complex*16 :: cpipix, cpipiy, cpipiz
    integer*4 :: iwmin, iwmax
    real*8 :: crf
    ! long range correction
    real*8, dimension(:), allocatable :: flr_l, flr_w, flr_lja, flr_ljc, flr_bua, flr_buc, flr_bur
    real*8, dimension(:), allocatable :: ulr_l, ulr_w, flrcor
    ! external force
    integer*4 :: iemin, iemax

    ! average energy of water
    real*8 :: uavex, uavextot
    integer*4, dimension(:), allocatable :: firstbin


    contains

!
!======================== Molecular forces ==================================
!
    subroutine force_init
        implicit none
        integer*4 :: iw

        spii= 2.0/sqrt(M_PI)
        alfa = 3.09/rcut
        facq = PH_4epi*PH_e*PH_e*PH_Na/(1000.0*PH_A2m)

        totpres = 0

        iw = nmove/nproc
        iwmin = iw*me + 1
        iwmax = iw*(me + 1)
        if (me == nproc-1) iwmax = nmove

        call force_allocate
        call force_fouri_init
        call force_real_init
        call force_intra_init
!        if(lfext) call force_extern_init

        print *, 'qsum', sum(qs)
        print *, 'nwall', nwall, ngogo, npull, nstatic

        call force_zero

        return
    end subroutine force_init

    ! zero accumulators
    subroutine force_zero
        implicit none

        uavextot = 0.0

        return
    end subroutine force_zero

    !
    subroutine force_get
        implicit none

        fx = 0.0
        fy = 0.0
        fz = 0.0
        fwsx = 0.0
        fwsy = 0.0
        fwsz = 0.0

        vir = 0.0
        elj = 0.0                       ! LJ energy
        eqq  = 0.0                      ! coulombic energy
        uewa = 0.0
 
!        print *, 'f'
!        print *, 'sum1', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        call force_real
        !write(*,'(i1,3(1x,f10.3))') 1,sum(abs(fx)), sum(abs(fy)), sum(abs(fz)) 
!        print *, 'g'
!        print *, 'sum2', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        call force_fouri
        !write(*,'(i1,3(1x,f10.3))') 3,sum(abs(fx)), sum(abs(fy)), sum(abs(fz)) 
!        print *, 'i'
!        print *, 'sum3', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        call force_intra
!        print *, 'sum4', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        !write(*,'(i1,3(1x,f10.3))') 2,sum(abs(fx)), sum(abs(fy)), sum(abs(fz)) 
        !print *, 'h'
!        if (lfext) call force_extern
!        print *, 'j'
        call force_sum           ! gathers, sums, and broadcasts
!        print *, 'k'
        if (lchempot) call chempot

        return
    end subroutine force_get

    subroutine force_allocate
        implicit none

        allocate(fx(nmove), fy(nmove), fz(nmove))
        allocate(ftx(nriga+1:nmove), fty(nriga+1:nmove), ftz(nriga+1:nmove))
        allocate(fcx(nrig), fcy(nrig), fcz(nrig))
        allocate(tx(nrig), ty(nrig), tz(nrig))
        allocate(btx(nrig), bty(nrig), btz(nrig))

        allocate(elj(nijt), eqq(nijt))

        allocate(flr_l(nit), flr_w(nit))
        allocate(flr_lja(nit), flr_ljc(nit))
        allocate(flr_bua(nit), flr_buc(nit), flr_bur(nit))
        allocate(ulr_l(nit), ulr_w(nit), flrcor(nit))

        allocate(stat(2*nijt+10), bufs(2*nijt+10))

        allocate(firstbin(nmove))


        return
    end subroutine force_allocate
    !
    !====================== Forces computed in real space =================
    !
    subroutine force_real_init
        implicit none
        integer*4 :: it, jt, ijt, i
        real*8 :: rhowall, rholiq, flr, ul, xdens(TMX), sgn, ro

        rcutsq = rcut**2

        ! nearest neighbor list initialization
        call nnlist_init

        ! long range corrections

        ! site densities of liquid and solid phases
        rholiq = nmove/(lx*ly*lw)       ! liquid site densities
        rhowall = sum(nss(ntype, 1:nst(ntype)))/(2*cx*cy*cz)
        ! mole fractions of site types
        xdens = 0.0
        do i = 1, nmove
            xdens(itype(i)) = xdens(itype(i)) + 1.0
        end do
        xdens(1:nit-nst(ntype)) = xdens(1:nit-nst(ntype))/real(nmove)*rholiq
        xdens(nit-1) = real(nss(ntype, 1))/real(sum(nss(ntype, 1:2)))*rhowall
        xdens(nit) = real(nss(ntype, 2))/real(sum(nss(ntype, 1:2)))*rhowall

        ! lrcorr for wall
        flr_l = 0.0 ; flr_w = 0.0
        flr_lja = 0.0 ; flr_ljc = 0.0
        flr_bua = 0.0 ; flr_buc = 0.0 ; flr_bur = 1.0
        do it = 1, nit-nst(nwtype)                ! moving types
            do jt = 1, nit          ! wall types
                ijt = ijtype(it, jt)
                sgn = 1.0
                if (jt <= nit-nst(ntype)) sgn = -1.0
                select case (vdwnum(ijt))
                    case (0)
                        flr = 0.0
                        ul = 0.0
                    case (1)                          ! lennard-jones
                        flr = eps(ijt)*sig(ijt)*(sig(ijt)/(10.0*rcut**10) - 1.0/(4.0*rcut**4))
                        ul  = eps(ijt)*sig(ijt)*(sig(ijt)/(90.0*rcut**9) - 1.0/(12.0*rcut**3))

                        flr_lja(it) = flr_lja(it) + sgn*xdens(jt)*eps(ijt)*sig(ijt)**2/10.0
                        flr_ljc(it) = flr_ljc(it) + sgn*xdens(jt)*eps(ijt)*sig(ijt)/4.0
                    case (2)                          ! Buckingham
                        ro = 1.0/rhoi(ijt)
                        flr = ro*(ro + rcut)*aa(ijt)*exp(-rcut/ro) - cc(ijt)/(4.0*rcut**4)
                        ul  = ro*ro*(2.0*ro + rcut)*aa(ijt)*exp(-rcut/ro) - cc(ijt)/(12.0*rcut**3)

                        flr_bua(it) = flr_bua(it) + sgn*xdens(jt)*aa(ijt)*ro
                        flr_buc(it) = flr_buc(it) + sgn*xdens(jt)*cc(ijt)/4.0
                        flr_bur(it) = ro
                    case default                      ! no interaction
                        write(*, *) it, jt, ijt, vdwnum(ijt)
                        stop 'force_real_init(): wrong interaction type'
                end select
                if (sgn < 0.0) then
                   flr_l(it) = flr_l(it) + flr*xdens(jt)
                   ulr_l(it) = ulr_l(it) + ul*xdens(jt)
                else
                   flr_w(it) = flr_w(it) + flr*xdens(jt)
                   ulr_w(it) = ulr_w(it) + ul*xdens(jt)
                end if
            end do
        end do
        flr_lja = flr_lja*2.0*M_PI
        flr_ljc = flr_ljc*2.0*M_PI
        flr_bua = flr_bua*2.0*M_PI
        flr_buc = flr_buc*2.0*M_PI
        flr_l = flr_l*2.0*M_PI
        flr_w = flr_w*2.0*M_PI
        ulr_l = ulr_l*2.0*M_PI
        ulr_w = ulr_w*2.0*M_PI

        return
    end subroutine force_real_init
    !
    subroutine force_real
        implicit none
        integer*4 :: i, j, it, jt, k, ijt, nshl, ix, jx, jmin, sgn, ii, l
        real*8 :: aux, fij, xdist, ydist, zdist, rij, xi, yi, zi, qx, rsq, rsqi
        real*8 :: alfar, t, exp2a, erfx, r6, xio, yio, zio, dw, dwa
        real*8 :: layer
        logical :: luin
        include 'mpif.h'

        !ibcount = 0
      ! update of the neighbor list
      ! add movement of the wall
        aux = max(maxval((ocx0 - cx0)**2 + (ocy0 - cy0)**2 + (ocz0 - cz0)**2),          &
                  maxval((ox0(nriga+1:nmove) - x0(nriga+1:nmove))**2                    &
                       + (oy0(nriga+1:nmove) - y0(nriga+1:nmove))**2                    &
                       + (oz0(nriga+1:nmove) - z0(nriga+1:nmove))**2),                  &
                  dshsq)
        if (4.0*aux > skin*skin) call nnlist_make

!        aux = max(maxval((ox0(1:nriga:3) - x0(1:nriga:3))**2                            &
!                       + (oy0(1:nriga:3) - y0(1:nriga:3))**2                            &    
!                       + (oz0(1:nriga:3) - z0(1:nriga:3))**2),                          &
!                  maxval((ox0(nriga+1:nmove) - x0(nriga+1:nmove))**2                    &
!                       + (oy0(nriga+1:nmove) - y0(nriga+1:nmove))**2                    &
!                       + (oz0(nriga+1:nmove) - z0(nriga+1:nmove))**2))
!        if (4.0*aux > skin*skin) call nnlist_make


        uavex = 0.0
        luin = .true.

!        firstbin(1:nliqa) = 0
!        if (nsurf == 0) then
!            layer = lwh - 3.0
!            k = 0
!            do i = 1, nriga, 3
!                if (abs(z0(i)) > layer) then
!                    k = k + 1
!                    firstbin(i) = 1
!                end if
!            end do
! 
!            if (me == 0) then
!                if (k < 144) then ! only single processor
!                    print *, k, 'Too few molecules'
!                else if (k > 144) then
!                    print *, k, 'Too many molecules'
!                end if
!            end if
!        else
!            firstbin(nliqa+1:nmove) = 1
!        end if

      ! interactions between moving particles
!        print *, 'sumx', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        do k = 1, ipack
            i = list1(k)
            j = list2(k)
 
            xio = x0(i) - x0(j)
            xio = xio - lx*anint(xio*lxi)
            yio = y0(i) - y0(j)
            yio = yio - ly*anint(yio*lyi)
            zio = z0(i) - z0(j)
!            zio = zio - lz*anint(zio*lzi)
            rsq = xio*xio + yio*yio + zio*zio
            if (rsq > rcutsq) cycle
 
            it = itype(i)
            jt = itype(j)
            ijt = ijtype(it, jt)

            rsqi = 1.0/rsq
            rij = sqrt(rsq)

!                luin = .false.
!                if (firstbin(i)*firstbin(j) == 0) luin = .true.

            ! vdw
            select case (vdwnum(ijt))
                case (0)
                    fij = 0.0                     ! no vdw
                case (1)                          ! lennard-jones
                    r6 = sig(ijt)*rsqi**3
                    elj(ijt) = elj(ijt) + eps(ijt)*r6*(r6 - 1.0)
                        if (luin) uavex = uavex + eps(ijt)*r6*(r6 - 1.0)
                    fij = 6.0*rsqi*eps(ijt)*r6*(2.0*r6 - 1.0)
                case (2)                          ! Buckingham
                    r6 = cc(ijt)*rsqi**3
                    aux = aa(ijt)*exp(-rij*rhoi(ijt))
                    elj(ijt) = elj(ijt) + aux - r6
                        if (luin) uavex = uavex + aux - r6
                    fij = aux*rhoi(ijt)/rij - 6.0*rsqi*r6
                case default
                    write(*, *) i, j, ijt, vdwnum(ijt), vdwtyp(ijt)
                    stop 'force_real(): unknown interaction type (f-f)!'
            end select

            ! except pairs of molecules inside the solution
            ! coulombic 1-1
            aux = facq*qs(i)*qs(j)/rij
            alfar = alfa*rij
            t = 1.0/(1.0 + p*alfar)
            exp2a = exp(-alfar*alfar)
            erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp2a
            eqq(ijt) = eqq(ijt) + erfx*aux
                if (luin) uavex = uavex + erfx*aux
            fij = fij + (spii*alfar*exp2a + erfx)*aux*rsqi
            fx(i) = fx(i) + xio*fij
            fy(i) = fy(i) + yio*fij
            fz(i) = fz(i) + zio*fij
            fx(j) = fx(j) - xio*fij
            fy(j) = fy(j) - yio*fij
            fz(j) = fz(j) - zio*fij
            if (it == 3 .or. it == 4) then
                l = int((z0(i)+lwh)*rdeli + 0.5)
                efipos(it-2, l) = efipos(it-2, l) + zio*fij
            end if
            if (jt == 3 .or. jt == 4) then
                l = int((z0(j)+lwh)*rdeli + 0.5)
                efipos(jt-2, l) = efipos(jt-2, l) - zio*fij
            end if

!            l = int(rij*rdeli + 0.5)
!            ngc(ijt, l) = ngc(ijt, l) + 1 

            ! coulombic i /= 1 or j /= 1
            jmin = j + 1             ! don't cycle over the first pair already computed
            do ix = i, i+nsi(it)-1           ! coulombic forces and energy (and virial)
                xi = xio + rsx(ix)
                yi = yio + rsy(ix)
                zi = zio + rsz(ix)
                qx = facq*qs(ix)
                do jx = jmin, j+nsi(jt)-1           ! coulombic forces and energy (and virial)
                    xdist = xi - rsx(jx)
                    ydist = yi - rsy(jx)
                    zdist = zi - rsz(jx)
                    rsq = xdist*xdist + ydist*ydist + zdist*zdist
                    rsqi = 1.0/rsq
                    rij = sqrt(rsq)
                    aux = qx*qs(jx)/rij
                    alfar = alfa*rij
                    t = 1.0/(1.0 + p*alfar)
                    exp2a = exp(-alfar*alfar)
                    erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp2a
                    eqq(ijt) = eqq(ijt) + erfx*aux
                        if (luin) uavex = uavex + erfx*aux
                    fij = (spii*alfar*exp2a + erfx)*aux*rsqi
                    fx(ix) = fx(ix) + xdist*fij
                    fy(ix) = fy(ix) + ydist*fij
                    fz(ix) = fz(ix) + zdist*fij
                    fx(jx) = fx(jx) - xdist*fij
                    fy(jx) = fy(jx) - ydist*fij
                    fz(jx) = fz(jx) - zdist*fij

                    if (it == 3 .or. it == 4) then
                        l = int((z0(i)+lwh)*rdeli + 0.5)
                        efipos(it-2, l) = efipos(it-2, l) + zdist*fij
                    end if
                    if (jt == 3 .or. jt == 4) then
                        l = int((z0(j)+lwh)*rdeli + 0.5)
                        efipos(jt-2, l) = efipos(jt-2, l) - zdist*fij
                    end if

                    !l = int(rij*rdeli + 0.5)
                    !ii = ijtype(itype(ix), itype(jx))
                    !ngc(ii, l) = ngc(ii, l) + 1 
                end do
                jmin = j
            end do
        end do
        !write(*,'(i1,3(1x,f10.3))') 4,sum(abs(fx)), sum(abs(fy)), sum(abs(fz)) 
!        print *, 'sumy', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))


        do i = iwmin, iwmax
             if (i > nliqa) then            
                 sgn = int(sign(1.d0, z0(i)))
                 fwsz(sgn) = fwsz(sgn) + fz(i)
             end if
        end do
!        print *, 'sumz', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
 
    ! interactions between moving and frozen particles

        do k = 1, wpack
            i = listw1(k)
            j = listw2(k)

            xi = x0(i) - x0(j)
            xi = xi - lx*anint(xi*lxi)
            yi = y0(i) - y0(j)
            yi = yi - ly*anint(yi*lyi)
            zi = z0(i) - z0(j)
            !zi = zi - lz*anint(zi*lzi)
            rsq = xi*xi + yi*yi + zi*zi
            if (rsq > rcutsq) cycle
            sgn = int(sign(1.d0, z0(j)))      ! ???

            it = itype(i)
            ijt = ijtype(it, itype(j))

!                luin = .false.
!                if (firstbin(i) == 0) luin = .true.

            select case (vdwnum(ijt))
                case (0)                          ! no vdw
                    fij = 0.0
                case (1)                          ! lennard-jones
                    rsqi = 1.0/rsq
                    r6 = sig(ijt)*rsqi**3
                    elj(ijt) = elj(ijt) + eps(ijt)*r6*(r6 - 1.0)
                        if (luin) uavex = uavex + eps(ijt)*r6*(r6 - 1.0)
                    fij = 6.0*rsqi*eps(ijt)*r6*(2.0*r6 - 1.0)
                case (2)                          ! Buckingham
                    rsqi = 1.0/rsq
                    rij = sqrt(rsq)
                    r6 = cc(ijt)*rsqi**3
                    aux = aa(ijt)*exp(-rij*rhoi(ijt))
                    elj(ijt) = elj(ijt) + aux - r6
                        if (luin) uavex = uavex + aux - r6
                    fij = aux*rhoi(ijt)/rij - 6.0*rsqi*r6
                case default                      ! no interaction
                    stop 'force_real(): unknown interaction type (f-s)!'
            end select

            qx = facq*qs(j)
            do ix = i, i+nsi(it)-1           ! coulombic forces and energy (and virial)
                xdist = xi + rsx(ix)
                ydist = yi + rsy(ix)
                zdist = zi + rsz(ix)
                rsq = xdist*xdist + ydist*ydist + zdist*zdist

                rij = sqrt(rsq)
                rsqi = 1.0/rsq
                aux = qx*qs(ix)/rij
                alfar = alfa*rij
                t = 1.0/(1.0 + p*alfar)
                exp2a = exp(-alfar*alfar)
                erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp2a

                eqq(ijt) = eqq(ijt) + erfx*aux
                        if (luin) uavex = uavex + erfx*aux
                fij = fij + (spii*alfar*exp2a + erfx)*aux*rsqi

                fx(ix) = fx(ix) + xdist*fij
                fy(ix) = fy(ix) + ydist*fij
                fz(ix) = fz(ix) + zdist*fij
                if (it <= nltype) then
                fwsx(sgn) = fwsx(sgn) - xdist*fij
                fwsy(sgn) = fwsy(sgn) - ydist*fij
                fwsz(sgn) = fwsz(sgn) - zdist*fij
                end if
                if (it == 3 .or. it == 4) then
                    l = int((z0(i)+lwh)*rdeli + 0.5)
                    efipos(it-2, l) = efipos(it-2, l) + zdist*fij
                end if
                !print *, sgn, zdist*fij
                fij = 0.0

                l = int(rij*rdeli + 0.5)
                ii = ijtype(itype(ix), itype(j))
                ngc(ii, l) = ngc(ii, l) + 1 
        !        if (ii == 7 .and. rij < 2.375) ibcount = ibcount + 1
            end do
        end do
        !print *, 'ibcount', ibcount, sum(ngc(7, 1:119))
!        print *,'sum0', sum(abs(fx(1:nmove))), sum(abs(fy(1:nmove))), sum(abs(fz(1:nmove)))
        !print *, 1,fwsz(1)
        !print *, -1,fwsz(-1)
        !print *, nltype, fwsz(:)
!        print *, 'sumv', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        
        !write(*,'(i1,3(1x,f10.3))') 5,sum(abs(fx)), sum(abs(fy)), sum(abs(fz)) 

        ! long range correction to liquid (+surface) - wall interactions in z-direction
        ulr = 0.0
        do i = iwmin, iwmax     ! only my molecules
            it = itype(i)
            if (abs(qs(i)) < 0.7) cycle
            do j = -1, 1, 2
                aux = 0.0
                sgn = dble(j)
                dw = lwh - sgn*z0(i)
                if (dw < rcut) then
                    aux = aux + sgn*flr_w(it)
                    ulr = ulr + ulr_w(it)
!                    if (firstbin(i) == 0) uavex = uavex + ulr_w(it)
!                    uavex = uavex + ulr_w(it)
                else
                    aux = aux + sgn*flr_l(it)
                    ulr = ulr + ulr_l(it)

                    dwa = 1.0/(dw*dw)
                    aux = aux + sgn*(flr_lja(it)*dwa**5 - flr_ljc(it)*dwa**2) ! water rcut + wall 
                    ulr = ulr + (flr_lja(it)/(9.0*dw**9) - flr_ljc(it)/(3.0*dw**3))

                    t = flr_bur(it)
                    exp2a = flr_bua(it)*exp(-dw/t)
!                    print *, 't', t, dw, flr_bua(it)
                    aux = aux + sgn*((t + dw)*exp2a - flr_buc(it)*dwa**2)
                    ulr = ulr + (t*(2.0*t + dw)*exp2a - flr_buc(it)/(3.0*dw**3))
!                    if (firstbin(i) == 0) uavex = uavex + (t*(2.0*t + dw)*exp2a - flr_buc(it)/(3.0*dw**3))
!                    uavex = uavex + (t*(2.0*t + dw)*exp2a - flr_buc(it)/(3.0*dw**3))
                end if
                fz(i) = fz(i) - aux
                !fwsz(j) = fwsz(j) + aux
                if (it <= nltype) fwsz(j) = fwsz(j) + aux
            end do
        end do
        uavex = uavex + ulr
        if (.not.((abs(elj(8)) >= 0.0) .and. (abs(elj(8)) < 1e5))) then
            print *,'w', elj(8)
!            stop 'w'
        end if
!        print *, 'sumw', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
!        print *,'sum1', sum(abs(fx(1:nmove))), sum(abs(fy(1:nmove))), sum(abs(fz(1:nmove)))

        return
    end subroutine force_real

    subroutine force_intra_init 
        implicit none
        integer*4 :: ia

        ia = nang/nproc
        iamin = ia*me + 1
        iamax = ia*(me + 1)
        if (me == nproc-1) iamax = nang

        return
    end subroutine force_intra_init 
    subroutine force_intra 
        implicit none
        integer*4 :: ia, i, j, k, is, it, l, sgn
        real*8 :: ax, ay, az, bx, by, bz, fij, tmp, ai, bi, ab, cang, aa, bb

        do ia = iamin, iamax     ! go through all my angles
            it = iatype(ia)
            do is = 1, nangt(it)
                i = ang1(is, ia)
                j = ang2(is, ia)
                k = ang3(is, ia)

                ax = x0(k) - x0(j)          ! h - o
                ax = ax - lx*anint(ax*lxi)
                ay = y0(k) - y0(j)
                ay = ay - ly*anint(ay*lyi)
                az = z0(k) - z0(j)
                !az = az - lz*anint(az*lzi)

                bx = x0(j) - x0(i)          ! me - o
                bx = bx - lx*anint(bx*lxi)
                by = y0(j) - y0(i)
                by = by - ly*anint(by*lyi)
                bz = z0(j) - z0(i)
                !bz = bz - lz*anint(bz*lzi)

                ai = 1.0/sqrt(ax*ax + ay*ay + az*az) ! o-h
                bi = 1.0/sqrt(bx*bx + by*by + bz*bz)
!                if (it == 4) print *, j, k, 1.0/ai, 1.0/bi 
                ab = ax*bx + ay*by + az*bz
                cang = -ab*ai*bi
                fij = kt(it, is)*(acos(cang) - theta(it, is))/sqrt(1.0 - cang*cang)
                fij = fij*ai*bi

                sgn = int(sign(1.d0, z0(i)))      ! ???
                bb = ab*bi*bi ! m-o
                aa = ab*ai*ai
                tmp = aa*ax - bx
                fx(k) = fx(k) + fij*tmp         ! just for H
                fx(j) = fx(j) + fij*(bb*bx - ax - tmp)
!                fwsx(sgn) = fwsx(sgn) - fij*(bb*bx - ax)
                tmp = aa*ay - by
                fy(k) = fy(k) + fij*tmp         ! just for H
                fy(j) = fy(j) + fij*(bb*by - ay - tmp)
!                fwsy(sgn) = fwsy(sgn) - fij*(bb*by - ay)
                tmp = aa*az - bz 
                fz(k) = fz(k) + fij*tmp         ! just for H
                fz(j) = fz(j) + fij*(bb*bz - az - tmp)
!                fwsz(sgn) = fwsz(sgn) - fij*(bb*bz - az)

                !print *, 'i', fij*(bb*bz - az)

                l = int((cang + 1.0)*100) + 1
                if (l < 1) then
                    print *, 'angle', l, cang
                    print *, 'a1', i, j, k, cang
                    print *, 'ai', x0(i), y0(i), z0(i)
                    print *, 'aj', x0(j), y0(j), z0(j)
                    print *, 'ak', x0(k), y0(k), z0(k)
                end if
                ngang(it, l) = ngang(it, l) + 1

                ! forces along the constraint bonds (needed for pressure)
                ! H-O
!                tmp = -fx(k)*ax*ai
!                fx(k) = fx(k) - tmp
!                fx(j) = fx(j) + tmp
!                tmp = -fy(k)*ay*ai
!                fy(k) = fy(k) - tmp
!                fy(j) = fy(j) + tmp
!                tmp = -fz(k)*az*ai
!                fz(k) = fz(k) - tmp
!                fz(j) = fz(j) + tmp
!                ! O-Me
!                k = int(sign(1.d0, z0(i)))      ! ???
!                tmp = -fx(j)*bx*bi
!                fx(j) = fx(j) - tmp
!                fwsx(k) = fwsx(k) + tmp
!                tmp = -fy(j)*by*bi
!                fy(j) = fy(j) - tmp
!                fwsy(k) = fwsy(k) + tmp
!                tmp = -fz(j)*bz*bi
!                fz(j) = fz(j) - tmp
!                fwsz(k) = fwsz(k) + tmp
            end do
        end do
        return
    end subroutine force_intra 

    !
    !================== Ewald summation - reciprocal space ====================
    !
    subroutine force_fouri_init
        implicit none
        integer*4 :: iq, kl, km, kn, kmmin, knmin, it, is, js, i, j, kcount, kcount2, ix, vkmax, namax
        real*8 :: alf4i, fac, rl, rm, rn, rksq, rij, t, erfv, ksq, rkmaxsq, xi, yi, zi
        complex*16 :: Qws

        ! frequently used quantities
        pipix = 2.0*M_PI/Lx
        pipiy = 2.0*M_PI/Ly
        pipiz = 2.0*M_PI/Lz
        cpipix = (0.0, 1.0)*pipix
        cpipiy = (0.0, 1.0)*pipiy
        cpipiz = (0.0, 1.0)*pipiz
        alf4i = -1.0/(4.0*alfa**2)
        fac = 2.0*PH_e*PH_e*PH_Na/(PH_eps0*PH_A2m*1000.0*Lx*Ly*Lz)

        klmax = kmax
        klmax = nint(real(kmax)*lx/min(lx,ly,lz))
        kmmax = kmax
        kmmax = nint(real(kmax)*ly/min(lx,ly,lz))
        knmax = kmax
        knmax = nint(real(kmax)*lz/min(lx,ly,lz))
        rkmaxsq = (real(kmax)*max(pipix, pipiy, pipiz))**2 + 1.e-4  ! to avoid rounding errors for 5+1+1

        vkmax = int((4.0/3.0)*M_PI*real((kmax+1)*(2*kmmax+1)*(2*knmax+1)/4))
        if (nwall > ngogo) print *, 'force_fouri_init(): nwall > nmove'
        namax = max(nwall, ngogo)
        print *, 'namax', namax, nwall, ngogo
        allocate(el(namax,0:kmax), em(namax,-kmax:kmax), en(namax,-knmax:knmax))
        allocate(elm(namax), elmn(namax), ft(nmove))
        allocate(Qsumw1(vkmax), Qsumw2(vkmax), Ak(vkmax), myk(vkmax))
        allocate(Qsum(vkmax))

        ! moving vs. static wall atoms - define moving atoms' range
        n1min = nmove + 1
        n1max = nmove + nwall/2
        n2min = nmove + nwall/2 + 1
        n2max = nmove + nwall

        if (ngogo <= n1max) then
            n2max = n1max
            n1max = ngogo
        else
            n2max = ngogo
        end if
        print *, 'ngogo', ngogo, n1min, n1max, n2min, n2max

        ! the zeroth element (always the same)
        el(:, 0) = 1.0
        em(:, 0) = 1.0
        en(:, 0) = 1.0

        forall (i = 1:nwall)
            el(i, 1) = exp(cpipix*(x0(i+nmove)))
            em(i, 1) = exp(cpipiy*(y0(i+nmove)))
            en(i, 1) = exp(cpipiz*(z0(i+nmove)))
        end forall

        ! sum over the wall (particles that do not move) & precomputations
        do kl = 2, klmax
            forall (iq = 1:nwall) el(iq, kl) = el(iq, kl-1)*el(iq, 1)
        end do
        do km = 2, kmmax
            forall (iq = 1:nwall) em(iq, km) = em(iq, km-1)*em(iq, 1)
        end do
!        forall (km = 1:kmmax, iq = 1:nstatic) em(iq, -km) = conjg(em(iq, km))
        do km = 1, kmmax
            do iq = 1, nwall
                em(iq, -km) = conjg(em(iq, km))
            end do
        end do
        do kn = 2, knmax
            forall (iq = 1:nwall) en(iq, kn) = en(iq, kn-1)*en(iq, 1)
        end do
!        forall (kn = 1:knmax, iq = 1:nstatic) en(iq, -kn) = conjg(en(iq, kn))
        do kn = 1, knmax
            do iq = 1, nwall
                en(iq, -kn) = conjg(en(iq, kn))
            end do
        end do
 
        ewqWW = 0.0
        kcount = 0
        kcount2 = 0
        kmmin = 0
        knmin = 1
        do kl = 0, klmax
           rl = pipix*real(kl)
            do km = kmmin, kmmax
                rm = pipiy*real(km)
                forall (iq = 1:nwall) elm(iq) = qs(nmove + iq)*el(iq, kl)*em(iq, km)
                do kn = knmin, knmax
                    rn = pipiz*real(kn)
                    kcount = kcount + 1
                    if (kcount > vkmax) then
                        write(*, *) kcount, vkmax
                        stop 'force_fouri_init(): too large kcount!'
                    end if
                    myk(kcount) = -1
                    ksq = real(kl*kl + km*km) + real(kn*kn*kmax*kmax)/real(knmax*knmax)
!                    ksq = rl*rl + rm*rm + rn*rn
!                    if(ksq < rkmaxsq) then
                    if(ksq < real(kmax*kmax) + 2.01) then
                        kcount2 = kcount2 + 1
                        myk(kcount) = mod(kcount2, nproc)     ! precomputation
!                        if (kl == klmax .or. km == kmmax .or. kn == knmax) print *,  rl, rm, rn
                        rksq = rl*rl + rm*rm + rn*rn
                        Ak(kcount) = fac*exp(rksq*alf4i)/rksq ! precomputation

                        Qsumw1(kcount) = sum(elm(n1max+1-nmove:nwall/2)*en(n1max+1-nmove:nwall/2, kn)) ! precomputation
                        Qsumw2(kcount) = sum(elm(n2max+1-nmove:nwall)*en(n2max+1-nmove:nwall, kn)) ! precomputation

                        Qws = Qsumw1(kcount) + Qsumw2(kcount)
                        ewqWW = ewqWW + Ak(kcount)*Qws*conjg(Qws)
                    end if
                end do
                knmin = -knmax
            end do
            kmmin = -kmmax
        end do
        print *, 'nn', n1max+1-nmove, nwall/2
        print *, 'nn2', n2max+1-nmove, nwall
!        print *,'stat', nstatic
!        print *, n1max+1-nmove, nwall/2
!        print *, n2max+1-nmove, nwall
!        print *,'ng', ngogo,n2max+1-nmove, nwall
!        print *,'nn', n1min, n1max, n2min, n2max
!        print *,'QQ', sum(Qsumw1), sum(Qsumw2)
!        print *,'Q', sum(Qsumw1 + Qsumw2)

        ewqWW = 0.5*ewqWW
        if (kcount > vkmax) then
            print *, 'vkmax', vkmax, kcount
            stop 'too many kcounts'
        end if

        ! correction for self-term
        uself = -0.5*spii*facq*alfa*sum(qs(1:ntota)**2) ! atom selfinteraction
        ix = 0
        do it = 1, ntype        ! moltypes
            select case (moltypen(it))
                case (0)             ! rigmol selfinteraction
                    i = nt(it-1)        ! take the first molecule
                    do is = i + 1, i + ns(it)
                        do js = is + 1, i + ns(it)
                            rij = sqrt((x0(is) - x0(js))**2 + (y0(is) - y0(js))**2 + (z0(is) - z0(js))**2)
                            t = 1.0/(1.0 + p*alfa*rij)
                            erfv = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-(alfa*rij)**2)
                            uself = uself - nt(it)*facq*qs(is)*qs(js)*erfv/rij
                        end do
                    end do
                case (2)             ! constrained atoms selfinteraction
                    do iq = 1, nt(it)
                        ix = ix + 1
                        do is = 1, nconst(it)
                            i = constra(ix, is)
                            j = constrb(ix, is)
                            rij = length(it, is)
                            t = 1.0/(1.0 + p*alfa*rij)
                            erfv = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-(alfa*rij)**2)
                            uself = uself - facq*qs(i)*qs(j)*erfv/rij
                        end do
                    end do
                case default
            end select
        end do

        return
    end subroutine force_fouri_init
    !
    subroutine force_fouri
        implicit none
        integer*4 :: iq, kl, km, kn, kmmin, knmin, it, is, i, kcount
        real*8 :: rl, rm, rn
        real*8 :: dipz, Akx, Akx1, ewsurf
        complex*16 :: QsumLL, Qsurf


        ! el(0), em(0), and en(0) computed at the beginning in init sub.
        forall (i = 1:ngogo)
            el(i, 1) = exp(cpipix*x0(i))
            em(i, 1) = exp(cpipiy*y0(i))
            en(i, 1) = exp(cpipiz*z0(i))
        end forall

        do kl = 2, klmax
            forall (iq = 1:ngogo) el(iq, kl) = el(iq, kl-1)*el(iq, 1)
        end do
        do km = 2, kmmax
            forall (iq = 1:ngogo) em(iq, km) = em(iq, km-1)*em(iq, 1)
        end do
!        forall (km = 1:kmmax, iq = 1:ngogo) em(iq, -km) = conjg(em(iq, km))
        do km = 1, kmmax
            do iq = 1, ngogo
                em(iq, -km) = conjg(em(iq, km))
            end do
        end do
        do kn = 2, knmax
            forall (iq = 1:ngogo) en(iq, kn) = en(iq, kn-1)*en(iq, 1)
        end do
!        forall (kn = 1:knmax, iq = 1:ngogo) en(iq, -kn) = conjg(en(iq, kn))
        do kn = 1, knmax
            do iq = 1, ngogo
                en(iq, -kn) = conjg(en(iq, kn))
            end do
        end do
 
        ewsurf = 0.0 ! temporrary
        ewqLL = 0.0
        uewa = 0.0
        kcount = 0
        kmmin = 0
        knmin = 1
        do kl = 0, klmax                                  ! main k-space loop
            rl = pipix*real(kl)
            do km = kmmin, kmmax
                rm = pipiy*real(km)
                elm(1:ngogo) = qs(1:ngogo)*el(1:ngogo, kl)*em(1:ngogo, km)
                do kn = knmin, knmax
                    kcount = kcount + 1
                    if (myk(kcount) == me) then            ! parallel
                        rn = pipiz*real(kn)
                        elmn(1:ngogo) = elm(1:ngogo)*en(1:ngogo, kn)
                        QsumLL = sum(elmn(1:nmove))

                        Qsum(kcount) = QsumLL + Qsumw1(kcount) + Qsumw2(kcount) + sum(elmn(nmove+1:ngogo))
                        Akx = Ak(kcount)
                        ewqLL = ewqLL + Akx*QsumLL*conjg(QsumLL)
                        uewa = uewa + Akx*dble(Qsum(kcount)*conjg(Qsum(kcount)))

                        ! temporarry

!                        Qsurf = Qsumw1(kcount) + Qsumw2(kcount) + sum(elmn(nliqa+1:nmove))
!                        if (nsurf == 0) then
!                            do i = 1, nriga, 3
!                                !if (firstbin(i) == 1) Qsurf = Qsurf + elmn(i) + elmn(i+1) + elmn(i+2)
!                                Qsurf = Qsurf + elmn(i) + elmn(i+1) + elmn(i+2)
!                            end do
!                        end if
!                        ewsurf = ewsurf + Akx*Qsurf*conjg(Qsurf)

                        ft(1:nmove) = Akx*dimag(conjg(Qsum(kcount))*elmn(1:nmove))
                        fx(1:nmove) = fx(1:nmove) + ft(1:nmove)*rl
                        fy(1:nmove) = fy(1:nmove) + ft(1:nmove)*rm
                        fz(1:nmove) = fz(1:nmove) + ft(1:nmove)*rn
!                        Akx1 = Akx*dimag(conjg(Qsum(kcount))*(sum(elmn(n1min:n1max) + Qsumw1(kcount) + sum(elmn(nliqa+1:nliqa+nsurf/2)))))
!                        fwsx(1) = fwsx(1) + Akx1*rl
!                        fwsy(1) = fwsy(1) + Akx1*rm
!                        fwsz(1) = fwsz(1) + Akx1*rn
!                        Akx  = Akx*dimag(conjg(Qsum(kcount))*(sum(elmn(n2min:n2max) + Qsumw2(kcount) + sum(elmn(nliqa+nsurf/2+1:nmove)))))
!                        fwsx(-1) = fwsx(-1) + Akx*rl
!                        fwsy(-1) = fwsy(-1) + Akx*rm
!                        fwsz(-1) = fwsz(-1) + Akx*rn
                    end if
                end do
                knmin = -knmax
            end do
            kmmin = -kmmax
        end do
        ewqLL = 0.5*ewqLL
        uewa = 0.5*uewa
        ewqWL = uewa - ewqWW - ewqLL

        ewsurf = 0.5*ewsurf
        uavex = uavex + uewa - ewsurf! - 0.5*ewsurf

        ! 2D corrections
        dipz = sum(qs(1:nmove)*z0(1:nmove))
        crf = -2.0*M_PI*lzi*lxi*lyi
        urf = -crf*facq*dipz**2         ! energy of the 2D correction
        crf = 2.0*crf*facq*dipz         ! 2D correction to z-forces
        fz(iwmin:iwmax) = fz(iwmin:iwmax) + crf*qs(iwmin:iwmax)

        do i = nriga+1, nliqa
            it = itype(i)
            kl = int((z0(i)+lwh)*rdeli + 0.5)
            efipos(it-2, kl) = efipos(it-2, kl) + fz(i)
        end do
 
        return
    end subroutine force_fouri

    subroutine force_extern_init
        implicit none
        integer*4 :: ie

        !!! BE CAREFUL force applied to every atom, not molecule!
        ie = nang/nproc
        iemin = ie*me + 1
        iemax = ie*(me + 1)
        if (me == nproc-1) iemax = nliqa

        fxext = fxext
        fyext = fyext
        fzext = 0.0

        return
    end subroutine force_extern_init

    subroutine force_extern
        implicit none
        integer*4 :: i

        forall (i = iemin:iemax)
            fx(i) = fx(i) + fxext
            fy(i) = fy(i) + fyext
!            fz(i) = fz(i) + fzext
        end forall

        return
    end subroutine force_extern

    subroutine force_sum
        implicit none
        integer*4 :: i, j, ix, bmax
        integer*4 :: ierr
        real*8 :: aux

        integer*4 :: it, k, l, ig, ib, ib2
        real*8 :: rij, aux2, aux3
        include 'mpif.h'

        ! pressure correction for fixed bonds

!        aux = uavex + ulr
!        print *, 'start0', uavex


        aux = uavex
        call mpi_allreduce(aux, uavex, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        uavex = uavex + urf
        uavextot = uavextot + uavex
!        print *, 'start'
        fwsx(0) = sum(fx(1:nmove))
        fwsy(0) = sum(fy(1:nmove))
        fwsz(0) = sum(fz(1:nmove))

!        print *, 'suma', sum(abs(fx(nriga+1:nmove))), sum(abs(fy(nriga+1:nmove))), sum(abs(fz(nriga+1:nmove)))
        ftx = 0.0
        call mpi_allreduce(fx(nriga+1:nmove), ftx(nriga+1:nmove), nmove-nriga, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        fty = 0.0
        call mpi_allreduce(fy(nriga+1:nmove), fty(nriga+1:nmove), nmove-nriga, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        ftz = 0.0
        call mpi_allreduce(fz(nriga+1:nmove), ftz(nriga+1:nmove), nmove-nriga, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!        print *, 'sumb', sum(abs(ftx(nriga+1:nmove))), sum(abs(fty(nriga+1:nmove))), sum(abs(ftz(nriga+1:nmove)))

        ! rigmol forces & torques !!! asumes ns(1) = 3
        do i = 1, nrig
            j = 3*(i-1) + 1
            rsx(1:3) = x0(j:j+2) - cx0(i)
            rsy(1:3) = y0(j:j+2) - cy0(i)
            rsz(1:3) = z0(j:j+2) - cz0(i)
            btx(i) = sum(rsy(1:3)*fz(j:j+2) - rsz(1:3)*fy(j:j+2))
            bty(i) = sum(rsz(1:3)*fx(j:j+2) - rsx(1:3)*fz(j:j+2))
            btz(i) = sum(rsx(1:3)*fy(j:j+2) - rsy(1:3)*fx(j:j+2))
        end do
        tx = 0.0
        call mpi_allreduce(btx(1:nrig), tx(1:nrig), nrig, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        ty = 0.0
        call mpi_allreduce(bty(1:nrig), ty(1:nrig), nrig, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        tz = 0.0
        call mpi_allreduce(btz(1:nrig), tz(1:nrig), nrig, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        do i = 1, nrig
            j = 3*(i-1) + 1
            fx(i) = sum(fx(j:j+2))
            fy(i) = sum(fy(j:j+2))
            fz(i) = sum(fz(j:j+2))
        end do
        fcx = 0.0
        call mpi_allreduce(fx(1:nrig), fcx(1:nrig), nrig, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        fcy = 0.0
        call mpi_allreduce(fy(1:nrig), fcy(1:nrig), nrig, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        fcz = 0.0
        call mpi_allreduce(fz(1:nrig), fcz(1:nrig), nrig, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

        ! other statistics
        bufs = 0.0
        bufs(1:nijt) = elj(1:nijt)  ! energy
        bufs(nijt+1:2*nijt) = eqq(1:nijt)
        bmax = 2*nijt + 1
        bufs(bmax) = uewa
        bufs(bmax+1:bmax+3) = fwsx(-1:1)  ! wall force
        bufs(bmax+4:bmax+6) = fwsy(-1:1)  ! wall force
        bufs(bmax+7:bmax+9) = fwsz(-1:1)  ! wall force
        bmax = bmax + 9
!        print *, 'a', me, elj(1:nijt)
!        print *, 'b', me, bufs(1:nijt)
        call mpi_allreduce(bufs(1:bmax), stat(1:bmax), bmax, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!        if(me==0) print *, 'sb',me, stat(1:nijt)
        elj(1:nijt) = stat(1:nijt)
        eqq(1:nijt) = stat(nijt+1:2*nijt)
        bmax = bmax - 9
        uewa = stat(bmax)
        fwsx(-1:1) = stat(bmax+1:bmax+3)
        fwsy(-1:1) = stat(bmax+4:bmax+6)
        fwsz(-1:1) = stat(bmax+7:bmax+9)

!        uconf = sum(stat(1:bmax)) + uself + urf + ulr
        uconf = sum(stat(1:bmax)) + uself + urf
        if (.false.) then
            ftx(nriga+1:nmove) = fx(nriga+1:nmove)
            fty(nriga+1:nmove) = fy(nriga+1:nmove)
            ftz(nriga+1:nmove) = fz(nriga+1:nmove)

            fcx(1:nrig) = fx(1:nrig)
            fcy(1:nrig) = fy(1:nrig)
            fcz(1:nrig) = fz(1:nrig)

            do i = 1, nrig
                j = 3*(i-1) + 1
                rsx(1:3) = x0(j:j+2) - cx0(i)
                rsy(1:3) = y0(j:j+2) - cy0(i)
                rsz(1:3) = z0(j:j+2) - cz0(i)
                tx(i) = sum(rsy(1:3)*fz(j:j+2) - rsz(1:3)*fy(j:j+2))
                ty(i) = sum(rsz(1:3)*fx(j:j+2) - rsx(1:3)*fz(j:j+2))
                tz(i) = sum(rsx(1:3)*fy(j:j+2) - rsy(1:3)*fx(j:j+2))
            end do
            uconf = sum(elj(1:nijt)) + sum(eqq(1:nijt)) + uewa + uself + urf
!        print *, 'sumd', sum(abs(ftx(nriga+1:nmove))), sum(abs(fty(nriga+1:nmove))), sum(abs(ftz(nriga+1:nmove)))
        end if
!        print *, 'sumc', sum(abs(ftx(nriga+1:nmove))), sum(abs(fty(nriga+1:nmove))), sum(abs(ftz(nriga+1:nmove)))

        uconf = uconf - uself
        !print *, 'u', uself, urf, ulr, stat(1:bmax)
        press = (fwsz(1) - fwsz(-1))/(2.0*lx*ly)*1000.0/(PH_Na*1e-30)*1e-5
!        print *, 'press', press, real(nrig)*tempnow/(lx*ly*lw)*1000.0/(PH_Na*1e-30)*1e-5

        !if (nsurf > 0) call force_surfpress

        totpres = totpres + press


        fwxmeas = fwsx
        fwymeas = fwsy
        fwzmeas = fwsz

        if (me == 0) then
            do i = nriga+1, nliqa
                it = itype(i)
                l = int((z0(i)+lwh)*rdeli + 0.5)
                nefipos(it-2, l) = nefipos(it-2, l) + 1
                tefipos(it-2, l) = tefipos(it-2, l) + ftz(i)
            end do
        end if

!        print *, 'end'
!        print *, fwsx(0), fwsy(0), fwsz(0)
!        print *, fwsx(1), fwsy(1), fwsz(1)
!        print *, fwsx(-1), fwsy(-1), fwsz(-1)

!        aux = (sum(eqq(2:nijt) + elj(2:nijt)) + ewqWL)
!        print *, sum(stat(1:bmax)), bmax
!        print *, 'xlj', elj(1)/nrig
!        print *, 'xqq', (eqq(1)+ewqLL)/nrig
!        print *, 'xtot', (elj(1) + eqq(1) + ewqLL + (aux+ulr)/2.0 + uself)/nrig
        !print *, 'uself', uself
        !print *, 'elj',  sum(elj(1:nijt))
        !print *, 'eqq', sum(eqq(1:nijt))
        !print *, 'ewa', uewa, ewqWW, ewqWL, ewqLL
        !print *, 'ulr and urf', ulr, urf
!        print *, 'ewald', ewqWL/nrig
!        print *, 'waltot', aux/nrig
!        print *, 'ulr', ulr/nrig, urf/nrig
!        print *, 'press', press
!        print *, 'ew', uewa/nrig, ewqLL/nrig, ewqWL/nrig, ewqWW/nrig

        return
    end subroutine force_sum

    subroutine force_surfpress
        implicit none
        integer*4 :: i, j, it, ix, iy, is, ii, jj, kk, sgn
        real*8 :: ax, ay, az, au, fij, ai, xx, yy, zz
        

        print *,'a', fwsx(1), fwsy(1), fwsz(1)
        print *,'a', fwsx(-1), fwsy(-1), fwsz(-1)

        ix = 0
        do it = 1, ntype
            if (nconst(it) > 0) then
            do iy = 1, nt(it)
                ix = ix + 1

                is = nconst(it)
                ii = constra(is, ix)
                jj = constrb(is, ix)

                ax = x0(ii) - x0(jj)          ! h - o
                ax = ax - lx*anint(ax*lxi)
                ay = y0(ii) - y0(jj)
                ay = ay - ly*anint(ay*lyi)
                az = z0(ii) - z0(jj)
                ai = 1.0/sqrt(ax*ax + ay*ay + az*az)
                ax = ax*ai
                ay = ay*ai
                az = az*ai
                fij = ftx(jj)*ax + fty(jj)*ay + ftz(jj)*az
                !print *, ftz(ii), ftz(jj)
                !fij = 0.0
                xx = ftx(ii) + fij*ax
                yy = fty(ii) + fij*ay
                zz = ftz(ii) + fij*az

                if (nconst(it) == 2) then
                    ii = constra(2, ix)
                    jj = constrb(2, ix)

                    ax = x0(ii) - x0(jj)          ! h - o
                    ax = ax - lx*anint(ax*lxi)
                    ay = y0(ii) - y0(jj)
                    ay = ay - ly*anint(ay*lyi)
                    az = z0(ii) - z0(jj)
                    ai = 1.0/sqrt(ax*ax + ay*ay + az*az)
                    ax = ax*ai
                    ay = ay*ai
                    az = az*ai
                    fij = xx*ax + yy*ay + zz*az
                    xx = xx + fij*ax
                    yy = yy + fij*ay
                    zz = zz + fij*az
                else
                    ii = constra(2, ix)
                    jj = constrb(2, ix)
                    kk = constra(1, ix)
                    if (constrb(1, ix) /= jj) print *, 'constr', jj, constrb(1, ix)

                    ax = x0(ii) - x0(jj)          ! h - o
                    ax = ax - lx*anint(ax*lxi)
                    ay = y0(ii) - y0(jj)
                    ay = ay - ly*anint(ay*lyi)
                    az = z0(ii) - z0(jj)
                    au = x0(kk) - x0(jj)          ! h - o
                    ax = ax + au - lx*anint(au*lxi)
                    au = y0(kk) - y0(jj)
                    ay = ay + au - ly*anint(au*lyi)
                    az = az + z0(kk) - z0(jj)

                    ai = 1.0/sqrt(ax*ax + ay*ay + az*az)
                    ax = ax*ai
                    ay = ay*ai
                    az = az*ai
                    !print *, 'aa', ax, ay, az
                    fij = xx*ax + yy*ay + zz*az
                    xx = xx + fij*ax
                    yy = yy + fij*ay
                    zz = zz + fij*az
                end if

                sgn = int(sign(1.d0, z0(ii)))      ! ???
                fwsx(sgn) = fwsx(sgn) + xx
                fwsy(sgn) = fwsy(sgn) + yy
                fwsz(sgn) = fwsz(sgn) + zz

                !if (sgn == -1) then
                !print *, sgn
                !print *, xx, yy, zz
                !print *, fwsx(sgn), fwsy(sgn), fwsz(sgn)
                !end if

            end do
            end if
        end do
        print *,'b', fwsx(1), fwsy(1), fwsz(1)
        print *,'b', fwsx(-1), fwsy(-1), fwsz(-1)

        return
    end subroutine force_surfpress

    subroutine chempot
        implicit none
        real*8 :: h1r(3), h2r(3), or(3), mr(3)
        real*8 :: alph, cosa, sina, aco, asi, dum, ri
        real*8 :: uchem, uchemx, beta, rinsmin, rwinsmin, dw
        integer*4 :: i, j, jt, ijt, ix, jx, jmin, ic, pins
        real*8 :: xio, yio, zio, rsq, rsqi, rij, r6, aux, t, alfar, exp2a, erfx
        real*8 :: xi, yi, zi, qx, xdist, ydist, zdist, oxr, oyr, ozr

        integer*4 :: iq, kl, km, kn, kmmin, knmin, it, is, kcount, il, ir, k, kkmax, irot
        real*8 :: dipz, crf, Akx, Akx1, ewsurf
        real*8, dimension(3) :: qsw

        integer*8 :: ierr
        include 'mpif.h'

        aco = -0.333333
        asi = 0.942809
        rinsmin = 2.0**2
        rwinsmin = 2.35**2
        irot = 20
        qsw(1) = -0.8476
        qsw(2) =  0.4238
        qsw(3) =  0.4238

        beta = 1000.0/(PH_R*chtemper)

!        print *,'start', inave(3)
        do ic = 1, ichmax

            ! position
            oxr = (rnd(dum) - 0.5)*lx
            oyr = (rnd(dum) - 0.5)*ly
            do
                ozr = ((rnd(dum) - 0.5)*lwmax*2.0)
                if (abs(ozr) > lwmin) exit
            end do
            rins(1, 1) = oxr
            rins(1, 2) = oyr
            rins(1, 3) = ozr
            ! end position

            il = int((abs(ozr) - lwmin)*facbin) + 1

            if (il > ilmid) then
                irot = 20
            else
                irot = 1
            end if
            inave(il) = inave(il) + irot


            ! oxygen interactions

            i = 1
            uchem = 0.0
            pins = 1

            k = 0
            do j = 1, ntota
                if (itype(j) == 2) cycle
                zio = ozr - z0(j)
                if (abs(zio) > rcut) cycle
                xio = oxr - x0(j)
                xio = xio - lx*anint(xio*lxi)
                yio = oyr - y0(j)
                yio = yio - ly*anint(yio*lyi)
                rsq = xio*xio + yio*yio + zio*zio
                if (rsq > rcutsq) cycle

                it = 1
                jt = itype(j)
                ijt = ijtype(it, jt)

                if (rsq < rwinsmin) then
                    !print *, 'overlap', j, jt, sqrt(rsq), ozr
                    if (jt == 1) then
                        pins = 0
                        exit
                    else if (rsq < rinsmin) then
                        pins = 0
                        exit
                    end if
                end if
  
                rsqi = 1.0/rsq
                rij = sqrt(rsq)
  
                select case (vdwnum(ijt))
                    case (0)
                    case (1)                          ! lennard-jones
                        r6 = sig(ijt)*rsqi**3
                        uchem = uchem + eps(ijt)*r6*(r6 - 1.0)
                    case (2)                          ! Buckingham
                        r6 = cc(ijt)*rsqi**3
                        aux = aa(ijt)*exp(-rij*rhoi(ijt))
                        uchem = uchem + aux - r6
                    case default
                        write(*, *) i, j, ijt, vdwnum(ijt), vdwtyp(ijt)
                        stop 'force_real(): unknown interaction type (f-f)!'
                end select
 
             
                aux = facq*qsw(i)*qs(j)/rij
                alfar = alfa*rij
                t = 1.0/(1.0 + p*alfar)
                exp2a = exp(-alfar*alfar)
                erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp2a
                uchem = uchem + erfx*aux

                k = k + 1
                iarr(k) = j
                xarr(k) = xio
                yarr(k) = yio
                zarr(k) = zio
            end do
            if (pins == 0) exit
            kkmax = k

            uchemx = uchem

            do ir = 1, irot
            ! orientation
                do
                    h1r(1) = rnd(dum) - 0.5
                    h1r(2) = rnd(dum) - 0.5
                    h1r(3) = rnd(dum) - 0.5
                    ri = h1r(1)*h1r(1) + h1r(2)*h1r(2) + h1r(3)*h1r(3)
                    if (ri < 0.25 .and. ri > 0.001) exit
                end do
  
                ri = 1.0/sqrt(ri)
                h1r = h1r*ri
  
                ri = 1.0/sqrt(h1r(1)*h1r(1) + h1r(2)*h1r(2))
                h2r(1) = -h1r(2)*ri
                h2r(2) = h1r(1)*ri
                h2r(3) = 0.0
  
                mr(1) = (h2r(2)*h1r(3) - h2r(3)*h1r(2))
                mr(2) = (h2r(3)*h1r(1) - h2r(1)*h1r(3))
                mr(3) = (h2r(1)*h1r(2) - h2r(2)*h1r(1))
  
                alph = 2.0*M_PI*rnd(dum)
                cosa = cos(alph)
                sina = sin(alph)
  
                rins(2,:) = rins(1, :) + aco*h1r + asi*(cosa*h2r + sina*mr)
                rins(3,:) = rins(1, :) + h1r
            ! end orientation

                uchem = 0.0

                do k = 1, kkmax 
                    j = iarr(k)
                    xio = xarr(k)
                    yio = yarr(k)
                    zio = zarr(k)
            
                    jt = itype(j)
                    if (jt == 1) then
                        jmin = j + 1             ! don't cycle over the first pair already computed
                        do ix = 1, 3           ! coulombic forces and energy (and virial)
                            xi = xio + rins(ix, 1) - oxr
                            yi = yio + rins(ix, 2) - oyr
                            zi = zio + rins(ix, 3) - ozr
                            qx = facq*qsw(ix)
                            do jx = jmin, j+nsi(jt)-1           ! coulombic forces and energy (and virial)
                                xdist = xi - rsx(jx)
                                ydist = yi - rsy(jx)
                                zdist = zi - rsz(jx)
                                rsq = xdist*xdist + ydist*ydist + zdist*zdist
                                rsqi = 1.0/rsq
                                rij = sqrt(rsq)
                                aux = qx*qs(jx)/rij
                                alfar = alfa*rij
                                t = 1.0/(1.0 + p*alfar)
                                exp2a = exp(-alfar*alfar)
                                erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp2a
                                uchem = uchem + erfx*aux
                            end do
                            jmin = j
                        end do
                    else
                        do ix = 2, 3           ! coulombic forces and energy (and virial)
                            xi = xio + rins(ix, 1) - oxr
                            yi = yio + rins(ix, 2) - oyr
                            zi = zio + rins(ix, 3) - ozr
                            qx = facq*qsw(ix)
                            rsq = xi*xi + yi*yi + zi*zi
                            rsqi = 1.0/rsq
                            rij = sqrt(rsq)
                            aux = qx*qs(j)/rij
                            alfar = alfa*rij
                            t = 1.0/(1.0 + p*alfar)
                            exp2a = exp(-alfar*alfar)
                            erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp2a
                            uchem = uchem + erfx*aux
                        end do
                    end if

                end do

                uchem = uchem + uchemx
            ! End real part

            ! Ewald contribution

                forall (iq = 1:3)
                    el(iq, 1) = exp(cpipix*rins(iq, 1))
                    em(iq, 1) = exp(cpipiy*rins(iq, 2))
                    en(iq, 1) = exp(cpipiz*rins(iq, 3))
                end forall
                
                do kl = 2, klmax
                    forall (iq = 1:3) el(iq, kl) = el(iq, kl-1)*el(iq, 1)
                end do
                do km = 2, kmmax
                    forall (iq = 1:3) em(iq, km) = em(iq, km-1)*em(iq, 1)
                end do
!                forall (km = 1:kmmax, iq = 1:3) em(iq, -km) = conjg(em(iq, km))
                do km = 1, kmmax
                    do iq = 1, 3
                        em(iq, -km) = conjg(em(iq, km))
                    end do
                end do
                do kn = 2, knmax
                    forall (iq = 1:3) en(iq, kn) = en(iq, kn-1)*en(iq, 1)
                end do
!                forall (kn = 1:knmax, iq = 1:3) en(iq, -kn) = conjg(en(iq, kn))
                do kn = 1, knmax
                    do iq = 1, 3
                        en(iq, -kn) = conjg(en(iq, kn))
                    end do
                end do
             
                kcount = 0
                kmmin = 0
                knmin = 1
                do kl = 0, klmax                                  ! main k-space loop
                    do km = kmmin, kmmax
                        elm(1:3) = qsw(1:3)*el(1:3, kl)*em(1:3, km)
                        do kn = knmin, knmax
                            kcount = kcount + 1
                            uchem = uchem + 0.5*Ak(kcount)*dble(sum(elm(1:3)*en(1:3, kn))*conjg(Qsum(kcount)))
                        end do
                        knmin = -knmax
                    end do
                    kmmin = -kmmax
                end do
                
                uchem = uchem - crf*sum(qsw(1:3)*rins(1:3, 3))

                !end Ewald contribution
             
                do j = -1, 1, 2
                    dw = lwh - dble(j)*rins(1, 3)
                    if (dw < rcut) then
                        uchem = uchem + ulr_w(1)
                    else
                        uchem = uchem + ulr_l(1)
                        uchem = uchem + (flr_lja(1)/(9.0*dw**9) - flr_ljc(1)/(3.0*dw**3))
                        t = flr_bur(1)
                        exp2a = flr_bua(1)*exp(-dw/t)
                        uchem = uchem + (t*(2.0*t + dw)*exp2a - flr_buc(1)/(3.0*dw**3))
                    end if
                end do


                pave(il) = pave(il) + exp(-beta*uchem)      ! kJ/mol
             
            end do
            if (pins == 0) cycle
        end do

        pavesum = 0.0
        call mpi_allreduce(pave(1:lchmax), pavesum(1:lchmax), lchmax, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        inavesum = 0
        call mpi_allreduce(inave(1:lchmax), inavesum(1:lchmax), lchmax, mpi_integer8, mpi_sum, mpi_comm_world, ierr)

        return
    end subroutine chempot

end module force
