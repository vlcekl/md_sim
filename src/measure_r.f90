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
    real*8 :: rdel, rdeli, rfrst, rscnd
    integer*8 :: lmax

    logical :: l2d, ldiff, ldiff_r, lvisc, lhb

    ! correaltions functions
    integer*8, dimension(:, :), allocatable :: ngc
    integer*8, dimension(:, :, :), allocatable :: ng2d1, ng2d2
    integer*8, dimension(:, :), allocatable :: ngz, ngang
    real*8, dimension(:), allocatable :: ndipx, ndipy, ndipz, ndipsq, mux, muy, muz
    integer*8, dimension(:), allocatable :: nijat
    integer*4 :: igmin, igmax, lgxmax, lgymax, langmax
    real*8 :: dipnorm

    ! diffusion
    integer*8, dimension(:), allocatable :: iib
    real*8, dimension(:), allocatable :: diffx, diffy, diffz, xo, yo, zo
    real*8, dimension(:, :), allocatable :: dbinx, dbiny, dbinz, nbinz
    ! diffusion - rotation
    integer*8, dimension(:), allocatable :: sln
    real*8, dimension(:), allocatable :: xoo, yoo, zoo, sl1, sl2, sl3
    real*8, dimension(:, :), allocatable :: diff1, diff2, diff3
    integer*4 :: itima

    integer*8 :: rbin
    real*8 :: diftime, diftimax, timdif

    integer*4 :: vxmax
    real*8 :: dvxi
    real*8, dimension(-1:1) :: fwxmeas, fwymeas, fwzmeas, ftotx, ftoty, ftotz
    real*8, dimension(:, :), allocatable :: vxh, ivxh ! velocity profile

    real*8, dimension(:, :, :), allocatable :: rbonds
    integer*8, dimension(:, :, :), allocatable :: nbonds
    integer*8, dimension(:, :), allocatable :: laytyp, hblist
    integer*8, dimension(:), allocatable :: ibonds, iil

    real*8 :: rishell 
    integer*8, dimension(:), allocatable :: ibclust, wbclust, numib
    real*8, dimension(:, :), allocatable :: rbclust
    real*8, dimension(:), allocatable :: averij

    ! chemical potential
    real*8, dimension(3, 3) :: rins
    real*8 :: lwmin, lwmax, facbin, chtemper
    integer*8 :: ichmax, lchmax, ilmid
    real*8, dimension(:), allocatable :: pave, pavesum, xarr, yarr, zarr
    integer*8, dimension(:), allocatable :: inave, inavesum, iarr

    contains

    subroutine measure_init
        implicit none
        integer*4 :: ig


        itprn = nint(timprn/dt)  ! interval for measurement

        rdeli = 1.0/rdel

        ig = nmove/nproc
        igmin = ig*me + 1
        igmax = ig*(me + 1)
        if (me == nproc-1) igmax = nmove

        ! z-profile
        lmax = int(lw*rdeli + 0.5) + 1
        if (.true.) then
            !allocate(ngc(1+nstype-nltype, 0:lmax))
            allocate(ngc(nijt, 0:lmax))
        end if
        allocate(ngz(nit-2, 0:lmax), ndipx(0:lmax), ndipy(0:lmax), ndipz(0:lmax))
        allocate(ndipsq(0:lmax), mux(0:lmax), muy(0:lmax), muz(0:lmax))
        allocate(nijat(nijt))
        ! distance between oxygen and average of hydrogens (for dip. mom.)
        dipnorm = (xbody(1) - 0.5*(xbody(2) + xbody(3)))**2
        dipnorm = dipnorm + (ybody(1) - 0.5*(ybody(2) + ybody(3)))**2
        dipnorm = dipnorm + (zbody(1) - 0.5*(zbody(2) + zbody(3)))**2
        dipnorm = 1.0/(sqrt(dipnorm))

        ! flexi angles
        if (nsurf > 0) then
            langmax = 201
            allocate(ngang(nit-2, 0:langmax))
        end if

        ! diffusion
        if (ldiff) then
            allocate(xo(igmin:igmax), yo(igmin:igmax), zo(igmin:igmax))
            allocate(diffx(igmin:igmax), diffy(igmin:igmax), diffz(igmin:igmax))
            diffx = 0.0 ; diffy = 0.0 ; diffz = 0.0

            allocate(iib(nmove))

            allocate(dbinx(rbin, nit-2), dbiny(rbin, nit-2), dbinz(rbin, nit-2))
            dbinx = 0.0 ; dbiny = 0.0 ; dbinz = 0.0
            allocate(nbinz(rbin, nit-2))
            nbinz = 0.0

            ldiff_r = .false.
            ldiff_r = .true.
            if (ldiff_r) then
                itima = int(timdif/dt)/10 + 1
                allocate(diff1(rbin, itima), diff2(rbin, itima), diff3(rbin, itima))
                allocate(sl1(rbin), sl2(rbin), sl3(rbin), sln(rbin))
                allocate(xoo(igmin:igmax), yoo(igmin:igmax), zoo(igmin:igmax))
            end if
        end if

        ! 2d-profiles
        if (l2d) then
            lgxmax = 4.0*cx*rdeli + 1
            lgymax = 2.0*cy*rdeli + 1
            allocate(ng2d1(nit-2, 0:lgxmax, 0:lgymax), ng2d2(nit-2, 0:lgxmax, 0:lgymax))
        end if

        if (lhb) then
            allocate(laytyp(16, 2), hblist(1:6, nmol), nbonds(1:6, 1:6, 1:2), rbonds(1:6, 1:6, 1:2))
            allocate(ibonds(1:6), iil(1:6))
        end if

        if (lhb) then
            allocate(ibclust(0:8), wbclust(0:8), numib(6), averij(6), rbclust(0:8, 4))
        end if 

        ! viscosity
        lvisc = .false.
        if (lcouett .or. lfext) lvisc = .true.
        if (lvisc) then
            dvxi = 1.0/(lw/40.0)
            vxmax = int(lw*dvxi)
            allocate(vxh(nit-2, 0:vxmax), ivxh(nit-2, 0:vxmax))
        end if

        if (lchempot) then
            call rnd_init(14135621769d0+real(2*me*me + 7*me + 3))
            facbin = lwmin
            lwmin = lwh - lwmax  ! CAREFUL, this shwitches min and max 
            lwmax = lwh - facbin  ! from wall distance to center distance
            facbin = 1.0/(lwmax-lwmin)*float(lchmax)
            if (lchempot) then
                allocate(iarr(ntota), xarr(ntota), yarr(ntota), zarr(ntota))
                allocate(pave(lchmax), pavesum(lchmax))
                allocate(inave(lchmax), inavesum(lchmax))
            end if
            pave = 0.0
            inave = 0
            if (temper > 0.0) then
                chtemper = temper
            else
                chtemper = 300.0
            end if
        end if


        call measure_zero

        return
    end subroutine measure_init

    subroutine measure_zero
        implicit none
        integer*4 :: ig, ib, i

        ftotx = 0.0
        ftoty = 0.0
        ftotz = 0.0

        ngc = 0
        ngz = 0
        ndipx = 0.0
        ndipy = 0.0
        ndipz = 0.0
        ndipsq = 0.0
        if (l2d) then
            ng2d1 = 0
            ng2d2 = 0
        end if
        if (nsurf > 0) ngang = 0

        if (ldiff) then
            xo(igmin:igmax) = x0(igmin:igmax)
            yo(igmin:igmax) = y0(igmin:igmax)
            zo(igmin:igmax) = z0(igmin:igmax)
            diffx = 0.0 ; diffy = 0.0 ; diffz = 0.0
            dbinx = 0.0 ; dbiny = 0.0 ; dbinz = 0.0
            nbinz = 0.0
 
            iib = 0 
            if (rbin /= 3) then
                do ig = igmin, igmax
                    iib(ig) = int((z0(ig) + lwh)/(lw/real(rbin))) + 1
                end do
            else
                do ig = igmin, igmax
                    if (lwh - abs(z0(ig)) < rfrst) then
                        iib(ig) = 1
                    else if (lwh - abs(z0(ig)) < rscnd) then
                        iib(ig) = 2
                    else
                        iib(ig) = 3
                    end if
                end do
            end if
 
            diftime = 0.0

            if (ldiff_r) then
                do i = igmin, igmax
                    if (itype(i) == 2) then
                        ib = (i - 1)/3 * 3 + 1
                        print *, 'rot', i, ib
                        xoo(i) = x0(i) - x0(ib)
                        yoo(i) = y0(i) - y0(ib)
                        zoo(i) = z0(i) - z0(ib)
                    end if
                end do
                diff1 = 0.0 ; diff2 = 0.0 ; diff3 = 0.0
            end if
        end if


        if (lvisc) then
            vxh = 0.0
            ivxh = 0.0
        end if

        if (lhb) then
            laytyp = 0
            hblist = 0
            nbonds = 0
            rbonds = 0.0
            ibonds = 0
            iil = 0
        end if
        if (lhb) then
            ibclust = 0
            wbclust = 0
            rbclust = 0.0
        end if

        if (lchempot) then
            pave = 0.0
            inave = 0
        end if

        return
    end subroutine measure_zero


    subroutine measure_do
        implicit none
        integer*4 :: i, it, k, l, ig, ib, ib2, itcount
        real*8 :: rij, dipmz, dfac
        integer*4, save :: iti = 0
        integer*8, save :: nhbold

        
        ftotx = ftotx + fwxmeas
        ftoty = ftoty + fwymeas
        ftotz = ftotz + fwzmeas

        ! liquid - solid (g(r), orientation, charge distribution, ...)
        if (.true.) then
        mux = 0.0
        muy = 0.0
        muz = 0.0
        do i = igmin, igmax
            it = itype(i)
            rij = z0(i) + lwh
            l = int(rij*rdeli + 0.5)
            ngz(it, l) = ngz(it, l) + 1

            if (it == 1) then
                ! cos(theta) - angle of water dipole towards the surface normal
                muz(l) = muz(l) + (z0(i) - 0.5*(z0(i+1) + z0(i+2)))
                mux(l) = mux(l) + (x0(i) - 0.5*(x0(i+1) + x0(i+2)))
                muy(l) = muy(l) + (y0(i) - 0.5*(y0(i+1) + y0(i+2)))
            end if
 
            rij = lwh - abs(z0(i))
            if (l2d .and. rij < rscnd) then
                l = int(0.5*modulo(x0(i)+lxh + 4.0*cx, 4.0*cx)*rdeli + 0.5)
                k = int(0.5*modulo(y0(i)+lyh + 2.0*cy, 2.0*cy)*rdeli + 0.5)
                if (l < 0 .or. k < 0 .or. l > lgxmax .or. k > lgymax) print *, 'big', l, k
                if (rij < rfrst) then
                    ng2d1(it, l, k) = ng2d1(it, l, k) + 1
                else
                    ng2d2(it, l, k) = ng2d2(it, l, k) + 1
                end if
            end if
        end do
        ndipx = ndipx + mux
        ndipy = ndipy + muy
        ndipz = ndipz + muz
        ndipsq = ndipsq + mux**2 + muy**2 + muz**2
        end if

        if (lhb) then
            call measure_hb
!            print *, 'nhb', real(sum(nbonds)-nhbold)
            nhbold = sum(nbonds)
!            call measure_ib
        end if
        
        ! diffusion
        if (ldiff) then

            do i = igmin, igmax
                rij = x0(i) - xo(i)
                diffx(i) = diffx(i) + rij
                xo(i) = x0(i)
                rij = y0(i) - yo(i)
                diffy(i) = diffy(i) + rij
                yo(i) = y0(i)
                rij = z0(i) - zo(i)
                diffz(i) = diffz(i) + rij
                zo(i) = z0(i)
            end do

            if (ldiff_r) then
                dfac = 1.0
                iti = iti + 1
                if (mod(iti, 10) == 0) then
                    itcount = iti/10
                    sl1 = 0.0 ; sl2 = 0.0 ; sl3 = 0.0 ; sln = 0
                    do i = igmin, igmax
 
                        it = itype(i)
                        if (it /= 2) cycle
                        k = (i - 1)/3 * 3 + 1
                        ib = iib(k)
 
                        rij =       (x0(i) - x0(k))*xoo(i)
                        rij = rij + (y0(i) - y0(k))*yoo(i)
                        rij = rij + (z0(i) - z0(k))*zoo(i)
                        rij = rij*dfac
                        sl1(ib) = sl1(ib) + rij
                        sl2(ib) = sl2(ib) + 0.5*(3.0*rij*rij - 1.0)
                        sl3(ib) = sl3(ib) + 0.5*(5.0*rij*rij*rij - 3.0*rij)
                        sln(ib) = sln(ib) + 1
 
                    end do
                    do ib = 1, rbin
                        if (sln(ib) == 0) then
                            rij = 0.0
                        else
                            rij = 1.0/real(sln(ib))
                        end if
                        diff1(ib, itcount) = diff1(ib, itcount) + sl1(ib)*rij
                        diff2(ib, itcount) = diff2(ib, itcount) + sl2(ib)*rij
                        diff3(ib, itcount) = diff3(ib, itcount) + sl3(ib)*rij
                    end do
                end if
            end if

            if (diftime >= timdif) then

                do i = igmin, igmax
                    ib = iib(i)
                    if (rbin /= 3) then
                        ib2 = int((z0(i) + lwh)/(lw/real(rbin))) + 1
                    else
                            if (lwh - abs(z0(i)) < rfrst) then
                                ib2 = 1
                            else if (lwh - abs(z0(i)) < rscnd) then
                                ib2 = 2
                            else
                                ib2 = 3
                            end if
                    end if
 
                    it = itype(i)
 
                    if (ib == ib2) then
                        dbinx(ib, it) = dbinx(ib, it) + diffx(i)**2
                        dbiny(ib, it) = dbiny(ib, it) + diffy(i)**2
                        dbinz(ib, it) = dbinz(ib, it) + diffz(i)**2
                        nbinz(ib, it) = nbinz(ib, it) + 1.0
                    else
                        iib(i) = ib2
                        dbinx(ib, it) = dbinx(ib, it) + 0.5*diffx(i)**2
                        dbiny(ib, it) = dbiny(ib, it) + 0.5*diffy(i)**2
                        dbinz(ib, it) = dbinz(ib, it) + 0.5*diffz(i)**2
                        nbinz(ib, it) = nbinz(ib, it) + 0.5
      
                        dbinx(ib2, it) = dbinx(ib2, it) + 0.5*diffx(i)**2
                        dbiny(ib2, it) = dbiny(ib2, it) + 0.5*diffy(i)**2
                        dbinz(ib2, it) = dbinz(ib2, it) + 0.5*diffz(i)**2
                        nbinz(ib2, it) = nbinz(ib2, it) + 0.5
                    end if
                end do
 
                diftime = 0.0
                diffx = 0.0 ; diffy = 0.0 ; diffz = 0.0


                if (ldiff_r) then
                    do i = igmin, igmax
                        if (itype(i) == 2) then
                            k = (i - 1)/3 * 3 + 1
                            xoo(i) = (x0(i) - x0(k))
                            yoo(i) = (y0(i) - y0(k))
                            zoo(i) = (z0(i) - z0(k))
                        end if
                    end do
                    iti = 0
                end if

            else
                diftime = diftime + dt
            end if
        end if


        ! velocity profile
        if (lvisc) then
            do i = igmin, igmax
                it = itype(i)
                select case (it)
                    case(1)
                        k = (i-1)/3 + 1
                        l = int((cz0(k) + lwh)*dvxi)
                        vxh(it, l) = vxh(it, l) + cx1(k)
                        ivxh(it, l) = ivxh(it, l) + 1.0
                    case(3:)
                        l = int((z0(i) + lwh)*dvxi)
                        vxh(it, l) = vxh(it, l) + x1(i)
                        ivxh(it, l) = ivxh(it, l) + 1.0
                end select
            end do
        end if

        return
    end subroutine measure_do

    subroutine measure_output
        implicit none
        integer*4 :: l, it, jt, k, ib, ist, jst
        real*8 :: r, fac, facv, rx, ry, facv1, facv2, rh2o, r1, r2, x1, x2, p1, p2

        ! gz
!        print *, 'a'
        if (.true.) then
            call measure_gather_int2d(ngc)
        end if
        if (.true.) then
            call measure_gather_int2d(ngz)
            call measure_gather_dble1d(ndipx)
            call measure_gather_dble1d(ndipy)
            call measure_gather_dble1d(ndipz)
            call measure_gather_dble1d(ndipsq)
        end if

        if (lhb) then
            call measure_gather_int3d(nbonds)
            call measure_gather_dble3d(rbonds)

            call measure_gather_int1d(ibclust)
            call measure_gather_int1d(wbclust)
            call measure_gather_dble2d(rbclust)
        end if

        if (l2d) then
            call measure_gather_int3d(ng2d1)
            call measure_gather_int3d(ng2d2)
        end if

        if (nsurf > 0) then
            call measure_gather_int2d(ngang)
        end if

        if (ldiff) then
            forall (ib = 1:rbin, it = 1:nit-2, nbinz(ib, it) /= 0.0)
                dbinx(ib, it) = dbinx(ib, it)/(2.0*timdif*nbinz(ib, it))
                dbiny(ib, it) = dbiny(ib, it)/(2.0*timdif*nbinz(ib, it))
                dbinz(ib, it) = dbinz(ib, it)/(2.0*timdif*nbinz(ib, it))
            end forall
            call measure_gather_dble2d(dbinx)
            call measure_gather_dble2d(dbiny)
            call measure_gather_dble2d(dbinz)
            if (ldiff_r) then
                call measure_gather_dble2d(diff1)
                call measure_gather_dble2d(diff2)
                call measure_gather_dble2d(diff3)
            end if
            diff1 = diff1/real(nproc*int((timrun/timdif)))
            diff2 = diff2/real(nproc*int((timrun/timdif)))
            diff3 = diff3/real(nproc*int((timrun/timdif)))
        end if

!        print *, 'f'
        if (lvisc) then
            call measure_gather_dble2d(vxh)
            call measure_gather_dble2d(ivxh)
        end if

    if (me == 0) then
        write(*, *) '# wl1', ftotx(1)/real(itrun), ftoty(1)/real(itrun), ftotz(1)/real(itrun)
        write(*, *) '# wl2', ftotx(-1)/real(itrun), ftoty(-1)/real(itrun), ftotz(-1)/real(itrun)
        write(*, *) '# liq', ftotx(0)/real(itrun), ftoty(0)/real(itrun), ftotz(0)/real(itrun)
    end if
        if (me == 0) then
            open(4, file = filmeas, status='unknown')

            write(4, *) '# Total Forces'
            write(4, *) '# wl1', ftotx(1)/real(itrun), ftoty(1)/real(itrun), ftotz(1)/real(itrun)
            write(4, *) '# wl2', ftotx(-1)/real(itrun), ftoty(-1)/real(itrun), ftotz(-1)/real(itrun)
            write(4, *) '# liq', ftotx(0)/real(itrun), ftoty(0)/real(itrun), ftotz(0)/real(itrun)

         
            if (.true.) then
            fac = 3*lx*ly*lw/(2*M_PI*real(itrun)) !*real(nmol*nmol))
            lmax = idint(max(lxh, lyh, lwh)*rdeli) + 1

            k = 0
            do it = 1, ntype
                do ist = 1, nst(it)
                    do jt = it, ntype
                        do jst = 1, nst(jt)
                            k = k + 1
                            nijat(k) = nt(it)*nss(it, ist)*nt(jt)*nss(jt, jst)
                            if(it == ntype) nijat(k) = nijat(k)/4*ns(it)/3
                            if(jt == ntype) nijat(k) = nijat(k)/4*ns(jt)/3
                            print *, k, fac, fac/real(nijat(k)), nijat(k)
                        end do
                    end do
                end do
            end do

            if (sum(nijat) /= ntota) print *, 'wrong sum of atoms', sum(nijat), ntota

            do l = 1, lmax
                r = (real(l)+0.5)*rdel
                facv = fac/((r + rdel)**3 - r**3)
                write(4, 10) r, (facv*real(ngc(it, l))/real(nijat(it)), it = 1, nijt)
            end do

            write(4, *) '#'

            do l = 1, nijt
                k = int(2.375/rdel) + 1
                write(4, *) l, real(sum(ngc(l, 1:k)))/real(itrun)
            end do
!            do l = 1, lmax
!                r = (real(l)+0.5)*rdel
!                if (r < 3.000) then
!                    write(4, 10) r, (real(ngc(it, l))/real(itrun), it = 1, nijt)
!                end if
!            end do


            rh2o = 0.0333606866885  ! 1.0 g/cm^3 density of pure water
            rh2o = 0.033328051141044067 ! 0.997 g/cm^3
            facv = 1.0/((real(itrun)*rdel*lx*ly*rh2o))
            lmax = idint(lw*rdeli) + 1
            where (ngz(1,:) > 0) ndipsq = qs(1)**2*(ndipsq - (ndipx**2 + ndipy**2 + ndipz**2))/real(ngz(1,:))**2
            where (ngz(1,:) > 0) ndipz(:) = ndipz(:)*dipnorm/real(ngz(1,:))
            write(4, *) '# 1d'
            do l = 1, lmax
                r = (real(l)+0.5)*rdel
                write(4, 30) r - lwh, ndipz(l), ndipsq(l), (facv*ngz(it, l), it = 1, nit-2)
            end do
            end if

            if (l2d) then
                write(4, *) '# 2d'
                do it = 1, nit-2
                    facv1 = 10000.0/real(sum(ng2d1(it,:,:)))
                    facv2 = 10000.0/real(sum(ng2d2(it,:,:)))
                    do l = 1, int(4.0*cx*rdeli + 0.5)
                        do k = 1, int(2.0*cy*rdeli + 0.5)
                            rx = (real(l)+0.5)*rdel*2.0
                            ry = (real(k)+0.5)*rdel*2.0
                            write(4, 50) rx, ry, facv1*real(ng2d1(it, l, k)), facv2*real(ng2d2(it, l, k))
                        end do
                        write(4, *)
                    end do
                    write(4, *)
                end do
            end if

            if (nsurf > 0) then
                write(4, *) '# angle'
                fac = 1.0/real(itrun)
                do l = 1, 200
                    r = real(l-1)/100 - 1.0
                    write(4, 60) r, (fac*ngang(it, l)/real(nt(it)), it = 2, ntype-1)
                end do
            end if

            if (lhb) then
                write(4, *) '# Hydrogen Bonds'
                write(4, *) real(sum(nbonds))/real(itrun)
                write(4, *) '#  layers'
                write(4, 80) (ibonds(ib)/(itrun/2), ib = 1, 6)
                write(4, *) '# l1 l2  b1h b1o b2o b2h' 
                do ib = 1, 16
                    p1 = real(nbonds(laytyp(ib, 1), laytyp(ib, 2), 1))
                    p2 = real(nbonds(laytyp(ib, 1), laytyp(ib, 2), 2))
                    if (p1 + p2 > 0.0 .and. sum(ibonds(laytyp(ib, 1:2))) > 0) then
                        r1 = rbonds(laytyp(ib, 1), laytyp(ib, 2), 1)
                        r2 = rbonds(laytyp(ib, 1), laytyp(ib, 2), 2)
                        x1 = real(ibonds(laytyp(ib, 1)))
                        x2 = real(ibonds(laytyp(ib, 2)))
                        write(4, 90) laytyp(ib, 1), laytyp(ib, 2), p1/x1, p2/x1, p1/x2, p2/x2, r1/p1, r2/p2
                    end if
                end do
            end if

            if (lhb) then
                write(4, *) '# Ion Clusers'
                do ib = 0, 8
                    if (ibclust(ib) /= 0) write(4, 80) ib, real(ibclust(ib))/real(itrun/2), real(wbclust(ib))/real(ibclust(ib))
                end do
                write(4, *) '# Bondlenghts'
                do ib = 0, 8
                    if (ibclust(ib) /= 0) then
                        rbclust(ib, 1:4) = rbclust(ib, 1:4)/real(ibclust(ib))
                        write(4, 20) ib, (rbclust(ib, l), l = 1, 4)
                    end if
                end do
            end if

            if (ldiff) then
                ! convert from A^2/dt to cm^2/s*10^-5
                fac = (1e-8)**2/(1e-15)*1e5  
                write(4, *) '# Diffusion in x'
                do ib = 1, rbin
                    r = (real(ib)-0.5)*(lw/real(rbin)) - lwh
                    write(4, 70) r, (fac*dbinx(ib, it), it = 1, nit-2)
                end do 
                write(4, *) '# Diffusion in y'
                do ib = 1, rbin
                    r = (real(ib)-0.5)*(lw/real(rbin)) - lwh
                    write(4, 70) r, (fac*dbiny(ib, it), it = 1, nit-2)
                end do 
                write(4, *) '# Diffusion in z'
                do ib = 1, rbin
                    r = (real(ib)-0.5)*(lw/real(rbin)) - lwh
                    write(4, 70) r, (fac*dbinz(ib, it), it = 1, nit-2)
                end do 
                if (ldiff_r) then
                    write(4, *) '# Rotational diffusion l = 1'
                    do l = 1, itima
                        write(4, 70) (diff1(ib, l), ib = 1, rbin)
                    end do 

                    write(4, *) '# Rotational diffusion l = 2'
                    do l = 1, itima
                        write(4, 70) (diff2(ib, l), ib = 1, rbin)
                    end do 

                    write(4, *) '# Rotational diffusion l = 3'
                    do l = 1, itima
                        write(4, 70) (diff3(ib, l), ib = 1, rbin)
                    end do 

                end if
            end if

            if (lvisc) then
                write(4, *) '# visco'
                do l = 0, vxmax-1
                    r = real(l)/dvxi
                    write(4, *) r + 0.5/dvxi, (vxh(it, l)/ivxh(it, l), it = 1, nit-2)
                end do

            end if

            if (lchempot) then
                write(4, *) '# total ch.p.: ', -log(sum(pavesum(1:lchmax))/real(sum(inavesum(1:lchmax))))*PH_R*chtemper/1000.0
!                do l = lchmax, 1, -1
!                    write(4, *) lwh - (lwmin + (real(l)-0.5)*(lwmax-lwmin)/real(lchmax)),  -log(pavesum(l)/real(inavesum(l)))*PH_R*chtemper/1000.0, pavesum(l)/real(inavesum(l))
!                end do
            end if

            close(4, status='keep')
        end if

        call measure_deallocate

        return
10      format(f7.3,5(1x,f11.6))
30      format(f7.3,11(1x,f11.6))
50      format(f7.3,1x,f7.3,2(1x,f11.6))
60      format(f5.2,6(1x,f11.6))
70      format(12(1x,f11.6))
80      format(i2,2(1x,f7.3))
90      format(2(i2,1x),6(f11.6,1x))
20      format(i2,3(1x,f11.6))
    end subroutine measure_output

    subroutine measure_hb
        implicit none
        real*8, save :: cspi4
        integer*4, save :: cont = 0
        integer*4 :: iif, iis, iit, ii4
        integer*4 :: i, j, il, jl, ih, ihmax, k, ii, jj
        real*8 :: rij, xio, yio, zio, hi, rb

        cont = cont + 1

        ! make lists
        if (cont == 1) then        ! list of wall atoms

            cspi4 = cos(30.0/180.0*M_PI)

            laytyp(1:16, 1) = ((/1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6/))
            laytyp(1:16, 2) = ((/2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 6/))

            iis = 0
            iit = 0
            do i = nmove+1, ntota
                rij = abs(z0(i)) - lwh
                if (rij < -1.0) then
                    iis = iis + 1
                    hblist(1, iis) = i
                end if
            end do

            do i = nliqa+1, nmove
                if (qs(i) > 0) cycle
                if (itype(i) == itype(nliqa + 1)) then
                    iit = iit + 1
                    hblist(2, iit) = i
                else
                    iis = iis + 1
                    hblist(1, iis) = i
                end if
            end do

            iil(1) = iis
            iil(2) = iit

        end if
        if (mod(cont, 2) /= 0) return
!        if (mod(cont, 5) /= 0) return

        iif = 0
        iis = 0
        iit = 0
        ii4 = 0
        hblist(3:6, :) = 0
        do i = 1, nriga
            if (itype(i) == 2) cycle
            rij = lwh - abs(z0(i))
            if (rij < rscnd + 3.5) then
                if (rij < rscnd) then
                    if (rij < rfrst) then
                        iif = iif + 1
                        hblist(3, iif) = i
                    else
                        iis = iis + 1
                        hblist(4, iis) = i
                    end if
                else
                    iit = iit + 1
                    hblist(5, iit) = i
                end if
            else
                ii4 = ii4 + 1
                hblist(6, ii4) = i
            end if
        end do
        iil(3) = iif
        iil(4) = iis
        iil(5) = iit
        iil(6) = ii4

        ibonds = ibonds + iil

        ! check bonds

        do k = 1, 16
            if (me == mod(k, nproc)) then
                il = laytyp(k, 1)
                jl = laytyp(k, 2)
                do ii = 1, iil(il)
                    i = hblist(il, ii)
                    do jj = 1, iil(jl)
                        j = hblist(jl, jj)
                        if (z0(i)*z0(j) < 0.0 .and. k /= 16) cycle
                        xio = x0(i) - x0(j)
                        xio = xio - lx*anint(xio*lxi)
                        yio = y0(i) - y0(j)
                        yio = yio - ly*anint(yio*lyi)
                        zio = z0(i) - z0(j)
                        rij= xio*xio + yio*yio + zio*zio
                        if (rij < 3.5*3.5 .and. rij > 1e-5) then  ! possible bond
                            rij = sqrt(rij)
                            ihmax = 1
                            if(itype(i) == 1) ihmax = 2
                            if(abs(qs(i+1)) > 0.7) ihmax = 0 ! no hydrogen follows
                            do ih = 1, ihmax
                                hi =      (x0(i+ih) - x0(i))*xio
                                hi = hi + (y0(i+ih) - y0(i))*yio
                                hi = hi + (z0(i+ih) - z0(i))*zio
                                if (-hi/rij >= cspi4) then
                                    nbonds(il, jl, 1) = nbonds(il, jl, 1) + 1
                                    hi = x0(i+ih) - x0(j)
                                    rb = (hi - lx*anint(hi*lxi))**2
                                    hi = y0(i+ih) - y0(j)
                                    rb = rb + (hi - ly*anint(hi*lyi))**2
                                    rb = rb + (z0(i+ih) - z0(j))**2
                                    rbonds(il, jl, 1) = rbonds(il, jl, 1) + sqrt(rb)
                                end if
                            end do
                            ihmax = 1
                            if(itype(j) == 1) ihmax = 2
                            do ih = 1, ihmax
                                hi =      (x0(j+ih) - x0(j))*xio
                                hi = hi + (y0(j+ih) - y0(j))*yio
                                hi = hi + (z0(j+ih) - z0(j))*zio
                                if (hi/rij >= cspi4) then
                                    nbonds(il, jl, 2) = nbonds(il, jl, 2) + 1
                                    hi = x0(j+ih) - x0(i)
                                    rb = (hi - lx*anint(hi*lxi))**2
                                    hi = y0(j+ih) - y0(i)
                                    rb = rb + (hi - ly*anint(hi*lyi))**2
                                    rb = rb + (z0(j+ih) - z0(i))**2
                                    rbonds(il, jl, 2) = rbonds(il, jl, 2) + sqrt(rb)
                                end if
                            end do
                        end if
                    end do
                end do
            end if
        end do

        return
    end subroutine measure_hb

    subroutine measure_ib
        implicit none
        integer*4, save :: cont = 0
        integer*4 :: i, j, jl, k, jj, ilmax, ilmin, ibn, wbn
        real*8 :: rij, xio, yio, zio, rii
        character*40 :: fil
        character*40, save :: fil2

        cont = cont + 1
        if (cont == 1) then
            write(fil, '(i3)') me
            fil2 = adjustl(fil)
            i = len_trim(fil2)
            fil2(i+1:i+4) = 'icfg'
            open(4, file = fil2, status = 'new')
            close(4, status = 'keep')
        end if
        
        if (mod(cont, 2) /= 0) return

        if (mod(cont, 100) == 0) then
            open(4, file = fil2, access = 'append', status = 'old')
            write(4, *) '#', cont 
        end if

! print *, me, 'ib'
        k = 0
        do i = nriga+1, nliqa
            if (qs(i) < 0.0) exit
            k = k + 1
            if (mod(k, nproc) /= me) cycle
            ilmin = 1
            ilmax = 6
            rii = lwh - abs(z0(i))
            if (rii < rscnd) ilmax = 6
            if (rii > rscnd) ilmin = 1
            numib = 0
            averij = 0.0
            do jl = ilmin, ilmax
                do jj = 1, iil(jl)
                    j = hblist(jl, jj)
                    if (z0(i)*z0(j) < 0.0 .and. jl /= 6) cycle
                    xio = x0(i) - x0(j)
                    xio = xio - lx*anint(xio*lxi)
                    yio = y0(i) - y0(j)
                    yio = yio - ly*anint(yio*lyi)
                    zio = z0(i) - z0(j)
                    rij= xio*xio + yio*yio + zio*zio
                    if (rij < rishell*rishell) then
                        averij(jl) = averij(jl) + sqrt(rij)
                        numib(jl) = numib(jl) + 1
                    end if
                end do
            end do
            ibn = numib(1)*3 + numib(2) + numib(3)
            if (ibn > 8) stop 'too many ion neighbors'
            wbn = sum(numib(4:6))
!print *, ibn, wbn
            ibclust(ibn) = ibclust(ibn) + 1
            wbclust(ibn) = wbclust(ibn) + wbn
            if(numib(1) /= 0) rbclust(ibn, 1) = rbclust(ibn, 1) + averij(1)/real(numib(1))
            if(numib(2) /= 0) rbclust(ibn, 2) = rbclust(ibn, 2) + averij(2)/real(numib(2))
            if(numib(3) /= 0) rbclust(ibn, 2) = rbclust(ibn, 2) + averij(3)/real(numib(3))
            if(wbn /= 0) rbclust(ibn, 3) = rbclust(ibn, 3) + sum(averij(4:6))/real(wbn)
            rbclust(ibn, 4) = rbclust(ibn, 4) + rii

            if (mod(cont, 100) == 0) then
                write(4, *) i, x0(i), y0(i), z0(i), ibn, wbn
            end if

        end do

        if (mod(cont, 100) == 0) then
            k = 0
            do i = nriga+1, nliqa
                if (qs(i) > 0.0) cycle
                k = k + 1
                if (mod(k, nproc) /= me) cycle
                write(4, *) i, x0(i), y0(i), z0(i)
            end do
            close(4, status='keep')
        end if

        return
    end subroutine measure_ib

    subroutine measure_idump(fil) 
        implicit none
        character*40 :: fil
        print *, 'a', fil

        call config_quatsite(1,nrig)

        print *, 'b', fil
        open(1, file=fil, form='unformatted')
            write(1) x0(1:nmove), y0(1:nmove), z0(1:nmove)
        close(1, status='keep')

        print *, 'c'
        return
    end subroutine measure_idump



    subroutine measure_gather_int1d(iarr)
        implicit none
        integer*8 :: siz, ierr
        integer*8, dimension(:) :: iarr
        integer*8, dimension(:), allocatable :: ibufo
        include 'mpif.h'
        
        siz = size(iarr)
        allocate(ibufo(siz))
        ibufo = 0.0
        call mpi_allreduce(iarr, ibufo, siz, mpi_integer8, mpi_sum, mpi_comm_world, ierr)
        iarr = ibufo
        deallocate(ibufo)

        return
    end subroutine measure_gather_int1d

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

    subroutine measure_gather_int3d(iarr)
        implicit none
        integer*8 :: ierr, siz
        integer*8, dimension(:, :, :) :: iarr
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
    end subroutine measure_gather_int3d

    subroutine measure_gather_dble1d(iarr)
        implicit none
        integer*8 :: siz, ierr
        real*8, dimension(:) :: iarr
        real*8, dimension(:), allocatable :: ibufo
        include 'mpif.h'
        
        siz = size(iarr)
        allocate(ibufo(siz))
        ibufo = 0.0
        call mpi_allreduce(iarr, ibufo, siz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        iarr = ibufo
        deallocate(ibufo)

        return
    end subroutine measure_gather_dble1d

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

    subroutine measure_gather_dble3d(iarr)
        implicit none
        integer*8 :: ierr, siz
        real*8, dimension(:, :, :) :: iarr
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
    end subroutine measure_gather_dble3d

    subroutine measure_deallocate
        implicit none

        deallocate(ngz)
        if (l2d) deallocate(ng2d1, ng2d2)
        if (nsurf > 0) deallocate(ngang)
        if (lhb) then
            deallocate(laytyp)
            deallocate(hblist)
            deallocate(nbonds)
            deallocate(ibonds, iil)
        end if
        if (ldiff) then
            deallocate(xo, yo, zo)
            deallocate(diffx, diffy, diffz)
            deallocate(dbinx, dbiny, dbinz)
            deallocate(nbinz)
            if (ldiff_r) then
                deallocate(diff1, diff2, diff3)
                deallocate(sl1, sl2, sl3, sln)
                deallocate(xoo, yoo, zoo)
            end if
        end if
        if (lvisc) deallocate(vxh, ivxh)

        return
    end subroutine measure_deallocate

end module measure
