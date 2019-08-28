module nnlist
    use config

    integer*4, dimension(:), allocatable :: list1, list2, list3, list4
    integer*4, dimension(:), allocatable :: listw1, listw2, listw3, listw4
    integer*4 :: mypair, wmypair, ipack, wpack
    integer*4 :: icount = 0
    real*8, dimension(:), allocatable :: ox0, oy0, oz0, ocx0, ocy0, ocz0           ! old positions (for nn_list)
    real*8 :: skin, rskinsq, dshsq
    integer*4 :: nnmax, nnwmax

    contains

    subroutine nnlist_init
        implicit none
        integer*4 :: i, j, nrow, ncol, myrow, mycol, jmin, jmax, inmin, inmax, icnt

        allocate(ocx0(nrig), ocy0(nrig), ocz0(nrig))
!        allocate(ox0(1:nmove), oy0(1:nmove),oz0(1:nmove))
        allocate(ox0(nriga+1:nmove), oy0(nriga+1:nmove),oz0(nriga+1:nmove))

!        rskinsq = (rcut + skin)**2
        rskinsq = rcut + skin
        ! correction for the differece between O atom and COM
        rskinsq = rskinsq + 2.0*sqrt((x0(1)-cx0(1))**2 + (y0(1)-cy0(1))**2 + (z0(1)-cz0(1))**2)
        rskinsq = rskinsq*rskinsq

        dshsq = dshx**2 + dshy**2 + dshz**2  ! max movement of the wall (Couette flow)

    ! moving - moving
        nnmax = (nrig + nmove - nriga)*(nrig + nmove - nriga + 1)/(2*nproc) + 1
        allocate(list3(nnmax), list4(nnmax))

        print *, 'nn', me, nmove, nliqa, nnmax, nrig, nriga
        icnt = 0
        mypair = 0
        do i = 1, nmove-1
            if (nsi(itype(i)) == 0) cycle
            do j = i+1, nmove
                if (nsi(itype(j)) == 0) cycle
                if (i > nliqa .and. j > nliqa ) then
                    if (nnlist_constr(i, j)) cycle
                end if
                if (mod(icnt, nproc) == me) then
                    mypair = mypair + 1
                    if (mypair > nnmax) then
                        write(*,*) 'mypair-nnmax', me, mypair, nnmax, icnt, mod(icnt, nproc)
                        stop 'mypair too large'
                    end if
                    list3(mypair) = i
                    list4(mypair) = j
                end if
                icnt = icnt + 1
            end do
        end do
        write(*,*) 'mypair-nnmax', me, mypair, nnmax
        nnmax = nnmax/1
        allocate(list1(nnmax), list2(nnmax))

    ! moving - frozen
        nnwmax = (nrig + nmove - nriga)*nwall/nproc + 1
        allocate(listw3(nnwmax), listw4(nnwmax))

        icnt = 0
        wmypair = 0
        do i = 1, nmove
            if (nsi(itype(i)) == 0) cycle
            do j = nmove+1, ntota
                if (i > nliqa) then
                    if (nnlist_constr(i, j)) cycle
                end if
                if (mod(icnt, nproc) == me) then
                    wmypair = wmypair + 1
                    listw3(wmypair) = i
                    listw4(wmypair) = j
                end if
                icnt = icnt + 1
            end do
        end do
        if (wmypair > nnwmax) then
            write(*,*) me, wmypair, nnwmax
            stop 'nnlist_init(): nnlist, iw too big!'
        end if
        nnwmax = nnwmax/1
        allocate(listw1(nnwmax), listw2(nnwmax))

        call nnlist_make

        return
    end subroutine nnlist_init

    subroutine nnlist_make
        implicit none
        integer*4 :: k, left, i, iw
        real*8 :: xijs, yijs, zijs

        icount = icount + 1

        forall (i = 1:nrig)
!            cx0(i) = cx0(i) - lx*aint(cx0(i)*lxhi)
            ocx0(i) = cx0(i)
!            cy0(i) = cy0(i) - ly*aint(cy0(i)*lyhi)
            ocy0(i) = cy0(i)
            ocz0(i) = cz0(i)
        end forall
!        forall (i = 1:nriga:3)
!            x0(i) = x0(i) - lx*aint(x0(i)*lxhi)
!            ox0(i) = x0(i)
!            y0(i) = y0(i) - ly*aint(y0(i)*lyhi)
!            oy0(i) = y0(i)
!            oz0(i) = z0(i)
!        end forall
        forall (i = nriga+1:nmove)
!            x0(i) = x0(i) - lx*aint(x0(i)*lxhi)
            ox0(i) = x0(i)
!            y0(i) = y0(i) - ly*aint(y0(i)*lyhi)
            oy0(i) = y0(i)
            oz0(i) = z0(i)
        end forall

        ! liquid - liquid
        left = 0
        do k = 1, mypair                
            xijs = x0(list3(k)) - x0(list4(k))
            xijs = xijs - lx*anint(xijs*lxi)
            yijs = y0(list3(k)) - y0(list4(k))
            yijs = yijs - ly*anint(yijs*lyi)
            zijs = z0(list3(k)) - z0(list4(k))
!            zijs = zijs - lz*anint(zijs*lzi)
            if(xijs*xijs + yijs*yijs + zijs*zijs < rskinsq) then
                left = left + 1
                list1(left) = list3(k)
                list2(left) = list4(k)
            end if
        end do
        if (left > nnmax) then
            write(*, *) me, nnmax, left
            stop 'nnlist_make(): nnlist too big!'
        end if
        ipack = left

        ! liquid (+surface) - solid
        left = 0
        do k = 1, wmypair
            xijs = x0(listw3(k)) - x0(listw4(k))
            xijs = xijs - lx*anint(xijs*lxi)
            yijs = y0(listw3(k)) - y0(listw4(k))
            yijs = yijs - ly*anint(yijs*lyi)
            zijs = z0(listw3(k)) - z0(listw4(k))
!            zijs = zijs - lz*anint(zijs*lzi)
            if(xijs*xijs + yijs*yijs + zijs*zijs < rskinsq) then
                left = left + 1
                listw1(left) = listw3(k)
                listw2(left) = listw4(k)
            end if
        end do
        if (left > nnwmax) then
            write(*, *) me, nnwmax, left
            stop 'nnlist_make(): wnnlist too big!'
        end if
        wpack = left

        return
    end subroutine nnlist_make

    ! test if the pair i-j is constrained
    logical function nnlist_constr(i, j)
        implicit none
        integer*4 :: i, j, it, ix, iy, is, ii, jj, kk

!        ix = 0
!        do it = 1, ntype
!            if (nangt(it) > 0) then
!            do iy = 1, nt(it)
!                ix = ix + 1
!                    do is = 1, nconst(it)
!                        ii = ang1(is, ix)
!                        jj = ang2(is, ix)
!                        kk = ang3(is, ix)
!                        if ((i == ii .or. i == jj .or. i == kk) .and. (j == ii .or. j == jj .or. j == kk)) then
!                            nnlist_constr = .true.
!                            return
!                        end if
!                    end do
!                end do
!            end if
!        end do

        ix = 0
        do it = 1, ntype
            if (nconst(it) > 0) then
            do iy = 1, nt(it)
                ix = ix + 1
                    do is = 1, nconst(it)
                        ii = constra(is, ix)
                        jj = constrb(is, ix)
                        if ((i == ii .and. j == jj) .or. (i == jj .and. j == ii)) then
                            nnlist_constr = .true.
                            return
                        end if
                    end do
                end do
            end if
        end do

        nnlist_constr = .false.
        return
    end function nnlist_constr

end module nnlist

