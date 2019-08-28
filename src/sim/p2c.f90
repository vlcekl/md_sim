!---------------------------------------------------
!       Converts pdb file into cfg file
!---------------------------------------------------

implicit none
character*40 :: str
integer*4 :: i, j, k, nwat, nsurf, nli, ncl, nsi, noo, nmove
integer*4 :: nsistart, noostart, nsurfstart, nwatstart, ntype
real*8, dimension(:), allocatable :: x0, y0, z0, x1, y1, z1, mi
integer*8, dimension(:), allocatable :: isi
real*8 :: lx, ly, lz, temper
real*8 :: dummy, start, sumpx0, sumpy0, sumpz0, sumvsq, star, PH_R, dt

nwat = 258
nli = 0
ncl = 0
nsi = 1072
noo = 2117
nsurf = 54
nmove = nwat*3 + nsurf*2 + nli + ncl + nsi + noo
ntype = 5

nsistart = 0
noostart = nsi-nsurf
nsurfstart = noostart + noo
nwatstart = nsurfstart + 3*nsurf


temper = 300.0

lx = 38.154
ly = 38.154
lz = 38.154

allocate(x0(nmove), y0(nmove), z0(nmove))
allocate(x1(nmove), y1(nmove), z1(nmove), mi(nmove))
x0 = 0. ; y0 = 0. ; z0 = 0.
x1 = 0. ; y1 = 0. ; z1 = 0.
allocate(isi(nmove))
isi = 0

open(2, file = 'y.pdb', status = 'old')

    do i = 1, nmove
        read(2, *) str, j, str, j, x0(i), y0(i), z0(i), str
        if (x0(i) > lx/2.0) x0(i) = x0(i) - lx
        if (y0(i) > ly/2.0) y0(i) = y0(i) - ly
        if (z0(i) > lz/2.0) z0(i) = z0(i) - lz
    end do

close(2, status = 'keep')


open(1, file = 'y.cfg', status = 'unknown')

    write(1, *)'# SIMULATION BOX SIZE'
    write(1, *) lx, ly, lz
    write(1, *) ''
    write(1, *)'# NUMBER OF PARTICLE TYPES'
    write(1, *) ntype
    write(1, *) ''
    write(1, *)'# PARTICLE COORDINATES'

    k = 0
    ! water
    write(1, *) 'spce ', 'rigmol', nwat
    do i = nwatstart + 1, nwatstart + 3*nwat, 3
        k = k + 1
        write(1, *) k, 1, 2, 2
        write(1, *)  x0(i), y0(i), z0(i)
        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.

        write(1, *) ' -0.87675800E+00 -0.44312204E+00  0.34852716E-01 -0.18364030E+00 '
        write(1, *) '  0.41033240E-02 -0.62047122E-02 -0.57192250E-03 -0.47272334E-02 '
        write(1, *) '  0.87855904E-04 -0.74536571E-04 -0.23239829E-04 -0.31416875E-04 '
        write(1, *) ' -0.24818947E-05  0.33148425E-05  0.90744827E-06  0.96187011E-05 '
        write(1, *) ' -0.98686539E-07  0.19854743E-07  0.11440026E-06  0.16099893E-06 '

        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.
        write(1, *)  0., 0., 0.
    end do
    write(1, *) ''
    k = 3*k

    write(1, *) 'Li ', 'atmol', nli
    write(1, *) ''

    write(1, *) 'Cl ', 'atmol', ncl
    write(1, *) ''

    ! oh
    write(1, *) 'tOH ', 'surfmol', nsurf
    do i = nsurfstart + 1, nsurfstart + 3*nsurf, 3
        k = k + 1
        write(1, *) k, 1
        write(1, *)  x0(i+1), y0(i+1), z0(i+1)
        write(1, *)  x0(i+1), y0(i+1), z0(i+1)
        k = k + 1
        write(1, *) k, 2
        write(1, *)  x0(i+2), y0(i+2), z0(i+2)
        write(1, *)  x0(i+2), y0(i+2), z0(i+2)
    end do
    write(1, *) ''

    ! si
    write(1, *) 'rutile ', 'crystal', 1
    do i = nsistart + 1, nsistart + nsi-nsurf, 1
        k = k + 1
        write(1, *) k, 1
        write(1, *)  x0(i), y0(i), z0(i)
    end do
    j = 0
    do i = nsurfstart + 1, nsurfstart + 3*nsurf, 3
        k = k + 1
        write(1, *) k, 1
        write(1, *)  x0(i), y0(i), z0(i)
        j = j + 1
        isi(j) = k
    end do
    do i = noostart + 1, noostart + noo, 1
        k = k + 1
        write(1, *) k, 2
        write(1, *)  x0(i), y0(i), z0(i)
    end do
    write(1, *) ''
    ! thermostat
    write(1, *)'# THERMOSTAT'
    write(1, *)  0., 0., 0., 0., 0.
    write(1, *)  0., 0., 0., 0., 0.

close(1, status = 'keep')

open(1, file = 'y.top', status = 'unknown')
    write(1, *)'# Atom connections, specific interactions'

    write(1, *) 1, 'spce ', 'rigmol', nwat
    write(1, *)'sites:', 3

    write(1, *) 2, 'Li ', 'atmol', nli
    write(1, *)'no constraints'
    write(1, *) 3, 'Cl ', 'atmol', ncl
    write(1, *)'no constraints'

    write(1, *) 4, 'tOH ', 'surfmol', nsurf
    write(1, *)'constraints:', 2
    j = 0
    do i = 3*nwat + nli + ncl + 1, 3*nwat + nli + ncl + 2*nsurf, 2
        j = j + 1
        write(1, *) isi(j), i
        write(1, *) i, i+1
    end do
    write(1, *)'angles:', 1
    j = 0
    do i = 3*nwat + nli + ncl + 1, 3*nwat + nli + ncl + 2*nsurf, 2
        j = j + 1
        write(1, *) isi(j), i, i+1
    end do
    write(1, *) 5, 'rutile ', 'crystal', 1

close(1, status = 'keep')

deallocate(x0, y0, z0)

stop
10  format(A6,I5,A1, A4, A14,  3f8.3)
20  format(A6,I5)
30  format(A6)
end

!  end of readcfg.f90 
