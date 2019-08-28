module sysdef
    use global
    use comm
    implicit none

    integer*4, parameter :: AMX = 12     ! number of relaxed atom types
    integer*4, parameter :: TMXSQ = (TMX+1)*TMX/2

    ! time
    real*8 :: dt

    ! molecules
    character*40 :: filfld, filtop
    character*20, dimension(TMX) :: moltype, molname
    integer*8, dimension(TMX) :: moltypen
    integer*8 :: ntype, nt(0:TMX)
    integer*8 :: nltype, nstype
    integer*8 :: nrig, nliq, nmol

    ! atoms
    character*2, dimension(TMX, SMX) :: sitnam
    integer*8 :: nit, nijt
    integer*8 :: nriga, nliqa, nsurf, nmove, nwall, ntota, ngogo
    integer*8 :: ncons, nang
    integer*8 :: nshift, npull, nstatic
    integer*8, dimension(TMX) :: nst, ns, nsi  !nst, ns - moltypes, nsi - attypes
    integer*8, dimension(TMX, SMX) :: nss
    real*8, dimension(TMX, SMX) :: ms, qms           ! site mass, charge, & LJ pars. 
    real*8, dimension(TMX, SMX) :: psx, psy, psz    ! site coordinates in principal axis frame
    real*8, dimension(AMX) :: qshift

    ! general (physical) interactions (export)
!    integer*8, dimension(NMX) :: itype
    integer*8, dimension(:), allocatable :: itype, walltype
    integer*8, dimension(TMX, TMX) :: ijtype
    character*8, dimension(TMXSQ) :: vdwtyp
    integer*8, dimension(TMXSQ) :: vdwnum
    real*8, dimension(TMXSQ) :: eps, sig, aa, rhoi, cc
    integer*8, dimension(1000) :: iatype, ictype
    real*8 :: rcut

    ! specific (chemical) interactions and constraints (export)
    integer*8, dimension(TMX) :: nconst, nangt
    real*8, dimension(TMX, SMX) :: length, kt, theta
    integer*8, dimension(SMX, 1000) :: constra, constrb, ang1, ang2, ang3
    real*8 :: cx, cy, cz, lw, lwh

    ! derived quantities (practical => export all)
!    real*8, dimension(NMX) :: mi, mii, qs         ! molecular masses, inverse mass & charge
    real*8, dimension(:), allocatable :: mi, mii, qs         ! molecular masses, inverse mass & charge
    real*8, dimension(TMX) :: pmix, pmiy, pmiz    ! principal moments of inertia
    real*8, dimension(SMX) :: xbody, ybody, zbody ! water site coordinates
    real*8 :: mw, inxx, inyy, inzz                ! water mass & moments of inertia

    logical :: lchempot
    logical :: lfext, lcouett
    real*8 :: fxext, fyext, fzext
    real*8 :: dshx, dshy, dshz
    integer*8 :: nwtype

    contains

    subroutine sysdef_zero 
        implicit none

        ! zero
        nst = 0 ; nss = 0 ; ns = 0 ; nsi = 0
        ms = 0.0
        psx = 0.0; psy = 0.0; psz = 0.0
        length = 0.0 ; kt = 0.0 ; theta = 0.0
        nconst = 0 ; nangt = 0
        moltypen = 0

        eps = 0.0 ; sig = 0.0
        aa = 0.0 ; rhoi = 0.0 ; cc = 0.0
        
        return
    end subroutine sysdef_zero 

    !========================================================================
    !      Reads molecular parameters from the field file
    !========================================================================
    subroutine sysdef_input
        implicit none
        character*80 :: str
        character*20 :: mname, mtype
        integer*4 :: it, is, ijt, ic, ia, ntx, nsx, i, ix, jx
        integer*8 :: nrtype
        real*8 :: six, epx

        if (me == 0) then
        open(4, file = filfld, status = 'old')
            ! molecules
            read(4, *) str  ! molecules and atoms
            read(4, *) ntype
            do it = 1, ntype
                read(4, *) str  ! # Molecule
                read(4, *) ix, molname(it), moltype(it), nt(it)
                if (ix /= it) stop 'sysdef_input(): wrong molecule number!'

                read(4, *) str   ! sites
                read(4, *) nst(it)
                do is = 1, nst(it)
                    read(4, *) sitnam(it, is), nss(it, is)
                    read(4, *) ms(it, is), qms(it, is)
                end do
                ns(it) = sum(nss(it, 1:nst(it)))

                select case(moltype(it))
                    case ('rigmol')
                        read(4, *) str          ! geometry
                        do is = 1, ns(it)
                            read(4, *) psx(it, is), psy(it, is), psz(it, is)
                        end do
                    case ('atmol')
                    case ('surfmol')
                        ! constraints
                        read(4, *) str
                        read(4, *) nconst(it)
                        do is = 1, nconst(it)
                            read(4, *) length(it, is)
                        end do
                        ! angles
                        read(4, *) str
                        read(4,*) nangt(it)
                        do is = 1, nangt(it)
                            read(4, *) kt(it, is), theta(it, is)
                        end do
                    case ('crystal')
                        read(4, *) ns(it)
                        read(4, *) cx, cy, cz
                        read(4, *) lw
                        read(4, *) nshift
                        do is = 1, nshift
                            read(4, *) qshift(is)
                        end do
                    case default
                        write(*,*) moltype(it)
                        stop 'sysdef_input(): Unknown moltype!'
                end select
            end do

            ! vdw interactions (non-coulombic)
            read(4, *) str
            ntx = sum(nst(1:ntype))
            do i = 1, ntx*(ntx+1)/2
                read(4, *) ijt
                read(4, *) ix, str, str
                read(4, *) jx, str, str
                read(4, *) vdwtyp(ijt)
                select case (vdwtyp(ijt))
                    case ('zero')
                    case ('12-6')
                        read(4, *) sig(ijt), eps(ijt)
                        six = (sig(ijt)/eps(ijt))**(1.0/6.0)
                        epx = eps(ijt)**2/(4.0*sig(ijt))
                        sig(ijt) = six
                        eps(ijt) = epx
                        vdwtyp(ijt) = 'LeJo'
                    case ('LeJo')
                        read(4, *) sig(ijt), eps(ijt)
                    case ('Buck')
                        read(4, *) aa(ijt), rhoi(ijt), cc(ijt)
                    case default
                        write(*, *) vdwtyp(ijt)
                        stop 'sysdef_input(): Unknown type of interaction!'
                end select
                ijtype(ix, jx) = ijt
            end do

        close(4)

        ! topofile
        open(1, file=filtop, status='old')
            read(1, *) str
            ia = 0
            ic = 0
            do it = 1, ntype
                read(1, *) ix, mname, mtype, ntx
                if ((ix /= it) .or. (mname .ne. molname(it)) .or. (mtype .ne. moltype(it)) .or. (ntx /= nt(it))) then
                    write(*, *) ix, it
                    write(*, *) mname, molname(it)
                    write(*, *) mtype, moltype(it)
                    write(*, *) ntx, nt(it)
                    stop 'sysdef_input(): Error in topofile!'
                end if
                select case (mtype)
                    case ('rigmol')
                        read(1, *) str, nsx
                        if (nsx .ne. ns(it)) stop 'sysdef_input(): Error in topofile (sites)!'
                    case ('atmol')
                        read(1, *) str
                    case ('surfmol')
                        read(1, *) str, nsx
                        if (nsx .ne. nconst(it)) stop 'sysdef_input(): Error in topofile (nconst)!'
                        do i = 1, ntx
                            ic = ic + 1
                            do is = 1, nsx
                                ictype(ic) = it
                                read(1, *) constra(is, ic), constrb(is, ic)
                            end do
                        end do
                        read(1, *) str, nsx
                        if (nsx .ne. nangt(it)) stop 'sysdef_input(): Error in topofile (nangt)!'
                        if (nsx > 0) then
                            do i = 1, ntx
                                ia = ia + 1
                                do is = 1, nsx
                                    iatype(ia) = it
                                    read(1, *) ang1(is, ia), ang2(is, ia), ang3(is, ia)
                                end do
                            end do
                        end if
                end select
            end do
        close(1, status='keep')
        ncons = ic
        nang = ia
        if (ncons > 1000 .or. nang > 1000) stop 'sysdef_input(): too large number of constraints or angles ic or ia'

        ! moltypes
        nrtype = 0
        nltype = 0
        nstype = 0
        nwtype = 0
        do it = 1, ntype
            select case (moltype(it))
                case ('rigmol')
                    moltypen(it) = 0
                    nrtype = it
                    nltype = it
!                    nstype = it
!                    nwtype = it
                case ('atmol')
                    moltypen(it) = 1
                    nltype = it
!                    nstype = it
!                    nwtype = it
                case ('surfmol')
                    moltypen(it) = 2
                    nstype = it
!                    nwtype = it
                case ('crystal')
                    moltypen(it) = 3
                    nwtype = it
                case default
                    stop 'sysdef_init(): Unknown moltype!'
            end select
        end do
        end if
!        if (nstype /= ntype-1 .or. nwtype /= ntype) stop 'sysdef_init(): Wrong molype numbers!'
        call bcast_int(ntype)
        call bcast_int(nrtype)
        call bcast_int(nltype)
        call bcast_int(nstype)
        call bcast_int(nwtype)
        call bcast_int(nang)
        call bcast_int(ncons)
        call bcast_int1d(nst)
        call bcast_int1d(nt)
        call bcast_int1d(ns)
        call bcast_int1d(moltypen)
        call bcast_int1d(nconst)
        call bcast_int1d(nangt)
        call bcast_int1d(nangt)
        call bcast_dbl(cx)
        call bcast_dbl(cy)
        call bcast_dbl(cz)
        call bcast_dbl(lw)
        call bcast_int(nshift)
        call bcast_dbl1d(qshift)
        call bcast_dbl1d(sig)
        call bcast_dbl1d(eps)
        call bcast_dbl1d(aa)
        call bcast_dbl1d(rhoi)
        call bcast_dbl1d(cc)
        call bcast_int2d(nss)
        call bcast_int2d(ijtype)
        call bcast_dbl2d(ms)
        call bcast_dbl2d(qms)
        call bcast_dbl2d(psx)
        call bcast_dbl2d(psy)
        call bcast_dbl2d(psz)
        call bcast_dbl2d(length)
        !call bcast_dbl2d(uharm)
        call bcast_dbl2d(kt)
        call bcast_dbl2d(theta)
        

        call bcast_int1d(ictype)
        call bcast_int1d(iatype)
        call bcast_int2d(constra)
        call bcast_int2d(constrb)
        call bcast_int2d(ang1)
        call bcast_int2d(ang2)
        call bcast_int2d(ang3)

        ! # of molecules
        nt(0) = 0
        nrig = sum(nt(1:nrtype))
        nliq = sum(nt(1:nltype))
        nmol = sum(nt(1:ntype))
        ! # of atoms
        nriga = sum(ns(1:nrtype)*nt(1:nrtype))
        nliqa = sum(ns(1:nltype)*nt(1:nltype))
        nsurf = sum(ns(nltype+1:nstype)*nt(nltype+1:nstype))
        nmove = nliqa + nsurf
        nwall = ns(nwtype)*nt(nwtype)
        print *, 'nwall', nwall
        npull = 0
        if (lcouett) npull = nwall/2
        ngogo = nmove + npull
        nstatic = nwall - npull
        ntota = nmove + nwall
        if (nt(nwtype) > 1) stop 'sysdef_init(): wrong number of crystals'
        if (nstatic < 0) stop 'sysdef_init(): wrong number of npull'
        ! # of vdw interactions
        nit = sum(nst(1:ntype))
        nijt = nit*(nit+1)/2

        if (me == 0) then
        ! interaction types
            do ijt = 1, nijt
                select case (vdwtyp(ijt))
                    case ('zero')
                        vdwnum(ijt) = 0
                    case ('LeJo')
                        vdwnum(ijt) = 1
                    case ('Buck')
                        vdwnum(ijt) = 2
                    case default
                        vdwnum(ijt) = -1
                end select
            end do
        end if

        call bcast_int1d(vdwnum)

        return
    end subroutine sysdef_input

    subroutine sysdef_alloc
        implicit none

        allocate(itype(ntota), qs(ntota), mi(ntota), mii(ntota))
        itype = 0 ; qs = 0.0 ; mi = 0.0 ; mii = 0.0
        allocate(walltype(nmove+1:ntota))
        walltype = 0

        return
    end subroutine sysdef_alloc
    !========================================================================
    !              Initialisation of all molecular interactions
    !        (derived & special quantities + some useful parameters)
    !========================================================================
    subroutine sysdef_init
        implicit none
        integer*4 :: it, jt, is, ijt, i, iit, ix, ist, jjt
        real*8 :: v_unit


        ijt = 0
        do it = 1, nit
            do jt = it, nit
                ijt = ijt + 1
                ijtype(it, jt) = ijt
                ijtype(jt, it) = ijt
            end do
        end do
        if (ijt /= nijt) stop 'sysdef_init(): something rotten in the interaction types!'

        eps = 4.0*eps
        sig = sig**6
        forall (it = 1:nijt, rhoi(it) /= 0.0) rhoi(it) = 1.0/rhoi(it)

        theta = theta*M_PI/180.0
        ! atom type and mass and charge
        v_unit = PH_A2m/(dt*1.0e-15) ! velocity: A/timestep->SI/1000
        ms = ms*v_unit**2/1000.0

        i = 0    ! atom
        jjt = 0  ! interaction type
        do it = 1, max(nltype, nstype)
            do ix = 1, nt(it)
                do ist = 1, nst(it)
                    iit = jjt + ist
                    if (moltypen(it) == 0) then
                        nsi(iit) = 0
                        if (ist == 1) nsi(iit) = ns(it)
                    else
                        nsi(iit) = 1
                    end if
                    do is = 1, nss(it, ist)
                        i = i + 1
                        itype(i) = iit
                        mi(i) = ms(it, ist)
                        qs(i) = qms(it, ist)
                    end do
                end do
            end do
            jjt = jjt + nst(it)
        end do
        if (i /= nmove) stop 'sysdef_init(): i /= nmove !'
        if (jjt + nst(nwtype) /= nit) stop 'sysdef_init(): wrong number of interaction types !'

        do i = nmove+1, ntota    ! crystal site types from config file
            is = itype(i)
            if (is > 0) then
                qs(i) = qms(ntype, is)
            else
                qs(i) = qshift(-is)
                if (-is <= nss(ntype, 1)) then
                    is = 1
                else
                    is = 2
                end if
            end if
            itype(i) = jjt + is ! remap surface types to vdw types
        end do
!        do i = 1, ntota
!            print *, i, itype(i)
!        end do

        forall (i = 1:nmove) mii(i) = 1.0/mi(i)
        forall (i = nmove+1:ntota) mii(i) = 0.0    ! frozen sites = infinite mass

        ! rigid bodies
        ! suppose that the tenzor is already diagonal
        pmix = 0.0 ; pmiy = 0.0 ; pmiz = 0.0
        do it = 1, ntype
            if (moltypen(it) == 0) then
                ix = sum(nt(1:it-1))
                ! geometry
                psx(it,1:ns(it)) = psx(it,1:ns(it)) - sum(psx(it,1:ns(it))*mi(ix+1:ix+ns(it)))/sum(mi(ix+1:ix+ns(it)))     ! coordinates relative to COM
                psy(it,1:ns(it)) = psy(it,1:ns(it)) - sum(psy(it,1:ns(it))*mi(ix+1:ix+ns(it)))/sum(mi(ix+1:ix+ns(it)))
                psz(it,1:ns(it)) = psz(it,1:ns(it)) - sum(psz(it,1:ns(it))*mi(ix+1:ix+ns(it)))/sum(mi(ix+1:ix+ns(it)))
                pmix(it) = sum(mi(ix+1:ix+ns(it))*(psy(it,1:ns(it))**2 + psz(it,1:ns(it))**2))
                pmiy(it) = sum(mi(ix+1:ix+ns(it))*(psz(it,1:ns(it))**2 + psx(it,1:ns(it))**2))
                pmiz(it) = sum(mi(ix+1:ix+ns(it))*(psx(it,1:ns(it))**2 + psy(it,1:ns(it))**2))
            end if
        end do

        ! specific - only for water
        mw = sum(nss(1, 1:nst(1))*ms(1, 1:nst(1)))   ! mass of a water molecule
        forall (is = 1:ns(1))        ! water site principal axis coordinates
            xbody(is) = psx(1, is)
            ybody(is) = psy(1, is)
            zbody(is) = psz(1, is)
        end forall
        inxx = pmix(1)               ! water moments of inertia
        inyy = pmiy(1)
        inzz = pmiz(1)

        lwh = 0.5*lw

        return
    end subroutine sysdef_init

    subroutine sysdef_output
        use global
        implicit none
        integer*4 :: it, jt, ijt, ix, jx, is, js, ia, ic, i, jsmin
 
        
        ! FIELD FILE
        open(1, file=filfld, status='unknown')
 
            write(1, *) '# MOLECULES'
            write(1, '(i3)') ntype
            do it = 1, ntype
 
                write(1, *) '# Molecule'
                write(1, '(i2,1x,a6,1x,a8,1x,i6)') it,  molname(it), moltype(it), nt(it)
 
                write(1, *) '# sites'
                write(1, '(i3)') nst(it)
                do is = 1, nst(it)
                    write(1, '(a2,1x,i3)') sitnam(it, is), nss(it, is)
                    write(1, '(es13.5,1x,es13.5)') ms(it, is), qms(it, is)
                end do
 
                select case(moltype(it))
                    case ('rigmol')
                        write(1, *) '# geometry'
                        ns(it) = sum(nss(it, 1:nst(it)))
                        do is = 1, ns(it)
                            write(1, '(3(es13.5,1x))') psx(it, is), psy(it, is), psz(it, is)
                        end do
                    case ('atmol')
                    case ('surfmol')
                        write(1, *) '# constraints'
                        write(1, '(i3)') nconst(it)
                        do is = 1, nconst(it)
                            write(1, '(es13.5)') length(it, is)
                        end do
                        write(1, *) '# angles'
                        write(1, '(i3)') nangt(it)
                        do is = 1, nangt(it)
                            write(1, '(2(es13.5,1x))') kt(it, is), theta(it, is)
                        end do
                    case ('crystal')
                        ! surface layers
                        write(1, '(i5)') ns(it)
                        write(1, '(3(es13.5,1x))') cx, cy, cz
                        write(1, '(es13.5)') lw
                        write(1, '(i3)') nshift
                        do is = 1, nshift
                            write(1, '(es13.5)') qshift(is)
                        end do
                    case default
                        write(*, *) moltype(it)
                        stop 'field_out(): Unknown moltype!'
                end select
            end do
 
            ! vdw interactions (non-coulombic)

            sig = sig**(1.0/6.0) ! from internal form
            eps = eps/4.0
            forall (it = 1:nijt, rhoi(it) /= 0.0) rhoi(it) = 1.0/rhoi(it)

                            
            write(1, *) '# VDW INTERACTIONS'
            ijt = 0
            ix = 0
            do it = 1, ntype
                do is = 1, nst(it)
                    ix = ix + 1
                    jx = ix - 1
                    do jt = it, ntype
                        if (it == jt) then
                            jsmin = is
                        else
                            jsmin = 1
                        end if
                        do js = jsmin, nst(jt)
                            jx = jx + 1
                            ijt = ijt + 1
                            write(1, '(i3)') ijt
                            write(1, '(i3,1x,a6,1x,a2)') ix, molname(it), sitnam(it, is)
                            write(1, '(i3,1x,a6,1x,a2)') jx, molname(jt), sitnam(jt, js)
                            write(1, '(a4)') vdwtyp(ijt)

                            select case (vdwtyp(ijt))
                                case ('zero')
                                case ('LeJo')
                                    write(1, '(2(es13.5,1x))') sig(ijt), eps(ijt)
                                case ('Buck')
                                    write(1, '(3(es13.5,1x))') aa(ijt), rhoi(ijt), cc(ijt)
                                case ('spec')
                                    write(1, *) 'insert'
                                case default
                                    write(*, *) it, jt, ijt, vdwtyp(ijt)
                                    stop 'field_out(): Unknown vdwtype!'
                            end select
                            write(1, *)
                        end do
                    end do
                end do
            end do

            sig = sig**6    ! back to internal form
            eps = 4.0*eps
            forall (it = 1:nijt, rhoi(it) /= 0.0) rhoi(it) = 1.0/rhoi(it)
 
        close(1, status='keep')
 
        ! topofile
        open(1, file=filtop)
            write(1, *) '# Atom connections, specific interactions'
            ia = 0
            ic = 0
            do it = 1, ntype
                write(1, *) it, molname(it), moltype(it), nt(it)
                select case (moltype(it))
                    case ('rigmol')
                        write(1, *) 'sites:', ns(it)
                    case ('atmol')
                        write(1, *) 'no connections'
                    case ('surfmol')
                        write(1, *) 'constraints:', nconst(it)
                        do i = 1, nt(it)
                            ic = ic + 1
                            ictype(ic) = it
                            do is = 1, nconst(it)
                                write(1, *) constra(is, ic), constrb(is, ic)
                            end do
                        end do
                        write(1, *) 'angles:', nangt(it)
                        do i = 1, nt(it)
                            ia = ia + 1
                            iatype(ia) = it
                            do is = 1, nangt(it)
                                write(1, *) ang1(is, ia), ang2(is, ia), ang3(is, ia)
                            end do
                        end do
                end select
            end do
 
        close(1, status='keep')
 
        return
    end subroutine sysdef_output
        
end module sysdef

