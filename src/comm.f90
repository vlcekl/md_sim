!
!  File name: comm.f90
!  Date:      2009/04/10 23:45
!  Author:    Lukas Vlcek
! 
module comm

    contains

    subroutine bcast_intx(iarr)
        implicit none
        integer*8 :: ierr
        integer*8, intent(INOUT) :: iarr
        include 'mpif.h'

        call mpi_bcast(iarr, 1, mpi_integer8, 1, mpi_comm_world, ierr)

        return
    end subroutine bcast_intx

    subroutine bcast_int1dx(iarr)
        implicit none
        integer*8 :: ierr, siz
        integer*8, dimension(:) :: iarr
        integer*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_integer8, 1, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_int1dx

    subroutine bcast_dblx(iarr)
        implicit none
        integer*8 :: ierr
        real*8, intent(INOUT) :: iarr
        include 'mpif.h'

        call mpi_bcast(iarr, 1, mpi_double_precision, 1, mpi_comm_world, ierr)

        return
    end subroutine bcast_dblx

    subroutine bcast_dbl1dx(iarr)
        implicit none
        integer*8 :: ierr, siz
        real*8, dimension(:) :: iarr
        real*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_double_precision, 1, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_dbl1dx

    subroutine bcast_dbl2dx(iarr)
        implicit none
        integer*8 :: ierr, siz
        real*8, dimension(:,:) :: iarr
        real*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_double_precision, 1, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_dbl2dx

    subroutine bcast_int2dx(iarr)
        implicit none
        integer*8 :: ierr, siz
        integer*8, dimension(:,:) :: iarr
        integer*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_integer8, 1, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_int2dx

    subroutine bcast_log(iarr)
        implicit none
        logical :: ierr
        logical, intent(INOUT) :: iarr
        include 'mpif.h'

        call mpi_bcast(iarr, 1, mpi_logical, 0, mpi_comm_world, ierr)

        return
    end subroutine bcast_log

    subroutine bcast_int(iarr)
        implicit none
        integer*8 :: ierr
        integer*8, intent(INOUT) :: iarr
        include 'mpif.h'

        call mpi_bcast(iarr, 1, mpi_integer8, 0, mpi_comm_world, ierr)

        return
    end subroutine bcast_int

    subroutine bcast_dbl(iarr)
        implicit none
        integer*8 :: ierr
        real*8, intent(INOUT) :: iarr
        include 'mpif.h'

        call mpi_bcast(iarr, 1, mpi_double_precision, 0, mpi_comm_world, ierr)

        return
    end subroutine bcast_dbl

    subroutine bcast_int1d(iarr)
        implicit none
        integer*8 :: ierr, siz
        integer*8, dimension(:) :: iarr
        integer*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_integer8, 0, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_int1d

    subroutine bcast_dbl1d(iarr)
        implicit none
        integer*8 :: ierr, siz
        real*8, dimension(:) :: iarr
        real*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_double_precision, 0, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_dbl1d

    subroutine bcast_int2d(iarr)
        implicit none
        integer*8 :: ierr, siz
        integer*8, dimension(:,:) :: iarr
        integer*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_integer8, 0, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_int2d

    subroutine bcast_dbl2d(iarr)
        implicit none
        integer*8 :: ierr, siz
        real*8, dimension(:,:) :: iarr
        real*8, dimension(:), allocatable :: ibufi
        include 'mpif.h'

        siz = size(iarr)
        allocate(ibufi(siz))
        ibufi = 0 
        ibufi = reshape(iarr, shape(ibufi))
        call mpi_bcast(ibufi, siz, mpi_double_precision, 0, mpi_comm_world, ierr)
        iarr = reshape(ibufi, shape(iarr))
        deallocate(ibufi)

        return
    end subroutine bcast_dbl2d

end module comm

!  end of comm.f90 
