module procedures
  implicit none
contains

  subroutine linked_list(data, boxsize, ngrid, ll, lirst, rgrid)
    implicit none
    integer*8 :: ndata
    integer*8 :: i, ipx, ipy, ipz
    integer*8, intent(in) :: ngrid
    real*8, intent(in) :: boxsize
    real*8, intent(out) :: rgrid
    real*8, dimension(:,:), intent(in) :: data
    integer*8, dimension(:,:,:), allocatable, intent(out) :: lirst
    integer*8, dimension(:), allocatable, intent(out) :: ll

    ndata = size(data, dim=2)
    allocate(lirst(ngrid, ngrid, ngrid))
    allocate(ll(ndata))
    rgrid = boxsize / real(ngrid)
    
    lirst = 0
    ll = 0
    do i = 1, ndata
      ipx = int(data(1, i) / rgrid + 1.)
      ipy = int(data(2, i) / rgrid + 1.)
      ipz = int(data(3, i) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i

    end do

    do i = 1, ndata
      ipx = int(data(1, i) / rgrid + 1.)
      ipy = int(data(2, i) / rgrid + 1.)
      ipz = int(data(3, i) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
      endif
    end do

  end subroutine linked_list

  subroutine read_unformatted(input_filename, data, weight, ndata, has_velocity)
    implicit none
    integer*8 :: nrows, ncols
    character(len=500), intent(in) :: input_filename
    integer*8, intent(out) :: ndata
    real*8, allocatable, dimension(:,:), intent(out) :: data
    real*8, allocatable, dimension(:), intent(out) :: weight
    logical, intent(out) :: has_velocity

    has_velocity = .false.

    open(20, file=input_filename, status='old', form='unformatted')
    read(20) nrows
    read(20) ncols
    allocate(data(ncols, nrows))
    allocate(weight(nrows))
    read(20) data
    close(20)
    ndata = nrows
    weight = 1.0 ! default weights
    if (ncols .ge. 4) then
      if (ncols .eq. 4) weight = data(4, :)
      if (ncols .eq. 6) has_velocity = .true.
      if (ncols .eq. 7) then
        weight = data(7, :)
        has_velocity = .true.
      end if
    end if
      
  end subroutine read_unformatted

  subroutine binning(rmin, rmax, nbin, bin, bin_edges, rwidth)
    implicit none

    integer*8 :: i
    integer*8, intent(in) :: nbin
    real*8, intent(in) :: rmin, rmax
    real*8, intent(out) :: rwidth
    real*8, allocatable, dimension(:), intent(out) :: bin, bin_edges

    allocate(bin(nbin))
    allocate(bin_edges(nbin + 1))

    rwidth = (rmax - rmin) / nbin
    do i = 1, nbin + 1
      bin_edges(i) = rmin + (i - 1) * rwidth
    end do
    do i = 1, nbin
      bin(i) = bin_edges(i + 1) - rwidth / 2.
    end do

  end subroutine binning

end module procedures