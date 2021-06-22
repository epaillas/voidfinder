program mean_radial_velocity_vs_r
  use OMP_LIB
  implicit none
  
  real*8 :: rgrid, boxsize
  real*8 :: disx, disy, disz, dis, dis2, vr
  real*8 :: velx, vely, velz
  real*8 :: rwidth, dim1_max, dim1_min, dim1_min2, dim1_max2
  
  integer*8 :: ntracers, ncentres, dim1_nbin, rind
  integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*8 :: indx, indy, indz, nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid
  integer*4 :: nthreads
  
  integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, dimension(3) :: r, vel
  real*8, allocatable, dimension(:,:)  :: tracers, centres
  real*8, dimension(:), allocatable :: DD
  real*8, dimension(:), allocatable :: VV, mean_vr
  real*8, dimension(:), allocatable :: rbin, rbin_edges
  real*8, dimension(:, :), allocatable :: DD_i, VV_i

  logical :: has_velocity = .false.
  logical :: debug = .true.
  
  character(20), external :: str
  character(len=500) :: data_filename, data_filename_2, output_filename, nthreads_char
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char, box_char
  
  if (iargc() .ne. 9) then
      write(*,*) 'Some arguments are missing.'
      write(*,*) '1) data_filename'
      write(*,*) '2) data_filename_2'
      write(*,*) '3) output_filename'
      write(*,*) '4) boxsize'
      write(*,*) '5) dim1_min'
      write(*,*) '6) dim1_max'
      write(*,*) '7) dim1_nbin'
      write(*,*) '8) ngrid'
      write(*,*) '9) nthreads'
      write(*,*) ''
      stop
    end if
    
  call getarg(1, data_filename)
  call getarg(2, data_filename_2)
  call getarg(3, output_filename)
  call getarg(4, box_char)
  call getarg(5, dim1_min_char)
  call getarg(6, dim1_max_char)
  call getarg(7, dim1_nbin_char)
  call getarg(8, ngrid_char)
  call getarg(9, nthreads_char)
  
  read(box_char, *) boxsize
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(ngrid_char, *) ngrid
  read(nthreads_char, *) nthreads
  
  write(*,*) '-----------------------'
  write(*,*) 'Running mean_radial_velocity_vs_r.exe'
  write(*,*) 'input parameters:'
  write(*,*) ''
  write(*, *) 'data_filename: ', trim(data_filename)
  write(*, *) 'data_filename_2: ', trim(data_filename_2)
  write(*, *) 'boxsize: ', trim(box_char)
  write(*, *) 'output_filename: ', trim(output_filename)
  write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
  write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
  write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
  write(*, *) 'ngrid: ', trim(ngrid_char)
  write(*, *) 'nthreads: ', trim(nthreads_char)
  write(*,*) ''

  open(10, file=data_filename, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(tracers(ncols, nrows))
  read(10) tracers
  close(10)
  ntracers = nrows
  write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
  write(*,*) 'pos(min), pos(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))

  open(11, file=data_filename_2, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(centres(ncols, nrows))
  read(11) centres
  close(11)
  ncentres = nrows
  if (ncols .ge. 6) then
      has_velocity = .true.
      write(*,*) 'Centres file includes velocity. Pairwise velocity will be computed.'
  else
      write(*,*) 'Centres file is missing velocity. Centre motion will be ignored.'
  end if
  write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)

  allocate(rbin(dim1_nbin))
  allocate(rbin_edges(dim1_nbin + 1))
  allocate(DD(dim1_nbin))
  allocate(VV(dim1_nbin))
  allocate(mean_vr(dim1_nbin))
  allocate(DD_i(ncentres, dim1_nbin))
  allocate(VV_i(ncentres, dim1_nbin))
  
  rwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    rbin_edges(i) = dim1_min+(i-1)*rwidth
  end do
  do i = 1, dim1_nbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
  end do
  
  ! Construct linked list for tracers
  write(*,*) ''
  write(*,*) 'Constructing linked list...'
  allocate(lirst(ngrid, ngrid, ngrid))
  allocate(nlirst(ngrid, ngrid, ngrid))
  allocate(ll(ntracers))
  rgrid = (boxsize) / real(ngrid)
  
  lirst = 0
  ll = 0
  
  do i = 1, ntracers
    indx = int((tracers(1, i)) / rgrid + 1.)
    indy = int((tracers(2, i)) / rgrid + 1.)
    indz = int((tracers(3, i)) / rgrid + 1.)
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)lirst(indx,indy,indz)=i
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)nlirst(indx,indy,indz) = &
    nlirst(indx, indy, indz) + 1
  end do
  
  do i = 1, ntracers
    indx = int((tracers(1, i))/ rgrid + 1.)
    indy = int((tracers(2, i))/ rgrid + 1.)
    indz = int((tracers(3, i))/ rgrid + 1.)
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    &indz.gt.0.and.indz.le.ngrid) then
      ll(lirst(indx,indy,indz)) = i
      lirst(indx,indy,indz) = i
    endif
  end do
  
  write(*,*) 'Linked list successfully constructed'
  write(*,*) ''
  write(*,*) 'Starting loop over tracers...'
  
  DD = 0
  VV = 0
  DD_i = 0
  VV_i = 0
  mean_vr = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  ndif = int(dim1_max / rgrid + 1.)

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
  !$OMP ix, iy, iz, ix2, iy2, iz2, disx, disy, disz, dis, dis2, rind, &
  !$OMP vel, vr, velx, vely, velz)
  do i = 1, ncentres
    ipx = int(centres(1, i) / rgrid + 1.)
    ipy = int(centres(2, i) / rgrid + 1.)
    ipz = int(centres(3, i) / rgrid + 1.)

    do ix = ipx - ndif, ipx + ndif
      do iy = ipy - ndif, ipy + ndif
        do iz = ipz - ndif, ipz + ndif
          if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif + 1)**2) cycle
  
          ix2 = ix
          iy2 = iy
          iz2 = iz
  
          if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
          if (ix2 .lt. 1) ix2 = ix2 + ngrid
          if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
          if (iy2 .lt. 1) iy2 = iy2 + ngrid
          if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
          if (iz2 .lt. 1) iz2 = iz2 + ngrid
  
          ii = lirst(ix2,iy2,iz2)
          if(ii.ne.0) then
            do
              ii = ll(ii)
              disx = tracers(1, ii) - centres(1, i)
              disy = tracers(2, ii) - centres(2, i)
              disz = tracers(3, ii) - centres(3, i)

              if (disx .lt. -boxsize/2) disx = disx + boxsize
              if (disx .gt. boxsize/2) disx = disx - boxsize
              if (disy .lt. -boxsize/2) disy = disy + boxsize
              if (disy .gt. boxsize/2) disy = disy - boxsize
              if (disz .lt. -boxsize/2) disz = disz + boxsize
              if (disz .gt. boxsize/2) disz = disz - boxsize

              dis2 = disx ** 2 + disy **2 + disz ** 2

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)

                if (has_velocity) then
                  velx = tracers(4, ii) - centres(4, i)
                  vely = tracers(5, ii) - centres(5, i)
                  velz = tracers(6, ii) - centres(6, i)
                else
                  velx = tracers(4, ii)
                  vely = tracers(5, ii)
                  velz = tracers(6, ii)
                end if

                vr = (velx * disx + vely * disy + velz * disz) / dis

                rind = int((dis - dim1_min) / rwidth + 1)
                DD_i(i, rind) = DD_i(i, rind) + 1
                VV_i(i, rind) = VV_i(i, rind) + vr
              end if
  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  do i = 1, dim1_nbin
    DD(i) = SUM(DD_i(:, i))
    VV(i) = SUM(VV_i(:, i))
    mean_vr(i) = VV(i) / DD(i)
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_filename, status='replace')
  do i = 1, dim1_nbin
    write(12, fmt='(2f15.5)') rbin(i), mean_vr(i)
  end do

  end program mean_radial_velocity_vs_r
  