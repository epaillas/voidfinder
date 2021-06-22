program tpcf
  implicit none
  
  real*8 :: rgrid_x, rgrid_y, rgrid_z, vol, rhomean
  real*8 :: boxsize, boxsize_x, boxsize_y, boxsize_z
  real*8 :: disx, disy, disz, dis, dis2
  real*8 :: rwidth, dim1_max, dim1_min, dim1_max2, dim1_min2
  real*8 :: pi = 4.*atan(1.)
  real*8 :: qperp, qpara
  
  integer*8 :: ntracers, ncentres, dim1_nbin, rind
  integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*8 :: indx, indy, indz, nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif_x, ndif_y, ndif_z
  integer*8 :: ngrid
  integer*8 :: end, beginning, rate
  integer*4 :: argstat1, argstat2
  
  integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, allocatable, dimension(:,:)  :: tracers, centres
  real*8, dimension(:), allocatable :: DD, delta
  real*8, dimension(:), allocatable :: rbin, rbin_edges

  character(20), external :: str
  character(len=500) :: data_filename, data_filename_2, output_filename
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char, box_char
  character(len=10) :: qperp_char, qpara_char
  
  if (iargc() .lt. 8) then
      write(*,*) 'Some arguments are missing.'
      write(*,*) '1) data_filename'
      write(*,*) '2) data_filename_2'
      write(*,*) '3) output_filename'
      write(*,*) '4) boxsize'
      write(*,*) '5) dim1_min'
      write(*,*) '6) dim1_max'
      write(*,*) '7) dim1_nbin'
      write(*,*) '8) ngrid'
      write(*,*) '9) qperp (optional)'
      write(*,*) '10) qpara (optional)'
      write(*,*) ''
      stop
    end if

    call system_clock(beginning, rate)
    
  call get_command_argument(number=1, value=data_filename)
  call get_command_argument(number=2, value=data_filename_2)
  call get_command_argument(number=3, value=output_filename)
  call get_command_argument(number=4, value=box_char)
  call get_command_argument(number=5, value=dim1_min_char)
  call get_command_argument(number=6, value=dim1_max_char)
  call get_command_argument(number=7, value=dim1_nbin_char)
  call get_command_argument(number=8, value=ngrid_char)
  call get_command_argument(number=9, value=qperp_char, status=argstat1)
  call get_command_argument(number=10, value=qpara_char, status=argstat2)
  
  read(box_char, *) boxsize
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(ngrid_char, *) ngrid

  if (argstat1 == 0 .and. argstat2 == 0) then
    read(qperp_char, *) qperp
    read(qpara_char, *) qpara
  else
    qperp = 1.0
    qpara = 1.0
  end if
  
  write(*,*) '-----------------------'
  write(*,*) 'Running tpcf.exe'
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
  write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)

  ! Account for potential geometrical distortions
  boxsize_x = boxsize / qperp
  boxsize_y = boxsize / qperp
  boxsize_z = boxsize / qpara

  tracers(1,:) = tracers(1,:) / qperp
  tracers(2,:) = tracers(2,:) / qperp
  tracers(3,:) = tracers(3,:) / qpara

  centres(1,:) = centres(1,:) / qperp
  centres(2,:) = centres(2,:) / qperp
  centres(3,:) = centres(3,:) / qpara


  if (qperp .ne. 1.0 .or. qpara .ne. 1.0) then
    write(*,*) 'Positions have been shifted due to geometrical distortions'
    write(*,*) 'qperp, qpara: ', qperp, qpara
    write(*,*) 'boxsize_x: ', boxsize_x
    write(*,*) 'boxsize_y: ', boxsize_y
    write(*,*) 'boxsize_z: ', boxsize_z
    write(*,*) 'tracers(min), tracers(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))
    write(*,*) 'centres(min), centres(max) = ', minval(centres(1,:)), maxval(centres(1,:))
  end if

  allocate(rbin(dim1_nbin))
  allocate(rbin_edges(dim1_nbin + 1))
  allocate(DD(dim1_nbin))
  allocate(delta(dim1_nbin))
  
  rwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    rbin_edges(i) = dim1_min+(i-1)*rwidth
  end do
  do i = 1, dim1_nbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
  end do
  
  ! Mean density inside the box
  rhomean = ntracers / (boxsize_x * boxsize_y * boxsize_z)
  
  ! Construct linked list for tracers
  write(*,*) ''
  write(*,*) 'Constructing linked list...'
  allocate(lirst(ngrid, ngrid, ngrid))
  allocate(nlirst(ngrid, ngrid, ngrid))
  allocate(ll(ntracers))
  rgrid_x = (boxsize_x) / real(ngrid)
  rgrid_y = (boxsize_y) / real(ngrid)
  rgrid_z = (boxsize_z) / real(ngrid)
  
  lirst = 0
  ll = 0
  
  do i = 1, ntracers
    indx = int((tracers(1, i)) / rgrid_x + 1.)
    indy = int((tracers(2, i)) / rgrid_y + 1.)
    indz = int((tracers(3, i)) / rgrid_z + 1.)
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)lirst(indx,indy,indz)=i
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)nlirst(indx,indy,indz) = &
    nlirst(indx, indy, indz) + 1
  end do
  
  do i = 1, ntracers
    indx = int((tracers(1, i))/ rgrid_x + 1.)
    indy = int((tracers(2, i))/ rgrid_y + 1.)
    indz = int((tracers(3, i))/ rgrid_z + 1.)
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
  delta = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  
  do i = 1, ncentres
    ipx = int(centres(1, i) / rgrid_x + 1.)
    ipy = int(centres(2, i) / rgrid_y + 1.)
    ipz = int(centres(3, i) / rgrid_z + 1.)

    ndif_x = int(dim1_max / rgrid_x + 1.)
    ndif_y = int(dim1_max / rgrid_y + 1.)
    ndif_z = int(dim1_max / rgrid_z + 1.)
  
    do ix = ipx - ndif_x, ipx + ndif_x
      do iy = ipy - ndif_y, ipy + ndif_y
        do iz = ipz - ndif_z, ipz + ndif_z
          !if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif+ 1)**2) cycle
  
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

              if (disx .lt. -boxsize_x/2) disx = disx + boxsize_x
              if (disx .gt. boxsize_x/2) disx = disx - boxsize_x
              if (disy .lt. -boxsize_y/2) disy = disy + boxsize_y
              if (disy .gt. boxsize_y/2) disy = disy - boxsize_y
              if (disz .lt. -boxsize_z/2) disz = disz + boxsize_z
              if (disz .gt. boxsize_z/2) disz = disz - boxsize_z

              dis2 = disx ** 2 + disy ** 2 + disz ** 2

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) / rwidth + 1)
                DD(rind) = DD(rind) + 1
              end if
  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do

  do i = 1, dim1_nbin
    vol = 4./3 * pi * (rbin_edges(i + 1) ** 3 - rbin_edges(i) ** 3)
    delta(i) = DD(i) / (vol * rhomean * ncentres) - 1
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_filename, status='replace')
  do i = 1, dim1_nbin
    write(12, fmt='(2f15.5)') rbin(i), delta(i)
  end do

  call system_clock(end)
  print *, "elapsed time: ", real(end - beginning) / real(rate)

  end program tpcf
  