program tpcf
  use procedures
  use OMP_LIB
  implicit none
  
  real*8 :: rgrid, disx, disy, disz, dis, dis2, gridmin, gridmax
  real*8 :: rwidth, dim1_max, dim1_min, dim1_max2, dim1_min2
  
  integer*8 :: ng, nr, dim1_nbin, rind
  integer*8 :: i, ii, ix, iy, iz
  integer*8 :: ngrid, ipx, ipy, ipz, ndif
  integer*8 :: end, beginning, rate
  integer*4 :: nthreads
  
  integer*8, dimension(:, :, :), allocatable :: lirst_data, lirst_random
  integer*8, dimension(:), allocatable :: ll_data, ll_random
  
  real*8, allocatable, dimension(:,:)  :: data, random
  real*8, dimension(:), allocatable :: DD, DR, RR, delta
  real*8, dimension(:), allocatable :: weight_data, weight_random
  real*8, dimension(:), allocatable :: rbin, rbin_edges
  real*8, dimension(:, :), allocatable :: DD_i, DR_i, RR_i

  character(20), external :: str
  character(len=500) :: data_filename, output_filename, random_filename
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char
  character(len=10) :: nthreads_char, gridmin_char, gridmax_char
  character(len=10) :: estimator

  logical :: debug = .true.
  
  if (debug) then
    if (iargc() .lt. 11) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) data_filename'
        write(*,*) '2) random_filename'
        write(*,*) '3) output_filename'
        write(*,*) '4) dim1_min'
        write(*,*) '5) dim1_max'
        write(*,*) '6) dim1_nbin'
        write(*,*) '7) ngrid'
        write(*,*) '8) gridmin'
        write(*,*) '9) gridmax'
        write(*,*) '10) estimator'
        write(*,*) '11) nthreads'
        write(*,*) ''
        stop
      end if
    end if

  call system_clock(beginning, rate)
    
  call get_command_argument(number=1, value=data_filename)
  call get_command_argument(number=2, value=random_filename)
  call get_command_argument(number=3, value=output_filename)
  call get_command_argument(number=4, value=dim1_min_char)
  call get_command_argument(number=5, value=dim1_max_char)
  call get_command_argument(number=6, value=dim1_nbin_char)
  call get_command_argument(number=7, value=ngrid_char)
  call get_command_argument(number=8, value=gridmin_char)
  call get_command_argument(number=9, value=gridmax_char)
  call get_command_argument(number=10, value=estimator)
  call get_command_argument(number=11, value=nthreads_char)
  
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(ngrid_char, *) ngrid
  read(gridmin_char, *) gridmin
  read(gridmax_char, *) gridmax
  read(nthreads_char, *) nthreads

  if (debug) then
    write(*,*) '-----------------------'
    write(*,*) 'Running tpcf.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'data_filename: ', trim(data_filename)
    write(*, *) 'random_filename: ', trim(random_filename)
    write(*, *) 'output_filename: ', trim(output_filename)
    write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
    write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
    write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
    write(*, *) 'gridmin: ', trim(gridmin_char), ' Mpc'
    write(*, *) 'gridmax: ', trim(gridmax_char), ' Mpc'
    write(*, *) 'ngrid: ', trim(ngrid_char)
    write(*, *) 'estimator: ', trim(estimator)
    write(*, *) 'nthreads: ', trim(nthreads_char)
    write(*,*) ''
  end if

  call read_unformatted(data_filename, data, weight_data, ng)
  call read_unformatted(random_filename, random, weight_random, nr)

  if (debug) then
    write(*,*) 'ndata dim: ', size(data, dim=1), size(data, dim=2)
    write(*,*) 'data(min, max) = ', minval(data(:,:)), maxval(data(:,:))
    write(*,*) 'weight_data(min, max) = ', minval(weight_data), maxval(weight_data)
    write(*,*) 'nrandom dim: ', size(random, dim=1), size(random, dim=2)
    write(*,*) 'random(min), random(max) = ', minval(random(:,:)), maxval(random(:,:))
    write(*,*) 'weight_random(min, max) = ', minval(weight_random), maxval(weight_random)
  end if

  ! construct linked lists for data and random
  allocate(ll_data(ng))
  allocate(ll_random(nr))
  allocate(lirst_data(ngrid, ngrid, ngrid))
  allocate(lirst_random(ngrid, ngrid, ngrid))
  call linked_list(data, ngrid, gridmin, gridmax, ll_data, lirst_data, rgrid)
  call linked_list(random, ngrid, gridmin, gridmax, ll_random, lirst_random, rgrid)

  allocate(rbin(dim1_nbin))
  allocate(rbin_edges(dim1_nbin + 1))
  allocate(DD(dim1_nbin))
  allocate(DD_i(ng, dim1_nbin))
  allocate(DR(dim1_nbin))
  allocate(DR_i(ng, dim1_nbin))
  allocate(delta(dim1_nbin))
  if (estimator .eq. 'LS') then
    allocate(RR(dim1_nbin))
    allocate(RR_i(nr, dim1_nbin))
  end if

  rwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    rbin_edges(i) = dim1_min+(i-1)*rwidth
  end do
  do i = 1, dim1_nbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
  end do

  DD = 0
  DD_i = 0
  DR = 0
  DR_i = 0
  delta = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  ndif = int(dim1_max / rgrid + 1.)
  if (estimator .eq. 'LS') then
    RR_i = 0
    RR = 0
  end if

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if

  ! Loop over data data
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
  !$OMP ix, iy, iz, disx, disy, disz, dis, dis2, rind)
  do i = 1, ng

    ipx = int((data(1, i) - gridmin) / rgrid + 1.)
    ipy = int((data(2, i) - gridmin) / rgrid + 1.)
    ipz = int((data(3, i) - gridmin) / rgrid + 1.)

    do ix = ipx - ndif, ipx + ndif, 1
      do iy = ipy - ndif, ipy + ndif, 1
        do iz = ipz - ndif, ipz + ndif, 1 
          if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif+ 1)**2) cycle

          ii = lirst_data(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_data(ii)
              disx = data(1, ii) - data(1, i)
              disy = data(2, ii) - data(2, i)
              disz = data(3, ii) - data(3, i)

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) / rwidth + 1)
                DD_i(i, rind) = DD_i(i, rind) + weight_data(i) * weight_data(ii)
              end if

              if(ii .eq. lirst_data(ix, iy, iz)) exit

            end do
          end if

          ii = lirst_random(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_random(ii)
              disx = random(1, ii) - data(1, i)
              disy = random(2, ii) - data(2, i)
              disz = random(3, ii) - data(3, i)

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) / rwidth + 1)
                DR_i(i, rind) = DR_i(i, rind) + weight_data(i) * weight_random(ii)
              end if

              if (ii .eq. lirst_random(ix, iy, iz)) exit

            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  if (estimator .eq. 'LS') then
    ! Loop over random
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
    !$OMP ix, iy, iz, disx, disy, disz, dis, dis2, rind)
    do i = 1, nr

      ipx = int((random(1, i) - gridmin) / rgrid + 1.)
      ipy = int((random(2, i) - gridmin) / rgrid + 1.)
      ipz = int((random(3, i) - gridmin) / rgrid + 1.)

      do ix = ipx - ndif, ipx + ndif, 1
        do iy = ipy - ndif, ipy + ndif, 1
          do iz = ipz - ndif, ipz + ndif, 1 
            if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif + 1)**2) cycle

            ii = lirst_random(ix, iy, iz)
            if (ii .ne. 0) then
              do
                ii = ll_random(ii)
                disx = random(1, ii) - random(1, i)
                disy = random(2, ii) - random(2, i)
                disz = random(3, ii) - random(3, i)

                dis2 = disx * disx + disy * disy + disz * disz

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  dis = sqrt(dis2)
                  rind = int((dis - dim1_min) / rwidth + 1)
                  RR_i(i, rind) = RR_i(i, rind) + weight_random(i) * weight_random(ii)
                end if

                if(ii .eq. lirst_random(ix, iy, iz)) exit

              end do
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if

  ! Add up pair counts
  do i = 1, dim1_nbin
    DD(i) = SUM(DD_i(:, i))
    DR(i) = SUM(DR_i(:, i))
    if (estimator .eq. 'LS') then
      RR(i) = SUM(RR_i(:, i))
    end if
  end do

  ! Normalize pair counts
  DD = DD * 1. / (SUM(weight_data) ** 2)
  DR = DR * 1. / (SUM(weight_data) * SUM(weight_random))
  if (estimator .eq. 'LS') then
    RR = RR * 1. / (SUM(weight_random) ** 2)
  end if

  ! Calculate density contrast
  if (estimator .eq. 'DP') then
    delta = (DD / DR) - 1
  else if (estimator .eq. 'LS') then
    delta = (DD - 2 * DR + RR) / RR
  else
    write(*,*) 'Estimator for the correlation function was not recognized.'
    stop
  end if

  if (debug) then
    write(*,*) ''
    write(*,*) 'Calculation finished. Writing output...'
  end if

  open(12, file=output_filename, status='replace')
  do i = 1, dim1_nbin
    if (estimator .eq. 'LS') then
      write(12, fmt='(5E15.5)') rbin(i), delta(i), SUM(DD_i(:, i)),&
      & SUM(DR_i(:, i)), SUM(RR_i(:, i))
    else
      write(12, fmt='(4E15.5)') rbin(i), delta(i), SUM(DD_i(:, i)),&
      & SUM(DR_i(:, i))
    end if
  end do

  call system_clock(end)
  if (debug) then
    print *, "elapsed time: ", real(end - beginning) / real(rate)
  end if

  end program tpcf
    
