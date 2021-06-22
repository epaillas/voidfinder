program tpccf_monopole_2d
  use procedures
  use OMP_LIB
  implicit none
  
  real*8 :: rgrid, disx, disy, dis, dis2, gridmin, gridmax
  real*8 :: rwidth, dim1_max, dim1_min, dim1_max2, dim1_min2
  real*8 :: irwidth
  
  integer*8 :: ng1, ng2, nr1, nr2, rind
  integer*8 :: i, ii, ix, iy
  integer*8 :: dim1_nbin
  integer*8 :: ngrid, ipx, ipy, ndif
  integer*8 :: end, beginning, rate
  integer*4 :: nthreads

  integer*8, dimension(:, :), allocatable :: lirst_data2, lirst_random2
  integer*8, dimension(:), allocatable :: ll_data2, ll_random2
  
  real*8, allocatable, dimension(:,:)  :: data1, data2, random1, random2
  real*8, dimension(:), allocatable :: D1D2, D1R2, R1D2, R1R2, xi_r
  real*8, dimension(:), allocatable :: weight_data1, weight_data2
  real*8, dimension(:), allocatable :: weight_random1, weight_random2
  real*8, dimension(:), allocatable :: rbin, rbin_edges

  character(20), external :: str
  character(len=500) :: data_filename1, data_filename2, output_filename
  character(len=500) :: random_filename1, random_filename2
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char
  character(len=10) :: nthreads_char, gridmin_char, gridmax_char, ngrid_char
  character(len=2) :: estimator

  logical :: debug = .true.
  
  if (debug) then
    if (iargc() .lt. 13) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1)  data_filename1'
        write(*,*) '2)  data_filename2'
        write(*,*) '3)  random_filename1'
        write(*,*) '4)  random_filename2'
        write(*,*) '5)  output_filename'
        write(*,*) '6)  dim1_min'
        write(*,*) '7)  dim1_max'
        write(*,*) '8)  dim1_nbin'
        write(*,*) '9) ngrid'
        write(*,*) '10) gridmin'
        write(*,*) '11) gridmax'
        write(*,*) '12) estimator'
        write(*,*) '13) nthreads'
        write(*,*) ''
        stop
      end if
    end if

  call system_clock(beginning, rate)
    
  call get_command_argument(number=1, value=data_filename1)
  call get_command_argument(number=2, value=data_filename2)
  call get_command_argument(number=3, value=random_filename1)
  call get_command_argument(number=4, value=random_filename2)
  call get_command_argument(number=5, value=output_filename)
  call get_command_argument(number=6, value=dim1_min_char)
  call get_command_argument(number=7, value=dim1_max_char)
  call get_command_argument(number=8, value=dim1_nbin_char)
  call get_command_argument(number=9, value=ngrid_char)
  call get_command_argument(number=10, value=gridmin_char)
  call get_command_argument(number=11, value=gridmax_char)
  call get_command_argument(number=12, value=estimator)
  call get_command_argument(number=13, value=nthreads_char)
  
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(ngrid_char, *) ngrid
  read(gridmin_char, *) gridmin
  read(gridmax_char, *) gridmax
  read(nthreads_char, *) nthreads

  if (debug) then
    write(*,*) '-----------------------'
    write(*,*) 'Running tpccf_monopole_2d.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'data_filename1: ', trim(data_filename1)
    write(*, *) 'data_filename2: ', trim(data_filename2)
    write(*, *) 'random_filename1: ', trim(random_filename1)
    write(*, *) 'random_filename2: ', trim(random_filename2)
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

  call read_unformatted(data_filename1, data1, weight_data1, ng1)
  call read_unformatted(data_filename2, data2, weight_data2, ng2)
  call read_unformatted(random_filename2, random2, weight_random2, nr2)
  if (estimator .eq. 'LS') then
    call read_unformatted(random_filename1, random1, weight_random1, nr1)
  end if

  ! construct linked lists for data and random
  allocate(ll_data2(ng2))
  allocate(ll_random2(nr2))
  allocate(lirst_data2(ngrid, ngrid))
  allocate(lirst_random2(ngrid, ngrid))
  call linked_list_2d(data2, ngrid, gridmin, gridmax, ll_data2, lirst_data2, rgrid)
  call linked_list_2d(random2, ngrid, gridmin, gridmax, ll_random2, lirst_random2, rgrid)

  allocate(D1D2(dim1_nbin))
  allocate(D1R2(dim1_nbin))
  allocate(xi_r(dim1_nbin))
  if (estimator .eq. 'LS') then
    allocate(R1D2(dim1_nbin))
    allocate(R1R2(dim1_nbin))
  end if

  call binning(dim1_min, dim1_max, dim1_nbin, rbin, rbin_edges, rwidth)

  D1D2 = 0
  D1R2 = 0
  xi_r = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  ndif = int(dim1_max / rgrid + 1.)
  irwidth = 1 / rwidth
   if (estimator .eq. 'LS') then
    R1D2 = 0
    R1R2 = 0
  end if

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if

  ! Loop over data 1
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, &
  !$OMP ix, iy, disx, disy, dis, dis2, rind) &
  !$OMP REDUCTION(+:D1D2, D1R2)
  do i = 1, ng1

    ipx = int((data1(1, i) - gridmin) / rgrid + 1.)
    ipy = int((data1(2, i) - gridmin) / rgrid + 1.)

    do ix = ipx - ndif, ipx + ndif, 1
      do iy = ipy - ndif, ipy + ndif, 1
          if ((ix - ipx)**2 + (iy - ipy)**2 .gt. (ndif+ 1)**2) cycle

          ii = lirst_data2(ix, iy)
          if (ii .ne. 0) then
            do
              ii = ll_data2(ii)
              disx = data2(1, ii) - data1(1, i)
              disy = data2(2, ii) - data1(2, i)

              dis2 = disx * disx + disy * disy

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) * irwidth + 1)
                D1D2(rind) = D1D2(rind) + weight_data1(i) * weight_data2(ii)
              end if

              if(ii .eq. lirst_data2(ix, iy)) exit

            end do
          end if

          ii = lirst_random2(ix, iy)
          if (ii .ne. 0) then
            do
              ii = ll_random2(ii)
              disx = random2(1, ii) - data1(1, i)
              disy = random2(2, ii) - data1(2, i)

              dis2 = disx * disx + disy * disy

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) * irwidth + 1)
                D1R2(rind) = D1R2(rind) + weight_data1(i) * weight_random2(ii)
              end if

              if (ii .eq. lirst_random2(ix, iy)) exit

            end do
          end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  if (estimator .eq. 'LS') then
    ! Loop over randoms # 1
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, &
    !$OMP ix, iy, disx, disy, dis, dis2, rind) &
    !$OMP REDUCTION(+:R1R2, R1D2)
    do i = 1, nr1

      ipx = int((random1(1, i) - gridmin) / rgrid + 1.)
      ipy = int((random1(2, i) - gridmin) / rgrid + 1.)

      do ix = ipx - ndif, ipx + ndif, 1
        do iy = ipy - ndif, ipy + ndif, 1
            if ((ix - ipx)**2 + (iy - ipy)**2 .gt. (ndif + 1)**2) cycle

            ii = lirst_random2(ix, iy)
            if (ii .ne. 0) then
              do
                ii = ll_random2(ii)
                disx = random2(1, ii) - random1(1, i)
                disy = random2(2, ii) - random1(2, i)

                dis2 = disx * disx + disy * disy

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  dis = sqrt(dis2)
                  rind = int((dis - dim1_min) * irwidth + 1)
                  R1R2(rind) = R1R2(rind) + weight_random1(i) * weight_random2(ii)
                end if

                  if(ii .eq. lirst_random2(ix, iy)) exit

              end do
            end if

            ii = lirst_data2(ix, iy)
            if (ii .ne. 0) then
              do
                ii = ll_data2(ii)
                disx = data2(1, ii) - random1(1, i)
                disy = data2(2, ii) - random1(2, i)

                dis2 = disx * disx + disy * disy

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  dis = sqrt(dis2)
                  rind = int((dis - dim1_min) * irwidth + 1)
                  R1D2(rind) = R1D2(rind) + weight_random1(i) * weight_data2(ii)
                end if

                  if(ii .eq. lirst_random2(ix, iy)) exit

              end do
            end if


        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if

  ! Normalize pair counts
  D1D2 = D1D2 * 1. / (SUM(weight_data1) * SUM(weight_data2))
  D1R2 = D1R2 * 1. / (SUM(weight_data1) * SUM(weight_random2))
  if (estimator .eq. 'LS') then
    R1R2 = R1R2 * 1. / (SUM(weight_random1) * SUM(weight_random2))
    R1D2 = R1D2 * 1. / (SUM(weight_random1) * SUM(weight_data2))
  end if

  ! Calculate density contrast
  if (estimator .eq. 'DP') then
    xi_r = (D1D2 / D1R2) - 1
  else if (estimator .eq. 'LS') then
    xi_r = (D1D2 - D1R2 - R1D2 + R1R2) / R1R2
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
        write(12, fmt='(7E15.5)') rbin(i), xi_r(i), D1D2(i),&
        & D1R2(i), R1D2(i), R1R2(i)
      else
        write(12, fmt='(5E15.5)') rbin(i), xi_r(i), D1D2(i),&
        & D1R2(i)
      end if
  end do

  call system_clock(end)
  if (debug) then
    print *, "elapsed time: ", real(end - beginning) / real(rate)
  end if

  end program tpccf_monopole_2d
    
