PROGRAM recentering
    use OMP_LIB
    use procedures
    implicit none

    real*8 :: boxsize, delta, rgrid, mean_density, nden
    real*8 :: px, py, pz, disx, disy, disz, dis, dis2
    real*8 :: pxr, pyr, pzr, rvoidr, box2
    real*8 :: rvoid, rwidth, rvoidmax, rvoidmax2
    real*8 :: rnd, rnd_phi, rnd_theta, rnd_rvoid, rnd_vol
    real*8 :: rnd_px, rnd_py, rnd_pz, rnd_ng, rnd_nden
    real*8 :: pi = 4.*atan(1.)

    integer*8 :: ng, nc, nv, rind, stuck, nrows, ncols
    integer*8 :: id, ierr, process_num, iargc, filenumber
    integer*8 :: i, j, k, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: ipx, ipy, ipz, ndif, ngrid, use_weights
    integer*8, parameter ::  nrbin = 1000, nrc = 128
    integer*4 :: nthreads

    integer*8, dimension(:,:,:), allocatable :: lirst
    integer*8, dimension(:), allocatable :: ll

    real*8, allocatable, dimension(:,:)  :: tracers, voids
    real*8, allocatable, dimension(:) :: weight_tracers
    real*8, dimension(nrbin) :: rbin, cum_rbin

    character(len=500) :: input_tracers, input_centres, output_voids
    character(len=10) :: box_char, rvoidmax_char, delta_char, ngrid_char
    character(len=10) :: nthreads_char, use_weights_char
    character(len=1)  :: creturn = achar(13)

    if (iargc() .ne. 9) then
    write(*,*) 'recentering.exe: some parameters are missing.'
    write(*,*) ''
    write(*,*) '1) input_tracers'
    write(*,*) '2) input_centres'
    write(*,*) '3) output_voids'
    write(*,*) '4) boxsize'
    write(*,*) '5) density_threshold'
    write(*,*) '6) rvoidmax'
    write(*,*) '7) ngrid'
    write(*,*) '8) nthreads'
    write(*,*) '9) use_weights'
    stop
    end if

    call getarg(1, input_tracers)
    call getarg(2, input_centres)
    call getarg(3, output_voids)
    call getarg(4, box_char)
    call getarg(5, delta_char)
    call getarg(6, rvoidmax_char)
    call getarg(7, ngrid_char)
    call getarg(8, nthreads_char)
    call getarg(9, use_weights_char)

    read(box_char, *) boxsize
    read(rvoidmax_char, *) rvoidmax
    read(delta_char, *) delta
    read(ngrid_char, *) ngrid
    read(nthreads_char, *) nthreads
    read(use_weights_char, *) use_weights


    write(*,*) '-----------------------'
    write(*,*) 'Running recentering.exe'
    write(*,*) 'Input parameters:'
    write(*,*) ''
    write(*,*) 'input_tracers: ', trim(input_tracers)
    write(*,*) 'input_centres: ', trim(input_centres)
    write(*,*) 'output_voids: ', trim(output_voids)
    write(*,*) 'boxsize: ', trim(box_char), ' Mpc/h'
    write(*,*) 'rvoidmax: ', trim(rvoidmax_char), ' Mpc/h'
    write(*,*) 'density_threshold: ', trim(delta_char), ' * mean_density'
    write(*,*) 'random_centres: ', nrc
    write(*,*) 'ngrid: ', ngrid, ' Mpc/h'
    write(*,*) ''

    if (use_weights == 1) then
        call read_catalogue_type6(input_tracers, tracers, weight_tracers, ng)
    else
        call read_catalogue_type5(input_tracers, tracers, weight_tracers, ng)
    end if 

    call count_rows_ascii(input_centres, nc)

    allocate(voids(6, nc))
    allocate(ll(ng))
    allocate(lirst(ngrid, ngrid, ngrid))
    call linked_list(tracers, boxsize, ngrid, ll, lirst, rgrid)

    mean_density = ng * 1./(boxsize ** 3)
    rgrid = boxsize / ngrid
    ndif = int(rvoidmax / rgrid + 1.)
    rwidth = rvoidmax / nrbin
    rvoidmax2 = rvoidmax ** 2
    box2 = boxsize / 2
    nv = 0
    voids = 0
    open(11, file=input_centres, status='old')

    call OMP_SET_NUM_THREADS(nthreads)
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, rbin, cum_rbin, &
  !$OMP ix, iy, iz, ipx, ipy, ipz, ix2, iy2, iz2, dis2, dis, rvoid, &
  !$OMP rnd_vol, nden, ng, rnd_px, rnd_py, rnd_pz, pxr, pyr, pzr, rvoidr, &
  !$OMP rnd_phi, rnd_theta, px, py, pz, stuck) REDUCTION(+:nv)
  do i = 1, nc
    read(11,*) px, py, pz, rvoid, ng, nden

    pxr = px
    pyr = py
    pzr = pz
    rvoidr = rvoid

      stuck = 0
      do j = 1, nrc
        rbin = 0
        cum_rbin = 0

        call random_number(rnd)
        rnd_phi = rnd * 2 * pi
        call random_number(rnd)
        rnd_theta = acos(rnd * 2 - 1)

        rnd_px = rvoid/4 * sin(rnd_theta) * cos(rnd_phi) + px
        rnd_py = rvoid/4 * sin(rnd_theta) * sin(rnd_phi) + py
        rnd_pz = rvoid/4 * cos(rnd_theta) + pz

        if (sqrt((rnd_px-pxr)**2 + (rnd_py-pyr)**2 + (rnd_pz-pzr)**2 )&
        & .gt. rvoidr) cycle

        if (rnd_px .lt. 0) rnd_px = rnd_px + boxsize
        if (rnd_px .gt. boxsize) rnd_px = rnd_px - boxsize
        if (rnd_py .lt. 0) rnd_py = rnd_py + boxsize
        if (rnd_py .gt. boxsize) rnd_py = rnd_py - boxsize
        if (rnd_pz .lt. 0) rnd_pz = rnd_pz + boxsize
        if (rnd_pz .gt. boxsize) rnd_pz = rnd_pz - boxsize

        ipx = int(rnd_px / rgrid + 1.)
        ipy = int(rnd_py / rgrid + 1.)
        ipz = int(rnd_pz / rgrid + 1.)

        do ix = ipx - ndif, ipx + ndif, 1
          do iy = ipy - ndif, ipy + ndif, 1
            do iz = ipz - ndif, ipz + ndif, 1

              if (real((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2) .gt. (ndif + 1)**2) cycle

              ix2 = ix
              iy2 = iy
              iz2 = iz

              if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
              if (ix2 .lt. 1) ix2 = ix2 + ngrid
              if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
              if (iy2 .lt. 1) iy2 = iy2 + ngrid
              if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
              if (iz2 .lt. 1) iz2 = iz2 + ngrid

              ii = lirst(ix2, iy2, iz2)
              if (ii .ne. 0) then
                do
                  ii = ll(ii)

                  disx = tracers(1, ii) - rnd_px
                  disy = tracers(2, ii) - rnd_py
                  disz = tracers(3, ii) - rnd_pz

                  if (disx .lt. -box2) disx = disx + boxsize
                  if (disx .gt. box2) disx = disx - boxsize
                  if (disy .lt. -box2) disy = disy + boxsize
                  if (disy .gt. box2) disy = disy - boxsize
                  if (disz .lt. -box2) disz = disz + boxsize
                  if (disz .gt. box2) disz = disz - boxsize

                  dis2 = disx * disx + disy * disy + disz * disz

                  if (dis2 .lt. rvoidmax2) then
                    dis = sqrt(dis2)
                    rind = int(dis / rwidth + 1)
                    rbin(rind) = rbin(rind) + 1
                  end if

                  if (ii .eq. lirst(ix2, iy2, iz2)) exit
                end do
              end if

            end do
          end do
        end do

        stuck = stuck + 1

        cum_rbin(1) = rbin(1)
        do ii = 2, nrbin
          cum_rbin(ii) =  cum_rbin(ii - 1) + rbin(ii)
        end do

        do ii = nrbin, 1, -1
          rnd_rvoid = rwidth * ii
          rnd_vol = 4./3 * pi * rvoid ** 3
          rnd_ng = int(cum_rbin(ii))
          rnd_nden = cum_rbin(ii) / rnd_vol
          if (rnd_nden .lt. delta * mean_density .and. rnd_rvoid .gt. rvoid&
          & .and. rnd_ng .gt. 0) then
            rvoid = rnd_rvoid
            px = rnd_px
            py = rnd_py
            pz = rnd_pz
            ng = rnd_ng
            nden = rnd_nden / mean_density
            stuck = 0
            exit
          end if
        end do

        if (stuck .gt. 64) exit

      end do

      voids(1, i) = px
      voids(2, i) = py
      voids(3, i) = pz
      voids(4, i) = rvoid
      voids(5, i) = ng
      voids(6, i) = nden / mean_density

  end do
  !$OMP END PARALLEL DO

  deallocate(tracers)

  open(12, file=output_voids, status='unknown')
  do i = 1, nc
      write(12, '(4F10.3, 1I10, 1F10.3)') voids(1, i), voids(2, i),&
      voids(3, i), voids(4, i), int(voids(5, i)), voids(6, i)
  end do

end PROGRAM recentering
