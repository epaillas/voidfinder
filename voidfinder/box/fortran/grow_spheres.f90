PROGRAM grow_spheres
    use OMP_LIB 
    use procedures
    implicit none

    real*8 :: boxsize, delta, rgrid, mean_density, nden
    real*8 :: px, py, pz, disx, disy, disz, dis, dis2
    real*8 :: rvoid, rwidth, rvoidmax, rvoidmax2, vol
    real*8 :: pi = 4.*atan(1.), box2

    integer*8 :: ng, nc, nv, rind, nrows, ncols
    integer*8 :: id, ierr, process_num, iargc, filenumber
    integer*8 :: i, j, k, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: ipx, ipy, ipz, ndif, ngrid, use_weights
    integer*8, parameter :: nrbin = 1000
    integer*4 :: nthreads

    integer*8, dimension(:,:,:), allocatable :: lirst
    integer*8, dimension(:), allocatable :: ll

    real*8, allocatable, dimension(:,:)  :: tracers, centres, voids
    real*8, allocatable, dimension(:) :: weight_tracers, weight_centres
    real*8, dimension(nrbin) :: rbin, cum_rbin

    character(len=500) :: tracers_filename, centres_filename, output_filename
    character(len=10) :: box_char, rvoidmax_char, delta_char, ngrid_char
    character(len=20) :: nthreads_char, use_weights_char, tracers_fileformat
    character(len=1)  :: creturn = achar(13)


    if (iargc() .ne. 10) then
        write(*,*) 'grow_spheres.exe: some parameters are missing.'
        write(*,*) ''
        write(*,*) '1) tracers_filename'
        write(*,*) '2) centres_filename'
        write(*,*) '3) output_filename'
        write(*,*) '4) boxsize'
        write(*,*) '5) density_threshold'
        write(*,*) '6) rvoidmax'
        write(*,*) '7) ngrid'
        write(*,*) '8) nthreads'
        write(*,*) '9) use_weights'
        write(*,*) '10) tracers_fileformat'
        stop
    end if

    call getarg(1, tracers_filename)
    call getarg(2, centres_filename)
    call getarg(3, output_filename)
    call getarg(4, box_char)
    call getarg(5, delta_char)
    call getarg(6, rvoidmax_char)
    call getarg(7, ngrid_char)
    call getarg(8, nthreads_char)
    call getarg(9, use_weights_char)
    call getarg(10, tracers_fileformat)

    read(box_char, *) boxsize
    read(rvoidmax_char, *) rvoidmax
    read(delta_char, *) delta
    read(ngrid_char, *) ngrid
    read(nthreads_char, *) nthreads
    read(use_weights_char, *) use_weights


    write(*,*) '-----------------------'
    write(*,*) 'Running grow_spheres.exe'
    write(*,*) 'Input parameters:'
    write(*,*) ''
    write(*,*) 'tracers_filename: ', trim(tracers_filename)
    write(*,*) 'centres_filename: ', trim(centres_filename)
    write(*,*) 'output_filename: ', trim(output_filename)
    write(*,*) 'boxsize: ', trim(box_char), ' Mpc/h'
    write(*,*) 'rvoidmax: ', trim(rvoidmax_char), ' Mpc/h'
    write(*,*) 'density_threshold: ', trim(delta_char), ' * mean_density'
    write(*,*) 'ngrid: ', ngrid
    write(*,*) 'nthreads: ', nthreads
    write(*,*) 'use_weights: ', use_weights
    write(*,*) 'tracers_fileformat : ', tracers_fileformat

    ! read tracers file
    if (trim(tracers_fileformat) == 'ascii') then
        if (use_weights == 1) then
            call read_catalogue_type2(tracers_filename, tracers, weight_tracers, ng)
        else
            call read_catalogue_type1(tracers_filename, tracers, weight_tracers, ng)
        end if
    else
        if (use_weights == 1) then
            call read_catalogue_type6(tracers_filename, tracers, weight_tracers, ng)
        else
            call read_catalogue_type5(tracers_filename, tracers, weight_tracers, ng)
        end if
    end if

    ! read centres file
    call read_catalogue_type5(centres_filename, centres, weight_centres, nc)

    ! create linked list
    allocate(ll(ng))
    allocate(lirst(ngrid, ngrid, ngrid))
    call linked_list(tracers, boxsize, ngrid, ll, lirst, rgrid)

    ! define useful variables
    allocate(voids(6, nc))
    mean_density = sum(weight_tracers) * 1./(boxsize ** 3)
    ndif = int(rvoidmax / rgrid + 1)
    rwidth = rvoidmax / nrbin
    rvoidmax2 = rvoidmax ** 2
    box2 = boxsize / 2
    nv = 0
    voids = 0

    ! call to OpenMP
    call OMP_SET_NUM_THREADS(nthreads)
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, rbin, cum_rbin,&
    !$OMP ix, iy, iz, ipx, ipy, ipz, ix2, iy2, iz2, dis2, dis, rvoid, vol,&
    !$OMP nden, ng) REDUCTION(+:nv) 
    do i = 1, nc

        rbin = 0
        cum_rbin = 0

        ipx = int(centres(1, i) / rgrid + 1)
        ipy = int(centres(2, i) / rgrid + 1)
        ipz = int(centres(3, i) / rgrid + 1)


        do ix = ipx - ndif, ipx + ndif, 1
            do iy = ipy - ndif, ipy + ndif, 1
                do iz = ipz - ndif, ipz + ndif, 1

                    if (((ix - ipx) ** 2 + (iy - ipy)** 2 + (iz - ipz) ** 2) .gt. (ndif + 1)**2) cycle

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

                            disx = tracers(1, ii) - centres(1, i)
                            disy = tracers(2, ii) - centres(2, i)
                            disz = tracers(3, ii) - centres(3, i)

                            if (disx .lt. -box2) disx = disx + boxsize
                            if (disx .gt. box2) disx = disx - boxsize
                            if (disy .lt. -box2) disy = disy + boxsize
                            if (disy .gt. box2) disy = disy - boxsize
                            if (disz .lt. -box2) disz = disz + boxsize
                            if (disz .gt. box2) disz = disz - boxsize

                            dis2 = disx ** 2 + disy ** 2 + disz ** 2

                            if (dis2 .lt. rvoidmax2) then
                                dis = sqrt(dis2)
                                rind = int(dis / rwidth + 1)
                                rbin(rind) = rbin(rind) + weight_tracers(ii)
                            end if

                            if (ii .eq. lirst(ix2, iy2, iz2)) exit
                        end do
                    end if
                end do
            end do
        end do

        cum_rbin(1) = rbin(1)
        do ii = 2, nrbin
            cum_rbin(ii) =  cum_rbin(ii - 1) + rbin(ii)
        end do

        do ii = nrbin, 1, -1
            rvoid = rwidth * ii
            vol = 4./3 * pi * rvoid ** 3
            nden = cum_rbin(ii) / vol
            ng = int(cum_rbin(ii))
            if (nden .lt. delta * mean_density .and. ng .gt. 0) then
                nv = nv + 1
                voids(1, i) = centres(1, i)
                voids(2, i) = centres(2, i)
                voids(3, i) = centres(3, i)
                voids(4, i) = rvoid
                voids(5, i) = ng
                voids(6, i) = nden / mean_density
                exit
            end if
        end do
    end do
    !$OMP END PARALLEL DO

    deallocate(tracers)
    deallocate(centres)

    open(12, file=output_filename, status='unknown')
    do i = 1, nc
        if (voids(4, i) .ne. 0) then
            write(12, '(4F10.3, 1I10, 1F10.3)') voids(1, i), voids(2, i),&
            voids(3, i), voids(4, i), int(voids(5, i)), voids(6, i)
        end if
    end do

    write(*,*) ''
    write(*,*) nv, ' voids found.'


end PROGRAM grow_spheres
