subroutine test()
  print*, 'fcodes test ... Passed'
end subroutine test

subroutine resolution_grid(uc,mode,maxbin,nx,ny,nz,nbin,res_arr,bin_idx,s_grid)
  implicit none
  real*8, parameter :: PI = 3.141592653589793
  integer, intent(in) :: mode, maxbin,nx,ny,nz
  real, dimension(6),intent(in) :: uc
  integer,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: bin_idx
  real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: s_grid
  real, dimension(0:maxbin-1),intent(out) :: res_arr
  integer, intent(out) :: nbin
  ! locals
  integer, dimension(3) :: nxyz
  real       :: low_res,high_res,resol,tmp_val,tmp_min,val,start,finish
  real       :: r(3),s1(3),step(3)
  integer    :: i,j,k,n,xyzmin(3),xyzmax(3),hkl(3),sloc,ibin,mnloc
  logical    :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  if(debug) print*, 'fcodes_fast...'
  bin_idx = -100
  s_grid = 0.0
  n = 0
  r = 0.0; s1 = 0.0
  res_arr = 0.0
  step = 0.0
  mnloc = -100
  xyzmin = 0; xyzmax = 0; hkl = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1), xyzmax(2), 0
  if(debug) print*, 'unit cell = ', uc
  call get_resol(uc,real(xyzmax(1)),0.0,0.0,r(1))
  call get_resol(uc,0.0,real(xyzmax(2)),0.0,r(2))
  call get_resol(uc,0.0,0.0,real(xyzmax(3)),r(3))
  if(debug) print*,'a-max, b-max, c-max = ', r
  !
  sloc = minloc(r,1)
  hkl = 0
  do i = 1, 3
     if(sloc == i) hkl(i) = sloc/sloc
  end do

  do i = 0, xyzmax(sloc)-1
     !step = (i + 1.5) * hkl
     step = (i + 2.5) * hkl
     call get_resol(uc,step(1),step(2),step(3),resol)
     if(debug) print*, i,step(1),step(2),step(3),resol
     !print*, i,step(1),step(2),step(3),resol
     print*, i,resol
     res_arr(i) = resol
     nbin = i + 1
  end do
  print*, 'nbin=', nbin
  high_res = res_arr(nbin-1)
  call get_resol(uc,0.0,0.0,0.0,low_res)
  !print*,"Low res=",low_res,"High res=",high_res ,'A'

  print*, 'Creating resolution grid. Please wait...'

  ! Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           call get_resol(uc,real(i),real(j),real(k),resol)
           s_grid(i,j,k) = 1.0/resol
           if(k/=xyzmin(3) .and. j/=xyzmin(2) .and. i/=xyzmin(1))then
              s_grid(-i,-j,-k) = s_grid(i,j,k)
           end if
           if(resol < high_res .or. resol > low_res) cycle
!!$           ! Find the matching bin to resol
!!$           mnloc = minloc(sqrt((res_arr - resol)**2), DIM=1) - 1
!!$           !mnloc  = 1
!!$           !if(mnloc < 0 .or. mnloc > nbin-1) cycle
!!$           bin_idx(i,j,k) = mnloc
!!$           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
!!$           bin_idx(-i,-j,-k) = mnloc
           ! Find the matching bin to resol
           do ibin = 0, nbin - 1
              val = sqrt((res_arr(ibin) - resol)**2)
              if(ibin == 0)then
                 tmp_val = val; tmp_min = val
                 mnloc = ibin 
              else
                 tmp_val = val
                 if(tmp_val < tmp_min)then
                    tmp_min = val
                    mnloc = ibin
                 end if
              end if
           end do
           bin_idx(i,j,k) = mnloc
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           bin_idx(-i,-j,-k) = mnloc
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine resolution_grid


subroutine resolution_grid_from_given_resarr(uc,res_arr,mode,nbin,nx,ny,nz,bin_idx,s_grid)
   implicit none
   integer, intent(in) :: mode,nbin,nx,ny,nz
   real, dimension(6),intent(in) :: uc
   real, dimension(0:nbin-1),intent(in) :: res_arr
   integer,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: bin_idx
   real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: s_grid
   ! locals
   real, dimension(6) :: ucstar
   integer, dimension(3) :: nxyz
   real       :: low_res,high_res,resol,tmp_val,tmp_min,val,start,finish
   integer    :: i,j,k,mnloc,xyzmin(3),xyzmax(3),ibin
   logical    :: debug
   !
   debug         = .FALSE.
   if(mode == 1) debug = .TRUE.
   call cpu_time(start)
   !
   if(debug) print*, 'fcodes_fast...'
   bin_idx = -100
   s_grid = 0.0
   xyzmin = 0; xyzmax = 0
   nxyz = (/ nx, ny, nz /)
   xyzmin(1) = int(-nxyz(1)/2)
   xyzmin(2) = int(-nxyz(2)/2)
   xyzmin(3) = int(-nxyz(3)/2)
   xyzmax    = -(xyzmin+1)
   !
   print*, 'nbin=', nbin
   high_res = res_arr(nbin-1)
   !call get_reciprocal_basis(uc, ucstar)
   call get_resol(uc,0.0,0.0,0.0,low_res)
   !call get_resol2(ucstar,0.0,0.0,0.0,low_res)
   print*,"Low res=",low_res,"High res=",high_res ,'A'
   print*, 'Creating resolution grid. Please wait...'
   ! Friedel's Law
   do i=xyzmin(1), xyzmax(1)
      do j=xyzmin(2), xyzmax(2)
         do k=xyzmin(3), 0 !xyzmax(3)
            call get_resol(uc,real(i),real(j),real(k),resol)
            !call get_resol2(ucstar,real(i),real(j),real(k),resol)
            s_grid(i,j,k) = 1.0/resol
            if(k/=xyzmin(3) .and. j/=xyzmin(2) .and. i/=xyzmin(1))then
               s_grid(-i,-j,-k) = s_grid(i,j,k)
            end if
            if(resol < high_res .or. resol > low_res) cycle
            ! Find the matching bin to resol
            do ibin = 0, nbin - 1
               val = sqrt((res_arr(ibin) - resol)**2)
               if(ibin == 0)then
                  tmp_val = val; tmp_min = val
                  mnloc = ibin 
               else
                  tmp_val = val
                  if(tmp_val < tmp_min)then
                     tmp_min = val
                     mnloc = ibin
                  end if
               end if
            end do
            bin_idx(i,j,k) = mnloc
            if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
            bin_idx(-i,-j,-k) = mnloc
         end do
      end do
   end do
   call cpu_time(finish)
   if(debug) print*, 'time for calculation(s) = ', finish-start
 end subroutine resolution_grid_from_given_resarr


subroutine resol_grid_em(uc,mode,maxbin,nx,ny,nz,nbin,res_arr,bin_idx,s_grid)
  implicit none
  real*8, parameter :: PI = 3.141592653589793
  integer, intent(in) :: mode, maxbin,nx,ny,nz
  real, dimension(6),intent(in) :: uc
  integer,dimension(-nz/2:(nz-2)/2,-ny/2:(ny-2)/2,-nx/2:(nx-2)/2),intent(out) :: bin_idx
  real, dimension(-nz/2:(nz-2)/2,-ny/2:(ny-2)/2,-nx/2:(nx-2)/2),intent(out) :: s_grid
  real, dimension(0:maxbin-1),intent(out) :: res_arr
  integer, intent(out) :: nbin
  ! locals
  integer, dimension(3) :: nxyz
  real, dimension(0:maxbin-1) :: bin_arr
  real       :: low_res,high_res,resol,start,finish
  real       :: r(3),s1(3),step(3)
  integer    :: i,j,k,n,xyzmin(3),xyzmax(3),hkl(3),sloc,ibin,mnloc,ib
  logical    :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  if(debug) print*, 'fcodes_fast...'
  bin_idx = -100
  s_grid = 0.0
  n = 0
  r = 0.0; s1 = 0.0
  res_arr = 0.0
  step = 0.0
  mnloc = -100
  xyzmin = 0; xyzmax = 0; hkl = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1), xyzmax(2), 0
  if(debug) print*, 'unit cell = ', uc
  call get_resol(uc,real(xyzmax(1)),0.0,0.0,r(1))
  call get_resol(uc,0.0,real(xyzmax(2)),0.0,r(2))
  call get_resol(uc,0.0,0.0,real(xyzmax(3)),r(3))
  if(debug) print*,'a-max, b-max, c-max = ', r
  !
  sloc = minloc(r,1)
  hkl = 0
  do i = 1, 3
     if(sloc == i) hkl(i) = sloc/sloc
  end do

  do i = 0, xyzmax(sloc)-1
     !step = (i + 1.5) * hkl
     step = (i + 2.5) * hkl
     call get_resol(uc,step(1),step(2),step(3),resol)
     if(debug) print*, i,step(1),step(2),step(3),resol
     !print*, i,step(1),step(2),step(3),resol
     print*, i,resol
     res_arr(i) = resol
     bin_arr(i) = i + 2.5
     nbin = i + 1
  end do
  print*, 'nbin=', nbin
  high_res = res_arr(nbin-1)
  call get_resol(uc,0.0,0.0,0.0,low_res)
  !print*,"Low res=",low_res,"High res=",high_res ,'A'

  print*, 'Creating resolution grid. Please wait...'

  ! Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        outer: do k=xyzmin(3), 0 !xyzmax(3)
           do ib=0, nbin-1
              if(sqrt(real(k)**2 + real(j)**2 + real(i)**2) .le. bin_arr(ib))then
                 bin_idx(k,j,i) = ib
                 s_grid(k,j,i) = 1.0/res_arr(ib)
                 if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1))then
                    cycle outer
                 end if
                 bin_idx(-k,-j,-i) = ib
                 s_grid(-k,-j,-i) = s_grid(k,j,i)
                 cycle outer
              end if
           end do
        end do outer
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine resol_grid_em

subroutine resolution_grid_full(uc,highres,mode,maxbin,nx,ny,nz,resol_grid,s_grid,mask)
  implicit none
  integer, intent(in) :: mode,maxbin,nx,ny,nz
  real, intent(in) :: highres
  real, dimension(6),intent(in) :: uc
  real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: s_grid
  real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: resol_grid
  integer, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: mask
  ! locals
  integer, dimension(3) :: nxyz
  real       :: resol,start,finish
  integer    :: i,j,k,xyzmin(3),xyzmax(3)
  logical    :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  if(debug) print*, 'fcodes_fast...'
  resol_grid = 0.0
  s_grid = 0.0
  mask = 0
  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax
  if(debug) print*, 'unit cell = ', uc

  print*, 'Creating resolution grid. Please wait...'
  print*, 'High resolution cutoff: ', highres, 'A'

  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           call get_resol(uc,real(i),real(j),real(k),resol)
           if(resol < highres) cycle
           mask(i,j,k) = 1
           resol_grid(i,j,k) = resol
           s_grid(i,j,k) = 1.0/resol
           if(k == xyzmin(3) .or. j==xyzmin(2) .or. i==xyzmin(1)) cycle
           mask(-i,-j,-k) = 1
           resol_grid(-i,-j,-k) = resol_grid(i,j,k)
           s_grid(-i,-j,-k) = s_grid(i,j,k)
        end do
     end do
  end do
  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine resolution_grid_full

subroutine make_resarr(uc,maxbin,res_arr,nbin,firststep)
  implicit none
  integer, intent(in) :: maxbin
  real, dimension(6),intent(in) :: uc
  real, intent(in),optional :: firststep
  real, dimension(0:maxbin-1),intent(out) :: res_arr
  integer, intent(out) :: nbin
  ! locals
  integer :: i
  real :: fstep, step, resol

  res_arr = 0.0
  print*, 'first step: ', firststep
  ! current F2PY do not handle fortran optional arguments properly
  ! optional args are always present
  if(present(firststep))then
     if((firststep <= 0.0).or.(firststep > maxbin))then
        fstep = 2.5
     else
        fstep = firststep
     end if
  end if
  print*, 'fstep: ', fstep
  do i = 0, maxbin-2 ! ignore last line
     step = (i + fstep)
     call get_resol(uc,step,0.0,0.0,resol)
     print*, i,step,resol
     res_arr(i) = resol
     nbin = i + 1
  end do
  return
end subroutine make_resarr

subroutine conv3d_to_1d(f3d,uc,nx,ny,nz,mode,f1d,resol1d)
  ! this subroutine converts 3d grid data into 1d. only the
  ! hemispphereis taken.
  implicit none
  integer, intent(in) :: mode,nx,ny,nz
  real, dimension(6),intent(in) :: uc
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: f3d
  complex*16, dimension(nx*ny*(nz+2)/2),intent(out) :: f1d
  real, dimension(nx*ny*(nz+2)/2),intent(out) :: resol1d
  ! locals
  integer, dimension(nx*ny*(nz+2)/2) :: h, k, l
  real, dimension(nx*ny*(nz+2)/2) :: s2
  integer, dimension(3) :: nxyz
  real       :: astar, bstar, cstar
  real       :: resol,start,finish
  integer    :: i1,i2,i3,j,xyzmin(3),xyzmax(3)
  integer    :: nxh,nxh1,idxzero
  logical    :: debug
  
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1), xyzmax(2), 0
  if(debug) print*, 'unit cell = ', uc

  print*, 'Creating resolution array. Please wait...'

  j = 0
  do i1=xyzmin(1), xyzmax(1)
     do i2=xyzmin(2), xyzmax(2)
        do i3=xyzmin(3), 0
           j = j + 1
           h(j) = i1
           k(j) = i2
           l(j) = i3
           f1d(j) = f3d(i1,i2,i3)
        end do
     end do
  end do
  astar = 1.0 / uc(1)
  bstar = 1.0 / uc(2)
  cstar = 1.0 / uc(3)
  s2 = (h * astar)**2 + (k * bstar)**2 + (l * cstar)**2
  nxh = int(nx/2); nxh1 = nxh+1; idxzero = nxh1 * (nxh * nx + nxh1)
  s2(idxzero) = astar 
  !s2 = s2 + merge(1.0e-8, 0.0, s2 == 0.0)
  resol1d = 1.0 / sqrt(s2)
  call cpu_time(finish)
  print*, 'time for resol calculation: ', finish-start
  return
end subroutine conv3d_to_1d

subroutine calc_halffsc(hf1,hf2,bin_idx,nbin,mode,nx,ny,nz, &
     Fo,Eo,bin_stats,bin_arr_count)
  implicit none
  real*8, parameter :: PI = 3.141592653589793
  integer, intent(in) :: nbin,mode,nx,ny,nz
  complex*16, dimension(-nz/2:(nz-2)/2, -ny/2:(ny-2)/2, -nx/2:(nx-2)/2),intent(in)  :: hf1,hf2
  complex*16, dimension(-nz/2:(nz-2)/2, -ny/2:(ny-2)/2, -nx/2:(nx-2)/2),intent(out) :: Fo,Eo
  integer, dimension(-nz/2:(nz-2)/2, -ny/2:(ny-2)/2, -nx/2:(nx-2)/2),intent(in) :: bin_idx
  real*8, dimension(0:nbin-1,0:3),intent(out) :: bin_stats
  integer, dimension(0:nbin-1),intent(out) :: bin_arr_count
  ! locals
  integer, dimension(3) :: nxyz
  real*8, dimension(0:nbin-1, 0:17) :: bindata
  complex*16  :: fdiff
  real*8  :: A,B,A1,A2,B1,B2
  real :: start,finish
  integer :: i,j,k,xyzmin(3),xyzmax(3),ibin,indx
  logical :: debug,make_all_zero
  !
  debug         = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  !bin_stats[:,0-3] - bin_noise_var,bin_sgnl_var,bin_total_var,bin_fsc
  !bindata[:,0-5] - A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
  !bindata[:,6-9] - A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
  !bindata[:,10-13] - A_sum,B_sum,AA_sum,BB_sum
  !bindata[:,14-16] - F1_var, F2_var, F1F2_covar
  !bindata[:,17] - bin_arr_fdiff

  Fo = dcmplx(0.0d0, 0.0d0)
  Eo = dcmplx(0.0d0, 0.0d0)
  bin_stats = 0.0
  bindata = 0.0
  bin_arr_count = 0
  xyzmin = 0; xyzmax = 0
  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'Use only hemisphere data'
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1),xyzmax(2),0
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0!xyzmax(3)
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           indx = bin_idx(k,j,i)
           if(indx < 0 .or. indx > nbin-1) cycle
           bin_arr_count(indx) = bin_arr_count(indx) + 1
           fdiff = hf1(k,j,i) - hf2(k,j,i)
           bindata(indx,17) = bindata(indx,17) + real(fdiff * conjg(fdiff))
           ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
           A1 = real(hf1(k,j,i));  A2 = real(hf2(k,j,i))
           B1 = aimag(hf1(k,j,i)); B2 =  aimag(hf2(k,j,i))
           Fo(k,j,i) = (hf1(k,j,i) + hf2(k,j,i))/2.0                
           A = real(Fo(k,j,i)); B = aimag(Fo(k,j,i))
           bindata(indx,0)  = bindata(indx,0)  + A1
           bindata(indx,1)  = bindata(indx,1)  + B1
           bindata(indx,2)  = bindata(indx,2)  + A2
           bindata(indx,3)  = bindata(indx,3)  + B2
           bindata(indx,4)  = bindata(indx,4)  + A1*A2
           bindata(indx,5)  = bindata(indx,5)  + B1*B2
           bindata(indx,6)  = bindata(indx,6)  + A1*A1
           bindata(indx,7)  = bindata(indx,7)  + B1*B1
           bindata(indx,8)  = bindata(indx,8)  + A2*A2
           bindata(indx,9)  = bindata(indx,9)  + B2*B2
           bindata(indx,10) = bindata(indx,10) + A
           bindata(indx,11) = bindata(indx,11) + B
           bindata(indx,12) = bindata(indx,12) + A*A
           bindata(indx,13) = bindata(indx,13) + B*B
           Fo(-k,-j,-i) = conjg(Fo(k,j,i))
        end do
     end do
  end do
  if(debug) print*, 'bin_arr_count=', sum(bin_arr_count)
  do ibin=0, nbin-1 !to make compatible with python arrays
     bin_stats(ibin,0) = bindata(ibin,17) / (bin_arr_count(ibin) * 4)
     bin_stats(ibin,2) = (bindata(ibin,12) + bindata(ibin,13))/bin_arr_count(ibin) - &
          ((bindata(ibin,10)/bin_arr_count(ibin))**2 + &
          (bindata(ibin,11)/bin_arr_count(ibin))**2)
     bindata(ibin,14)  = (bindata(ibin,6) + bindata(ibin,7))/bin_arr_count(ibin) - &
          ((bindata(ibin,0)/bin_arr_count(ibin))**2 + &
          (bindata(ibin,1)/bin_arr_count(ibin))**2)
     bindata(ibin,15)  = (bindata(ibin,8) + bindata(ibin,9))/bin_arr_count(ibin) - &
          ((bindata(ibin,2)/bin_arr_count(ibin))**2 + &
          (bindata(ibin,3)/bin_arr_count(ibin))**2)
     bindata(ibin,16)  = (bindata(ibin,4) + bindata(ibin,5)) / bin_arr_count(ibin) - &
          (bindata(ibin,0) / bin_arr_count(ibin) * &
          bindata(ibin,2) / bin_arr_count(ibin) + &
          bindata(ibin,1) / bin_arr_count(ibin) * &
          bindata(ibin,3) / bin_arr_count(ibin))
     bin_stats(ibin,1) = bindata(ibin,16)
     bin_stats(ibin,3) = bin_stats(ibin,1) / (sqrt(bindata(ibin,14)) * sqrt(bindata(ibin,15)))
     if(debug)then
        print*,ibin,bin_stats(ibin,0),bin_stats(ibin,1), &
             bin_stats(ibin,2),bin_stats(ibin,3),bin_arr_count(ibin)
     end if
  end do
  ! Calculate normalized structure factors
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           indx = bin_idx(k,j,i)
           if(indx < 0 .or. indx > nbin-1) cycle
           if(bin_stats(indx,1) <= 0.0) cycle ! singal var cannot be negative
           Eo(k,j,i) = Fo(k,j,i)/sqrt(bin_stats(indx,2))
           Eo(-k,-j,-i) = conjg(Eo(k,j,i))
        end do
     end do
  end do
  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine calc_halffsc

subroutine calc_fsc(hf1,hf2,bin_idx,nbin,mode,binstats,bin_arr_count,nx,ny,nz)
  implicit none
  integer,intent(in) :: nbin,mode,nx,ny,nz
  integer, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: bin_idx
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: hf1,hf2
  real*8, dimension(0:nbin-1,0:1), intent(out) :: binstats
  integer, dimension(0:nbin-1),intent(out) :: bin_arr_count
  real*8, dimension(0:nbin-1,0:11) :: bindata
  real*8    :: A1,A2,B1,B2
  integer   :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,ibin,indx
  real      :: start, finish
  logical   :: debug, make_all_zero 
  !
  debug = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)
  !binstats[:,0:1] - F1F2_covar,bin_fsc
  !bindata[:,0:5] - A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
  !bindata[:,6:9] - A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
  !bindata[:,10:11] - F1_var,F2_var
  !
  binstats = 0.0
  bindata = 0.0
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'

  bin_arr_count = 0
  if(debug) print*, 'using hemisphere data...'
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, 0 !zmax
           indx = bin_idx(k,j,i)
           if(indx < 0 .or. indx > nbin-1)then
              cycle
           else
              bin_arr_count(indx) = bin_arr_count(indx) + 1
              ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
              A1 = real(hf1(k,j,i));  A2 = real(hf2(k,j,i))
              B1 = aimag(hf1(k,j,i)); B2 = aimag(hf2(k,j,i))
              bindata(indx,0) = bindata(indx,0) + A1
              bindata(indx,1) = bindata(indx,1) + B1
              bindata(indx,2) = bindata(indx,2) + A2
              bindata(indx,3) = bindata(indx,3) + B2
              bindata(indx,4) = bindata(indx,4) + A1*A2
              bindata(indx,5) = bindata(indx,5) + B1*B2
              bindata(indx,6) = bindata(indx,6) + A1*A1
              bindata(indx,7) = bindata(indx,7) + B1*B1
              bindata(indx,8) = bindata(indx,8) + A2*A2
              bindata(indx,9) = bindata(indx,9) + B2*B2
           end if
        end do
     end do
  end do
  if(debug) print*,'ibin F1F2_covar(ibin) F1_var(ibin) F2_var(ibin) bin_fsc(ibin) bin_reflex_count'
  do ibin=0, nbin-1 !to make compatible with python arrays
     binstats(ibin,0) = (bindata(ibin,4) + bindata(ibin,5)) / bin_arr_count(ibin) - &
          (bindata(ibin,0) / bin_arr_count(ibin) * &
          bindata(ibin,2) / bin_arr_count(ibin) + &
          bindata(ibin,1) / bin_arr_count(ibin) * &
          bindata(ibin,3) / bin_arr_count(ibin))
     bindata(ibin,10) = (bindata(ibin,6) + bindata(ibin,7))/bin_arr_count(ibin) - &
          ((bindata(ibin,0)/bin_arr_count(ibin))**2 + &
          (bindata(ibin,1)/bin_arr_count(ibin))**2)
     bindata(ibin,11) = (bindata(ibin,8) + bindata(ibin,9))/bin_arr_count(ibin) - &
          ((bindata(ibin,2)/bin_arr_count(ibin))**2 + &
          (bindata(ibin,3)/bin_arr_count(ibin))**2)
     binstats(ibin,1) = binstats(ibin,0) / (sqrt(bindata(ibin,10)) * sqrt(bindata(ibin,11)))
     if(debug)then
        print*,ibin,binstats(ibin,0),bindata(ibin,10),bindata(ibin,11), &
             binstats(ibin,1),bin_arr_count(ibin), &
             bindata(ibin,8), bindata(ibin,9), bindata(ibin,2), bindata(ibin,3)
     end if
  end do
  call cpu_time(finish)
  if(debug) print*, 'time for loop = ', finish-start
  return
end subroutine calc_fsc

subroutine calc_fsc_using_halfmaps(hf1,hf2,bin_idx,nbin,mode,nx,ny,nz, &
     Fo,Eo,bin_noise_var,bin_sgnl_var,bin_total_var,bin_fsc,bin_arr_count)
  implicit none
  real*8,    parameter :: PI = 3.141592653589793

  integer,                intent(in) :: nbin,mode,nx,ny,nz
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: hf1,hf2
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: Fo,Eo
  integer,   dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in) :: bin_idx
  real*8,    dimension(0:nbin-1),intent(out) :: bin_sgnl_var,bin_noise_var,bin_total_var,bin_fsc
  integer,   dimension(0:nbin-1),intent(out) :: bin_arr_count
  ! locals
  integer,   dimension(3)          :: nxyz

  real*8,    dimension(0:nbin-1) :: A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
  real*8,    dimension(0:nbin-1) :: A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
  real*8,    dimension(0:nbin-1) :: bin_arr_fdiff,A_sum,B_sum,AA_sum,BB_sum
  real*8,    dimension(0:nbin-1) :: F1_var, F2_var, F1F2_covar
  !
  complex*16  :: fdiff
  real*8     :: A,B,A1,A2,B1,B2,bin_sigvar,denominator
  real       :: start,finish
  integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin
  logical    :: debug,make_all_zero
  !
  debug         = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  Fo = dcmplx(0.0d0, 0.0d0)
  Eo = dcmplx(0.0d0, 0.0d0)

  bin_arr_fdiff = 0.0
  bin_sigvar    = 0.0
  bin_noise_var = 0.0
  bin_sgnl_var  = 0.0
  bin_total_var = 0.0

  F1F2_covar = 0.0
  F1_var = 0.0
  F2_var = 0.0
  bin_total_var = 0.0
  bin_fsc = 0.0

  A_sum = 0.0
  B_sum = 0.0
  AA_sum = 0.0
  BB_sum = 0.0

  A1_sum = 0.0; A2_sum = 0.0
  B1_sum = 0.0; B2_sum = 0.0
  A1A2_sum = 0.0; B1B2_sum = 0.0
  A1A1_sum = 0.0; B1B1_sum = 0.0
  A2A2_sum = 0.0; B2B2_sum = 0.0

  bin_arr_count = 0
  xyzmin = 0; xyzmax = 0
  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'Use only hemisphere data'
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1),xyzmax(2),0

  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0!xyzmax(3)
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           !New calculation of total-var and noise var
           fdiff = hf1(i,j,k) - hf2(i,j,k)
           bin_arr_fdiff(bin_idx(i,j,k)) = bin_arr_fdiff(bin_idx(i,j,k)) + real(fdiff * conjg(fdiff))
           Fo(i,j,k) = (hf1(i,j,k) + hf2(i,j,k))/2.0                
           A = real(Fo(i,j,k)); B = aimag(Fo(i,j,k))
           A_sum(bin_idx(i,j,k)) = A_sum(bin_idx(i,j,k)) + A
           AA_sum(bin_idx(i,j,k)) = AA_sum(bin_idx(i,j,k)) + A*A
           B_sum(bin_idx(i,j,k)) = B_sum(bin_idx(i,j,k)) + B
           BB_sum(bin_idx(i,j,k)) = BB_sum(bin_idx(i,j,k)) + B*B
           ! end of new calculation

           ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
           A1 = real(hf1(i,j,k));  A2 = real(hf2(i,j,k))
           B1 = aimag(hf1(i,j,k)); B2 = aimag(hf2(i,j,k))
           A1_sum(bin_idx(i,j,k)) = A1_sum(bin_idx(i,j,k)) + A1
           A2_sum(bin_idx(i,j,k)) = A2_sum(bin_idx(i,j,k)) + A2
           B1_sum(bin_idx(i,j,k)) = B1_sum(bin_idx(i,j,k)) + B1
           B2_sum(bin_idx(i,j,k)) = B2_sum(bin_idx(i,j,k)) + B2

           A1A2_sum(bin_idx(i,j,k)) = A1A2_sum(bin_idx(i,j,k)) + A1 * A2
           B1B2_sum(bin_idx(i,j,k)) = B1B2_sum(bin_idx(i,j,k)) + B1 * B2

           A1A1_sum(bin_idx(i,j,k)) = A1A1_sum(bin_idx(i,j,k)) + A1 * A1
           B1B1_sum(bin_idx(i,j,k)) = B1B1_sum(bin_idx(i,j,k)) + B1 * B1

           A2A2_sum(bin_idx(i,j,k)) = A2A2_sum(bin_idx(i,j,k)) + A2 * A2
           B2B2_sum(bin_idx(i,j,k)) = B2B2_sum(bin_idx(i,j,k)) + B2 * B2
           !if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           Fo(-i,-j,-k) = conjg(Fo(i,j,k))
           !if(-k <= xyzmax(3) .or. -j <= xyzmax(2) .or. -i <= xyzmax(1))then   
           !   Fo(-i,-j,-k) = conjg(Fo(i,j,k))
           !end if
        end do
     end do
  end do

  if(debug) print*, 'bin_arr_count=', sum(bin_arr_count)

  do ibin=0, nbin-1 !to make compatible with python arrays
     bin_noise_var(ibin) = bin_arr_fdiff(ibin) / (bin_arr_count(ibin) * 4)
     bin_total_var(ibin) = (AA_sum(ibin) + BB_sum(ibin))/bin_arr_count(ibin) &
          - ((A_sum(ibin)/bin_arr_count(ibin))**2 + (B_sum(ibin)/bin_arr_count(ibin))**2)

     F1F2_covar(ibin) = (A1A2_sum(ibin) + B1B2_sum(ibin)) / bin_arr_count(ibin) - &
          (A1_sum(ibin) / bin_arr_count(ibin) * A2_sum(ibin) / bin_arr_count(ibin) + &
          B1_sum(ibin) / bin_arr_count(ibin) * B2_sum(ibin) / bin_arr_count(ibin))

     F1_var(ibin) = (A1A1_sum(ibin) + B1B1_sum(ibin))/bin_arr_count(ibin) - &
          ((A1_sum(ibin)/bin_arr_count(ibin))**2 + (B1_sum(ibin)/bin_arr_count(ibin))**2)
     F2_var(ibin) = (A2A2_sum(ibin) + B2B2_sum(ibin))/bin_arr_count(ibin) - &
          ((A2_sum(ibin)/bin_arr_count(ibin))**2 + (B2_sum(ibin)/bin_arr_count(ibin))**2)

     bin_sgnl_var(ibin) = F1F2_covar(ibin)
     denominator = (sqrt(F1_var(ibin)) * sqrt(F2_var(ibin)))
     bin_fsc(ibin) = F1F2_covar(ibin) / denominator
     !print*,ibin,bin_total_var(ibin),denominator, denominator-bin_total_var(ibin)
     if(debug)then
        print*,ibin,bin_noise_var(ibin),bin_sgnl_var(ibin), &
             bin_total_var(ibin),bin_fsc(ibin),bin_arr_count(ibin)
     end if
  end do

  ! Calculate normalized structure factors
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           if(bin_sgnl_var(bin_idx(i,j,k)) <= 0.0) cycle ! singal var cannot be negative
           Eo(i,j,k) = Fo(i,j,k)/sqrt(bin_total_var(bin_idx(i,j,k)))
           !if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           Eo(-i,-j,-k) = conjg(Eo(i,j,k))
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine calc_fsc_using_halfmaps

!!$subroutine get_normalized_sf(hf1,hf2,bin_idx,nbin,mode,nx,ny,nz, &
!!$     Eo,bin_noise_var,bin_sgnl_var,bin_total_var,bin_fsc,bin_arr_count)
!!$  implicit none
!!$  real*8, parameter :: PI = 3.141592653589793
!!$
!!$  integer, intent(in) :: nbin,mode,nx,ny,nz
!!$  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: hf1,hf2
!!$  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2, 3),intent(out) :: Eo
!!$  integer, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in) :: bin_idx
!!$  real*8, dimension(0:nbin-1),intent(out) :: bin_sgnl_var,bin_noise_var,bin_total_var,bin_fsc
!!$  integer, dimension(0:nbin-1),intent(out) :: bin_arr_count
!!$  ! locals
!!$  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2) :: Fo
!!$  integer, dimension(3) :: nxyz
!!$  real*8, dimension(0:nbin-1) :: A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
!!$  real*8, dimension(0:nbin-1) :: A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
!!$  real*8, dimension(0:nbin-1) :: bin_arr_fdiff,A_sum,B_sum,AA_sum,BB_sum
!!$  real*8, dimension(0:nbin-1) :: F1_var, F2_var, F1F2_covar
!!$  !
!!$  complex*16  :: fdiff
!!$  real*8     :: A,B,A1,A2,B1,B2,bin_sigvar,denominator
!!$  real       :: start,finish
!!$  integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin
!!$  logical    :: debug,make_all_zero
!!$  !
!!$  debug         = .FALSE.
!!$  make_all_zero = .FALSE.
!!$  if(mode == 1) debug = .TRUE.
!!$  call cpu_time(start)
!!$
!!$  Fo = dcmplx(0.0d0, 0.0d0)
!!$  Eo = dcmplx(0.0d0, 0.0d0)
!!$
!!$  bin_arr_fdiff = 0.0
!!$  bin_sigvar    = 0.0
!!$  bin_noise_var = 0.0
!!$  bin_sgnl_var  = 0.0
!!$  bin_total_var = 0.0
!!$
!!$  F1F2_covar = 0.0
!!$  F1_var = 0.0
!!$  F2_var = 0.0
!!$  bin_total_var = 0.0
!!$  bin_fsc = 0.0
!!$
!!$  A_sum = 0.0
!!$  B_sum = 0.0
!!$  AA_sum = 0.0
!!$  BB_sum = 0.0
!!$
!!$  A1_sum = 0.0; A2_sum = 0.0
!!$  B1_sum = 0.0; B2_sum = 0.0
!!$  A1A2_sum = 0.0; B1B2_sum = 0.0
!!$  A1A1_sum = 0.0; B1B1_sum = 0.0
!!$  A2A2_sum = 0.0; B2B2_sum = 0.0
!!$
!!$  bin_arr_count = 0
!!$  xyzmin = 0; xyzmax = 0
!!$  nxyz = (/ nx, ny, nz /)
!!$
!!$  xyzmin(1) = int(-nxyz(1)/2)
!!$  xyzmin(2) = int(-nxyz(2)/2)
!!$  xyzmin(3) = int(-nxyz(3)/2)
!!$  xyzmax    = -(xyzmin+1)
!!$  if(debug) print*, 'Use only hemisphere data'
!!$  if(debug) print*, 'xyzmin = ', xyzmin
!!$  if(debug) print*, 'xyzmax = ', xyzmax(1),xyzmax(2),0
!!$
!!$  do i=xyzmin(1), xyzmax(1)
!!$     do j=xyzmin(2), xyzmax(2)
!!$        do k=xyzmin(3), 0
!!$           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
!!$           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
!!$           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
!!$           fdiff = hf1(i,j,k) - hf2(i,j,k)
!!$           bin_arr_fdiff(bin_idx(i,j,k)) = bin_arr_fdiff(bin_idx(i,j,k)) + real(fdiff * conjg(fdiff))
!!$           Fo(i,j,k) = (hf1(i,j,k) + hf2(i,j,k))/2.0                
!!$           A = real(Fo(i,j,k)); B = aimag(Fo(i,j,k))
!!$           A_sum(bin_idx(i,j,k)) = A_sum(bin_idx(i,j,k)) + A
!!$           AA_sum(bin_idx(i,j,k)) = AA_sum(bin_idx(i,j,k)) + A*A
!!$           B_sum(bin_idx(i,j,k)) = B_sum(bin_idx(i,j,k)) + B
!!$           BB_sum(bin_idx(i,j,k)) = BB_sum(bin_idx(i,j,k)) + B*B
!!$
!!$           ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
!!$           A1 = real(hf1(i,j,k));  A2 = real(hf2(i,j,k))
!!$           B1 = aimag(hf1(i,j,k)); B2 = aimag(hf2(i,j,k))
!!$           A1_sum(bin_idx(i,j,k)) = A1_sum(bin_idx(i,j,k)) + A1
!!$           A2_sum(bin_idx(i,j,k)) = A2_sum(bin_idx(i,j,k)) + A2
!!$           B1_sum(bin_idx(i,j,k)) = B1_sum(bin_idx(i,j,k)) + B1
!!$           B2_sum(bin_idx(i,j,k)) = B2_sum(bin_idx(i,j,k)) + B2
!!$
!!$           A1A2_sum(bin_idx(i,j,k)) = A1A2_sum(bin_idx(i,j,k)) + A1 * A2
!!$           B1B2_sum(bin_idx(i,j,k)) = B1B2_sum(bin_idx(i,j,k)) + B1 * B2
!!$
!!$           A1A1_sum(bin_idx(i,j,k)) = A1A1_sum(bin_idx(i,j,k)) + A1 * A1
!!$           B1B1_sum(bin_idx(i,j,k)) = B1B1_sum(bin_idx(i,j,k)) + B1 * B1
!!$
!!$           A2A2_sum(bin_idx(i,j,k)) = A2A2_sum(bin_idx(i,j,k)) + A2 * A2
!!$           B2B2_sum(bin_idx(i,j,k)) = B2B2_sum(bin_idx(i,j,k)) + B2 * B2
!!$           Fo(-i,-j,-k) = conjg(Fo(i,j,k))
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  if(debug) print*, 'bin_arr_count=', sum(bin_arr_count)
!!$
!!$  do ibin=0, nbin-1 !to make compatible with python arrays
!!$     bin_noise_var(ibin) = bin_arr_fdiff(ibin) / (bin_arr_count(ibin) * 4)
!!$     bin_total_var(ibin) = (AA_sum(ibin) + BB_sum(ibin))/bin_arr_count(ibin) &
!!$          - ((A_sum(ibin)/bin_arr_count(ibin))**2 + (B_sum(ibin)/bin_arr_count(ibin))**2)
!!$
!!$     F1F2_covar(ibin) = (A1A2_sum(ibin) + B1B2_sum(ibin)) / bin_arr_count(ibin) - &
!!$          (A1_sum(ibin) / bin_arr_count(ibin) * A2_sum(ibin) / bin_arr_count(ibin) + &
!!$          B1_sum(ibin) / bin_arr_count(ibin) * B2_sum(ibin) / bin_arr_count(ibin))
!!$
!!$     F1_var(ibin) = (A1A1_sum(ibin) + B1B1_sum(ibin))/bin_arr_count(ibin) - &
!!$          ((A1_sum(ibin)/bin_arr_count(ibin))**2 + (B1_sum(ibin)/bin_arr_count(ibin))**2)
!!$     F2_var(ibin) = (A2A2_sum(ibin) + B2B2_sum(ibin))/bin_arr_count(ibin) - &
!!$          ((A2_sum(ibin)/bin_arr_count(ibin))**2 + (B2_sum(ibin)/bin_arr_count(ibin))**2)
!!$
!!$     bin_sgnl_var(ibin) = F1F2_covar(ibin)
!!$     denominator = (sqrt(F1_var(ibin)) * sqrt(F2_var(ibin)))
!!$     bin_fsc(ibin) = F1F2_covar(ibin) / denominator
!!$     !print*,ibin,bin_total_var(ibin),denominator, denominator-bin_total_var(ibin)
!!$     if(debug)then
!!$        print*,ibin,bin_noise_var(ibin),bin_sgnl_var(ibin), &
!!$             bin_total_var(ibin),bin_fsc(ibin),bin_arr_count(ibin)
!!$     end if
!!$  end do
!!$
!!$  ! Calculate normalized structure factors
!!$  do i=xyzmin(1), xyzmax(1)
!!$     do j=xyzmin(2), xyzmax(2)
!!$        do k=xyzmin(3), 0
!!$           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
!!$           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
!!$           Eo(i,j,k,1) = hf1(i,j,k) / sqrt(F1_var(bin_idx(i,j,k)))
!!$           Eo(-i,-j,-k,1) = conjg(Eo(i,j,k,1))
!!$           Eo(i,j,k,2) = hf2(i,j,k) / sqrt(F2_var(bin_idx(i,j,k)))
!!$           Eo(-i,-j,-k,2) = conjg(Eo(i,j,k,2))
!!$           if(bin_sgnl_var(bin_idx(i,j,k)) <= 0.0) cycle ! singal var cannot be negative
!!$           Eo(i,j,k,3) = Fo(i,j,k)/sqrt(bin_total_var(bin_idx(i,j,k)))
!!$           Eo(-i,-j,-k,3) = conjg(Eo(i,j,k,3))
!!$        end do
!!$     end do
!!$  end do
!!$  call cpu_time(finish)
!!$  if(debug) print*, 'time for calculation(s) = ', finish-start
!!$end subroutine get_normalized_sf

subroutine get_normalized_sf(hf1,hf2,bin_idx,nbin,mode,nx,ny,nz, &
     Eo,bin_noise_var,bin_sgnl_var,bin_total_var,bin_fsc,bin_arr_count)
  implicit none
  real*8, parameter :: PI = 3.141592653589793

  integer, intent(in) :: nbin,mode,nx,ny,nz
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: hf1,hf2
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2, 3),intent(out) :: Eo
  integer, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in) :: bin_idx
  real*8, dimension(0:nbin-1),intent(out) :: bin_sgnl_var,bin_noise_var,bin_total_var,bin_fsc
  integer, dimension(0:nbin-1),intent(out) :: bin_arr_count
  ! locals
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2) :: Fo
  integer, dimension(3) :: nxyz
  real*8, dimension(0:nbin-1) :: A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
  real*8, dimension(0:nbin-1) :: A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
  real*8, dimension(0:nbin-1) :: bin_arr_fdiff,A_sum,B_sum,AA_sum,BB_sum
  real*8, dimension(0:nbin-1) :: F1_var, F2_var, F1F2_covar
  !
  complex*16  :: fdiff
  real*8     :: A,B,A1,A2,B1,B2,bin_sigvar,denominator
  real       :: start,finish
  integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin
  logical    :: debug,make_all_zero
  !
  debug         = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  Fo = dcmplx(0.0d0, 0.0d0)
  Eo = dcmplx(0.0d0, 0.0d0)

  bin_arr_fdiff = 0.0
  bin_sigvar    = 0.0
  bin_noise_var = 0.0
  bin_sgnl_var  = 0.0
  bin_total_var = 0.0

  F1F2_covar = 0.0
  F1_var = 0.0
  F2_var = 0.0
  bin_total_var = 0.0
  bin_fsc = 0.0

  A_sum = 0.0
  B_sum = 0.0
  AA_sum = 0.0
  BB_sum = 0.0

  A1_sum = 0.0; A2_sum = 0.0
  B1_sum = 0.0; B2_sum = 0.0
  A1A2_sum = 0.0; B1B2_sum = 0.0
  A1A1_sum = 0.0; B1B1_sum = 0.0
  A2A2_sum = 0.0; B2B2_sum = 0.0

  bin_arr_count = 0
  xyzmin = 0; xyzmax = 0
  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'Use only hemisphere data'
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1),xyzmax(2),0

  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           fdiff = hf1(i,j,k) - hf2(i,j,k)
           bin_arr_fdiff(bin_idx(i,j,k)) = bin_arr_fdiff(bin_idx(i,j,k)) + real(fdiff * conjg(fdiff))
           Fo(i,j,k) = (hf1(i,j,k) + hf2(i,j,k))/2.0                
           A = real(Fo(i,j,k)); B = aimag(Fo(i,j,k))
           A_sum(bin_idx(i,j,k)) = A_sum(bin_idx(i,j,k)) + A
           AA_sum(bin_idx(i,j,k)) = AA_sum(bin_idx(i,j,k)) + A*A
           B_sum(bin_idx(i,j,k)) = B_sum(bin_idx(i,j,k)) + B
           BB_sum(bin_idx(i,j,k)) = BB_sum(bin_idx(i,j,k)) + B*B

           ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
           A1 = real(hf1(i,j,k));  A2 = real(hf2(i,j,k))
           B1 = aimag(hf1(i,j,k)); B2 = aimag(hf2(i,j,k))
           A1_sum(bin_idx(i,j,k)) = A1_sum(bin_idx(i,j,k)) + A1
           A2_sum(bin_idx(i,j,k)) = A2_sum(bin_idx(i,j,k)) + A2
           B1_sum(bin_idx(i,j,k)) = B1_sum(bin_idx(i,j,k)) + B1
           B2_sum(bin_idx(i,j,k)) = B2_sum(bin_idx(i,j,k)) + B2

           A1A2_sum(bin_idx(i,j,k)) = A1A2_sum(bin_idx(i,j,k)) + A1 * A2
           B1B2_sum(bin_idx(i,j,k)) = B1B2_sum(bin_idx(i,j,k)) + B1 * B2

           A1A1_sum(bin_idx(i,j,k)) = A1A1_sum(bin_idx(i,j,k)) + A1 * A1
           B1B1_sum(bin_idx(i,j,k)) = B1B1_sum(bin_idx(i,j,k)) + B1 * B1

           A2A2_sum(bin_idx(i,j,k)) = A2A2_sum(bin_idx(i,j,k)) + A2 * A2
           B2B2_sum(bin_idx(i,j,k)) = B2B2_sum(bin_idx(i,j,k)) + B2 * B2
           Fo(-i,-j,-k) = conjg(Fo(i,j,k))
        end do
     end do
  end do

  if(debug) print*, 'bin_arr_count=', sum(bin_arr_count)

  do ibin=0, nbin-1 !to make compatible with python arrays
     bin_noise_var(ibin) = bin_arr_fdiff(ibin) / (bin_arr_count(ibin) * 4)
     bin_total_var(ibin) = (AA_sum(ibin) + BB_sum(ibin))/bin_arr_count(ibin) &
          - ((A_sum(ibin)/bin_arr_count(ibin))**2 + (B_sum(ibin)/bin_arr_count(ibin))**2)

     F1F2_covar(ibin) = (A1A2_sum(ibin) + B1B2_sum(ibin)) / bin_arr_count(ibin) - &
          (A1_sum(ibin) / bin_arr_count(ibin) * A2_sum(ibin) / bin_arr_count(ibin) + &
          B1_sum(ibin) / bin_arr_count(ibin) * B2_sum(ibin) / bin_arr_count(ibin))

     F1_var(ibin) = (A1A1_sum(ibin) + B1B1_sum(ibin))/bin_arr_count(ibin) - &
          ((A1_sum(ibin)/bin_arr_count(ibin))**2 + (B1_sum(ibin)/bin_arr_count(ibin))**2)
     F2_var(ibin) = (A2A2_sum(ibin) + B2B2_sum(ibin))/bin_arr_count(ibin) - &
          ((A2_sum(ibin)/bin_arr_count(ibin))**2 + (B2_sum(ibin)/bin_arr_count(ibin))**2)

     bin_sgnl_var(ibin) = F1F2_covar(ibin)
     denominator = (sqrt(F1_var(ibin)) * sqrt(F2_var(ibin)))
     bin_fsc(ibin) = F1F2_covar(ibin) / denominator
     !print*,ibin,bin_total_var(ibin),denominator, denominator-bin_total_var(ibin)
     if(debug)then
        print*,ibin,bin_noise_var(ibin),bin_sgnl_var(ibin), &
             bin_total_var(ibin),bin_fsc(ibin),bin_arr_count(ibin)
     end if
  end do

  ! Calculate normalized structure factors
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           Eo(i,j,k,1) = hf1(i,j,k) / sqrt(bin_total_var(bin_idx(i,j,k)))
           Eo(-i,-j,-k,1) = conjg(Eo(i,j,k,1))
           Eo(i,j,k,2) = hf2(i,j,k) / sqrt(bin_total_var(bin_idx(i,j,k)))
           Eo(-i,-j,-k,2) = conjg(Eo(i,j,k,2))
           !if(bin_sgnl_var(bin_idx(i,j,k)) <= 0.0) cycle ! singal var cannot be negative
           Eo(i,j,k,3) = Fo(i,j,k)/sqrt(bin_total_var(bin_idx(i,j,k)))
           Eo(-i,-j,-k,3) = conjg(Eo(i,j,k,3))
        end do
     end do
  end do
  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine get_normalized_sf


!!$subroutine calc_covar_and_fsc_betwn_anytwomaps(hf1,hf2,bin_idx,nbin,mode,&
!!$     F1F2_covar,bin_fsc,bin_arr_count,nx,ny,nz)
!!$  implicit none
!!$  integer,intent(in) :: nbin,mode,nx,ny,nz
!!$  integer,  dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: bin_idx
!!$  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: hf1,hf2
!!$  real*8,   dimension(0:nbin-1), intent(out) :: F1F2_covar,bin_fsc
!!$  integer,  dimension(0:nbin-1),intent(out) :: bin_arr_count
!!$  !
!!$  real*8,   dimension(0:nbin-1) :: F1_var,F2_var,A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
!!$  real*8,   dimension(0:nbin-1) :: A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
!!$  real*8    :: A1,A2,B1,B2
!!$  integer   :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,ibin
!!$  real      :: start, finish
!!$  logical   :: debug, make_all_zero 
!!$  !
!!$  debug = .FALSE.
!!$  make_all_zero = .FALSE.
!!$  if(mode == 1) debug = .TRUE.
!!$  call cpu_time(start)
!!$
!!$  F1F2_covar = 0.0
!!$  F1_var = 0.0
!!$  F2_var = 0.0
!!$  bin_fsc = 0.0
!!$
!!$  A1_sum = 0.0; A2_sum = 0.0
!!$  B1_sum = 0.0; B2_sum = 0.0
!!$  A1A2_sum = 0.0; B1B2_sum = 0.0
!!$  A1A1_sum = 0.0; B1B1_sum = 0.0
!!$  A2A2_sum = 0.0; B2B2_sum = 0.0
!!$
!!$  xmin = int(-nx/2); xmax = -(xmin+1)
!!$  ymin = int(-ny/2); ymax = -(ymin+1)
!!$  zmin = int(-nz/2); zmax = -(zmin+1)
!!$  !if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,zmax,']'
!!$  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'
!!$
!!$  bin_arr_count = 0
!!$  if(debug) print*, 'using hemisphere data...'
!!$  do i=xmin, xmax
!!$     do j=ymin, ymax
!!$        do k=zmin, 0 !zmax
!!$           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1)then
!!$              cycle
!!$           else
!!$              bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
!!$              ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
!!$              A1 = real(hf1(i,j,k));  A2 = real(hf2(i,j,k))
!!$              B1 = aimag(hf1(i,j,k)); B2 = aimag(hf2(i,j,k))
!!$              A1_sum(bin_idx(i,j,k)) = A1_sum(bin_idx(i,j,k)) + A1
!!$              A2_sum(bin_idx(i,j,k)) = A2_sum(bin_idx(i,j,k)) + A2
!!$              B1_sum(bin_idx(i,j,k)) = B1_sum(bin_idx(i,j,k)) + B1
!!$              B2_sum(bin_idx(i,j,k)) = B2_sum(bin_idx(i,j,k)) + B2
!!$
!!$              A1A2_sum(bin_idx(i,j,k)) = A1A2_sum(bin_idx(i,j,k)) + A1 * A2
!!$              B1B2_sum(bin_idx(i,j,k)) = B1B2_sum(bin_idx(i,j,k)) + B1 * B2
!!$
!!$              A1A1_sum(bin_idx(i,j,k)) = A1A1_sum(bin_idx(i,j,k)) + A1 * A1
!!$              B1B1_sum(bin_idx(i,j,k)) = B1B1_sum(bin_idx(i,j,k)) + B1 * B1
!!$
!!$              A2A2_sum(bin_idx(i,j,k)) = A2A2_sum(bin_idx(i,j,k)) + A2 * A2
!!$              B2B2_sum(bin_idx(i,j,k)) = B2B2_sum(bin_idx(i,j,k)) + B2 * B2
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  if(debug) print*,'ibin F1F2_covar(ibin) F1_var(ibin) F2_var(ibin) bin_fsc(ibin) bin_reflex_count'
!!$  do ibin=0, nbin-1 !to make compatible with python arrays
!!$
!!$     F1F2_covar(ibin) = (A1A2_sum(ibin) + B1B2_sum(ibin)) / bin_arr_count(ibin) - &
!!$          (A1_sum(ibin) / bin_arr_count(ibin) * A2_sum(ibin) / bin_arr_count(ibin) + &
!!$          B1_sum(ibin) / bin_arr_count(ibin) * B2_sum(ibin) / bin_arr_count(ibin))
!!$
!!$     F1_var(ibin) = (A1A1_sum(ibin) + B1B1_sum(ibin))/bin_arr_count(ibin) - &
!!$          ((A1_sum(ibin)/bin_arr_count(ibin))**2 + (B1_sum(ibin)/bin_arr_count(ibin))**2)
!!$     F2_var(ibin) = (A2A2_sum(ibin) + B2B2_sum(ibin))/bin_arr_count(ibin) - &
!!$          ((A2_sum(ibin)/bin_arr_count(ibin))**2 + (B2_sum(ibin)/bin_arr_count(ibin))**2)
!!$
!!$     !bin_sgnl_var(ibin) = F1F2_covar(ibin)
!!$     bin_fsc(ibin) = F1F2_covar(ibin) / (sqrt(F1_var(ibin)) * sqrt(F2_var(ibin)))
!!$
!!$     if(debug)then
!!$        print*,ibin,F1F2_covar(ibin),F1_var(ibin),F2_var(ibin),bin_fsc(ibin),bin_arr_count(ibin)
!!$     end if
!!$  end do
!!$
!!$  call cpu_time(finish)
!!$  if(debug) print*, 'time for loop = ', finish-start
!!$  return
!!$end subroutine calc_covar_and_fsc_betwn_anytwomaps


subroutine calc_covar_and_fsc_betwn_anytwomaps(hf1,hf2,bin_idx,nbin,mode,&
     F1F2_covar,bin_fsc,bin_arr_count,nx,ny,nz)
  implicit none
  integer,intent(in) :: nbin,mode,nx,ny,nz
  integer,  dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: bin_idx
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: hf1,hf2
  real*8,   dimension(0:nbin-1), intent(out) :: F1F2_covar,bin_fsc
  integer,  dimension(0:nbin-1),intent(out) :: bin_arr_count
  !
  real*8,   dimension(0:nbin-1) :: F1_var,F2_var,A1_sum,B1_sum,A2_sum,B2_sum,A1A2_sum,B1B2_sum
  real*8,   dimension(0:nbin-1) :: A1A1_sum,B1B1_sum,A2A2_sum,B2B2_sum
  real*8    :: A1,A2,B1,B2
  integer   :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,ibin, indx
  real      :: start, finish
  logical   :: debug, make_all_zero 
  !
  debug = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  F1F2_covar = 0.0
  F1_var = 0.0
  F2_var = 0.0
  bin_fsc = 0.0

  A1_sum = 0.0; A2_sum = 0.0
  B1_sum = 0.0; B2_sum = 0.0
  A1A2_sum = 0.0; B1B2_sum = 0.0
  A1A1_sum = 0.0; B1B1_sum = 0.0
  A2A2_sum = 0.0; B2B2_sum = 0.0

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  !if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,zmax,']'
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'

  bin_arr_count = 0
  if(debug) print*, 'using hemisphere data...'
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, 0 !zmax
           indx = bin_idx(k,j,i)
           if(indx < 0 .or. indx > nbin-1) cycle

           bin_arr_count(indx) = bin_arr_count(indx) + 1
           ! correspondence hf1 : A1 + iB1 ; hf2 = A2 + iB2
           A1 = real(hf1(k,j,i));  A2 = real(hf2(k,j,i))
           B1 = aimag(hf1(k,j,i)); B2 = aimag(hf2(k,j,i))
           A1_sum(indx) = A1_sum(indx) + A1
           A2_sum(indx) = A2_sum(indx) + A2
           B1_sum(indx) = B1_sum(indx) + B1
           B2_sum(indx) = B2_sum(indx) + B2
           A1A2_sum(indx) = A1A2_sum(indx) + A1 * A2
           B1B2_sum(indx) = B1B2_sum(indx) + B1 * B2
           A1A1_sum(indx) = A1A1_sum(indx) + A1 * A1
           B1B1_sum(indx) = B1B1_sum(indx) + B1 * B1
           A2A2_sum(indx) = A2A2_sum(indx) + A2 * A2
           B2B2_sum(indx) = B2B2_sum(indx) + B2 * B2

        end do
     end do
  end do
  if(debug) print*,'ibin F1F2_covar(ibin) F1_var(ibin) F2_var(ibin) bin_fsc(ibin) bin_reflex_count'
  do ibin=0, nbin-1 !to make compatible with python arrays

     F1F2_covar(ibin) = (A1A2_sum(ibin) + B1B2_sum(ibin)) / bin_arr_count(ibin) - &
          (A1_sum(ibin) / bin_arr_count(ibin) * A2_sum(ibin) / bin_arr_count(ibin) + &
          B1_sum(ibin) / bin_arr_count(ibin) * B2_sum(ibin) / bin_arr_count(ibin))

     F1_var(ibin) = (A1A1_sum(ibin) + B1B1_sum(ibin))/bin_arr_count(ibin) - &
          ((A1_sum(ibin)/bin_arr_count(ibin))**2 + (B1_sum(ibin)/bin_arr_count(ibin))**2)
     F2_var(ibin) = (A2A2_sum(ibin) + B2B2_sum(ibin))/bin_arr_count(ibin) - &
          ((A2_sum(ibin)/bin_arr_count(ibin))**2 + (B2_sum(ibin)/bin_arr_count(ibin))**2)

     !bin_sgnl_var(ibin) = F1F2_covar(ibin)
     bin_fsc(ibin) = F1F2_covar(ibin) / (sqrt(F1_var(ibin)) * sqrt(F2_var(ibin)))

     if(debug)then
        print*,ibin,F1F2_covar(ibin),F1_var(ibin),F2_var(ibin),bin_fsc(ibin),bin_arr_count(ibin)
     end if
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for loop = ', finish-start
  return
end subroutine calc_covar_and_fsc_betwn_anytwomaps



subroutine read_into_grid(bin_idx,bin_fsc,nbin,nx,ny,nz,fsc_weighted_grid)
  implicit none
  integer, intent(in) :: nbin,nx,ny,nz
  real*8,  dimension(0:nbin-1),intent(in) :: bin_fsc
  integer, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: bin_idx
  real*8,  dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: fsc_weighted_grid
  ! locals
  integer,   dimension(3) :: nxyz
  integer    :: i,j,k,xyzmin(3),xyzmax(3),indx
  !
  xyzmin = 0; xyzmax = 0
  fsc_weighted_grid = 0.0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  ! using Friedel's law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0!xyzmax(3)
           indx = bin_idx(k,j,i)
           if(indx < 0 .or. indx > nbin-1) cycle
           fsc_weighted_grid(k,j,i) = bin_fsc(indx)
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           fsc_weighted_grid(-k,-j,-i) = bin_fsc(indx)
        end do
     end do
  end do

end subroutine read_into_grid


subroutine read_into_grid2(bin_idx,bin_fsc,nbin,nx,ny,nz,wgrid)
  implicit none
  integer, intent(in) :: nbin,nx,ny,nz
  real*8,  dimension(0:nbin-1,2),intent(in) :: bin_fsc
  integer, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: bin_idx
  real*8,  dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2,2),intent(out) :: wgrid
  ! locals
  integer,   dimension(3) :: nxyz
  integer    :: i,j,k,xyzmin(3),xyzmax(3),indx,n
  !
  xyzmin = 0; xyzmax = 0
  wgrid = 0.0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  ! using Friedel's law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0!xyzmax(3)
           do n=1, 2
              indx = bin_idx(k,j,i)
              if(indx < 0 .or. indx > nbin-1) cycle
              wgrid(k,j,i,n) = bin_fsc(indx,n)
              if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
              wgrid(-k,-j,-i,n) = bin_fsc(indx,n)
           end do
        end do
     end do
  end do

end subroutine read_into_grid2

subroutine get_reciprocal_basis(uc,ucstar)
   implicit none
   real,dimension(6),intent(in) :: uc
   real,dimension(6),intent(out) :: ucstar
   real :: vol,sa,sb,sc
   vol = uc(1) * uc(2) * uc(3)
   ucstar(1) = uc(2)*uc(3)/vol
   ucstar(2) = uc(1)*uc(3)/vol
   ucstar(3) = uc(1)*uc(2)/vol
end subroutine get_reciprocal_basis

subroutine get_resol2(ucstar,h,k,l,resol)
   implicit none
   real,dimension(6),intent(in) :: ucstar
   real,intent(in) :: h,k,l
   real :: sa,sb,sc,s2,tmp
   real,intent(out) :: resol
   !
   sa = ucstar(1); sb = ucstar(2); sc = ucstar(3)
   s2 = ((h*sa)**2 + (k*sb)**2 + (l*sc)**2)/4.0
   if(s2 == 0.0) s2 = 1.0e-10 ! F(000) resolution hard coded
   tmp = sqrt(s2)
   resol = 1.0/(2.0*tmp)
   return
end subroutine get_resol2

subroutine get_resol(uc,h,k,l,resol)
  implicit none
  real,dimension(6),intent(in) :: uc
  !integer,intent(in) :: h,k,l
  real,intent(in) :: h,k,l
  real :: a,b,c,vol,sa,sb,sc,s2,tmp
  real,intent(out) :: resol
  !
  a = uc(1); b = uc(2); c = uc(3)
  !print*, a,b,c
  vol = a*b*c
  sa = b*c/vol; sb = a*c/vol; sc = a*b/vol
  s2 = ((h*sa)**2 + (k*sb)**2 + (l*sc)**2)/4.0
  if(s2 == 0.0) s2 = 1.0e-10 ! F(000) resolution hard coded
  tmp = sqrt(s2)
  resol = 1.0/(2.0*tmp)
  return
end subroutine get_resol

subroutine get_st(nx,ny,nz,t,st,s1,s2,s3)
  implicit none
  real*8, parameter :: PI = 3.141592653589793
  integer,intent(in) :: nx,ny,nz
  real*8,dimension(3),intent(in) :: t
  real*8 :: sv(3)
  complex*16 :: xj
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: st
  integer,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: s1,s2,s3
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax

  xj = dcmplx(0.0d0,1.0d0)
  st = dcmplx(0.0d0,1.0d0)
  sv = 0.0

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, zmax
           s1(i,j,k) = i
           s2(i,j,k) = j
           s3(i,j,k) = k
           sv(1) = i
           sv(2) = j
           sv(3) = k
           st(i,j,k) = exp(2.0d0 * PI * xj * dot_product(t,sv))
        end do
     end do
  end do
  return
end subroutine get_st

subroutine get_xyz(uc,nx,ny,nz,xyz)
  implicit none
  integer,intent(in) :: nx,ny,nz
  real*8, intent(in) :: uc(6)
  real*8,dimension(-nz/2:nz/2-1,-ny/2:ny/2-1,-nx/2:nx/2-1,3),intent(out) :: xyz
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax
  xyz = 0.0d0
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, zmax
           xyz(k,j,i,1) = (k/nz) * uc(3)
           xyz(k,j,i,2) = (j/ny) * uc(2)
           xyz(k,j,i,3) = (i/nx) * uc(1)
        end do
     end do
  end do
  return  
end subroutine get_xyz



subroutine fsc_weight_calculation(fsc_weighted_grid,F1,F2,bin_idx,nbin,mode,nx,ny,nz)
  implicit none
  integer,intent(in) :: nbin,mode,nx,ny,nz
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: F1,F2
  integer,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx 
  real*8,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: fsc_weighted_grid
  real*8,dimension(0:nbin-1) :: bin_fsc, F1F2_covar
  integer,dimension(0:nbin-1) :: bin_arr_count
  integer :: i,j,k,indx,xmin,xmax,ymin,ymax,zmin,zmax
  real :: start, finish
  logical :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'

  bin_fsc  = 0.0
  fsc_weighted_grid = 0.0

  call calc_covar_and_fsc_betwn_anytwomaps(F1,F2,bin_idx,nbin,mode,&
       F1F2_covar,bin_fsc,bin_arr_count,nx,ny,nz)

  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, 0
           indx = bin_idx(i,j,k)
           if(indx < 0 .or. indx > nbin-1) cycle
           fsc_weighted_grid(i,j,k) = bin_fsc(indx)/(1.0 - bin_fsc(indx)**2)
           if(k == zmin .or. j == ymin .or. i == xmin) cycle
           fsc_weighted_grid(-i,-j,-k) = fsc_weighted_grid(i,j,k)
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for fsc and weights calculation = ', finish-start
end subroutine fsc_weight_calculation

subroutine calc_avg_maps(all_maps,bin_idx,s_grid,wgt,Bf_arr,uc,nbin,nmaps,nbf,mode,nx,ny,nz,avgmaps_all)
  implicit none
  !
  integer,  intent(in) :: nbin,nmaps,nbf,mode,nx,ny,nz
  real,     dimension(6),                                           intent(in) :: uc
  real,     dimension(nbf),                                         intent(in) :: Bf_arr
  real*8,   dimension(0:nmaps-1,0:nmaps-1,0:nbin-1),                intent(in) :: wgt
  !real*8,   dimension(0:nbin-1),                                    intent(in) :: res_arr
  integer,  dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx
  real,     dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: s_grid
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1),intent(in) :: all_maps
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1,nbf),intent(out) :: avgmaps_all
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,imap,jmap,ibf
  real    :: resol,s,start, finish, lowres, highres
  real    :: Bfac(nbf)
  logical :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.

  avgmaps_all = dcmplx(0.0d0, 0.0d0)

  call cpu_time(start)
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,zmax,']'
  !if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'

  Bfac = -1.0 * Bf_arr

  call get_resol(uc,0.0,0.0,0.0,lowres)
  call get_resol(uc,real(xmin),real(ymin),real(zmin),highres)
  print*,'low high resol= ', lowres, highres

  print*, 'Calculating average maps...'
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, zmax
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1)then
              !avgmaps_all(i,j,k) = 0.0
              cycle
           else
              !s = 1.0 / res_arr(bin_idx(i,j,k)) ! Assuming that resolution never becomes zero
              !call get_resol(uc,real(i),real(j),real(k),resol)
              !if(resol == 0.0)then
              !   print*, i,j,k
              !end if

              !s = 1.0 / resol
              s = s_grid(i,j,k)

              do imap=0, nmaps-1
                 do jmap=0, nmaps-1
                    do ibf = 1, nbf
                       avgmaps_all(i,j,k,imap,ibf) = avgmaps_all(i,j,k,imap,ibf) + &
                            all_maps(i,j,k,jmap) * wgt(imap,jmap,bin_idx(i,j,k)) * &
                            exp(-1.0 * (Bfac(ibf)/4.0) * s**2) !B-factor sharpening/blurring
                       !if(k==zmin .or. j==ymin .or. i==xmin) cycle
                       !avgmaps_all(-i,-j,-k,imap,ibf) = conjg(avgmaps_all(i,j,k,imap,ibf))
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for avgerage maps calculation = ', finish-start  
end subroutine calc_avg_maps

! This subroutine is not used at the moment
subroutine calc_avg_maps_3d(all_maps,bin_idx,Smat,Fmat,Tinv,Bf_arr,uc,nbin, &
     nmaps,nbf,mode,nx,ny,nz,avgmaps_all)
  implicit none
  !
  integer,  intent(in) :: nbin,nmaps,nbf,mode,nx,ny,nz
  real,     dimension(6),                                           intent(in) :: uc
  real,     dimension(nbf),                                         intent(in) :: Bf_arr
  real*8,   dimension(0:nmaps-1,0:nmaps-1,0:nbin-1),                intent(in) :: Smat,Tinv
  integer,  dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx
  real,     dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1,0:nmaps-1),intent(in) :: Fmat
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1),intent(in) :: all_maps
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1,nbf),intent(out) :: avgmaps_all
  ! locals
  real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1,0:nmaps-1) :: F_dot_Tinv
  real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:nmaps-1,0:nmaps-1) :: S_dot_FdotTinv
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,imap,jmap,ibf
  real    :: resol,s,start, finish, lowres, highres
  real    :: Bfac(nbf)
  logical :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.

  avgmaps_all = dcmplx(0.0d0, 0.0d0)
  Bfac = -1.0 * Bf_arr
  F_dot_Tinv = 0.0
  S_dot_FdotTinv = 0.0

  call cpu_time(start)
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,zmax,']'



  call get_resol(uc,0.0,0.0,0.0,lowres)
  call get_resol(uc,real(xmin),real(ymin),real(zmin),highres)
  print*,'low high resol= ', lowres, highres

  print*, 'Calculating average maps...'
  ! New code
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, zmax
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1)then
              !avgmaps_all(i,j,k) = 0.0
              cycle
           else
              call get_resol(uc,real(i),real(j),real(k),resol)
              if(resol == 0.0)then
                 print*, i,j,k
              end if
              !
              do imap=0, nmaps-1
                 do jmap=0, nmaps-1
                    F_dot_Tinv(i,j,k,imap,jmap) = F_dot_Tinv(i,j,k,imap,jmap) + &
                         Fmat(i,j,k,imap,jmap) * Tinv(jmap,imap,bin_idx(i,j,k))
                 end do
              end do
              do imap=0, nmaps-1
                 do jmap=0, nmaps-1
                    S_dot_FdotTinv(i,j,k,imap,jmap) = S_dot_FdotTinv(i,j,k,imap,jmap) + &
                         Smat(imap,jmap,bin_idx(i,j,k)) * F_dot_Tinv(i,j,k,jmap,imap) 
                 end do
              end do
              !
              s = 1.0 / resol
              do imap=0, nmaps-1
                 do jmap=0, nmaps-1
                    do ibf = 1, nbf
                       avgmaps_all(i,j,k,imap,ibf) = avgmaps_all(i,j,k,imap,ibf) + &
                            S_dot_FdotTinv(i,j,k,imap,jmap) * all_maps(i,j,k,jmap) * &
                            exp(-1.0 * (Bfac(ibf)/4.0) * s**2) !B-factor sharpening/blurring
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do
  ! End new code

  call cpu_time(finish)
  if(debug) print*, 'time for fsc and weights calculation = ', finish-start  

end subroutine calc_avg_maps_3d

subroutine apply_bfactor_to_map(mapin,Bf_arr,uc,nx,ny,nz,nbf,mode,all_mapout)
  implicit none
  !
  integer,  intent(in) :: nbf,mode,nx,ny,nz
  real,     dimension(6),                                           intent(in) :: uc
  real,     dimension(nbf),                                         intent(in) :: Bf_arr
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: mapin
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nbf),intent(out) :: all_mapout
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,ibf
  real    :: resol,s,start, finish!, lowres, highres
  real    :: Bfac(nbf)
  logical :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.

  all_mapout = dcmplx(0.0d0, 0.0d0)
  Bfac = -1.0 * Bf_arr

  call cpu_time(start)
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,zmax,']'
  !if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'

  !call get_resol(uc,0.0,0.0,0.0,lowres)
  !call get_resol(uc,real(xmin),real(ymin),real(zmin),highres)
  !print*,'low high resol= ', lowres, highres

  print*, 'Applying B factors to map...'
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, zmax
           call get_resol(uc,real(i),real(j),real(k),resol)
           if(resol == 0.0)then
              print*, i,j,k
              print*, 'Resol cannot be zero or negative!'
              stop
           end if

           s = 1.0 / resol
           do ibf = 1, nbf
              all_mapout(i,j,k,ibf) = all_mapout(i,j,k,ibf) + &
                   mapin(i,j,k) * exp(-1.0 * (Bfac(ibf)/4.0) * s**2) !B-factor sharpening/blurring
              !if(k==zmin .or. j==ymin .or. i==xmin) cycle
              !all_mapout(-i,-j,-k,ibf) = conjg(all_mapout(i,j,k,ibf))
           end do
        end do
     end do
  end do
  call cpu_time(finish)
  if(debug) print*, 'time for map blurring/sharpening = ', finish-start 

end subroutine apply_bfactor_to_map


subroutine tricubic(RM,F,FRS,nc,mode,nx,ny,nz)
  implicit none
  real*8,dimension(3,3),intent(in):: RM
  integer,intent(in):: nx,ny,nz,mode,nc
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nc),intent(in):: F
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nc),intent(out):: FRS
  integer :: x1(3,-1:2)
  real*8 :: x(3),xd(3),s(3)
  real*8 :: ul,ul2
  real*8 :: vl(4)
  !real :: high_res,resol
  integer :: i,j,h,k,l,nmin,ncopies
  logical :: debug
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic
  complex*16 :: fl(-1:2,-1:2,-1:2,nc),esz(-1:2,-1:2,nc),esyz(-1:2,nc)
  complex*16 :: esxyz(nc)
  !
  FRS = dcmplx(0.0d0, 0.0d0)
  x = 0.0d0
  xd = 0.0d0
  fl = dcmplx(0.0d0, 0.0d0)
  esz = dcmplx(0.0d0, 0.0d0)
  esyz = dcmplx(0.0d0, 0.0d0)
  esxyz = dcmplx(0.0d0, 0.0d0)
  
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.
  !   Body
  ncopies = nc
  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2
  nmin = min(nx,ny,nz)

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  if(debug) write(*,*) nxyz,nxyzmn,nxyzmx

  do l = zmin, zmax
     do k = ymin, ymax
        do h = xmin, 0!xmax
           s(1) = real(h)
           s(2) = real(k)
           s(3) = real(l)
           x = matmul(transpose(RM),s)
           do i = 1, 3
              x1(i,0) = floor(x(i))
              xd(i) = x(i) - real(x1(i,0))
              if(abs(xd(i)).gt.1.0) then
                 print*, 'Something is wrong ',xd(i)
                 stop
              endif
              x1(i,1)  = x1(i,0) + 1
              x1(i,2)  = x1(i,0) + 2
              x1(i,-1) = x1(i,0) - 1
           end do
           !
           !  Careful here: we may get to the outside of the array
           do i = 1,3
              do j= -1,2
                 x1(i,j) = min(nxyzmx(i),max(nxyzmn(i),x1(i,j)))
              enddo
           enddo
           do ic=1, ncopies
              fl(-1:2,-1:2,-1:2,ic) = F(x1(1,-1:2),x1(2,-1:2),x1(3,-1:2),ic)
           end do

           !
           !  Alternattive implementation
           !  along z
           ul = xd(3)
           ul2 = ul*ul
           vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
           vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
           vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
           vl(4) = ul2*(ul-1.0d0)
           vl = 0.5d0*vl
           do ic=1,ncopies
              do j=-1,2
                 do i=-1,2
                    esz(i,j,ic) = dot_product(vl,fl(i,j,-1:2,ic))
                 enddo
              enddo
           end do
           ul = xd(2)
           ul2 = ul*ul
           vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
           vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
           vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
           vl(4) = ul2*(ul-1.0d0)
           vl = 0.5d0*vl
           do ic=1,ncopies
              do i=-1,2
                 esyz(i,ic) = dot_product(vl,esz(i,-1:2,ic))
              enddo
           end do
           ul = xd(1)
           ul2 = ul*ul
           vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
           vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
           vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
           vl(4) = ul2*(ul-1.0d0)
           vl = 0.5d0*vl
           do ic=1,ncopies
              esxyz(ic) =  dot_product(vl,esyz(-1:2,ic))
              FRS(h,k,l,ic) = esxyz(ic)
              if((h == xmin).or.(k == ymin).or.(l == zmin)) cycle
              FRS(-h,-k,-l,ic) = conjg(esxyz(ic))
           end do
        end do
     end do
  end do
  return
end subroutine tricubic

subroutine mtz2_3d(h,k,l,f,nobs,nx,ny,nz,f3d)
  implicit none
  integer, intent(in)  :: nx,ny,nz,nobs
  real, dimension(nobs), intent(in)  :: h,k,l
  complex*16, dimension(nobs),intent(in)  :: f
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: f3d
  complex*16, dimension(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2) :: f3d2
  integer :: i,i1,i2,i3
  !
  f3d = dcmplx(0.0d0, 0.0d0)
  f3d2 = dcmplx(0.0d0, 0.0d0)
  do i = 1, nobs
     i1 = int(h(i))
     i2 = int(k(i))
     i3 = int(l(i))
     f3d2(-i3,-i2,-i1)    = f(i) !to comply with MRC
     f3d2(i3,i2,i1) = conjg(f(i))    
  end do
  f3d = f3d2(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2)
  return
end subroutine mtz2_3d

subroutine prepare_hkl(hf1,bin_idx,cbin,nx,ny,nz,mode,h,k,l,ampli,phase)
  implicit none
  real*8,    parameter :: PI = 3.141592653589793
  integer, intent(in) :: cbin,mode,nx,ny,nz
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: hf1
  integer, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: bin_idx
  !integer, dimension(nx*ny*(ny+2)/2),intent(out) :: h,k,l
  !real*8, dimension(nx*ny*(ny+2)/2),intent(out) :: ampli, phase
  integer, dimension(4*cbin*cbin*cbin),intent(out) :: h,k,l
  real*8, dimension(4*cbin*cbin*cbin),intent(out) :: ampli, phase
  ! locals
  integer,   dimension(3)          :: nxyz
  integer    :: xyzmin(3),xyzmax(3)
  integer    :: i1,i2,i3,j
  real*8       :: hf1_real, hf1_imag
  logical    :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  !call cpu_time(start)

  xyzmin = 0; xyzmax = 0
  nxyz = (/ nx, ny, nz /)
  h = 0; k = 0; l = 0
  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax

  print*, 'cbin= ', cbin
  j = 0
  ! using Friedel's Law
  do i1=xyzmin(1), xyzmax(1)
     do i2=xyzmin(2), xyzmax(2)
        do i3=xyzmin(3), 0 !xyzmax(3)
           if(bin_idx(i1,i2,i3) < 0 .or. bin_idx(i1,i2,i3) > cbin) cycle
           j = j + 1
           h(j) = i1
           k(j) = i2
           l(j) = i3
           hf1_real = real(hf1(i1,i2,i3))
           hf1_imag = aimag(hf1(i1,i2,i3))
           ampli(j) = sqrt(hf1_real**2 + hf1_imag**2)
           phase(j) = atan2(hf1_imag,hf1_real)*180.0d0/PI
        end do
     end do
  end do
  return
end subroutine prepare_hkl

subroutine prepare_hkl2(f_arr,nx,ny,nz,nmap,mode,hkl,ampli,phase)
  implicit none
  real*8,    parameter :: PI = 3.141592653589793
  integer, intent(in) :: mode,nx,ny,nz,nmap
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nmap),intent(in)  :: f_arr
  integer, dimension(nx*ny*(ny+2)/2,3),intent(out) :: hkl
  real*8, dimension(nx*ny*(ny+2)/2,nmap),intent(out) :: ampli
  real*8, dimension(nx*ny*(ny+2)/2),intent(out) :: phase
  ! locals
  integer,   dimension(3)          :: nxyz
  integer    :: xyzmin(3),xyzmax(3)
  integer    :: i1,i2,i3,j,im
  real*8       :: f_real, f_imag
  logical    :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  !call cpu_time(start)

  xyzmin = 0; xyzmax = 0
  hkl = 0
  ampli = 0.0d0
  phase = 0.0d0
  nxyz = (/ nx, ny, nz /)
  
  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax

  j = 0
  ! using Friedel's Law
  do i1=xyzmin(1), xyzmax(1)
     do i2=xyzmin(2), xyzmax(2)
        do i3=xyzmin(3), 0 !xyzmax(3)
           !if(i3 < 0) cycle
           j = j + 1
           hkl(j,1) = i1
           hkl(j,2) = i2
           hkl(j,3) = i3
           do im=1, nmap
              f_real = real(f_arr(i1,i2,i3,im))
              f_imag = aimag(f_arr(i1,i2,i3,im))
              ampli(j,im) = sqrt(f_real**2 + f_imag**2)
              if(im==1) phase(j) = atan2(f_imag,f_real)*180.0d0/PI
           end do
        end do
     end do
  end do
  return
end subroutine prepare_hkl2

subroutine prepare_hkl_bfac(s_grid,f1,f2,Bfac,nx,ny,nz,nbf,mode,h,k,l,ampli,noise,phase)
  implicit none
  real*8, parameter :: PI = 3.141592653589793

  integer, intent(in) :: mode,nx,ny,nz,nbf
  real, dimension(nbf), intent(in) :: Bfac
  real, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: s_grid
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in)  :: f1,f2
  integer, dimension(nx*ny*(ny+2)/2),intent(out) :: h,k,l
  real*8, dimension(nx*ny*(ny+2)/2,nbf),intent(out) :: ampli,noise
  real*8, dimension(nx*ny*(ny+2)/2),intent(out) :: phase
  ! locals
  integer, dimension(3) :: nxyz
  integer :: xyzmin(3),xyzmax(3)
  integer :: i1,i2,i3,j,ibf
  real*8  :: f1_real,f1_imag,f2_real,f2_imag
  real :: s
  logical :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  !call cpu_time(start)

  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  h = 0; k = 0; l = 0
  
  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1),xyzmax(1),0

  j = 0
  ! using Friedel's Law
  do i1=xyzmin(1), xyzmax(1)
     do i2=xyzmin(2), xyzmax(2)
        do i3=xyzmin(3), 0!xyzmax(3)
           !if(i3 < 0) cycle
           !call get_resol(uc,real(i1),real(i2),real(i3),resol)
           !s = 1.0 / resol
           s = s_grid(i1,i2,i3)
           j = j + 1
           h(j) = i1
           k(j) = i2
           l(j) = i3
           f1_real = real(f1(i1,i2,i3))
           f1_imag = aimag(f1(i1,i2,i3))
           f2_real = real(f2(i1,i2,i3))
           f2_imag = aimag(f2(i1,i2,i3))
           phase(j) = atan2(f1_imag,f1_real)*180.0d0/PI
           do ibf=1, nbf
              if(ibf == 1)then
                 ampli(j,1) = sqrt(f1_real**2 + f1_imag**2)
                 noise(j,1) = sqrt(f2_real**2 + f2_imag**2)
               else
                  ampli(j,ibf) = ampli(j,1) * exp((Bfac(ibf)/4.0) * s**2)
                  noise(j,ibf) = noise(j,1) * exp((Bfac(ibf)/4.0) * s**2)
               end if
           end do
        end do
     end do
  end do
  return
end subroutine prepare_hkl_bfac

subroutine add_random_phase_beyond(F_ori,F_all_random,bin_idx,rand_idx,nx,ny,nz,F_beyond_random)
  implicit none
  real*8,    parameter :: PI = 3.141592653589793
  integer,   intent(in) :: nx,ny,nz,rand_idx
  integer,  dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: F_ori,F_all_random
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: F_beyond_random
  ! locals
  integer,dimension(3) :: nxyz
  real       :: start,finish
  integer    :: xyzmin(3),xyzmax(3)
  integer    :: i1,i2,i3

  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  F_beyond_random = dcmplx(0.0d0, 0.0d0)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  print*, 'xyzmin = ', xyzmin
  print*, 'xyzmax = ', xyzmax(1), xyzmax(3), 0

  call cpu_time(start)
  do i1=xyzmin(1), xyzmax(1)
     do i2=xyzmin(2), xyzmax(2)
        do i3=xyzmin(3), 0
           if(bin_idx(i3,i2,i1) <= rand_idx)then
              F_beyond_random(i3,i2,i1) = F_ori(i3,i2,i1)
              if(i3 == xyzmin(3) .or. i2 == xyzmin(2) .or. i1 == xyzmin(1)) cycle
              F_beyond_random(-i3,-i2,-i1) = conjg(F_ori(i3,i2,i1))
           else
              F_beyond_random(i3,i2,i1) = F_all_random(i3,i2,i1)
              if(i3 == xyzmin(3) .or. i2 == xyzmin(2) .or. i1 == xyzmin(1)) cycle
              F_beyond_random(-i3,-i2,-i1) = conjg(F_beyond_random(i3,i2,i1))
           end if
        end do
     end do
  end do
  call cpu_time(finish)
  print*, 'time for phase randomisation looping(s) = ', finish-start
  return
end subroutine add_random_phase_beyond

subroutine cutmap(fin,bin_idx,smax,mode,nbin,nx,ny,nz,fout)
  implicit none
  integer,   intent(in) :: smax,mode,nbin,nx,ny,nz
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: fin
  integer,   dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: bin_idx
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: fout
  ! locals
  integer,   dimension(3) :: nxyz
  integer,   dimension(0:nbin-1) :: bin_arr_count
  !
  real       :: start,finish
  integer    :: i,j,k,n,xyzmin(3),xyzmax(3)
  logical    :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  fout = dcmplx(0.0d0, 0.0d0)

  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax
  if(debug) print*, 'nbin=', nbin

  ! Using Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           n = n + 1
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           if(bin_idx(i,j,k) > smax) cycle
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           fout(i,j,k) = fin(i,j,k)
           if(k==xyzmin(3) .or. j==xyzmin(2) .or. i==xyzmin(1)) cycle
           fout(-i,-j,-k) = conjg(fout(i,j,k))
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for Eo calculation(s) = ', finish-start
end subroutine cutmap

subroutine cutmap_arr(fin,bin_idx,smax,mode,nbin,nx,ny,nz,n,fout)
  implicit none
  integer,   intent(in) :: smax,mode,nbin,nx,ny,nz,n
  complex*16, dimension(n, -nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: fin
  integer,   dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: bin_idx
  complex*16, dimension(n, -nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: fout
  ! locals
  integer,   dimension(3) :: nxyz
  real       :: start,finish
  integer    :: i,j,k,l,xyzmin(3),xyzmax(3)
  logical    :: debug
  !
  debug         = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)
  fout = dcmplx(0.0d0, 0.0d0)
  xyzmin = 0; xyzmax = 0
  nxyz = (/ nx, ny, nz /)
  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax
  if(debug) print*, 'nbin=', nbin

  ! Using Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0
           do l=1, n
              if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
              if(bin_idx(i,j,k) > smax) cycle
              fout(l,i,j,k) = fin(l,i,j,k)
              if(k==xyzmin(3) .or. j==xyzmin(2) .or. i==xyzmin(1)) cycle
              fout(l,-i,-j,-k) = conjg(fout(l,i,j,k))
           end do
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for Eo calculation(s) = ', finish-start
end subroutine cutmap_arr


subroutine trilinear(RM,F,FRS,ncopies,mode,nx,ny,nz)
  implicit none
  real*8,intent(in):: RM(3,3)
  integer,intent(in):: nx,ny,nz,mode,ncopies
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,ncopies),intent(in):: F
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,ncopies),intent(out):: FRS
  ! locals
  integer :: x0(3),x1(3)
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  real*8 :: x(3),xd(3),s(3)
  complex*16 :: c000,c001,c010,c011,c100,c101,c110,c111,c00,c01,c10,c11,c0,c1,c
  integer :: h,k,l,i
  integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic


  FRS = dcmplx(0.0d0, 0.0d0)
  x = 0.0d0
  xd = 0.0d0

  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  !write(*,*) nxyz,nxyzmn,nxyzmx

  do l = zmin, zmax
     do k = ymin, ymax
        outer: do h = xmin, 0!xmax
           s(1) = h
           s(2) = k
           s(3) = l
           x = matmul(transpose(RM),s)
           do i = 1, 3
              !x(i)  = dot_product(RM(:,i),s) ! Note that RM is now transposed
              x0(i) = floor(x(i))
              x1(i) = x0(i) + 1
              if((nxyzmx(i) < x0(i)) .or. (x0(i) < nxyzmn(i)) &
                   .or. (nxyzmx(i) < x1(i)) .or. (x1(i) < nxyzmn(i)))then
                 cycle outer
              end if
              xd(i) = (x(i)-real(x0(i)))!/(x1(i)-x0(i))
              if(abs(xd(i)).gt.1.0) then
                 print*, 'Something is wrong ',xd(i)
                 stop
              endif
           end do
           !  Careful here: we may get to the outside of the array
           do i = 1,3
              x1(i) = min(nxyzmx(i),max(nxyzmn(i),x1(i)))
           enddo
           do ic = 1, ncopies
              c000 = F(x0(1),x0(2),x0(3),ic)
              c001 = F(x0(1),x0(2),x1(3),ic)
              c010 = F(x0(1),x1(2),x0(3),ic)
              c011 = F(x0(1),x1(2),x1(3),ic)
              c100 = F(x1(1),x0(2),x0(3),ic)
              c101 = F(x1(1),x0(2),x1(3),ic)
              c110 = F(x1(1),x1(2),x0(3),ic)
              c111 = F(x1(1),x1(2),x1(3),ic)
              ! Interpolation along x direction
              c00 = c000*(1.0d0-xd(1)) + c100*xd(1)
              c01 = c001*(1.0d0-xd(1)) + c101*xd(1)
              c10 = c010*(1.0d0-xd(1)) + c110*xd(1)
              c11 = c011*(1.0d0-xd(1)) + c111*xd(1)
              ! Interpolation along y direction
              c0 = c00*(1.0d0-xd(2)) + c10*xd(2)
              c1 = c01*(1.0d0-xd(2)) + c11*xd(2)
              ! Interpolation along z direction
              c = c0*(1.0d0-xd(3)) + c1*xd(3)
              FRS(h,k,l,ic) = c
              if((h == xmin).or.(k == ymin).or.(l == zmin)) cycle
              FRS(-h,-k,-l,ic) = conjg(c)
           end do
        end do outer
     end do
  end do
  return
end subroutine trilinear


subroutine trilinear_new(RM,uc,F,FRS,highres,nx,ny,nz)
  implicit none
  real*8,intent(in):: RM(3,3)
  real,intent(in) :: uc(6), highres
  integer,intent(in):: nx,ny,nz
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in):: F
  complex*16,dimension(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2),intent(out):: FRS
  complex*16,dimension(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2) :: F2
  !integer :: x0(3),x1(3)
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  !real*8 :: x(3),xd(3),s(3),xdr(3)
  complex*16 :: c000,c001,c010,c011,c100,c101,c110,c111,c00,c01,c10,c11,c0,c1,cc
  integer :: h,k,l,i,indx
  integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic
  integer :: x01,x02,x03,x11,x12,x13
  real*8 :: x1,x2,x3,xd1,xd2,xd3,xdr1,xdr2,xdr3
  real :: a,b,c, resol, denominator, tol, highres_inv
  real      :: start, finish, st2,end2
  logical :: cubic

  cubic = .false.
  tol = 1e-5

  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  a = uc(1)
  b = uc(2)
  c = uc(3)

  highres_inv = (1.0/highres)**2
  if((a-b) < tol .and. (a-c) < tol)then
     denominator = a**2
     cubic = .true.
  end if

  call cpu_time(start)
  do l = zmin, zmax
     do k = ymin, ymax
        do h = xmin, 0!xmax
           if(cubic)then
              resol = (h**2 + k**2 + l**2)/denominator
           else
              resol = ((h/a)**2 + (k/b)**2 + (l/c)**2)
           end if
           if(resol > highres_inv) cycle
           x1 = RM(1,1) * h + RM(2,1) * k + RM(3,1) * l
           x2 = RM(1,2) * h + RM(2,2) * k + RM(3,2) * l
           x3 = RM(1,3) * h + RM(2,3) * k + RM(3,3) * l
           x01 = floor(x1)
           x02 = floor(x2)
           x03 = floor(x3)
           x11 = x01 + 1
           x12 = x02 + 1
           x13 = x03 + 1
           x11 = min(nxyzmx(1),max(nxyzmn(1),x11))
           x12 = min(nxyzmx(2),max(nxyzmn(2),x12))
           x13 = min(nxyzmx(3),max(nxyzmn(3),x13))
           xd1 = x1 - real(x01)
           xd2 = x2 - real(x02)
           xd3 = x3 - real(x03)
           xdr1 = 1.0d0 - xd1
           xdr2 = 1.0d0 - xd2
           xdr3 = 1.0d0 - xd3
           FRS(h,k,l) = ((F(x01,x02,x03)*xdr1 + F(x11,x02,x03)*xd1) * xdr2 &
                + (F(x01,x12,x03)*xdr1 + F(x11,x12,x03)*xd1) * xd2) * xdr3 &
                + ((F(x01,x02,x13)*xdr1 + F(x11,x02,x13)*xd1) * xdr2 &
                + (F(x01,x12,x13)*xdr1 + F(x11,x12,x13)*xd1) * xd2) * xd3
           FRS(-h,-k,-l) = conjg(FRS(h,k,l))
        end do
     end do
  end do
  call cpu_time(finish)
  print*, 'time for loop = ', finish-start
  return
end subroutine trilinear_new


subroutine trilinear2(F,bin_idx,RM,nbin,ncopies,mode,nx,ny,nz,FRS)
  implicit none
  real*8,intent(in):: RM(3,3)
  integer,intent(in):: nx,ny,nz,mode,ncopies,nbin
  integer, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: bin_idx
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,ncopies),intent(in):: F
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,ncopies),intent(out):: FRS
  ! locals
  integer :: x0(3),x1(3)
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  real*8 :: x(3),xd(3),s(3),xdr(3)
  complex*16 :: c000,c001,c010,c011,c100,c101,c110,c111,c00,c01,c10,c11,c0,c1,c
  integer :: h,k,l,i,indx
  integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic
  real      :: start, finish


  FRS = dcmplx(0.0d0, 0.0d0)
  x = 0.0d0
  xd = 0.0d0

  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  !write(*,*) nxyz,nxyzmn,nxyzmx
  call cpu_time(start)
  do l = zmin, zmax
     do k = ymin, ymax
        outer: do h = xmin, 0!xmax
           indx = bin_idx(h,k,l)
           if(indx < 0 .or. indx > nbin-1) cycle
           s(1) = h
           s(2) = k
           s(3) = l
           x = matmul(transpose(RM),s)
           do i = 1, 3
              !x(i)  = dot_product(RM(:,i),s) ! Note that RM is now transposed
              x0(i) = floor(x(i))
              x1(i) = x0(i) + 1
              if((nxyzmx(i) < x0(i)) .or. (x0(i) < nxyzmn(i)) &
                   .or. (nxyzmx(i) < x1(i)) .or. (x1(i) < nxyzmn(i)))then
                 cycle outer
              end if
              xd(i) = (x(i)-real(x0(i)))!/(x1(i)-x0(i))
              if(abs(xd(i)).gt.1.0) then
                 print*, 'Something is wrong ',xd(i)
                 stop
              endif
           end do
           !  Careful here: we may get to the outside of the array
           do i = 1,3
              x1(i) = min(nxyzmx(i),max(nxyzmn(i),x1(i)))
           enddo
           xdr = 1.0d0 - xd
           do ic = 1, ncopies
              c000 = F(x0(1),x0(2),x0(3),ic)
              c001 = F(x0(1),x0(2),x1(3),ic)
              c010 = F(x0(1),x1(2),x0(3),ic)
              c011 = F(x0(1),x1(2),x1(3),ic)
              c100 = F(x1(1),x0(2),x0(3),ic)
              c101 = F(x1(1),x0(2),x1(3),ic)
              c110 = F(x1(1),x1(2),x0(3),ic)
              c111 = F(x1(1),x1(2),x1(3),ic)
              ! Interpolation along x direction
              !c00 = c000*(1.0d0-xd(1)) + c100*xd(1)
              !c01 = c001*(1.0d0-xd(1)) + c101*xd(1)
              !c10 = c010*(1.0d0-xd(1)) + c110*xd(1)
              !c11 = c011*(1.0d0-xd(1)) + c111*xd(1)
              c00 = c000*xdr(1) + c100*xd(1)
              c01 = c001*xdr(1) + c101*xd(1)
              c10 = c010*xdr(1) + c110*xd(1)
              c11 = c011*xdr(1) + c111*xd(1)
              ! Interpolation along y direction
              !c0 = c00*(1.0d0-xd(2)) + c10*xd(2)
              !c1 = c01*(1.0d0-xd(2)) + c11*xd(2)
              c0 = c00*xdr(2) + c10*xd(2)
              c1 = c01*xdr(2) + c11*xd(2)
              ! Interpolation along z direction
              !c = c0*(1.0d0-xd(3)) + c1*xd(3)
              c = c0*xdr(3) + c1*xd(3)
              FRS(h,k,l,ic) = c
              if((h == xmin).or.(k == ymin).or.(l == zmin)) cycle
              FRS(-h,-k,-l,ic) = conjg(c)
           end do
        end do outer
     end do
  end do
  call cpu_time(finish)
  !print*, 'time for loop = ', finish-start
  return
end subroutine trilinear2


subroutine trilinearn(F,RM,ncopies,mode,nx,ny,nz,FRS)
  implicit none
  integer,intent(in):: nx,ny,nz,mode,ncopies
  real*8,intent(in):: RM(ncopies,3,3)
  complex*16,dimension(ncopies,-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in):: F
  complex*16,dimension(ncopies,-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out):: FRS
  ! locals
  integer :: x0(3),x1(3)
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  real*8 :: x(3),xd(3),s(3)
  complex*16 :: c000,c001,c010,c011,c100,c101,c110,c111,c00,c01,c10,c11,c0,c1,c
  integer :: h,k,l,i
  integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic
  real      :: start, finish

  FRS = dcmplx(0.0d0, 0.0d0)
  x = 0.0d0
  xd = 0.0d0

  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  call cpu_time(start)
  do l = zmin, zmax
     do k = ymin, ymax
        outer: do h = xmin, 0
           s(1) = h
           s(2) = k
           s(3) = l
           do ic = 1, ncopies
              x = matmul(transpose(RM(ic,:,:)),s)
              do i = 1, 3
                 x0(i) = floor(x(i))
                 x1(i) = x0(i) + 1
                 if((nxyzmx(i) < x0(i)) .or. (x0(i) < nxyzmn(i)) &
                      .or. (nxyzmx(i) < x1(i)) .or. (x1(i) < nxyzmn(i)))then
                    cycle outer
                 end if
                 xd(i) = (x(i)-real(x0(i)))
                 if(abs(xd(i)).gt.1.0) then
                    print*, 'Something is wrong ',xd(i)
                    stop
                 endif
              end do
              do i = 1,3
                 x1(i) = min(nxyzmx(i),max(nxyzmn(i),x1(i)))
              enddo
              c000 = F(ic,x0(1),x0(2),x0(3))
              c001 = F(ic,x0(1),x0(2),x1(3))
              c010 = F(ic,x0(1),x1(2),x0(3))
              c011 = F(ic,x0(1),x1(2),x1(3))
              c100 = F(ic,x1(1),x0(2),x0(3))
              c101 = F(ic,x1(1),x0(2),x1(3))
              c110 = F(ic,x1(1),x1(2),x0(3))
              c111 = F(ic,x1(1),x1(2),x1(3))
              ! Interpolation along x direction
              c00 = c000*(1.0d0-xd(1)) + c100*xd(1)
              c01 = c001*(1.0d0-xd(1)) + c101*xd(1)
              c10 = c010*(1.0d0-xd(1)) + c110*xd(1)
              c11 = c011*(1.0d0-xd(1)) + c111*xd(1)
              ! Interpolation along y direction
              c0 = c00*(1.0d0-xd(2)) + c10*xd(2)
              c1 = c01*(1.0d0-xd(2)) + c11*xd(2)
              ! Interpolation along z direction
              c = c0*(1.0d0-xd(3)) + c1*xd(3)
              FRS(ic,h,k,l) = c
              if((h == xmin).or.(k == ymin).or.(l == zmin)) cycle
              FRS(ic,-h,-k,-l) = conjg(c)
           end do
        end do outer
     end do
  end do
  call cpu_time(finish)
  !print*, 'time for loop = ', finish-start
  return
end subroutine trilinearn

                 
subroutine trilinear_map(RM,arr1,arr2,nx,ny,nz,mode)
  implicit none
  real*8,dimension(3,3),intent(in):: RM
  integer,intent(in):: nx,ny,nz,mode
  !real*8,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in):: arr1
  !real*8,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out):: arr2
  real*8,dimension(-nz/2:(nz-2)/2,-ny/2:(ny-2)/2,-nx/2:(nx-2)/2),intent(in):: arr1
  real*8,dimension(-nz/2:(nz-2)/2,-ny/2:(ny-2)/2,-nx/2:(nx-2)/2),intent(out):: arr2
  ! locals
  logical :: debug
  integer :: x0(3),x1(3)
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  real*8 :: x(3),xd(3),s(3)
  real*8 :: c000,c001,c010,c011,c100,c101,c110,c111,c00,c01,c10,c11,c0,c1,c
  integer :: h,k,l,i
  integer :: xmin,xmax,ymin,ymax,zmin,zmax

  debug = .FALSE.
  if(mode == 1) debug = .TRUE.

  arr2 = 0.0d0
  x = 0.0d0

  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  if(debug) write(*,*) nxyz,nxyzmn,nxyzmx

  do h = zmin, zmax
     do k = ymin, ymax
        outer: do l = xmin, xmax
           s(1) = real(h) 
           s(2) = real(k) 
           s(3) = real(l)
           x = matmul(transpose(RM),s)

           do i = 1, 3
              x0(i) = floor(x(i))
              !x0(i) = min(nxyzmx(i),max(nxyzmn(i),x0(i)))
              x1(i) = x0(i) + 1
              !x1(i) = min(nxyzmx(i),max(nxyzmn(i),x1(i)))
              if((nxyzmx(i) < x0(i)) .or. (x0(i) < nxyzmn(i)) &
                   .or. (nxyzmx(i) < x1(i)) .or. (x1(i) < nxyzmn(i)))then
                 cycle outer
              end if
              xd(i) = (x(i)-real(x0(i)))!/(x1(i)-x0(i))
              if(abs(xd(i)).gt.1.0) then
                 print*, 'Something is wrong ',xd(i)
                 stop
              endif
           end do

           !  Careful here: we may get to the outside of the array
           do i = 1,3
              x1(i) = min(nxyzmx(i),max(nxyzmn(i),x1(i)))
           enddo
           !if(k==65) print*, h, k, l, x0, x1
           !c000 = arr1(x0(1),x0(2),x0(3))
           !c001 = arr1(x0(1),x0(2),x1(3))
           !c010 = arr1(x0(1),x1(2),x0(3))
           !c011 = arr1(x0(1),x1(2),x1(3))
           !c100 = arr1(x1(1),x0(2),x0(3))
           !c101 = arr1(x1(1),x0(2),x1(3))
           !c110 = arr1(x1(1),x1(2),x0(3))
           !c111 = arr1(x1(1),x1(2),x1(3))
           c000 = arr1(x0(3),x0(2),x0(1))
           c001 = arr1(x1(3),x0(2),x0(1))
           c010 = arr1(x0(3),x1(2),x0(1))
           c011 = arr1(x1(3),x1(2),x0(1))
           c100 = arr1(x0(3),x0(2),x1(1))
           c101 = arr1(x1(3),x0(2),x1(1))
           c110 = arr1(x0(3),x1(2),x1(1))
           c111 = arr1(x1(3),x1(2),x1(1))
           ! Interpolation along x direction
           c00 = c000*(1.0d0-xd(1)) + c100*xd(1)
           c01 = c001*(1.0d0-xd(1)) + c101*xd(1)
           c10 = c010*(1.0d0-xd(1)) + c110*xd(1)
           c11 = c011*(1.0d0-xd(1)) + c111*xd(1)
           ! Interpolation along y direction
           c0 = c00*(1.0d0-xd(2)) + c10*xd(2)
           c1 = c01*(1.0d0-xd(2)) + c11*xd(2)
           ! Interpolation along z direction
           c = c0*(1.0d0-xd(3)) + c1*xd(3)
           !arr2(h,k,l) = c
           arr2(l,k,h) = c ! 16.09.2020
           !print*, h, k, l, arr1(h,k,l), arr2(h,k,l)
        end do outer
     end do
  end do
  return
end subroutine trilinear_map

subroutine tricubic_map(RM,arr1,arr2,ncopies,mode,nx,ny,nz)
  implicit none
  real*8,dimension(3,3),intent(in):: RM
  integer,intent(in):: nx,ny,nz,mode,ncopies
  real*8,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,ncopies),intent(in):: arr1
  real*8,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,ncopies),intent(out):: arr2
  integer :: x1(3,-1:2)
  real*8 :: x(3),xd(3),s(3)
  real*8 :: ul,ul2
  real*8 :: vl(4)
  integer :: i,j,h,k,l,nmin
  logical :: debug
  integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
  integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic
  real*8 :: fl(-1:2,-1:2,-1:2,ncopies),esz(-1:2,-1:2,ncopies),esyz(-1:2,ncopies)
  real*8 :: esxyz(ncopies)
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.
  arr2 = 0.0d0
  fl = 0.0d0
  esz = 0.0d0
  esyz = 0.0d0
  esxyz = 0.0d0
  x1 = 0
  x = 0.0d0; xd = 0.0d0; s = 0.0d0; vl = 0.0d0
  !   Body
  nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
  nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
  nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2
  nmin = min(nx,ny,nz)

  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  if(debug) write(*,*) nxyz,nxyzmn,nxyzmx

  do l = zmin, zmax
     do k = ymin, ymax
        outer: do h = xmin, xmax
           s(1) = real(h)
           s(2) = real(k)
           s(3) = real(l)
           x = matmul(transpose(RM),s)
           do i = 1, 3
              x1(i,0) = floor(x(i))
              xd(i) = x(i) - real(x1(i,0))
              if(abs(xd(i)).gt.1.0) then
                 print*, 'Something is wrong ',xd(i)
                 stop
              endif
              x1(i,1)  = x1(i,0) + 1
              x1(i,2)  = x1(i,0) + 2
              x1(i,-1) = x1(i,0) - 1
              if((nxyzmx(i) < x1(i,0)) .or. (x1(i,0) < nxyzmn(i)) &
                   .or. (nxyzmx(i) < x1(i,1)) .or. (x1(i,1) < nxyzmn(i)) &
                   .or. (nxyzmx(i) < x1(i,2)) .or. (x1(i,2) < nxyzmn(i)) &
                   .or. (nxyzmx(i) < x1(i,-1)) .or. (x1(i,-1) < nxyzmn(i)))then
                 cycle outer
              end if
           end do
           !
           !  Careful here: we may get to the outside of the array
           !do i = 1,3
           !   do j= -1,2
           !      x1(i,j) = min(nxyzmx(i),max(nxyzmn(i),x1(i,j)))
           !   enddo
           !enddo
           do ic=1, ncopies
              fl(-1:2,-1:2,-1:2,ic) = arr1(x1(1,-1:2),x1(2,-1:2),x1(3,-1:2),ic)
              !fl(-1:2,-1:2,-1:2,ic) = arr1(x1(3,-1:2),x1(2,-1:2),x1(1,-1:2),ic)
           end do

           !
           !  Alternattive implementation
           !  along z
           ul = xd(3)
           ul2 = ul*ul
           vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
           vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
           vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
           vl(4) = ul2*(ul-1.0d0)
           vl = 0.5d0*vl
           do ic=1,ncopies
              do j=-1,2
                 do i=-1,2
                    esz(i,j,ic) = dot_product(vl,fl(i,j,-1:2,ic))
                 enddo
              enddo
           end do
           ul = xd(2)
           ul2 = ul*ul
           vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
           vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
           vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
           vl(4) = ul2*(ul-1.0d0)
           vl = 0.5d0*vl
           do ic=1,ncopies
              do i=-1,2
                 esyz(i,ic) = dot_product(vl,esz(i,-1:2,ic))
              enddo
           end do
           ul = xd(1)
           ul2 = ul*ul
           vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
           vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
           vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
           vl(4) = ul2*(ul-1.0d0)
           vl = 0.5d0*vl
           do ic=1,ncopies
              esxyz(ic) =  dot_product(vl,esyz(-1:2,ic))
              !arr2(h,k,l,ic) = esxyz(ic)
              arr2(l,k,h,ic) = esxyz(ic)
           end do
        end do outer
     end do
  end do
  return
end subroutine tricubic_map

subroutine average_intensity_in_bins(F1,F2,bin_idx,nbin,mode,nx,ny,nz,I1_mean,I2_mean)
  implicit none
  integer,intent(in) :: nbin,mode,nx,ny,nz
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: F1,F2
  integer,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx 
  real*8,dimension(0:nbin-1),intent(out) :: I1_mean, I2_mean
  real*8,dimension(0:nbin-1) :: I1_sum,I2_sum
  integer,dimension(0:nbin-1) :: bin_arr_count
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,ibin
  real :: start, finish
  logical :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)

  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'

  bin_arr_count = 0
  I1_sum = 0.0
  I2_sum = 0.0
  I1_mean = 0.0
  I2_mean = 0.0

  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, 0 !zmax
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           I1_sum(bin_idx(i,j,k)) = I1_sum(bin_idx(i,j,k)) + real(F1(i,j,k) * conjg(F1(i,j,k)))
           I2_sum(bin_idx(i,j,k)) = I2_sum(bin_idx(i,j,k)) + real(F2(i,j,k) * conjg(F2(i,j,k)))
        end do
     end do
  end do

  if(debug) print*,'ibin   < I1 >     < I2 >    No. reflx'
  do ibin=0, nbin-1 !to make compatible with python arrays
     I1_mean(ibin) = I1_sum(ibin)/bin_arr_count(ibin)
     I2_mean(ibin) = I2_sum(ibin)/bin_arr_count(ibin)
     if(debug)then
        print*,ibin,I1_mean(ibin),I2_mean(ibin),bin_arr_count(ibin)
     end if
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for calculation = ', finish-start
end subroutine average_intensity_in_bins

subroutine setbin(nbin,low_res,high_res,res_arr,binmode)
   implicit none
   integer,intent(in) :: nbin,binmode
   real,intent(in) :: low_res,high_res
   real,intent(out) :: res_arr(nbin)
   real :: stolmax3,stolmin3,stolinc,stol3max,stol2min,stol2max
   integer :: ibin
   !
   if(binmode == 1)then
      ! equal volume in reci. space
      print*, 'Equal volume in reciprocal space'
      stolmax3 = (1./(2.*high_res))**3
      stolmin3 = (1./(2.*low_res))**3
      stol2min = exp(log(stolmin3)*2./3.)
   else
      ! Linearly increasing volumes
      print*, 'Linearly increasing volume in reciprocal space'
      stolmax3 = (1./(2.*high_res))**2
      stolmin3 = (1./(2.*low_res))**2
      stol2min = stolmin3
   end if
   stolinc  = (stolmax3-stolmin3)/nbin
   do ibin = 1, nbin
      stol3max = stolinc * ibin + stolmin3
      if(binmode == 1)then
         stol2max = exp(log(stol3max)*2./3.)
         res_arr(ibin) = sqrt(1.0/stol2max)/2.0
      else
         res_arr(ibin) = sqrt(1.0/stol3max)/2.0
      end if
   end do
   return
 end subroutine setbin

subroutine binarr_1d(res_arr,resol_1d,bin_idx,nbin,nx,nref)
  implicit none
  integer,intent(in) :: nx,nbin
  real, dimension(0:nbin-1),intent(in) :: res_arr
  real, dimension(nx),intent(in) :: resol_1d
  integer, dimension(nx),intent(out) :: bin_idx
  integer, dimension(0:nbin-1),intent(out) :: nref

  real       :: tmp_val,tmp_min,val
  integer    :: i,ibin,mnloc
  logical    :: noval

  bin_idx = -100
  nref = 0
  noval = .false.
  do i=1, nx
     ! Find the matching bin to resol
     do ibin = 0, nbin - 1
        if (resol_1d(i) > res_arr(0) .or. &
             resol_1d(i) < res_arr(nbin-1))then
           noval = .true.
           exit
        end if
        noval = .false.
        val = sqrt((res_arr(ibin) - resol_1d(i))**2)
        if(ibin == 0)then
           tmp_val = val; tmp_min = val
           mnloc = ibin
        else
           tmp_val = val
           if(tmp_val < tmp_min)then
              tmp_min = val
              mnloc = ibin
           end if
        end if
     end do
     if(.not.(noval))then
        bin_idx(i) = mnloc
        nref(mnloc) = nref(mnloc) + 1
     end if

  end do
  do ibin=0, nbin-1
     print*, res_arr(ibin), nref(ibin)
  end do

end subroutine binarr_1d

subroutine bin_data_1d(intensity,binidx_1d,nbin,nx,nmax,bindat)
  implicit none
  integer, intent(in) :: nx,nbin,nmax
  real*8, dimension(nx),intent(in) :: intensity
  integer, dimension(nx),intent(in) :: binidx_1d
  real*8, dimension(0:nbin-1,nmax),intent(out) :: bindat
  integer, dimension(0:nbin-1) :: cntarr
  integer :: i
  real :: start, end
  !
  cntarr = 0
  !print*, nx, nmax, nbin
  call cpu_time(start)
  do i=1, nx
     if(binidx_1d(i) == -100) cycle
     cntarr(binidx_1d(i)) = cntarr(binidx_1d(i)) + 1
     bindat(binidx_1d(i),cntarr(binidx_1d(i))) = intensity(i)
  end do
  call cpu_time(end)
  print*, 'time for data binning: ', end-start, 's'
end subroutine bin_data_1d


subroutine bin_idx_with_given_nbin(uc,nbin,res_arr,bin_idx,nx,ny,nz,binmode)
  implicit none
  integer,intent(in) :: nbin,nx,ny,nz,binmode
  real,dimension(6),intent(in) :: uc
  integer,dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: bin_idx
  real,dimension(0:nbin-1),intent(out) :: res_arr
  !
  integer,dimension(3) :: nxyz
  real :: low_res,high_res,resol,tmp_val,tmp_min,val
  integer :: i,j,k,xyzmin(3),xyzmax(3),ibin,mnloc
  !
  call get_resol(uc,real(1),real(0),real(0),low_res)
  call get_resol(uc,real(nx/2),real(0),real(0),high_res)
  call setbin(nbin,low_res,high_res,res_arr,binmode)
  print*,"Low res=",low_res,"High res=",high_res

  bin_idx = -100
  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  print*, 'xyzmin = ', xyzmin
  print*, 'xyzmax = ', xyzmax(1), xyzmax(2), 0

  ! Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           call get_resol(uc,real(i),real(j),real(k),resol)
           if(resol < high_res .or. resol > low_res) cycle
           ! Find the matching bin to resol
           do ibin = 0, nbin - 1
              val = sqrt((res_arr(ibin) - resol)**2)
              if(ibin == 0)then
                 tmp_val = val; tmp_min = val
                 mnloc = ibin 
              else
                 tmp_val = val
                 if(tmp_val < tmp_min)then
                    tmp_min = val
                    mnloc = ibin
                 end if
              end if
           end do
           bin_idx(i,j,k) = mnloc
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           bin_idx(-i,-j,-k) = mnloc
        end do
     end do
  end do
  return
end subroutine bin_idx_with_given_nbin

subroutine apply_bin_scale(mapin,bin_idx,bin_scale,nbins,nmaps,nx,ny,nz,mode,all_mapout)
  implicit none
  !
  integer,intent(in) :: nmaps,mode,nbins,nx,ny,nz
  real,dimension(0:nbins-1),intent(in) :: bin_scale
  integer,dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in) :: bin_idx
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nmaps),intent(in) :: mapin
  complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nmaps),intent(out) :: all_mapout
  integer :: i,j,k,xmin,xmax,ymin,ymax,zmin,zmax,ibf
  real    :: start, finish
  logical :: debug
  !
  debug = .FALSE.
  if(mode == 1) debug = .TRUE.

  all_mapout = dcmplx(0.0d0, 0.0d0)

  call cpu_time(start)
  xmin = int(-nx/2); xmax = -(xmin+1)
  ymin = int(-ny/2); ymax = -(ymin+1)
  zmin = int(-nz/2); zmax = -(zmin+1)
  if(debug) print*, '[',xmin,xmax,'],[', ymin,ymax,'],[',zmin,0,']'


  print*, 'Applying bin-scale to maps...'
  do i=xmin, xmax
     do j=ymin, ymax
        do k=zmin, 0
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbins-1) cycle
           do ibf = 1, nmaps
              all_mapout(i,j,k,ibf) = mapin(i,j,k,ibf) * bin_scale(bin_idx(i,j,k))
              if(k==zmin .or. j==ymin .or. i==xmin) cycle
              all_mapout(-i,-j,-k,ibf) = conjg(all_mapout(i,j,k,ibf))
           end do
        end do
     end do
  end do
  call cpu_time(finish)
  if(debug) print*, 'time for map blurring/sharpening = ', finish-start 

end subroutine apply_bin_scale

subroutine calc_power_spectrum(Fo,bin_idx,nbin,mode,nx,ny,nz,bin_total_var)
  implicit none
  real*8,    parameter :: PI = 3.141592653589793

  integer,                intent(in) :: nbin,mode,nx,ny,nz
  complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: Fo
  integer,   dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in) :: bin_idx
  real*8,    dimension(0:nbin-1),intent(out) :: bin_total_var
  ! locals
  integer,   dimension(3)          :: nxyz
  integer,   dimension(0:nbin-1) :: bin_arr_count
  real*8,    dimension(0:nbin-1) :: A_sum,B_sum,AA_sum,BB_sum
  !
  real*8     :: A,B
  real       :: start,finish
  integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin
  logical    :: debug,make_all_zero
  !
  debug         = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  bin_total_var = 0.0


  A_sum = 0.0
  B_sum = 0.0
  AA_sum = 0.0
  BB_sum = 0.0


  bin_arr_count = 0
  xyzmin = 0; xyzmax = 0
  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'Use only hemisphere data'
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax(1),xyzmax(2),0

  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           !New calculation of total-var     
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           A = real(Fo(i,j,k)); B = aimag(Fo(i,j,k))
           A_sum(bin_idx(i,j,k)) = A_sum(bin_idx(i,j,k)) + A
           AA_sum(bin_idx(i,j,k)) = AA_sum(bin_idx(i,j,k)) + A*A
           B_sum(bin_idx(i,j,k)) = B_sum(bin_idx(i,j,k)) + B
           BB_sum(bin_idx(i,j,k)) = BB_sum(bin_idx(i,j,k)) + B*B
           ! end of new calculation
        end do
     end do
  end do

  do ibin=0, nbin-1 !to make compatible with python arrays
     bin_total_var(ibin) = (AA_sum(ibin) + BB_sum(ibin))/bin_arr_count(ibin) &
          - ((A_sum(ibin)/bin_arr_count(ibin))**2 + (B_sum(ibin)/bin_arr_count(ibin))**2)
  
     if(debug)then
        print*,ibin,bin_total_var(ibin),bin_arr_count(ibin)
     end if
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for calculation(s) = ', finish-start
end subroutine calc_power_spectrum

subroutine calc_derivatives(e0,e1,wgrid,w2grid,sv,dFRS,dRdq,xyz_sum,vol,nx,ny,nz,df,ddf)
  implicit none
  real*8, parameter :: PI = 3.141592653589793
  integer, intent(in) :: nx,ny,nz
  complex*16, dimension(nx,ny,nz),intent(in)  :: e0,e1
  complex*16, dimension(nx,ny,nz,3),intent(in)  :: dFRS
  real*8, dimension(nx,ny,nz),intent(in) :: wgrid,w2grid
  real*8, dimension(nx,ny,nz,3),intent(in) :: sv
  real*8, dimension(3,3,3),intent(in) :: dRdq
  real*8, dimension(6),intent(in):: xyz_sum
  real*8, intent(in) :: vol
  real*8, dimension(6),intent(out) :: df
  real*8, dimension(6,6),intent(out) :: ddf
  !locals
  integer :: i,j,k,l,n
  real*8 :: tp2
  complex*16 :: xj,tpi
  real*8, dimension(3,3) :: a,b
  real*8, dimension(nx,ny,nz) :: wfsc

  xj = dcmplx(0.0d0,1.0d0)
  a = 0.0d0
  b = 0.0d0
  df = 0.0d0
  ddf = 0.0d0

  tp2 = (2.0d0 * PI)**2
  tpi = (2.0d0 * PI * xj)
  ! translation derivatives
  do i=1,3
     df(i) = real(sum(wgrid * e0 * conjg(e1 * tpi * sv(:,:,:,i))))
     do j=1,3
        if(i==1 .or. (i>1 .and. j>=i))then
           ddf(i,j) = -tp2 * sum(w2grid * sv(:,:,:,i) * sv(:,:,:,j))
        else
           ddf(i,j) = ddf(j,i)
        end if
     end do
  end do
  !
  ! rotation derivatives
  do i=1,3
     a = 0.0d0
     do k=1,3
        do l=1,3
           if(k==1 .or. (k>1 .and. l>=k))then
              a(k,l) = sum(wgrid * real(conjg(e0) * (dFRS(:,:,:,k) * &
                   sv(:,:,:,l) * dRdq(i,k,l))))
           else
              a(k,l) = a(l,k)
           end if
        end do
     end do
     df(i+3) = sum(a)
  end do
  wfsc = wgrid * real(conjg(e0) * e1)
  do i=1,3
     do j=1,3
        if(i==1 .or. (i>1 .and. j>=i))then
           b = 0.0d0
           n = 0
           do k=1,3
              do l=1,3
                 if(k==1 .or. (k>1 .and. l>=k))then
                    n = n + 1
                    b(k,l) = (-tp2/vol) * xyz_sum(n) * &
                         sum(wfsc * sv(:,:,:,k) * sv(:,:,:,l) * &
                         dRdq(i,k,l) * dRdq(j,k,l))
                 else
                    b(k,l) = b(l,k)
                 end if
              end do
           end do
           ddf(i+3,j+3) = sum(b)
        else
           ddf(i+3,j+3) = ddf(j+3,i+3)
        end if
     end do
  end do
end subroutine calc_derivatives

subroutine calc_derivatives2(e0,e1,wgrid,w2grid,sv,dFRS,dRdq,dqda,xyz_sum,vol,nx,ny,nz,df,ddf)
   implicit none
   real*8, parameter :: PI = 3.141592653589793
   integer, intent(in) :: nx,ny,nz
   complex*16, dimension(nx,ny,nz),intent(in)  :: e0,e1
   complex*16, dimension(nx,ny,nz,3),intent(in)  :: dFRS
   real*8, dimension(nx,ny,nz),intent(in) :: wgrid,w2grid
   real*8, dimension(nx,ny,nz,3),intent(in) :: sv
   real*8, dimension(3,3,3),intent(in) :: dRdq
   real*8, dimension(6),intent(in):: xyz_sum
   real*8, dimension(3),intent(in):: dqda
   real*8, intent(in) :: vol
   real*8, dimension(6),intent(out) :: df
   real*8, dimension(6,6),intent(out) :: ddf
   !locals
   integer :: i,j,k,l,n
   real*8 :: tp2
   complex*16 :: xj,tpi
   real*8, dimension(3,3) :: a,b
   real*8, dimension(nx,ny,nz) :: wfsc
 
   xj = dcmplx(0.0d0,1.0d0)
   a = 0.0d0
   b = 0.0d0
   df = 0.0d0
   ddf = 0.0d0
 
   tp2 = (2.0d0 * PI)**2
   tpi = (2.0d0 * PI * xj)
   ! translation derivatives
   do i=1,3
      df(i) = real(sum(wgrid * e0 * conjg(e1 * tpi * sv(:,:,:,i))))
      do j=1,3
         if(i==1 .or. (i>1 .and. j>=i))then
            ddf(i,j) = -tp2 * sum(w2grid * sv(:,:,:,i) * sv(:,:,:,j))
         else
            ddf(i,j) = ddf(j,i)
         end if
      end do
   end do
   !
   ! rotation derivatives
   do i=1,3
      a = 0.0d0
      do k=1,3
         do l=1,3
            if(k==1 .or. (k>1 .and. l>=k))then
               a(k,l) = sum(wgrid * real(conjg(e0) * (dFRS(:,:,:,k) * &
                    sv(:,:,:,l) * dRdq(i,k,l) * dqda(i))))
            else
               a(k,l) = a(l,k)
            end if
         end do
      end do
      df(i+3) = sum(a)
   end do
   wfsc = wgrid * real(conjg(e0) * e1)
   do i=1,3
      do j=1,3
         if(i==1 .or. (i>1 .and. j>=i))then
            b = 0.0d0
            n = 0
            do k=1,3
               do l=1,3
                  if(k==1 .or. (k>1 .and. l>=k))then
                     n = n + 1
                     b(k,l) = (-tp2/vol) * xyz_sum(n) * &
                          sum(wfsc * sv(:,:,:,k) * sv(:,:,:,l) * &
                          dRdq(i,k,l) * dRdq(j,k,l) * dqda(j) * dqda(i))
                  else
                     b(k,l) = b(l,k)
                  end if
               end do
            end do
            ddf(i+3,j+3) = sum(b)
         else
            ddf(i+3,j+3) = ddf(j+3,i+3)
         end if
      end do
   end do
 end subroutine calc_derivatives2

subroutine differencemap(Fo,Fc,bin_idx,res_arr,smax,mode,nbin,nx,ny,nz, &
     diffmap)
  implicit none

  integer,   intent(in) :: mode,nbin,nx,ny,nz
  real,      intent(in) :: smax
  real,      dimension(0:nbin-1),intent(in) :: res_arr
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2), &
       intent(in) :: Fo, Fc
  integer,   dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),&
       intent(in) :: bin_idx
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:3),&
       intent(out) :: diffmap
  ! locals
  integer,   dimension(3)        :: nxyz
  integer,   dimension(0:nbin-1) :: bin_arr_count
  real*8,    dimension(0:nbin-1) :: FoFc_sum,FcFc_sum,FoFo_sum
  real*8,    dimension(0:nbin-1) :: scale_d
  !
  integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin,cbin
  real       :: start, finish
  logical    :: debug,make_all_zero
  real*8     :: d_tmp, d2_tmp

  !
  debug         = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  bin_arr_count = 0
  scale_d       = 0.0d0
  FoFo_sum      = 0.0d0
  FoFc_sum      = 0.0d0
  FcFc_sum      = 0.0d0
  diffmap       = dcmplx(0.0d0, 0.0d0)
  
  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax

  cbin =  minloc(sqrt((res_arr - smax)**2), DIM=1) - 1
  
  ! Not using Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), xyzmax(3)
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle  
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           FoFo_sum(bin_idx(i,j,k)) = FoFo_sum(bin_idx(i,j,k)) + &
                real(Fo(i,j,k) * conjg(Fo(i,j,k)))
           FoFc_sum(bin_idx(i,j,k)) = FoFc_sum(bin_idx(i,j,k)) + &
                real(Fo(i,j,k) * conjg(Fc(i,j,k)))
           FcFc_sum(bin_idx(i,j,k)) = FcFc_sum(bin_idx(i,j,k)) + &
                real(Fc(i,j,k) * conjg(Fc(i,j,k)))
        end do
     end do
  end do

  ! Estimate scale_d in resolution bins
  if(debug) print*,'ibin  resol  scale_d  sqrt(F1^2/F2^2)  bincount'
  do ibin=0, nbin-1
     d_tmp = FoFc_sum(ibin) / FcFc_sum(ibin)
     d2_tmp = sqrt(FoFo_sum(ibin) / FcFc_sum(ibin))
     if(d_tmp < 0.0) make_all_zero = .TRUE.
     if(make_all_zero)then
        scale_d(ibin) = 0.0
     else
        scale_d(ibin) = d_tmp
     end if
     print*,ibin, res_arr(ibin), scale_d(ibin), d2_tmp, bin_arr_count(ibin)
  end do
  
  ! Calculate difference map
  if(cbin > (nbin-1)) cbin = nbin-1
  print*
  print*, 'Resolution of Difference map: ', res_arr(cbin), ' A'
  print*, 'Calculating map. Please wait...'
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), xyzmax(3)
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           if(bin_idx(i,j,k) > cbin) cycle
           diffmap(i,j,k,0) = Fo(i,j,k) - scale_d(bin_idx(i,j,k)) * Fc(i,j,k)
           diffmap(i,j,k,1) = scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)
           diffmap(i,j,k,2) = Fo(i,j,k)
           diffmap(i,j,k,3) = scale_d(bin_idx(i,j,k)) * Fc(i,j,k)
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for difference map calculation(s) = ', finish-start
  return
end subroutine differencemap

subroutine diffmap_norm(Fo,Fc,bin_idx,res_arr,smax,mode,nbin,nx,ny,nz, &
     diffmap)
  implicit none

  integer,   intent(in) :: mode,nbin,nx,ny,nz
  real,      intent(in) :: smax
  real,      dimension(0:nbin-1),intent(in) :: res_arr
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2), &
       intent(in) :: Fo, Fc
  integer,   dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),&
       intent(in) :: bin_idx
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,0:3),&
       intent(out) :: diffmap
  ! locals
  complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2) :: Eo, Ec
  integer,   dimension(3)        :: nxyz
  integer,   dimension(0:nbin-1) :: bin_arr_count
  real*8,    dimension(0:nbin-1) :: FoFc_sum,FcFc_sum,FoFo_sum
  real*8,    dimension(0:nbin-1) :: fsc
  !
  integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin,cbin
  real       :: start, finish
  logical    :: debug, make_all_zero
  real*8     :: d_tmp, d2_tmp

  !
  debug         = .FALSE.
  make_all_zero = .FALSE.
  if(mode == 1) debug = .TRUE.
  call cpu_time(start)

  bin_arr_count = 0
  fsc           = 0.0d0
  FoFo_sum      = 0.0d0
  FoFc_sum      = 0.0d0
  FcFc_sum      = 0.0d0
  diffmap       = dcmplx(0.0d0, 0.0d0)
  Eo            = dcmplx(0.0d0, 0.0d0)
  Ec            = dcmplx(0.0d0, 0.0d0)
  
  xyzmin = 0; xyzmax = 0

  nxyz = (/ nx, ny, nz /)

  xyzmin(1) = int(-nxyz(1)/2)
  xyzmin(2) = int(-nxyz(2)/2)
  xyzmin(3) = int(-nxyz(3)/2)
  xyzmax    = -(xyzmin+1)
  if(debug) print*, 'xyzmin = ', xyzmin
  if(debug) print*, 'xyzmax = ', xyzmax

  cbin =  minloc(sqrt((res_arr - smax)**2), DIM=1) - 1
  
  ! Not using Friedel's Law
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), xyzmax(3)
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle  
           bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
           FoFo_sum(bin_idx(i,j,k)) = FoFo_sum(bin_idx(i,j,k)) + &
                real(Fo(i,j,k) * conjg(Fo(i,j,k)))
           FoFc_sum(bin_idx(i,j,k)) = FoFc_sum(bin_idx(i,j,k)) + &
                real(Fo(i,j,k) * conjg(Fc(i,j,k)))
           FcFc_sum(bin_idx(i,j,k)) = FcFc_sum(bin_idx(i,j,k)) + &
                real(Fc(i,j,k) * conjg(Fc(i,j,k)))
        end do
     end do
  end do

  ! Calculate normalized structure factors
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), 0 !xyzmax(3)
           if(k == xyzmin(3) .or. j == xyzmin(2) .or. i == xyzmin(1)) cycle
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           if(FoFc_sum(bin_idx(i,j,k)) <= 0.0) cycle ! singal var cannot be negative
           Eo(i,j,k) = Fo(i,j,k)/sqrt(FoFo_sum(bin_idx(i,j,k)))
           Ec(i,j,k) = Fc(i,j,k)/sqrt(FcFc_sum(bin_idx(i,j,k)))
           Eo(-i,-j,-k) = conjg(Eo(i,j,k))
           Ec(-i,-j,-k) = conjg(Ec(i,j,k))
        end do
     end do
  end do

  ! Set FSC to zero after covar bacame less than zero
  if(debug) print*,'ibin  resol  FoFo   FcFc   FoFc   FSC    bincount'
  do ibin=0, nbin-1
     if(FoFc_sum(ibin) < 0.0) make_all_zero = .TRUE.
     if(make_all_zero) FoFc_sum(ibin) = 0.0
     if(.not.(make_all_zero) .and. (FoFc_sum(ibin) > 0.0)) then
        fsc(ibin) = FoFc_sum(ibin) / &
           (sqrt(FcFc_sum(ibin)) * sqrt(FoFo_sum(ibin)))
     end if
     if(debug)then
        print*,ibin, res_arr(ibin), FoFo_sum(ibin),FcFc_sum(ibin), &
               FoFc_sum(ibin),fsc(ibin),bin_arr_count(ibin) 
     end if    
  end do


  !if(debug) print*,'ibin  resol  FoFo   FcFc   FoFc   FSC    bincount'
  !do ibin=0, nbin-1
  !   if(FcFc_sum(ibin) == 0.0 .or. FoFo_sum(ibin) == 0.0) cycle
  !   !if(FoFc_sum(ibin) == 0.0) cycle
  !   fsc(ibin) = FoFc_sum(ibin) / (sqrt(FcFc_sum(ibin)) * sqrt(FoFo_sum(ibin)))
  !   print*,ibin, res_arr(ibin), FoFo_sum(ibin),FcFc_sum(ibin), FoFc_sum(ibin), &
  !        fsc(ibin),bin_arr_count(ibin)
  !end do

  
  ! Calculate difference map
  if(cbin > (nbin-1)) cbin = nbin-1
  print*
  print*, 'Resolution of Difference map: ', res_arr(cbin), ' A'
  print*, 'Calculating map. Please wait...'
  do i=xyzmin(1), xyzmax(1)
     do j=xyzmin(2), xyzmax(2)
        do k=xyzmin(3), xyzmax(3)
           if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
           if(bin_idx(i,j,k) > cbin) cycle
           if(fsc(bin_idx(i,j,k)) <= 0.0) cycle
           diffmap(i,j,k,0) = (Eo(i,j,k) - Ec(i,j,k)) * sqrt(fsc(bin_idx(i,j,k)))
           diffmap(i,j,k,1) = (Ec(i,j,k) - Eo(i,j,k)) * sqrt(fsc(bin_idx(i,j,k)))
           diffmap(i,j,k,2) = Eo(i,j,k) * sqrt(fsc(bin_idx(i,j,k)))
           diffmap(i,j,k,3) = Ec(i,j,k) * sqrt(fsc(bin_idx(i,j,k)))
        end do
     end do
  end do

  call cpu_time(finish)
  if(debug) print*, 'time for difference map calculation(s) = ', finish-start
  return
end subroutine diffmap_norm

subroutine scalesigmafval(hf1,hf2,Fc,bin_idx,res_arr,nbin,mode,nx,ny,nz,Fo, &
   scale_d,sigma,bin_noise_var,fval)
implicit none

integer, intent(in) :: mode,nbin,nx,ny,nz
real, dimension(0:nbin-1) :: res_arr
complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: hf1,hf2,Fc
integer, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx
complex*16, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(out) :: Fo
real*8, dimension(0:nbin-1),intent(out) :: bin_noise_var,scale_d,sigma
real*8, intent(out) :: fval
! locals
integer,   dimension(3)          :: nxyz
integer,   dimension(0:nbin-1) :: bin_arr_count
real*8,    dimension(0:nbin-1) :: bin_arr_fdiff,FoFc_sum,FcFc_sum,FomDFc_sum,FoFo_sum
real*8,    dimension(0:nbin-1) :: DFcmFo_sum
real*8,    dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2) :: term1
!
complex*16  :: fdiff
integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin
real       :: start, finish
logical    :: debug,make_all_zero
!
real*8       :: d_tmp,sigma_tmp,cc,term2!,fval

!
debug         = .FALSE.
make_all_zero = .FALSE.
if(mode == 1) debug = .TRUE.
call cpu_time(start)

Fo = dcmplx(0.0d0, 0.0d0)


bin_arr_count = 0
bin_arr_fdiff = 0.0d0
bin_noise_var = 0.0d0
FoFo_sum      = 0.0d0
FoFc_sum      = 0.0d0
FcFc_sum      = 0.0d0
FomDFc_sum    = 0.0d0
DFcmFo_sum    = 0.0d0
term2         = 0.0d0
term1         = 0.0d0

xyzmin = 0; xyzmax = 0

nxyz = (/ nx, ny, nz /)

xyzmin(1) = int(-nxyz(1)/2)
xyzmin(2) = int(-nxyz(2)/2)
xyzmin(3) = int(-nxyz(3)/2)
xyzmax    = -(xyzmin+1)
if(debug) print*, 'xyzmin = ', xyzmin
if(debug) print*, 'xyzmax = ', xyzmax

! Use only hemisphere data
do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
       do k=xyzmin(3), 0!xyzmax(3) 
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
         bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
         fdiff = hf1(i,j,k) - hf2(i,j,k)
         bin_arr_fdiff(bin_idx(i,j,k)) = bin_arr_fdiff(bin_idx(i,j,k)) + real(fdiff * conjg(fdiff))
         Fo(i,j,k) = (hf1(i,j,k) + hf2(i,j,k))/2.0
         FoFo_sum(bin_idx(i,j,k)) = FoFo_sum(bin_idx(i,j,k)) + real(Fo(i,j,k) * conjg(Fo(i,j,k)))
         FoFc_sum(bin_idx(i,j,k)) = FoFc_sum(bin_idx(i,j,k)) + real(Fo(i,j,k) * conjg(Fc(i,j,k)))
         FcFc_sum(bin_idx(i,j,k)) = FcFc_sum(bin_idx(i,j,k)) + real(Fc(i,j,k) * conjg(Fc(i,j,k)))
      end do
   end do
end do

! Estimate scale_d in resolution bins
do ibin=0, nbin-1 !to make compatible with python arrays
   bin_noise_var(ibin) = bin_arr_fdiff(ibin) / (bin_arr_count(ibin) * 4)
   d_tmp = FoFc_sum(ibin) / FcFc_sum(ibin)
   cc = FoFc_sum(ibin) / sqrt(FoFo_sum(ibin) * FcFc_sum(ibin))
   if(d_tmp < 0.0) make_all_zero = .TRUE.
   if(make_all_zero)then
      scale_d(ibin) = 0.0
   else
      scale_d(ibin) = d_tmp
   end if
end do

bin_arr_count = 0
do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
      do k=xyzmin(3), 0!xyzmax(3)
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
         bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
         ! To calculate unexplained signal in resol. bins
         DFcmFo_sum(bin_idx(i,j,k)) = DFcmFo_sum(bin_idx(i,j,k)) + &
              real((scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)) * &
              conjg(scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)))
         ! Calculate functional value (likelihood)
         term1(i,j,k) = real((scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)) * &
              conjg(scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k))) / &
              (sigma(bin_idx(i,j,k)) + bin_noise_var(bin_idx(i,j,k)))
         term2 = term2 + log10(sigma(bin_idx(i,j,k)) + bin_noise_var(bin_idx(i,j,k)))
      end do
   end do
end do

if(debug) print*,'ibin  resol   noise_var  unexpl.signal   scale_d    bin_reflex_count'
make_all_zero = .FALSE.
do ibin=0, nbin-1
   sigma_tmp = DFcmFo_sum(ibin) / bin_arr_count(ibin) - bin_noise_var(ibin)
   if(sigma_tmp < 0.0) make_all_zero = .TRUE.
   if(make_all_zero)then
      sigma(ibin) = 0.0
   else
      sigma(ibin) = sigma_tmp
   end if
   if(debug)then
      print*,ibin,res_arr(ibin),bin_noise_var(ibin),sigma(ibin),scale_d(ibin),&
           bin_arr_count(ibin)
   end if
end do

do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
      do k=xyzmin(3), 0!xyzmax(3)
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle     
         ! Calculate functional value (likelihood)
         term1(i,j,k) = real((scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)) * &
              conjg(scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k))) / &
              (sigma(bin_idx(i,j,k)) + bin_noise_var(bin_idx(i,j,k)))
         term2 = term2 + log10(sigma(bin_idx(i,j,k)) + bin_noise_var(bin_idx(i,j,k)))
      end do
   end do
end do
fval = sum(term1) + term2
if(debug) print*, 'Functional (likelihood) value = ', fval
call cpu_time(finish)
if(debug) print*, 'time for scale, sigma and functional calculation(s) = ', finish-start
return
end subroutine scalesigmafval

subroutine scalesigmafval_full(Fo,Fc,bin_idx,res_arr,mode,nbin,nx,ny,nz,scale_d, &
   sigma,totalvar,fval)
implicit none

integer,   intent(in) :: mode,nbin,nx,ny,nz
real,      dimension(0:nbin-1),intent(in) :: res_arr
complex*8, dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: Fo,Fc
integer,   dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2),intent(in) :: bin_idx
real*8,    dimension(0:nbin-1),intent(out) :: scale_d,sigma,totalvar
real*8,                        intent(out) :: fval
! locals
integer,   dimension(3)        :: nxyz
integer,   dimension(0:nbin-1) :: bin_arr_count
real*8,    dimension(0:nbin-1) :: FoFo_sum,FoFc_sum,FcFc_sum,FomDFc_sum,noisevar,DFcmFo_sum
real*8,    dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2) :: term1
!
integer    :: i,j,k,xyzmin(3),xyzmax(3),ibin
real       :: start, finish
logical    :: debug,make_all_zero
real*8     :: d_tmp,sigma_tmp, cc, term2

!
debug         = .FALSE.
make_all_zero = .FALSE.
if(mode == 1) debug = .TRUE.
call cpu_time(start)

bin_arr_count = 0
scale_d       = 0.0d0
sigma         = 0.0d0
totalvar      = 0.0d0
noisevar      = 0.0d0
FoFo_sum      = 0.0d0
FoFc_sum      = 0.0d0
FcFc_sum      = 0.0d0
FomDFc_sum    = 0.0d0
DFcmFo_sum    = 0.0d0
term2         = 0.0d0
term1         = 0.0d0

xyzmin = 0; xyzmax = 0

nxyz = (/ nx, ny, nz /)

xyzmin(1) = int(-nxyz(1)/2)
xyzmin(2) = int(-nxyz(2)/2)
xyzmin(3) = int(-nxyz(3)/2)
xyzmax    = -(xyzmin+1)
if(debug) print*, 'xyzmin = ', xyzmin
if(debug) print*, 'xyzmax = ', xyzmax

! Use only hemisphere data
do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
      do k=xyzmin(3), 0!xyzmax(3)
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle  
         bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
         FoFo_sum(bin_idx(i,j,k)) = FoFo_sum(bin_idx(i,j,k)) + real(Fo(i,j,k) * conjg(Fo(i,j,k)))
         FoFc_sum(bin_idx(i,j,k)) = FoFc_sum(bin_idx(i,j,k)) + real(Fo(i,j,k) * conjg(Fc(i,j,k)))
         FcFc_sum(bin_idx(i,j,k)) = FcFc_sum(bin_idx(i,j,k)) + real(Fc(i,j,k) * conjg(Fc(i,j,k)))
      end do
   end do
end do

! Estimate scale_d in resolution bins
do ibin=0, nbin-1
   d_tmp = FoFc_sum(ibin) / FcFc_sum(ibin)
   if(d_tmp < 0.0) make_all_zero = .TRUE.
   if(make_all_zero)then
      scale_d(ibin) = 0.0
   else
      scale_d(ibin) = d_tmp
   end if
end do

bin_arr_count = 0
do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
      do k=xyzmin(3), 0!xyzmax(3)
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
         bin_arr_count(bin_idx(i,j,k)) = bin_arr_count(bin_idx(i,j,k)) + 1
         ! To calculate unexplained signal in resol. bins
         DFcmFo_sum(bin_idx(i,j,k)) = DFcmFo_sum(bin_idx(i,j,k)) + &
              real((scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)) * &
              conjg(scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)))
      end do
   end do
end do

! Calculate unexplained signal in resol. bins
if(debug) print*,'ibin  resol  noise_var  unexp.signalvar   scale_d    bincount'
make_all_zero = .FALSE.
do ibin=0, nbin-1
   sigma_tmp = DFcmFo_sum(ibin) / bin_arr_count(ibin) - noisevar(ibin)
   if(sigma_tmp < 0.0) make_all_zero = .TRUE.
   if(make_all_zero)then
      sigma(ibin) = 0.0
   else
      sigma(ibin) = sigma_tmp
   end if
   if(debug)then
      print*,ibin,res_arr(ibin),noisevar(ibin),sigma(ibin),scale_d(ibin), &
           bin_arr_count(ibin)
   end if
end do

do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
      do k=xyzmin(3), 0!xyzmax(3)
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
         ! Calculate functional value (likelihood)
         term1(i,j,k) = real((scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k)) * &
              conjg(scale_d(bin_idx(i,j,k)) * Fc(i,j,k) - Fo(i,j,k))) / &
              (sigma(bin_idx(i,j,k)) + noisevar(bin_idx(i,j,k)))
         term2 = term2 + log10(sigma(bin_idx(i,j,k)) + noisevar(bin_idx(i,j,k)))
      end do
   end do
end do
fval = sum(term1) + term2
if(debug) print*, 'Functional (likelihood) value = ', fval
call cpu_time(finish)
if(debug) print*, 'time for scale, sigma and functional calculation(s) = ', finish-start
return
end subroutine scalesigmafval_full

subroutine tricubic_zoom(scale,F,FRS,nc,mode,nx,ny,nz)
implicit none
real,intent(in) :: scale
integer,intent(in):: nx,ny,nz,mode,nc
complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nc),intent(in):: F
complex*16,dimension(-nx/2:(nx-2)/2,-ny/2:(ny-2)/2,-nz/2:(nz-2)/2,nc),intent(out):: FRS
integer :: x1(3,-1:2)
real*8 :: x(3),xd(3),s(3)
real*8 :: ul,ul2
real*8 :: vl(4)
!real :: high_res,resol
integer :: i,j,h,k,l,nmin,ncopies
logical :: debug
integer :: nxyz(3),nxyzmn(3),nxyzmx(3)
integer :: xmin,xmax,ymin,ymax,zmin,zmax,ic
complex*16 :: fl(-1:2,-1:2,-1:2,nc),esz(-1:2,-1:2,nc),esyz(-1:2,nc)
complex*16 :: esxyz(nc)
!
FRS = dcmplx(0.0d0, 0.0d0)
x = 0.0d0
xd = 0.0d0
fl = dcmplx(0.0d0, 0.0d0)
esz = dcmplx(0.0d0, 0.0d0)
esyz = dcmplx(0.0d0, 0.0d0)
esxyz = dcmplx(0.0d0, 0.0d0)

debug = .FALSE.
if(mode == 1) debug = .TRUE.
!   Body
ncopies = nc
nxyz(1) = nx; nxyz(2) = ny; nxyz(3) =nz
nxyzmn(1) = -nx/2; nxyzmn(2) = -ny/2; nxyzmn(3) = -nz/2
nxyzmx(1) = (nx-2)/2; nxyzmx(2) = (ny-2)/2; nxyzmx(3) = (nz-2)/2
nmin = min(nx,ny,nz)

xmin = int(-nx/2); xmax = -(xmin+1)
ymin = int(-ny/2); ymax = -(ymin+1)
zmin = int(-nz/2); zmax = -(zmin+1)

if(debug) write(*,*) nxyz,nxyzmn,nxyzmx

do l = zmin, zmax
   do k = ymin, ymax
      do h = xmin, 0!xmax
         s(1) = h
         s(2) = k
         s(3) = l
         !x = matmul(transpose(RM),s)
         x = scale * s
         do i = 1, 3
            x1(i,0) = floor(x(i))
            xd(i) = x(i) - real(x1(i,0))
            if(abs(xd(i)).gt.1.0) then
               print*, 'Something is wrong ',xd(i)
               stop
            endif
            x1(i,1)  = x1(i,0) + 1
            x1(i,2)  = x1(i,0) + 2
            x1(i,-1) = x1(i,0) - 1
         end do
         !
         !  Careful here: we may get to the outside of the array
         do i = 1,3
            do j= -1,2
               x1(i,j) = min(nxyzmx(i),max(nxyzmn(i),x1(i,j)))
            enddo
         enddo
         do ic=1, ncopies
            fl(-1:2,-1:2,-1:2,ic) = F(x1(1,-1:2),x1(2,-1:2),x1(3,-1:2),ic)
         end do

         !
         !  Alternattive implementation
         !  along z
         ul = xd(3)
         ul2 = ul*ul
         vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
         vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
         vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
         vl(4) = ul2*(ul-1.0d0)
         vl = 0.5d0*vl
         do ic=1,ncopies
            do j=-1,2
               do i=-1,2
                  esz(i,j,ic) = dot_product(vl,fl(i,j,-1:2,ic))
               enddo
            enddo
         end do
         ul = xd(2)
         ul2 = ul*ul
         vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
         vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
         vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
         vl(4) = ul2*(ul-1.0d0)
         vl = 0.5d0*vl
         do ic=1,ncopies
            do i=-1,2
               esyz(i,ic) = dot_product(vl,esz(i,-1:2,ic))
            enddo
         end do
         ul = xd(1)
         ul2 = ul*ul
         vl(1) = ul*((2.0d0-ul)*ul-1.0d0)
         vl(2) = ul2*(3.0d0*ul-5.0d0)+2.0d0
         vl(3) = ul*((4.0d0-3.0d0*ul)*ul+1.0d0)
         vl(4) = ul2*(ul-1.0d0)
         vl = 0.5d0*vl
         do ic=1,ncopies
            esxyz(ic) =  dot_product(vl,esyz(-1:2,ic))
            FRS(h,k,l,ic) = esxyz(ic)
            if((h == xmin).or.(k == ymin).or.(l == zmin)) cycle
            FRS(-h,-k,-l,ic) = conjg(esxyz(ic))
         end do
      end do
   end do
end do
return
end subroutine tricubic_zoom

subroutine ll_derivatives(Fo,Fc,bin_idx,D,totalvar,mode,nbin,nx,ny,nz,dll,ddll)
implicit none
integer,                     intent(in) :: mode,nbin,nx,ny,nz
real*8,      dimension(0:nbin-1),intent(in) :: D,totalvar
complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: Fo,Fc
integer,   dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(in)  :: bin_idx
complex*16, dimension(-nx/2:(nx-2)/2, -ny/2:(ny-2)/2, -nz/2:(nz-2)/2),intent(out) :: dll,ddll
!
integer    :: i,j,k,xyzmin(3),xyzmax(3),nxyz(3)

xyzmin = 0; xyzmax = 0; nxyz = 0

dll = dcmplx(0.0d0, 0.0d0)
ddll = dcmplx(0.0d0, 0.0d0)

nxyz = (/ nx, ny, nz /)

xyzmin(1) = int(-nxyz(1)/2)
xyzmin(2) = int(-nxyz(2)/2)
xyzmin(3) = int(-nxyz(3)/2)
xyzmax    = -(xyzmin+1)
!print*, D

do i=xyzmin(1), xyzmax(1)
   do j=xyzmin(2), xyzmax(2)
      do k=xyzmin(3), xyzmax(3)
         if(bin_idx(i,j,k) < 0 .or. bin_idx(i,j,k) > nbin-1) cycle
         if(totalvar(bin_idx(i,j,k)) <= 0.0) cycle
         
         dll(i,j,k) = ((2.0 * D(bin_idx(i,j,k)) * &
              (D(bin_idx(i,j,k))*Fc(i,j,k)) - Fo(i,j,k))) / &
              (totalvar(bin_idx(i,j,k)))

         ddll(i,j,k) =  (2.0 * D(bin_idx(i,j,k)) *  D(bin_idx(i,j,k))) / &
              (totalvar(bin_idx(i,j,k)))
         
      end do
   end do
end do
return
end subroutine ll_derivatives
