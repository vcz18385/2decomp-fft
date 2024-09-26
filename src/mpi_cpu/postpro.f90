!! SPDX-License-Identifier: BSD-3-Clause

! This is the main 2D pencil decomposition module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Coarser mesh support for statistic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_coarser_mesh_statS(i_skip, j_skip, k_skip, from1)

      implicit none

      integer, intent(IN) :: i_skip, j_skip, k_skip
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      integer, dimension(3) :: skip
      integer :: i
      logical, save :: first_call = .true.

      if (first_call) then
         first_call = .false.
         call decomp_2d_warning(__FILE__, __LINE__, 0, "init_coarser_mesh_stat is deprecated and will be removed")
      end if

      coarse_mesh_starts_from_1 = from1
      iskipS = i_skip
      jskipS = j_skip
      kskipS = k_skip

      skip(1) = iskipS
      skip(2) = jskipS
      skip(3) = kskipS

      do i = 1, 3
         if (from1) then
            xstS(i) = (xstart(i) + skip(i) - 1) / skip(i)
            if (mod(xstart(i) + skip(i) - 1, skip(i)) /= 0) xstS(i) = xstS(i) + 1
            xenS(i) = (xend(i) + skip(i) - 1) / skip(i)
         else
            xstS(i) = xstart(i) / skip(i)
            if (mod(xstart(i), skip(i)) /= 0) xstS(i) = xstS(i) + 1
            xenS(i) = xend(i) / skip(i)
         end if
         xszS(i) = xenS(i) - xstS(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            ystS(i) = (ystart(i) + skip(i) - 1) / skip(i)
            if (mod(ystart(i) + skip(i) - 1, skip(i)) /= 0) ystS(i) = ystS(i) + 1
            yenS(i) = (yend(i) + skip(i) - 1) / skip(i)
         else
            ystS(i) = ystart(i) / skip(i)
            if (mod(ystart(i), skip(i)) /= 0) ystS(i) = ystS(i) + 1
            yenS(i) = yend(i) / skip(i)
         end if
         yszS(i) = yenS(i) - ystS(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            zstS(i) = (zstart(i) + skip(i) - 1) / skip(i)
            if (mod(zstart(i) + skip(i) - 1, skip(i)) /= 0) zstS(i) = zstS(i) + 1
            zenS(i) = (zend(i) + skip(i) - 1) / skip(i)
         else
            zstS(i) = zstart(i) / skip(i)
            if (mod(zstart(i), skip(i)) /= 0) zstS(i) = zstS(i) + 1
            zenS(i) = zend(i) / skip(i)
         end if
         zszS(i) = zenS(i) - zstS(i) + 1
      end do

      return
   end subroutine init_coarser_mesh_statS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Coarser mesh support for visualization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_coarser_mesh_statV(i_skip, j_skip, k_skip, from1)

      implicit none

      integer, intent(IN) :: i_skip, j_skip, k_skip
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      integer, dimension(3) :: skip
      integer :: i
      logical, save :: first_call = .true.

      if (first_call) then
         first_call = .false.
         call decomp_2d_warning(__FILE__, __LINE__, 0, "init_coarser_mesh_stat is deprecated and will be removed")
      end if

      coarse_mesh_starts_from_1 = from1
      iskipV = i_skip
      jskipV = j_skip
      kskipV = k_skip

      skip(1) = iskipV
      skip(2) = jskipV
      skip(3) = kskipV

      do i = 1, 3
         if (from1) then
            xstV(i) = (xstart(i) + skip(i) - 1) / skip(i)
            if (mod(xstart(i) + skip(i) - 1, skip(i)) /= 0) xstV(i) = xstV(i) + 1
            xenV(i) = (xend(i) + skip(i) - 1) / skip(i)
         else
            xstV(i) = xstart(i) / skip(i)
            if (mod(xstart(i), skip(i)) /= 0) xstV(i) = xstV(i) + 1
            xenV(i) = xend(i) / skip(i)
         end if
         xszV(i) = xenV(i) - xstV(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            ystV(i) = (ystart(i) + skip(i) - 1) / skip(i)
            if (mod(ystart(i) + skip(i) - 1, skip(i)) /= 0) ystV(i) = ystV(i) + 1
            yenV(i) = (yend(i) + skip(i) - 1) / skip(i)
         else
            ystV(i) = ystart(i) / skip(i)
            if (mod(ystart(i), skip(i)) /= 0) ystV(i) = ystV(i) + 1
            yenV(i) = yend(i) / skip(i)
         end if
         yszV(i) = yenV(i) - ystV(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            zstV(i) = (zstart(i) + skip(i) - 1) / skip(i)
            if (mod(zstart(i) + skip(i) - 1, skip(i)) /= 0) zstV(i) = zstV(i) + 1
            zenV(i) = (zend(i) + skip(i) - 1) / skip(i)
         else
            zstV(i) = zstart(i) / skip(i)
            if (mod(zstart(i), skip(i)) /= 0) zstV(i) = zstV(i) + 1
            zenV(i) = zend(i) / skip(i)
         end if
         zszV(i) = zenV(i) - zstV(i) + 1
      end do

      return
   end subroutine init_coarser_mesh_statV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Coarser mesh support for probe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_coarser_mesh_statP(i_skip, j_skip, k_skip, from1)

      implicit none

      integer, intent(IN) :: i_skip, j_skip, k_skip
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      integer, dimension(3) :: skip
      integer :: i
      logical, save :: first_call = .true.

      if (first_call) then
         first_call = .false.
         call decomp_2d_warning(__FILE__, __LINE__, 0, "init_coarser_mesh_stat is deprecated and will be removed")
      end if

      coarse_mesh_starts_from_1 = from1
      iskipP = i_skip
      jskipP = j_skip
      kskipP = k_skip

      skip(1) = iskipP
      skip(2) = jskipP
      skip(3) = kskipP

      do i = 1, 3
         if (from1) then
            xstP(i) = (xstart(i) + skip(i) - 1) / skip(i)
            if (mod(xstart(i) + skip(i) - 1, skip(i)) /= 0) xstP(i) = xstP(i) + 1
            xenP(i) = (xend(i) + skip(i) - 1) / skip(i)
         else
            xstP(i) = xstart(i) / skip(i)
            if (mod(xstart(i), skip(i)) /= 0) xstP(i) = xstP(i) + 1
            xenP(i) = xend(i) / skip(i)
         end if
         xszP(i) = xenP(i) - xstP(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            ystP(i) = (ystart(i) + skip(i) - 1) / skip(i)
            if (mod(ystart(i) + skip(i) - 1, skip(i)) /= 0) ystP(i) = ystP(i) + 1
            yenP(i) = (yend(i) + skip(i) - 1) / skip(i)
         else
            ystP(i) = ystart(i) / skip(i)
            if (mod(ystart(i), skip(i)) /= 0) ystP(i) = ystP(i) + 1
            yenP(i) = yend(i) / skip(i)
         end if
         yszP(i) = yenP(i) - ystP(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            zstP(i) = (zstart(i) + skip(i) - 1) / skip(i)
            if (mod(zstart(i) + skip(i) - 1, skip(i)) /= 0) zstP(i) = zstP(i) + 1
            zenP(i) = (zend(i) + skip(i) - 1) / skip(i)
         else
            zstP(i) = zstart(i) / skip(i)
            if (mod(zstart(i), skip(i)) /= 0) zstP(i) = zstP(i) + 1
            zenP(i) = zend(i) / skip(i)
         end if
         zszP(i) = zenP(i) - zstP(i) + 1
      end do

      return
   end subroutine init_coarser_mesh_statP

   ! Copy data from a fine-resolution array to a coarse one for statistic
   subroutine fine_to_coarseS(ipencil, var_fine, var_coarse)

      implicit none

      real(mytype), dimension(:, :, :) :: var_fine
      real(mytype), dimension(:, :, :) :: var_coarse
      integer, intent(IN) :: ipencil

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer :: i, j, k
      logical, save :: first_call = .true.

      if (first_call) then
         first_call = .false.
         call decomp_2d_warning(__FILE__, __LINE__, 0, "fine_to_coarse is deprecated and will be removed")
      end if

      if (ipencil == 1) then
         allocate (wk(xstS(1):xenS(1), xstS(2):xenS(2), xstS(3):xenS(3)))
         call alloc_x(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = xstS(3), xenS(3)
               do j = xstS(2), xenS(2)
                  do i = xstS(1), xenS(1)
                     wk(i, j, k) = wk2((i - 1) * iskipS + 1, (j - 1) * jskipS + 1, (k - 1) * kskipS + 1)
                  end do
               end do
            end do
         else
            do k = xstS(3), xenS(3)
               do j = xstS(2), xenS(2)
                  do i = xstS(1), xenS(1)
                     wk(i, j, k) = wk2(i * iskipS, j * jskipS, k * kskipS)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 2) then
         allocate (wk(ystS(1):yenS(1), ystS(2):yenS(2), ystS(3):yenS(3)))
         call alloc_y(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = ystS(3), yenS(3)
               do j = ystS(2), yenS(2)
                  do i = ystS(1), yenS(1)
                     wk(i, j, k) = wk2((i - 1) * iskipS + 1, (j - 1) * jskipS + 1, (k - 1) * kskipS + 1)
                  end do
               end do
            end do
         else
            do k = ystS(3), yenS(3)
               do j = ystS(2), yenS(2)
                  do i = ystS(1), yenS(1)
                     wk(i, j, k) = wk2(i * iskipS, j * jskipS, k * kskipS)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 3) then
         allocate (wk(zstS(1):zenS(1), zstS(2):zenS(2), zstS(3):zenS(3)))
         call alloc_z(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = zstS(3), zenS(3)
               do j = zstS(2), zenS(2)
                  do i = zstS(1), zenS(1)
                     wk(i, j, k) = wk2((i - 1) * iskipS + 1, (j - 1) * jskipS + 1, (k - 1) * kskipS + 1)
                  end do
               end do
            end do
         else
            do k = zstS(3), zenS(3)
               do j = zstS(2), zenS(2)
                  do i = zstS(1), zenS(1)
                     wk(i, j, k) = wk2(i * iskipS, j * jskipS, k * kskipS)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      end if

      deallocate (wk, wk2)

      return
   end subroutine fine_to_coarseS

   ! Copy data from a fine-resolution array to a coarse one for visualization
   subroutine fine_to_coarseV(ipencil, var_fine, var_coarse)

      implicit none

      real(mytype), dimension(:, :, :) :: var_fine
      real(mytype), dimension(:, :, :) :: var_coarse
      integer, intent(IN) :: ipencil

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer :: i, j, k
      logical, save :: first_call = .true.

      if (first_call) then
         first_call = .false.
         call decomp_2d_warning(__FILE__, __LINE__, 0, "fine_to_coarse is deprecated and will be removed")
      end if

      if (ipencil == 1) then
         allocate (wk(xstV(1):xenV(1), xstV(2):xenV(2), xstV(3):xenV(3)))
         call alloc_x(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = xstV(3), xenV(3)
               do j = xstV(2), xenV(2)
                  do i = xstV(1), xenV(1)
                     wk(i, j, k) = wk2((i - 1) * iskipV + 1, (j - 1) * jskipV + 1, (k - 1) * kskipV + 1)
                  end do
               end do
            end do
         else
            do k = xstV(3), xenV(3)
               do j = xstV(2), xenV(2)
                  do i = xstV(1), xenV(1)
                     wk(i, j, k) = wk2(i * iskipV, j * jskipV, k * kskipV)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 2) then
         allocate (wk(ystV(1):yenV(1), ystV(2):yenV(2), ystV(3):yenV(3)))
         call alloc_y(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = ystV(3), yenV(3)
               do j = ystV(2), yenV(2)
                  do i = ystV(1), yenV(1)
                     wk(i, j, k) = wk2((i - 1) * iskipV + 1, (j - 1) * jskipV + 1, (k - 1) * kskipV + 1)
                  end do
               end do
            end do
         else
            do k = ystV(3), yenV(3)
               do j = ystV(2), yenV(2)
                  do i = ystV(1), yenV(1)
                     wk(i, j, k) = wk2(i * iskipV, j * jskipV, k * kskipV)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 3) then
         allocate (wk(zstV(1):zenV(1), zstV(2):zenV(2), zstV(3):zenV(3)))
         call alloc_z(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = zstV(3), zenV(3)
               do j = zstV(2), zenV(2)
                  do i = zstV(1), zenV(1)
                     wk(i, j, k) = wk2((i - 1) * iskipV + 1, (j - 1) * jskipV + 1, (k - 1) * kskipV + 1)
                  end do
               end do
            end do
         else
            do k = zstV(3), zenV(3)
               do j = zstV(2), zenV(2)
                  do i = zstV(1), zenV(1)
                     wk(i, j, k) = wk2(i * iskipV, j * jskipV, k * kskipV)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      end if

      deallocate (wk, wk2)

      return
   end subroutine fine_to_coarseV

   ! Copy data from a fine-resolution array to a coarse one for probe
   subroutine fine_to_coarseP(ipencil, var_fine, var_coarse)

      implicit none

      real(mytype), dimension(:, :, :) :: var_fine
      real(mytype), dimension(:, :, :) :: var_coarse
      integer, intent(IN) :: ipencil

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer :: i, j, k
      logical, save :: first_call = .true.

      if (first_call) then
         first_call = .false.
         call decomp_2d_warning(__FILE__, __LINE__, 0, "fine_to_coarse is deprecated and will be removed")
      end if

      if (ipencil == 1) then
         allocate (wk(xstP(1):xenP(1), xstP(2):xenP(2), xstP(3):xenP(3)))
         call alloc_x(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = xstP(3), xenP(3)
               do j = xstP(2), xenP(2)
                  do i = xstP(1), xenP(1)
                     wk(i, j, k) = wk2((i - 1) * iskipP + 1, (j - 1) * jskipP + 1, (k - 1) * kskipP + 1)
                  end do
               end do
            end do
         else
            do k = xstP(3), xenP(3)
               do j = xstP(2), xenP(2)
                  do i = xstP(1), xenP(1)
                     wk(i, j, k) = wk2(i * iskipP, j * jskipP, k * kskipP)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 2) then
         allocate (wk(ystP(1):yenP(1), ystP(2):yenP(2), ystP(3):yenP(3)))
         call alloc_y(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = ystP(3), yenP(3)
               do j = ystP(2), yenP(2)
                  do i = ystP(1), yenP(1)
                     wk(i, j, k) = wk2((i - 1) * iskipP + 1, (j - 1) * jskipP + 1, (k - 1) * kskipP + 1)
                  end do
               end do
            end do
         else
            do k = ystP(3), yenP(3)
               do j = ystP(2), yenP(2)
                  do i = ystP(1), yenP(1)
                     wk(i, j, k) = wk2(i * iskipP, j * jskipP, k * kskipP)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 3) then
         allocate (wk(zstP(1):zenP(1), zstP(2):zenP(2), zstP(3):zenP(3)))
         call alloc_z(wk2, opt_global=.true.)
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = zstP(3), zenP(3)
               do j = zstP(2), zenP(2)
                  do i = zstP(1), zenP(1)
                     wk(i, j, k) = wk2((i - 1) * iskipP + 1, (j - 1) * jskipP + 1, (k - 1) * kskipP + 1)
                  end do
               end do
            end do
         else
            do k = zstP(3), zenP(3)
               do j = zstP(2), zenP(2)
                  do i = zstP(1), zenP(1)
                     wk(i, j, k) = wk2(i * iskipP, j * jskipP, k * kskipP)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      end if

      deallocate (wk, wk2)

      return
   end subroutine fine_to_coarseP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Find sub-domain information held by current processor
   !   INPUT:
   !     nx, ny, nz - global data dimension
   !     pdim(3)    - number of processor grid in each dimension,
   !                  valid values: 1 - distibute locally;
   !                                2 - distribute across p_row;
   !                                3 - distribute across p_col
   !   OUTPUT:
   !     lstart(3)  - starting index
   !     lend(3)    - ending index
   !     lsize(3)   - size of the sub-block (redundant)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

      implicit none

      integer, intent(IN) :: nx, ny, nz
      integer, dimension(3), intent(IN) :: pdim
      integer, dimension(3), intent(OUT) :: lstart, lend, lsize

      integer, allocatable, dimension(:) :: st, en, sz
      integer :: i, gsize

      do i = 1, 3

         if (i == 1) then
            gsize = nx
         else if (i == 2) then
            gsize = ny
         else if (i == 3) then
            gsize = nz
         end if

         if (pdim(i) == 1) then        ! all local
            lstart(i) = 1
            lend(i) = gsize
            lsize(i) = gsize
         elseif (pdim(i) == 2) then    ! distribute across dims(1)
            allocate (st(0:dims(1) - 1))
            allocate (en(0:dims(1) - 1))
            allocate (sz(0:dims(1) - 1))
            call distribute(gsize, dims(1), st, en, sz)
            lstart(i) = st(coord(1))
            lend(i) = en(coord(1))
            lsize(i) = sz(coord(1))
            deallocate (st, en, sz)
         elseif (pdim(i) == 3) then    ! distribute across dims(2)
            allocate (st(0:dims(2) - 1))
            allocate (en(0:dims(2) - 1))
            allocate (sz(0:dims(2) - 1))
            call distribute(gsize, dims(2), st, en, sz)
            lstart(i) = st(coord(2))
            lend(i) = en(coord(2))
            lsize(i) = sz(coord(2))
            deallocate (st, en, sz)
         end if

      end do
      return

   end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   - distibutes grid points in one dimension
   !   - handles uneven distribution properly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine distribute(data1, proc, st, en, sz)

      implicit none
      ! data1 -- data size in any dimension to be partitioned
      ! proc  -- number of processors in that dimension
      ! st    -- array of starting index
      ! en    -- array of ending index
      ! sz    -- array of local size  (redundent)
      integer data1, proc, st(0:proc - 1), en(0:proc - 1), sz(0:proc - 1)
      integer i, size1, nl, nu

      size1 = data1 / proc
      nu = data1 - size1 * proc
      nl = proc - nu
      st(0) = 1
      sz(0) = size1
      en(0) = size1
      do i = 1, nl - 1
         st(i) = st(i - 1) + size1
         sz(i) = size1
         en(i) = en(i - 1) + size1
      end do
      size1 = size1 + 1
      do i = nl, proc - 1
         st(i) = en(i - 1) + 1
         sz(i) = size1
         en(i) = en(i - 1) + size1
      end do

      ! Safety checks
      if (en(proc - 1) /= data1) &
         call decomp_2d_abort(__FILE__, __LINE__, en(proc - 1), "Invalid distribution.")
      if (sz(proc - 1) /= (data1 - st(proc - 1) + 1)) &
         call decomp_2d_abort(__FILE__, __LINE__, sz(proc - 1), "Invalid distribution.")

   end subroutine distribute
