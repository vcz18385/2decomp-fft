!! SPDX-License-Identifier: BSD-3-Clause

! This file contains the routines that transpose data from X to Y pencil
submodule(decomp_2d) d2d_transpose_x_to_y

   use decomp_2d_constants, only: mytype

   implicit none

contains

   module subroutine transpose_x_to_y_real_short_mpi_cpu(src, dst)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_x_to_y(src, dst, decomp_main)

   end subroutine transpose_x_to_y_real_short_mpi_cpu

   module subroutine transpose_x_to_y_real_long_mpi_cpu(src, dst, decomp)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      if (decomp_profiler_transpose) call decomp_profiler_start("transp_x_y_r")

      if (dims(1) == 1) then
         dst = src
      else
         call transpose_x_to_y_real_mpi_cpu(src, dst, decomp, work1_r, work2_r)
      end if

      if (decomp_profiler_transpose) call decomp_profiler_end("transp_x_y_r")

   end subroutine transpose_x_to_y_real_long_mpi_cpu

   subroutine transpose_x_to_y_real_mpi_cpu(src, dst, decomp, wk1, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(out) :: wk1, wk2

      integer :: s1, s2, s3, d1, d2, d3
      integer :: ierror

      s1 = SIZE(src, 1)
      s2 = SIZE(src, 2)
      s3 = SIZE(src, 3)
      d1 = SIZE(dst, 1)
      d2 = SIZE(dst, 2)
      d3 = SIZE(dst, 3)

      ! rearrange source array as send buffer
      call mem_split_xy_real_mpi_cpu(src, s1, s2, s3, wk1, dims(1), &
                             decomp%x1dist, decomp)

      ! define receive buffer
      ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
      call MPI_ALLTOALL(wk1, decomp%x1count, real_type, &
                        wk2, decomp%y1count, real_type, &
                        DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")

#else
      call MPI_ALLTOALLV(wk1, decomp%x1cnts, decomp%x1disp, real_type, &
                         wk2, decomp%y1cnts, decomp%y1disp, real_type, &
                         DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

      ! rearrange receive buffer
      call mem_merge_xy_real_mpi_cpu(wk2, d1, d2, d3, dst, dims(1), &
                             decomp%y1dist, decomp)

   end subroutine transpose_x_to_y_real_mpi_cpu

   module subroutine transpose_x_to_y_complex_short_mpi_cpu(src, dst)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_x_to_y(src, dst, decomp_main)

   end subroutine transpose_x_to_y_complex_short_mpi_cpu

   module subroutine transpose_x_to_y_complex_long_mpi_cpu(src, dst, decomp)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      if (decomp_profiler_transpose) call decomp_profiler_start("transp_x_y_c")

      if (dims(1) == 1) then
         dst = src
      else
         call transpose_x_to_y_complex_mpi_cpu(src, dst, decomp, work1_c, work2_c)
      end if

      if (decomp_profiler_transpose) call decomp_profiler_end("transp_x_y_c")

   end subroutine transpose_x_to_y_complex_long_mpi_cpu

   subroutine transpose_x_to_y_complex_mpi_cpu(src, dst, decomp, wk1, wk2)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      complex(mytype), dimension(:), intent(OUT) :: wk1, wk2

      integer :: s1, s2, s3, d1, d2, d3
      integer :: ierror

      s1 = SIZE(src, 1)
      s2 = SIZE(src, 2)
      s3 = SIZE(src, 3)
      d1 = SIZE(dst, 1)
      d2 = SIZE(dst, 2)
      d3 = SIZE(dst, 3)

      ! rearrange source array as send buffer
      call mem_split_xy_complex_mpi_cpu(src, s1, s2, s3, wk1, dims(1), &
                                decomp%x1dist, decomp)

      ! define receive buffer
      ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
      call MPI_ALLTOALL(wk1, decomp%x1count, complex_type, &
                        wk2, decomp%y1count, complex_type, &
                        DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")

#else
      call MPI_ALLTOALLV(wk1, decomp%x1cnts, decomp%x1disp, complex_type, &
                         wk2, decomp%y1cnts, decomp%y1disp, complex_type, &
                         DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

      ! rearrange receive buffer
      call mem_merge_xy_complex_mpi_cpu(wk2, d1, d2, d3, dst, dims(1), &
                                decomp%y1dist, decomp)

   end subroutine transpose_x_to_y_complex_mpi_cpu

   subroutine mem_split_xy_real_mpi_cpu(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      real(mytype), dimension(n1, n2, n3), intent(IN) :: in
      real(mytype), dimension(*), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%x1count + 1
#else
         pos = decomp%x1disp(m) + 1
#endif

         do k = 1, n3
            do j = 1, n2
               do i = i1, i2
                  out(pos) = in(i, j, k)
                  pos = pos + 1
               end do
            end do
         end do
      end do

   end subroutine mem_split_xy_real_mpi_cpu

   subroutine mem_split_xy_complex_mpi_cpu(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(n1, n2, n3), intent(IN) :: in
      complex(mytype), dimension(*), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%x1count + 1
#else
         pos = decomp%x1disp(m) + 1
#endif

         do k = 1, n3
            do j = 1, n2
               do i = i1, i2
                  out(pos) = in(i, j, k)
                  pos = pos + 1
               end do
            end do
         end do
      end do

   end subroutine mem_split_xy_complex_mpi_cpu

   subroutine mem_merge_xy_real_mpi_cpu(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      real(mytype), dimension(*), intent(IN) :: in
      real(mytype), dimension(n1, n2, n3), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%y1count + 1
#else
         pos = decomp%y1disp(m) + 1
#endif

         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  out(i, j, k) = in(pos)
                  pos = pos + 1
               end do
            end do
         end do
      end do

   end subroutine mem_merge_xy_real_mpi_cpu

   subroutine mem_merge_xy_complex_mpi_cpu(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(*), intent(IN) :: in
      complex(mytype), dimension(n1, n2, n3), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%y1count + 1
#else
         pos = decomp%y1disp(m) + 1
#endif

         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  out(i, j, k) = in(pos)
                  pos = pos + 1
               end do
            end do
         end do
      end do

   end subroutine mem_merge_xy_complex_mpi_cpu

end submodule d2d_transpose_x_to_y
