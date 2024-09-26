!! SPDX-License-Identifier: BSD-3-Clause

! This is the main 2D pencil decomposition module

module decomp_2d

   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use factor
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_profiler

   implicit none

   private        ! Make everything private unless declared public

   ! Default parameter opt_global in alloc subroutines
   logical, parameter :: default_opt_global = .false.

   ! some key global variables
   integer, save, public :: nx_global, ny_global, nz_global  ! global size

   ! parameters for 2D Cartesian topology
   integer, save, dimension(2) :: dims, coord
   integer, save, public :: DECOMP_2D_COMM_CART_X = MPI_COMM_NULL
   integer, save, public :: DECOMP_2D_COMM_CART_Y = MPI_COMM_NULL
   integer, save, public :: DECOMP_2D_COMM_CART_Z = MPI_COMM_NULL
   integer, save :: DECOMP_2D_COMM_ROW = MPI_COMM_NULL
   integer, save :: DECOMP_2D_COMM_COL = MPI_COMM_NULL

   ! define neighboring blocks (to be used in halo-cell support)
   !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
   ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom
   integer, save, dimension(3, 6) :: neighbour

   ! flags for periodic condition in three dimensions
   logical, save :: periodic_x, periodic_y, periodic_z

   !
   ! Output for the log can be changed by the external code before calling decomp_2d_init
   !
   !    0 => No log output
   !    1 => Master rank log output to stdout
   !    2 => Master rank log output to the file "decomp_2d_setup.log"
   !    3 => All ranks log output to a dedicated file
   !
   ! The default value is 2 (3 for debug builds)
   !
#ifdef DEBUG
   integer, public, save :: decomp_log = D2D_LOG_TOFILE_FULL
#else
   integer, public, save :: decomp_log = D2D_LOG_TOFILE
#endif

   !
   ! Debug level can be changed by the external code before calling decomp_2d_init
   !
   ! The environment variable "DECOMP_2D_DEBUG" can be used to change the debug level
   !
   ! Debug checks are performed only when the preprocessor variable DEBUG is defined
   !
#ifdef DEBUG
   integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_INFO
#else
   integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_OFF
#endif

   ! derived type to store decomposition info for a given global data size
   TYPE, public :: DECOMP_INFO
      ! staring/ending index and size of data held by current processor
      integer, dimension(3) :: xst, xen, xsz  ! x-pencil
      integer, dimension(3) :: yst, yen, ysz  ! y-pencil
      integer, dimension(3) :: zst, zen, zsz  ! z-pencil

      ! in addition to local information, processors also need to know
      ! some global information for global communications to work

      ! how each dimension is distributed along pencils
      integer, allocatable, dimension(:) :: &
         x1dist, y1dist, y2dist, z2dist

      ! send/receive buffer counts and displacements for MPI_ALLTOALLV
      integer, allocatable, dimension(:) :: &
         x1cnts, y1cnts, y2cnts, z2cnts
      integer, allocatable, dimension(:) :: &
         x1disp, y1disp, y2disp, z2disp

#ifdef EVEN
      ! buffer counts for MPI_ALLTOALL for padded-alltoall
      integer :: x1count, y1count, y2count, z2count
      ! evenly distributed data
      logical :: even
#endif

   END TYPE DECOMP_INFO

   ! main (default) decomposition information for global size nx*ny*nz
   TYPE(DECOMP_INFO), target, save, public :: decomp_main

   ! staring/ending index and size of data held by current processor
   ! duplicate 'decomp_main', needed by apps to define data structure
   integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
   integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
   integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

   ! These are the buffers used by MPI_ALLTOALL(V) calls
   integer, save :: decomp_buf_size = 0
   ! Shared real/complex buffers
   real(mytype), target, allocatable, dimension(:) :: work1, work2
   ! Real/complex pointers to CPU buffers
   real(mytype), pointer, contiguous, dimension(:) :: work1_r, work2_r
   complex(mytype), pointer, contiguous, dimension(:) :: work1_c, work2_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! To define smaller arrays using every several mesh points
   integer, save, dimension(3), public :: xszS, yszS, zszS, xstS, ystS, zstS, xenS, yenS, zenS
   integer, save, dimension(3), public :: xszV, yszV, zszV, xstV, ystV, zstV, xenV, yenV, zenV
   integer, save, dimension(3), public :: xszP, yszP, zszP, xstP, ystP, zstP, xenP, yenP, zenP
   logical, save :: coarse_mesh_starts_from_1
   integer, save :: iskipS, jskipS, kskipS
   integer, save :: iskipV, jskipV, kskipV
   integer, save :: iskipP, jskipP, kskipP

   ! public user routines
   public :: decomp_2d_init, decomp_2d_finalize, &
             transpose_x_to_y, transpose_y_to_z, &
             transpose_z_to_y, transpose_y_to_x, &
             decomp_info_init, decomp_info_finalize, partition, &
             decomp_info_print, &
             init_coarser_mesh_statS, fine_to_coarseS, &
             init_coarser_mesh_statV, fine_to_coarseV, &
             init_coarser_mesh_statP, fine_to_coarseP, &
             alloc_x, alloc_y, alloc_z, &
             update_halo, &
             get_decomp_info, &
             get_decomp_dims, &
             d2d_log_is_active, &
             d2d_log_get_unit, &
             d2d_log_close_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! These are routines to perform global data transpositions
   !
   !   Four combinations are available, enough to cover all situations
   !    - transpose_x_to_y (X-pencil --> Y-pencil)
   !    - transpose_y_to_z (Y-pencil --> Z-pencil)
   !    - transpose_z_to_y (Z-pencil --> Y-pencil)
   !    - transpose_y_to_x (Y-pencil --> X-pencil)
   !
   !   Generic interface provided here to support multiple data types
   !    - real and complex types supported through generic interface
   !    - single/double precision supported through pre-processing
   !       * see 'mytype' variable at the beginning
   !    - an optional argument can be supplied to transpose data whose
   !      global size is not the default nx*ny*nz
   !       * as the case in fft r2c/c2r interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface decomp_2d_init
      module procedure decomp_2d_init_ref
   end interface decomp_2d_init

   interface decomp_2d_finalize
      module procedure decomp_2d_finalize_ref
   end interface decomp_2d_finalize

   interface transpose_x_to_y
#if defined(_CPU)
      module subroutine transpose_x_to_y_real_long_mpi_cpu(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_x_to_y_real_long_mpi_cpu
      module subroutine transpose_x_to_y_real_short_mpi_cpu(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_x_to_y_real_short_mpi_cpu
      module subroutine transpose_x_to_y_complex_long_mpi_cpu(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_x_to_y_complex_long_mpi_cpu
      module subroutine transpose_x_to_y_complex_short_mpi_cpu(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_x_to_y_complex_short_mpi_cpu
#endif
#if defined(_GPU)
      module subroutine transpose_x_to_y_real_long_mpi_gpu(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_x_to_y_real_long_mpi_gpu
      module subroutine transpose_x_to_y_real_short_mpi_gpu(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_x_to_y_real_short_mpi_gpu
      module subroutine transpose_x_to_y_complex_long_mpi_gpu(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_x_to_y_complex_long_mpi_gpu
      module subroutine transpose_x_to_y_complex_short_mpi_gpu(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_x_to_y_complex_short_mpi_gpu
#endif
   end interface transpose_x_to_y


   interface transpose_y_to_z
      module subroutine transpose_y_to_z_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_z_real_long
      module subroutine transpose_y_to_z_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_z_real_short
      module subroutine transpose_y_to_z_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_z_complex_long
      module subroutine transpose_y_to_z_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_z_complex_short
   end interface transpose_y_to_z

   interface transpose_z_to_y
      module subroutine transpose_z_to_y_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_z_to_y_real_long
      module subroutine transpose_z_to_y_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_z_to_y_real_short
      module subroutine transpose_z_to_y_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_z_to_y_complex_long
      module subroutine transpose_z_to_y_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_z_to_y_complex_short
   end interface transpose_z_to_y

   interface transpose_y_to_x
      module subroutine transpose_y_to_x_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_x_real_long
      module subroutine transpose_y_to_x_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_x_real_short
      module subroutine transpose_y_to_x_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_x_complex_long
      module subroutine transpose_y_to_x_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_x_complex_short
   end interface transpose_y_to_x

   interface update_halo
      module procedure update_halo_real
      module procedure update_halo_real_short
      module procedure update_halo_complex
      module procedure update_halo_complex_short
   end interface update_halo

   interface alloc_x
      module procedure alloc_x_freal
      module procedure alloc_x_freal_short
      module procedure alloc_x_dreal
      module procedure alloc_x_dreal_short
      module procedure alloc_x_fcplx
      module procedure alloc_x_fcplx_short
      module procedure alloc_x_dcplx
      module procedure alloc_x_dcplx_short
      module procedure alloc_x_ints
      module procedure alloc_x_ints_short
      module procedure alloc_x_logs
      module procedure alloc_x_logs_short
   end interface alloc_x

   interface alloc_y
      module procedure alloc_y_freal
      module procedure alloc_y_freal_short
      module procedure alloc_y_dreal
      module procedure alloc_y_dreal_short
      module procedure alloc_y_fcplx
      module procedure alloc_y_fcplx_short
      module procedure alloc_y_dcplx
      module procedure alloc_y_dcplx_short
      module procedure alloc_y_ints
      module procedure alloc_y_ints_short
      module procedure alloc_y_logs
      module procedure alloc_y_logs_short
   end interface alloc_y

   interface alloc_z
      module procedure alloc_z_freal
      module procedure alloc_z_freal_short
      module procedure alloc_z_dreal
      module procedure alloc_z_dreal_short
      module procedure alloc_z_fcplx
      module procedure alloc_z_fcplx_short
      module procedure alloc_z_dcplx
      module procedure alloc_z_dcplx_short
      module procedure alloc_z_ints
      module procedure alloc_z_ints_short
      module procedure alloc_z_logs
      module procedure alloc_z_logs_short
   end interface alloc_z

   interface

      module function d2d_log_is_active()
         logical :: d2d_log_is_active
      end function d2d_log_is_active

      module function d2d_log_get_unit()
         integer :: d2d_log_get_unit
      end function d2d_log_get_unit

      module subroutine d2d_log_close_unit(io_unit)
         integer, intent(in) :: io_unit
      end subroutine d2d_log_close_unit

      module subroutine d2d_log(given_io_unit)
         integer, intent(in), optional :: given_io_unit
      end subroutine d2d_log

      module subroutine decomp_info_print(d2d, io_unit, d2dname)
         type(decomp_info), intent(in) :: d2d
         integer, intent(in) :: io_unit
         character(len=*), intent(in) :: d2dname
      end subroutine decomp_info_print

   end interface

contains

#include "decomp_2d_init_fin.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Return the default decomposition object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! FIXME avoid a copy and return a pointer to decomp_main
   ! TODO list the external codes using this subroutine
   subroutine get_decomp_info(decomp)

      implicit none

      ! FIXME TYPE(DECOMP_INFO), pointer :: decomp
      TYPE(DECOMP_INFO), intent(OUT) :: decomp

      ! FIXME decomp => decomp_main
      decomp = decomp_main

      return
   end subroutine get_decomp_info

   !
   ! Return the 2D processor grid
   !
   function get_decomp_dims()

      implicit none

      integer, dimension(2) :: get_decomp_dims

      get_decomp_dims = dims

   end function get_decomp_dims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advanced Interface allowing applications to define globle domain of
   ! any size, distribute it, and then transpose data among pencils.
   !  - generate 2D decomposition details as defined in DECOMP_INFO
   !  - the default global data size is nx*ny*nz
   !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
   !  - multiple global sizes can co-exist in one application, each
   !    using its own DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_info_init(nx, ny, nz, decomp)

      use, intrinsic:: iso_c_binding, only: c_f_pointer, c_loc

      implicit none

      integer, intent(IN) :: nx, ny, nz
      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      integer :: buf_size, status, errorcode

      ! verify the global size can actually be distributed as pencils
      if (nx_global < dims(1) .or. ny_global < dims(1) .or. ny_global < dims(2) .or. nz_global < dims(2)) then
         errorcode = 6
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Invalid 2D processor grid. '// &
                              'Make sure that min(nx,ny) >= p_row and '// &
                              'min(ny,nz) >= p_col')
      end if

      ! distribute mesh points
      allocate (decomp%x1dist(0:dims(1) - 1), decomp%y1dist(0:dims(1) - 1), &
                decomp%y2dist(0:dims(2) - 1), decomp%z2dist(0:dims(2) - 1))
      call get_dist(nx, ny, nz, decomp)

      ! generate partition information - starting/ending index etc.
      call partition(nx, ny, nz, (/1, 2, 3/), &
                     decomp%xst, decomp%xen, decomp%xsz)
      call partition(nx, ny, nz, (/2, 1, 3/), &
                     decomp%yst, decomp%yen, decomp%ysz)
      call partition(nx, ny, nz, (/2, 3, 1/), &
                     decomp%zst, decomp%zen, decomp%zsz)

      ! prepare send/receive buffer displacement and count for ALLTOALL(V)
      allocate (decomp%x1cnts(0:dims(1) - 1), decomp%y1cnts(0:dims(1) - 1), &
                decomp%y2cnts(0:dims(2) - 1), decomp%z2cnts(0:dims(2) - 1))
      allocate (decomp%x1disp(0:dims(1) - 1), decomp%y1disp(0:dims(1) - 1), &
                decomp%y2disp(0:dims(2) - 1), decomp%z2disp(0:dims(2) - 1))
      call prepare_buffer(decomp)

      ! allocate memory for the MPI_ALLTOALL(V) buffers
      ! define the buffers globally for performance reason

      buf_size = max(decomp%xsz(1) * decomp%xsz(2) * decomp%xsz(3), &
                     max(decomp%ysz(1) * decomp%ysz(2) * decomp%ysz(3), &
                         decomp%zsz(1) * decomp%zsz(2) * decomp%zsz(3)))

#ifdef EVEN
      ! padded alltoall optimisation may need larger buffer space
      buf_size = max(buf_size, &
                     max(decomp%x1count * dims(1), decomp%y2count * dims(2)))
      ! evenly distributed data ?
      if (mod(nx, dims(1)) == 0 .and. mod(ny, dims(1)) == 0 .and. &
          mod(ny, dims(2)) == 0 .and. mod(nz, dims(2)) == 0) then
         decomp%even = .true.
      else
         decomp%even = .false.
      end if
#endif

      ! check if additional memory is required
      if (buf_size > decomp_buf_size) then
         decomp_buf_size = buf_size
         if (associated(work1_r)) nullify (work1_r)
         if (associated(work2_r)) nullify (work2_r)
         if (associated(work1_c)) nullify (work1_c)
         if (associated(work2_c)) nullify (work2_c)
         if (allocated(work1)) deallocate (work1)
         if (allocated(work2)) deallocate (work2)
         allocate (work1(2 * buf_size), STAT=status)
         if (status /= 0) then
            errorcode = 2
            call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                 'Out of memory when allocating 2DECOMP workspace')
         end if
         allocate (work2(2 * buf_size), STAT=status)
         if (status /= 0) then
            errorcode = 2
            call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                 'Out of memory when allocating 2DECOMP workspace')
         end if
         call c_f_pointer(c_loc(work1), work1_r, [buf_size])
         call c_f_pointer(c_loc(work2), work2_r, [buf_size])
         call c_f_pointer(c_loc(work1), work1_c, [buf_size])
         call c_f_pointer(c_loc(work2), work2_c, [buf_size])
      end if

   end subroutine decomp_info_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Release memory associated with a DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_info_finalize(decomp)

      implicit none

      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      if (allocated(decomp%x1dist)) deallocate (decomp%x1dist)
      if (allocated(decomp%y1dist)) deallocate (decomp%y1dist)
      if (allocated(decomp%y2dist)) deallocate (decomp%y2dist)
      if (allocated(decomp%z2dist)) deallocate (decomp%z2dist)
      if (allocated(decomp%x1cnts)) deallocate (decomp%x1cnts)
      if (allocated(decomp%y1cnts)) deallocate (decomp%y1cnts)
      if (allocated(decomp%y2cnts)) deallocate (decomp%y2cnts)
      if (allocated(decomp%z2cnts)) deallocate (decomp%z2cnts)
      if (allocated(decomp%x1disp)) deallocate (decomp%x1disp)
      if (allocated(decomp%y1disp)) deallocate (decomp%y1disp)
      if (allocated(decomp%y2disp)) deallocate (decomp%y2disp)
      if (allocated(decomp%z2disp)) deallocate (decomp%z2disp)

      return
   end subroutine decomp_info_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  Define how each dimension is distributed across processors
   !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
   !    such global information is required locally at MPI_ALLTOALLV time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_dist(nx, ny, nz, decomp)

      implicit none

      integer, intent(IN) :: nx, ny, nz
      TYPE(DECOMP_INFO), intent(INOUT) :: decomp
      integer, allocatable, dimension(:) :: st, en

      allocate (st(0:dims(1) - 1))
      allocate (en(0:dims(1) - 1))
      call distribute(nx, dims(1), st, en, decomp%x1dist)
      call distribute(ny, dims(1), st, en, decomp%y1dist)
      deallocate (st, en)

      allocate (st(0:dims(2) - 1))
      allocate (en(0:dims(2) - 1))
      call distribute(ny, dims(2), st, en, decomp%y2dist)
      call distribute(nz, dims(2), st, en, decomp%z2dist)
      deallocate (st, en)

      return
   end subroutine get_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine prepare_buffer(decomp)

      implicit none

      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      integer :: i

      !LG : AJOUTS "bidons" pour eviter un plantage en -O3 avec gcc9.3
      !       * la fonction sortait des valeurs 'aleatoires'
      !         et le calcul plantait dans MPI_ALLTOALLV
      !       * pas de plantage en O2

      if (nrank == 0) then
         open (newunit=i, file='temp.dat', form='unformatted')
         write (i) decomp%x1dist, decomp%y1dist, decomp%y2dist, decomp%z2dist, &
            decomp%xsz, decomp%ysz, decomp%zsz
         close (i, status='delete')
      end if

      ! MPI_ALLTOALLV buffer information

      do i = 0, dims(1) - 1
         decomp%x1cnts(i) = decomp%x1dist(i) * decomp%xsz(2) * decomp%xsz(3)
         decomp%y1cnts(i) = decomp%ysz(1) * decomp%y1dist(i) * decomp%ysz(3)
         if (i == 0) then
            decomp%x1disp(i) = 0  ! displacement is 0-based index
            decomp%y1disp(i) = 0
         else
            decomp%x1disp(i) = decomp%x1disp(i - 1) + decomp%x1cnts(i - 1)
            decomp%y1disp(i) = decomp%y1disp(i - 1) + decomp%y1cnts(i - 1)
         end if
      end do

      do i = 0, dims(2) - 1
         decomp%y2cnts(i) = decomp%ysz(1) * decomp%y2dist(i) * decomp%ysz(3)
         decomp%z2cnts(i) = decomp%zsz(1) * decomp%zsz(2) * decomp%z2dist(i)
         if (i == 0) then
            decomp%y2disp(i) = 0  ! displacement is 0-based index
            decomp%z2disp(i) = 0
         else
            decomp%y2disp(i) = decomp%y2disp(i - 1) + decomp%y2cnts(i - 1)
            decomp%z2disp(i) = decomp%z2disp(i - 1) + decomp%z2cnts(i - 1)
         end if
      end do

      ! MPI_ALLTOALL buffer information
#ifdef EVEN
      ! For evenly distributed data, following is an easier implementation.
      ! But it should be covered by the more general formulation below.
      !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
      !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1)
      !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
      !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

      ! For unevenly distributed data, pad smaller messages. Note the
      ! last blocks along pencils always get assigned more mesh points
      ! for X <=> Y transposes
      decomp%x1count = decomp%x1dist(dims(1) - 1) * &
                       decomp%y1dist(dims(1) - 1) * decomp%xsz(3)
      decomp%y1count = decomp%x1count
      ! for Y <=> Z transposes
      decomp%y2count = decomp%y2dist(dims(2) - 1) * &
                       decomp%z2dist(dims(2) - 1) * decomp%zsz(1)
      decomp%z2count = decomp%y2count
#endif

   end subroutine prepare_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Postprocessing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "postpro.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"

end module decomp_2d
