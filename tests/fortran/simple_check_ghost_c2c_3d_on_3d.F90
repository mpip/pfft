
      program test

      implicit none

      include "mpif.h"
#include "fftw3.f"
      include "pfft.f"
      
      integer np(3), myrank, ierror, comm_cart_3d
      integer(ptrdiff_t_kind) :: l, m
      integer(ptrdiff_t_kind) :: n(3)
      integer(ptrdiff_t_kind) :: alloc_local_gc, alloc_local
      integer(ptrdiff_t_kind) :: gc_below(3), gc_above(3)
      integer(ptrdiff_t_kind) :: local_ngc(3), local_gc_start(3)
      integer(ptrdiff_t_kind) :: local_ni(3), local_i_start(3)
      integer(ptrdiff_t_kind) :: local_no(3), local_o_start(3)
      integer(8) gcplan
      complex(8), allocatable ::  data_in(:)
      real(8) error
      
      n = (/ 8,8,8 /)
      np = (/ 2,2,2 /)
      gc_below = (/ 0,0,0 /)
      gc_above = (/ 0,0,8 /)

!     Initialize MPI and PFFT
      call MPI_Init(ierror)
      call dpfft_init();
      
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
      
!     Create two-dimensional process grid of
!     size np(1) x np(2), if possible
      call dpfft_create_procmesh(ierror, 3, MPI_COMM_WORLD, np,&
           & comm_cart_3d)
      if (ierror .ne. 0) then
        if(myrank .eq. 0) then
          write(*,*) "Error: This test file only works with 8 processes"
        endif
        call MPI_Finalize(ierror)
        call exit(1)
      endif
      
!     Get parameters of data distribution
      call dpfft_local_size_dft_3d(alloc_local, n, &
           & comm_cart_3d, PFFT_TRANSPOSED_NONE, &
           & local_ni, local_i_start, local_no, local_o_start)

!     Get parameters of ghostcell distribution
      call dpfft_local_size_gc_3d(alloc_local_gc, &
           & local_ni, local_i_start, alloc_local, &
           & gc_below, gc_above, &
           & local_ngc, local_gc_start)

!     Allocate memory
      allocate(data_in(alloc_local_gc))
!       data_in  = reshape(data_in,  (/ local_ni(1), local_ni(2), local_ni(3) /))

!     Plan parallel ghost cell send
      call dpfft_plan_cgc_3d(gcplan, n, gc_below, gc_above, &
           & data_in, comm_cart_3d, PFFT_GC_TRANSPOSED)
      
!     Initialize input with random numbers
      call dpfft_init_input_c2c_3d(n, local_ni, local_i_start, &
           & data_in)

!     Check gcell input
      call dpfft_apr_complex_3d(data_in, local_ni, local_i_start, &
           & "gcell input"//CHAR(0), comm_cart_3d)

!     Execute parallel ghost cell send
      call dpfft_exchange(gcplan)

!     Check gcell output 
      call dpfft_apr_complex_3d( &
           & data_in, local_ngc, local_gc_start, &
           & "exchanged gcells"//CHAR(0), comm_cart_3d)

!     Execute adjoint parallel ghost cell send
      call dpfft_reduce(gcplan)

!     Check input
      call dpfft_apr_complex_3d( &
           & data_in, local_no, local_o_start, &
           & "reduced gcells"//CHAR(0), comm_cart_3d)

!     Scale data
      m=1
      do l=1,local_ni(1) * local_ni(2) * local_ni(3)
        data_in(l) = data_in(l) / 3
        m = m+1
      enddo

!     Print error of back transformed data
      call dpfft_check_output_c2c_3d(error, n, local_ni, local_i_start, &
     &      data_in, comm_cart_3d)
      
      if(myrank .eq. 0) then
        write(*,*) "Error after one gcell exchange and reduce&
            & of size n=(", n(1), n(2), n(3), "),"
        write(*,*) "gc_below = (", gc_below(1), gc_below(2), gc_below(3), "),"
        write(*,*) "gc_above = (", gc_above(1), gc_above(2), gc_above(3), "):"
        write(*,*) "maxerror = ", error 
      endif

!     free mem and finalize
      call dpfft_destroy_gcplan(gcplan)
      call MPI_Comm_free(comm_cart_3d, ierror)
      deallocate(data_in)
  
!     Finalize MPI
      call dpfft_cleanup()
      call MPI_Finalize(ierror)
      end
      
