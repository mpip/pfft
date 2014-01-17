
      program test

      implicit none

      include "mpif.h"
#include "fftw3.f"
      include "pfft.f"

      integer np(2), myrank, ierror, comm_cart_2d
      integer kinds_forw(3), kinds_back(3)
      integer(ptrdiff_t_kind) :: l, m
      integer(ptrdiff_t_kind) :: n(3), ni(3), no(3), ln(3)
      integer(ptrdiff_t_kind) :: alloc_local_forw, alloc_local_back, alloc_local
      integer(ptrdiff_t_kind) :: howmany
      integer(ptrdiff_t_kind) :: local_ni(3), local_i_start(3)
      integer(ptrdiff_t_kind) :: local_n(3), local_start(3)
      integer(ptrdiff_t_kind) :: local_no(3), local_o_start(3)
      integer(8) plan_forw, plan_back
      real(8), allocatable ::  data_in(:)
      real(8), allocatable ::  data_out(:)
      real(8) error

!     Set size of FFT and process mesh
      ni = (/ 16,16,16 /)
      n  = (/ 29,27,31 /)
      do l=1,3
        no(l) = ni(l)
      enddo
      howmany = 1

      np = (/ 2,2 /)

!     Set PFFT kinds of 1d R2R trafos
      kinds_forw(1) = PFFT_REDFT00; kinds_back(1) = PFFT_REDFT00;
      kinds_forw(2) = PFFT_REDFT01; kinds_back(2) = PFFT_REDFT10;
      kinds_forw(3) = PFFT_RODFT00; kinds_back(3) = PFFT_RODFT00;

!     Set logical DFT sizes corresponding to FFTW manual:
!     for REDFT00 N=2*(n-1), for RODFT00 N=2*(n+1), otherwise N=2*n
      ln(1) = 2*(n(1)-1);
      ln(2) = 2*n(2);
      ln(3) = 2*(n(3)+1); 

!     Initialize MPI and PFFT
      call MPI_Init(ierror)
      call dpfft_init();
      
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
      
!     Create two-dimensional process grid of
!     size np(1) x np(2), if possible
      call dpfft_create_procmesh_2d(ierror, MPI_COMM_WORLD, &
     &     np(1), np(2), comm_cart_2d)
      if (ierror .ne. 0) then
        if(myrank .eq. 0) then
          write(*,*) "Error: This test file only works with 4 processes"
        endif
        call MPI_Finalize(ierror)
        call exit(1)
      endif

!     Get parameters of data distribution
      call dpfft_local_size_many_r2r( &
     &     alloc_local_forw, 3, n, ni, n, howmany, &
     &     PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, &
     &     comm_cart_2d, PFFT_TRANSPOSED_NONE, &
     &     local_ni, local_i_start, local_n, local_start);

      call dpfft_local_size_many_r2r( &
     &     alloc_local_back, 3, n, n, no, howmany, &
     &     PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, &
     &     comm_cart_2d, PFFT_TRANSPOSED_NONE, &
     &     local_n, local_start, local_no, local_o_start);


!     Allocate enough memory for both trafos
      alloc_local = alloc_local_forw
      if(alloc_local_back .gt. alloc_local_forw) then
        alloc_local = alloc_local_back
      endif
      allocate(data_in(alloc_local))
      allocate(data_out(alloc_local))

!     Plan parallel forward FFT
      call dpfft_plan_many_r2r( &
     &     plan_forw, 3, n, ni, n, howmany, &
     &     PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, &
     &     data_in, data_out, comm_cart_2d, &
     &     kinds_forw, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_DESTROY_INPUT)
      
!     Plan parallel backward FFT
      call dpfft_plan_many_r2r( &
     &     plan_back, 3, n, n, no, howmany, &
     &     PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, &
     &     data_out, data_in, comm_cart_2d, &
     &     kinds_back, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_DESTROY_INPUT)

!     Initialize input with random numbers
      call dpfft_init_input_real_3d(ni, local_ni, local_i_start, &
     &     data_in)

!     Execute parallel forward FFT
      call dpfft_execute(plan_forw)

!     Execute parallel backward FFT
      call dpfft_execute(plan_back)

!     Scale data
      m=1
      do l=1,local_ni(1) * local_ni(2) * local_ni(3)
        data_in(l) = data_in(l) / (ln(1)*ln(2)*ln(3))
        m = m+1
      enddo

!     Print error of back transformed data
      call dpfft_check_output_real_3d(error, ni, local_ni, local_i_start, &
     &      data_in, comm_cart_2d)
      
      if(myrank .eq. 0) then
        write(*,*) "Error after one forward and backward&
            & trafo of size n=(", n(1), n(2), n(3), "):"
        write(*,*) "maxerror = ", error 
      endif

!     free mem and finalize
      call dpfft_destroy_plan(plan_forw)
      call dpfft_destroy_plan(plan_back)
      call MPI_Comm_free(comm_cart_2d, ierror)
      deallocate(data_out)
      deallocate(data_in)
  
!     Finalize MPI
      call MPI_Finalize(ierror)
      end

