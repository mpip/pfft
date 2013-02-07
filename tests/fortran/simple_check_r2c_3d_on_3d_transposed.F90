

      program test

      implicit none

      include "mpif.h"
#include "fftw3.f"
      include "pfft.f"

      integer np(3), myrank, ierror, comm_cart_3d
      integer(ptrdiff_t_kind) :: l, m
      integer(ptrdiff_t_kind) :: n(3)
      integer(ptrdiff_t_kind) :: alloc_local
      integer(ptrdiff_t_kind) :: local_ni(3), local_i_start(3)
      integer(ptrdiff_t_kind) :: local_no(3), local_o_start(3)
      integer(8) plan_forw, plan_back
      real(8), allocatable ::  data_in(:)
      complex(8), allocatable ::  data_out(:)
      real(8) error

      n = (/ 29,27,31 /)
      np = (/ 2,2,2 /)

!     Initialize MPI and PFFT
      call MPI_Init(ierror)
      call dpfft_init();
      
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
      
!     Create two-dimensional process grid of
!     size np(1) x np(2), if possible
      call dpfft_create_procmesh(ierror, 3, MPI_COMM_WORLD, &
     &     np, comm_cart_3d)
      if (ierror .ne. 0) then
        if(myrank .eq. 0) then
          write(*,*) "Error: This test file only works with 8 processes"
        endif
        call MPI_Finalize(ierror)
        call exit(1)
      endif

!     Get parameters of data distribution
      call dpfft_local_size_dft_r2c_3d( &
     &     alloc_local, n, comm_cart_3d, PFFT_TRANSPOSED_OUT, &
     &     local_ni, local_i_start, local_no, local_o_start);

!     Allocate memory
      allocate(data_in(2*alloc_local))
      allocate(data_out(alloc_local))

!     Plan parallel forward FFT
      call dpfft_plan_dft_r2c_3d(plan_forw, n, data_in, data_out, comm_cart_3d, &
     &     PFFT_FORWARD, PFFT_TRANSPOSED_OUT + PFFT_MEASURE + PFFT_DESTROY_INPUT)
      
!     Plan parallel backward FFT
      call dpfft_plan_dft_c2r_3d(plan_back, n, data_out, data_in, comm_cart_3d, &
     &     PFFT_BACKWARD, PFFT_TRANSPOSED_IN + PFFT_MEASURE + PFFT_DESTROY_INPUT)

!     Initialize input with random numbers
      call dpfft_init_input_r2c_3d(n, local_ni, local_i_start, &
     &     data_in)

!     Initialize input with random numbers
!      m = 0
!      do k2=0,n(2)
!        do k1=0,n(1)
!          do k0=0,n(0)
!            in(m) = 1000.0/(m+1)
!            m = m+1
!          enddo
!          do k0=n(0), n(0)/2+1
!            in(m) = 0
!            m= m+1
!          enddo
!        enddo
!      enddo


!     Execute parallel forward FFT
      call dpfft_execute(plan_forw)

!     Execute parallel backward FFT
      call dpfft_execute(plan_back)


!     Scale data
      m=1
      do l=1, local_ni(1) * local_ni(2) * local_ni(3) 
        data_in(l) = data_in(l) / (n(1)*n(2)*n(3))
        m = m+1
      enddo

!     Print error of back transformed data
      call dpfft_check_output_c2r_3d(error, n, local_ni, local_i_start, &
     &      data_in, comm_cart_3d)
      
      if(myrank .eq. 0) then
        write(*,*) "Error after one forward and backward&
            & trafo of size n=(", n(1), n(2), n(3), "):"
        write(*,*) "maxerror = ", error 
      endif

!     free mem and finalize
      call dpfft_destroy_plan(plan_forw)
      call dpfft_destroy_plan(plan_back)
      call MPI_Comm_free(comm_cart_3d, ierror)
      deallocate(data_out)
      deallocate(data_in)
  
!     Finalize MPI
      call MPI_Finalize(ierror)
      end













  
  
