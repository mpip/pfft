program main
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  include "fftw3-mpi.f03"
  include "pfft.f03"

  integer np(2)
  integer(C_INTPTR_T) :: n(3)
  integer(C_INTPTR_T) :: alloc_local
  integer(C_INTPTR_T) :: local_ni(3), local_i_start(3)
  integer(C_INTPTR_T) :: local_no(3), local_o_start(3)
  double precision err
  complex(C_DOUBLE_COMPLEX), pointer :: in(:,:,:), out(:,:,:)
  type(C_PTR) :: plan_forw, plan_back, cin, cout
  integer comm_cart_2d

  integer myrank, ierror

  n  = [31,27,29]
  np = [2,2]

  ! Initialize MPI and PFFT
  call MPI_Init(ierror)
  call pfft_init()
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)

  ! Create two-dimensional process grid of size np(1) x np(2), if possible
  ierror =  pfft_create_procmesh_2d(MPI_COMM_WORLD, np(1), np(2), comm_cart_2d)
  if (ierror .ne. 0) then
    if(myrank .eq. 0) then
      write(*,*) "Error: This test file only works with ", np(1)*np(2), " processes"
    endif
    call MPI_Finalize(ierror)
    call exit(1)
  endif

  ! Get parameters of data distribution
  alloc_local = pfft_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE, &
      local_ni, local_i_start, local_no, local_o_start);

  ! Allocate memory
  cin  = pfft_alloc_complex(alloc_local)
  cout = pfft_alloc_complex(alloc_local)

  ! Convert data pointers to Fortran format
  call c_f_pointer(cin,  in,  [local_ni(3), local_ni(2), local_ni(1)])
  call c_f_pointer(cout, out, [local_no(3), local_no(2), local_no(1)])

  ! Plan parallel forward FFT
  plan_forw = pfft_plan_dft_3d( &
      n, in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_DESTROY_INPUT)

  ! Plan parallel backward FFT
  plan_back = pfft_plan_dft_3d( &
      n, out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_DESTROY_INPUT)

  ! Initialize input with random numbers
  call pfft_init_input_c2c_3d(n, local_ni, local_i_start, &
      in)

  ! Execute parallel forward FFT
  call pfft_execute(plan_forw)
  
  ! Execute parallel backward FFT
  call pfft_execute(plan_back)

  ! Scale data
  in = in / (n(1)*n(2)*n(3))
  
  ! Print error of back transformed data
  err = pfft_check_output_c2c_3d(n, local_ni, local_i_start, in, comm_cart_2d)
  if(myrank .eq. 0) then
    write(*,*) "Error after one forward and backward trafo of size n=(", n(1), ", ", n(2), ", ", n(3), "):"
    write(*,*) "maxerror = ", err
  endif 

  ! Free mem and finalize
  call pfft_destroy_plan(plan_forw)
  call pfft_destroy_plan(plan_back)
  call MPI_Comm_free(comm_cart_2d, ierror)
  call pfft_free(cin)
  call pfft_free(cout)
  call MPI_Finalize(ierror)
end program main

