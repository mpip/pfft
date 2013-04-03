! This program tests the over-/under-sampled PFFT.
! We first compute an oversampled forward PFFT and second
! revert it with an undersampled backward PFFT.
! This should give back the initial input values
! (up to the typical scaling factor)

program main
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  include "fftw3-mpi.f03"
  include "pfft.f03"

  integer np(3)
  integer(C_INTPTR_T) :: n(4), ni(4), no(4)
  integer(C_INTPTR_T) :: alloc_local_forw, alloc_local_back, alloc_local
  integer(C_INTPTR_T) :: howmany
  integer(C_INTPTR_T) :: local_ni(4), local_i_start(4)
  integer(C_INTPTR_T) :: local_n(4), local_start(4)
  integer(C_INTPTR_T) :: local_no(4), local_o_start(4)
  integer(C_INTPTR_T) :: blocks(4)
  complex(C_DOUBLE_COMPLEX), pointer :: in(:,:,:,:), out(:,:,:,:)
  type(C_PTR) :: plan_forw, plan_back, cin, cout
  integer comm_cart_3d
  integer myrank, ierror
  double precision err

  ! Set size of FFT and process mesh
  ni = [ 8, 8, 8, 8]
  n  = [13,14,19,17]
  no = [ 8, 8, 8, 8]
  np = [2,2,2]
  howmany = 1
  blocks = [PFFT_DEFAULT_BLOCK,PFFT_DEFAULT_BLOCK,PFFT_DEFAULT_BLOCK,PFFT_DEFAULT_BLOCK]

  ! Initialize MPI and PFFT
  call MPI_Init(ierror)
  call pfft_init()
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)

  ! Create three-dimensional process grid of size np(1) x np(2) x np(3), if possible
  ierror =  pfft_create_procmesh(3, MPI_COMM_WORLD, np, comm_cart_3d)
  if (ierror .ne. 0) then
    if(myrank .eq. 0) then
      write(*,*) "Error: This test file only works with ", np(1)*np(2)*np(3), " processes"
    endif
    call MPI_Finalize(ierror)
    call exit(1)
  endif

  ! Get parameters of data distribution
  alloc_local_forw = pfft_local_size_many_dft(4, n, ni, n, howmany, &
      blocks, blocks, &
      comm_cart_3d, PFFT_TRANSPOSED_NONE, &
      local_ni, local_i_start, local_n, local_start)

  alloc_local_back = pfft_local_size_many_dft(4, n, n, no, howmany, &
      blocks, blocks, &
      comm_cart_3d, PFFT_TRANSPOSED_NONE, &
      local_n, local_start, local_no, local_o_start)

  ! Allocate enough memory for both trafos
  alloc_local = max(alloc_local_forw,alloc_local_back)
  cin  = pfft_alloc_complex(alloc_local)
  cout = pfft_alloc_complex(alloc_local)

  ! Convert data pointers to Fortran format
  call c_f_pointer(cin,  in,  [local_ni(4),local_ni(3), local_ni(2), local_ni(1)])
  call c_f_pointer(cout, out, [local_no(4),local_no(3), local_no(2), local_no(1)])

  ! Plan parallel forward FFT
  plan_forw = pfft_plan_many_dft( &
      4, n, ni, n, howmany, blocks, blocks, &
      in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_DESTROY_INPUT)

  ! Plan parallel backward FFT
  plan_back = pfft_plan_many_dft( &
      4, n, n, no, howmany, blocks, blocks, &
      out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_DESTROY_INPUT)

  ! Initialize input with random numbers
  call pfft_init_input_c2c(4, ni, local_ni, local_i_start, &
      in)

  ! Execute parallel forward FFT
  call pfft_execute(plan_forw)
  
  ! Execute parallel backward FFT
  call pfft_execute(plan_back)

  ! Scale data
  in = in / (n(1)*n(2)*n(3)*n(4))

  ! Print error of back transformed data
  err = pfft_check_output_c2c(4, ni, local_ni, local_i_start, in, comm_cart_3d)
  if(myrank .eq. 0) then
    write(*,*) "Error after one forward and backward trafo of size n=(", n(1), ", ", n(2), ", ", n(3), ", ", n(4), "):"
    write(*,*) "maxerror = ", err
  endif 

  ! Free mem and finalize
  call pfft_destroy_plan(plan_forw)
  call pfft_destroy_plan(plan_back)
  call MPI_Comm_free(comm_cart_3d, ierror)
  call pfft_free(cin)
  call pfft_free(cout)
  call MPI_Finalize(ierror)
end program main

