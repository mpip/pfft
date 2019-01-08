! a simple example adapted from simple_check_ousam_r2c_4d_on_3d.F03

program main
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  include "fftw3-mpi.f03"
  include "pfft.f03"

  integer np(2)
  integer(C_INTPTR_T) :: n(3), nr(3), nc(3)
  integer(C_INTPTR_T) :: alloc_local
  integer(C_INTPTR_T) :: howmany
  integer(C_INTPTR_T) :: local_nr(3), local_r_start(3)
  integer(C_INTPTR_T) :: local_nc(3), local_c_start(3)
  integer(C_INTPTR_T) :: blocks(3)
  real(C_DOUBLE), pointer :: dr(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: dc(:,:,:)
  type(C_PTR) :: plan_c2r, plan_r2c, cdc, cdr
  integer comm_cart
  integer myrank, ierror
  double precision err

  ! Pruned size (complex arrays)
  nc = [ 3, 3, 9]

  ! padded size (transform and the real array)
  n  = nc + 5
  nr = n

  ! processor grid
  np = [ 2, 1]
  howmany = 1
  blocks = [PFFT_DEFAULT_BLOCK,PFFT_DEFAULT_BLOCK,PFFT_DEFAULT_BLOCK]

  ! Initialize MPI and PFFT
  call MPI_Init(ierror)
  call pfft_init()
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)

  ! Create process grid 
  ierror =  pfft_create_procmesh(2, MPI_COMM_WORLD, np, comm_cart)
  if (ierror .ne. 0) then
    if(myrank .eq. 0) then
      write(*,*) "Error: This test file only works with ", np(1)*np(2), " processes"
    endif
    call MPI_Finalize(ierror)
    call exit(1)
  endif

  ! Get parameters of data distribution
  alloc_local = pfft_local_size_many_dft_r2c(3, n, nr, nc, howmany, &
      blocks, blocks, &
      comm_cart, PFFT_TRANSPOSED_NONE, &
      local_nr, local_r_start, local_nc, local_c_start)
  
  ! Allocate memory

  cdr  = pfft_alloc_real( local_nr(1)*local_nr(2)*local_nr(3) )
  cdc  = pfft_alloc_complex(alloc_local)

  ! Convert data pointers to Fortran format
  call c_f_pointer(cdr, dr, [local_nr(3), local_nr(2), local_nr(1)])
  call c_f_pointer(cdc, dc, [local_nc(3), local_nc(2), local_nc(1)])

  ! Plan parallel R2C FFT
  plan_r2c = pfft_plan_many_dft_r2c( &
       3, n, nr, nc, howmany, blocks, blocks, &
       dr, dc, comm_cart, PFFT_FORWARD, PFFT_TRANSPOSED_NONE + PFFT_MEASURE )

  ! Plan parallel C2R FFT
  plan_c2r = pfft_plan_many_dft_c2r( &
       3, n, nc, nr, howmany, blocks, blocks, &
       dc, dr, comm_cart, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE + PFFT_MEASURE )

  ! Initialize complex input with random numbers
  call pfft_init_input_complex_hermitian( 3, nc, local_nc, local_c_start, dc)

  ! Execute parallel backward FFT
  call pfft_execute(plan_c2r)
  
  ! Execute parallel forward FFT
  call pfft_execute(plan_r2c)
  
  ! Scale data
  dc = dc / (n(1)*n(2)*n(3))

    ! Print error of back transformed data
  err = pfft_check_output_complex_hermitian(3, nc, local_nc, local_c_start, dc, comm_cart)
  if(myrank .eq. 0) then
    write(*,*) "Error after one forward and backward trafo of size n=(", n(1), ", ", n(2), ", ", n(3), ", ",  "):"
    write(*,*) "maxerror = ", err
  endif 

  ! Free mem and finalize
  call pfft_destroy_plan(plan_r2c)
  call pfft_destroy_plan(plan_c2r)
  call MPI_Comm_free(comm_cart, ierror)
  call pfft_free(cdr)
  call pfft_free(cdc)
  call MPI_Finalize(ierror)
end program main

