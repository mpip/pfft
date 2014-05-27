#!/bin/sh -e

test_files=""

# test_files="$test_files adv_check_ghost_c2c"
# test_files="$test_files manual_min_c2c"
test_files="$test_files minimal_check_c2c"
test_files="$test_files minimal_check_c2c_transposed"
# test_files="$test_files serial_c2c"
test_files="$test_files simple_check_c2c"
test_files_3d="$test_files_3d simple_check_c2c_3d_on_3d"
test_files_3d="$test_files_3d simple_check_c2c_3d_on_3d_transposed"
test_files="$test_files simple_check_c2c_4d"
test_files_3d="$test_files_3d simple_check_c2c_4d_on_3d"
test_files_3d="$test_files_3d simple_check_c2c_4d_on_3d_transposed"
test_files="$test_files simple_check_c2c_4d_transposed"
test_files="$test_files simple_check_c2c_transposed"
test_files="$test_files simple_check_c2r_c2c"
test_files="$test_files simple_check_c2r_c2c_ousam_shifted"
test_files="$test_files simple_check_c2r_c2c_shifted"
test_files="$test_files simple_check_ghost_c2c"
test_files_3d="$test_files_3d simple_check_ghost_c2c_3d_on_3d"
test_files="$test_files simple_check_ghost_c2r"
test_files="$test_files simple_check_ousam_c2c"
test_files="$test_files simple_check_ousam_c2c_4d"
test_files_3d="$test_files_3d simple_check_ousam_c2c_4d_on_3d"
test_files_3d="$test_files_3d simple_check_ousam_c2c_4d_on_3d_transposed"
test_files="$test_files simple_check_ousam_c2c_4d_transposed"
test_files="$test_files simple_check_ousam_c2c_transposed"
test_files="$test_files simple_check_ousam_c2r"
test_files="$test_files simple_check_ousam_r2c"
test_files="$test_files simple_check_ousam_r2c_4d"
test_files_3d="$test_files_3d simple_check_ousam_r2c_4d_on_3d"
test_files_3d="$test_files_3d simple_check_ousam_r2c_4d_on_3d_transposed"
test_files="$test_files simple_check_ousam_r2c_4d_transposed"
test_files="$test_files simple_check_ousam_r2c_transposed"
test_files="$test_files simple_check_ousam_r2r"
test_files="$test_files simple_check_ousam_r2r_transposed"
test_files="$test_files simple_check_r2c"
test_files_3d="$test_files_3d simple_check_r2c_3d_on_3d"
test_files_3d="$test_files_3d simple_check_r2c_3d_on_3d_transposed"
test_files="$test_files simple_check_r2c_4d"
test_files_3d="$test_files_3d simple_check_r2c_4d_on_3d"
test_files_3d="$test_files_3d simple_check_r2c_4d_on_3d_transposed"
test_files="$test_files simple_check_r2c_4d_transposed"
test_files="$test_files simple_check_r2c_transposed"
test_files="$test_files simple_check_r2r"
test_files_3d="$test_files_3d simple_check_r2r_3d_on_3d"
test_files_3d="$test_files_3d simple_check_r2r_3d_on_3d_transposed"
test_files="$test_files simple_check_r2r_4d"
test_files_3d="$test_files_3d simple_check_r2r_4d_on_3d"
test_files_3d="$test_files_3d simple_check_r2r_4d_on_3d_transposed"
test_files="$test_files simple_check_r2r_4d_transposed"
test_files="$test_files simple_check_r2r_transposed"

 
## Run tests with 2d procmesh
for name in $test_files; do
  echo "## mpirun -np 4 ${name}:"
  mpirun -np 4 ./${name}
done

## Run tests with 2d procmesh
for name in $test_files_3d; do
  echo "## mpirun -np 8 ${name}:"
  mpirun -np 8 ./${name}
done
