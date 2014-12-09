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
# test_files="$test_files simple_check_c2r_c2c"
# test_files="$test_files simple_check_c2r_c2c_ousam_shifted"
# test_files="$test_files simple_check_c2r_c2c_shifted"
test_files="$test_files simple_check_ghost_c2c"
test_files_3d="$test_files_3d simple_check_ghost_c2c_3d_on_3d"
test_files="$test_files simple_check_ghost_r2c_input"
test_files="$test_files simple_check_ghost_r2c_input_padded"
test_files="$test_files simple_check_ghost_r2c_output"
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

test_files="$test_files simple_check_ousam_r2c_padded"

test_files="$test_files simple_check_c2c_newarray"
test_files="$test_files simple_check_c2c_transposed_newarray"
test_files="$test_files simple_check_r2c_newarray"
test_files="$test_files simple_check_r2c_transposed_newarray"
test_files="$test_files simple_check_r2c_padded_newarray"
test_files="$test_files simple_check_r2c_padded_transposed_newarray"
test_files="$test_files simple_check_ousam_c2c_4d_newarray"
test_files="$test_files simple_check_ousam_c2c_4d_transposed_newarray"
test_files="$test_files simple_check_ousam_r2c_4d_newarray"
test_files="$test_files simple_check_ousam_r2c_4d_transposed_newarray"

test_files="$test_files simple_check_c2c_inplace"
test_files="$test_files simple_check_c2c_transposed_inplace"


tol="1e-12"

failed="" 
## Run tests with 2d procmesh
for name in $test_files; do
  echo "## mpirun -np 4 ${name}:"
  res=$(mpirun -np 4 ./${name} | grep -i "error")
  echo $res
  err=$(echo $res | sed -n 's/.*maxerror = \([^;]*\);/\1/p')
  works=$(awk "BEGIN{print ($err < $tol)}")
  if [ $works -eq 0 ]; then
    failed="$failed $name"
  fi
done

## Run tests with 2d procmesh
failed_3d=""
for name in $test_files_3d; do
  echo "## mpirun -np 8 ${name}:"
  res=$(mpirun -np 8 ./${name} | grep -i "error")
  echo $res
  err=$(echo $res | sed -n 's/.*maxerror = \([^;]*\);/\1/p')
  works=$(awk "BEGIN{print ($err < $tol)}")
  if [ $works -eq 0 ]; then
    failed_3d="$failed_3d $name"
  fi
done

if test -z $failed && test -z $failed_3d; then
  echo -e  "\nAll tests reached accuracy level $tol."
else
  echo -e  "\nThe following tests did not reach accuracy level $tol:"
fi

for name in $failed; do
  echo "mpirun -np 4 $name"
done 
for name in $failed_3d; do
  echo "mpirun -np 8 $name"
done 
