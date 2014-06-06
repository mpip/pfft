#!/bin/sh -e

FFTWDIR=$HOME/local/fftw-3.3-debug
SRCDIR=../../..

# test_files="$test_files adv_check_ghost_c2c"
# test_files="$test_files manual_min_c2c"
test_files="$test_files minimal_check_c2c"
test_files="$test_files minimal_check_c2c_transposed"
# test_files="$test_files serial_c2c"
test_files="$test_files simple_check_c2c"
test_files="$test_files simple_check_c2c_3d_on_3d"
test_files="$test_files simple_check_c2c_3d_on_3d_transposed"
test_files="$test_files simple_check_c2c_4d"
test_files="$test_files simple_check_c2c_4d_on_3d"
test_files="$test_files simple_check_c2c_4d_on_3d_transposed"
test_files="$test_files simple_check_c2c_4d_transposed"
test_files="$test_files simple_check_c2c_transposed"
test_files="$test_files simple_check_c2r_c2c"
test_files="$test_files simple_check_c2r_c2c_ousam_shifted"
test_files="$test_files simple_check_c2r_c2c_shifted"
test_files="$test_files simple_check_ghost_c2c"
test_files="$test_files simple_check_ghost_c2c_3d_on_3d"
test_files="$test_files simple_check_ghost_r2c_input"
test_files="$test_files simple_check_ghost_r2c_output"
test_files="$test_files simple_check_ousam_c2c"
test_files="$test_files simple_check_ousam_c2c_4d"
test_files="$test_files simple_check_ousam_c2c_4d_on_3d"
test_files="$test_files simple_check_ousam_c2c_4d_on_3d_transposed"
test_files="$test_files simple_check_ousam_c2c_4d_transposed"
test_files="$test_files simple_check_ousam_c2c_transposed"
test_files="$test_files simple_check_ousam_c2r"
test_files="$test_files simple_check_ousam_r2c"
test_files="$test_files simple_check_ousam_r2c_4d"
test_files="$test_files simple_check_ousam_r2c_4d_on_3d"
test_files="$test_files simple_check_ousam_r2c_4d_on_3d_transposed"
test_files="$test_files simple_check_ousam_r2c_4d_transposed"
test_files="$test_files simple_check_ousam_r2c_transposed"
test_files="$test_files simple_check_ousam_r2r"
test_files="$test_files simple_check_ousam_r2r_transposed"
test_files="$test_files simple_check_r2c"
test_files="$test_files simple_check_r2c_3d_on_3d"
test_files="$test_files simple_check_r2c_3d_on_3d_transposed"
test_files="$test_files simple_check_r2c_4d"
test_files="$test_files simple_check_r2c_4d_on_3d"
test_files="$test_files simple_check_r2c_4d_on_3d_transposed"
test_files="$test_files simple_check_r2c_4d_transposed"
test_files="$test_files simple_check_r2c_transposed"
test_files="$test_files simple_check_r2r"
test_files="$test_files simple_check_r2r_3d_on_3d"
test_files="$test_files simple_check_r2r_3d_on_3d_transposed"
test_files="$test_files simple_check_r2r_4d"
test_files="$test_files simple_check_r2r_4d_on_3d"
test_files="$test_files simple_check_r2r_4d_on_3d_transposed"
test_files="$test_files simple_check_r2r_4d_transposed"
test_files="$test_files simple_check_r2r_transposed"

test_files="$test_files simple_check_ousam_r2c_padded"

PFFTLIB=../.libs/lib*pfft*.a

cd .. && make && cd - 
for name in $test_files; do
  mpicc $SRCDIR/tests/${name}.c -o ${name} \
    -I$FFTWDIR/include -I$SRCDIR/api -std=gnu99 -g -Wall \
    $PFFTLIB $FFTWDIR/lib/libfcs_fftw3.a
done
