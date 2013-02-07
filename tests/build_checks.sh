#!/bin/sh -e

FFTWDIR=$HOME/local/fftw-3.3-debug
SRCDIR=../../..

## Extra checks:
test_files="$test_files minimal_check_c2c minimal_check_c2c_transposed"

## 3D:
test_files="$test_files simple_check_c2c simple_check_c2c_transposed"
test_files="$test_files simple_check_r2c simple_check_r2c_transposed"
test_files="$test_files simple_check_r2r simple_check_r2r_transposed"

## 4D:
test_files="$test_files simple_check_c2c_4d simple_check_c2c_4d_transposed"
test_files="$test_files simple_check_r2c_4d simple_check_r2c_4d_transposed"
test_files="$test_files simple_check_r2r_4d simple_check_r2r_4d_transposed"

## 3D on 3D procmesh:
test_files="$test_files simple_check_c2c_3d_on_3d simple_check_c2c_3d_on_3d_transposed"
test_files="$test_files simple_check_r2c_3d_on_3d simple_check_r2c_3d_on_3d_transposed"
test_files="$test_files simple_check_r2r_3d_on_3d simple_check_r2r_3d_on_3d_transposed"

## 4D on 3D procmesh:
test_files="$test_files simple_check_c2c_4d_on_3d simple_check_c2c_4d_on_3d_transposed"
test_files="$test_files simple_check_r2c_4d_on_3d simple_check_r2c_4d_on_3d_transposed"
test_files="$test_files simple_check_r2r_4d_on_3d simple_check_r2r_4d_on_3d_transposed"



## 3D Oversampled:
test_files="$test_files simple_check_ousam_c2c simple_check_ousam_c2c_transposed"
test_files="$test_files simple_check_ousam_r2c simple_check_ousam_r2c_transposed"

## 4D Oversampled:
test_files="$test_files simple_check_ousam_r2c_4d simple_check_ousam_r2c_4d_transposed"
test_files="$test_files simple_check_ousam_c2c_4d simple_check_ousam_c2c_4d_transposed"

## 4D on 3D procmesh, Oversampled:
test_files="$test_files simple_check_ousam_c2c_4d_on_3d simple_check_ousam_c2c_4d_on_3d_transposed"
test_files="$test_files simple_check_ousam_r2c_4d_on_3d simple_check_ousam_r2c_4d_on_3d_transposed"

## Ghost cell send
test_files="$test_files simple_check_ghost_c2c"
test_files="$test_files simple_check_ghost_c2c_3d_on_3d"


PFFTLIB=../.libs/lib*pfft*.a

cd .. && make && cd - 
for name in $test_files; do
  mpicc $SRCDIR/tests/${name}.c -o ${name} \
    -I$FFTWDIR/include -I$SRCDIR/api -std=gnu99 -g -Wall \
    $PFFTLIB $FFTWDIR/lib/libfcs_fftw3.a
done
