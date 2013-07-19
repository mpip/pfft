#!/bin/sh -e

# change version numbers before new release
OLD_FFTW_VERSION=3.3.2
NEW_FFTW_VERSION=3.3.3

OLD_PFFT_VERSION=1.0.5-alpha
NEW_PFFT_VERSION=1.0.6-alpha

sed -i "s/\(^FFTW_VERSION=\)$OLD_FFTW_VERSION/\1$NEW_FFTW_VERSION/" *.sh
sed -i "s/\(^PFFT_VERSION=\)$OLD_PFFT_VERSION/\1$NEW_PFFT_VERSION/" *.sh

export BUILD_RELEASE_CALLED_FROM_MAKE_RELEASE="yes"
sh build_release_gcc.sh

