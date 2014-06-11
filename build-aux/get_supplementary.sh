# Assume we want to install them below $HOME/local.
myprefix=${HOME}/local
tempdir="./supplementary.tmp"

# Show where the auto-tools will be installed to
echo "supplementary libraries will be installed to: $myprefix"

# Create a temporary directory to store the tar-balls
echo "creating temporary directory"
mkdir ${tempdir}
cd ${tempdir}

GSL=gsl-1.9

M4=m4-1.4.16
AUTOCONF=autoconf-2.68
AUTOMAKE=automake-1.12.3
LIBTOOL=libtool-2.4.2

# Do the following in a scratch directory.
echo "downloading supplementary libraries"
wget ftp://ftp.gnu.org/gnu/gsl/${GSL}.tar.gz
echo "extracting supplementary libraries"
gzip -dc ${GSL}.tar.gz | tar xvf -
echo "installing supplementary libraries to $myprefix"
cd ${GSL}
./configure -C --prefix=$myprefix && make && make install
export PATH=${myprefix}:${PATH}

# Clean up the temporary directory and downloaded tar-balls
cd ../..
#echo "deleting temporary directory `pwd`/${tempdir}"
#rm -rf ${tempdir}

# Notification to update the PATH variable
echo "Please add $myprefix/lib to your LD_LIBRARY_PATH variable"
echo "and configure ScaFaCoS using "
echo "      CFLAGS=-I$myprefix/include ../scafacos.git/configure ..."
