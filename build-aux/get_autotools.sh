# Assume we want to install them below $HOME/local.
myprefix=${HOME}/local
tempdir="./autotools.tmp"

# Show where the auto-tools will be installed to
echo "auto-tools will be installed to: $myprefix"

# Create a temporary directory to store the tar-balls
echo "creating temporary directory"
mkdir ${tempdir}
cd ${tempdir}

M4=m4-1.4.16
AUTOCONF=autoconf-2.68
AUTOMAKE=automake-1.12.3
LIBTOOL=libtool-2.4.2

# Do the following in a scratch directory.
echo "downloading auto-tools"
wget http://ftp.gnu.org/gnu/m4/${M4}.tar.gz
wget http://ftp.gnu.org/gnu/autoconf/${AUTOCONF}.tar.gz
wget http://ftp.gnu.org/gnu/automake/${AUTOMAKE}.tar.gz
wget http://ftp.gnu.org/gnu/libtool/${LIBTOOL}.tar.gz
echo "extracting auto-tools"
gzip -dc ${M4}.tar.gz | tar xvf -
gzip -dc ${AUTOCONF}.tar.gz | tar xvf -
gzip -dc ${AUTOMAKE}.tar.gz | tar xvf -
gzip -dc ${LIBTOOL}.tar.gz | tar xvf -
echo "installing auto-tools to $myprefix"
cd ${M4}
./configure -C --prefix=$myprefix && make && make install
export PATH=${myprefix}:${PATH}
cd ../${AUTOCONF}
./configure -C --prefix=$myprefix && make && make install
cd ../${AUTOMAKE}
./configure -C --prefix=$myprefix && make && make install
cd ../${LIBTOOL}
./configure -C --prefix=$myprefix && make && make install

# Clean up the temporary directory and downloaded tar-balls
cd ../..
#echo "deleting temporary directory `pwd`/${tempdir}"
#rm -rf ${tempdir}

# Notification to update the PATH variable
echo "$myprefix/bin was added to your PATH variable to use the newly installed auto-tools"
