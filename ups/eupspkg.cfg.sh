
config()
{
	function finish {
		cat config.log
  	}
 	trap finish EXIT
	options="--prefix=$PREFIX --disable-threads --with-lapack-stub "
	options=$options"--enable-plplot=no "
	# The following if statement allows the package to be built if
	# someone chooses to use their system libraries for fftw in place of
	# the package provided through eups. It checks if the eups fftw
	# environment variable is set, and if it is the configuration is
	# set appropriately.
	if [ "$FFTW_DIR" ]; then
		options=$options"--with-fftw-incdir=$FFTW_DIR/include "
		options=$options"--with-fftw-libdir=$FFTW_DIR/lib "
	elif [ "$CONDA_PREFIX" ]; then
                options=$options"--enable-openblas "
		options=$options"--with-fftw-incdir=${CONDA_PREFIX}/include "
		options=$options"--with-fftw-libdir=${CONDA_PREFIX}/lib "
		options=$options"--with-openblas-incdir=${CONDA_PREFIX}/include "
		options=$options"--with-openblas-libdir=${CONDA_PREFIX}/lib "
        fi
	./configure ${options}
}

install()
{
    default_install
    # Install headers under 'include' directory
    # in case the user wants to delete 'src'.
    mkdir -p $PREFIX/include
    cp config.h $PREFIX/include/
    cp src/*.h $PREFIX/include/
    mkdir -p $PREFIX/include/wcs
    cp src/wcs/*.h $PREFIX/include/wcs
    mkdir -p $PREFIX/include/levmar/
    cp src/levmar/*.h $PREFIX/include/levmar
    mkdir -p $PREFIX/include/fits
    cp src/fits/*.h $PREFIX/include/fits
}

