
config()
{
	function finish {
		cat config.log
  	}
  	trap finish EXIT
	cd lapack_functions
	# Run the lapack_functions build script. This is a wrapper
	# around scons which checks and removes any sconsflags which
	# may be set but which break the build
	python build.py
	cd ../
	cp lapack_functions/lib/*lapack* lib/
	if [ -z $PSFEX_DIR ]; then
		PSFEX_DIR=$(pwd)
	fi
	options="--prefix=$PREFIX --disable-threads --with-clapack=$PSFEX_DIR/lapack_functions "
	options=$options"--enable-plplot=no "
	# The following if statement allows the package to be built if
	# someone chooses to use their system libraries for fftw in place of
	# the package provided through eups. It checks if the eups fftw
	# environment variable is set, and if it is the configuration is
	# set appropriately.
	if [ "$FFTW_DIR" ]; then
		options=$options"--with-fftw-incdir=$FFTW_DIR/include "
		options=$options"--with-fftw-libdir=$FFTW_DIR/lib "
	fi
	./configure $options
}

build()
{
	make 
	scons prefix=$PREFIX version=$VERSION
}

install()
{
	clean_old_install
	make install
	scons install prefix=$PREFIX version=$VERSION
}

