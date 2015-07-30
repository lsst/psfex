
config(){
	cd lapack_functions
	scons -Q
	cd ../

	options=''
	options=$options"--prefix=$PREFIX "
	options=$options"--disable-threads "
	options=$options"--with-clapack=$PSFEX_DIR/lapack_functions "
	if [ "$FFTW_DIR" ]; then
		options=$options"--with-fftw-incdir=$FFTW_DIR/include "
		options=$options"--with-fftw-libdir=$FFTW_DIR/lib "
	fi
	#echo $options
	pwd
	./configure $options
}

build(){
	make CFLAGS=-g
	scons opt=3 prefix=$PREFIX version=$VERSION
}

install(){
	make install
	scons opt=3 install prefix=$PREFIX version=$VERSION
}

