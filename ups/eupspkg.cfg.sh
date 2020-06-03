
config()
{
	function finish {
		cat config.log
  	}
 	trap finish EXIT
	options="--prefix=$PREFIX --disable-threads --with-lapack-stub "
	options=$options"--enable-plplot=no "
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

