#if !defined(WRITE_FITS_IMAGE)
#define WRITE_FITS_IMAGE 1
void write_fits_image(const char *filename, float *data, int nx, int ny, float xc, float yc);
#endif
