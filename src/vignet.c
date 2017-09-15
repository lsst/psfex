/*
*				vignet.c
*
* Manipulate image rasters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<assert.h>
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"vignet.h"


/****** vignet_resample ******************************************************
PROTO	int	vignet_resample(float *pix1, int w1, int h1,
		float *pix2, int w2, int h2, double dx, double dy, float step2,
		float stepi)
PURPOSE	Scale and shift a small image through sinc interpolation, with
	adjustable spatial wavelength cut-off. Image parts which lie outside
	boundaries are set to 0.

INPUT	Input raster,
	input raster width,
	input raster height,
	output raster,
	output raster width,
	output raster height,
	shift in x,
	shift in y,
	output pixel scale.	
OUTPUT	RETURN_ERROR if the images do not overlap, RETURN_OK otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	13/09/2010
 ***/

int
vignet_resample(const float *s_pix1, const int w1, const int h1, /* input */
		float *s_pix2, const int w2, const int h2,       /* output */
		const double dx, const double dy,
		const float step2,
		float stepi)
{
   static float	*stat_s_pix2 = NULL;	/* save an old version of s_pix2.  Scary! */

   if (stepi <= 0.0) {
      stepi = 1.0;
   }
   const double dstepi = 1.0/stepi;
   const double mx1 = w1/2;		/* Im1 center x-coord*/
   const double mx2 = w2/2;		/* Im2 center x-coord*/
   double xs1 = mx1 + dx - mx2*step2;	/* Im1 start x-coord */

   if ((int)xs1 >= w1) {
      return RETURN_ERROR;
   }
   int ixs2 = 0;			/* Int part of Im2 start x-coord */
   if (xs1 < 0.0) {
      const int dix2 = 1 - xs1/step2;
      /*-- Simply leave here if the images do not overlap in x */
      if (dix2 >= w2) {
	 return RETURN_ERROR;
      }
      ixs2 += dix2;
      xs1 += dix2*step2;
   }
   int nx2 = (w1 - 1 - xs1)/step2 + 1;	/* nb of interpolated Im2 pixels along x */
   const int ix2 = w2 - ixs2;
   if (nx2 > ix2) {
      nx2 = ix2;
   }
   if (nx2 <= 0) {
      return RETURN_ERROR;
   }
   const double my1 = h1/2;		/* Im1 center y-coord */
   const double my2 = h2/2;		/* Im2 center y-coord */
   double ys1 = my1 + dy - my2*step2;	/* Im1 start y-coord */
   if ((int)ys1 >= h1) {
      return RETURN_ERROR;
   }
   int iys2 = 0;			/* Int part of Im2 start y-coord */
   if (ys1 < 0.0) {
      const int diy2 = 1 - ys1/step2;
      /*-- Simply leave here if the images do not overlap in y */
      if (diy2 >= h2) {
	 return RETURN_ERROR;
      }
      iys2 += diy2;
      ys1 += diy2*step2;
   }
   int ny2 = (h1 - 1 - ys1)/step2 + 1;	/* nb of interpolated Im2 pixels along y */
   const int iy2 = h2 - iys2;
   if (ny2 > iy2) {
      ny2 = iy2;
   }
   if (ny2 <= 0) {
      return RETURN_ERROR;
   }
   
   /* Set the yrange for the x-resampling with some margin for interpolation */
   int iys1a = ys1;			   /* Int part of Im1 start y-coord with margin */
   const int hmh = (INTERPW/2)/dstepi + 2; /* Interpolant start */
   const int interph = 2*hmh;
   const int hmw = (INTERPW/2)/dstepi + 2;
   const int interpw =  2*hmw;
   iys1a -= hmh;
   if (iys1a < 0) {
      iys1a = 0;
   }
   int ny1 = (ys1 + ny2*step2) + interpw - hmh;	/* Interpolated Im1 y size */
   if (ny1 > h1) {				/* with margin */
      ny1 = h1;
   }
   /* Express everything relative to the effective Im1 start (with margin) */
   ny1 -= iys1a;
   ys1 -= (double)iys1a;

   /* Initialize destination buffer to zero if s_pix2 != NULL */
   if (!s_pix2) {
      assert(stat_s_pix2 != NULL);
      s_pix2 = stat_s_pix2;
   } else {
      memset(s_pix2, 0, (size_t)(w2*h2)*sizeof(float));
      stat_s_pix2 = s_pix2;
   }

   /* Allocate 2-D arrays for data */
   const float **pix1; QMALLOC(pix1, const float *, h1);
   for (int k = 0; k < h1; k++) {
      pix1[k] = &s_pix1[k*w1];
   }

   float **pix2; QMALLOC(pix2, float *, h2);
   for (int k = 0; k < h2; k++) {
      pix2[k] = &s_pix2[k*w2];
   }

   /* Allocate interpolant stuff for the x direction */
   double *s_mask; QMALLOC(s_mask, double, nx2*interpw); /* Storage for interpolation masks */
   double **mask; QMALLOC(mask, double *, nx2); /* interpolation masks */
   for (int k = 0; k < nx2; k++) {
      mask[k] = &s_mask[k*interpw];
   }

   int *nmask;   QMALLOC(nmask, int, nx2);	     /* Interpolation mask sizes */
   int *start;   QMALLOC(start, int, nx2);	     /* Int part of Im1 conv starts */

   /* Compute the local interpolant and data starting points in x */
   double x1 = xs1;
   for (int j = 0; j < nx2; j++, x1 += step2) {
      const int ix1 = x1;
      int ix = ix1 - hmw;
      double dxm = (ix1 - x1 - hmw)*dstepi;/* starting point in the interp. func */
      int n;
      if (ix < 0) {
	 n = interpw + ix;
	 dxm -= (double)ix*dstepi;
	 ix = 0;
      } else {
	 n = interpw;
      }
      int t = w1 - ix;
      if (n > t) {
	 n = t;
      }
      start[j] = ix;
      nmask[j] = n;
      double norm = 0.0;
      double x = dxm;
      for (int i = 0; i < n; i++, x += dstepi) {
	 double pval = INTERPF(x);
	 mask[j][i] = pval;
	 norm += pval;
      }
      norm = (norm > 0.0) ? 1.0/norm : dstepi;
      for (int i = 0; i < n; i++) {
	 mask[j][i] *= norm;
      }
   }

   float *s_pix_tmp; QCALLOC(s_pix_tmp, float, nx2*ny1); /* Intermediary frame-buffer */
   float **pix_tmp; QCALLOC(pix_tmp, float *, ny1);
   for (int k = 0; k < ny1; k++) {
      pix_tmp[k] = &s_pix_tmp[k*nx2];
   }

   /* Make the interpolation in x (this includes transposition) */
   for (int i = 0; i < ny1; i++) {
      for (int j = 0; j < nx2; j++) {
	 float val = 0.0; 
	 for (int k = 0; k < nmask[j]; k++) {
	    const double pval = pix1[iys1a + i][start[j] + k];
	    val += pval*mask[j][k];
	 }
	 pix_tmp[i][j] = val;
      }
   }
  
   /* Reallocate interpolant stuff for the y direction */
   QREALLOC(s_mask, double, ny2*interph);	/* Interpolation masks */
   QREALLOC(mask, double *, ny2*interph);	/* Interpolation masks */
   for (int k = 0; k < ny2; k++) {
      mask[k] = &s_mask[k*interph];
   }
   
   QREALLOC(nmask, int, ny2);		/* Interpolation mask sizes */
   QREALLOC(start, int, ny2);		/* Int part of Im1 conv starts */
   
   /* Compute the local interpolant and data starting points in y */
   double y1 = ys1;
   for (int j = 0; j < ny2; j++, y1 += step2) {
      const int iy1 = y1;
      int iy = iy1 - hmh;
      double dym = (iy1 - y1 - hmh)*dstepi;/* starting point in the interp. func */
      int n;
      if (iy < 0) {
	 n = interph + iy;
	 dym -= (double)iy*dstepi;
	 iy = 0;
      } else {
	 n = interph;
      }
      const int t = ny1 - iy;
      if (n > t) {
	 n = t;
      }
      start[j] = iy;
      nmask[j] = n;
      double norm = 0.0;
      double y = dym;
      for (int k = 0; k < n; k++, y += dstepi) {
	 const double pval = INTERPF(y);
	 mask[j][k] = pval;
	 norm += pval;
      }
      norm = (norm > 0.0) ? 1.0/norm : dstepi;
      for (int k = 0; k < n; k++) {
	 mask[j][k] *= norm;
      }
   }
   
   /* Make the interpolation in y  and transpose once again */
   for (int i = 0; i < nx2; i++) {
      for (int j = 0; j < ny2; j++) {
	 float val = 0.0; 
	 for (int k = 0; k < nmask[j]; k++) {
	    val += pix_tmp[start[j] + k][i]*mask[j][k];
	 }
	 pix2[iys2 + j][ixs2 + i] = val;
      }
   }
   
   /* Free memory */
   free(s_pix_tmp);
   free(pix_tmp);
   free(pix1);
   free(pix2);
   free(s_mask);
   free(mask);
   free(nmask);
   free(start);
   
   return RETURN_OK;
}
/****** vignet_resample_pixel ******************************************************
PROTO	int	vignet_resample(float *pix1, int w1, int h1,
		float *pix2, int w2, int h2, double dx, double dy, float step2,
		float stepi)
PURPOSE	Scale and shift a small image consisting of delta functions
	Image parts which lie outside boundaries are set to 0.

INPUT	Input raster,
	input raster width,
	input raster height,
	output raster,
	output raster width,
	output raster height,
	shift in x,
	shift in y,
	output pixel scale.	
OUTPUT	RETURN_ERROR if the images do not overlap, RETURN_OK otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	13/09/2010
 ***/

int
vignet_resample_pixel(const float *s_pix1, const int w1, const int h1, /* input */
		      float *s_pix2, const int w2, const int h2,       /* output */
		      const double dx, const double dy,
		      const float step2,
		      float stepi)
{
   static float	*stat_s_pix2 = NULL;	/* save an old version of s_pix2.  Scary! */

   if (stepi <= 0.0) {
      stepi = 1.0;
   }
   const double dstepi = 1.0/stepi;
   const double mx1 = w1/2;		/* Im1 center x-coord*/
   const double mx2 = w2/2;		/* Im2 center x-coord*/
   double xs1 = mx1 + dx - mx2*step2;	/* Im1 start x-coord */

   if ((int)xs1 >= w1) {
      return RETURN_ERROR;
   }
   int ixs2 = 0;			/* Int part of Im2 start x-coord */
   if (xs1 < 0.0) {
      const int dix2 = 1 - xs1/step2;
      /*-- Simply leave here if the images do not overlap in x */
      if (dix2 >= w2) {
	 return RETURN_ERROR;
      }
      ixs2 += dix2;
      xs1 += dix2*step2;
   }
   int nx2 = (w1 - 1 - xs1)/step2 + 1;	/* nb of interpolated Im2 pixels along x */
   const int ix2 = w2 - ixs2;
   if (nx2 > ix2) {
      nx2 = ix2;
   }
   if (nx2 <= 0) {
      return RETURN_ERROR;
   }
   const double my1 = h1/2;		/* Im1 center y-coord */
   const double my2 = h2/2;		/* Im2 center y-coord */
   double ys1 = my1 + dy - my2*step2;	/* Im1 start y-coord */
   if ((int)ys1 >= h1) {
      return RETURN_ERROR;
   }
   int iys2 = 0;			/* Int part of Im2 start y-coord */
   if (ys1 < 0.0) {
      const int diy2 = 1 - ys1/step2;
      /*-- Simply leave here if the images do not overlap in y */
      if (diy2 >= h2) {
	 return RETURN_ERROR;
      }
      iys2 += diy2;
      ys1 += diy2*step2;
   }
   int ny2 = (h1 - 1 - ys1)/step2 + 1;	/* nb of interpolated Im2 pixels along y */
   const int iy2 = h2 - iys2;
   if (ny2 > iy2) {
      ny2 = iy2;
   }
   if (ny2 <= 0) {
      return RETURN_ERROR;
   }
   
   /* Set the yrange for the x-resampling with some margin for interpolation */
   int iys1a = ys1;			   /* Int part of Im1 start y-coord with margin */
   const int hmh = (INTERPW/2)/dstepi + 2; /* Interpolant start */
   const int interph = 2*hmh;
   const int hmw = (INTERPW/2)/dstepi + 2;
   const int interpw =  2*hmw;
   iys1a -= hmh;
   if (iys1a < 0) {
      iys1a = 0;
   }
   int ny1 = (ys1 + ny2*step2) + interpw - hmh;	/* Interpolated Im1 y size */
   if (ny1 > h1) {				/* with margin */
      ny1 = h1;
   }
   /* Express everything relative to the effective Im1 start (with margin) */
   ny1 -= iys1a;
   ys1 -= (double)iys1a;

   /* Initialize destination buffer to zero if s_pix2 != NULL */
   if (!s_pix2) {
      assert(stat_s_pix2 != NULL);
      s_pix2 = stat_s_pix2;
   } else {
      memset(s_pix2, 0, (size_t)(w2*h2)*sizeof(float));
      stat_s_pix2 = s_pix2;
   }

   /* Allocate 2-D arrays for data */
   const float **pix1; QMALLOC(pix1, const float *, h1);
   for (int k = 0; k < h1; k++) {
      pix1[k] = &s_pix1[k*w1];
   }

   float **pix2; QMALLOC(pix2, float *, h2);
   for (int k = 0; k < h2; k++) {
      pix2[k] = &s_pix2[k*w2];
   }

   /* Allocate interpolant stuff for the x direction */
   double *s_mask; QMALLOC(s_mask, double, nx2*interpw); /* Storage for interpolation masks */
   double **mask; QMALLOC(mask, double *, nx2); /* interpolation masks */
   for (int k = 0; k < nx2; k++) {
      mask[k] = &s_mask[k*interpw];
   }

   int *nmask;   QMALLOC(nmask, int, nx2);	     /* Interpolation mask sizes */
   int *start;   QMALLOC(start, int, nx2);	     /* Int part of Im1 conv starts */

   /* Compute the local interpolant and data starting points in x */
   double x1 = xs1;
   for (int j = 0; j < nx2; j++, x1 += step2) {
      const int ix1 = x1;
      int ix = ix1 - hmw;
      double dxm = (ix1 - x1 - hmw)*dstepi;/* starting point in the interp. func */
      int n;
      if (ix < 0) {
	 n = interpw + ix;
	 dxm -= (double)ix*dstepi;
	 ix = 0;
      } else {
	 n = interpw;
      }
      int t = w1 - ix;
      if (n > t) {
	 n = t;
      }
      start[j] = ix;
      nmask[j] = n;
      //double norm = 0.0;
      double x = dxm;
#define	INTERPF_LINEAR_DOWN(x)	((fabs(x) > step2) ? 0.0 : step2 - fabs(x))

      for (int i = 0; i < n; i++, x += dstepi) {
	 double pval = INTERPF_LINEAR_DOWN(x);
	 mask[j][i] = pval;
	 //norm += pval;
      }
#if 0
      norm = (norm > 0.0) ? 1.0/norm : dstepi;
      for (int i = 0; i < n; i++) {
	 mask[j][i] *= norm;
      }
#endif
   }

   float *s_pix_tmp; QCALLOC(s_pix_tmp, float, nx2*ny1); /* Intermediary frame-buffer */
   float **pix_tmp; QCALLOC(pix_tmp, float *, ny1);
   for (int k = 0; k < ny1; k++) {
      pix_tmp[k] = &s_pix_tmp[k*nx2];
   }

   /* Make the interpolation in x (this includes transposition) */
   for (int i = 0; i < ny1; i++) {
      for (int j = 0; j < nx2; j++) {
	 float val = 0.0; 
	 for (int k = 0; k < nmask[j]; k++) {
	    const double pval = pix1[iys1a + i][start[j] + k];
	    val += pval*mask[j][k];
	 }
	 pix_tmp[i][j] = val;
      }
   }
  
   /* Reallocate interpolant stuff for the y direction */
   QREALLOC(s_mask, double, ny2*interph);	/* Interpolation masks */
   QREALLOC(mask, double *, ny2*interph);	/* Interpolation masks */
   for (int k = 0; k < ny2; k++) {
      mask[k] = &s_mask[k*interph];
   }
   
   QREALLOC(nmask, int, ny2);		/* Interpolation mask sizes */
   QREALLOC(start, int, ny2);		/* Int part of Im1 conv starts */
   
   /* Compute the local interpolant and data starting points in y */
   double y1 = ys1;
   for (int j = 0; j < ny2; j++, y1 += step2) {
      const int iy1 = y1;
      int iy = iy1 - hmh;
      double dym = (iy1 - y1 - hmh)*dstepi;/* starting point in the interp. func */
      int n;
      if (iy < 0) {
	 n = interph + iy;
	 dym -= (double)iy*dstepi;
	 iy = 0;
      } else {
	 n = interph;
      }
      const int t = ny1 - iy;
      if (n > t) {
	 n = t;
      }
      start[j] = iy;
      nmask[j] = n;
      //double norm = 0.0;
      double y = dym;
      for (int k = 0; k < n; k++, y += dstepi) {
	 const double pval = INTERPF_LINEAR_DOWN(y);
	 mask[j][k] = pval;
	 //norm += pval;
      }
      #if 0
      norm = (norm > 0.0) ? 1.0/norm : dstepi;
      for (int k = 0; k < n; k++) {
	 mask[j][k] *= norm;
      }
      #endif
   }
   
   /* Make the interpolation in y  and transpose once again */
   for (int i = 0; i < nx2; i++) {
      for (int j = 0; j < ny2; j++) {
	 float val = 0.0; 
	 for (int k = 0; k < nmask[j]; k++) {
	    val += pix_tmp[start[j] + k][i]*mask[j][k];
	 }
	 pix2[iys2 + j][ixs2 + i] = val;
      }
   }
   
   /* Free memory */
   free(s_pix_tmp);
   free(pix_tmp);
   free(pix1);
   free(pix2);
   free(s_mask);
   free(mask);
   free(nmask);
   free(start);
   
   return RETURN_OK;
}

/******************************** vignet_copy ********************************/
/*
Copy a small part of the image. Image parts which lie outside boundaries are
set to 0.
*/
int     vignet_copy(float *pix1, int w1, int h1,
		float *pix2, int w2, int h2, int idx, int idy, vigopenum vigop)
  {
   int          x,y, xmin,ymin, nx,ny, off1,off2;

  if (vigop==VIGNET_CPY)
/*-- First put the pix2 background to zero */
    memset(pix2, 0, (size_t)(w2*h2)*sizeof(float));

/* Set the image boundaries */
  ymin = h2/2+idy-h1/2;
  if ((ny=h2-ymin)>h1)
    ny = h1;
  else if (ny<=0)
    return RETURN_ERROR;
  if (ymin<0)
    {
    pix1 -= ymin*w1;
    ny += ymin;
    }
  else
    pix2 += ymin*w2;

  xmin = w2/2+idx-w1/2;
  if ((nx=w2-xmin)>w1)
    nx = w1;
  else if (nx<=0)
    return RETURN_ERROR;
  if (xmin<0)
    {
    pix1 -= xmin;
    nx += xmin;
    }
  else
    pix2 += xmin;

/* Offsets */
  off1 = w1-nx;
  off2 = w2-nx;
/* Copy the right pixels to the destination */
  switch(vigop)
    {
    case VIGNET_CPY:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) = *(pix1++);
      break;

    case VIGNET_ADD:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) += *(pix1++);
      break;

    case VIGNET_SUB:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) -= *(pix1++);
      break;

    case VIGNET_MUL:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) *= *(pix1++);
      break;

    case VIGNET_DIV:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--; pix2++)
          if (*pix1)
            *pix2 /= *(pix1++);
          else
            *pix2 = (*pix2>0.0)? BIG:-BIG;

    default:
      error(EXIT_FAILURE, "*Internal Error*: unknown operation in ",
			"vignet_copy()");
      
      break;
    }

  return RETURN_OK;
  }


 /**************************** vignet_aperflux******************************/
/*
Compute the total flux within a circular aperture.
*/
float	vignet_aperflux(float *ima, float *var, int w, int h,
			float dxc, float dyc, float aper,
			float gain, float backnoise, float *fluxvar)

  {
   float		*imat,*vart,
			r2, raper,raper2, rintlim,rintlim2,rextlim2, mx,my,
			dx,dx1,dy,dy2, pix,pvar,invbacknoise2, invgain, area,
			offsetx,offsety,scalex,scaley,scale2, locarea, vthresh;
   double		tv, sigtv;
   int			x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy,
			gainflag = 0;
   long			pos;

  pvar = backnoise*backnoise;
  invbacknoise2 = pvar>0.0? 1.0/pvar : 0.0;
  invgain = gain>0.0? 1.0/gain : 0.0;
/* Integration radius */
  raper = aper/2.0;
  raper2 = raper*raper;
/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
  rintlim = raper - 0.75;
  rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
  rextlim2 = (raper + 0.75)*(raper + 0.75);
  tv = sigtv = area = 0.0;
  scaley = scalex = 1.0/APER_OVERSAMP;
  scale2 = scalex*scaley;
  offsetx = 0.5*(scalex-1.0);
  offsety = 0.5*(scaley-1.0);
  vthresh = BIG/2.0;
  mx = dxc + (float)(w/2);
  my = dyc + (float)(h/2);

  xmin = (int)(mx-raper+0.499999);
  xmax = (int)(mx+raper+1.499999);
  ymin = (int)(my-raper+0.499999);
  ymax = (int)(my+raper+1.499999);
  if (xmin < 0 || xmax > w || ymin < 0 || ymax > h)
    return 0.0;

  for (y=ymin; y<ymax; y++)
    {
    imat = ima + (pos = y*w + xmin);
    vart = var + pos;
    for (x=xmin; x<xmax; x++, imat++, vart++)
      {
      dx = x - mx;
      dy = y - my;
      if ((r2=dx*dx+dy*dy) < rextlim2)
        {
        if (r2> rintlim2)
          {
          dx += offsetx;
          dy += offsety;
          locarea = 0.0;
          for (sy=APER_OVERSAMP; sy--; dy+=scaley)
            {
            dx1 = dx;
            dy2 = dy*dy;
            for (sx=APER_OVERSAMP; sx--; dx1+=scalex)
              if (dx1*dx1+dy2<raper2)
                locarea += scale2;
            }
          }
        else
          locarea = 1.0;
        area += locarea;
/*------ Here begin tests for pixel and/or weight overflows. Things are a */
/*------ bit intricated to have it running as fast as possible in the most */
/*------ common cases */
        if ((pix=*imat)<=-BIG || (var && (pvar=*vart)>=vthresh))
          {
          if ((x2=(int)(2*mx+0.49999-x))>=0 && x2<w
		&& (y2=(int)(2*my+0.49999-y))>=0 && y2<h
		&& (pix=*(imat + (pos = y2*w + x2)))>-BIG)
            {
            if (var)
              {
              pvar = *(var + pos);
              if (pvar>=vthresh)
                pix = pvar = 0.0;
              }
            }
          else
            {
            pix = 0.0;
            if (var)
              pvar = 0.0;
            }
          }
        tv += locarea*pix;
        sigtv += locarea*pvar;
        if (gainflag && pix>0.0 && gain>0.0)
          sigtv += pix*invgain*pvar*invbacknoise2;
        }
      }
    }

  if (!gainflag && tv>0.0)
    sigtv += tv*invgain;

  if (fluxvar)
    *fluxvar = sqrt(sigtv);

  return tv;
  }



