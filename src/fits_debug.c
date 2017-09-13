#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fits_debug.h"

#define FITS_SIZE 2880

/*****************************************************************************/
/*
 * Write a string value
 */
static int
write_card_s(fd, ncard, record, keyword, val, commnt)
FILE *fd;
int ncard;
char *record;
char *keyword;
char *val;
char *commnt;
{
   char *card = &record[80*ncard];
   char value[20];			/* blank-padded val, if needed */

   if(strlen(val) < 8) {		/* FITS requires at least 8 chars */
      sprintf(value, "%-8s", val);
      val = value;
   }

   if(*keyword == '\0' || !strcmp(keyword,"COMMENT") || !strcmp(keyword,"END")
		       || !strcmp(keyword,"HISTORY")) {
      if(commnt[0] != '\0') {
	 fprintf(stderr, "You can't add a comment to a COMMENT, END, or HISTORY card\n");
      }
      sprintf(card,"%-8.8s%-72s",keyword,val);
   } else {
      sprintf(card,"%-8.8s= '%s' %c%-*s",keyword,val,
	      (commnt[0] == '\0' ? ' ' : '/'),
	      (int)(80 - 14 - strlen(val)), commnt);
   }
/*
 * Write record if full
 */
   if(++ncard == 36) {
      if(fwrite(record, FITS_SIZE, 1, fd) != 1) {
	 fprintf(stderr, "Cannot write header record");
      }
      ncard = 0;
   }
   
   return(ncard);
}   

/*****************************************************************************/
/*
 * Write an integer (%d)
 */
static int
write_card_d(fd, ncard, record, keyword, val, commnt)
FILE *fd;
int ncard;
char *record;
char *keyword;
int val;
char *commnt;
{
   char *card = &record[80*ncard];

   sprintf(card,"%-8.8s= %20d %c%-48s",keyword,val,
	   		(commnt[0] == '\0' ? ' ' : '/'),commnt);
/*
 * Write record if full
 */
   if(++ncard == 36) {
      if(fwrite(record, FITS_SIZE, 1, fd) != 1) {
	 fprintf(stderr, "Cannot write header record\n");
      }
      ncard = 0;
   }
   
   return(ncard);
}   

/*****************************************************************************/
/*
 * Write a logical value
 */
static int
write_card_l(fd, ncard, record, keyword, val, commnt)
FILE *fd;
int ncard;
char *record;
char *keyword;
char *val;
char *commnt;
{
   char *card = &record[80*ncard];

   if(strcmp(val,"T") != 0 && strcmp(val,"F") != 0) {
      fprintf(stderr, "Invalid logical %s for keyword %s\n",val,keyword);
      val = "?";
   }
   sprintf(card,"%-8.8s= %20s %c%-48s",keyword,val,
	   		(commnt[0] == '\0' ? ' ' : '/'),commnt);
/*
 * Write record if full
 */
   if(++ncard == 36) {
      if(fwrite(record, FITS_SIZE, 1, fd) != 1) {
	 fprintf(stderr, "Cannot write header record\n");
      }
      ncard = 0;
   }
   
   return(ncard);
}   

/*
 * Byte swap ABCDABCD -> DCBADCBA in place (e.g. sun <--> vax)
 */
static void
swap_4(char *arr,			// array to swap
       const int n)			// number of bytes
{
   if (n%4 != 0) {
      fprintf(stderr, "Attempt to byte swap non-multiple of 4 bytes: %d", n);
   }

    for (char *end = arr + n; arr < end; arr += 4) {
        char t = arr[0];
        arr[0] = arr[3];
        arr[3] = t;
        t = arr[1];
        arr[1] = arr[2];
        arr[2] = t;
    }
}

void
write_fits_image(const char *filename,
		 float *data,
		 int nx,
		 int ny)
{
   FILE *fd;				/* file descriptor */
   
   if((fd = fopen(filename,"w")) == NULL) {
      fprintf(stderr, "I'm afraid that I can't open %s\n", filename);
      return;
   }

   char record[FITS_SIZE + 1];		/* write buffer */
   const int swap_bytes = 1;		/* e.g. X86 */

   int ncard = 0;			/* number of cards written */
   ncard = write_card_l(fd,ncard,record,"SIMPLE","T","");
   ncard = write_card_d(fd,ncard,record,"BITPIX",-32,"");
   ncard = write_card_d(fd,ncard,record,"NAXIS",2,"");
   ncard = write_card_d(fd,ncard,record,"NAXIS1",nx,"");
   ncard = write_card_d(fd,ncard,record,"NAXIS2",ny,"");
   ncard = write_card_s(fd,ncard,record,"END","","");
      
   while(ncard != 0) {
      ncard = write_card_s(fd,ncard,record,"","","");
   }

   float *buff = NULL;			/* buffer for writing data rows */
   if(swap_bytes) {
      if((buff = (float *)malloc(nx*sizeof(float))) == NULL) {
	 fprintf(stderr, "Can't allocate storage for buff\n");
	 fclose(fd);
	 return;
      }
   }

   int i = 0;
   for(; i < ny; i++) {
      if(swap_bytes) {
	 memcpy(buff, &data[i*nx], nx*sizeof(float));
	 swap_4((char *)buff, nx*sizeof(float));
      } else {
	 buff = &data[i*nx];
      }

      if(fwrite(buff, nx*sizeof(float), 1, fd) != 1) {
	 fprintf(stderr, "Error writing %dth row of image\n", i + 1);
	 break;
      }
   }
   if(i == ny) {		/* no write errors */
      int nbyte = nx*ny*sizeof(float);
      if(nbyte%FITS_SIZE != 0) {
	 nbyte = FITS_SIZE - nbyte%FITS_SIZE;
	 memset(record, ' ', nbyte);
	 if(fwrite(record, nbyte, 1, fd) != 1) {
	    fprintf(stderr, "error padding file to multiple of fits block size\n");
	 }	       
      }
   }

   if(buff) {
      free(buff);
   }
   fclose(fd);
}
