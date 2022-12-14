/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: oftware Author:  Kevin L. Bolling
           Updated by Jeff Hinrichs
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
//#include <malloc.h>
#include <math.h>

#define FLOAT32 4
#define FLOAT64 8
#define INT32 4
#define INT64 8

#define ASCII_F 0
#define IEEE_F  1

long n_nodes, n_cells, n_surf, n_faces, i, str_ncells;
static int charsize_out = 8;
short structflag = 0;
FILE *fp;
char tmpname[9];
int ieerror;
int typei = 1;
int typef = 2;
int typed = 8;
int num0 = 0;
int num1 = 1;
int iflag64 = 0, rflag64 = 0;
int filetype;

void write_ascii_float(long num, float *val);
void write_ascii_double(long num, double *val);
void write_ascii_int(long num, int *val);
void write_ascii_long(long num, long *val);


/* --------------------------------------------------------- */

void gmvwrite_openfile(char filenam[])
{
   fp = fopen(filenam, "w");

   /*  Write header.  */
   strcpy(tmpname,"gmvinput");
   fwrite(tmpname,sizeof(char),8,fp);
   strcpy(tmpname, "ieee    ");
   fwrite(tmpname,sizeof(char),8,fp);
   filetype = IEEE_F;
}

/* --------------------------------------------------------- */

void gmvwrite_openfile_ascii(char filenam[])
{
   fp = fopen(filenam, "w");

   /*  Write header.  */
   fprintf(fp,"gmvinput ascii\n");
   filetype = ASCII_F;
   charsize_out = 32;
}

/* --------------------------------------------------------- */

void gmvwrite_openfile_ir(char filenam[], int isize, int rsize)
{

   /*  Check if this machine can address 64bit integers (isize = 8).  */
   if (isize == 8)
     {
      if (sizeof(long) < 8)
        {
         printf("gmvwrite error, cannot write 8 byte integers on this machine.\n");
         exit(0);
        }
     }

   fp = fopen(filenam, "w");

   /*  Write header.  */
   strcpy(tmpname, "gmvinput");
   fwrite(tmpname, sizeof(char), 8, fp);
   if (isize == 4 && rsize == 8)
     {
      strcpy(tmpname, "ieeei4r8");
      fwrite(tmpname, sizeof(char), 8, fp);
      rflag64 = 1;
     }
   else if (isize == 8 && rsize == 4)
     {
      strcpy(tmpname, "ieeei8r4");
      fwrite(tmpname, sizeof(char), 8, fp);
      iflag64 = 1;
     }
   else if (isize == 8 && rsize == 8)
     {
      strcpy(tmpname, "ieeei8r8");
       fwrite(tmpname, sizeof(char), 8, fp);
      iflag64 = 1;
      rflag64 = 1;
     }
   else  /* just do the normal "ieee    " thing */
     {
      strcpy(tmpname, "ieee    ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   filetype = IEEE_F;
}

/* --------------------------------------------------------- */

void gmvwrite_openfile_cxir(char filenam[], int isize, int rsize)
{

   /*  Check if this machine can address 64bit integers (isize = 8).  */
   if (isize == 8)
     {
      if (sizeof(long) < 8)
        {
         printf("gmvwrite error, cannot write 8 byte integers on this machine.\n");
         exit(0);
        }
     }

   fp = fopen(filenam, "w");
   charsize_out = 32;

   /*  Write header.  */
   strcpy(tmpname, "gmvinput");
   fwrite(tmpname, sizeof(char), 8, fp);
   if (isize == 4 && rsize == 8)
     {
      strcpy(tmpname, "iecxi4r8");
      fwrite(tmpname, sizeof(char), 8, fp);
      rflag64 = 1;
     }
   else if (isize == 8 && rsize == 4)
     {
      strcpy(tmpname, "iecxi8r4");
      fwrite(tmpname, sizeof(char), 8, fp);
      iflag64 = 1;
     }
   else if (isize == 8 && rsize == 8)
     {
      strcpy(tmpname, "iecxi8r8");
       fwrite(tmpname, sizeof(char), 8, fp);
      iflag64 = 1;
      rflag64 = 1;
     }
   else  /* just do the normal "iecxi4r4" thing */
     {
      strcpy(tmpname, "iecxi4r4");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   filetype = IEEE_F;
}

/* --------------------------------------------------------- */

void gmvwrite_openfile_ir_ascii(char filenam[], int isize, int rsize)
{

   /*  Check if this machine can address 64bit integers (isize = 8).  */
   if (isize == 8)
     {
      if (sizeof(long) < 8)
        {
         printf("gmvwrite error, cannot write 8 byte integers on this machine.\n");
         exit(0);
        }
     }

   fp = fopen(filenam, "w");
   charsize_out = 32;

   /* write header */
   fprintf(fp,"gmvinput ascii\n");

   if (isize == 4 && rsize == 8)
     {
      rflag64 = 1;
     }
   else if (isize == 8 && rsize == 4)
     {
      iflag64 = 1;
     }
   else if (isize == 8 && rsize == 8)
     {
      iflag64 = 1;
      rflag64 = 1;
     }
   filetype = ASCII_F;
}

/* --------------------------------------------------------- */

void gmvwrite_closefile(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endgmv  ");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"endgmv\n");
   fclose(fp);
}

/* --------------------------------------------------------- */

void gmvwrite_nodes_fromfile(char *filename, long nndes)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("nodes   fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"nodes   fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);

   n_nodes = nndes;

   free(tmpbuf);
   return;
 }

/* --------------------------------------------------------- */

void gmvwrite_node_data(void *nndes, void *x, void *y, void *z)
 {
   int tmpnndes;
   float *tempx=NULL, *tempy=NULL, *tempz=NULL;
   double *tempx64=NULL, *tempy64=NULL, *tempz64=NULL;
   long tmpnndes64, longcount, totnodes;

   strcpy(tmpname, "nodes   ");

   if (iflag64)
      totnodes = tmpnndes64 = *((long *) nndes);
   else
      totnodes = tmpnndes = *((int *) nndes);

   if (rflag64)
     {
      tempx64 = (double *) malloc(sizeof(double) * totnodes);
      tempy64 = (double *) malloc(sizeof(double) * totnodes);
      tempz64 = (double *) malloc(sizeof(double) * totnodes);
     }
   else
     {
      tempx = (float*) malloc(sizeof(float) * totnodes);
      tempy = (float*) malloc(sizeof(float) * totnodes);
      tempz = (float*) malloc(sizeof(float) * totnodes);
     }
    if (rflag64)
      {
       for (longcount = 0; longcount < totnodes; longcount++)
	 {
	  tempx64[longcount] = *((double *) x + longcount);
          tempy64[longcount] = *((double *) y + longcount);
	  tempz64[longcount] = *((double *) z + longcount);
         }
      }
    else
      {
       for (longcount = 0; longcount < totnodes; longcount++)
         {
          tempx[longcount] = *((float *) x + longcount);
          tempy[longcount] = *((float *) y + longcount);
          tempz[longcount] = *((float *) z + longcount);
	 }
      }

   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   else
      fprintf(fp,"nodes  ");

   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnndes64, INT64, 1, fp);
      else
         fprintf(fp,"%ld\n",tmpnndes64);
      n_nodes = tmpnndes64;
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnndes, INT32, 1, fp);
      else
         fprintf(fp,"%d\n",tmpnndes);
      n_nodes = tmpnndes;
     }
   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx64, FLOAT64, n_nodes, fp);
         fwrite(tempy64, FLOAT64, n_nodes, fp);
         fwrite(tempz64, FLOAT64, n_nodes, fp);
        }
      else
        {
         write_ascii_double(n_nodes, tempx64);
         write_ascii_double(n_nodes, tempy64);
         write_ascii_double(n_nodes, tempz64);
        }
      free(tempx64), free(tempy64), free(tempz64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx, FLOAT32, n_nodes, fp);
         fwrite(tempy, FLOAT32, n_nodes, fp);
         fwrite(tempz, FLOAT32, n_nodes, fp);
        }
      else
        {
         write_ascii_float(n_nodes, tempx);
         write_ascii_float(n_nodes, tempy);
         write_ascii_float(n_nodes, tempz);
        }
      free(tempx), free(tempy), free(tempz);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_node_data_struct(void *nxv, void *nyv, void *nzv,
			       void *x,  void *y, void *z)
 {
   int alt, count;
   int tmpnxv, tmpnyv, tmpnzv;
   long tmpnxv64, tmpnyv64, tmpnzv64, longcount, ix, iy, iz, icx, icy, icz;
   float *tempx=NULL, *tempy=NULL, *tempz=NULL;
   double *tempx64=NULL, *tempy64=NULL, *tempz64=NULL;

   strcpy(tmpname, "nodes   ");
   alt =-1;

   structflag = 1;

   if (iflag64)
     {
      tmpnxv64 = *((long *) nxv);
      tmpnyv64 = *((long *) nyv);
      tmpnzv64 = *((long *) nzv);
     }
   else
     {
      tmpnxv = *((int *) nxv);
      tmpnyv = *((int *) nyv);
      tmpnzv = *((int *) nzv);
     }

   if (rflag64)
   {
     if (iflag64)
     {
       tempx64 = (double *) malloc(sizeof(double) * tmpnxv64);
       tempy64 = (double *) malloc(sizeof(double) * tmpnyv64);
       tempz64 = (double *) malloc(sizeof(double) * tmpnzv64);
     }
     else
     {
       tempx64 = (double *) malloc(sizeof(double) * tmpnxv);
       tempy64 = (double *) malloc(sizeof(double) * tmpnyv);
       tempz64 = (double *) malloc(sizeof(double) * tmpnzv);
     }
   }
   else
   {
     if (iflag64)
     {
       tempx = (float*)malloc(sizeof(float)*tmpnxv64);
       tempy = (float*)malloc(sizeof(float)*tmpnyv64);
       tempz = (float*)malloc(sizeof(float)*tmpnzv64);
     }
     else
     {
       tempx = (float*)malloc(sizeof(float)*tmpnxv);
       tempy = (float*)malloc(sizeof(float)*tmpnyv);
       tempz = (float*)malloc(sizeof(float)*tmpnzv);
     }
   }

   if (rflag64)
     {
      if (iflag64)
        {
	 for (longcount = 0; longcount < tmpnxv64; longcount++)
	    tempx64[longcount] = *((double *) x + longcount);
	 for (longcount = 0; longcount < tmpnyv64; longcount++)
	    tempy64[longcount] = *((double *) y + longcount);
	 for (longcount = 0; longcount < tmpnzv64; longcount++)
	    tempz64[longcount] = *((double *) z + longcount);
	}
      else
	{
	 for (count = 0; count < tmpnxv; count++)
	    tempx64[count] = *((double *) x + count);
	 for (count = 0; count < tmpnyv; count++)
	    tempy64[count] = *((double *) y + count);
	 for (count = 0; count < tmpnzv; count++)
	    tempz64[count] = *((double *) z + count);
	}
     }
   else
     {
      if (iflag64)
	{
         for(longcount = 0; longcount < tmpnxv64; longcount++)
            tempx[longcount] = *((float *) x + longcount);
         for(longcount = 0; longcount < tmpnyv64; longcount++)
	    tempy[longcount] = *((float *) y + longcount);
	 for(longcount = 0; longcount < tmpnzv64; longcount++)
            tempz[longcount] = *((float *) z + longcount);
	}
      else
	{
         for(count = 0; count < tmpnxv; count ++)
            tempx[count] = *((float *) x + count);
         for(count = 0; count < tmpnyv; count ++)
	    tempy[count] = *((float *) y + count);
	 for(count = 0; count < tmpnzv; count ++)
            tempz[count] = *((float *) z + count);
	}
     }

   if (iflag64)
     {
      ix = tmpnxv64;
      iy = tmpnyv64;
      iz = tmpnzv64;
     }
   else
     {
      ix = tmpnxv;
      iy = tmpnyv;
      iz = tmpnzv;
     }

   /*  Write node header.  */
   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(&alt,INT32,1,fp);
      if (iflag64)
        {
         fwrite(&tmpnxv64, INT64, 1, fp);
         fwrite(&tmpnyv64, INT64, 1, fp);
         fwrite(&tmpnzv64, INT64, 1, fp);
        }
      else
        {
         fwrite(&tmpnxv, INT32, 1, fp);
         fwrite(&tmpnyv, INT32, 1, fp);
         fwrite(&tmpnzv, INT32, 1, fp);
        }
     }
   else
     {
      fprintf(fp,"nodes %d",alt);
      if (iflag64)
        {
         fprintf(fp," %ld %ld %ld\n",tmpnxv64,tmpnyv64,tmpnzv64);
        }
      else
        {
         fprintf(fp," %d %d %d\n",tmpnxv,tmpnyv,tmpnzv);
        }
     }

   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx64, FLOAT64, ix, fp);
         fwrite(tempy64, FLOAT64, iy, fp);
         fwrite(tempz64, FLOAT64, iz, fp);
        }
      else
        {
         write_ascii_double(ix, tempx64);
         write_ascii_double(iy, tempy64);
         write_ascii_double(iz, tempz64);
        }
      free(tempx64), free(tempy64), free(tempz64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx, FLOAT32, ix, fp);
         fwrite(tempy, FLOAT32, iy, fp);
         fwrite(tempz, FLOAT32, iz, fp);
        }
      else
        {
         write_ascii_float(ix, tempx);
         write_ascii_float(iy, tempy);
         write_ascii_float(iz, tempz);
        }
     }

   n_nodes = ix * iy * iz;
   icx = ix - 1;
   icy = iy - 1;
   icz = iz - 1;
   if (icx < 1) icx = 1;
   if (icy < 1) icy = 1;
   if (icz < 1) icz = 1;
   str_ncells = icx * icy * icz;
 }

/* --------------------------------------------------------- */

void gmvwrite_node_data_lstruct(void *nxv, void *nyv, void *nzv,
			        void *x, void *y, void *z)
{
  int alt, tmpnxv, tmpnyv, tmpnzv;
  long tmpnxv64, tmpnyv64, tmpnzv64, nndes, count, ix, iy, iz, icx, icy, icz;
  float *tempx=NULL, *tempy=NULL, *tempz=NULL;
  double *tempx64=NULL, *tempy64=NULL, *tempz64=NULL;

   strcpy(tmpname, "nodes   ");
   alt = -2;

   structflag = 1;

   if (iflag64)
     {
      tmpnxv64 = *((long *) nxv);
      tmpnyv64 = *((long *) nyv);
      tmpnzv64 = *((long *) nzv);
     }
   else
     {
      tmpnxv = *((int *) nxv);
      tmpnyv = *((int *) nyv);
      tmpnzv = *((int *) nzv);
     }

  if (iflag64)
   {
    nndes = tmpnxv64 * tmpnyv64 * tmpnzv64;
    str_ncells = (tmpnxv64-1) * (tmpnyv64-1) * (tmpnzv64-1);
   }
  else
   {
    nndes = tmpnxv * tmpnyv * tmpnzv;
    str_ncells = (tmpnxv-1) * (tmpnyv-1) * (tmpnzv-1);
   }
  if (rflag64)
    {
     tempx64 = (double *) malloc(sizeof(double) * nndes);
     tempy64 = (double *) malloc(sizeof(double) * nndes);
     tempz64 = (double *) malloc(sizeof(double) * nndes);
    }
  else
    {
     tempx = (float*)malloc(sizeof(float)*nndes);
     tempy = (float*)malloc(sizeof(float)*nndes);
     tempz = (float*)malloc(sizeof(float)*nndes);
    }

  if (rflag64)
     for (count = 0; count < nndes; count++)
       {
	tempx64[count] = *((double *) x + count);
        tempy64[count] = *((double *) y + count);
	tempz64[count] = *((double *) z + count);
       }
     else
       for (count = 0; count < nndes; count++)
         {
          tempx[count] = *((float *) x + count);
          tempy[count] = *((float *) y + count);
          tempz[count] = *((float *) z + count);
         }

   if (iflag64)
     {
      ix = tmpnxv64;
      iy = tmpnyv64;
      iz = tmpnzv64;
     }
   else
     {
      ix = tmpnxv;
      iy = tmpnyv;
      iz = tmpnzv;
     }

   /*  Write node x,y,z's for logicaly structured grids.  */
   n_nodes = nndes;
   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(&alt,INT32,1,fp);
     }
   else
      fprintf(fp,"nodes %d",alt);

   if (iflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tmpnxv64, INT64, 1, fp);
         fwrite(&tmpnyv64, INT64, 1, fp);
         fwrite(&tmpnzv64, INT64, 1, fp);
        }
      else
        {
         fprintf(fp," %ld %ld %ld\n",tmpnxv64,tmpnyv64,tmpnzv64);
        }
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tmpnxv, INT32, 1, fp);
         fwrite(&tmpnyv, INT32, 1, fp);
         fwrite(&tmpnzv, INT32, 1, fp);
        }
      else
        {
         fprintf(fp," %d %d %d\n",tmpnxv,tmpnyv,tmpnzv);
        }
     }
   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx64, FLOAT64, nndes, fp);
         fwrite(tempy64, FLOAT64, nndes, fp);
         fwrite(tempz64, FLOAT64, nndes, fp);
        }
      else
        {
         write_ascii_double(nndes, tempx64);
         write_ascii_double(nndes, tempy64);
         write_ascii_double(nndes, tempz64);
        }
      free(tempx64), free(tempy64), free(tempz64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx,FLOAT32,nndes,fp);
         fwrite(tempy,FLOAT32,nndes,fp);
         fwrite(tempz,FLOAT32,nndes,fp);
        }
      else
        {
         write_ascii_float(nndes, tempx);
         write_ascii_float(nndes, tempy);
         write_ascii_float(nndes, tempz);
        }
      free(tempx), free(tempy), free(tempz);
     }

   n_nodes = ix * iy * iz;
   icx = ix - 1;
   icy = iy - 1;
   icz = iz - 1;
   if (icx < 1) icx = 1;
   if (icy < 1) icy = 1;
   if (icz < 1) icz = 1;
   str_ncells = icx * icy * icz;
}

/* --------------------------------------------------------- */

void gmvwrite_node_data_amr(int nxc, int nyc, int nzc, void *x0, void *y0,
		            void *z0, void *dx, void *dy, void *dz)
{
  int tnxc, tnyc, tnzc, neg3 = -3;
  float tx0, ty0, tz0, tdx, tdy, tdz;
  double tx064, ty064, tz064, tdx64, tdy64, tdz64;

  n_nodes = nxc * nyc * nzc;

  tnxc = nxc;  tnyc = nyc;  tnzc = nzc;
  if (rflag64)
    {
     tx064 = *((double *) x0);
     ty064 = *((double *) y0);
     tz064 = *((double *) z0);
     tdx64 = *((double *) dx);
     tdy64 = *((double *) dy);
     tdz64 = *((double *) dz);
    }
  else
    {
     tx0 = *((float *) x0);
     ty0 = *((float *) y0);
     tz0 = *((float *) z0);
     tdx = *((float *) dx);
     tdy = *((float *) dy);
     tdz = *((float *) dz);
    }

   /*  Write node header data.  */
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "nodes   ");
      fwrite(tmpname, sizeof(char), 8, fp);
      fwrite(&neg3, INT32, 1, fp);
      fwrite(&tnxc, INT32, 1, fp);
      fwrite(&tnyc, INT32, 1, fp);
      fwrite(&tnzc, INT32, 1, fp);
     }
   else
     {
      fprintf(fp,"%s %d %d %d %d\n", tmpname, neg3, tnxc, tnyc, tnzc);
     }
   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tx064, FLOAT64, 1, fp);
         fwrite(&ty064, FLOAT64, 1, fp);
         fwrite(&tz064, FLOAT64, 1, fp);
         fwrite(&tdx64, FLOAT64, 1, fp);
         fwrite(&tdy64, FLOAT64, 1, fp);
         fwrite(&tdz64, FLOAT64, 1, fp);
        }
      else
        {
         fprintf(fp,"%lg %lg %lg %lg %lg %lg\n", tx064, ty064, tz064,
                 tdx64, tdy64, tdz64);
        }
     }
  else
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tx0, FLOAT32, 1, fp);
         fwrite(&ty0, FLOAT32, 1, fp);
         fwrite(&tz0, FLOAT32, 1, fp);
         fwrite(&tdx, FLOAT32, 1, fp);
         fwrite(&tdy, FLOAT32, 1, fp);
         fwrite(&tdz, FLOAT32, 1, fp);
        }
      else
        {
         fprintf(fp,"%g %g %g %g %g %g\n", tx0, ty0, tz0, tdx, tdy, tdz);
        }
     }
}

/* --------------------------------------------------------- */

void gmvwrite_nodev_fromfile(char *filename, long nndes)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("nodev   fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"nodev   fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
   fprintf(fp,"%s\n",tmpbuf);

   n_nodes = nndes;

   free(tmpbuf);
   return;
 }

/* --------------------------------------------------------- */

void gmvwrite_nodev_data(void *nndes, void *x, void *y, void *z)
 {
   int tmpnndes, i, j;
   float *tempx=NULL, *tempy=NULL, *tempz=NULL, *tempxyz=NULL;
   double *tempx64=NULL, *tempy64=NULL, *tempz64=NULL, *tempxyz64=NULL;
   long tmpnndes64, longcount, totnodes;

   strcpy(tmpname, "nodev   ");

   if (iflag64)
      totnodes = tmpnndes64 = *((long *) nndes);
   else
      totnodes = tmpnndes = *((int *) nndes);

   if (rflag64)
     {
      tempx64 = (double *) malloc(sizeof(double) * totnodes);
      tempy64 = (double *) malloc(sizeof(double) * totnodes);
      tempz64 = (double *) malloc(sizeof(double) * totnodes);
     }
   else
     {
      tempx = (float*) malloc(sizeof(float) * totnodes);
      tempy = (float*) malloc(sizeof(float) * totnodes);
      tempz = (float*) malloc(sizeof(float) * totnodes);
     }

    if (rflag64)
      {
       for (longcount = 0; longcount < totnodes; longcount++)
	 {
	  tempx64[longcount] = *((double *) x + longcount);
	  tempy64[longcount] = *((double *) y + longcount);
	  tempz64[longcount] = *((double *) z + longcount);
	 }
      }
    else
      {
       for (longcount = 0; longcount < totnodes; longcount++)
         {
          tempx[longcount] = *((float *) x + longcount);
          tempy[longcount] = *((float *) y + longcount);
          tempz[longcount] = *((float *) z + longcount);
	 }
      }

   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   else
      fprintf(fp,"%s",tmpname);

   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnndes64, INT64, 1, fp);
      else
         fprintf(fp,"%ld\n",tmpnndes64);
      n_nodes = tmpnndes64;
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnndes, INT32, 1, fp);
      else
         fprintf(fp,"%d\n",tmpnndes);
      n_nodes = tmpnndes;
     }

   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         tempxyz64 = (double *) malloc(sizeof(double) * totnodes * 3);
         j = 0;
         for (i = 0; i < totnodes; i++)
           {
            tempxyz64[j] = tempx64[i];  j++;
            tempxyz64[j] = tempy64[i];  j++;
            tempxyz64[j] = tempz64[i];  j++;
           }
         fwrite(tempxyz64, FLOAT64, 3*totnodes, fp);
         free(tempxyz64);
        }
      else
        {
         for (i = 0; i < totnodes; i++)
           {
            fprintf(fp,"%lg %lg %lg \n",tempx64[i], tempy64[i], tempz64[i]);
           }
        }
      free(tempx64), free(tempy64), free(tempz64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         tempxyz = (float *) malloc(sizeof(float) * totnodes * 3);
         j = 0;
         for (i = 0; i < totnodes; i++)
           {
            tempxyz[j] = tempx[i];  j++;
            tempxyz[j] = tempy[i];  j++;
            tempxyz[j] = tempz[i];  j++;
           }
         fwrite(tempxyz, FLOAT32, 3*totnodes, fp);
         free(tempxyz);
        }
      else
        {
         for (i = 0; i < totnodes; i++)
           {
            fprintf(fp,"%lg %lg %lg \n",tempx[i], tempy[i], tempz[i]);
           }
        }
      free(tempx), free(tempy), free(tempz);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_nodev_data_lstruct(void *nxv, void *nyv, void *nzv,
			        void *x, void *y, void *z)
{
  int alt, tmpnxv, tmpnyv, tmpnzv, i, j;
  long tmpnxv64, tmpnyv64, tmpnzv64, nndes, count;
  float *tempx=NULL, *tempy=NULL, *tempz=NULL, *tempxyz=NULL;
  double *tempx64=NULL, *tempy64=NULL, *tempz64=NULL, *tempxyz64=NULL;

   strcpy(tmpname, "nodev   ");
   alt = -2;

   structflag = 1;

   if (iflag64)
     {
      tmpnxv64 = *((long *) nxv);
      tmpnyv64 = *((long *) nyv);
      tmpnzv64 = *((long *) nzv);
     }
   else
     {
      tmpnxv = *((int *) nxv);
      tmpnyv = *((int *) nyv);
      tmpnzv = *((int *) nzv);
     }

   if (iflag64)
     {
      nndes = tmpnxv64 * tmpnyv64 * tmpnzv64;
      str_ncells = (tmpnxv64-1) * (tmpnyv64-1) * (tmpnzv64-1);
     }
   else
     {
      nndes = tmpnxv * tmpnyv * tmpnzv;
      str_ncells = (tmpnxv-1) * (tmpnyv-1) * (tmpnzv-1);
     }
   if (rflag64)
     {
      tempx64 = (double *) malloc(sizeof(double) * nndes);
      tempy64 = (double *) malloc(sizeof(double) * nndes);
      tempz64 = (double *) malloc(sizeof(double) * nndes);
     }
   else
     {
      tempx = (float*)malloc(sizeof(float)*nndes);
      tempy = (float*)malloc(sizeof(float)*nndes);
      tempz = (float*)malloc(sizeof(float)*nndes);
     }

   if (rflag64)
      for (count = 0; count < nndes; count++)
        {
 	 tempx64[count] = *((double *) x + count);
	 tempy64[count] = *((double *) y + count);
	 tempz64[count] = *((double *) z + count);
        }
   else
      for (count = 0; count < nndes; count++)
        {
         tempx[count] = *((float *) x + count);
         tempy[count] = *((float *) y + count);
         tempz[count] = *((float *) z + count);
        }

   /*  Write node x,y,z's for logicaly structured grids.  */
   n_nodes = nndes;
   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(&alt,INT32,1,fp);
     }
   else
      fprintf(fp,"nodes %d",alt);

   if (iflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tmpnxv64, INT64, 1, fp);
         fwrite(&tmpnyv64, INT64, 1, fp);
         fwrite(&tmpnzv64, INT64, 1, fp);
        }
      else
         fprintf(fp," %ld %ld %ld\n",tmpnxv64,tmpnyv64,tmpnzv64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tmpnxv, INT32, 1, fp);
         fwrite(&tmpnyv, INT32, 1, fp);
         fwrite(&tmpnzv, INT32, 1, fp);
        }
      else
         fprintf(fp," %d %d %d\n",tmpnxv,tmpnyv,tmpnzv);
     }
   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         tempxyz64 = (double *) malloc(sizeof(double) * nndes);
         j = 0;
         for (i = 0; i < nndes; i++)
           {
            tempxyz64[j] = tempx64[i];  j++;
            tempxyz64[j] = tempy64[i];  j++;
            tempxyz64[j] = tempz64[i];  j++;
           }
         fwrite(tempxyz64, FLOAT64, 3*nndes, fp);
         free(tempxyz64);
        }
      else
        {
         for (i = 0; i < nndes; i++)
           {
            fprintf(fp,"%lg %lg %lg \n",tempx64[i], tempy64[i], tempz64[i]);
           }
        }
      free(tempx64), free(tempy64), free(tempz64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         tempxyz = (float *) malloc(sizeof(float) * nndes);
         j = 0;
         for (i = 0; i < nndes; i++)
           {
            tempxyz[j] = tempx[i];  j++;
            tempxyz[j] = tempy[i];  j++;
            tempxyz[j] = tempz[i];  j++;
           }
         fwrite(tempxyz, FLOAT32, 3*nndes, fp);
         free(tempxyz);
        }
      else
        {
         for (i = 0; i < nndes; i++)
          {
           fprintf(fp,"%lg %lg %lg \n",tempx[i],tempy[i],tempz[i]);
          }
        }
      free(tempx), free(tempy), free(tempz);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_cells_amr(void *numcells, void *numtop, void *daughters)
{
  int i, tmpnumcells, tmpnumtop, *tmpdaughters=NULL;
  long ilong, tmpnumcells64, tmpnumtop64, *tmpdaughters64=NULL;

   if (iflag64)
     {
      tmpnumcells64 = *((long *) numcells);
      tmpnumtop64 = *((long *) numtop);
      tmpdaughters64 = (long *) malloc(sizeof(long) * tmpnumcells64);
      for (ilong = 0; ilong < tmpnumcells64; ilong++)
	tmpdaughters64[ilong] = *((long *) daughters + ilong);
      n_cells = tmpnumtop64;
     }
   else
     {
      tmpnumcells = *((int *) numcells);
      tmpnumtop = *((int *) numtop);
      tmpdaughters = (int *) malloc(sizeof(int) * tmpnumcells);
      for (i = 0; i < tmpnumcells; i++)
        tmpdaughters[i] = *((int *) daughters + i);
      n_cells = tmpnumtop;
     }

  strcpy(tmpname, "cells   ");
  if (filetype == IEEE_F)
     fwrite(tmpname, sizeof(char), 8, fp);
  else
     fprintf(fp,"cells ");
  if (iflag64)
    {
     if (filetype == IEEE_F)
       {
        fwrite(&tmpnumcells64, INT64, 1, fp);
        fwrite(&tmpnumtop64, INT64, 1, fp);
        fwrite(tmpdaughters64, INT64, tmpnumcells64, fp);
       }
     else
       {
        fprintf(fp,"%ld %ld\n", tmpnumcells64, tmpnumtop64);
        for (i = 0; i < tmpnumcells64; i++)
          {
           fprintf(fp,"%ld \n",tmpdaughters64[i]);
          }
       }
     free(tmpdaughters64);
    }
  else
    {
     if (filetype == IEEE_F)
       {
        fwrite(&tmpnumcells, INT32, 1, fp);
        fwrite(&tmpnumtop, INT32, 1, fp);
        fwrite(tmpdaughters, INT32, tmpnumcells, fp);
       }
     else
       {
        fprintf(fp,"%d %d\n", tmpnumcells, tmpnumtop);
        for (i = 0; i < tmpnumcells; i++)
          {
           fprintf(fp,"%d \n",tmpdaughters[i]);
          }
       }
     free(tmpdaughters);
    }
}

/* --------------------------------------------------------- */

void gmvwrite_cells_fromfile(char *filename, long nclls)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("cells   fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"cells   fromfile \"%s\"",filename);
   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);

   n_cells = nclls;
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_cell_header(void *ncells)
{
  int tmpncells;
  long tmpncells64;

  strcpy(tmpname, "cells   ");

  if (iflag64)
     tmpncells64 = *((long *) ncells);
  else
     tmpncells = *((int *) ncells);

  if (filetype == IEEE_F)
     fwrite(tmpname,sizeof(char),8,fp);
  if (iflag64)
    {
     if (filetype == IEEE_F)
        fwrite(&tmpncells64, INT64, 1, fp);
     else
        fprintf(fp,"cells %ld\n",tmpncells64);
     n_cells = tmpncells64;
    }
  else
    {
     if (filetype == IEEE_F)
        fwrite(&tmpncells,INT32,1,fp);
     else
        fprintf(fp,"cells %d\n",tmpncells);
     n_cells = tmpncells;
    }

   if (structflag) n_cells = str_ncells;
}

/* --------------------------------------------------------- */

void gmvwrite_cell_type(char cell_type[], int nverts, void *nodes)
 {
   int tmpnverts, count;
   int *tempnodes=NULL, *nodesptr=NULL;
   long *tempnodes64=NULL, *nodesptr64=NULL;

   tmpnverts = nverts;

   if (iflag64)
     {
      tempnodes64 = (long *) malloc(sizeof(long) * nverts);
      nodesptr64 = (long *) nodes;
     }
   else
     {
      tempnodes = (int *) malloc(sizeof(int)*nverts);
      nodesptr = (int *) nodes;
     }

   if (iflag64)
      for (count = 0; count < nverts; count++)
	  tempnodes64[count] = nodesptr64[count];
   else
      for (count = 0;count < nverts;count++)
	  tempnodes[count] = nodesptr[count];

   if (filetype == IEEE_F)
     {
      fwrite(cell_type,sizeof(char),8,fp);
      fwrite(&tmpnverts,INT32,1,fp);
     }
   else
      fprintf(fp,"%s %d",cell_type,tmpnverts);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempnodes64,INT64,nverts,fp);
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp," %ld",tempnodes64[i]);
         fprintf(fp,"\n");
        }
      free(tempnodes64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempnodes, INT32, nverts, fp);
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp," %d",tempnodes[i]);
         fprintf(fp,"\n");
        }
      free(tempnodes);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_general_cell_type(char cell_type[], int nverts[], int nfaces,
			        void *nodeids)
 {
   int sumverts=0, count, *tempnverts, tmpnfaces, *tempnodeids=NULL;
   long *tempnodeids64=NULL, *nodesptr64=NULL;
   int *nodesptr=NULL;

   tmpnfaces = nfaces;

   tempnverts = (int*)malloc(sizeof(int)*nfaces);

   for (count = 0;count < nfaces;count++)
       tempnverts[count] = nverts[count];

   for(count = 0; count < nfaces; count++)
     sumverts += nverts[count];

   if (iflag64)
     {
      tempnodeids64 = (long *) malloc(sizeof(long) * sumverts);
      nodesptr64 = (long *) nodeids;
     }
   else
     {
      tempnodeids = (int *) malloc(sizeof(int)*sumverts);
      nodesptr = (int *) nodeids;
     }

   if (iflag64)
      for (count = 0;count < sumverts;count++)
          tempnodeids64[count] = nodesptr64[count];
   else
      for (count = 0; count < sumverts; count++)
	  tempnodeids[count] = nodesptr[count];

   if (filetype == IEEE_F)
     {
      fwrite(cell_type,sizeof(char),8,fp);
      fwrite(&tmpnfaces,INT32,1,fp);
      fwrite(tempnverts,INT32,nfaces,fp);
     }
   else
     {
      fprintf(fp,"%s %d\n",cell_type,tmpnfaces);
      for (i = 0; i < nfaces; i++)
         fprintf(fp," %d",tempnverts[i]);
      fprintf(fp,"\n");
     }
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempnodeids64,INT64,sumverts,fp);
      else
        {
         for (i = 0; i < sumverts; i++)
            fprintf(fp," %ld",tempnodeids64[i]);
         fprintf(fp,"\n");
        }
      free(tempnodeids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempnodeids, INT32, sumverts, fp);
      else
        {
         for (i = 0; i < sumverts; i++)
            fprintf(fp," %d",tempnodeids[i]);
         fprintf(fp,"\n");
        }
      free(tempnodeids);
     }
   free(tempnverts);
 }

/* --------------------------------------------------------- */

void gmvwrite_faces_fromfile(char *filename, long nclls)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("faces   fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"faces   fromfile \"%s\"",filename);
   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);

   n_cells = nclls;
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_face_header(void *nfaces, void *ncells)
{
  int tmpnfaces, tmpncells;
  long tmpnfaces64, tmpncells64;

  strcpy(tmpname, "faces   ");

  if (iflag64)
    {
     tmpnfaces64 = *((long *) nfaces);
     tmpncells64 = *((long *) ncells);
    }
  else
    {
     tmpnfaces = *((int *) nfaces);
     tmpncells = *((int *) ncells);
    }

   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   if (iflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tmpnfaces64, INT64, 1, fp);
         fwrite(&tmpncells64, INT64, 1, fp);
        }
      else
         fprintf(fp,"faces %ld %ld\n",tmpnfaces64,tmpncells64);
      n_faces = tmpnfaces64;
      n_cells = tmpncells64;
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tmpnfaces,INT32,1,fp);
         fwrite(&tmpncells,INT32,1,fp);
        }
      else
        fprintf(fp,"faces %d %d\n",tmpnfaces,tmpncells);
      n_faces = tmpnfaces;
      n_cells = tmpncells;
     }
}

/* --------------------------------------------------------- */

void gmvwrite_face_data(int nverts, void *nodeids, void *cellid1,
                   void *cellid2)
{
   int  tempnverts, *tempnodeids=NULL, tmpcellid1, tmpcellid2;
   long *tempnodeids64=NULL, *nodesptr64=NULL, tmpcellid164, tmpcellid264;
   int *nodesptr=NULL, count;

   tempnverts = nverts;

   /*  Save face nodes.  */
   if (iflag64)
     {
      tempnodeids64 = (long *) malloc(sizeof(long) * nverts);
      nodesptr64 = (long *) nodeids;
     }
   else
     {
      tempnodeids = (int *) malloc(sizeof(int) * nverts);
      nodesptr = (int *) nodeids;
     }

   if (iflag64)
      for (count = 0;count < nverts;count++)
         tempnodeids64[count] = nodesptr64[count];
   else
       for (count = 0; count < nverts; count++)
	  tempnodeids[count] = nodesptr[count];

   /*  Save face cells.  */
   if (iflag64)
     {
      tmpcellid164 = *((long *) cellid1);
      tmpcellid264 = *((long *) cellid2);
     }
   else
     {
      tmpcellid1 = *((int *) cellid1);
      tmpcellid2 = *((int *) cellid2);
     }

   if (filetype == IEEE_F)
      fwrite(&tempnverts,INT32,1,fp);
   else
      fprintf(fp,"%d",tempnverts);
   if (iflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempnodeids64,INT64,nverts,fp);
         fwrite(&tmpcellid164,INT64,1,fp);
         fwrite(&tmpcellid264,INT64,1,fp);
        }
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp," %ld",tempnodeids64[i]);
         fprintf(fp," %ld %ld\n",tmpcellid164,tmpcellid264);
        }
      free(tempnodeids64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempnodeids, INT32, nverts, fp);
         fwrite(&tmpcellid1, INT32, 1, fp);
         fwrite(&tmpcellid2, INT32, 1, fp);
        }
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp," %d",tempnodeids[i]);
         fprintf(fp," %d %d\n",tmpcellid1,tmpcellid2);
        }
      free(tempnodeids);
     }

 }

/* --------------------------------------------------------- */

void gmvwrite_vfaces_fromfile(char *filename, long nclls)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("vfaces  fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"vfaces  fromfile \"%s\"",filename);
   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);

   n_cells = nclls;
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_vface_header(void *nfaces)
{
  int tmpnfaces;
  long tmpnfaces64;

  strcpy(tmpname, "vfaces  ");

  if (iflag64)
     tmpnfaces64 = *((long *) nfaces);
  else
     tmpnfaces = *((int *) nfaces);

   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnfaces64, INT64, 1, fp);
      else
         fprintf(fp,"vfaces %ld\n",tmpnfaces64);
      n_faces = tmpnfaces64;
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnfaces,INT32,1,fp);
      else
         fprintf(fp,"vfaces %d\n",tmpnfaces);
      n_faces = tmpnfaces;
     }
}

/* --------------------------------------------------------- */

void gmvwrite_vface_data(int nverts, int facepe, void *oppface,
                         int oppfacepe, void *cellid, void *nodeids)
{
   int  tempnverts, *tempnodeids=NULL, tempcellid, tempfacepe, tempoppface,
        tempoppfacepe;
   long *tempnodeids64=NULL, *nodesptr64=NULL, tempcellid64, tempoppface64;
   int *nodesptr=NULL, count;

   tempnverts = nverts;
   tempfacepe = facepe;
   tempoppfacepe = oppfacepe;

   /*  Save face nodes.  */
   if (iflag64)
     {
      tempoppface64 = *((long *) oppface);
      tempcellid = *((long *) cellid);
      tempnodeids64 = (long *) malloc(sizeof(long) * nverts);
      nodesptr64 = (long *) nodeids;
     }
   else
     {
      tempoppface = *((int *) oppface);
      tempcellid = *((int *) cellid);
      tempnodeids = (int *) malloc(sizeof(int) * nverts);
      nodesptr = (int *) nodeids;
     }

   if (iflag64)
      for (count = 0;count < nverts;count++)
         tempnodeids64[count] = nodesptr64[count];
   else
      for (count = 0; count < nverts; count++)
	 tempnodeids[count] = nodesptr[count];

   /*  Save face cells.  */
   if (iflag64)
     {
      tempcellid64 = *((long *) cellid);
     }
   else
     {
      tempcellid = *((int *) cellid);
     }

   if (filetype == IEEE_F)
     {
      fwrite(&tempnverts,INT32,1,fp);
      fwrite(&tempfacepe,INT32,1,fp);
     }
   else
      fprintf(fp," %d %d",tempnverts,tempfacepe);
   if (iflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tempoppface64,INT64,1,fp);
         fwrite(&tempoppfacepe,INT32,1,fp);
         fwrite(&tempcellid64,INT64,1,fp);
         fwrite(tempnodeids64,INT64,nverts,fp);
        }
      else
        {
         fprintf(fp," %ld %d %ld ",tempoppface64,tempoppfacepe,tempcellid64);
         for (i = 0; i < nverts; i++)
            fprintf(fp," %ld",tempnodeids64[i]);
         fprintf(fp,"\n");
        }
      free(tempnodeids64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(&tempoppface,INT32,1,fp);
         fwrite(&tempoppfacepe,INT32,1,fp);
         fwrite(&tempcellid, INT32, 1, fp);
         fwrite(tempnodeids, INT32, nverts, fp);
        }
      else
        {
         fprintf(fp," %d %d %d ",tempoppface,tempoppfacepe,tempcellid);
         for (i = 0; i < nverts; i++)
            fprintf(fp," %d",tempnodeids[i]);
         fprintf(fp,"\n");
        }
      free(tempnodeids);
     }
 }


/* --------------------------------------------------------- */

void gmvwrite_xfaces_fromfile(char *filename, long nfces, long nclls)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("xfaces  fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"xfaces  fromfile \"%s\"",filename);
   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);

   n_faces = nfces;
   n_cells = nclls;
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_xface_header(void *nfaces)
{
  int tmpnfaces;
  long tmpnfaces64;

   strcpy(tmpname, "xfaces  ");
   if (iflag64)
      tmpnfaces64 = *((long *) nfaces);
   else
      tmpnfaces = *((int *) nfaces);
   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnfaces64, INT64, 1, fp);
      else
         fprintf(fp,"xfaces %ld\n",tmpnfaces64);
      n_faces = tmpnfaces64;
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnfaces,INT32,1,fp);
      else
         fprintf(fp,"xfaces %d\n",tmpnfaces);
      n_faces = tmpnfaces;
     }
}

/* --------------------------------------------------------- */

void gmvwrite_xface_data(long totverts, void *nverts, void *nodeids,
                         void *cellid, void *oppface, void *facepe,
                         void *oppfacepe)
{
  int  *tempface=NULL, *tempnodeids=NULL, *tempptr=NULL;
  long  *tempface64=NULL, *tempnodeids64=NULL, *tempptr64=NULL, nc, tempnverts;
  int count;

   /*  Allocate output arrays.  */
   tempnverts = totverts;
   if (iflag64)
     {
      tempface64 = (long *)malloc(sizeof(long) * n_faces);
      tempnodeids64 = (long *)malloc(sizeof(long) * tempnverts);
     }
   else
     {
      tempface = (int *)malloc(sizeof(int) * n_faces);
      tempnodeids = (int *)malloc(sizeof(int) * tempnverts);
     }

   /*  Save and write numverts array and the face nodes array.  */
   if (iflag64)
     {
      tempptr64 = nverts;
      for (count = 0; count < n_faces; count++)
         tempface64[count] = tempptr64[count];
      tempptr64 = nodeids;
      for (count = 0; count < totverts; count++)
         tempnodeids64[count] = tempptr64[count];
     }
   else
     {
      tempptr = nverts;
      for (count = 0; count < n_faces; count++)
         tempface[count] = tempptr[count];
      tempptr = nodeids;
      for (count = 0; count < totverts; count++)
	 tempnodeids[count] = tempptr[count];
     }
   if (iflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempface64,INT64,n_faces,fp);
         fwrite(tempnodeids64,INT64,totverts,fp);
        }
      else
        {
         write_ascii_long(n_faces, tempface64);
         write_ascii_long(totverts, tempnodeids64);
        }
      free(tempnodeids64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempface,INT32,n_faces,fp);
         fwrite(tempnodeids,INT32,totverts,fp);
        }
      else
        {
         write_ascii_int(n_faces, tempface);
         write_ascii_int(totverts, tempnodeids);
        }
      free(tempnodeids);
     }

   /*  Save and write face cells, also calculate the  */
   /*  number of cells from the largest cell number.  */
   nc = 0;
   if (iflag64)
     {
      tempptr64 = cellid;
      for (count = 0; count < n_faces; count++)
        {
         tempface64[count] = tempptr64[count];
         if (tempface64[count] > nc) nc = tempface64[count];
        }
     }
   else
     {
      tempptr = cellid;
      for (count = 0; count < n_faces; count++)
        {
         tempface[count] = tempptr[count];
         if (tempface[count] > nc) nc = tempface[count];
        }
     }
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempface64,INT64,n_faces,fp);
      else
         write_ascii_long(n_faces, tempface64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempface,INT32,n_faces,fp);
      else
         write_ascii_int(n_faces, tempface);
     }
   n_cells = nc;

   /*  Save and write opposite face data.  */
   if (iflag64)
     {
      tempptr64 = oppface;
      for (count = 0; count < n_faces; count++)
        {
         tempface64[count] = tempptr64[count];
        }
     }
   else
     {
      tempptr = oppface;
      for (count = 0; count < n_faces; count++)
        {
         tempface[count] = tempptr[count];
        }
     }
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempface64,INT64,n_faces,fp);
      else
         write_ascii_long(n_faces, tempface64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempface,INT32,n_faces,fp);
      else
         write_ascii_int(n_faces, tempface);
     }

   /*  Save and write face pe data.  */
   if (iflag64)
     {
      tempptr64 = facepe;
      for (count = 0; count < n_faces; count++)
        {
         tempface64[count] = tempptr64[count];
        }
     }
   else
     {
      tempptr = facepe;
      for (count = 0; count < n_faces; count++)
        {
         tempface[count] = tempptr[count];
        }
     }
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempface64,INT64,n_faces,fp);
      else
         write_ascii_long(n_faces, tempface64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempface,INT32,n_faces,fp);
      else
        write_ascii_int(n_faces, tempface);
     }

   /*  Save and write opposite face pe data.  */
   if (iflag64)
     {
      tempptr64 = oppfacepe;
      for (count = 0; count < n_faces; count++)
        {
         tempface64[count] = tempptr64[count];
        }
     }
   else
     {
      tempptr = oppfacepe;
      for (count = 0; count < n_faces; count++)
        {
         tempface[count] = tempptr[count];
        }
     }
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempface64,INT64,n_faces,fp);
      else
         write_ascii_long(n_faces, tempface64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempface,INT32,n_faces,fp);
      else
         write_ascii_int(n_faces, tempface);
     }

   if (iflag64)
      free(tempface64);
   else
      free(tempface);
 }

/* --------------------------------------------------------- */

void gmvwrite_material_fromfile(char *filename)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("materialfromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"materialfromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s",tmpbuf);

   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_material_header(int nmats, int data_type)
 {
  int tmpnmats, tmpdata_type;

  strcpy(tmpname, "material");
  tmpnmats = nmats;
  tmpdata_type = data_type;

   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(&tmpnmats,INT32,1,fp);
      fwrite(&tmpdata_type,INT32,1,fp);
     }
   else
      fprintf(fp,"%s %d %d\n",tmpname,tmpnmats,tmpdata_type);
 }

/* --------------------------------------------------------- */

void gmvwrite_material_name(char matname[])
 {
   if (filetype == IEEE_F)
      fwrite(matname,sizeof(char),charsize_out,fp);
   else
      fprintf(fp,"%s\n",matname);
 }

/* --------------------------------------------------------- */

void gmvwrite_material_ids(int matids[], int data_type)
 {
  int *tempmatids=NULL;
  long npts=0, count;

   if (data_type == 0)
     {
      npts = n_cells;
     }
   else if (data_type == 1)
     {
      npts = n_nodes;
     }
   tempmatids = (int*)malloc(sizeof(int)*npts);

   for (count = 0; count <npts; count++)
     {
      tempmatids[count] = matids[count];
     }

   if (filetype == IEEE_F)
     fwrite(tempmatids,INT32,npts,fp);
   else
     {
      write_ascii_int(npts, tempmatids);
    }
   free(tempmatids);
 }

/* --------------------------------------------------------- */

void gmvwrite_velocity_data(int data_type, void *u, void *v, void *w)
 {
  int tmpdata_type;
  long npts=0, i;
  float *tempu=NULL, *tempv=NULL, *tempw=NULL;
  double *tempu64=NULL, *tempv64=NULL, *tempw64=NULL;
  strcpy(tmpname, "velocity");

   tmpdata_type = data_type;

  if(data_type == 0)
   {
    npts = n_cells;
   }
  else if(data_type == 1)
   {
    npts = n_nodes;
   }
  else if(data_type == 2)
   {
    npts = n_faces;
   }

  if (rflag64)
    {
     tempu64 = (double*)malloc(sizeof(double)*npts);
     tempv64 = (double*)malloc(sizeof(double)*npts);
     tempw64 = (double*)malloc(sizeof(double)*npts);
    }
  else
    {
     tempu = (float*)malloc(sizeof(float)*npts);
     tempv = (float*)malloc(sizeof(float)*npts);
     tempw = (float*)malloc(sizeof(float)*npts);
    }

  if (rflag64)
    {
     for (i = 0; i < npts; i++)
       {
        tempu64[i] = *((double *) u + i);
        tempv64[i] = *((double *) v + i);
        tempw64[i] = *((double *) w + i);
       }
    }
  else
    {
     for (i = 0; i < npts; i++)
       {
        tempu[i] = *((float *) u + i);
        tempv[i] = *((float *) v + i);
        tempw[i] = *((float *) w + i);
       }
    }

   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(&tmpdata_type,INT32,1,fp);
     }
   else
      fprintf(fp,"velocity %d\n",tmpdata_type);

   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempu64,FLOAT64,npts,fp);
         fwrite(tempv64,FLOAT64,npts,fp);
         fwrite(tempw64,FLOAT64,npts,fp);
        }
      else
        {
         write_ascii_double(npts, tempu64);
         write_ascii_double(npts, tempv64);
         write_ascii_double(npts, tempw64);
        }
      free(tempu64); free(tempv64); free(tempw64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempu,FLOAT32,npts,fp);
         fwrite(tempv,FLOAT32,npts,fp);
         fwrite(tempw,FLOAT32,npts,fp);
        }
      else
        {
         write_ascii_float(npts, tempu);
         write_ascii_float(npts, tempv);
         write_ascii_float(npts, tempw);
        }
      free(tempu); free(tempv); free(tempw);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_variable_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "variable");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"variable \n");
}

/* --------------------------------------------------------- */

void gmvwrite_variable_name_data(int data_type, char varname[], void *vids)
 {
  int tmpdata_type;
  long count, npts=0;
  float *tempvids=NULL;
  double *tempvids64=NULL;

  tmpdata_type = data_type;

  if(data_type == 0)
   {
    npts = n_cells;
   }
  else if(data_type == 1)
   {
    npts = n_nodes;
   }
  else if(data_type == 2)
   {
    npts = n_faces;
   }

  if (rflag64)
   {
    tempvids64 = (double*)malloc(sizeof(double)*npts);
   }
  else
   {
    tempvids = (float*)malloc(sizeof(float)*npts);
   }

  if (rflag64)
    {
     for (count = 0; count < npts; count++)
       {
        tempvids64[count] = *((double *) vids + count);
       }
    }
  else
    {
     for (count = 0; count < npts; count++)
       {
        tempvids[count] = *((float *) vids + count);
       }
    }

   if (filetype == IEEE_F)
     {
      fwrite(varname, sizeof(char), charsize_out, fp);
      fwrite(&tmpdata_type, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d\n",varname,tmpdata_type);
   if (rflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempvids64, FLOAT64, npts, fp);
      else
        {
         write_ascii_double(npts, tempvids64);
        }
      free(tempvids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempvids, FLOAT32, npts, fp);
      else
        {
         write_ascii_float(npts, tempvids);
        }
      free(tempvids);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_variable_endvars(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endvars ");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"endvars \n");
}

/* --------------------------------------------------------- */

void gmvwrite_flag_fromfile(char *filename)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("flags   fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"flags   fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);

   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_flag_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "flags   ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"flags \n");
}

/* --------------------------------------------------------- */

void gmvwrite_flag_name(char flagname[], int numtypes, int data_type)
{
  int tmpdata_type, tmpnumtypes;

   tmpdata_type = data_type;
   tmpnumtypes = numtypes;

   if (filetype == IEEE_F)
     {
      fwrite(flagname, sizeof(char), charsize_out, fp);
      fwrite(&tmpnumtypes, INT32, 1, fp);
      fwrite(&tmpdata_type, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d %d\n",flagname,tmpnumtypes,tmpdata_type);
}

/* --------------------------------------------------------- */

void gmvwrite_flag_subname(char subname[])
{
   if (filetype == IEEE_F)
      fwrite(subname, sizeof(char), charsize_out, fp);
   else
      fprintf(fp,"%s \n",subname);
}

/* --------------------------------------------------------- */

void gmvwrite_flag_data(int data_type, int flag_data[])
{
  int *tmpdata;
  long numthings, i;

  if (data_type)
    numthings = n_nodes;
  else
    numthings = n_cells;
  tmpdata = (int *) malloc(sizeof(int) * numthings);
  for (i = 0; i < numthings; i++)
     tmpdata[i] = flag_data[i];

   if (filetype == IEEE_F)
      fwrite(tmpdata, INT32, numthings, fp);
   else
     {
      write_ascii_int(numthings, tmpdata);
     }
   free(tmpdata);
}

/* --------------------------------------------------------- */

void gmvwrite_flag_endflag(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endflag ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"endflag \n");
}

/* --------------------------------------------------------- */

void gmvwrite_polygons_fromfile(char *filename)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("polygonsfromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"polygonsfromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);

   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_polygons_header(void)
{
   strcpy(tmpname, "polygons");
   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   else
   fprintf(fp,"polygons \n");
}

/* --------------------------------------------------------- */

void gmvwrite_polygons_data(int nverts, int matnum, void *x, void *y, void *z)
 {
   float tx[3000], ty[3000], tz[3000];
   double tx64[3000], ty64[3000], tz64[3000];
   int tmpnverts, tmpmatnum, count;

   tmpnverts = nverts;
   tmpmatnum = matnum;
   if (rflag64)
     {
      for (count = 0; count <  nverts; count++)
        {
         tx64[count] = *((double *) x + count);
         ty64[count] = *((double *) y + count);
         tz64[count] = *((double *) z + count);
        }
     }
   else
     {
      for (count = 0; count <  nverts; count++)
        {
         tx[count] = *((float *) x + count);
         ty[count] = *((float *) y + count);
         tz[count] = *((float *) z + count);
        }
     }

   if (filetype == IEEE_F)
     {
      fwrite(&tmpmatnum, INT32, 1, fp);
      fwrite(&tmpnverts, INT32, 1, fp);
     }
   else
      fprintf(fp,"%d %d\n",tmpmatnum,tmpnverts);
   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tx64,FLOAT64,nverts,fp);
         fwrite(ty64,FLOAT64,nverts,fp);
         fwrite(tz64,FLOAT64,nverts,fp);
        }
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp,"%lg ",tx64[i]);
         fprintf(fp,"\n");
         for (i = 0; i < nverts; i++)
            fprintf(fp,"%lg ",ty64[i]);
         fprintf(fp,"\n");
         for (i = 0; i < nverts; i++)
            fprintf(fp,"%lg ",tz64[i]);
         fprintf(fp,"\n");
        }
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tx,FLOAT32,nverts,fp);
         fwrite(ty,FLOAT32,nverts,fp);
         fwrite(tz,FLOAT32,nverts,fp);
        }
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp,"%g ",tx[i]);
         fprintf(fp,"\n");
         for (i = 0; i < nverts; i++)
            fprintf(fp,"%g ",ty[i]);
         fprintf(fp,"\n");
         for (i = 0; i < nverts; i++)
            fprintf(fp,"%g ",tz[i]);
         fprintf(fp,"\n");
        }
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_polygons_endpoly(void)
{
   strcpy(tmpname, "endpoly ");
   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   else
      fprintf(fp,"endpoly \n");
}

/* --------------------------------------------------------- */

void gmvwrite_tracers_header(int ntracers, void *x, void *y, void *z)
{
  int tmpntracers, count;
  float *tempx=NULL, *tempy=NULL, *tempz=NULL;
  double *tempx64=NULL, *tempy64=NULL, *tempz64=NULL;

  if (rflag64)
    {
     tempx64 = (double*)malloc(sizeof(double)*ntracers);
     tempy64 = (double*)malloc(sizeof(double)*ntracers);
     tempz64 = (double*)malloc(sizeof(double)*ntracers);
    }
  else
    {
     tempx = (float*)malloc(sizeof(float)*ntracers);
     tempy = (float*)malloc(sizeof(float)*ntracers);
     tempz = (float*)malloc(sizeof(float)*ntracers);
    }

  strcpy(tmpname, "tracers ");

  tmpntracers = ntracers;
  if (rflag64)
    {
     for (count = 0; count< ntracers; count++)
       {
        tempx64[count] = *((double *) x + count);
        tempy64[count] = *((double *) y + count);
        tempz64[count] = *((double *) z + count);
       }
    }
  else
    {
     for (count = 0; count< ntracers; count++)
       {
        tempx[count] = *((float *) x + count);
        tempy[count] = *((float *) y + count);
        tempz[count] = *((float *) z + count);
       }
    }

   if (filetype == IEEE_F)
     {
      fwrite(tmpname, sizeof(char), 8, fp);
      fwrite(&tmpntracers, INT32, 1, fp);
     }
   else
      fprintf(fp,"tracers %d\n",tmpntracers);
   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx64, FLOAT64, ntracers, fp);
         fwrite(tempy64, FLOAT64, ntracers, fp);
         fwrite(tempz64, FLOAT64, ntracers, fp);
        }
      else
        {
         write_ascii_double((long)ntracers, tempx64);
         write_ascii_double((long)ntracers, tempy64);
         write_ascii_double((long)ntracers, tempz64);
        }
      free(tempx64); free(tempy64); free(tempz64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempx, FLOAT32, ntracers, fp);
         fwrite(tempy, FLOAT32, ntracers, fp);
         fwrite(tempz, FLOAT32, ntracers, fp);
        }
      else
        {
         write_ascii_float((long)ntracers, tempx);
         write_ascii_float((long)ntracers, tempy);
         write_ascii_float((long)ntracers, tempz);
        }
      free(tempx); free(tempy); free(tempz);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_tracers_name_data(int ntracers, char tracername[], void *data)
 {
  float *tmpdata=NULL;
  double *tmpdata64=NULL;
  int count;

  if (rflag64)
     tmpdata64 = (double*)malloc(sizeof(double) * ntracers);
  else
     tmpdata = (float*)malloc(sizeof(float) * ntracers);

  if (rflag64)
    {
     for (count = 0; count < ntracers; count++)
       {
        tmpdata64[count] = *((double *) data + count);
       }
    }
  else
    {
     for (count = 0; count < ntracers; count++)
      {
       tmpdata[count] = *((float *) data + count);
      }
    }

   if (filetype == IEEE_F)
      fwrite(tracername,sizeof(char),charsize_out,fp);
   else
      fprintf(fp,"%s\n",tracername);
   if (rflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata64,FLOAT64,ntracers,fp);
      else
         write_ascii_double((long)ntracers, tmpdata64);
      free(tmpdata64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata,FLOAT32,ntracers,fp);
      else
         write_ascii_float((long)ntracers, tmpdata);
      free(tmpdata);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_tracers_endtrace(void)
{
   strcpy(tmpname, "endtrace");
   if (filetype == IEEE_F)
      fwrite(tmpname, sizeof(char), 8, fp);
   else
      fprintf(fp,"endtrace \n");
}

/* --------------------------------------------------------- */

void gmvwrite_probtime(double ptime)
{
  float tmpptime;
  double tmpptime64;
  strcpy(tmpname, "probtime");

  if (rflag64)
     tmpptime64 = ptime;
  else
     tmpptime = ptime;

   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      if (rflag64)
         fwrite(&tmpptime64,FLOAT64,1,fp);
      else
         fwrite(&tmpptime,FLOAT32,1,fp);
     }
   else
     {
      if (rflag64)
         fprintf(fp,"probtime %lg\n",tmpptime64);
      else
         fprintf(fp,"probtime %g\n",tmpptime);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_cycleno(int cyclenum)
{
  int tmpcyclenum;
  strcpy(tmpname, "cycleno ");
  tmpcyclenum = cyclenum;

   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(&tmpcyclenum,INT32,1,fp);
     }
   else
      fprintf(fp,"cycleno %d\n",tmpcyclenum);
}

/* --------------------------------------------------------- */

void gmvwrite_nodeids(void *nodeids)
{
  int i, *tmpnodeids=NULL;
  long ilong, *tmpnodeids64=NULL;

  if (iflag64)
    {
     tmpnodeids64 = (long *) malloc(sizeof(long) * n_nodes);
     for (ilong = 0; ilong < n_nodes; ilong++)
	tmpnodeids64[ilong] = *((long *) nodeids + ilong);
    }
 else
    {
     tmpnodeids = (int *) malloc(sizeof(int) * n_nodes);
     for (i = 0; i < n_nodes; i++)
        tmpnodeids[i] = *((int *) nodeids + i);
    }

   strcpy(tmpname, "nodeids ");
   if (filetype == IEEE_F)
      fwrite(tmpname, sizeof(char), 8, fp);
   else
      fprintf(fp,"nodeids \n");
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpnodeids64, INT64, n_nodes, fp);
      else
         write_ascii_long(n_nodes, tmpnodeids64);
      free(tmpnodeids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpnodeids, INT32, n_nodes, fp);
      else
         write_ascii_int(n_nodes, tmpnodeids);
      free(tmpnodeids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_nodeids_fromfile(char *filename)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("nodeids fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"nodeids fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_cellids(void *cellids)
{
  int i, *tmpcellids=NULL;
  long ilong, *tmpcellids64=NULL;

  if (iflag64)
    {
     tmpcellids64 = (long *) malloc(sizeof(long) * n_cells);
     for (ilong = 0; ilong < n_cells; ilong++)
	tmpcellids64[ilong] = *((long *) cellids + ilong);
    }
 else
    {
     tmpcellids = (int *) malloc(sizeof(int) * n_cells);
     for (i = 0; i < n_cells; i++)
        tmpcellids[i] = *((int *) cellids + i);
    }

   strcpy(tmpname, "cellids ");
   if (filetype == IEEE_F)
      fwrite(tmpname, sizeof(char), 8, fp);
   else
      fprintf(fp,"cellids \n");
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpcellids64, INT64, n_cells, fp);
      else
         write_ascii_long(n_cells, tmpcellids64);
      free(tmpcellids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpcellids, INT32, n_cells, fp);
      else
         write_ascii_int(n_cells, tmpcellids);
      free(tmpcellids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_cellids_fromfile(char *filename)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("cellids fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"cellids fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_surface_header(void *nsurf)
{
  int tmpnsurf;
  long tmpnsurf64;

  strcpy(tmpname, "surface ");
  if (iflag64)
    {
     tmpnsurf64 = *((long *) nsurf);
    }
  else
    {
     tmpnsurf = *((int *) nsurf);
    }

   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnsurf64, INT64, 1, fp);
      else
         fprintf(fp,"surface %ld\n",tmpnsurf64);
      n_surf = tmpnsurf64;
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(&tmpnsurf,INT32,1,fp);
      else
         fprintf(fp,"surface %d\n",tmpnsurf);
      n_surf = tmpnsurf;
     }
}

/* --------------------------------------------------------- */

void gmvwrite_surface_data(int nverts, void *nodeids)
{
   int  tempnverts, *tempnodeids=NULL;
   long *tempnodeids64=NULL, *nodesptr64=NULL;
   int *nodesptr=NULL, count;

   tempnverts = nverts;

   /*  Save surface nodes.  */
   if (iflag64)
     {
      tempnodeids64 = (long *) malloc(sizeof(long) * nverts);
      nodesptr64 = (long *) nodeids;
     }
   else
     {
      tempnodeids = (int *) malloc(sizeof(int) * nverts);
      nodesptr = (int *) nodeids;
     }

   if (iflag64)
      for (count = 0;count < nverts;count++)
         tempnodeids64[count] = nodesptr64[count];
   else
      for (count = 0; count < nverts; count++)
	 tempnodeids[count] = nodesptr[count];

   if (filetype == IEEE_F)
      fwrite(&tempnverts,INT32,1,fp);
   else
      fprintf(fp,"%d",tempnverts);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempnodeids64,INT64,nverts,fp);
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp," %ld",tempnodeids64[i]);
         fprintf(fp,"\n");
        }
      free(tempnodeids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempnodeids, INT32, nverts, fp);
      else
        {
         for (i = 0; i < nverts; i++)
            fprintf(fp," %d",tempnodeids[i]);
         fprintf(fp,"\n");
        }
      free(tempnodeids);
     }

}

/* --------------------------------------------------------- */

void gmvwrite_surface_fromfile(char *filename, long nsrf)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("surface fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"surface fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);

   free(tmpbuf);

   n_surf = nsrf;
   return;
 }

/* --------------------------------------------------------- */

void gmvwrite_surfmats(int matids[])
{
  int *tempmatids, count;

  strcpy(tmpname, "surfmats");

  tempmatids = (int*)malloc(sizeof(int)*n_surf);

  for (count = 0; count <n_surf; count++)
    {
     tempmatids[count] = matids[count];
    }

   if (filetype == IEEE_F)
     {
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(tempmatids,INT32,n_surf,fp);
     }
   else
     {
      fprintf(fp,"surfmats \n");
      write_ascii_int(n_surf, tempmatids);
     }
   free(tempmatids);
}

/* --------------------------------------------------------- */

void gmvwrite_surfvel( void *u, void *v, void *w)
{
  long i;
  float *tempu=NULL, *tempv=NULL, *tempw=NULL;
  double *tempu64=NULL, *tempv64=NULL, *tempw64=NULL;

   strcpy(tmpname, "surfvel ");

   if (rflag64)
     {
      tempu64 = (double*)malloc(sizeof(double)*n_surf);
      tempv64 = (double*)malloc(sizeof(double)*n_surf);
      tempw64 = (double*)malloc(sizeof(double)*n_surf);
     }
   else
     {
      tempu = (float*)malloc(sizeof(float)*n_surf);
      tempv = (float*)malloc(sizeof(float)*n_surf);
      tempw = (float*)malloc(sizeof(float)*n_surf);
     }

  if (rflag64)
    {
     for (i = 0; i < n_surf; i++)
       {
        tempu64[i] = *((double *) u + i);
        tempv64[i] = *((double *) v + i);
        tempw64[i] = *((double *) w + i);
       }
    }
  else
    {
     for (i = 0; i < n_surf; i++)
       {
        tempu[i] = *((float *) u + i);
        tempv[i] = *((float *) v + i);
        tempw[i] = *((float *) w + i);
       }
    }

   if (filetype == IEEE_F)
      fwrite(tmpname,sizeof(char),8,fp);
   else
      fprintf(fp,"surfvel \n");

   if (rflag64)
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempu64,FLOAT64,n_surf,fp);
         fwrite(tempv64,FLOAT64,n_surf,fp);
         fwrite(tempw64,FLOAT64,n_surf,fp);
        }
      else
        {
         write_ascii_double(n_surf, tempu64);
         write_ascii_double(n_surf, tempv64);
         write_ascii_double(n_surf, tempw64);
        }
      free(tempu64); free(tempv64); free(tempw64);
     }
   else
     {
      if (filetype == IEEE_F)
        {
         fwrite(tempu,FLOAT32,n_surf,fp);
         fwrite(tempv,FLOAT32,n_surf,fp);
         fwrite(tempw,FLOAT32,n_surf,fp);
        }
      else
        {
         write_ascii_float(n_surf, tempu);
         write_ascii_float(n_surf, tempv);
         write_ascii_float(n_surf, tempw);
        }
      free(tempu); free(tempv); free(tempw);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_surfvars_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "surfvars");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"surfvars \n");
}

/* --------------------------------------------------------- */

void gmvwrite_surfvars_name_data(char varname[], void *vids)
{
  long count;
  float *tempvids=NULL;
  double *tempvids64=NULL;

  if (rflag64)
   {
    tempvids64 = (double*)malloc(sizeof(double)*n_surf);
   }
  else
   {
    tempvids = (float*)malloc(sizeof(float)*n_surf);
   }

  if (rflag64)
    {
     for (count = 0; count < n_surf; count++)
       {
        tempvids64[count] = *((double *) vids + count);
       }
    }
  else
    {
     for (count = 0; count < n_surf; count++)
       {
        tempvids[count] = *((float *) vids + count);
       }
    }

   if (filetype == IEEE_F)
     fwrite(varname, sizeof(char), charsize_out, fp);
   else
      fprintf(fp,"%s\n",varname);

   if (rflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempvids64, FLOAT64, n_surf, fp);
      else
         write_ascii_double(n_surf, tempvids64);
      free(tempvids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempvids, FLOAT32, n_surf, fp);
      else
         write_ascii_float(n_surf, tempvids);
      free(tempvids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_surfvars_endsvar(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endsvar ");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"endsvar \n");
}

/* --------------------------------------------------------- */

void gmvwrite_surfflag_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "surfflag");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"surfflag \n");
}

/* --------------------------------------------------------- */

void gmvwrite_surfflag_name(char flagname[], int numtypes)
{
  int tmpnumtypes;

  tmpnumtypes = numtypes;

   if (filetype == IEEE_F)
     {
      fwrite(flagname, sizeof(char), charsize_out, fp);
      fwrite(&tmpnumtypes, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d\n",flagname,tmpnumtypes);
}

/* --------------------------------------------------------- */

void gmvwrite_surfflag_subname(char subname[])
{
   if (filetype == IEEE_F)
      fwrite(subname, sizeof(char), charsize_out, fp);
   else
      fprintf(fp,"%s\n",subname);
}

/* --------------------------------------------------------- */

void gmvwrite_surfflag_data(int flag_data[])
{
  int *tmpdata, i;

  tmpdata = (int *) malloc(sizeof(int) * n_surf);
  for (i = 0; i < n_surf; i++)
     tmpdata[i] = flag_data[i];

   if (filetype == IEEE_F)
      fwrite(tmpdata, INT32, n_surf, fp);
   else
      write_ascii_int(n_surf, tmpdata);
  free(tmpdata);
}

/* --------------------------------------------------------- */

void gmvwrite_surfflag_endflag(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endsflag");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"endsflag\n");
}

/* --------------------------------------------------------- */


void gmvwrite_units_fromfile(char *filename)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("units   fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"units   fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);

   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_units_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "units   ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"units   \n");
}

/* --------------------------------------------------------- */

void gmvwrite_units_typehdr(int data_type, int numtypes)
{
  int tmpdata_type, tmpnumtypes;

  tmpdata_type = data_type;
  tmpnumtypes = numtypes;

  if (tmpdata_type == 0)
     strcpy(tmpname, "cells   ");
  else if (tmpdata_type == 1)
     strcpy(tmpname, "nodes   ");
  else if (tmpdata_type == 2)
     strcpy(tmpname, "faces   ");

   if (filetype == IEEE_F)
     {
      fwrite(tmpname, sizeof(char), 8, fp);
      fwrite(&tmpnumtypes, sizeof(int), 1, fp);
     }
   else
      fprintf(fp,"%s %d\n",tmpname,tmpnumtypes);
}

/* --------------------------------------------------------- */

void gmvwrite_units_name(char fldname[], char unitname[])
{
   if (filetype == IEEE_F)
     {
      fwrite(fldname, sizeof(char), charsize_out, fp);
      fwrite(unitname, sizeof(char), 16, fp);
     }
   else
      fprintf(fp," %s %s\n",fldname,unitname);
}

/* --------------------------------------------------------- */

void gmvwrite_units_endunit(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endunit ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"endunit \n");
}

/* --------------------------------------------------------- */

void gmvwrite_vinfo_header(void)
 {
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "vinfo   ");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"vinfo \n");
 }

/* --------------------------------------------------------- */

void gmvwrite_vinfo_name_data(int nelem, int nlines, char varname[], void *vids)
{
  int tmpnelem, tmpnlines;
  long count, npts, j, jj;
  float *tempvids=NULL;
  double *tempvids64=NULL;

  tmpnelem = nelem;
  tmpnlines = nlines;

  npts = nelem * nlines;

  if (rflag64)
    {
     tempvids64 = (double*)malloc(sizeof(double)*npts);
    }
  else
    {
     tempvids = (float*)malloc(sizeof(float)*npts);
    }

  if (rflag64)
    {
     for (count = 0; count < npts; count++)
       {
        tempvids64[count] = *((double *) vids + count);
       }
    }
  else
    {
     for (count = 0; count < npts; count++)
       {
        tempvids[count] = *((float *) vids + count);
       }
    }

   if (filetype == IEEE_F)
     {
      fwrite(varname, sizeof(char), charsize_out, fp);
      fwrite(&tmpnelem, INT32, 1, fp);
      fwrite(&tmpnlines, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d %d\n",varname,tmpnelem,tmpnlines);
   if (rflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempvids64, FLOAT64, npts, fp);
      else
        {
         jj = 0;
         for (j = 0; j < tmpnlines; j++)
           {
            for (i = 0; i < tmpnelem; i++)
              {
               fprintf(fp,"%lg ",tempvids64[jj]);
               jj++;
              }
            fprintf(fp,"\n");
           }
        }
      free(tempvids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempvids, FLOAT32, npts, fp);
      else
        {
         jj = 0;
         for (j = 0; j < tmpnlines; j++)
           {
            for (i = 0; i < tmpnelem; i++)
              {
               fprintf(fp,"%g ",tempvids[jj]);
               jj++;
              }
            fprintf(fp,"\n");
           }
        }
      free(tempvids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_vinfo_endvinfo(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endvinfo");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"endvinfo\n");
}


/* --------------------------------------------------------- */

void gmvwrite_traceids(int ntracers, void *traceids)
{
  int i, *tmptraceids=NULL;
  long ilong, *tmptraceids64=NULL;

  if (iflag64)
    {
     tmptraceids64 = (long *) malloc(sizeof(long) * ntracers);
     for (ilong = 0; ilong < ntracers; ilong++)
	tmptraceids64[ilong] = *((long *) traceids + ilong);
    }
  else
    {
     tmptraceids = (int *) malloc(sizeof(int) * ntracers);
     for (i = 0; i < ntracers; i++)
        tmptraceids[i] = *((int *) traceids + i);
    }

   strcpy(tmpname, "traceids");
   if (filetype == IEEE_F)
      fwrite(tmpname, sizeof(char), 8, fp);
   else
      fprintf(fp,"%s |n",tmpname);

   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmptraceids64, INT64, ntracers, fp);
      else
         write_ascii_long((long)ntracers, tmptraceids64);
      free(tmptraceids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmptraceids, INT32, ntracers, fp);
      else
         write_ascii_int((long)ntracers, tmptraceids);
      free(tmptraceids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_traceids_fromfile(char *filename)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("traceidsfromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"traceidsfromfile \"%s\"",filename);
   fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_faceids(void *faceids)
{
  int i, *tempfaceids=NULL;
  long ilong, *tempfaceids64=NULL;

  if (iflag64)
    {
     tempfaceids64 = (long *) malloc(sizeof(long) * n_faces);
     for (ilong = 0; ilong < n_faces; ilong++)
	tempfaceids64[ilong] = *((long *) faceids + ilong);
    }
  else
    {
     tempfaceids = (int *) malloc(sizeof(int) * n_faces);
     for (i = 0; i < n_faces; i++)
        tempfaceids[i] = *((int *) faceids + i);
    }

   strcpy(tmpname, "faceids ");
   if (filetype == IEEE_F)
      fwrite(tmpname, sizeof(char), 8, fp);
   else
      fprintf(fp,"%s |n",tmpname);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tempfaceids64, INT64, n_faces, fp);
      else
         write_ascii_long(n_faces, tempfaceids64);
      free(tempfaceids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tempfaceids, INT32, n_faces, fp);
      else
         write_ascii_int(n_faces, tempfaceids);
      free(tempfaceids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_faceids_fromfile(char *filename)
{
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("faceids fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"faceids fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);
   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_group_fromfile(char *filename)
 {
   char *tmpbuf;

   tmpbuf = (char *)malloc((strlen(filename) + strlen("groups  fromfile ") + 3)
            * sizeof(char));
   sprintf(tmpbuf,"groups  fromfile \"%s\"",filename);

   if (filetype == IEEE_F)
      fwrite(tmpbuf,sizeof(char),strlen(tmpbuf),fp);
   else
      fprintf(fp,"%s\n",tmpbuf);

   free(tmpbuf);
   return;
}

/* --------------------------------------------------------- */

void gmvwrite_group_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "groups  ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"groups \n");
}

/* --------------------------------------------------------- */

void gmvwrite_group_data(char groupname[], int data_type, int numgrp,
     void *group_data)
{
  int tmpdata_type, tmpnumgrp, *tmpdata=NULL, *groupptr=NULL, count;
  long *tmpdata64=NULL, *groupptr64=NULL;

  tmpdata_type = data_type;
  tmpnumgrp = numgrp;

   /*  Save group data elements.  */
   if (iflag64)
     {
      tmpdata64 = (long *) malloc(sizeof(long) * tmpnumgrp);
      groupptr64 = group_data;
     }
   else
     {
      tmpdata = (int *) malloc(sizeof(int) * tmpnumgrp);
      groupptr = group_data;
     }

   if (iflag64)
      for (count = 0;count < tmpnumgrp;count++)
         tmpdata64[count] = groupptr64[count];
   else
      for (count = 0; count < tmpnumgrp; count++)
	 tmpdata[count] = groupptr[count];

   if (filetype == IEEE_F)
     {
      fwrite(groupname, sizeof(char), charsize_out, fp);
      fwrite(&tmpdata_type, INT32, 1, fp);
      fwrite(&tmpnumgrp, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d %d ",groupname,tmpdata_type,tmpnumgrp);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata64, INT64, tmpnumgrp, fp);
      else
        {
         for (i = 0; i < tmpnumgrp; i++)
            fprintf(fp,"%ld ",tmpdata64[i]);
         fprintf(fp,"\n");
        }
      free(tmpdata64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata, INT32, tmpnumgrp, fp);
      else
        {
         for (i = 0; i < tmpnumgrp; i++)
            fprintf(fp,"%d ",tmpdata[i]);
         fprintf(fp,"\n");
        }
      free(tmpdata);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_group_endgroup(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endgrp  ");
      fwrite(tmpname, sizeof(char), 8, fp);
     }
   else
     fprintf(fp,"endgrp \n");
}

/* --------------------------------------------------------- */

void gmvwrite_surfids(void *surfids)
{
  int i, *tmpsurfids=NULL;
  long ilong, *tmpsurfids64=NULL;

  if (iflag64)
    {
     tmpsurfids64 = (long *) malloc(sizeof(long) * n_surf);
     for (ilong = 0; ilong < n_surf; ilong++)
	tmpsurfids64[ilong] = *((long *) surfids + ilong);
    }
  else
    {
     tmpsurfids = (int *) malloc(sizeof(int) * n_surf);
     for (i = 0; i < n_surf; i++)
        tmpsurfids[i] = *((int *) surfids + i);
    }

   strcpy(tmpname, "surfids ");
   if (filetype == IEEE_F)
      fwrite(tmpname, sizeof(char), 8, fp);
   else
      fprintf(fp,"%s |n",tmpname);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpsurfids64, INT64, n_surf, fp);
      else
         write_ascii_long(n_surf, tmpsurfids64);
      free(tmpsurfids64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpsurfids, INT32, n_surf, fp);
      else
         write_ascii_int(n_surf, tmpsurfids);
      free(tmpsurfids);
     }
}

/* --------------------------------------------------------- */

void gmvwrite_codename(char codename[])
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "codename");
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(codename, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"codename %s\n",codename);
}

/* --------------------------------------------------------- */

void gmvwrite_codever(char codever[])
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "codever ");
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(codever, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"codever %s\n",codever);
}

/* --------------------------------------------------------- */

void gmvwrite_simdate(char simdate[])
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "simdate ");
      fwrite(tmpname,sizeof(char),8,fp);
      fwrite(simdate, sizeof(char), 8, fp);
     }
   else
      fprintf(fp,"simdate %s\n",simdate);
}

/* --------------------------------------------------------- */

void gmvwrite_subvars_header(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "subvars ");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"subvars \n");
}

/* --------------------------------------------------------- */

void gmvwrite_subvars_name_data(int data_type, int numelem, char varname[],
                                void *vids, void *vdata)
 {
  int tmpdata_type, tmpnumelem, *tmpdata=NULL, *subvarids=NULL;
  long count, npts;
  long *tmpdata64=NULL, *subvarids64=NULL;
  float *tmpvdata=NULL;
  double *tmpvdata64=NULL;

  tmpdata_type = data_type;
  tmpnumelem = numelem;

  npts = numelem;

   /*  Save subvars element ids.  */
   if (iflag64)
     {
      tmpdata64 = (long *) malloc(sizeof(long) * tmpnumelem);
      subvarids64 = vids;
     }
   else
     {
      tmpdata = (int *) malloc(sizeof(int) * tmpnumelem);
      subvarids = vids;
     }

   if (iflag64)
      for (count = 0;count < tmpnumelem;count++)
         tmpdata64[count] = subvarids64[count];
   else
      for (count = 0; count < tmpnumelem; count++)
	 tmpdata[count] = subvarids[count];

   /*  Save subvars field data.  */
  if (rflag64)
   {
    tmpvdata64 = (double*)malloc(sizeof(double)*npts);
   }
  else
   {
    tmpvdata = (float*)malloc(sizeof(float)*npts);
   }

  if (rflag64)
    {
     for (count = 0; count < npts; count++)
       {
        tmpvdata64[count] = *((double *) vdata + count);
       }
    }
  else
    {
     for (count = 0; count < npts; count++)
       {
        tmpvdata[count] = *((float *) vdata + count);
       }
    }

   /*  Write out the data.  */
   if (filetype == IEEE_F)
     {
      fwrite(varname, sizeof(char), charsize_out, fp);
      fwrite(&tmpdata_type, INT32, 1, fp);
      fwrite(&tmpnumelem, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d %d\n",varname,tmpdata_type,tmpnumelem);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata64, INT64, tmpnumelem, fp);
      else
         write_ascii_long(npts, tmpdata64);
      free(tmpdata64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata, INT32, tmpnumelem, fp);
      else
         write_ascii_int(npts, tmpdata);
      free(tmpdata);
     }
   if (rflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpvdata64, FLOAT64, npts, fp);
      else
         write_ascii_double(npts, tmpvdata64);
      free(tmpvdata64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpvdata, FLOAT32, npts, fp);
      else
         write_ascii_float(npts, tmpvdata);
      free(tmpvdata);
     }
 }

/* --------------------------------------------------------- */

void gmvwrite_subvars_endsubv(void)
{
   if (filetype == IEEE_F)
     {
      strcpy(tmpname, "endsubv ");
      fwrite(tmpname,sizeof(char),8,fp);
     }
   else
      fprintf(fp,"endsubv \n");
}

/* --------------------------------------------------------- */

void gmvwrite_ghosts(int data_type, int numghst, void *ghost_data)
{
  int tmpdata_type, tmpnumghst, *tmpdata=NULL, *ghostptr=NULL, count;
  long *tmpdata64=NULL, *ghostptr64=NULL;

  tmpdata_type = data_type;
  tmpnumghst = numghst;

   strcpy(tmpname, "groups  ");

   /*  Save ghost data elements.  */
   if (iflag64)
     {
      tmpdata64 = (long *) malloc(sizeof(long) * tmpnumghst);
      ghostptr64 = ghost_data;
     }
   else
     {
      tmpdata = (int *) malloc(sizeof(int) * tmpnumghst);
      ghostptr = ghost_data;
     }

   if (iflag64)
      for (count = 0; count < tmpnumghst; count++)
         tmpdata64[count] = ghostptr64[count];
   else
      for (count = 0; count < tmpnumghst; count++)
	 tmpdata[count] = ghostptr[count];

   if (filetype == IEEE_F)
     {
      fwrite(tmpname, sizeof(char), charsize_out, fp);
      fwrite(&tmpdata_type, INT32, 1, fp);
      fwrite(&tmpnumghst, INT32, 1, fp);
     }
   else
      fprintf(fp,"%s %d %d ",tmpname,tmpdata_type,tmpnumghst);
   if (iflag64)
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata64, INT64, tmpnumghst, fp);
      else
        {
         for (i = 0; i < tmpnumghst; i++)
            fprintf(fp,"%ld ",tmpdata64[i]);
         fprintf(fp,"\n");
        }
      free(tmpdata64);
     }
   else
     {
      if (filetype == IEEE_F)
         fwrite(tmpdata, INT32, tmpnumghst, fp);
      else
        {
         for (i = 0; i < tmpnumghst; i++)
            fprintf(fp,"%d ",tmpdata[i]);
         fprintf(fp,"\n");
        }
      free(tmpdata);
     }
}

/* --------------------------------------------------------- */

void write_ascii_float(long num, float *val)
{
  long j;
  int count;

   count = 0;
   for (j = 0; j < num; j++)
     {
      fprintf(fp,"%g ",val[j]);
      count++;
      if (count >= 10)
        {
         fprintf(fp,"\n");
         count = 0;
        }
     }
   if (count != 0)
      fprintf(fp,"\n");
}

/* --------------------------------------------------------- */

void write_ascii_double(long num, double *val)
{
  long j;
  int count;

   count = 0;
   for (j = 0; j < num; j++)
     {
      fprintf(fp,"%lg ",val[j]);
      count++;
      if (count >= 10)
        {
         fprintf(fp,"\n");
         count = 0;
        }
     }
   if (count != 0)
      fprintf(fp,"\n");
}

/* --------------------------------------------------------- */

void write_ascii_int(long num, int *val)
{
  long j;
  int count;

   count = 0;
   for (j = 0; j < num; j++)
     {
      fprintf(fp,"%d ",val[j]);
      count++;
      if (count >= 10)
        {
         fprintf(fp,"\n");
         count = 0;
        }
     }
   if (count != 0)
      fprintf(fp,"\n");
}

/* --------------------------------------------------------- */

void write_ascii_long(long num, long *val)
{
  long j;
  int count;

   count = 0;
   for (j = 0; j < num; j++)
     {
      fprintf(fp,"%ld ",val[j]);
      count++;
      if (count >= 10)
        {
         fprintf(fp,"\n");
         count = 0;
        }
     }
   if (count != 0)
      fprintf(fp,"\n");
}

/* --------------------------------------------------------- */
