/* weight.cpp: handle setting of weights for observations

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "watdefs.h"
#include "weight.h"
#include "date.h"

#define WEIGHT_RECORD struct weight_record

WEIGHT_RECORD
   {
   double jd1, jd2, weight;
   double mag_weight;
   int mag1, mag2;
   char mpc_code[5];
   };

static int n_weight_recs;
static WEIGHT_RECORD *weight_recs;

int debug_printf( const char *format, ...);                /* runge.cpp */

static int parse_weight_record( WEIGHT_RECORD *w, const char *buff)
{
   int rval = 0;

   if( *buff == ' ')
      {
      int i;

      memcpy( w->mpc_code, buff + 1, 3);
      w->mpc_code[3] = '\0';
      w->jd1 = 0.;
      w->jd2 = 3000000.;
      w->mag1 = -1;
      w->mag2 = 300;
      for( i = 0; i < 2; i++)
         if( buff[i * 11 + 6] != ' ')
            {
            const double jd = (double)dmy_to_day( atoi( buff + i * 11 + 14),
                                                  atoi( buff + i * 11 + 11),
                                                  atoi( buff + i * 11 + 6),
                                                  CALENDAR_JULIAN_GREGORIAN);

            if( i)
               w->jd2 = jd;
            else
               w->jd1 = jd;
            }
      w->weight = atof( buff + 38);
      if( buff[41] == 'm' && buff[42] == ':')
         w->mag_weight = atof( buff + 43);
      else
         w->mag_weight = 1.;
      rval = 1;
      }
   return( rval);
}

int load_up_weight_records( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");

   if( ifile)
      {
      int i, j = 0;

      for( i = 0; i < 2; i++)
         {
         char buff[120];
         WEIGHT_RECORD w;

         fseek( ifile, 0L, SEEK_SET);
         while( fgets( buff, sizeof( buff), ifile))
            if( parse_weight_record( &w, buff))
               {
               if( !i)
                  n_weight_recs++;
               else
                  weight_recs[j++] = w;
               }
         if( !i)
            weight_recs = (WEIGHT_RECORD *)calloc( n_weight_recs,
                                             sizeof( WEIGHT_RECORD));
         if( !weight_recs)
            debug_printf( "%d weight recs not alloced\n", n_weight_recs);
         }
      fclose( ifile);
      }
   return( n_weight_recs);
}

void free_weight_recs( void)
{
   if( weight_recs)
      {
      free( weight_recs);
      weight_recs = NULL;
      }
   n_weight_recs = 0;
}

double get_observation_weight( const double jd, const int mag_in_tenths,
                  const char *mpc_code, double *mag_weight)
{
   int i;

   if( mag_weight)
      *mag_weight = 1.;
   for( i = 0; i < n_weight_recs; i++)
      {
      WEIGHT_RECORD *w = weight_recs + i;

      if( !memcmp( mpc_code, w->mpc_code, 3)
                              || !memcmp( "   ", w->mpc_code, 3))
         if( jd > w->jd1 && jd < w->jd2)
            if( mag_in_tenths > w->mag1 && mag_in_tenths < w->mag2)
               {
               if( mag_weight)
                  *mag_weight = 1.;
               return( w->weight);
               }
      }
   return( 1.);         /* shouldn't actually reach here */
}
