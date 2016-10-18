/* pluto.cpp: functions for Pluto coordinates from Meeus theory

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

#include <math.h>
#include <stdlib.h>
#include "watdefs.h"
#include "lunar.h"

#define PI 3.1415926535897932384626433
#define N_COEFFS 36
#define COEFFS struct coeffs

COEFFS
   {
   char j, s, p, dummy_to_maintain_alignment;
   short lon_a, lon_b, lat_a, lat_b, rad_a, rad_b;
   };

/*
Lifted from p 247, Meeus, Astro Algorithms.  Note that this provides
excellent accuracy for years 1885-2099,  but is not really intended
for use outside that range. */

int DLL_FUNC calc_pluto_loc( const void FAR *data, double DLLPTR *loc,
                                          const double t, const long precision)
{
   double lat, lon, r, j, s, p, cosine, sine, arg;
   int i, prec = (int)precision;
   COEFFS FAR *tptr;
   long FAR *long_coeffs = (long FAR *)((char FAR *)data + 58610U);
   COEFFS FAR *coeffs = (COEFFS FAR *)(long_coeffs + 42);

                     /* assume t in julian centuries from J2000.0 */
   j =  34.35 + 3034.9057 * t;       /* jupiter's mean longitude */
   s =  50.08 + 1222.1138 * t;       /* saturn's mean longitude */
   p = 238.96 +  144.9600 * t;       /* pluto's mean longitude */
   j *= PI / 180.;
   s *= PI / 180.;
   p *= PI / 180.;
   lon = 238.956785 + 144.96 * t;
   lat = -3.908202;
   r = 407.247248;      /* temporarily in tenths of AUs; fixed at the end */
   for( i = 0; i < 7; i++)
      {
      long FAR *ltptr;

      if( i == 6)
         arg = j - p;
      else
         arg = (double)(i + 1) * p;
      cosine = cos( arg) * 1.e-6;
      sine = sin( arg) * 1.e-6;
      ltptr = long_coeffs + (i * 6);
      lon += (double)ltptr[0] * sine + (double)ltptr[1] * cosine;
      lat += (double)ltptr[2] * sine + (double)ltptr[3] * cosine;
      r   += (double)ltptr[4] * sine + (double)ltptr[5] * cosine;
      }
   tptr = coeffs;
   for( i = 0; i < N_COEFFS; i++, tptr++)
      if( abs( tptr->lon_a) > prec || abs( tptr->lon_b) > prec ||
          abs( tptr->lat_a) > prec || abs( tptr->lat_b) > prec ||
          abs( tptr->rad_a) > prec || abs( tptr->rad_b) > prec)
         {
         if( !tptr->j)
            arg = 0.;
         else
            arg = ((tptr->j == 1) ? j : j * (double)tptr->j);
         if( tptr->s < 0)
            arg -= (tptr->s == -1) ? s : s + s;
         if( tptr->s > 0)
            arg += (tptr->s ==  1) ? s : s + s;
         if( tptr->p)
            arg += p * (double)tptr->p;
         cosine = cos( arg) * 1.e-6;
         sine = sin( arg) * 1.e-6;
         lon += sine * (double)tptr->lon_a + cosine * (double)tptr->lon_b;
         lat += sine * (double)tptr->lat_a + cosine * (double)tptr->lat_b;
         r   += sine * (double)tptr->rad_a + cosine * (double)tptr->rad_b;
         }
   *loc++ = lon * PI / 180.;     /* cvt to radians */
   *loc++ = lat * PI / 180.;
   *loc++ = r / 10.;    /* convert back to units of AUs */
   return( 0);
}
