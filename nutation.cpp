/* nutation.cpp: implements 1980 IAU nutation theory

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

#ifdef TEST_PROGRAM
#include <stdio.h>
#endif
#include <math.h>
#include "watdefs.h"
#include "lunar.h"

/* 21 Sep 2003:  James Miller and Mark Huss pointed out an error in the
coefficients for the 52rd term,  now shown (correctly) as 0, 0, 3, 2, 2.
This could lead to an error in the nutation in longitude of up to
.0006 degrees, or about 2.5 arcseconds.  In order to fix this,  I had
to make some small improvements to "de-gibberish-ize" the code,  at
least a little. */

#define N_NUTATION_COEFFS 62
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define CVT (PI / 180.)
#define NUTATION_COMPACT( D, M, Mp, F, Om)   \
          ((D+3) * 7*7*7*7 + (M+3) * 7*7*7 + (Mp+3) * 7*7 + (F+3) * 7 + Om+3)

      /* nutation formula comes from p 132-5, Meeus,  Astro Algor */
      /* input is time in julian centuries from 2000. */
      /* *d_lon is nutation (delta phi) in arcseconds */
      /* *d_obliq is nutation (delta epsilon) in arcseconds */
      /* Either pointer can be NULL,  in which case that value is  */
      /* not computed.  (I added this because sometimes,  you want */
      /* only d_lon or d_obliq;  in such cases,  computing _both_  */
      /* is a waste of perfectly good CPU time) */
      /*    Note that the coefficients are "compacted" using the above */
      /* NUTATION_COMPACT macro.  Had the coeffs been stored as five   */
      /* bytes per line instead of compacted into a 16-bit int,  three */
      /* bytes per line would have been wasted,  for a total of 186    */
      /* bytes.  Back in the early 1990s,  when this code was written, */
      /* that actually mattered.  It's obviously not something I'd do  */
      /* today.  But I'm reluctant to modify code that works Just Fine. */

int DLL_FUNC nutation( const double t, double DLLPTR *d_lon,
                                       double DLLPTR *d_obliq)
{
   static double linear_part[5] = {445267.111480, 35999.050340,
            477198.867398, 483202.017538, -1934.136261 };
   static long coeffs[5 * 3] = { 29785036L, -19142L,  189474L,
                                 35752772L, - 1603L, -300000L,
                                 13496298L,  86972L, 56250L,
                                  9327191L, -36825L, 327270L,
                                 12504452L,  20708L, 450000L };
   static short args[3 * N_NUTATION_COEFFS + 1] = {
         NUTATION_COMPACT(-2, 0, 0, 2, 2),-13187,5736,    /* 01 */
         NUTATION_COMPACT( 0, 0, 0, 2, 2), -2274, 977,    /* 02 */
         NUTATION_COMPACT( 0, 0, 0, 0, 2),  2062,-895,    /* 03 */
         NUTATION_COMPACT( 0, 1, 0, 0, 0),  1426,  54,    /* 04 */
         NUTATION_COMPACT( 0, 0, 1, 0, 0),   712,  -7,    /* 05 */
         NUTATION_COMPACT(-2, 1, 0, 2, 2),  -517, 224,    /* 06 */
         NUTATION_COMPACT( 0, 0, 0, 2, 1),  -386, 200,    /* 07 */
         NUTATION_COMPACT( 0, 0, 1, 2, 2),  -301, 129,    /* 08 */
         NUTATION_COMPACT(-2,-1, 0, 2, 2),   217, -95,    /* 09 */
         NUTATION_COMPACT(-2, 0, 1, 0, 0),  -158,   0,    /* 10 */
         NUTATION_COMPACT(-2, 0, 0, 2, 1),   129, -70,    /* 11 */
         NUTATION_COMPACT( 0, 0,-1, 2, 2),   123, -53,    /* 12 */
         NUTATION_COMPACT( 2, 0, 0, 0, 0),    63,   0,    /* 13 */
         NUTATION_COMPACT( 0, 0, 1, 0, 1),    63, -33,    /* 14 */
         NUTATION_COMPACT( 2, 0,-1, 2, 2),   -59,  26,    /* 15 */
         NUTATION_COMPACT( 0, 0,-1, 0, 1),   -58,  32,    /* 16 */
         NUTATION_COMPACT( 0, 0, 1, 2, 1),   -51,  27,    /* 17 */
         NUTATION_COMPACT(-2, 0, 2, 0, 0),    48,   0,    /* 18 */
         NUTATION_COMPACT( 0, 0,-2, 2, 1),    46, -24,    /* 19 */
         NUTATION_COMPACT( 2, 0, 0, 2, 2),   -38,  16,    /* 20 */
         NUTATION_COMPACT( 0, 0, 2, 2, 2),   -31,  13,    /* 21 */
         NUTATION_COMPACT( 0, 0, 2, 0, 0),    29,   0,    /* 22 */
         NUTATION_COMPACT(-2, 0, 1, 2, 2),    29, -12,    /* 23 */
         NUTATION_COMPACT( 0, 0, 0, 2, 0),    26,   0,    /* 24 */
         NUTATION_COMPACT(-2, 0, 0, 2, 0),   -22,   0,    /* 25 */
         NUTATION_COMPACT( 0, 0,-1, 2, 1),    21, -10,    /* 26 */
         NUTATION_COMPACT( 0, 2, 0, 0, 0),    17,   0,    /* 27 */
         NUTATION_COMPACT( 2, 0,-1, 0, 1),    16,  -8,    /* 28 */
         NUTATION_COMPACT(-2, 2, 0, 2, 2),   -16,   7,    /* 29 */
         NUTATION_COMPACT( 0, 1, 0, 0, 1),   -15,   9,    /* 30 */
         NUTATION_COMPACT(-2, 0, 1, 0, 1),   -13,   7,    /* 31 */
         NUTATION_COMPACT( 0,-1, 0, 0, 1),   -12,   6,    /* 32 */
         NUTATION_COMPACT( 0, 0, 2,-2, 0),    11,   0,    /* 33 */
         NUTATION_COMPACT( 2, 0,-1, 2, 1),   -10,   5,    /* 34 */
         NUTATION_COMPACT( 2, 0, 1, 2, 2),    -8,   3,    /* 35 */
         NUTATION_COMPACT( 0, 1, 0, 2, 2),     7,  -3,    /* 36 */
         NUTATION_COMPACT(-2, 1, 1, 0, 0),    -7,   0,    /* 37 */
         NUTATION_COMPACT( 0,-1, 0, 2, 2),    -7,   3,    /* 38 */
         NUTATION_COMPACT( 2, 0, 0, 2, 1),    -7,   3,    /* 39 */
         NUTATION_COMPACT( 2, 0, 1, 0, 0),     6,   0,    /* 40 */
         NUTATION_COMPACT(-2, 0, 2, 2, 2),     6,  -3,    /* 41 */
         NUTATION_COMPACT(-2, 0, 1, 2, 1),     6,  -3,    /* 42 */
         NUTATION_COMPACT( 2, 0,-2, 0, 1),    -6,   3,    /* 43 */
         NUTATION_COMPACT( 2, 0, 0, 0, 1),    -6,   3,    /* 44 */
         NUTATION_COMPACT( 0,-1, 1, 0, 0),     5,   0,    /* 45 */
         NUTATION_COMPACT(-2,-1, 0, 2, 1),    -5,   3,    /* 46 */
         NUTATION_COMPACT(-2, 0, 0, 0, 1),    -5,   3,    /* 47 */
         NUTATION_COMPACT( 0, 0, 2, 2, 1),    -5,   3,    /* 48 */
         NUTATION_COMPACT( 0, 0,-2, 2, 2),    -3,   0,    /* 49 */
         NUTATION_COMPACT(-2, 0, 2, 0, 1),     4,   0,    /* 50 */
         NUTATION_COMPACT(-2, 1, 0, 2, 1),     4,   0,    /* 51 */
         NUTATION_COMPACT( 0, 0, 3, 2, 2),    -3,   0,    /* 52 */
         NUTATION_COMPACT( 2,-1,-1, 2, 2),    -3,   0,    /* 53 */
         NUTATION_COMPACT( 0,-1, 1, 2, 2),    -3,   0,    /* 54 */
         NUTATION_COMPACT( 2,-1, 0, 2, 2),    -3,   0,    /* 55 */
         NUTATION_COMPACT(-1,-1, 1, 0, 0),    -3,   0,    /* 56 */
         NUTATION_COMPACT(-1, 0, 1, 0, 0),    -4,   0,    /* 57 */
         NUTATION_COMPACT(-2, 1, 0, 0, 0),    -4,   0,    /* 58 */
         NUTATION_COMPACT( 0, 0, 1,-2, 0),     4,   0,    /* 59 */
         NUTATION_COMPACT( 1, 0, 0, 0, 0),    -4,   0,    /* 60 */
         NUTATION_COMPACT( 0, 1, 1, 0, 0),    -3,   0,    /* 61 */
         NUTATION_COMPACT( 0, 0, 1, 2, 0),     3,   0,    /* 62 */
         0 };
   static char time_dependent[16 + 9] = {
    -16, -2, 2, -34, 1, 12, -4, 0, -5, 0, 1, 0, 0, 1, 0, -1,
    -31, -5, 5, -1, 0, -6, 0, -1, 3 };

   double terms[5];
   double t2;
   int i;

   t2 = t * t;
   for( i = 0; i < 5; i++)
      {
      terms[i] = linear_part[i] * t + (double)coeffs[i * 3] / 100000.;
      terms[i] += t2 * (double)coeffs[i * 3 + 1] * 1.e-7;
      terms[i] += t2 * t / (double)coeffs[i * 3 + 2];
      terms[i] *= PI / 180.;
      }

            /* The largest term in d_lon won't fit into a short int. */
            /* Everything else does,  luckily. */
   if( d_lon)
      *d_lon = (-171996. - 174.2 * t) * sin( terms[4]);
   if( d_obliq)
      *d_obliq = (92025. + 8.9 * t) * cos( terms[4]);

   for( i = 0; args[i]; i += 3)
      {
      double total_arg = 0., coeff;
      int j, mult = args[i];

      for( j = 4; j >= 0; j--, mult /= 7)
         if( mult % 7 != 3)
            total_arg += (double)( mult % 7 - 3) * terms[j];

      coeff = (double)args[i + 1];
      if( i < 16 && time_dependent[i])
         coeff += (double)time_dependent[i] * t / 10.;
      if( i == 26 || i == 28)
         coeff += (double)(27 - i) * t / 10.;
      if( d_lon)
         *d_lon += coeff * sin( total_arg);

      if( args[i + 2])
         {
         coeff = (double)args[i + 2];
         if( i < 9 && time_dependent[i + 16])
            coeff += (double)time_dependent[i + 16] * t / 10;
         if( d_obliq)
            *d_obliq += coeff * cos( total_arg);
         }
      }
   if( d_lon)
      *d_lon *= .0001;
   if( d_obliq)
      *d_obliq *= .0001;
   return( 0);
}

/* To make the test program,  use */

/* cl -W4 -Ox -AT -I\myincl -DTEST_PROGRAM nutation.cpp */

#ifdef TEST_PROGRAM

void main( int argc, char **argv)
{
   double d_lon, d_obliq;
   double t = (atof( argv[1]) - 2451545.) / 36525.;

   nutation( t, &d_lon, &d_obliq);
   printf( "%lf %lf\n", d_lon, d_obliq);
}
#endif
