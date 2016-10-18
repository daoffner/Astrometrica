/* precess.cpp: functions for computing Earth precession

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

/* This implements the 'precise precession in ecliptic coords' method from
Meeus' _Astronomical Algorithms_,  which in turn comes from _Connaissance
des Temps_ pour 1984,  p. XXX and XL.  I'd previously implemented ecliptic
precession (in 'precess.cpp') by making the _equatorial_ precession matrix,
then rotating it by the obliquity.  This works,  but it occurred to me
that precession "ought to" be done in ecliptic coordinates.  In that frame,
precession is mostly just a rotation around the ecliptic pole;  doing it
in equatorial coordinates is a poor fit to the actual physics.

   That led me to examine the actual formulae,  and sure enough,  the cubic
terms for equatorial precession become huge as one goes a hundred centuries
or so into the past or future.  (If you animate the earth's polar motion
using this method,  you see the earth's pole do most of a complete
precession circle around the ecliptic pole;  then it just goes wandering
off as the cubic term diverges.)  In the ecliptic frame,  the same terms are
almost negligible.

   So... the "right" thing to do is the exact opposite of what I initially
did:  one should just implement ecliptic precession,  then rotate _that_
matrix by the obliquity.  The new setup_precession() does exactly this,  and
the results are indeed more reasonable than those from the earlier code.  (I
don't actually know how _accurate_ they are,  merely that they don't lead as
rapidly to having the earth tumble end over end.)

   Incidentally,  the previous formulae are preserved in 'testprec.cpp'
for testing/comparison purposes.

    To create a generalized matrix to precess from t1 to t2,  we create
a matrix to go from 2000 to t1;  invert it; then multiply it by a
precession matrix that goes from 2000 to t2.  (With some checks to see
if t1=2000 or t2=2000,  both common cases;  when these happen -- which
is most of the time -- we can just create one matrix, and possibly
invert it.)

   There's also some code at the end which I used to make sure that the
result from this code matched the tried,  true,  previous code (with
small deviations for dates near the present.)   */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"         /* for obliquity( ) prototype */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078

#include <stdlib.h>

static int setup_ecliptic_precession_from_j2000( double DLLPTR *matrix, const double year)
{
   const double t = (year - 2000.) / 100.;
   const double S2R = (PI / 180.) / 3600.;   /* converts arcSeconds to Radians */
   const double eta = t * (47.0029 * S2R + (-.03302 * S2R + 6.e-5 * S2R * t) * t);
   const double pie = 174.876384 * PI / 180. -
           t * (869.8089 * S2R - .03536 * S2R * t);
   const double p = t * (5029.0966 * S2R + (1.11113 * S2R - 6.e-5 * S2R * t) * t);

   set_identity_matrix( matrix);
#ifdef UNNECESSARY_MATH
   spin_matrix( matrix, matrix + 3, -pie);
#else       /* can get the same result without as much math: */
   matrix[0] = matrix[4] = cos( pie);
   matrix[1] = sin( pie);
   matrix[3] = -matrix[1];
#endif
   spin_matrix( matrix + 3, matrix + 6, -eta);
   spin_matrix( matrix + 3, matrix, -p);
   spin_matrix( matrix, matrix + 3, pie);
   return( 0);
}

#define SEMIRANDOM_GARBAGE1  314.8145501e+12
#define SEMIRANDOM_GARBAGE2  -9.19001473e-08

int DLL_FUNC setup_ecliptic_precession( double DLLPTR *matrix,
                     const double year_from, const double year_to)
{
   int rval;
         /* It's pretty common to precess a few zillion data points.  So   */
         /* it helps to cache the most recently computed precession matrix */
         /* so that repeated calls don't result in repeated computation.   */
   static double prev_year_from = SEMIRANDOM_GARBAGE1;
   static double prev_year_to   = SEMIRANDOM_GARBAGE2;
   static double prev_matrix[9];

   if( fabs( year_from - year_to) < 1.e-5)   /* dates sensibly equal; */
      {                                      /* avoid pointless math */
      set_identity_matrix( matrix);
      return( 0);
      }

   if( year_from == prev_year_from && year_to == prev_year_to)
      {
      memcpy( matrix, prev_matrix, 9 * sizeof( double));
      return( 0);
      }
               /* Similarly,  it's common to precess first  */
               /* in one direction,  then the other:        */
   if( year_from == prev_year_to && year_to == prev_year_from)
      {
      memcpy( matrix, prev_matrix, 9 * sizeof( double));
      invert_orthonormal_matrix( matrix);
      return( 0);
      }
               /* OK,  this is definitely something new: */
   if( year_from == 2000.)
      rval = setup_ecliptic_precession_from_j2000( matrix, year_to);
   else
      {
      rval = setup_ecliptic_precession_from_j2000( matrix, year_from);
      invert_orthonormal_matrix( matrix);
      if( year_to != 2000.)
         {
         double product[9], tmatrix[9];
         int i, j;

         setup_ecliptic_precession_from_j2000( tmatrix, year_to);
         for( i = 0; i < 3; i++)
            for( j = 0; j < 3; j++)
               product[j + i * 3] = matrix[i * 3    ] * tmatrix[j]
                                  + matrix[i * 3 + 1] * tmatrix[j + 3]
                                  + matrix[i * 3 + 2] * tmatrix[j + 6];
         memcpy( matrix, product, 9 * sizeof( double));
         }
      }
               /* Store matrix for likely subsequent use: */
   memcpy( prev_matrix, matrix, 9 * sizeof( double));
   prev_year_from = year_from;
   prev_year_to = year_to;
   return( rval);
}

int DLL_FUNC setup_precession( double DLLPTR *matrix, const double year_from, const double year_to)
{
   const double obliquity1 = mean_obliquity( (year_from - 2000.) / 100.);
   const double obliquity2 = mean_obliquity( (year_to - 2000.) / 100.);

   setup_ecliptic_precession( matrix, year_from, year_to);
   pre_spin_matrix( matrix + 1, matrix + 2, obliquity1);
   spin_matrix( matrix + 3, matrix + 6, obliquity2);
   return( 0);
}

static const double sin_obliq_2000 = .397777156;
static const double cos_obliq_2000 = .917482062;

void DLL_FUNC equatorial_to_ecliptic( double *vect)
{
   double temp;

   temp    = vect[2] * cos_obliq_2000 - vect[1] * sin_obliq_2000;
   vect[1] = vect[1] * cos_obliq_2000 + vect[2] * sin_obliq_2000;
   vect[2] = temp;
}

void DLL_FUNC ecliptic_to_equatorial( double *vect)
{
   double temp;

   temp    = vect[2] * cos_obliq_2000 + vect[1] * sin_obliq_2000;
   vect[1] = vect[1] * cos_obliq_2000 - vect[2] * sin_obliq_2000;
   vect[2] = temp;
}

int DLL_FUNC precess_vector( const double DLLPTR *matrix,
                                      const double DLLPTR *v1,
                                      double DLLPTR *v2)
{
   int i = 3;

   while( i--)
      {
      *v2++ = matrix[0] * v1[0] + matrix[1] * v1[1] + matrix[2] * v1[2];
      matrix += 3;
      }
   return( 0);
}

int DLL_FUNC deprecess_vector( const double DLLPTR *matrix,
                                      const double DLLPTR *v1,
                                      double DLLPTR *v2)
{
   int i = 3;

   while( i--)
      {
      *v2++ = matrix[0] * v1[0] + matrix[3] * v1[1] + matrix[6] * v1[2];
      matrix++;
      }
   return( 0);
}

int DLL_FUNC precess_ra_dec( const double DLLPTR *matrix,
                        double DLLPTR *p_out,
                        const double DLLPTR *p_in, int backward)
{
   double v1[3], v2[3];
   const double old_ra = p_in[0];

   v1[0] = cos( p_in[0]) * cos( p_in[1]);
   v1[1] = sin( p_in[0]) * cos( p_in[1]);
   v1[2] =                 sin( p_in[1]);
   if( backward)
      deprecess_vector( matrix, v1, v2);
   else
      precess_vector( matrix, v1, v2);
   if( v2[1] || v2[0])
      p_out[0] = atan2( v2[1], v2[0]);
   else
      p_out[0] = 0.;
   p_out[1] = asine( v2[2]);
   while( p_out[0] - old_ra > PI)
      p_out[0] -= PI * 2.;
   while( p_out[0] - old_ra <-PI)
      p_out[0] += PI * 2.;
   return( 0);
}
