/* lunar2.cpp: functions for modest-precision lunar coords

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

/* Implements a simplified lunar ephemeris via the method described
in Meeus' _Astronomical Algorithms_.  The actual series coefficients
are stored in 'vsop.bin',  which must be read into a buffer before
you call these functions.  At the time I wrote all this -- early
1990s -- memory was a scarce resource,  so I packed bytes as much
as possible.  I suppose all this would still be a good idea on some
embedded systems.  It's probably not the way I would do things now.
The code is not easy to follow.  But it _does_ all work,  and runs
fast and has a small footprint.     */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "watdefs.h"
#include "lunar.h"

#define Lp  fund[0]
#define D   fund[1]
#define M   fund[2]
#define Mp  fund[3]
#define F   fund[4]
#define A1  fund[5]
#define A2  fund[6]
#define A3  fund[7]
#define T   fund[8]

#define N_TERM1 60
#define TERM1 struct term1

/* 28 Jul 2011:  modified these structures so they're packed (they
   are read from disk files) and that 32-bit ints are specified.
   Obviously,  they'll still break on non-Intel order hardware.  */

#pragma pack( 1)
TERM1
   {
   char d, m, mp, f;
   int32_t sl, sr;
   };
#pragma pack( )

#define N_TERM2 60
#define TERM2 struct term2

#pragma pack( 1)
TERM2
   {
   char d, m, mp, f;
   int32_t sb;
   };
#pragma pack( )

#define PI 3.1415926535897932384626433
#define CVT (PI / 180.)

int DLL_FUNC lunar_fundamentals( const void FAR *data, const double t,
                                        double DLLPTR *fund)
{
   int i, j;
   const double FAR *tptr = (const double FAR *)((const char FAR *)data + 60554U);
   double tpow;

   for( i = 0; i < 5; i++)
      {
      fund[i] = *tptr++;
      tpow = t;
      for( j = 4; j; j--, tpow *= t)
         fund[i] += tpow * (*tptr++);
      }

   A1 = 119.75 + 131.849 * t;
   A2 =  53.09 + 479264.290 * t;
   A3 = 313.45 + 481266.484 * t;
   T = t;
   for( i = 0; i < N_FUND - 1; i++)      /* convert to radians */
      {
      fund[i] = fmod( fund[i], 360.);
      if( fund[i] < 0.) fund[i] += 360.;
      fund[i] *= CVT;
      }
   return( 0);
}

int DLL_FUNC lunar_lon_and_dist( const void FAR *data, const double DLLPTR *fund,
                 double DLLPTR *lon, double DLLPTR *r, const long precision)
{
   int i, j;
   const TERM1 FAR *term1 = (const TERM1 FAR *)((const char FAR *)data + 59354U);
   const TERM1 FAR *tptr;
   double sl = 0., sr = 0., e;

   tptr = term1;
   e = 1. - .002516 * T - .0000074 * T * T;
   for( i = N_TERM1; i; i--, tptr++)
      if( labs( tptr->sl) > precision || labs( tptr->sr) > precision)
         {
         double arg, term;

         switch( tptr->d)
            {
            case  1:   arg = D;     break;
            case -1:   arg =-D;     break;
            case  2:   arg = D+D;   break;
            case -2:   arg =-D-D;   break;
            case  0:   arg = 0.;    break;
            default:   arg = (double)tptr->d * D;  break;
            }
         switch( tptr->m)
            {
            case  1:   arg += M;     break;
            case -1:   arg -= M;     break;
            case  2:   arg += M+M;   break;
            case -2:   arg -= M+M;  break;
            case  0:           ;    break;
            default:   arg += (double)tptr->m * M;  break;
            }
         switch( tptr->mp)
            {
            case  1:   arg += Mp;      break;
            case -1:   arg -= Mp;      break;
            case  2:   arg += Mp+Mp;   break;
            case -2:   arg -= Mp+Mp;   break;
            case  0:           ;       break;
            default:   arg += (double)tptr->mp * Mp;  break;
            }
         switch( tptr->f)
            {
            case  1:   arg += F;     break;
            case -1:   arg -= F;     break;
            case  2:   arg += F+F;   break;
            case -2:   arg -= F+F;  break;
            case  0:           ;    break;
            default:   arg += (double)tptr->f * F;  break;
            }
         if( tptr->sl)
            {
            term = (double)tptr->sl * sin( arg);
            for( j = abs( tptr->m); j; j--)
               term *= e;
            sl += term;
            }
         if( tptr->sr)
            {
            term = (double)tptr->sr * cos( arg);
            for( j = abs( tptr->m); j; j--)
               term *= e;
            sr += term;
            }
         }
   if( precision < 3959L)
      sl += 3958. * sin( A1) + 1962. * sin( Lp - F) + 318. * sin( A2);
   *lon = (Lp * 180. / PI) + sl * 1.e-6;
   while( *lon < 0.)
      *lon += 360.;
   while( *lon > 360.)
      *lon -= 360.;
   *r = 385000.56 + sr / 1000.;
   return( 0);
}

double DLL_FUNC lunar_lat( const void FAR *data, const double DLLPTR *fund,
                                           const long precision)
{
   int i, j;
   const TERM2 FAR *tptr;
   const TERM2 FAR *term2 = (const TERM2 FAR *)((const char FAR *)data + 60074U);
   double rval = 0., e;

   tptr = term2;
   e = 1. - .002516 * T - .0000074 * T * T;
   for( i = N_TERM2; i; i--, tptr++)
      if( labs( tptr->sb) > precision)
         {
         double arg, term;

         switch( tptr->d)
            {
            case  1:   arg = D;     break;
            case -1:   arg =-D;     break;
            case  2:   arg = D+D;   break;
            case -2:   arg =-D-D;   break;
            case  0:   arg = 0.;    break;
            default:   arg = (double)tptr->d * D;  break;
            }
         switch( tptr->m)
            {
            case  1:   arg += M;     break;
            case -1:   arg -= M;     break;
            case  2:   arg += M+M;   break;
            case -2:   arg -= M+M;  break;
            case  0:           ;    break;
            default:   arg += (double)tptr->m * M;  break;
            }
         switch( tptr->mp)
            {
            case  1:   arg += Mp;      break;
            case -1:   arg -= Mp;      break;
            case  2:   arg += Mp+Mp;   break;
            case -2:   arg -= Mp+Mp;   break;
            case  0:           ;       break;
            default:   arg += (double)tptr->mp * Mp;  break;
            }
         switch( tptr->f)
            {
            case  1:   arg += F;     break;
            case -1:   arg -= F;     break;
            case  2:   arg += F+F;   break;
            case -2:   arg -= F+F;  break;
            case  0:           ;    break;
            default:   arg += (double)tptr->f * F;  break;
            }
         term = (double)tptr->sb * sin( arg);
         for( j = abs( tptr->m); j; j--)
            term *= e;
         rval += term;
         }
   if( precision < 2236L)
      rval += -2235. * sin( Lp) + 382. * sin( A3) + 175. * sin( A1 - F) +
               175. * sin( A1 + F) + 127. * sin(Lp - Mp) - 115. * sin(Lp+Mp);
   return( rval * 1.e-6);
}
