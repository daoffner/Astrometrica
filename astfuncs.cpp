/* astfuncs.cpp: functions for asteroid/comet two-body ephems

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
#include <stddef.h>
#include <stdint.h>
#include "watdefs.h"
#include "comets.h"

#define TEN_MILLION 1.e+7
#define HUND_MILLION 1.e+8
#define LARGE_AXIS      ((uint32_t)(2100000000U))
#define VERY_LARGE_AXIS ((uint32_t)(3150000000U))
#define PI 3.141592653589793238462643383279502884197169399375105
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define SQRT_2 1.414213562
#define THRESH 1.e-8
#define MIN_THRESH 1.e-15
#define CUBE_ROOT( X)  (exp( log( X) / 3.))

static double kepler( const double ecc, double mean_anom);
void setup_orbit_vectors( ELEMENTS DLLPTR *e);  /* astfuncs.cpp */
void comet_posn_part_ii( const ELEMENTS DLLPTR *elem, const double t,
                                    double DLLPTR *loc, double DLLPTR *vel);

/* Asteroid elements on the Guide CD-ROM are stored in a compressed format,
where a full set of elements consumes six long integers = 24 bytes.  The
grisly details are described in the file \COMPRESS\ASTEROID.DOC on the
Guide CD-ROM.        */

int DLL_FUNC setup_elems_from_ast_file( ELEMENTS DLLPTR *class_elem,
              const uint32_t DLLPTR *elem, const double t_epoch)
{
   double mean_anomaly;

   mean_anomaly           = (PI/180.) * (double)elem[0] / TEN_MILLION;
   class_elem->asc_node   = (PI/180.) * (double)elem[3] / TEN_MILLION;
   class_elem->arg_per    = (PI/180.) * (double)elem[4] / TEN_MILLION;
   class_elem->incl       = (PI/180.) * (double)elem[5] / TEN_MILLION;

   class_elem->major_axis = (double)elem[1] / HUND_MILLION;
               /* Originally,  the above line was sufficient;  it handled  */
               /* major axes out to 42.9 AU.  When more distant asteroids  */
               /* were discovered,  the code was modified to work out to   */
               /* 84 AU;  then a final modification was made to allow it   */
               /* to work out to infinity... which is,  admittedly,  what  */
               /* I should have done right from the beginning.             */
   if( elem[1] > VERY_LARGE_AXIS)   /* kludge to accommodate large axes */
      {
      double tval = 4. - class_elem->major_axis / 10.5;

      class_elem->major_axis = 63. / tval;
      }
   else if( elem[1] > LARGE_AXIS)   /* kludge to accommodate large axes */
      class_elem->major_axis = class_elem->major_axis * 4. - 63.;

   class_elem->ecc        = (double)elem[2] / HUND_MILLION;
   class_elem->q = class_elem->major_axis * (1. - class_elem->ecc);
   class_elem->epoch = t_epoch;
   class_elem->mean_anomaly = mean_anomaly;
   derive_quantities( class_elem, SOLAR_GM);
   class_elem->perih_time = t_epoch - mean_anomaly * class_elem->t0;
   class_elem->central_obj = 0;     /* all asteroids orbit the sun */
   return( 0);
}

/* We want to have,  in the ELEMENTS structure,  the ratio of the minor
and major axes;  the longitude of perihelion;  and a unit vector,
"sideways",  that lies in the plane of the orbit and points at right angles
to the direction of perihelion. */

void setup_orbit_vectors( ELEMENTS DLLPTR *e)
{
   const double sin_incl = sin( e->incl), cos_incl = cos( e->incl);
   double FAR *vec;
   double vec_len;
   double up[3];
   int i;

   e->minor_to_major = sqrt( fabs( 1. - e->ecc * e->ecc));
   e->lon_per = e->asc_node + atan2( sin( e->arg_per) * cos_incl,
                                       cos( e->arg_per));
   vec = e->perih_vec;

#ifdef DIV_BY_ZERO_IF_INCL_EQ_90
   vec[0] = cos( e->lon_per);
   vec[1] = sin( e->lon_per);
   vec[2] = (sin_incl / cos_incl) * sin( e->lon_per - e->asc_node);
   vec_len = sqrt( 1. + vec[2] * vec[2]);
#else
   vec[0] = cos( e->lon_per) * cos_incl;
   vec[1] = sin( e->lon_per) * cos_incl;
   vec[2] = sin_incl * sin( e->lon_per - e->asc_node);
   vec_len = sqrt( cos_incl * cos_incl + vec[2] * vec[2]);
   if( cos_incl < 0.)      /* for retrograde cases, make sure */
      vec_len *= -1.;      /* 'vec' has correct orientation   */
#endif
   for( i = 0; i < 3; i++)
      vec[i] /= vec_len;
            /* 'up' is a vector perpendicular to the plane of the orbit */
   up[0] =  sin( e->asc_node) * sin_incl;
   up[1] = -cos( e->asc_node) * sin_incl;
   up[2] = cos_incl;

   e->sideways[0] = up[1] * vec[2] - up[2] * vec[1];
   e->sideways[1] = up[2] * vec[0] - up[0] * vec[2];
   e->sideways[2] = up[0] * vec[1] - up[1] * vec[0];
}

void DLL_FUNC derive_quantities( ELEMENTS DLLPTR *e, const double gm)
{
   if( e->ecc != 1.)    /* for non-parabolic orbits: */
      {
      e->major_axis = e->q / fabs(1. - e->ecc);
      e->t0 = e->major_axis * sqrt( e->major_axis / gm);
      }
   else
      {
      e->w0 = (3. / SQRT_2) / (e->q * sqrt( e->q / gm));
      e->major_axis = e->t0 = 0.;
      }
  setup_orbit_vectors( e);
}

/* 'asinh' = 'arc-hyperbolic sine.'  Most compilers now implement this.
I _think_ the Microsoft compiler is the only one with this problem. */

#ifdef _MSC_VER
static double asinh( const double z)
{
   return( log( z + sqrt( z * z + 1.)));
}
#endif

/* If the eccentricity is very close to parabolic,  and the eccentric
anomaly is quite low,  you can get an unfortunate situation where
roundoff error keeps you from converging.  Consider the just-barely-
elliptical case,  where in Kepler's equation,

M = E - e sin( E)

   E and e sin( E) can be almost identical quantities.  To
around this,  near_parabolic( ) computes E - e sin( E) by expanding
the sine function as a power series:

E - e sin( E) = E - e( E - E^3/3! + E^5/5! - ...)
= (1-e)E + e( -E^3/3! + E^5/5! - ...)

   It's a little bit expensive to do this,  and you only need do it
quite rarely.  (I only encountered the problem because I had orbits
that were supposed to be 'pure parabolic',  but due to roundoff,
they had e = 1+/- epsilon,  with epsilon _very_ small.)  So 'near_parabolic'
is only called if we've gone seven iterations without converging. */

static double near_parabolic( const double ecc_anom, const double e)
{
   const double anom2 = (e > 1. ? ecc_anom * ecc_anom : -ecc_anom * ecc_anom);
   double term = e * anom2 * ecc_anom / 6.;
   double rval = (1. - e) * ecc_anom - term;
   int n = 4;

   while( fabs( term) > 1e-15)
      {
      term *= anom2 / (double)(n * (n + 1));
      rval -= term;
      n += 2;
      }
   return( rval);
}

/* For a full description of this function,  see KEPLER.HTM on the Guide
Web site,  http://www.projectpluto.com.  There was a long thread about
solutions to Kepler's equation on sci.astro.amateur,  and I decided to
go into excruciating detail as to how it's done below. */

#define MAX_ITERATIONS 7

static double kepler( const double ecc, double mean_anom)
{
   double curr, err, thresh, offset = 0.;
   double delta_curr = 1.;
   int is_negative = 0, n_iter = 0;

   if( !mean_anom)
      return( 0.);

   if( ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
      {
      curr = atan2( sin( mean_anom), cos( mean_anom) - ecc);
            /* two correction steps,  and we're done */
      for( n_iter = 2; n_iter; n_iter--)
         {
         err = curr - ecc * sin( curr) - mean_anom;
         curr -= err / (1. - ecc * cos( curr));
         }
      return( curr);
      }

   if( ecc < 1.)
      if( mean_anom < -PI || mean_anom > PI)
         {
         double tmod = fmod( mean_anom, PI * 2.);

         if( tmod > PI)             /* bring mean anom within -pi to +pi */
            tmod -= 2. * PI;
         else if( tmod < -PI)
            tmod += 2. * PI;
         offset = mean_anom - tmod;
         mean_anom = tmod;
         }

   if( mean_anom < 0.)
      {
      mean_anom = -mean_anom;
      is_negative = 1;
      }

   curr = mean_anom;
   thresh = THRESH * fabs( 1. - ecc);
               /* Due to roundoff error,  there's no way we can hope to */
               /* get below a certain minimum threshhold anyway:        */
   if( thresh < MIN_THRESH)
      thresh = MIN_THRESH;
   if( (ecc > .8 && mean_anom < PI / 3.) || ecc > 1.)    /* up to 60 degrees */
      {
      double trial = mean_anom / fabs( 1. - ecc);

      if( trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
         {
         if( mean_anom < PI)
            trial = CUBE_ROOT( 6. * mean_anom);
         else        /* hyperbolic w/ 5th & higher-order terms predominant */
            trial = asinh( mean_anom / ecc);
         }
      curr = trial;
      }
   if( ecc > 1. && mean_anom > 4.)    /* hyperbolic, large-mean-anomaly case */
      curr = log( mean_anom);

   if( ecc < 1.)
      while( fabs( delta_curr) > thresh)
         {
         if( n_iter++ > MAX_ITERATIONS)
            err = near_parabolic( curr, ecc) - mean_anom;
         else
            err = curr - ecc * sin( curr) - mean_anom;
         delta_curr = -err / (1. - ecc * cos( curr));
         curr += delta_curr;
         }
   else
      while( fabs( delta_curr) > thresh)
         {
         if( n_iter++ > MAX_ITERATIONS)
            err = -near_parabolic( curr, ecc) - mean_anom;
         else
            err = ecc * sinh( curr) - curr - mean_anom;
         delta_curr = -err / (ecc * cosh( curr) - 1.);
         curr += delta_curr;
         }
   return( is_negative ? offset - curr : offset + curr);
}

void comet_posn_part_ii( const ELEMENTS DLLPTR *elem, const double t,
                                    double DLLPTR *loc, double DLLPTR *vel)
{
   double true_anom, r, x, y, r0;

   if( elem->ecc == 1.)    /* parabolic */
      {
      double g = elem->w0 * t * .5;

      y = CUBE_ROOT( g + sqrt( g * g + 1.));
      true_anom = 2. * atan( y - 1. / y);
      }
   else           /* got the mean anomaly;  compute eccentric,  then true */
      {
      double ecc_anom;

      ecc_anom = kepler( elem->ecc, elem->mean_anomaly);
      if( elem->ecc > 1.)     /* hyperbolic case */
         {
         x = (elem->ecc - cosh( ecc_anom));
         y = sinh( ecc_anom);
         }
      else           /* elliptical case */
         {
         x = (cos( ecc_anom) - elem->ecc);
         y =  sin( ecc_anom);
         }
      y *= elem->minor_to_major;
      true_anom = atan2( y, x);
      }

   r0 = elem->q * (1. + elem->ecc);
   r = r0 / (1. + elem->ecc * cos( true_anom));
   x = r * cos( true_anom);
   y = r * sin( true_anom);
   loc[0] = elem->perih_vec[0] * x + elem->sideways[0] * y;
   loc[1] = elem->perih_vec[1] * x + elem->sideways[1] * y;
   loc[2] = elem->perih_vec[2] * x + elem->sideways[2] * y;
   loc[3] = r;
   if( vel && (elem->angular_momentum != 0.))
      {
      double angular_component = elem->angular_momentum / (r * r);
      double radial_component = elem->ecc * sin( true_anom) *
                                elem->angular_momentum / (r * r0);
      double x1 = x * radial_component - y * angular_component;
      double y1 = y * radial_component + x * angular_component;
      int i;

      for( i = 0; i < 3; i++)
         vel[i] = elem->perih_vec[i] * x1 + elem->sideways[i] * y1;
      }
}

int DLL_FUNC comet_posn_and_vel( ELEMENTS DLLPTR *elem, double t,
                  double DLLPTR *loc, double DLLPTR *vel)
{
   t -= elem->perih_time;
   if( elem->ecc != 1.)    /* not parabolic */
      {
      t /= elem->t0;
      if( elem->ecc < 1.)     /* elliptical case;  throw out extra orbits */
         {                    /* to fit mean anom between -PI and PI */
         t = fmod( t, PI * 2.);
         if( t < -PI) t += 2. * PI;
         if( t >  PI) t -= 2. * PI;
         }
      elem->mean_anomaly = t;
      }
   comet_posn_part_ii( elem, t, loc, vel);
   return( 0);
}

int DLL_FUNC comet_posn( ELEMENTS DLLPTR *elem, double t, double DLLPTR *loc)
{
   return( comet_posn_and_vel( elem, t, loc, NULL));
}
