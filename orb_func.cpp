/* orb_func.cpp: basic orbital element/numerical integration funcs

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
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "lsquare.h"
#include "date.h"
#include "afuncs.h"
#include "monte0.h"

#ifndef _MSC_VER
   #define CONSOLE
#endif

int perturbers = 0;
int integration_method = 0;
extern int debug_level;

      /* In the console version of Find_Orb,  the following two functions */
      /* get remapped to Curses functions.  In the non-interactive one,   */
      /* they're mapped to 'do-nothings'.  See fo.cpp & find_orb.cpp.     */
void refresh_console( void);
void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes);

#define AUTOMATIC_PERTURBERS  1
#define J2000 2451545.
#define JD_TO_YEAR( jd)  (2000. + ((jd) - J2000) / 365.25)
#define MAX_CONSTRAINTS 5
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define SQRT_2 1.41421356

const double SRP1AU = 2.3e-7;
   /* "Solar radiation pressure at 1 AU",  in kg*AU^3 / (m^2*d^2),       */
   /* from a private communication from Steve Chesley.  This means       */
   /* that if you had a one-kilogram object showing one square           */
   /* meter of surface area to the sun,  and it was one AU from the      */
   /* sun,  and it absorbed all the solar radiation (and re-radiated     */
   /* it isotropically;  i.e.,  the re-radiated photons didn't cause     */
   /* any net force),  then that object would accelerate away from       */
   /* the sun at 2.3e-7 AU/day^2.                                        */
   /*                                                                    */
   /*   One can derive SRP1AU  from basic principles.   The 'solar       */
   /* constant' is C = 1367.6 AU^2*W/m^2; i.e.,  if you set up a         */
   /* one square meter solar panel with 100% efficiency one AU from      */
   /* the Sun, pointed straight at the sun, it would generate            */
   /* 1367.6 watts. To get the above number, one uses                    */
   /*                                                                    */
   /* SRP1AU = C * d^2 / (meters_per_AU * c)                             */
   /*                                                                    */
   /*    C = 1367.6 AU^2*W/m^2 = 1367.6 AU^2*kg/s^3                      */
   /*    d = 86400 seconds/day                                           */
   /*    meters_per_AU = 149597870691 m/AU                               */
   /*    c = 299792458 m/s                                               */
   /*                                                                    */
   /*    The result indicates that if the solar panel in question        */
   /* had a mass of one kilogram,  it would accelerate away from the     */
   /* sun at 2.27636e-7 AU/day^2.  I think Steve rounded off with a      */
   /* one-percent error because that's a decent match to the accuracy    */
   /* you can expect with real-world objects that reflect and re-radiate */
   /* photons,  instead of politely absorbing them and then re-radiating */
   /* them isotropically.                                                */

int n_extra_params = 0, setting_outside_of_arc = 1;
double solar_pressure[3], uncertainty_parameter = 99.;

int debug_printf( const char *format, ...);
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
             const double r1, const double angle_param);   /* orb_func.cpp */
int set_locs( const double *orbit, double t0, OBSERVE FAR *obs, int n_obs);
int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect);         /* runge.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
static int evaluate_limited_orbit( const double *orbit,
                    const int planet_orbiting, const double epoch,
                    const char *limited_orbit, double *constraints);
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
int find_relative_orbit( const double jd, const double *ivect,
               ELEMENTS *elements, const int ref_planet);
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
static inline void look_for_best_subarc( const OBSERVE FAR *obs,
       const int n_obs, const double max_arc_len, int *start, int *end);
int check_for_perturbers( const double t_cen, const double *vect); /* sm_vsop*/
int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                int *idx1, int *idx2);      /* elem_out.c */
void set_distance( OBSERVE FAR *obs, double r);             /* orb_func.c */
double find_r_given_solar_r( const OBSERVE FAR *obs, const double solar_r);

void set_distance( OBSERVE FAR *obs, double r)
{
   int i;

   obs->r = r;
   for( i = 0; i < 3; i++)
      obs->obj_posn[i] = obs->obs_posn[i] + r * obs->vect[i];
   obs->solar_r = vector3_length( obs->obj_posn);
}

double find_r_given_solar_r( const OBSERVE FAR *obs, const double solar_r)
{
   double r_dot_v = 0., r_dot_r = 0., b, c, discr, rval = -1.;
   int i;

   for( i = 0; i < 3; i++)
      {
      r_dot_r += obs->obs_posn[i] * obs->obs_posn[i];
      r_dot_v += obs->obs_posn[i] * obs->vect[i];
      }
   b = 2. * r_dot_v;
   c = r_dot_r - solar_r * solar_r;
   discr = b * b - 4 * c;
   if( discr > 0.)
      rval = (-b + sqrt( discr)) / 2.;
   return( rval);
}

int calc_derivatives( const double jd, const double *ival, double *oval,
                           const int reference_planet);
double take_rk_step( const double jd, ELEMENTS *ref_orbit,
                 const double *ival, double *ovals,
                 const int n_vals, const double step);      /* runge.cpp */
double take_pd89_step( const double jd, ELEMENTS *ref_orbit,
                 const double *ival, double *ovals,
                 const int n_vals, const double step);      /* runge.cpp */
void compute_ref_state( ELEMENTS *ref_orbit, double *ref_state,
                                          const double jd);
int symplectic_6( double jd, ELEMENTS *ref_orbit, double *vect,
                                          const double dt);
static int is_unreasonable_orbit( const double *orbit);     /* orb_func.cpp */

double integration_tolerance = 1.e-12;

char *runtime_message;
int show_runtime_messages = 1;

// static int reference_planet = 0;
static int perturbers_automatically_found;

static int reset_auto_perturbers( const double jd, const double *orbit)
{
   const int perturbing_planet = check_for_perturbers(
                                (jd - J2000) / 36525., orbit);

   perturbers = AUTOMATIC_PERTURBERS | (1 << perturbing_planet);
   if( perturbing_planet == 3)    /* add in the moon,  too: */
      perturbers |= (1 << 10);
   if( perturbing_planet == 10)    /* or vice versa:        */
      {
      debug_printf( "Moon is chief perturber\n");
      perturbers |= (1 << 3);
      }
   perturbers_automatically_found |= perturbers;
   return( perturbing_planet);
}


#define STEP_INCREMENT 2
#define ENCKE_CHANGE_LIMIT .01

int integrate_orbit( double *orbit, const double t0, const double t1)
{
   static double stepsize = 2.;
   const double min_stepsize = 1e-6;
   const double chicken = .9;
   int reset_of_elements_needed = 1;
   const double step_increase = chicken * integration_tolerance
                 / pow( STEP_INCREMENT, (integration_method ? 9. : 5.));
   const int use_cowell = atoi( get_environment_ptr( "COWELL"));
   double t = t0;
#ifdef CONSOLE
   static time_t real_time = (time_t)0;
   double prev_t = t, last_err = 0.;
#endif
   int n_rejects = 0, rval;
   int n_steps = 0, saved_perturbers = perturbers, prev_n_steps = 0;
   int going_backward = (t1 < t0);
   static int n_changes;
   ELEMENTS ref_orbit;

   rval = is_unreasonable_orbit( orbit);
   if( rval)
      {
      debug_printf( "Unreasonable %d\n", rval);
      return( -1);
      }
   ref_orbit.central_obj = -1;
   stepsize = fabs( stepsize);
   if( going_backward)
      stepsize = -stepsize;
   while( t != t1 && !rval)
      {
      double delta_t, new_t = ceil( (t - .5) / stepsize + .5) * stepsize + .5;

      if( perturbers & AUTOMATIC_PERTURBERS)
         reset_auto_perturbers( t, orbit);
      if( reset_of_elements_needed || !(n_steps % 50))
         if( use_cowell != 1)
            {
            extern int best_fit_planet;

            find_relative_orbit( t, orbit, &ref_orbit, best_fit_planet);
            reset_of_elements_needed = 0;
            if( !perturbers && !use_cowell)  /* for unperturbed cases,  do  */
               {                    /* a two-body soln and go home */
               double ref_state[9];

               compute_ref_state( &ref_orbit, ref_state, t1);
               rval = is_unreasonable_orbit( orbit);
               if( rval)
                  debug_printf( "Unreasonable %d after two-body soln\n",
                                 rval);
               memcpy( orbit, ref_state, 6 * sizeof( double));
               return( 0);
               }
            }
      n_steps++;
#ifdef CONSOLE
      if( time( NULL) != real_time && show_runtime_messages)
         {
         char buff[80];
         extern int best_fit_planet, n_posns_cached;
         extern double best_fit_planet_dist;
#define TEST_PLANET_CACHING_HASH_FUNCTION
#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
         extern long total_n_searches, total_n_probes, max_probes_required;
#endif

         if( runtime_message)
            move_add_nstr( 9, 10, runtime_message, -1);
         sprintf( buff, "t = %.5lf; %.5lf to %.5lf; step %.4e   ",
                 JD_TO_YEAR( t), JD_TO_YEAR( t0), JD_TO_YEAR( t1), stepsize);
         if( prev_n_steps)    /* i.e.,  not our first time through here */
            sprintf( buff + strlen( buff), "%d step/sec  ",
                        n_steps - prev_n_steps);
         move_add_nstr( 10, 10, buff, -1);
         prev_n_steps = n_steps;
         real_time = time( NULL);
         sprintf( buff, " %02ld:%02ld:%02ld; %lf; %d cached   ",
                     (real_time / 3600) % 24L,
                     (real_time / 60) % 60, real_time % 60, t - prev_t,
                     n_posns_cached);
         prev_t = t;
         move_add_nstr( 11, 10, buff, -1);
         sprintf( buff, "%d steps; %d rejected; center %d, ",
                           n_steps, n_rejects, best_fit_planet);
         format_dist_in_buff( buff + strlen( buff), best_fit_planet_dist);
         strcat( buff, "  ");
         move_add_nstr( 12, 10, buff, -1);
         sprintf( buff, "last err: %.3e/%.3e  n changes: %d  ",
                        last_err, step_increase, n_changes);
         move_add_nstr( 13, 10, buff, -1);
         if( use_cowell != 1)
            {
            sprintf( buff, "e = %.5lf; q = ", ref_orbit.ecc);
            format_dist_in_buff( buff + strlen( buff), ref_orbit.q);
            strcat( buff, "     ");
            move_add_nstr( 18, 10, buff, -1);
            }
         sprintf( buff, "Pos: %11.6lf %11.6lf %11.6lf",
                     orbit[0], orbit[1], orbit[2]);
         move_add_nstr( 14, 10, buff, -1);
         sprintf( buff, "Vel: %11.6lf %11.6lf %11.6lf",
                     orbit[3], orbit[4], orbit[5]);
         move_add_nstr( 15, 10, buff, -1);
#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
         if( total_n_searches)
            {
            sprintf( buff, "%ld searches; avg %.2lf max %ld     ",
                            total_n_searches,
                            (double)total_n_probes / (double)total_n_searches,
                            max_probes_required);
            move_add_nstr( 16, 10, buff, -1);
            }
#endif
         refresh_console( );
         }
#endif

               /* Make sure we don't step completely past */
               /* the time t1 we want to stop at!         */
      if( (!going_backward && new_t > t1) || (going_backward && new_t < t1))
         new_t = t1;
      delta_t = new_t - t;

      switch( integration_method)
         {
         case 1:
            symplectic_6( t, &ref_orbit, orbit, delta_t);
            break;
         default:
            {
            double new_vals[6];
            const double err = (integration_method ?
                   take_pd89_step( t, &ref_orbit, orbit, new_vals, 6, delta_t) :
                   take_rk_step( t, &ref_orbit, orbit, new_vals, 6, delta_t));

            if( !stepsize)
               exit( 0);
            if( err < integration_tolerance
                        || fabs( stepsize) < min_stepsize)  /* it's good! */
               {
               memcpy( orbit, new_vals, 6 * sizeof( double));
               if( err < step_increase)
                  if( fabs( delta_t - stepsize) < fabs( stepsize * .01))
                     {
                     n_changes++;
                     stepsize *= STEP_INCREMENT;
                     }
               }
            else           /* failed:  try again with a smaller step */
               {
               n_rejects++;
               new_t = t;
               stepsize /= STEP_INCREMENT;
               reset_of_elements_needed = 1;
               }
#ifdef CONSOLE
            last_err = err;
#endif
            }
            break;
         }
      t = new_t;
      rval = is_unreasonable_orbit( orbit);
      if( rval)
         {
         debug_printf( "Unreasonable %d at %.5lg (%.5lf to %.5lf)\n",
                rval, t - t0, t0, t1);
         debug_printf( "Stepsize %lg\n", stepsize);
         }
      }
   perturbers = saved_perturbers;
   return( rval);
}

/* At times,  the orbits generated by 'full steps' or Herget or other methods
   are completely unreasonable.  The exact definition of 'unreasonable'
   is pretty darn fuzzy.  The following function says that if at the epoch,
   an object is more than one light-year from the Sun (about 60000 AU) or
   is travelling at faster than 5% of the speed of light,  the orbit is
   "unreasonable".
      I'm sure this definition could be tightened a lot without any real
   trouble.  The fastest natural objects in the solar system are comets
   impacting the Sun,  which they do at the escape speed of about 600 km/s,
   or .2% of the speed of light.  (With the limit being a lot tighter than
   that as one gets away from the sun;  for example,  near the earth's
   orbit,  nothing goes much faster than about 70 km/s.)  And the furthest
   objects seen orbiting the sun are of the order of 100 AU away,  as
   opposed to the 60000 AU distance given below.   */

static const double max_reasonable_dist = 1.0 * 365.25 * AU_PER_DAY;
static const double max_reasonable_speed = AU_PER_DAY * .05;

static int is_unreasonable_orbit( const double *orbit)
{
   int rval;

   if( vector3_length( orbit) > max_reasonable_dist ||
       vector3_length( orbit + 3) > max_reasonable_speed)
      rval = -1;
   else
      rval = 0;
   return( rval);
}

static void light_time_lag( const double *orbit, const double *observer, double *result)
{
   double r = 0., delta, dt, solar_r2 = 0., afact;
   int i;

   for( i = 0; i < 3; i++)
      {
      delta = orbit[i] - observer[i];
      r += delta * delta;
      solar_r2 += orbit[i] * orbit[i];
      }
   r = sqrt( r);
   dt = -r / AU_PER_DAY;
   afact = 1. - SOLAR_GM * .5 * dt * dt / (sqrt( solar_r2) * solar_r2);
   for( i = 0; i < 3; i++)
      result[i] = afact * orbit[i] + dt * orbit[i + 3];
}

static void set_solar_r( OBSERVE FAR *ob)
{
   ob->solar_r = vector3_length( ob->obj_posn);
}

/* In setting the locations for each observation,  a little bit of
trickery is used to improve performance.  First (pass = 0),  we
set observations made _before_ t0,  integrating _backwards_.  Then
we set the observations made _after_ t0,  integrating _forwards_.
A small point,  but unfortunately,  it does make the code a little
less intuitively obvious. */

int set_locs( const double *orbit, const double t0, OBSERVE FAR *obs,
                       const int n_obs)
{
   int i, pass;

   for( i = 0; i < n_obs && obs[i].jd < t0; i++)
      ;

   for( pass = 0; pass < 2; pass++)
      {
      int j = (pass ? i - 1 : i);
      double curr_orbit[6], tvals[6];
      double curr_t = t0;

      memcpy( curr_orbit, orbit, 6 * sizeof( double));
      while( j < n_obs && j >= 0)
         {
         OBSERVE FAR *optr = obs + j;

         integrate_orbit( curr_orbit, curr_t, optr->jd);
         curr_t = optr->jd;
         FMEMCPY( tvals, optr->obs_posn, 3 * sizeof( double));
         light_time_lag( curr_orbit, tvals, tvals + 3);
         FMEMCPY( optr->obj_posn, tvals + 3, 3 * sizeof( double));
         FMEMCPY( optr->obj_vel, curr_orbit + 3, 3 * sizeof( double));
         j += (pass ? -1 : 1);
         }
      }

            /* We've now set the object heliocentric positions and */
            /* velocities,  in ecliptic J2000,  for each observation */
            /* time.  Now let's go back and find observer-centric */
            /* computed RA/decs and distances to the object at those */
            /* times. */
   for( i = 0; i < n_obs; i++)
      {
      double loc[3], ra, dec, temp, r = 0.;
      static const double sin_obliq_2000 = .397777156;
      static const double cos_obliq_2000 = .917482062;
      int j;

      for( j = 0; j < 3; j++)
         {
         loc[j] = obs[i].obj_posn[j] - obs[i].obs_posn[j];
         r += loc[j] * loc[j];
         }
      r = sqrt( r);
      obs[i].r = r;
      temp = loc[1] * cos_obliq_2000 - loc[2] * sin_obliq_2000;
      loc[2] = loc[2] * cos_obliq_2000 + loc[1] * sin_obliq_2000;
      loc[1] = temp;
      ra = atan2( loc[1], loc[0]);
      if( r > 100000. || r <= 0.)
         debug_printf( "???? bad r: %lf %lf %lf\n",
                  loc[0], loc[1], loc[2]);
      if( r)
         dec = asine( loc[2] / r);
      else
         dec = 0.;
      while( ra - obs[i].ra > PI)
         ra -= 2. * PI;
      while( ra - obs[i].ra < -PI)
         ra += 2. * PI;
      obs[i].computed_ra = ra;
      obs[i].computed_dec = dec;
      set_solar_r( obs + i);
      }
   return( 0);
}

/* 2010 Nov 4:  revised so that RMS residuals are computed as
root-mean-square of RA and dec treated separately,  meaning the
previous results needed to be multipled by sqrt(.5)  (Gareth
Williams kindly steered me the right way on this).          */

double compute_rms( const OBSERVE FAR *obs, const int n_obs)
{
   double rval = 0.;
   int i, n_included = 0;

   for( i = n_obs; i; i--, obs++)
      if( obs->is_included)
         {
         const double d_dec = obs->computed_dec - obs->dec;
         const double d_ra  = (obs->computed_ra  - obs->ra) * cos( obs->computed_dec);

         rval += d_dec * d_dec + d_ra * d_ra;
         n_included++;
         }
   rval = sqrt( rval / (double)(n_included * 2));
   return( 3600. * (180. / PI) * rval);
}

static double eval_3x3_determinant( const double *a, const double *b,
                         const double *c)
{
   return( a[0] * (b[1] * c[2] - c[1] * b[2])
         + b[0] * (c[1] * a[2] - a[1] * c[2])
         + c[0] * (a[1] * b[2] - b[1] * a[2]));
}

/* 'find_transfer_orbit' finds the state vector that can get an object
from the location/time described at 'obs1' to that described at 'obs2'.
It does this with the logic given in 'herget.htm#sund_xplns'. */
/* 24 Jun 2007:  realized 'accel0' wasn't used any more and removed it. */
/* 2011 May 8:  Realized this code is responsible for much of the time
   it takes to initially compute an orbit.  One simple speed-up is to
   take advantage of the fact that in the method of Herget,  we're
   often tweaking a "nearby" orbit,  i.e.,  we already have an
   approximate orbit;  using this as our starting estimate ought to
   result in faster convergence.       */

#define XFER_TOO_FAR_AWAY          -1
#define XFER_TOO_FAST              -2
#define XFER_TOO_MANY_ITERATIONS   -3
#define XFER_OK                     0

clock_t t_transfer;

static int find_transfer_orbit( double *orbit, OBSERVE FAR *obs1,
                OBSERVE FAR *obs2,
                const int already_have_approximate_orbit)
{
   double r = 0.;
   const double jd1 = obs1->jd - obs1->r / AU_PER_DAY;
   const double jd2 = obs2->jd - obs2->r / AU_PER_DAY;
   const double delta_t = jd2 - jd1;
   double deriv[6];
   double orbit2[6];
   double diff_squared = 999.;
            /* Iterate until the error is less than 1e-8 of the object */
            /* distance.  This corresponds to an error of about .002 arcsec */
            /* in the actual position (with the '+ .01' allowing for */
            /* some margin for very close objects,  such as artsats). */
            /*   Assume a maximum of ten iterations,  just to be safe */
   const double target_diff = 1.e-8 * (obs2->r + .01);
   int i, max_iterations = 10;
   clock_t t0 = clock( );

   if( obs1->r > max_reasonable_dist || obs2->r > max_reasonable_dist)
      {
      debug_printf( "Bad xfer: %lf %lf\n", obs1->r, obs2->r);
      return( XFER_TOO_FAR_AWAY);
      }

   set_distance( obs1, obs1->r);
   set_distance( obs2, obs2->r);

   if( !already_have_approximate_orbit)
      {
      double speed_squared = 0., speed;
      const int saved_perturbers = perturbers;
      const double mid_jd = (jd1 + jd2) / 2.;

      for( i = 0; i < 3; i++)
         {
                     /* find midpoint between obs1 and obs2... */
         orbit[i] = (obs2->obj_posn[i] + obs1->obj_posn[i]) / 2.;
                     /* and the average velocity connecting them: */
         orbit[i + 3] = (obs2->obj_posn[i] - obs1->obj_posn[i]) / delta_t;
         speed_squared += orbit[i + 3] * orbit[i + 3];
         }

      if( speed_squared > max_reasonable_speed * max_reasonable_speed)
         return( XFER_TOO_FAST);
                     /* speed_squared is in (AU/day)^2, speed in km/second: */
      speed = sqrt( speed_squared) * AU_IN_KM / 86400.;
      if( debug_level > 3)
         debug_printf( "Speed = %lf km/s\n", speed);
      if( perturbers & AUTOMATIC_PERTURBERS)
         reset_auto_perturbers( mid_jd, orbit);
      calc_derivatives( mid_jd, orbit, deriv, -1);
      perturbers = saved_perturbers;
      for( i = 0; i < 3; i++)
         orbit[3 + i] -= deriv[i + 3] * delta_t / 2.;
      }
   for( i = 0; i < 3; i++)
      orbit[i] = obs1->obj_posn[i];

   for( i = 0; i < 3 && diff_squared > target_diff * target_diff; i++)
      {
      int j;

      memcpy( orbit2, orbit, 6 * sizeof( double));
      integrate_orbit( orbit2, jd1, jd2);
      diff_squared = 0;
      for( j = 0; j < 3; j++)
         {
         const double delta = orbit2[j] - obs2->obj_posn[j];

         orbit[j + 3] -=  delta / delta_t;
         diff_squared += delta * delta;
         }
      if( debug_level > 3)
         debug_printf( "Initial pass %d: %.3lg\n", i, diff_squared);
      }

   while( diff_squared > target_diff * target_diff && max_iterations--)
      {
      double delta[4][3], discr;
      int pass;

      for( pass = 0; pass < 4; pass++)
         {
         const double h = 1.e-5;

         memcpy( orbit2, orbit, 6 * sizeof( double));
         if( pass)
            orbit2[(pass - 1) + 3] += h;
         integrate_orbit( orbit2, jd1, jd2);
         for( i = 0; i < 3; i++)
            orbit2[i] -= obs2->obj_posn[i];
         memcpy( delta[pass], orbit2, 3 * sizeof( double));
         if( pass)
            for( i = 0; i < 3; i++)
               delta[pass][i] = (delta[pass][i] - delta[0][i]) / h;
         }
      discr = eval_3x3_determinant( delta[1], delta[2], delta[3]);
      if( discr)
         {
         orbit[3] -= eval_3x3_determinant( delta[0], delta[2], delta[3]) / discr;
         orbit[4] -= eval_3x3_determinant( delta[1], delta[0], delta[3]) / discr;
         orbit[5] -= eval_3x3_determinant( delta[1], delta[2], delta[0]) / discr;
         }
      diff_squared = 0.;
      for( i = 0; i < 3; i++)
         diff_squared += delta[0][i] * delta[0][i];
      if( debug_level > 3)
         debug_printf( "Transfer orbit: %d, %.3lg; discr=%lg\n", max_iterations,
                        diff_squared, discr);
      }
   if( max_iterations <= 0)
      return( XFER_TOO_MANY_ITERATIONS);
               /* adjust for light-time lag: */
   r = -SOLAR_GM / (obs1->solar_r * obs1->solar_r * obs1->solar_r);
   for( i = 0; i < 3; i++)
      {
      const double time_lag = obs1->r / AU_PER_DAY;

      orbit[i] += time_lag * (orbit[3 + i] + orbit[i] * r * time_lag / 2.);
      orbit[i + 3] += orbit[i] * r * time_lag;
      }
   t_transfer += clock( ) - t0;
   return( XFER_OK);
}

int drop_excluded_obs( OBSERVE *obs, int *n_obs)
{
   int rval = 0;

   while( *n_obs && !obs->is_included)    /* skip excluded obs at */
      {                                   /* start of arc */
      rval++;
      obs++;
      (*n_obs)--;
      }
   while( *n_obs && !obs[*n_obs - 1].is_included)   /* skip excluded obs at */
      (*n_obs)--;                                   /* end of arc */
   return( rval);
}

/* 'find_trial_orbit' tries to find an orbit linking the first and last
(included) observations in 'obs',  with the distance to the first object
being r1 AU and a radial velocity defined by 'angle_param'.

   First,  excluded observations at the beginning and end are skipped over
by advancing the obs pointer and/or decrementing n_obs.  If there are at
least two observations remaining,  we set the first observation to be at
distance r1.  Next,  we compute the distance in space between that point
and the ray defined by the second observation,  and we also compute the
maximum distance the object could go (assuming a parabolic orbit with
the object at escape velocity) during the time 'dt' between the two
observations.

   If the object can't get to the ray from the first observation without
going faster than escape velocity,  then we return an error code of -2.
This basically says,  "There are no non-hyperbolic orbits that satisfy
these two observations with the value you gave for r1.  Sorry."

   If the object _can_ get there,  we compute an orbit with an assumed
radial velocity linking the two observations.  If angle_param == -1,
the orbit will be one in which the object is at escape speed,  going
_toward_ us.  If angle_param == 1,  it'll be at escape speed,  going
_away_ from us.  In between,  you'll get assorted elliptical orbits.
(If angle_param > 1 or less than -1,  you'll get an hyperbolic orbit.)
*/

int find_sr_ranges( double *ranges, const double *q1, const double *p1,
                                    const double *q2, const double *p2,
                                    const double gm, const double dt);

int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                 const double r1, const double angle_param)
{
   int rval = 0;

   uncertainty_parameter = 99.;
   obs += drop_excluded_obs( obs, &n_obs);
   if( n_obs < 2)       /* need at least two valid obs */
      rval = -1;
   else
      {
      double r2 = 0., dist2 = 0., escape_dist2;
      OBSERVE FAR *endptr = obs + n_obs - 1;
      const double dt = endptr->jd - obs->jd;
      int i;

      set_distance( obs, r1);
      for( i = 0; i < 3; i++)
         r2 += (obs->obj_posn[i] - endptr->obs_posn[i]) * endptr->vect[i];
      set_distance( endptr, r2);
      for( i = 0; i < 3; i++)
         {
         const double delta = endptr->obj_posn[i] - obs->obj_posn[i];

         dist2 += delta * delta;
         }
               /* Escape velocity for an object r AU from the sun would be */
               /* sqrt( 2 * SOLAR_GM / r).  We'll use the square of this:  */
      escape_dist2 = 2. * SOLAR_GM / obs->solar_r;
               /* ...and multiply by dt squared to get the square of the   */
               /* distance the object would travel if it's at the esc speed: */
      escape_dist2 *= dt * dt;
      if( escape_dist2 < dist2)     /* only hyperbolic orbits exist */
         rval = -2;
      else
         {
         set_distance( endptr, r2 + angle_param * sqrt( escape_dist2 - dist2));
         find_transfer_orbit( orbit, obs, endptr, 0);
         set_locs( orbit, obs->jd, obs, n_obs);
         if( n_obs > 2)
            adjust_herget_results( obs, n_obs, orbit);
         }
      }
   return( rval);
}

static double haltonize( int idx, const int radix)
{
   int divisor = 1, rval = 0;

   while( idx)
      {
      rval = (rval * radix) + idx % radix;
      idx /= radix;
      divisor *= radix;
      }
   return( (double)rval / (double)divisor);
}

double sr_min_r = 0., sr_max_r = 0.;

int find_nth_sr_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                            const int orbit_number)
{
   int rval = 0;

   uncertainty_parameter = 99.;
   obs += drop_excluded_obs( obs, &n_obs);
   if( n_obs < 2)       /* need at least two valid obs */
      {
      debug_printf( "SR fail: n_obs = %d\n", n_obs);
      rval = -1;
      }
   else
      {
      const double rand1 = haltonize( orbit_number + 1, 2);
      const double rand2 = haltonize( orbit_number + 1, 3);
      double dist = 0.;
      int i;
      static int n_roots;
      static double roots[10], total_dist;

      if( !orbit_number)
         {
         if( sr_max_r)        /* user override of SR ranges */
            {
            roots[0] = sr_min_r;
            roots[1] = sr_max_r;
            n_roots = 1;
            }
         else        /* "normal" determination of possible ranges */
            n_roots = find_sr_ranges( roots,
                        obs[n_obs - 1].obs_posn, obs[n_obs - 1].vect,
                        obs[    0    ].obs_posn, obs[    0    ].vect,
                        SOLAR_GM, obs->jd - obs[n_obs - 1].jd);
         if( debug_level > 1)
            {
            debug_printf( "%d ranges\n", n_roots);
            for( i = 0; i < n_roots * 2; i++)
               debug_printf( "Root %d: %lf\n", i, roots[i]);
            }
         total_dist = 0.;
         for( i = 0; i < n_roots; i += 2)
            total_dist += roots[i + 1] - roots[i];
         if( debug_level > 1)
            debug_printf( "Total dist: %lf\n", total_dist);
         }
      for( i = 0; i < n_roots; i += 2)
         {
         const double d2 = dist + roots[i + 1] - roots[i];

         if( d2 < total_dist * rand1)
            dist = d2;
         else
            {
            dist = roots[i] + (total_dist * rand1 - dist);
            i = n_roots;      /* break out of loop */
            }
         }
      find_trial_orbit( orbit, obs, n_obs, dist, 2. * rand2 - 1.);
      }
   return( rval);
}

static OBSERVE *get_real_arc( OBSERVE *obs, int *n_obs,
                              int *n_real_obs)
{
   int idx1, idx2;

   *n_real_obs = get_idx1_and_idx2( *n_obs, obs, &idx1, &idx2);
   *n_obs = idx2 - idx1 + 1;
   return( obs + idx1);
}

#define RAD2SEC (180. * 3600. / PI)

void set_obs_vect( OBSERVE FAR *obs);        /* mpc_obs.h */

int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit)
{
   int i, n_found, rval = -1;
   double avg_x2 = 0., avg_x = 0., avg_xt = 0., avg_t = 0.;
   double avg_y2 = 0., avg_y = 0., avg_yt = 0., avg_t2 = 0.;

   obs = get_real_arc( obs, &n_obs, &n_found);
   if( n_found < 3)   /* must have at least three obs */
      return( -2);
   if( is_unreasonable_orbit( orbit))
      return( -3);
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         const double dx = obs[i].computed_ra - obs[i].ra;
         const double dy = obs[i].computed_dec - obs[i].dec;
         const double dt = obs[i].jd - obs[0].jd;

         avg_x2 += dx * dx;
         avg_y2 += dy * dy;
         avg_xt += dx * dt;
         avg_yt += dy * dt;
         avg_x += dx;
         avg_y += dy;
         avg_t += dt;
         avg_t2 += dt * dt;
         }

// if( n_found)   /* we already know n_found >= 3 */
      {
      double determ;

      avg_x2 /= (double)n_found;
      avg_y2 /= (double)n_found;
      avg_xt /= (double)n_found;
      avg_yt /= (double)n_found;
      avg_x /= (double)n_found;
      avg_y /= (double)n_found;
      avg_t /= (double)n_found;
      avg_t2 /= (double)n_found;
      determ = avg_t2 - avg_t * avg_t;
      if( determ)
         {
         double ax = (avg_xt - avg_x * avg_t) / determ;
         double bx = (avg_t2 * avg_x - avg_t * avg_xt) / determ;
         double ay = (avg_yt - avg_y * avg_t) / determ;
         double by = (avg_t2 * avg_y - avg_t * avg_yt) / determ;
         OBSERVE start = obs[0];
         OBSERVE end = obs[n_obs - 1];

         start.ra  = start.computed_ra - bx;
         end.ra    = end.computed_ra - (bx + (end.jd - start.jd) * ax);
         start.dec = start.computed_dec - by;
         end.dec   = end.computed_dec - (by + (end.jd - start.jd) * ay);
         set_obs_vect( &start);
         set_obs_vect( &end);
         set_distance( &start, start.r);
         set_distance( &end, end.r);
         uncertainty_parameter = 99.;
         find_transfer_orbit( orbit, &start, &end, 1);
         set_locs( orbit, start.jd, obs, n_obs);
         rval = 0;
         }
      }
   return( rval);
}

int herget_method( OBSERVE FAR *obs, int n_obs, double r1, double r2,
         double *orbit, double *d_r1, double *d_r2, const char *limited_orbit)
{
   double *xresid, *yresid, *dx_dr1, *dx_dr2, *dy_dr1, *dy_dr2;
   double delta, a = 0., b = 0., c = 0., e = 0., f = 0., determ;
   int i, n_real_obs;
   const int using_pseudo_vaisala = (r1 < 0.);
   double orbit2[6];
   double orbit_offset[6], *constraint = NULL;
   int planet_orbiting = 0, n_constraints = 0;
   char tstr[80];
   OBSERVE temp_obs1, temp_obs2;

   if( limited_orbit && !*limited_orbit)
      limited_orbit = NULL;

   obs = get_real_arc( obs, &n_obs, &n_real_obs);
   if( n_real_obs < 2)        /* should never happen */
      {
      debug_printf( "??? n_obs %d, n_real_obs = %d\n",
                     n_obs, n_real_obs);
      return( -1);
      }

   temp_obs1 = obs[0];
            /* Look "ahead" up to 100 days: */
   for( i = n_obs - 1; i && obs[i].jd - obs[0].jd > 100.; i--)
      ;
   temp_obs2 = obs[i];
   uncertainty_parameter = 99.;
   if( using_pseudo_vaisala)
      {
      r2 = find_r_given_solar_r( &temp_obs2, -r1);
      r1 = find_r_given_solar_r( &temp_obs1, -r1);

      if( debug_level > 7)
         debug_printf( "r1 = %lf; r2 = %lf\n",
                  r1, r2);
      if( r1 < 0. || r2 < 0.)
         return( 1);
      if( d_r1)
         *d_r1 = r1;
      if( d_r2)
         *d_r2 = r2;
      }
   set_distance( &temp_obs1, r1);
   set_distance( &temp_obs2, r2);
   runtime_message = tstr;
   if( using_pseudo_vaisala)
      sprintf( tstr, "Vaisala %lf\n", obs->solar_r);
   else
      strcpy( tstr, "H/xfer orbit (1)");
               /* Compute the trial orbit in the local orbit2 array.  That */
               /* way,  if we find it's completely stupid,  we've not      */
               /* done anything to the plain old 'orbit' vector,  which    */
               /* may already contain a valid state vector:                */
   if( find_transfer_orbit( orbit2, &temp_obs1, &temp_obs2, 0)
                   || is_unreasonable_orbit( orbit2))
      {
      runtime_message = NULL;
      return( -2);
      }
               /* But now that we know it's a good result,  let's copy:     */
   memcpy( orbit, orbit2, 6 * sizeof( double));
   strcpy( tstr, using_pseudo_vaisala ? "Vaisala set_locs" : "H/set_locs (1)");
   set_locs( orbit, temp_obs1.jd, obs, n_obs);
   runtime_message = NULL;
   if( !d_r1 || !d_r2 || using_pseudo_vaisala)
      return( 0);

   *d_r1 = *d_r2 = 0.;
   if( n_real_obs == 2)
      return( 0);       /* good as done... */

   xresid = (double *)calloc( 6 * (n_obs + MAX_CONSTRAINTS), sizeof( double));
   if( !xresid)
      return( -3);
   yresid = xresid + n_obs + MAX_CONSTRAINTS;
   dx_dr1 = yresid + n_obs + MAX_CONSTRAINTS;
   dx_dr2 = dx_dr1 + n_obs + MAX_CONSTRAINTS;
   dy_dr1 = dx_dr2 + n_obs + MAX_CONSTRAINTS;
   dy_dr2 = dy_dr1 + n_obs + MAX_CONSTRAINTS;

   if( limited_orbit)
      {
      planet_orbiting = find_best_fit_planet( temp_obs1.jd, orbit, orbit2);
      for( i = 0; i < 6; i++)
         orbit_offset[i] = orbit[i] - orbit2[i];
      constraint = xresid + n_obs;
      n_constraints = evaluate_limited_orbit( orbit2, planet_orbiting,
                                       temp_obs1.jd, limited_orbit, constraint);
      }

   for( i = 0; i < n_obs; i++)
      {
      xresid[i] = obs[i].computed_ra - obs[i].ra;
      yresid[i] = obs[i].computed_dec - obs[i].dec;
      }

   delta = r1 / 10000.;
   if( delta > .1) delta = .1;
   set_distance( &temp_obs1, r1 + delta);
   runtime_message = tstr;
   strcpy( tstr, "H/xfer orbit (2)");
   memcpy( orbit2, orbit, 6 * sizeof( double));
   find_transfer_orbit( orbit2, &temp_obs1, &temp_obs2, 1);
   strcpy( tstr, "H/set_locs (2)");
   set_locs( orbit2, temp_obs1.jd, obs, n_obs);
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double dx = obs[i].computed_ra - obs[i].ra;
         double dy = obs[i].computed_dec - obs[i].dec;

         dx_dr1[i] = (dx - xresid[i]) / delta;
         dy_dr1[i] = (dy - yresid[i]) / delta;
         }

   if( limited_orbit)
      {
      double constraint2[MAX_CONSTRAINTS];

      for( i = 0; i < 6; i++)
         orbit2[i] -= orbit_offset[i];
      evaluate_limited_orbit( orbit2, planet_orbiting, temp_obs1.jd,
                                limited_orbit, constraint2);
      for( i = 0; i < n_constraints; i++)
         {
         dx_dr1[n_obs + i] = (constraint2[i] - constraint[i]) / delta;
         dy_dr1[n_obs + i] = 0.;
         }
      }

   delta = r2 / 10000.;
   if( delta > .1) delta = .1;
   set_distance( &temp_obs1, r1);
   set_distance( &temp_obs2, r2 + delta);
   strcpy( tstr, "H/xfer orbit (3)");
   memcpy( orbit2, orbit, 6 * sizeof( double));
   find_transfer_orbit( orbit2, &temp_obs1, &temp_obs2, 1);
   strcpy( tstr, "H/set_locs (3)");
   set_locs( orbit2, temp_obs1.jd, obs, n_obs);
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double dx = obs[i].computed_ra - obs[i].ra;
         double dy = obs[i].computed_dec - obs[i].dec;

         dx_dr2[i] = (dx - xresid[i]) / delta;
         dy_dr2[i] = (dy - yresid[i]) / delta;
         }

   if( limited_orbit)
      {
      double constraint2[MAX_CONSTRAINTS];

      for( i = 0; i < 6; i++)
         orbit2[i] -= orbit_offset[i];
      evaluate_limited_orbit( orbit2, planet_orbiting, temp_obs1.jd,
                                limited_orbit, constraint2);
      for( i = 0; i < n_constraints; i++)
         {
         dx_dr2[n_obs + i] = (constraint2[i] - constraint[i]) / delta;
         dy_dr2[n_obs + i] = 0.;
         }
      }

                    /* OK,  we now have all values & derivatives needed... */
   for( i = 0; i < n_obs + n_constraints; i++)
      if( obs[i].is_included || i >= n_obs)
         {
         c += xresid[i] * dx_dr1[i] + yresid[i] * dy_dr1[i];
         a += dx_dr1[i] * dx_dr1[i] + dy_dr1[i] * dy_dr1[i];
         b += dx_dr1[i] * dx_dr2[i] + dy_dr1[i] * dy_dr2[i];

         f += xresid[i] * dx_dr2[i] + yresid[i] * dy_dr2[i];
         /*  d = b;  */
         e += dx_dr2[i] * dx_dr2[i] + dy_dr2[i] * dy_dr2[i];
         }

   free( xresid);
   determ = a * e - b * b;
   *d_r1 = -(e * c - b * f) / determ;
   *d_r2 = -(a * f - c * b) / determ;
   if( *d_r1 > r1 / 3.) *d_r1 = r1 / 3.;
   if( *d_r1 <-r1 / 3.) *d_r1 = -r1 / 3.;
   if( *d_r2 > r2 / 3.) *d_r2 = r2 / 3.;
   if( *d_r2 <-r2 / 3.) *d_r2 = -r2 / 3.;
   runtime_message = NULL;
   return( 0);
}

static int setup_parabolic( const double *iparams, double *orbit)
{
   double r = 0., escape_vel;
   int i;

   for( i = 0; i < 3; i++)
      {
      orbit[i] = iparams[i];
      r += iparams[i] * iparams[i];
      }
   r = sqrt( r);
   escape_vel = SQRT_2 * GAUSS_K / sqrt( r);
   orbit[3] = escape_vel * cos( iparams[3]) * cos( iparams[4]);
   orbit[4] = escape_vel * sin( iparams[3]) * cos( iparams[4]);
   orbit[5] = escape_vel *                    sin( iparams[4]);
   return( 0);
}

void improve_parabolic( OBSERVE FAR *obs, int n_obs, double *orbit,
                                                              double epoch)
{
   void *lsquare = lsquare_init( 5);
   double *xresids = (double *)calloc( 2 * n_obs + 10 * n_obs, sizeof( double));
   double *yresids = xresids + n_obs;
   double *slopes = yresids + n_obs;
   double params2[5], params[5], differences[5];
   const double delta_val = .0001;
   double v2 = 0.;
   int i, j;

   uncertainty_parameter = 99.;
   memcpy( params, orbit, 3 * sizeof( double));
   params[3] = atan2( orbit[4], orbit[3]);
   for( i = 3; i < 6; i++)
      v2 += orbit[i] * orbit[i];
   params[4] = asin( orbit[5] / sqrt( v2));
   setup_parabolic( params, orbit);

   set_locs( orbit, epoch, obs, n_obs);
   for( i = 0; i < n_obs; i++)
      {
      xresids[i] = obs[i].computed_ra - obs[i].ra;
      yresids[i] = obs[i].computed_dec - obs[i].dec;
      }

   for( i = 0; i < 5; i++)
      {
      memcpy( params2, params, 5 * sizeof( double));
      params2[i] = params[i] - delta_val;
      setup_parabolic( params2, orbit);
      set_locs( orbit, epoch, obs, n_obs);
      for( j = 0; j < n_obs; j++)
         {
         slopes[i + j * 10] = obs[j].computed_ra;
         slopes[i + j * 10 + 5] = obs[j].computed_dec;
         }

      params2[i] = params[i] + delta_val;
      setup_parabolic( params2, orbit);
      set_locs( orbit, epoch, obs, n_obs);
      for( j = 0; j < n_obs; j++)
         {
         slopes[i + j * 10] -= obs[j].computed_ra;
         slopes[i + j * 10 + 5] -= obs[j].computed_dec;
         slopes[i + j * 10] /= 2. * delta_val;
         slopes[i + j * 10 + 5] /= 2. * delta_val;
         }
      }

   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         lsquare_add_observation( lsquare, xresids[i], 1., slopes + i * 10);
         lsquare_add_observation( lsquare, yresids[i], 1., slopes + i * 10 + 5);
         }

   free( xresids);
   lsquare_solve( lsquare, differences);
   lsquare_free( lsquare);

   for( i = 0; i < 5; i++)
      params[i] += differences[i];
   setup_parabolic( params, orbit);
   set_locs( orbit, epoch, obs, n_obs);
}

static int evaluate_limited_orbit( const double *orbit,
                    const int planet_orbiting, const double epoch,
                    const char *limited_orbit, double *constraints)
{
   int rval = 0;

   if( limited_orbit)
      {
      ELEMENTS elem;
      extern double planet_mass[];

      calc_classical_elements( &elem, orbit, epoch, 1,
                              SOLAR_GM * planet_mass[planet_orbiting]);
      while( *limited_orbit && limited_orbit[1] == '=')
         {
         double value = atof( limited_orbit + 2);
         char variable = *limited_orbit;

         while( *limited_orbit && *limited_orbit != ',')
            limited_orbit++;
         if( limited_orbit[-1] == 'k')
            value /= AU_IN_KM;
         switch( variable)
            {
            case 'q':
               constraints[rval] = elem.major_axis * (1. - elem.ecc) - value;
               rval++;
               break;
            case 'Q':
               constraints[rval] = elem.major_axis * (1. + elem.ecc) - value;
               rval++;
               break;
            case 'e':
               if( value)
                  constraints[rval++] = elem.ecc - value;
               else        /* handle e=0 (circular) orbits separately: */
                  {
                  constraints[rval++] = orbit[0] * orbit[3] +
                              orbit[1] * orbit[4] + orbit[2] * orbit[5];
                  constraints[rval++] = vector3_length( orbit)
                               - elem.major_axis;
                  }
               break;
            case 'i':
               constraints[rval++] = elem.incl * 180. / PI - value;
               break;
            case 'P':         /* convert to a major axis */
               if( limited_orbit[-1] == 'd')    /* convert from days to yrs */
                  value /= 365.25;
               if( limited_orbit[-1] == 'h')    /* convert from hrs to yrs */
                  value /= 365.25 * 24.;
               if( limited_orbit[-1] == 'm')    /* convert from mins to yrs */
                  value /= 365.25 * 1440.;
               value = pow( value * sqrt( planet_mass[planet_orbiting]),
                                                     2. / 3.);
                           /* then _don't_ break; just fall through to the */
                           /* major axis code:                             */
            case 'a':
               constraints[rval++] = 1. / elem.major_axis - 1. / value;
               break;
            case 'n':
               value = 360. / value;         /* now value = period in days */
               value /= 365.25;              /* now value = period in years */
               value = pow( value * sqrt( planet_mass[planet_orbiting]),
                                                     2. / 3.);
               constraints[rval++] = 1. / elem.major_axis - 1. / value;
               break;
            case 'A':            /* area/mass ratio */
               if( n_extra_params == 1)
                  constraints[rval++] =
                             solar_pressure[0] * SOLAR_GM / SRP1AU - value;
               break;
            case 'K':
               {
               extern double comet_magnitude_slope_param;

               comet_magnitude_slope_param = value;
               }
               break;
            case 'G':
               {
               extern double asteroid_magnitude_slope_param;

               asteroid_magnitude_slope_param = value;
               }
               break;
            }
         if( *limited_orbit == ',')
            limited_orbit++;
         }
      }
   return( rval);
}

#define MAX_N_PARAMS 9

double one_sigma_eigenvect[MAX_N_PARAMS];
void jacobi_eigenvalues( double *a, const int size, double *eigenvals,
                        double *eigenvects);       /* eigen.cpp */

void **calloc_double_dimension_array( const int x, const int y, size_t obj_size)
{
   void **rval = (void **)calloc( x * sizeof( void *) + x * y * obj_size, 1);
   int i;

   rval[0] = (void *)(rval + x);
   for( i = 1; i < x; i++)
      rval[i] = (void *)( (char *)rval[0] + i * y * obj_size);
   return( rval);
}

static double dot_prod( const double *a, const double *b)
{
   return( a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static double dotted_dist( OBSERVE FAR *obs)
{
   return( dot_prod( obs->vect, obs->obj_posn) - dot_prod( obs->vect, obs->obs_posn));
}

/* Describing what 'full_improvement()' does requires an entire separate
file of commentary: see 'full.txt'. */

int full_improvement( OBSERVE FAR *obs, int n_obs, double *orbit,
                 const double epoch, const char *limited_orbit,
                 int sigmas_requested)
{
   const int n_params = (!limited_orbit || memcmp( limited_orbit, "np=", 3)
                ? 6 + n_extra_params : atoi( limited_orbit + 3));
   void *lsquare = lsquare_init( n_params);
   double FAR *xresids;
   double FAR *yresids;
   double FAR *slopes;
   double constraint_slope[MAX_CONSTRAINTS * MAX_N_PARAMS];
   double element_slopes[MAX_N_PARAMS][MONTE_N_ENTRIES];
   double elements_in_array[MONTE_N_ENTRIES];
   double orbit2[6], differences[10];
   double original_orbit[6], original_params[3];
   double orbit_offset[6], working_epoch = epoch;
   double orbit_at_working_epoch[6];
   static double **unit_vectors = NULL;
   double **new_unit_vectors = NULL;
   double constraint[MAX_CONSTRAINTS];
   double sigma_squared = 0.;       /* see Danby, p. 243, (7.5.20) */
   double scale_factor = 1.;
#ifdef REMOVE_HALF_STEPS
   double before_rms, after_rms;
#endif
   int planet_orbiting, n_constraints = 0;
   int i, j, n_skipped_obs = 0, err_code = 0;
   int n_included_observations = 0;
   const int n_total_obs = n_obs;
   const char *covariance_filename = "covar.txt";
   char tstr[80];
   ELEMENTS elem;
   const int use_symmetric_derivatives =
                      atoi( get_environment_ptr( "SYMMETRIC"));
   const int showing_deltas_in_debug_file =
                      atoi( get_environment_ptr( "DEBUG_DELTAS"));
   const double r_mult = 1e+2;
   extern double planet_mass[];
   const double sigma_override = atof( get_environment_ptr( "NO_SIGMAS"));

   if( !unit_vectors)         /* initialize to an identity matrix */
      {
      unit_vectors = (double **)calloc_double_dimension_array(
                     MAX_N_PARAMS, MAX_N_PARAMS, sizeof( double));
      for( i = 0; i < MAX_N_PARAMS; i++)
         unit_vectors[i][i] = 1.;
      }
   if( get_idx1_and_idx2( n_obs, obs, &i, &j) < 3)
      return( -1);
   if( is_unreasonable_orbit( orbit))
      return( -2);
               /* We save the input orbit;  if there's an error,  we can */
               /* restore it:         */
   memcpy( original_orbit, orbit, 6 * sizeof( double));
   memcpy( original_params, solar_pressure, 3 * sizeof( double));
   sprintf( tstr, "full improvement: %lf  ", JD_TO_YEAR( epoch));
   runtime_message = tstr;
   for( i = 0; i < n_obs; i++)
      {
      obs[i].computed_ra  = obs[i].ra;
      obs[i].computed_dec = obs[i].dec;
      }

               /* Drop unincluded observations from the start of the arc: */
   while( n_obs && !obs->is_included)
      {
      obs++;
      n_obs--;
      n_skipped_obs++;
      }
               /* ... and then from the end of the arc: */
   while( n_obs && !obs[n_obs - 1].is_included)
      n_obs--;
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         n_included_observations++;
   if( n_included_observations < 3)
      return( -3);
   if( limited_orbit && !*limited_orbit)
      limited_orbit = NULL;
   if( n_params < 6)
      limited_orbit = NULL;

                  /* Computing sigmas/covariances with an epoch well outside */
                  /* the observed arc can be troublesome.  It's slower and   */
                  /* can cause numerical instability.  Setting a suitable    */
                  /* 'sigma_override' value to tell Find_Orb not to compute  */
                  /* covariances if the epoch is more than a certain number  */
                  /* of days outside the arc can help to avoid all that.     */
   if( sigma_override)
      if( epoch < obs[0].jd - sigma_override ||
                 epoch > obs[n_obs - 1].jd + sigma_override)
         sigmas_requested = NO_ORBIT_SIGMAS_REQUESTED;

   if( sigmas_requested == HELIOCENTRIC_SIGMAS_ONLY)
      {                          /* sigmas computed relative to the sun */
      planet_orbiting = 0;       /* even if the object is orbiting a planet */
      for( i = 0; i < 6; i++)
        orbit2[i] = orbit[i];
      }
   else
      planet_orbiting = find_best_fit_planet( epoch, orbit, orbit2);

   for( i = 0; i < 6; i++)
      orbit_offset[i] = orbit[i] - orbit2[i];

   if( limited_orbit)
      n_constraints = evaluate_limited_orbit( orbit2, planet_orbiting,
                                       epoch, limited_orbit, constraint);

               /* evaluate elements of 'orbit2',  then put */
               /* into an array form: */
   calc_classical_elements( &elem, orbit2, epoch, 1,
                              SOLAR_GM * planet_mass[planet_orbiting]);
   put_orbital_elements_in_array_form( &elem, elements_in_array);

      /* If the epoch is outside the actual arc of observations, we'll do a
      lot of integrating between the epoch and the actual arc of data. We
      can speed things up by using a "working epoch", placed within the
      actual arc.  Except we can't do that if we're after covariance info. */

   uncertainty_parameter = 99.;
   if( sigmas_requested == NO_ORBIT_SIGMAS_REQUESTED)
      {
      if( working_epoch < obs[0].jd)
         working_epoch = obs[0].jd;
      if( working_epoch > obs[n_obs - 1].jd)
         working_epoch = obs[n_obs - 1].jd;
      integrate_orbit( orbit, epoch, working_epoch);
      }
   memcpy( orbit_at_working_epoch, orbit, 6 * sizeof( double));

   xresids = (double FAR *)FCALLOC( (2 + 2 * n_params) * n_obs + n_params, sizeof( double));
   yresids = xresids + n_obs;
   slopes = yresids + n_obs;

   sprintf( tstr, "fi/setting locs: %lf  ", JD_TO_YEAR( working_epoch));
   set_locs( orbit, working_epoch, obs, n_obs);
#ifdef REMOVE_HALF_STEPS
   before_rms = compute_rms( obs, n_obs);
#endif
   if( limited_orbit && *limited_orbit == 'R')
      constraint[n_constraints++] =
              r_mult * (dotted_dist( obs + n_obs - 1) - atof( limited_orbit + 2));


   sprintf( tstr, "fi/locs set: %lf  ", JD_TO_YEAR( working_epoch));
   for( i = 0; i < n_obs; i++)
      {
      xresids[i] = obs[i].computed_ra  - obs[i].ra;
      yresids[i] = obs[i].computed_dec - obs[i].dec;
      }

   for( i = 0; i < n_params; i++)
      {
      static double delta_vals[9] =
                 { 1e-10, 1e-10, 1e-10, 1e-8, 1e-8, 1e-8, .000001, .000001, .0001 };
//               { 1e-8, 1e-8, 1e-8, 1e-10, 1e-10, 1e-10, .000001, .000001, .0001 };
//               { 2e-8, 2e-8, 2e-8, 1e-10, 1e-10, 1e-10, .0001, .001, .001 };
      const double delta_val = delta_vals[i];
      double worst_error_squared = 0, worst_error_in_arcseconds;

      memcpy( orbit2, orbit, 6 * sizeof( double));
      if( i < 6)           /* adjust position/velocity */
         for( j = 0; j < 6; j++)
            orbit2[j] -= unit_vectors[i][j] * delta_val;
      else                 /* adjust non-gravitational param */
         solar_pressure[i - 6] -= delta_val;

      sprintf( tstr, "Evaluating %d of %d    ", i + 1, n_params);
      set_locs( orbit2, working_epoch, obs, n_obs);
      for( j = 0; j < n_obs; j++)
         {
         slopes[i +   2 * j     * n_params] = obs[j].computed_ra;
         slopes[i + (2 * j + 1) * n_params] = obs[j].computed_dec;
         }

                   /* For computing constrained orbits and element sigmas, */
                   /* we need the orbit at the "real" epoch.  If we don't */
                   /* need sigmas or constraints,  we can save ourselves  */
                   /* some work by skipping the following.                */
      if( sigmas_requested != NO_ORBIT_SIGMAS_REQUESTED || limited_orbit)
         {
         double rel_orbit[6];

         memcpy( rel_orbit, orbit2, 6 * sizeof( double));
         integrate_orbit( rel_orbit, working_epoch, epoch);
         for( j = 0; j < 6; j++)
            rel_orbit[j] -= orbit_offset[j];
                  /* evaluate elements of 'rel_orbit',  then  put */
                  /* into an array form: */
         calc_classical_elements( &elem, rel_orbit, epoch, 1,
                              SOLAR_GM * planet_mass[planet_orbiting]);
         put_orbital_elements_in_array_form( &elem, element_slopes[i]);
         for( j = 0; j < MONTE_N_ENTRIES; j++)
            {
            element_slopes[i][j] -= elements_in_array[j];
            element_slopes[i][j] /= delta_val;
            }
         if( limited_orbit)
            {
            double constraint2[MAX_CONSTRAINTS];

            evaluate_limited_orbit( rel_orbit, planet_orbiting, working_epoch,
                                      limited_orbit, constraint2);
            for( j = 0; j < n_constraints; j++)
               constraint_slope[i + j * n_params] =
                       (constraint2[j] - constraint[j]) / delta_val;
            if( *limited_orbit == 'R')
               {
               const double tconstraint =
                     r_mult * (dotted_dist( obs + n_obs - 1) - atof( limited_orbit + 2));

               constraint_slope[i] = (constraint[0] - tconstraint) / delta_val;
               }
            }
         }
      if( i >= 6)
         solar_pressure[i - 6] += 2. * delta_val;
      if( use_symmetric_derivatives)
         {
         for( j = 0; j < 6; j++)
            orbit2[j] = 2. * orbit[j] - orbit2[j];
         sprintf( tstr, "Evaluating %d of %d rev   ", i + 1, n_params);
         set_locs( orbit2, working_epoch, obs, n_obs);
         }
      else
         for( j = 0; j < n_obs; j++)
            {
            const double orig_computed_ra  = obs[j].ra  + xresids[j];
            const double orig_computed_dec = obs[j].dec + yresids[j];

            obs[j].computed_ra  = 2. * orig_computed_ra  - obs[j].computed_ra;
            obs[j].computed_dec = 2. * orig_computed_dec - obs[j].computed_dec;
            }

      for( j = 0; j < n_obs; j++)
         {
         double *slope_ptr = slopes + i + 2 * j * n_params;
         double error_squared;

         slope_ptr[0] -= obs[j].computed_ra;
         if( slope_ptr[0] > PI)
            slope_ptr[0] -= PI + PI;
         if( slope_ptr[0] < -PI)
            slope_ptr[0] += PI + PI;
         slope_ptr[0] *= cos( obs[j].dec);
         slope_ptr[n_params] -= obs[j].computed_dec;
         error_squared = slope_ptr[0] * slope_ptr[0] +
                     slope_ptr[n_params] * slope_ptr[n_params];
         if( worst_error_squared < error_squared)
            worst_error_squared = error_squared;
         slope_ptr[0]        /= 2. * delta_val;
         slope_ptr[n_params] /= 2. * delta_val;
         }
      worst_error_in_arcseconds =
                 sqrt( worst_error_squared) * 3600. * (180. / PI);
      if( showing_deltas_in_debug_file)
         debug_printf( "Total change on param %d: %lf arcseconds; delta %.3e\n",
            i, worst_error_in_arcseconds, delta_val);
                     /* Attempt to keep the error at 1.5 arcsecs: */
      delta_vals[i] *= 1.5 / worst_error_in_arcseconds;
      if( i >= 6)         /* put solar pressure back where it was: */
         solar_pressure[i - 6] -= delta_val;
      }

// runtime_message = NULL;
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double loc_vals[22];
         const double xresid = xresids[i] * cos( obs[i].dec);
         const double yresid = yresids[i];
         const double one_arcsec = (PI / 180.) / 3600.;

         FMEMCPY( loc_vals, slopes + i * 2 * n_params,
                                         2 * n_params * sizeof( double));
         lsquare_add_observation( lsquare, xresid,
                                          obs[i].weight, loc_vals);
         lsquare_add_observation( lsquare, yresid,
                                          obs[i].weight, loc_vals + n_params);
         sigma_squared +=  obs[i].weight * obs[i].weight
                  * (xresid * xresid + yresid * yresid + one_arcsec * one_arcsec);
         }
   i = n_included_observations * 2 - n_params;
   if( i > 0)
      sigma_squared /= (double)i;

   if( limited_orbit)
      for( j = 0; j < n_constraints; j++)
         lsquare_add_observation( lsquare, constraint[j], 1.,
                                            constraint_slope + j * n_params);

   if( *covariance_filename && sigmas_requested != NO_ORBIT_SIGMAS_REQUESTED)
      {
      FILE *ofile = fopen( covariance_filename, "wb");
      double *matrix = lsquare_covariance_matrix( lsquare);
      double *wtw = lsquare_wtw_matrix( lsquare);
      double eigenvals[10], eigenvects[100], element_sigmas[MONTE_N_ENTRIES];
      int pass;
      char tbuff[20];

      setvbuf( ofile, NULL, _IONBF, 0);
      fprintf( ofile, "Orbit: %.7lf %.7lf %.7lf %.7lf %.7lf %.7lf\nepoch JD %.5lf (%.5lf)\n",
               orbit[0], orbit[1], orbit[2],
               orbit[3], orbit[4], orbit[5], working_epoch, epoch);
      fprintf( ofile, "Unit vectors:\n");
      for( i = 0; i < 6; i++)
         for( j = 0; j < 6; j++)
            {
            put_double_in_buff( tbuff, unit_vectors[i][j]);
            fprintf( ofile, "%s%s", tbuff, (j == 5 ? "\n" : " "));
            }
//    jacobi_eigenvalues( matrix, n_params, eigenvals, eigenvects);
      jacobi_eigenvalues( wtw, n_params, eigenvals, eigenvects);

      for( i = 0; i < n_params; i++)
         eigenvals[i] /= sigma_squared;
      for( i = 0; i < n_params * n_params; i++)
         matrix[i] *= sigma_squared;      /* Danby, p 243, (7.5.21) */
      for( i = 0; i < n_obs; i++)
         if( obs[i].is_included)
            {
            fprintf( ofile, "%4d: ", i);
            for( j = 0; j < n_params * 2; j++)
               {
               put_double_in_buff( tbuff, slopes[i * 2 * n_params + j]);
               fprintf( ofile, "%s", tbuff);
               }
            fprintf( ofile, "\n");
            }
      for( pass = (matrix ? 0 : 2); pass < 4; pass++)
         {
         const char *titles[] = {  "Covariance:\n", "Correlation:\n",
                                               "WtW:\n", "Eigenvectors:\n" };

         fprintf( ofile, "%s", titles[pass]);
         for( i = 0; i < n_params; i++)
            {
            for( j = 0; j < n_params; j++)
               {
               double oval;

               if( pass == 0 || pass == 1)
                  {
                  oval = matrix[i + j * n_params];
                  if( pass == 1)      /* normalize to get correlation, not covar */
                     oval /= sqrt( matrix[i * (n_params + 1)] * matrix[j * (n_params + 1)]);
                  }
               else if( pass == 2)
                  oval = wtw[i + j * n_params];
               else        /* if( pass == 3) */
                  oval = eigenvects[j + i * n_params];
               if( pass == 1 || pass == 3)      /* correlation or eigenvects */
                  sprintf( tbuff, "%10.6lf", oval);  /* values are -1 to 1 */
               else                             /* covar/WtW values can be */
                  put_double_in_buff( tbuff, oval);   /* huge or tiny */
               fprintf( ofile, "%s", tbuff);
               }
            fprintf( ofile, "\n");
            }
         }
      fprintf( ofile, "Eigenvalues:\n");
      for( i = 0; i < n_params; i++)
         {
         put_double_in_buff( tbuff, eigenvals[i]);
         fprintf( ofile, "%s", tbuff);
         }
      for( i = 0; i < n_params; i++)
         one_sigma_eigenvect[i] = 0;
      if( n_params >= 6)
         {
         int k;

         fprintf( ofile, "\n\nNew unit vectors: \n");
         new_unit_vectors = (double **)calloc_double_dimension_array(
                                n_params, n_params, sizeof( double));
         for( k = 0; k < n_params; k++)
            for( i = 0; i < n_params; i++)
               {
               for( j = 0; j < n_params; j++)
                  new_unit_vectors[k][i] +=
                        unit_vectors[j][i] * eigenvects[j + k * n_params];
               put_double_in_buff( tbuff, new_unit_vectors[k][i]);
               fprintf( ofile, "%s%s", tbuff, (i == n_params - 1 ? "\n" : ""));
               }

         for( i = 0; i < n_params; i++)
            one_sigma_eigenvect[i] = new_unit_vectors[0][i] / sqrt( eigenvals[0]);
         fprintf( ofile, "\nOne-sigma eigenvect:\n");
         for( i = 0; i < n_params; i++)
            {
            put_double_in_buff( tbuff, one_sigma_eigenvect[i]);
            fprintf( ofile, "%s", tbuff);
            }
         }
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         double sigma_squared = 0.;
         static const char *text[MONTE_N_ENTRIES] = {
                           "Tp", "e", "q", "Q", "1/a", "i", "M",
                           "omega", "Omega" };

         for( j = 0; j < n_params; j++)
            {
            int k;
            double dot = 0.;

            for( k = 0; k < n_params; k++)
               dot += element_slopes[k][i] * matrix[k + j * n_params];
            sigma_squared += dot * element_slopes[j][i];
            }
         element_sigmas[i] = sqrt( sigma_squared);
         put_double_in_buff( tbuff, sqrt( sigma_squared));
         fprintf( ofile, "\n   %-6s %s", text[i], tbuff);
         fprintf( ofile, "\n  ");
         for( j = 0; j < n_params; j++)
            {
            put_double_in_buff( tbuff, element_slopes[j][i]);
            fprintf( ofile, "%s", tbuff);
            }
         }
      if( n_params == 7)      /* SRP included */
         fprintf( ofile, "\nsigma_AMR: %.6lf\n",
                   sqrt( matrix[n_params * 6 + 6]) * SOLAR_GM / SRP1AU);
      if( n_params == 8 || n_params == 9)
         {           /* comet A1, A2, maybe A3 included */
         char tbuff2[40];

         put_double_in_buff( tbuff, sqrt( matrix[n_params * 6 + 6]));
         put_double_in_buff( tbuff2, sqrt( matrix[n_params * 7 + 7]));
         fprintf( ofile, "\nsigma_A1, A2: %s %s\n", tbuff, tbuff2);
         if( n_params == 9)
            {
            put_double_in_buff( tbuff, sqrt( matrix[n_params * 8 + 8]));
            fprintf( ofile, "\nsigma_A3: %s\n", tbuff);
            }
         }

      fprintf( ofile, "\n\n");
      uncertainty_parameter = dump_monte_data_to_file( ofile, element_sigmas,
            elem.major_axis, elem.ecc, planet_orbiting);
      fclose( ofile);
      free( matrix);
      }
   FFREE( xresids);
   lsquare_solve( lsquare, differences);
   lsquare_free( lsquare);
   for( i = 0; i < 6 && i < n_params; i++)
      for( j = 0; j < 6; j++)
         {
         double max_difference = obs->r * .7, ratio;

         if( j > 2)     /* velocity component */
            max_difference /= (obs[n_obs - 1].jd - obs[0].jd) * .5;
         ratio = fabs( differences[i] * unit_vectors[i][j] / max_difference);
         if( ratio > scale_factor)
            scale_factor = ratio;
         }
   for( i = 0; i < n_params && !err_code; i++)
      {
      if( i < 6)      /* adjust position/velocity */
         for( j = 0; j < 6; j++)
            orbit[j] += unit_vectors[i][j] * differences[i] / scale_factor;
      if( i == 5)    /* is our new 'orbit' state vector reasonable?  */
         err_code = is_unreasonable_orbit( orbit);

      if( i >= 6)
         solar_pressure[i - 6] += differences[i];
      }
               /* If the orbit "blew up" or otherwise failed,  restore */
               /* the original version:  */
   if( err_code || is_unreasonable_orbit( orbit))
      {
      memcpy( orbit, original_orbit, 6 * sizeof( double));
      memcpy( solar_pressure, original_params, 3 * sizeof( double));
      debug_printf( "Failed full step: %d\n", err_code);
      }
   sprintf( tstr, "Final setting of orbit    ");
// debug_printf( "Before rms %.4lf\n", before_rms);
   i = 6;
   do
      {
      if( setting_outside_of_arc)
         set_locs( orbit, working_epoch, obs - n_skipped_obs, n_total_obs);
      else
         set_locs( orbit, working_epoch, obs, n_obs);
#ifdef REMOVE_HALF_STEPS
      after_rms = compute_rms( obs, n_obs);
      debug_printf( "'after' RMS %.4lf\n", after_rms);
      if( after_rms > before_rms * 1.1 && !limited_orbit)
         {
         debug_printf( "Trying half step\n");
//       printf( "Half stepping: %d: %.4lf to %.4lf\n", i, before_rms, after_rms);
         for( j = 0; j < 6; j++)
            orbit[j] = (orbit[j] + orbit_at_working_epoch[j]) * .5;
         i--;
         }
      else
#endif
         i = 0;
      }
      while( i);

   if( new_unit_vectors)
      {
#if 1
      for( i = 0; i < 6; i++)
         for( j = 0; j < 6; j++)
            unit_vectors[i][j] = new_unit_vectors[i][j];
#endif
      free( new_unit_vectors);
      }

   strcpy( tstr, "Return to correct epoch    ");
   if( !err_code)
      integrate_orbit( orbit, working_epoch, epoch);
   runtime_message = NULL;
   return( err_code);
}

static inline double dot_product( const double *a, const double *b)
{
   return( a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* 'score_orbit_arc' basically looks at a series of observations and
assigns a 'score',  in which having a long arc and no really long gaps
in it results in a higher (better) score.  It used to just look at how
long the arc was.  But having decent "fill" during the arc is important,
too.  If you're given two observations mere seconds apart,  plus one a week
later,  that's not as good as having three spaced a day apart each.  There
are other considerations,  too,  such as how much angle the observations
span,  and they really ought to be included.  But at least right now,
they aren't.  */

static inline double score_orbit_arc( const OBSERVE FAR *obs, const int n_obs)
{
   double longest_omission = 0.;
   double arc_length = obs[n_obs - 1].jd - obs[0].jd;
   int i;

   for( i = 0; i < n_obs - 1; i++)
      {
      const double gap = obs[i + 1].jd - obs[i].jd;

      if( longest_omission < gap)
         longest_omission = gap;
      }
   return( arc_length - longest_omission * .9999);
}

/* When looking for subarcs,  we want to look just for those that are less
than 45 degrees long.  Longer than that,  and the usual initial orbit
methods (Gauss,  Herget) may diverge.  (Not necessarily -- in fact,
they'll often still work with much longer arcs -- but 45 degrees means
you can be pretty confident that a decent orbit will be found.)  A quick
dot-product of the unit vectors,  compared to cos(45 degrees) =
sqrt(2) / 2.,  suffices to check the arc length. */

static inline void look_for_best_subarc( const OBSERVE FAR *obs,
       const int n_obs, const double max_arc_len, int *start, int *end)
{
   double best_score = 0., score;
   const double cos_45_deg = 1.414213 / 2.;
   int i, j;

   *start = *end = 0;
   for( i = j = 0; i < n_obs - 1; i++)
      {
      while( j < n_obs - 1 && obs[j + 1].jd - obs[i].jd <= max_arc_len
                  && dot_product( obs[i].vect, obs[j + 1].vect) > cos_45_deg)
         j++;
      score = score_orbit_arc( obs + i, j - i + 1);
      if( score > best_score)
         {
         best_score = score;
         *start = i;
         *end = j;
         }
      }
}

/* The following was created by plotting semimajor axes versus
inclinations as a scatterplot,  from the recent 'mpcorb.dat' file.
The result is a 30x60 graph,  with 30 vertical bins representing
2 degrees each in inclination (from incl=0 at the top to incl=60
at the bottom),  and 60 bins across representing .1 AU in semimajor
axis (from 0 AU to 6 AU).  The natural logarithm of the counts in
each bin are shown;  those ranged from zero to 9,  so I didn't
need to do any scaling.  This basically reflects the scatterplot
provided by MPC at http://www.cfa.harvard.edu/iau/plot/OrbEls51.gif
(except that the vertical axis is flipped).   */

/*           0 AU      1         2         3         4         5         6 */
static const char *incl_vs_a_scattergram[30] = {
/* incl=0*/ "000000012222232223333778877777776331002520000000001332000000",
            "000000112333233333334789988888886331202630100000002442000000",
            "000000112223333233334799888876776332111520000000002442000000",
            "000000022233333333333699888876776342111520010000002442000000",
            "000000011222233233333577787877886443011510000000002442000000",
/* 10 deg*/ "000000011222223323332357678868886443101410000010001442000000",
            "000000012223332232222246689866776331000410000000002441000000",
            "000000012222323222332124688766776331100310000000001340000000",
            "000000002222222223462222466655786331011200000000001341000000",
            "000000012122223223562233355544675332000100000000000431000000",
/* 20 deg*/ "000000011122232223562255456433565210010100000000001331000000",
            "000000111222222123662366555423565210100100000000000331000000",
            "000000011122223222562356444533464100000000000000001330000000",
            "000000012112212222442244355422464100000000000000000321000000",
            "000000011112122122231222156421353110000000000000000321000000",
/* 30 deg*/ "000000011112111122221111245410132000000000000000000230000000",
            "000000011110011212211111135410122000000000000000000220000000",
            "000000001111211221221120033320011100000000000000000120000000",
            "000000000112111211201110122310000000000000000000000110000000",
            "000000000011122212201110111110000000000000000000000000000000",
/* 40 deg*/ "000000001101111121111110101100001000000000000000000000000000",
            "000000000011010001101110000001000000000000000000000000000000",
            "000000000001100001101010010100000000000000000000000000000000",
            "000000000000011011100110010000000000000000000000000000000000",
            "000000000100100000000000000000000000000000000000000000000000",
/* 50 deg*/ "000000001000100100100000000000000000000000000000000000000000",
            "000000000000100011010100000000000000000000000000000000000000",
            "000000000000011000110000000000000000000000000000000000000000",
            "000000000000010000000011000000000000000000000000000000000000",
            "000000000000000000000000000000000000000000000000000000000000" };

/* evaluate_initial_orbit( ) is supposed to give a "score" of sorts,      */
/* describing just how likely this orbit seems to be,  given its rms      */
/* errors (lower is better);  eccentricity (same,  with highly hyperbolic */
/* solutions considered very unlikely and therefore bumping up the return */
/* value);  inclination (anything above .5 radians,  or about 30 degrees, */
/* is considered unlikely).  Also,  main-belt objects (based on a) are    */
/* slightly encouraged.                                                   */

double evaluate_initial_orbit( const OBSERVE FAR *obs,
                              const int n_obs, const double *orbit)
{
   extern double planet_mass[];
   double rms_err = compute_rms( obs, n_obs), rval, rel_orbit[6];
   ELEMENTS elem;
   int planet_orbiting = find_best_fit_planet( obs->jd,
                                  orbit, rel_orbit);

   calc_classical_elements( &elem, rel_orbit, obs[0].jd, 1,
                              SOLAR_GM * planet_mass[planet_orbiting]);
   rval = rms_err + elem.ecc / 2.;
   if( !planet_orbiting)                  /* for heliocentric orbits... */
      {
      const int xbin = (int)( elem.major_axis * 10.);
      const int ybin = (int)( elem.incl * (180. / PI) / 2.);

//    if( elem.ecc > .5)                     /* discourage eccentric orbits */
//       rval += (elem.ecc - .5) * 2.;
      if( elem.ecc > 1.01)                /* _strongly_ discourage hyperbolics */
         rval += (elem.ecc - 1.01) * 1000.;
      if( elem.incl > .5)                 /* gently discourage high-incl and */
         rval += (elem.incl - .5) * .2;   /* retrograde orbits */
      if( xbin >= 0 && ybin >= 0 && xbin < 60 && ybin < 30)
         rval -= (double)( incl_vs_a_scattergram[ybin][xbin] - '0') * .1;
      }
   return( rval);
}

static double attempt_improvements( double *orbit, OBSERVE *obs, const int n_obs)
{
   int method;
   double curr_score;

#ifdef CONSOLE
   if( show_runtime_messages)
      move_add_nstr( 14, 10, "Improving solution...        ", -1);
#endif
   set_locs( orbit, obs[0].jd, obs, n_obs);
   curr_score = evaluate_initial_orbit( obs, n_obs, orbit);
   for( method = 0; method < 2; method++)
      {
      int significant_improvement_occurred = 1;
      int iter = 0;
      int max_iter = 5;
      double temp_orbit[6];

                        /* We're willing to try the Herget,  then full  */
                        /* step methods, five times... _if_ they result */
                        /* in real improvement,  and we mandate two     */
                        /* tries with both Herget and full step  */
      while( significant_improvement_occurred && iter++ < max_iter)
         {
         double score;

         memcpy( temp_orbit, orbit, 6 * sizeof( double));
#ifdef CONSOLE
         if( show_runtime_messages)
            {
            char msg_buff[80];

            sprintf( msg_buff, "%s step: radii %lf, %lf",
                        (method ? "full" : "Herget"),
                        obs[0].r, obs[n_obs - 1].r);
            move_add_nstr( 14, 10, msg_buff, -1);
            }
#endif
         if( !method)         /* doing an Herget step */
            {
            double r1 = obs[0].r, r2 = obs[n_obs - 1].r;
            double d_r1, d_r2;

            herget_method( obs, n_obs, r1, r2, temp_orbit, &d_r1, &d_r2, NULL);
            r1 += d_r1;
            r2 += d_r2;
            adjust_herget_results( obs, n_obs, temp_orbit);
            if( debug_level > 3)
               debug_printf( "Adjusting Herget results\n");
            }
         else        /* doing a full step */
            full_improvement( obs, n_obs, temp_orbit, obs[0].jd, NULL,
                           NO_ORBIT_SIGMAS_REQUESTED);

         score = evaluate_initial_orbit( obs, n_obs, temp_orbit);
         if( debug_level > 2)
            debug_printf( "Method %d, run %d: score %lf\n",
                                    method, iter, score);
                     /* If we didn't improve things by a score of .1, */
                     /* we might as well stop:                        */
         if( score > curr_score - .1 && iter > 2)
            significant_improvement_occurred = 0;
                     /* If things totally fell apart, also stop:   */
         if( score > 50000.)
            significant_improvement_occurred = 0;
         if( score < curr_score)
            {
            memcpy( orbit, temp_orbit, 6 * sizeof( double));
            curr_score = score;
            if( debug_level > 2)
               debug_printf( "Improvement accepted\n");
            }
         }
      }
   return( curr_score);
}

#define INITIAL_ORBIT_NOT_YET_FOUND       -2
#define INITIAL_ORBIT_FAILED              -1
#define INITIAL_ORBIT_FOUND                0

double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit)
{
   int i, rval = INITIAL_ORBIT_NOT_YET_FOUND;
   int start = 0;
   double arclen;
#ifdef CONSOLE
   char msg_buff[80];
#endif

   if( debug_level)
      debug_printf( "initial_orbit(): %d obs;", n_obs);
   assert( orbit);
   for( i = 0; i < n_obs; i++)
      {
      obs[i].computed_ra  = obs[i].ra;
      obs[i].computed_dec = obs[i].dec;
      }

   while( n_obs && !obs->is_included)
      {
      obs++;
      n_obs--;
      }
   while( n_obs && !obs[n_obs - 1].is_included)
      n_obs--;
   if( debug_level)
      debug_printf( "  %d left\n", n_obs);
   if( n_obs < 2)
      return( -1);
   arclen = obs[n_obs - 1].jd - obs[0].jd;
   if( arclen > 730.)         /* two-year maximum */
      arclen = 730.;

   perturbers = AUTOMATIC_PERTURBERS;
   while( rval == INITIAL_ORBIT_NOT_YET_FOUND)
      {
      int end, n_subarc_obs;
      const double max_arg_length_for_vaisala = 30.;
      double best_score = 1e+300;
      double best_orbit[6];

      for( i = 0; i < 6; i++)
         best_orbit[i] = 0.;
      look_for_best_subarc( obs, n_obs, arclen, &start, &end);
      if( debug_level > 1)
         debug_printf( "  Current arc: %d to %d\n", start, end);
      if( start == end)       /* no such arc found */
         return( INITIAL_ORBIT_FAILED);
      for( i = 0; i < start; i++)
         obs[i].is_included = 0;
      for( i = start; i <= end; i++)
         obs[i].is_included = 1;
      for( i = end + 1; i < n_obs; i++)
         obs[i].is_included = 0;
      for( i = 0; i < n_obs; i++)      /* solely to ensure a non-zero r */
         obs[i].r = 1.;
      arclen = obs[end].jd - obs[start].jd;
      if( debug_level)
         debug_printf( "From %lf to %lf (%lf days)\n", obs[start].jd, obs[end].jd, arclen);
      n_subarc_obs = end - start + 1;

      if( n_subarc_obs >= 3)     /* at least three observations;  try Gauss */
         {
#ifdef CONSOLE
         if( show_runtime_messages)
            move_add_nstr( 14, 10, "In Gauss solution", -1);
#endif
         for( i = 0; i < 3 && rval == INITIAL_ORBIT_NOT_YET_FOUND; i++)
            {
            double epoch =
                convenient_gauss( obs + start, n_subarc_obs, orbit, 1., i);

            if( debug_level)
               debug_printf( "Gauss epoch: JD %lf\n", epoch);
            if( epoch)
               {
               double score;

               if( debug_level > 2)
                  debug_printf( "About to set locs\n");
               set_locs( orbit, epoch, obs + start, n_subarc_obs);
               score = evaluate_initial_orbit( obs + start, n_subarc_obs, orbit);
               if( debug_level > 2)
                  debug_printf( "Locations set; score %lf\n", score);
               if( score < 1000.)
                  {
                  integrate_orbit( orbit, epoch, obs[start].jd);
                  epoch = obs[start].jd;
                  score = attempt_improvements( orbit, obs + start, n_subarc_obs);
                  if( debug_level > 2)
                     debug_printf( "After improvements: score %lf\n", score);
                  }
               if( score < best_score)
                  {
                  best_score = score;
                  integrate_orbit( orbit, epoch, obs[start].jd);
                  memcpy( best_orbit, orbit, 6 * sizeof( double));
                  if( best_score < 2.)   /* good enough result: declare */
                     rval = INITIAL_ORBIT_FOUND;     /* victory/go home */
                  }
               }
            else        /* break out of Gauss loop */
               i = 3;
            }
         }           /* end of trying Gauss */
#ifdef CONSOLE
      if( show_runtime_messages)
         move_add_nstr( 14, 10, "Gauss done", -1);
#endif
      if( arclen < max_arg_length_for_vaisala && best_score > 2.)
         for( i = 0; i < 2; i++)
            {
            int orbit_looks_reasonable = 1;
            double pseudo_r;

            if( i)          /* dist from observer (second) pass:  some ad hoc */
               {            /* code that says,  "for long arcs,  start farther */
                            /* from the observer".                             */
               pseudo_r = arclen / 500.;
               if( pseudo_r < .00001)     /* 1e-5 AU = 1490 km */
                  pseudo_r = .00001;
               }
            else                  /* (first) Vaisala pass */
               pseudo_r = .1;

            while( pseudo_r < (i ? 5. : 100.) && orbit_looks_reasonable)
               {
               double pseudo_r_to_use;
               int herget_rval;
               double score;

               if( i)                          /* 2nd pass, dist from observer */
                  pseudo_r_to_use = pseudo_r;
               else                            /* 1st pass, dist from sun */
                  pseudo_r_to_use = -(1. + pseudo_r);
               herget_rval = herget_method( obs + start, n_subarc_obs,
                                    pseudo_r_to_use,
                                    pseudo_r_to_use,
                                    orbit, NULL, NULL, NULL);
               if( herget_rval < 0)    /* herget method failed */
                  score = 1.e+7;
               else if( herget_rval > 0)        /* vaisala method failed, */
                  score = 9e+5;                 /* but we should keep trying */
               else
                  {
                  adjust_herget_results( obs + start, n_subarc_obs, orbit);
                  score = evaluate_initial_orbit( obs + start, n_subarc_obs, orbit);
                  }
               if( debug_level > 2)
                  debug_printf( "%d, pseudo-r %lf: score %lf, herget rval %d\n",
                         i, pseudo_r, score, herget_rval);
               if( score < best_score)
                  {
                  best_score = score;
                  memcpy( best_orbit, orbit, 6 * sizeof( double));
                  }
#ifdef CONSOLE
               if( show_runtime_messages)
                  {
                  sprintf( msg_buff, "Method %d, r=%.4lf", i, pseudo_r);
                  move_add_nstr( 14, 10, msg_buff, -1);
                  }
#endif
               if( score > 5e+4)   /* usually means eccentricity > 100! */
                  {
                  orbit_looks_reasonable = 0;      /* should stop looking */
                  if( debug_level > 2)
                     debug_printf( "%d: Flipped out at %lf\n", i, pseudo_r);
                  }
               pseudo_r *= 1.2;
               }
            }
      if( best_score < 5. && rval == INITIAL_ORBIT_NOT_YET_FOUND)
         {           /* got a good orbit using Vaisala or Herget */
         rval = INITIAL_ORBIT_FOUND;
         memcpy( orbit, best_orbit, 6 * sizeof( double));
         if( n_obs > 2)
            attempt_improvements( orbit, obs + start, n_subarc_obs);
         }
      arclen *= .7;     /* If we failed,  try again with a shorter arc */
      }

   perturbers_automatically_found = 0;
   if( rval != INITIAL_ORBIT_FOUND)
      {
      start = 0;
      herget_method( obs, n_obs, 1., 1., orbit, NULL, NULL, NULL);
      }
   set_locs( orbit, obs[start].jd, obs, n_obs);
// perturbers = old_perturbers;        /* restore perturbers, if any */
   perturbers = perturbers_automatically_found & (~AUTOMATIC_PERTURBERS);
   return( obs[start].jd);    /* ...and return epoch = JD of first observation */
}

/* Examine 'peirce.cpp' for details on the following.  This is just
   the Rayleigh-distribution version,  not the root-finding one needed
   for a Gaussian distribution. */

double peirce_rayleigh_func( const int N, const int n, const int m)
{
   const double a = (double)(N - n - m) / (double)n;
   const double b = 2. / (double)(N - n);
   const double part_1 = (double)n / (double)(N - n);
   const double part_2 = (double)( N - n) / (double)N;
   const double log_P = (double)n * log( part_1) + (double)N * log( part_2);
   const double log_R = -.5;
   const double tval = exp( b * (log_P - (double)n * log_R));

   return( sqrt( 1. + a * (1 - tval)));
}
