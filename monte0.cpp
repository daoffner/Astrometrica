#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "afuncs.h"
#include "monte0.h"
#include "mt64.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

void put_orbital_elements_in_array_form( const ELEMENTS *elem,
                  double *output_array)
{
   output_array[MONTE_TP] = elem->perih_time;
   output_array[MONTE_ECC] = elem->ecc;
   output_array[MONTE_q] = elem->major_axis * (1. - elem->ecc);
   output_array[MONTE_Q] = elem->major_axis * (1. + elem->ecc);
   output_array[MONTE_INV_A] = 1. / elem->major_axis;
   output_array[MONTE_INCL] = elem->incl * 180. / PI;
   output_array[MONTE_MEAN_ANOM] = elem->mean_anomaly * 180. / PI;
   output_array[MONTE_ARG_PER] = elem->arg_per * 180. / PI;
   output_array[MONTE_ASC_NODE] = elem->asc_node * 180. / PI;
}

void add_monte_orbit( double *monte_data, const ELEMENTS *elem,
                  const int n_orbits)
{
   double tarr[MONTE_N_ENTRIES];
   double *offsets = monte_data + 2 * MONTE_N_ENTRIES;
   int i;

   put_orbital_elements_in_array_form( elem, tarr);
   if( !n_orbits)                /* initializing step */
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         offsets[i] = tarr[i];
         monte_data[i] = monte_data[i + MONTE_N_ENTRIES] = 0.;
         }
   else
      {
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         double delta = tarr[i] - offsets[i];

         if( i >= MONTE_INCL && i <= MONTE_ASC_NODE)
            {
            if( delta > 180.)            /* keep angular arguments in the */
               delta -= 360.;            /* proper range,  not wrapping   */
            else if( delta < -180.)      /* around at +/- 180 degrees     */
               delta += 360.;
            }
         monte_data[i] += delta;
         monte_data[i + MONTE_N_ENTRIES] += delta * delta;
         }
      }
}

void compute_monte_sigmas( double *sigmas, const double *monte_data,
                  const int n_orbits)
{
   int i;

   for( i = 0; i < MONTE_N_ENTRIES; i++)
      {
      const double avg_square = monte_data[i + MONTE_N_ENTRIES] / (double)n_orbits;
      const double avg_value = monte_data[i] / (double)n_orbits;

      sigmas[i] = sqrt( avg_square - avg_value * avg_value);
      }
}


/* Add some Gaussian noise to each RA/dec value.  Two Gaussian-distributed
random numbers are generated using the Box-Muller transform,  scaled
according to the observation weight;  this is   then added to the RA and
dec. the noise_in_arcseconds.  The original RA/decs are stored in an array;
calling 'remove_gaussian_noise_from_obs()' restores them, removing the
noise.

   2012 Jul 22:  Marco Micheli pointed out that 'r' should be scaled by
the observation weight.  It now is.  And yes,  it should have been that
way right from the beginning...

   2012 Feb 9:  switched to use of MT64 (Mersenne Twister,  64-bit version)
for generating pseudo-random numbers.  Default C-library PRNGs are sometimes
adequate and sometimes horrible.  MT64 is quite well thought of,  and should
enable us to assume that any problems are _not_ due to insufficiently random
numbers.  Also, if one uses the C-library PRNG,  you'll get different results
on different systems with different implementations.  */

double *add_gaussian_noise_to_obs( int n_obs, OBSERVE *obs,
                 const double noise_in_arcseconds)
{
   const double noise_in_radians = noise_in_arcseconds * PI / (180. * 3600.);
   double *stored_ra_decs = (double *)calloc( 2 * n_obs, sizeof( double));
   double *tptr = stored_ra_decs;
   static uint64_t *mt64_state = NULL;

   if( !mt64_state)
      {
      const uint64_t mt64_seed_value = (uint64_t)0x31415926;

      mt64_state = (uint64_t *)calloc( MT_STATE_SIZE, sizeof( uint64_t));
      init_mt64( mt64_seed_value, mt64_state);
      }
   while( n_obs--)
      {
      const double rt = log( 1. - mt64_double( mt64_state));
      const double r = sqrt( -2. * rt) / obs->weight;
      const double theta = 2. * PI * mt64_double( mt64_state);
      const double dx = r * cos( theta) * noise_in_radians;
      const double dy = r * sin( theta) * noise_in_radians;

      *tptr++ = obs->ra;
      *tptr++ = obs->dec;
      obs->ra += dx / cos( obs->dec);
      obs->dec += dy;
      obs++;
      }
   return( stored_ra_decs);
}

void remove_gaussian_noise_from_obs( int n_obs, OBSERVE *obs,
                           double *stored_ra_decs)
{
   double *tptr = stored_ra_decs;

   while( n_obs--)
      {
      obs->ra = *tptr++;
      obs->dec = *tptr++;
      obs++;
      }
   free( stored_ra_decs);
}

   /* For some time,  I displayed the covariance matrix and sigmas using the
      sprintf format specifier %10.3g.  That worked well,  except that values
      such as 1010 or 999999 were rendered as 1.01e+003 or 9.99e+005.  (999
      was left as 999,  and I've no problem with the use of scientific notation
      beyond six digits... purely a matter of personal preference;  I didn't
      want four,  five,  and six-digit numbers shown in SN.)  Also,  this
      code will remove unnecessary zeroes in the exponent,  so that 'e+004',
      for example,  becomes 'e+4'.           */

char *put_double_in_buff( char *buff, const double ival)
{
   if( fabs( ival) < 999.999 || fabs( ival) > 999999.)
      {
      sprintf( buff, "%10.3g", ival);
      if( buff[5] == 'e' && buff[7] == '0')
         {
         memmove( buff + 1, buff, 7);
         *buff = ' ';
         }
      if( buff[6] == 'e' && buff[8] == '0')
         {
         memmove( buff + 1, buff, 8);
         *buff = ' ';
         }
      }
   else
      sprintf( buff, "%10d", (int)ival);
   while( *buff == ' ')
      buff++;
   return( buff);
}

/* Just to explicate some of the error calculus below:

   We actually don't know sigma_a right off the bat,  since a is
not entirely defined for near-parabolic orbits.  Instead,  we have
sigma(1/a),  from which we can compute

sigma_a = sigma(1/a) * a^2

   Then,  since P_years = a^1.5,

sigma_P_years = sigma_a * 1.5 * sqrt(a)

   ...and,  since n = 360. / period in days = 360 / (days_per_year * P_years),

sigma_n = 360 * sigma_P_in_days / P_days^2
*/

double dump_monte_data_to_file( FILE *ofile, const double *sigmas,
            const double semimajor_axis, const double ecc,
            const int planet_orbiting)
{
   double uparam = 97.;   /* assume we won't get a "correct" uncertainty */
   static const char *text[MONTE_N_ENTRIES] = {
                           "Tp", "e", "q", "Q", "1/a", "i", "M",
                           "omega", "Omega" };
   static const char *units_text[MONTE_N_ENTRIES] = {
                           "days", "", "AU", "AU", "1/AU", "deg",
                           "deg", "deg", "deg" };
   const double sigma_a = sigmas[MONTE_INV_A] * semimajor_axis * semimajor_axis;
   int i;
   char tbuff[20], *tptr;

   fprintf( ofile, "Sigmas:\n");
   for( i = 0; i < MONTE_N_ENTRIES; i++)
      {
      double oval = sigmas[i];

      if( !strcmp( units_text[i], "deg"))    /* sigma's really in radians; */
         oval *= 180. / PI;                  /* show it in degrees  */
      put_double_in_buff( tbuff, sigmas[i]);
      fprintf( ofile, "sigma_%-5s%s %s",
                  text[i], tbuff, units_text[i]);
      if( !strcmp( units_text[i], "AU"))    /* show in km,  too */
         {
         tptr = put_double_in_buff( tbuff, sigmas[i] * AU_IN_KM);
         fprintf( ofile, " (%s km)", tptr);
         }
      fprintf( ofile, "\n");
      }
// if( semimajor_axis > sigma_a * .3)
   if( semimajor_axis > 0.)
      {
//    const double sigma_P_in_days = 365.25 * 1.5 * sigma_a *
//        sqrt( semimajor_axis / planet_mass[planet_orbiting]);
//    const double per_yrs = elem.t0 * (2. * PI / 365.25);
      extern double planet_mass[];
      const double per_yrs = semimajor_axis * sqrt( semimajor_axis / planet_mass[planet_orbiting]);
      const double days_per_year = 365.25;
      const double per_days = per_yrs * days_per_year;
//    const double sigma_P_in_days = days_per_year * 1.5 * sigma_a
//            * sqrt( semimajor_axis / planet_mass[planet_orbiting]);
      const double sigma_P_in_days = per_days * 1.5 * sigma_a / semimajor_axis;
      const double GAUSS_K = .01720209895;
      const double sigma_n = 360. * sigma_P_in_days / (per_days * per_days);
      const double runoff_coeff =
               3. * 3600. * (180. / PI) * GAUSS_K;
      const double runoff = (runoff_coeff / per_yrs) *
             (sigmas[0] * ecc + 10. * sigma_P_in_days / per_yrs);
      const double uparam_const = 1.49;  /* =ln(648000/9) */

      uparam = log( runoff) / uparam_const + 1.;

      if( semimajor_axis > sigma_a * .3)
         {
         put_double_in_buff( tbuff, sigma_n);
         fprintf( ofile, "sigma_n:   %s\n", tbuff);
         put_double_in_buff( tbuff, sigma_a);
         fprintf( ofile, "sigma_a:   %s AU", tbuff);
         tptr = put_double_in_buff( tbuff, sigma_a * AU_IN_KM);
         fprintf( ofile, " (%s km)\n", tptr);
         if( sigma_P_in_days < 999.)
            {
            put_double_in_buff( tbuff, sigma_P_in_days);
            fprintf( ofile, "sigma_P:   %s days\n", tbuff);
            }
         else
            fprintf( ofile, "sigma_P: %12.3lg years\n",
                        sigma_P_in_days / 365.25);
         }
      fprintf( ofile, "P = %.2lf years; U=%.1lf\n",
                  per_yrs, uparam);
      }
   return( uparam);
}

/* From http://www.minorplanetcenter.net/iau/info/UValue.html :

The U value is calculated in the following manner. First, calculate:

      RUNOFF = (dT * e + 10 / P * dP) * GAUSS_K * (180. / PI) / P * 3600 * 3

where dT is the uncertainty in the perihelion time (in days)
      e is the eccentricity
      P is the orbital period (in years)
      dP is the uncertainty in the orbital period (in days)
      ko is the Gaussian constant in degrees
         = 180 / pi * 0.01720209895
      3600 converts to seconds of arc
      3 is a empirical factor to make the formal errors more
         closely model reality
and   RUNOFF is the in-orbit longitude runoff in seconds of
         arc per decade
*/
