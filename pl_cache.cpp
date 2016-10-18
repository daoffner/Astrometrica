/* pl_cache.cpp: computes and caches planetary positions

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
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"
#include "jpleph.h"

const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int debug_printf( const char *format, ...);                /* runge.cpp */
extern int debug_level;

#define J2000 2451545.0
#define J0 (J2000 - 2000. * 365.25)

char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int generic_message_box( const char *message, const char *box_type);

int format_jpl_ephemeris_info( char *buff);           /* pl_cache.cpp */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end);
int compute_rough_planet_loc( const double t_cen, const int planet_idx,
                                          double *vect);    /* sm_vsop.cpp */

static int planet_posn_raw( const int planet_no, const double jd,
                            double *vect_2000)
{
   static void *ps_1996_data[10];
   const int jpl_center = 11;         /* default to heliocentric */
   int i, rval = 0;
   static const char *jpl_filename = NULL;
   static void *jpl_eph = NULL;

   if( !planet_no)            /* the sun */
      {
      vect_2000[0] = vect_2000[1] = vect_2000[2] = 0.;
      if( !jd && jpl_eph)       /* return version data: */
         {
         vect_2000[0] = (double)jpl_get_long( jpl_eph, JPL_EPHEM_EPHEMERIS_VERSION);
         vect_2000[1] = jpl_get_double( jpl_eph, JPL_EPHEM_START_JD);
         vect_2000[2] = jpl_get_double( jpl_eph, JPL_EPHEM_END_JD);
         }
      return( 0);
      }

   if( !jpl_filename)
      {
      FILE *ifile;

#if defined (_WIN32) || defined( __WATCOMC__)
      jpl_filename = get_environment_ptr( "JPL_FILENAME");
#else
      jpl_filename = get_environment_ptr( "LINUX_JPL_FILENAME");
#endif
      if( *jpl_filename)
         jpl_eph = jpl_init_ephemeris( jpl_filename, NULL, NULL);
      if( !jpl_eph)
         if( (ifile = fopen( "jpl_eph.txt", "rb")) != NULL)
            {
            char buff[100];

            while( !jpl_eph && fgets_trimmed( buff, sizeof( buff), ifile))
               if( *buff && *buff != ';')
                  jpl_eph = jpl_init_ephemeris( buff, NULL, NULL);
            if( debug_level)
               debug_printf( "Ephemeris file %s\n", buff);
            fclose( ifile);
            }
      if( debug_level && jpl_eph)
         {
         debug_printf( "\nEphemeris time span years %.3lf to %.3lf\n",
               (jpl_get_double( jpl_eph, JPL_EPHEM_START_JD) - J0) / 365.25,
               (jpl_get_double( jpl_eph, JPL_EPHEM_END_JD)   - J0) / 365.25);
         debug_printf( "Ephemeris version %d\n", jpl_get_long( jpl_eph, JPL_EPHEM_EPHEMERIS_VERSION));
         debug_printf( "Kernel size %d, record size %d, swap_bytes %d\n",
               jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_SIZE),
               jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_RECORD_SIZE),
               jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_SWAP_BYTES));
         debug_printf( "ncon = %d AU=%lf emrat = %lf\n",
               jpl_get_long( jpl_eph, JPL_EPHEM_N_CONSTANTS),
               jpl_get_double( jpl_eph, JPL_EPHEM_AU_IN_KM),
               jpl_get_double( jpl_eph, JPL_EPHEM_EARTH_MOON_RATIO));
         }
      }

   if( jpl_eph)
      {
      double state[6];            /* DE gives both posn & velocity */
      int failure_code;

      if( planet_no < 0)          /* flag to unload everything */
         {
         jpl_close_ephemeris( jpl_eph);
         jpl_eph = NULL;
         jpl_filename = NULL;
         return( 0);
         }
      else if( planet_no == 10)
         failure_code = jpl_pleph( jpl_eph, jd, 10, 3, state, 0);
      else
         failure_code = jpl_pleph( jpl_eph, jd,
                 (planet_no == 3) ? 13 : planet_no, jpl_center, state, 0);
      if( !failure_code)         /* we're done */
         {
         if( debug_level > 5)
            debug_printf( "JD %lf, planet %d: (%lf %lf %lf)\n",
                     jd, planet_no, state[0], state[1], state[2]);
         equatorial_to_ecliptic( state);
         memcpy( vect_2000, state, 3 * sizeof( double));
         return( 0);
         }
      else
         if( debug_level)
            debug_printf( "Failed: JD %lf, planet %d, code %d\n",
                           jd, planet_no, failure_code);
      }

   if( planet_no == 10)        /* the moon */
      {
      double tloc[4];

      if( !compute_elp_xyz( NULL, (jd - J2000) / 36525., 0., tloc))
         for( i = 0; i < 3; i++)
            vect_2000[i] = tloc[i] / AU_IN_KM;
      else
         {
         static int first_time = 1;

         rval = -3;
         if( first_time)
            generic_message_box( "\
The file 'elp82.dat' is missing.  Without\n\
this file,  lunar positions will be very inaccurate\n\
and orbits will be somewhat less accurate.  Either\n\
get the file,  or set up JPL ephemerides.  See\n\
http://www.projectpluto.com/find_orb.htm#de_eph\n\
for details on how to do this.", "o");
         first_time = 0;
         compute_rough_planet_loc( (jd - J2000) / 36525., 10, vect_2000);
         }
      return( rval);
      }

   if( planet_no < 0)          /* flag to unload everything */
      {
      for( i = 0; i < 10; i++)
         if( ps_1996_data[i])
            {
            unload_ps1996_series( ps_1996_data[i]);
            ps_1996_data[i] = NULL;
            }
      return( 0);
      }

   if( !ps_1996_data[planet_no])
      ps_1996_data[planet_no] = load_ps1996_series( NULL, jd, planet_no);

   if( !ps_1996_data[planet_no])
      rval = -1;
   else if( get_ps1996_position( jd, ps_1996_data[planet_no], vect_2000, 0))
      {
      unload_ps1996_series( ps_1996_data[planet_no]);
      ps_1996_data[planet_no] = load_ps1996_series( NULL, jd, planet_no);
      if( !ps_1996_data[planet_no])
         rval = -2;
      else if( get_ps1996_position( jd, ps_1996_data[planet_no], vect_2000, 0))
         rval = -3;
      }

   if( !rval)
      equatorial_to_ecliptic( vect_2000);
   else
      {
      static int first_time = 1;

      if( first_time)
         {
         generic_message_box( "\
The file 'ps_1996.dat' is missing.  Without it,\n\
planetary positions will be of low accuracy.  Either\n\
get the file,  or set up JPL ephemerides.  See\n\
http://www.projectpluto.com/find_orb.htm#de_eph\n\
for details on how to do this.", "o");
         debug_printf( "Loading ps_1996: rval %d, planet %d, JD %lf\n",
                  rval, planet_no, jd);
         }
      first_time = 0;
      if( planet_no > 0 && planet_no < 9)
         compute_rough_planet_loc( (jd - J2000) / 36525., planet_no, vect_2000);
      }

   return( rval);
}

#define POSN_CACHE struct posn_cache

POSN_CACHE
   {
   double jd;
   double vect[3];
   int8_t planet_no;
   };

int n_posns_cached = 0;

/* Hash the JD and planet number by putting them in a buffer,  then taking
each byte from said buffer,  XORing it with the current return value,  then
multiplying the current return value by a big prime.  The end result is
moduloed by the table size.   Seems to result in a suitably "pseudo-random"
hash table filling algorithm.  (This was borrowed from a hash algorithm
used in the Linux ext3 file system.  It's actually pretty slow,  but for
this use,  the time it adds is negligible compared to everything else
that's going on.  You may run into trouble if you use it elsewhere.)

   I've tried some alternative hash algorithms,  in particular because
of a weird compiler bug when compiling in Puppy Linux (and also reported
with some versions of g++ on Ubuntu 8.04).  For some reason,  setting rval
to be positive with the commented-out code failed;  rval remained negative,
resulting in a segmentation fault.  To evade this,  rval is now rendered
positive with 'rval &= 0x7fffffff' instead.   */

static long hash_function( const int planet_no, const double jd,
                     const long table_size)
{
#ifdef ALT_HASH_ALGORITHM
   long rval = planet_no * 314159 + (long)( (jd - 2400000.) * (jd - J2000));
#else
   const int bufflen = sizeof( double) + sizeof( int);
   unsigned char buff[sizeof( double) + sizeof( int)];
   long rval = 1;
   const long big_prime = 1234567891;
   int i;

   memcpy( buff, &jd, sizeof( double));
   memcpy( buff + sizeof( double), &planet_no, sizeof( int));
   for( i = 0; i < bufflen; i++)
      rval = (rval ^ buff[i]) * big_prime;
#endif
   rval &= 0x7fffffff;
   rval %= table_size;
// if( rval < 0)
//    rval += table_size;
   return( rval);
}

/* Computing planetary positions is somewhat expensive if we're using
JPL ephemerides,  and _very_ expensive if we aren't (the PS-1996 method
is used,  which involves lots of trig series).  And frequently,  we'll be
requesting the same data over and over (for example,  if we're integrating
over a particular time span repeatedly).  So it makes sense to cache the
planetary positions.

   Below,  this is done with a hash table in a very standard sort of way.
The planet_no and jd are hashed, we look in the 'cache' table,  we do a
quadratic search if there's an hash collision.  If we find the data,  we
return it.  If we don't find it,  we call the planet_posn_raw( ) (uncached)
function, and add the result to the cache, and _then_ return it.

   When the table is more than 80% full,  we double the table size,  dump
everything computed to date,  and start from scratch.  This is admittedly
mildly wasteful,  but I don't think the performance benefit of expanding
the cache and adding everything we've got back in would be worthwhile.
(It wouldn't be hard to do,  though.)

   If the cache gets above some limit (currently set to a million cached
positions,  or about 36 MBytes),  we stop growing the cache.  So if you
had a _really_ long integration,  the cache gets dumped,  rebuilt to the
same size,  dumped again,  built to the same size,  etc.

   Previously,  the data was stored using a balanced tree.  I don't know
what possessed me to do something that dumb.  (At the very least,  had
I keyed the tree using the above hash function,  entries to the tree would
have been nearly random,  and a plain old unbalanced tree would have worked
Just Fine.)

   The following three long ints keep track of the number of searches and
probes done,  and the "worst-case" maximum number of probes required,
just to check that the hash function is truly random enough.      */

#define TEST_PLANET_CACHING_HASH_FUNCTION

#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
long total_n_searches = 0, total_n_probes = 0, max_probes_required = 0;
#endif

int planet_posn( const int planet_no, const double jd, double *vect_2000)
{
   static POSN_CACHE *cache = NULL;
   static POSN_CACHE **hash_table = NULL;
   int loc, found_it = 0, rval = 0, n_probes = 0;
   static int table_size = 100000;

   if( !planet_no)            /* the sun */
      {
      vect_2000[0] = vect_2000[1] = vect_2000[2] = 0.;
      return( 0);
      }

   if( planet_no < 0 || n_posns_cached >= table_size * 4 / 5)
      {                                  /* flag to unload everything */
      if( planet_no < 0)
         planet_posn_raw( -1, 0., NULL);
      if( hash_table)
         free( hash_table);
      if( cache)
         free( cache);
      cache = NULL;
      hash_table = NULL;
      n_posns_cached = 0;
      if( planet_no < 0)
         return( 0);
      }

   if( !cache)
      {
      static int max_table_size;

      if( !max_table_size)
         {
         max_table_size = atoi( get_environment_ptr( "MAX_CACHE"));
         if( !max_table_size)
            max_table_size = 3000000;  /* Default cache size.  This is only */
         }                             /* about 108 MBytes,  so we could    */
      table_size *= 2;                 /* consider raising this limit      */
      if( table_size > max_table_size)
         table_size = max_table_size;
      cache = (POSN_CACHE *)calloc( table_size * 4 / 5, sizeof( POSN_CACHE));
      hash_table = (POSN_CACHE **)calloc( table_size, sizeof( POSN_CACHE *));
      }

   loc = hash_function( planet_no, jd, table_size);

               /* Try to find a table entry with the JD and planet number */
               /* that have been requested: */
   while( hash_table[loc] && !found_it)
      {
      n_probes++;
      found_it = (hash_table[loc]->planet_no == planet_no &&
                               hash_table[loc]->jd == jd);
      if( !found_it)
         loc = (loc + n_probes) % table_size;
      }

#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
   if( max_probes_required < n_probes)
      max_probes_required = n_probes;
   total_n_searches++;
   total_n_probes += n_probes;
#endif

               /* If we _didn't_ find such an entry,  then hash_table[loc] */
               /* must be null.  Set it to a valid entry,  then fill it    */
               /* with the correct planetary position.                     */
   if( !found_it)
      {
      hash_table[loc] = cache + n_posns_cached;
      rval = planet_posn_raw( planet_no, jd, hash_table[loc]->vect);
      hash_table[loc]->planet_no = planet_no;
      hash_table[loc]->jd = jd;
      n_posns_cached++;
      }
   memcpy( vect_2000, hash_table[loc]->vect, 3 * sizeof( double));
   return( rval);
}

      /* In the following,  we get the earth's position for a particular    */
      /* instant,  just to ensure that JPL ephemerides (if any) are loaded. */
      /* Then we call with planet = JD = 0,  which causes the info about    */
      /* the JPL ephemerides to be put into the 'state vector'.             */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end)
{
   double vect_2000[3];

   planet_posn_raw( 3, 2451545., vect_2000);
   planet_posn_raw( 0, 0., vect_2000);
   *de_version = (int)vect_2000[0];
   if( jd_start)
      *jd_start = (int)vect_2000[1];
   if( jd_end)
      *jd_end = (int)vect_2000[2];
   return( 0);
}

int format_jpl_ephemeris_info( char *buff)
{
   int de_version;
   double jd_start, jd_end;
   const double j2000 = 2451545.;

   get_jpl_ephemeris_info( &de_version, &jd_start, &jd_end);
   if( !de_version && !jd_start && !jd_end)
      {
      strcpy( buff, "\nNo JPL DE ephemeris file loaded; using (slower) PS1996 series\n");
      strcat( buff, "See http://www.projectpluto.com/find_orb.htm#de_eph for\n");
      strcat( buff, "info on how/why to use JPL DE ephemerides\n");
      }
   else
      sprintf( buff,
            "\nUsing DE-%d; covers years %.1lf to %.1lf\n", de_version,
            (jd_start - j2000) / 365.25 + 2000.,
            (jd_end   - j2000) / 365.25 + 2000.);
               /* Kludge to allow re-use of 'version' string */
               /* between 'about' and 'main' dialogues:      */
   return( de_version);
}
