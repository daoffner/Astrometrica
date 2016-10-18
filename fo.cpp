/* fo.cpp: main driver for non-interactive Find_Orb

Copyright (C) 2012, Project Pluto

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "watdefs.h"
#include "weight.h"
#include "afuncs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "date.h"
#include "monte0.h"

int debug_level = 0;

int sanity_test_observations( const char *filename);
int debug_printf( const char *format, ...);                /* runge.cpp */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */

/* In this non-interactive version of Find_Orb,  we just print out warning
messages such as "3 observations were made in daylight" or "couldn't find
thus-and-such file" and proceed as if nothing had happened: */

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
   printf( "%s\n", prompt);
   return( 0);
}

static void object_comment_text( char *buff, const OBJECT_INFO *id)
{
   sprintf( buff, "%d observations; ", id->n_obs);
   make_date_range_text( buff + strlen( buff), id->jd_start, id->jd_end);
}

static OBSERVE FAR *load_object( FILE *ifile, OBJECT_INFO *id,
                       double *curr_epoch, double *epoch_shown, double *orbit)
{
   extern int n_obs_actually_loaded;
   extern int perturbers;
   int got_vector;
   OBSERVE FAR *obs = load_observations( ifile, id->packed_desig,
                                                id->n_obs);

   if( debug_level || n_obs_actually_loaded != id->n_obs)
      {
      char buff[80];

      debug_printf( " %d observations loaded\n", n_obs_actually_loaded);
      make_date_range_text( buff, obs[0].jd,
                                    obs[n_obs_actually_loaded - 1].jd);
      id->n_obs = n_obs_actually_loaded;
      debug_printf( "%s\n", buff);
      }
   obj_desig_to_perturber( id->packed_desig);
   got_vector = fetch_previous_solution( obs, id->n_obs, orbit,
                            curr_epoch, &perturbers);
   *epoch_shown = (got_vector ? *curr_epoch : floor( *curr_epoch) + .5);
   return( obs);
}

/* In the (interactive) console Find_Orb,  these allow some functions in
orb_func.cpp to show info as orbits are being computed.  In this
non-interactive code,  they're mapped to do nothing. */

void refresh_console( void)
{
}

void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes)
{
}

int main( const int argc, const char **argv)
{
   char tbuff[300];
   int n_ids, i;
   OBJECT_INFO *ids;
   int heliocentric_only = 1;

   for( i = 1; i < argc; i++)       /* check to see if we're debugging: */
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'd':
               debug_level = atoi( argv[i] + 2);
               if( !debug_level)
                  debug_level = 1;
               debug_printf( "fo: debug_level = %d; %s %s\n",
                           debug_level, __DATE__, __TIME__);
               break;
            case 'c':
               {
               extern int combine_all_observations;

               combine_all_observations = 1;
               }
               break;
            case 'h':                     /* show planet-centric orbits */
               heliocentric_only = 0;
               break;
            case 'm':
               {
               extern int integration_method;

               integration_method = atoi( argv[i] + 2);
               }
               break;
            case 'a':
               {
               extern int separate_periodic_comet_apparitions;

               separate_periodic_comet_apparitions ^= 1;
               }
            case 's':
               sanity_test_observations( argv[1]);
               printf( "Sanity check complete\n");
               return( 0);
               break;
            default:
               printf( "Unknown command-line option '%s'\n", argv[i]);
               return( -1);
            }

   load_up_weight_records( "weight.txt");
   if( debug_level)
      debug_printf( "Default weighting table read\n");

   if( argc < 2)
      {
      printf( "'fo' needs the name of an input file of MPC-formatted\n");
      printf( "astrometry as a command-line argument.\n");
      return( -2);
      }

   ids = find_objects_in_file( argv[1], &n_ids, NULL);
   if( n_ids <= 0)
      {        /* no objects found,  or file not found */
      const char *err_msg;

      if( n_ids == -1)
         err_msg = "Couldn't locate the file '%s'";
      else
         err_msg = "No objects found in file '%s'";
      printf( err_msg, argv[1]);
      return( -1);
      }

// if( !interactive_mode)
      {
      FILE *ifile = fopen( argv[1], "rb");

      printf( "Processing %d objects\n", n_ids);
      for( i = 0; i < n_ids; i++)
         {
         const char *orbit_constraints = "";
         OBSERVE FAR *obs;
         const int n_obs = ids[i].n_obs;

         object_comment_text( tbuff, ids + i);
                  /* Abbreviate 'observations:' to 'obs:' */
         text_search_and_replace( tbuff, "ervations", "");
         printf( " %d: %s: %s", i + 1, ids[i].obj_name, tbuff);
         if( n_obs < 2)
            printf( "; skipping\n");
         else
            {
            extern int append_elements_to_element_file;
            extern char orbit_summary_text[];
            long file_offset = ids[i].file_offset - 40L;
            const int element_precision = 5;
            double epoch_shown, curr_epoch, orbit[6];

                /* Start a bit ahead of the actual data,  just in case */
                /* there's a #Weight: or similar command in there: */
            if( file_offset < 0L)
               file_offset = 0L;
            fseek( ifile, file_offset, SEEK_SET);
            obs = load_object( ifile, ids + i, &curr_epoch, &epoch_shown, orbit);

            if( i)
               append_elements_to_element_file = 1;
            create_obs_file( obs, n_obs, i);
            write_out_elements_to_file( orbit, curr_epoch, epoch_shown,
                  obs, n_obs, orbit_constraints, element_precision,
                  0, heliocentric_only);
            printf( "; %s\n", orbit_summary_text);
            }
         }
      fclose( ifile);
      }
                   {
                   extern clock_t t_transfer;

                   printf( "%.2lf seconds in transfer orbits\n",
                            (double)t_transfer / (double)CLOCKS_PER_SEC);
                   }
   return( 0);
}
