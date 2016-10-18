/* findorb.cpp: main driver for console Find_Orb

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

#ifdef _WIN32
#define PURE_WINDOWS_VERSION  1
#endif

#ifdef PURE_WINDOWS_VERSION
#include <windows.h>
#endif

/* "mycurses" is a minimal implementation of Curses for DOS,  containing
   just what's needed for a few of my DOS apps (including this one).  It's
   currently used only in the OpenWATCOM implementation,  and even there,
   it's not turned on by default. */

#ifdef USE_MYCURSES
   #include "mycurses.h"
   #include "bmouse.h"
   #include <conio.h>
BMOUSE global_bmouse;
#else
            /* If we're using PDCurses,  I currently use the DLL version.  */
   #define PDC_DLL_BUILD
   #ifdef MOUSE_MOVED
      #undef MOUSE_MOVED
   #endif
   #ifdef COLOR_MENU
      #undef COLOR_MENU
   #endif
   #include "curses.h"
#endif
      /* The 'usual' Curses library provided with Linux lacks a few things */
      /* that PDCurses and MyCurses have, such as definitions for ALT_A    */
      /* and such.  'curs_lin.h' fills in these gaps.   */
#ifndef ALT_A
   #include "curs_lin.h"
#endif

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

extern int perturbers;

#define KEY_MOUSE_MOVE 31000
#define KEY_TIMER      31001
#define AUTO_REPEATING 31002

#define CTRL(c) ((c) & 0x1f)

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

#define RESIDUAL_FORMAT_SHOW_DELTAS              64

int simplex_method( OBSERVE FAR *obs, int n_obs, double *orbit,
               const double r1, const double r2);          /* orb_func.cpp */
int superplex_method( OBSERVE FAR *obs, int n_obs, double *orbit);
static void show_a_file( const char *filename);
static void put_colored_text( const char *text, const int line_no,
               const int column, const int n_bytes, const int color);
int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
             const double r1, const double angle_param);   /* orb_func.cpp */
int find_nth_sr_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                            const int orbit_number);       /* orb_func.cpp */
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile);
int debug_printf( const char *format, ...);                /* runge.cpp */
static void get_mouse_data( int *mouse_x, int *mouse_y, int *mouse_z, unsigned long *button);
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                              /* ephem0.cpp */
int store_defaults( const int ephemeris_output_options);    /* elem_out.c */
int get_defaults( int *ephemeris_output_options);           /* elem_out.c */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
int sort_obs_by_date_and_remove_duplicates( OBSERVE *obs, int n_obs);
int create_b32_ephemeris( const char *filename, const double epoch,
                const double *orbit, const int n_steps,         /* b32_eph.c */
                const double ephem_step, const double jd_start);
void put_observer_data_in_text( const char FAR *mpc_code, char *buff);
                                                            /* mpc_obs.c */
int filter_obs( OBSERVE FAR *obs, const int n_obs,
                  const double max_residual_in_arcseconds);
int find_precovery_plates( const char *filename, const double *orbit,
                           double epoch_jd);          /* ephem0.cpp */
void set_statistical_ranging( const int new_using_sr);      /* elem_out.cpp */
int link_arcs( OBSERVE *obs, int n_obs, const double r1, const double r2);
int find_circular_orbits( OBSERVE FAR *obs1, OBSERVE FAR *obs2,
               double *orbit, const int desired_soln);   /* orb_fun2.cpp */
void set_up_observation( OBSERVE FAR *obs);               /* mpc_obs.cpp */
double peirce_rayleigh_func( const int N, const int n, const int m);
#ifdef _WIN32
int clipboard_to_file( const char *filename, const int append);      /* mpc_obs.cpp */
#endif

#define COLOR_BACKGROUND            1
#define COLOR_ORBITAL_ELEMENTS      2
#define COLOR_FINAL_LINE            3
#define COLOR_SELECTED_OBS          4
#define COLOR_HIGHLIT_BUTTON        5
#define COLOR_EXCLUDED_AND_SELECTED 8
#define COLOR_EXCLUDED_OBS          6
#define COLOR_OBS_INFO              7
#define COLOR_MESSAGE_TO_USER       8
#define COLOR_RESIDUAL_LEGEND       9
#define COLOR_MENU                 10
#define COLOR_DEFAULT_INQUIRY       9
#define COLOR_ATTENTION             7
#define COLOR_SCROLL_BAR           11


#define COLOR_GRAY                  9

#ifdef USE_MYCURSES
static int curses_kbhit( )
{
   return( kbhit( ) ? 0: ERR);
}
#else
static int curses_kbhit( )
{
   int c;

   nodelay( stdscr, TRUE);
   c = getch( );
   nodelay( stdscr, FALSE);
   if( c != ERR)     /* no key waiting */
      ungetch( c);
   return( c);
}
#endif

static int extended_getch( void)
{
#ifdef _WIN32
   int rval = getch( );

   if( !rval)
      rval = 256 + getch( );
#endif
#ifdef __WATCOMC__
   int rval = 0;

   while( !rval)
      {
      if( !curses_kbhit( ))
         {
         rval = getch( );
         if( !rval)
            rval = 256 + getch( );
         }
#ifdef USE_MYCURSES
      else
         {
         mouse_read( &global_bmouse);
         if( global_bmouse.released)
            rval = KEY_MOUSE;
         }
#endif
      if( !rval)
         napms( 200);
      }
#endif
#if !defined (_WIN32) && !defined( __WATCOMC__)
   int rval = getch( );

   if( rval == 27)
      {
      clock_t end_time = clock( ) + CLOCKS_PER_SEC / 2;

      nodelay( stdscr, TRUE);
      do
         {
         rval = getch( );
         }
         while( rval == ERR && clock( ) < end_time);
      nodelay( stdscr, FALSE);
      if( rval == ERR)    /* no second key found */
         rval = 27;       /* just a plain ol' Escape */
      else
         rval += (ALT_A - 'a');
      }
#endif
   return( rval);
}

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
   int i, j, rval, line, col, n_lines = 1, line_start = 0, box_size = 0;
   const int side_borders = 1;   /* leave a blank on either side */
   int real_width;
   char tbuff[200];
   chtype *buffered_screen;

   for( i = 0; prompt[i]; i++)
      if( prompt[i] == '\n')
         {
         if( box_size < i - line_start)
            box_size = i - line_start;
         line_start = i;
         n_lines++;
         }
   if( box_size < i - line_start)
      box_size = i - line_start;
   if( box_size > getmaxx( stdscr) - 2)
      box_size = getmaxx( stdscr) - 2;

   line = (getmaxy( stdscr) - n_lines) / 2;
   col = (getmaxx( stdscr) - box_size) / 2;
   real_width = side_borders * 2 + box_size;
   tbuff[real_width] = '\0';
         /* Store rectangle behind the 'inquiry box': */
   buffered_screen = (chtype *)calloc( n_lines * real_width,
                     sizeof( chtype));
   for( i = 0; i < n_lines; i++)
      for( j = 0; j < real_width; j++)
         buffered_screen[j + i * real_width] =
                  mvinch( line + i, col - side_borders + j);
   for( i = 0; prompt[i]; )
      {
      int n_spaces, color_to_use = color;

      for( j = i; prompt[j] && prompt[j] != '\n'; j++)
         ;
      memset( tbuff, ' ', side_borders);
      memcpy( tbuff + side_borders, prompt + i, j - i);
      n_spaces = box_size + side_borders - (j - i);
      if( n_spaces > 0)
         memset( tbuff + side_borders + j - i, ' ', n_spaces);
      if( !i)
         color_to_use |= 0x4000;     /* overline */
      if( !prompt[j] || !prompt[j + 1])
         color_to_use |= 0x2000;     /* underline */
      put_colored_text( tbuff, line, col - side_borders,
             real_width, color_to_use);
      i = j;
      if( prompt[i] == '\n')
         i++;
      line++;
      }
   if( buff)
      {
      memset( tbuff, ' ', real_width);
      put_colored_text( tbuff, line, col - side_borders,
             real_width, color);
      move( line, col);
      echo( );
      refresh( );
      rval = getnstr( buff, max_len);
      noecho( );
      }
   else
      {
      refresh( );
      rval = extended_getch( );
      }
   line -= n_lines;     /* put back to top of box */
   for( i = 0; i < n_lines; i++)
      mvaddchnstr( line + i, col - side_borders,
            buffered_screen + i * real_width, real_width);
   free( buffered_screen);
   return( rval);
}

/* In the (interactive) console Find_Orb,  these allow some functions
in orb_func.cpp to show info as orbits are being computed.  In the
non-interactive 'fo' code,  they're mapped to do nothing (see 'fo.cpp'). */

void refresh_console( void)
{
   refresh( );
}

void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes)
{
   mvaddnstr( col, row, msg, n_bytes);
}


static int extract_date( const char *buff, double *jd)
{
   int rval = 0;

               /* If the date seems spurious,  use 'now' as our zero point: */
   if( *jd < 10000.)
      {
      static const double jan_1970 = 2440587.5;
      const double jd_now = jan_1970 + time( NULL) / 86400.;

      *jd = jd_now;
      }
   *jd = get_time_from_string( *jd, buff,
                     FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
   rval = 2;
   if( *jd == 0.)
      rval = -1;
   return( rval);
}

/* Here's a simplified example of the use of the 'ephemeris_in_a_file'
   function... nothing fancy,  but it shows how it's used.  */

static char mpc_code[8];
static char ephemeris_start[80], ephemeris_step_size[80];
static int ephemeris_output_options, n_ephemeris_steps;

static void create_ephemeris( const double *orbit, const double epoch_jd,
         OBSERVE *obs, const int n_obs)
{
   int c = 1;
   char buff[1000];
   double jd_start = 0., jd_end = 0., step = 0.;

   while( c > 0)
      {
      int format_start;
      const int ephem_type = (ephemeris_output_options & 7);
      const char *ephem_type_strings[] = {
               "Observables",
               "State vectors",
               "Cartesian coord positions",
               "MPCORB output",
               "8-line elements",
               "Close approaches",
               NULL };

      jd_start = 0.;
      format_start = extract_date( ephemeris_start, &jd_start);
      step = get_step_size( ephemeris_step_size, NULL, NULL);
      if( format_start == 1 || format_start == 2)
         {
         if( step && format_start == 1)  /* time was relative to 'right now' */
            jd_start = floor( (jd_start - .5) / step) * step + .5;
         sprintf( buff, " (Ephem start: JD %.5lf = ", jd_start);
         full_ctime( buff + strlen( buff), jd_start,
                    FULL_CTIME_DAY_OF_WEEK_FIRST | CALENDAR_JULIAN_GREGORIAN);
         strcat( buff, ")\n");
         jd_end = jd_start + step * (double)n_ephemeris_steps;
         sprintf( buff + strlen( buff), " (Ephem end:   JD %.5lf = ", jd_end);
         full_ctime( buff + strlen( buff), jd_end,
                    FULL_CTIME_DAY_OF_WEEK_FIRST | CALENDAR_JULIAN_GREGORIAN);
         strcat( buff, ")\n");
         }
      else
         {
         strcpy( buff, "(Ephemeris starting time isn't valid)\n");
         jd_start = jd_end = 0.;
         }

      sprintf( buff + strlen( buff), "T  Ephem start: %s\n", ephemeris_start);
      sprintf( buff + strlen( buff), "N  Number steps: %d\n",
                                        n_ephemeris_steps);
      sprintf( buff + strlen( buff), "S  Step size: %s\n", ephemeris_step_size);
      sprintf( buff + strlen( buff), "L  Location: (%s) ", mpc_code);
      put_observer_data_in_text( mpc_code, buff + strlen( buff));
      strcat( buff, "\n");
      if( ephem_type == OPTION_OBSERVABLES)    /* for other tables,        */
         {                          /* these options are irrelevant:       */
         sprintf( buff + strlen( buff), "Z  Motion info in ephemerides: %s\n",
                  (ephemeris_output_options & OPTION_MOTION_OUTPUT) ? "On" : "Off");
         sprintf( buff + strlen( buff), "A  Alt/az info in ephemerides: %s\n",
                  (ephemeris_output_options & OPTION_ALT_AZ_OUTPUT) ? "On" : "Off");
         sprintf( buff + strlen( buff), "R  Radial velocity in ephemerides: %s\n",
                  (ephemeris_output_options & OPTION_RADIAL_VEL_OUTPUT) ? "On" : "Off");
         sprintf( buff + strlen( buff), "P  Phase angle in ephemerides: %s\n",
                  (ephemeris_output_options & OPTION_PHASE_ANGLE_OUTPUT) ? "On" : "Off");
         sprintf( buff + strlen( buff), "B  Phase angle bisector: %s\n",
                  (ephemeris_output_options & OPTION_PHASE_ANGLE_BISECTOR) ? "On" : "Off");
         sprintf( buff + strlen( buff), "H  Heliocentric ecliptic: %s\n",
                  (ephemeris_output_options & OPTION_HELIO_ECLIPTIC) ? "On" : "Off");
         sprintf( buff + strlen( buff), "X  Topocentric ecliptic: %s\n",
                  (ephemeris_output_options & OPTION_TOPO_ECLIPTIC) ? "On" : "Off");
         sprintf( buff + strlen( buff), "G  Ground track: %s\n",
                  (ephemeris_output_options & OPTION_GROUND_TRACK) ? "On" : "Off");
         if( ephemeris_output_options & OPTION_MOTION_OUTPUT)
            sprintf( buff + strlen( buff), "E  Separate motions: %s\n",
                  (ephemeris_output_options & OPTION_SEPARATE_MOTIONS) ? "On" : "Off");
         }
      sprintf( buff + strlen( buff), "C  %s\n", ephem_type_strings[ephem_type]);
      sprintf( buff + strlen( buff), "?  Help about making ephemerides\n");
      sprintf( buff + strlen( buff), "Hit 'm' to make ephemeris, 'q' to return to main display");
      c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
      switch( c)
         {
         case 't': case 'T':
            inquire( "Enter start of ephemeris (YYYY MM DD, or JD, or 'now'):",
                  ephemeris_start, sizeof( ephemeris_start), COLOR_MESSAGE_TO_USER);
            break;
         case 'n': case 'N':
            inquire( "Number of steps:", buff, sizeof( buff), COLOR_MESSAGE_TO_USER);
            if( atoi( buff) > 0)
               n_ephemeris_steps = atoi( buff);
            break;
         case 's': case 'S':
            inquire( "Enter step size in days: ",
                  ephemeris_step_size, sizeof( ephemeris_step_size), COLOR_MESSAGE_TO_USER);
            break;
         case 'l': case 'L':
            inquire( "Enter MPC code: ", buff, sizeof( buff), COLOR_MESSAGE_TO_USER);
            if( strlen( buff) == 3)
               strcpy( mpc_code, buff);
            break;
         case 'a': case 'A':
            ephemeris_output_options ^= OPTION_ALT_AZ_OUTPUT;
            break;
         case 'b': case 'B':
            ephemeris_output_options ^= OPTION_PHASE_ANGLE_BISECTOR;
            break;
         case 'h': case 'H':
            ephemeris_output_options ^= OPTION_HELIO_ECLIPTIC;
            break;
         case 'x': case 'X':
            ephemeris_output_options ^= OPTION_TOPO_ECLIPTIC;
            break;
         case 'e': case 'E':
            ephemeris_output_options ^= OPTION_SEPARATE_MOTIONS;
            break;
         case 'c': case 'C':
            if( ephem_type == OPTION_CLOSE_APPROACHES)  /* end of cycle: */
               ephemeris_output_options -= OPTION_CLOSE_APPROACHES;
            else                    /* not at end: move forward */
               ephemeris_output_options++;
            break;
         case 'g': case 'G':
            ephemeris_output_options ^= OPTION_GROUND_TRACK;
            break;
         case 'p': case 'P':
            ephemeris_output_options ^= OPTION_PHASE_ANGLE_OUTPUT;
            break;
         case 'r': case 'R':
            ephemeris_output_options ^= OPTION_RADIAL_VEL_OUTPUT;
            break;
         case 'z': case 'Z':
            ephemeris_output_options ^= OPTION_MOTION_OUTPUT;
            break;
         case 'm': case 'M':
            {
            const char *err_msg = NULL;
            const double min_jd = 10000.;  /* = year -4685,  i.e., */
                                           /* obviously wrong */

            if( jd_start < min_jd)
               err_msg = "You need to set a valid starting date!";
            else if( !n_ephemeris_steps)
               err_msg = "You need to set the number of ephemeris steps!";
            else if( !step)
               err_msg = "You need to set a valid step size!";
            else                 /* yes,  we can make an ephemeris */
               c = -2;
            if( err_msg)
               inquire( err_msg, NULL, 0, COLOR_FINAL_LINE);
            }
            break;
         case 'q': case 'Q':
            c = -1;
            break;
#ifdef KEY_EXIT
         case KEY_EXIT:
            exit( 0);
            break;
#endif
         default:
            show_a_file( "dosephem.txt");
            break;
         }
      }

   if( c == -2)         /* yes,  we're making an ephemeris */
      {
      if( !strcmp( mpc_code, "32b"))
         {
         inquire( ".b32 filename:",
                               buff, sizeof( buff), COLOR_DEFAULT_INQUIRY);
         if( *buff)
            {
            strcat( buff, ".b32");
            create_b32_ephemeris( buff, epoch_jd, orbit, n_ephemeris_steps,
                     atof( ephemeris_step_size), jd_start);    /* b32_eph.c */
            }
         }
      else
         {
         double lon, rho_sin_phi, rho_cos_phi;
         extern const char *ephemeris_filename;
         int planet_no = 3;

         planet_no = get_observer_data( mpc_code, buff, &lon,
                                           &rho_cos_phi, &rho_sin_phi);
         if( ephemeris_in_a_file( ephemeris_filename, orbit, obs, n_obs,
               planet_no, epoch_jd, jd_start, ephemeris_step_size, lon,
               rho_cos_phi, rho_sin_phi, n_ephemeris_steps,
               mpc_code, ephemeris_output_options))
            inquire( "Ephemeris generation failed!  Hit any key:", NULL, 0,
                              COLOR_MESSAGE_TO_USER);
         else
            show_a_file( ephemeris_filename);
         }
      }
}


static void object_comment_text( char *buff, const OBJECT_INFO *id)
{
   sprintf( buff, "%d observations; ", id->n_obs);
   make_date_range_text( buff + strlen( buff), id->jd_start, id->jd_end);
}

static int compare_by_last_obs_time = 0;

static int id_compare( const void *a, const void *b)
{
   const OBJECT_INFO *aptr = (const OBJECT_INFO *)a;
   const OBJECT_INFO *bptr = (const OBJECT_INFO *)b;
   int rval;

   if( compare_by_last_obs_time)
      rval = (aptr->jd_end > bptr->jd_end ? 1 : -1);
   else
      rval = strcmp( aptr->packed_desig, bptr->packed_desig);
   return( rval);
}

/* select_object_in_file( ) uses the find_objects_in_file( ) function to
   get a list of all the objects listed in a file of observations.  It
   prints out their IDs and asks the user to hit a key to select one.
   The name of that object is then copied into the 'obj_name' buffer.

   At present,  it's used only in main( ) when you first start findorb,
   so you can select the object for which you want an orbit.         */

int select_object_in_file( OBJECT_INFO *ids, const int n_ids)
{
   static int choice = 0, show_packed = 0;
   int rval = -1;

   qsort( ids, n_ids, sizeof( OBJECT_INFO), id_compare);
   if( ids && n_ids)
      {
      int i, curr_page = 0, err_message = 0, force_full_width_display = 0;

      clear( );
      while( rval == -1)
         {
         const int n_lines = getmaxy( stdscr) - 3;
         int column_width = 16;
         int c, n_cols = getmaxx( stdscr) / column_width;
         char buff[280];

         if( force_full_width_display)
            n_cols = 1;
         if( choice < 0)
            choice = 0;
         if( choice > n_ids - 1)
            choice = n_ids - 1;
         while( curr_page > choice)
            curr_page -= n_lines;
         while( curr_page + n_lines * n_cols <= choice)
            curr_page += n_lines;
         if( curr_page < 0)      /* ensure that we wrap around: */
            curr_page = 0;
         if( n_ids < n_lines * n_cols)
            n_cols = n_ids / n_lines + 1;
         column_width = getmaxx( stdscr) / n_cols;
//       if( column_width > 80)
//          column_width = 80;
         for( i = 0; i < n_lines * n_cols; i++)
            {
            char desig[181];
            int color = COLOR_BACKGROUND;

            if( i + curr_page < n_ids)
               {
               if( show_packed)
                  sprintf( desig, "'%s'", ids[i + curr_page].packed_desig);
               else
                  strcpy( desig, ids[i + curr_page].obj_name);
               if( i + curr_page == choice)
                  {
                  sprintf( buff, "Object %d of %d: %s",
                              choice + 1, n_ids, desig);
                  put_colored_text( buff, n_lines + 1, 0, -1,
                                                COLOR_SELECTED_OBS);
                  object_comment_text( buff, ids + choice);
                  put_colored_text( buff, n_lines + 2, 0, -1,
                                                COLOR_SELECTED_OBS);
                  color = COLOR_HIGHLIT_BUTTON;
                  }
               else
                  if( ids[i + curr_page].solution_exists)
                     color = COLOR_OBS_INFO;
               }
            else                        /* just want to erase space: */
               *desig = '\0';
            desig[column_width] = '\0';    /* trim to fit available space */
            sprintf( buff, "%-*s", column_width, desig);
            if( n_cols == 1 && i + curr_page < n_ids)
               object_comment_text( buff + 25, ids + i + curr_page);
            put_colored_text( buff, i % n_lines,
                   (i / n_lines) * column_width,
                   (n_cols == 1 ? -1 : column_width), color);
            }

         put_colored_text( "Use arrow keys to find your choice,  then hit the space bar or Enter",
                                 n_lines, 0, -1, COLOR_SELECTED_OBS);
         if( err_message)
            put_colored_text( "Not a valid choice",
                                 n_lines + 2, 0, -1, COLOR_FINAL_LINE);
         put_colored_text( "Quit", n_lines + 2, 75, 4, COLOR_HIGHLIT_BUTTON);
         put_colored_text( "Next", n_lines + 2, 70, 4, COLOR_HIGHLIT_BUTTON);
         put_colored_text( "Prev", n_lines + 2, 65, 4, COLOR_HIGHLIT_BUTTON);
         put_colored_text( "End", n_lines + 2, 61, 3, COLOR_HIGHLIT_BUTTON);
         put_colored_text( "Start", n_lines + 2, 55, 5, COLOR_HIGHLIT_BUTTON);
         refresh( );
         flushinp( );
         c = extended_getch( );
         err_message = 0;
         if( c == KEY_MOUSE)
            {
            int x, y, z;
            unsigned long button;

            get_mouse_data( &x, &y, &z, &button);
#ifdef BUTTON5_PRESSED
            if( button & BUTTON4_PRESSED)   /* actually 'wheel up' */
               c = KEY_UP;
            else if( button & BUTTON5_PRESSED)   /* actually 'wheel down' */
               c = KEY_DOWN;
            else
#endif
              if( y < n_lines)
               choice = curr_page + y + (x / column_width) * n_lines;
            else if( y == n_lines + 2)
               {
               if( x >= 75)
                  c = 27;          /* quit */
               else if( x >= 70)
                  c = KEY_NPAGE;   /* 'next page' */
               else if( x >= 65)
                  c = KEY_PPAGE;   /* 'prev page' */
               else if( x >= 61)
                  c = KEY_END;     /* end of list */
               else if( x >= 55)
                  c = KEY_HOME;    /* start of list */
               }
#ifndef USE_MYCURSES
            if( button & BUTTON1_DOUBLE_CLICKED)
               rval = choice;
#endif
            }
                     /* if a letter/number is hit,  look for an obj that */
                     /* starts with that letter/number: */
         if( c > ' ' && c <= 'z' && isalnum( c))
            for( i = 1; i < n_ids; i++)
               {
               const int loc = (i + choice) % n_ids;

               if( toupper( c) == toupper( ids[loc].obj_name[0]))
                  {
                  choice = loc;
                  i = n_ids;
                  }
               }
         else switch( c)
            {
            case 9:
               force_full_width_display ^= 1;
               break;
            case ';':
               compare_by_last_obs_time ^= 1;
               qsort( ids, n_ids, sizeof( OBJECT_INFO), id_compare);
               break;
            case ' ':
            case 13:
               rval = choice;
               break;
#ifdef KEY_C2
            case KEY_C2:
#endif
            case KEY_DOWN:
               choice++;
               break;
#ifdef KEY_A2
            case KEY_A2:
#endif
            case KEY_UP:
               choice--;
               break;
#ifdef KEY_B1
            case KEY_B1:
#endif
            case KEY_LEFT:
               choice -= n_lines;
               break;
#ifdef KEY_B3
            case KEY_B3:
#endif
            case KEY_RIGHT:
               choice += n_lines;
            break;
            case KEY_C3:         /* "PgDn" = lower right key in keypad */
            case KEY_NPAGE:
            case 'n':
               choice += n_lines * n_cols;
               break;
            case KEY_A3:         /* "PgUp" = upper right key in keypad */
            case KEY_PPAGE:
            case 'p':
               choice -= n_lines * n_cols;
               break;
            case KEY_C1:
            case KEY_END:
            case 'e':
               choice = n_ids - 1;
               break;
            case KEY_A1:
            case KEY_HOME:
            case 's':
               choice = 0;
               break;
            case ',':
               show_packed ^= 1;
               break;
            case KEY_MOUSE:      /* already handled above */
               break;
#ifdef KEY_RESIZE
            case KEY_RESIZE:
               resize_term( 0, 0);
               break;
#endif
            case 'q': case 'Q': case 27:
#ifdef KEY_EXIT
            case KEY_EXIT:
#endif
               rval = -2;
               break;
            }
         }
      if( debug_level > 3)
         debug_printf( "rval = %d; leaving select_object_in_file\n", rval);
      }
   return( rval);
}

void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */

static const char *command_text =
          "help Full Herget epheM Vaisa resiD Gauss constr New Quit";
static const char *command_remap =
          "???? ffff hhhhhh mmmmm vvvvv ddddd ggggg llllll nnn qqqqqqqqqqqq";

void show_basic_info( const OBSERVE FAR *obs, const int n_obs)
{
   char buff[80];
   double r1, r2;
   int i;

   get_r1_and_r2( n_obs, obs, &r1, &r2);    /* orb_func.cpp */
   strcpy( buff, "R1:");
   format_dist_in_buff( buff + 3, r1);
   put_colored_text( buff, 0, 0, 15, COLOR_BACKGROUND);

   strcpy( buff, "  R2:");
   format_dist_in_buff( buff + 5, r2);
   put_colored_text( buff, 0, 10, -1, COLOR_BACKGROUND);

   for( i = 0; command_text[i]; i++)
      if( command_text[i] != ' ')
         put_colored_text( (char *)command_text + i, 0, i + 24, 1,
                                 COLOR_MENU);
}

static const char *perturber_names[] = {
          "Mercury", "Venus", "Earth", "Mars", "Jupiter",
          "Saturn", "Uranus", "Neptune", "Pluto", "Moon", "Asteroids", NULL };

void show_perturbers( void)
{
   int i;

   for( i = 0; perturber_names[i]; i++)
      {
      int color = COLOR_BACKGROUND;
      const int shift_amt = (i == 10 ? 20 : i + 1);
      char buff[20];

      strcpy( buff, "(o)");
      if( (perturbers >> shift_amt) & 1)
         color = COLOR_HIGHLIT_BUTTON;
      else
         {
         buff[1] = (char)( '0' + (i + 1) % 10);
         if( i == 10)
            buff[1] = 'a';
         }
      put_colored_text( buff, 1, i * 7, 3, color);
      strcpy( buff, perturber_names[i]);
      strcpy( buff + 3, " ");
      put_colored_text( buff, 1, i * 7 + 3, strlen( buff), COLOR_BACKGROUND);
      }
               /* Clear to end of line: */
   put_colored_text( " ", 1, i * 7, -1, COLOR_BACKGROUND);
}

int max_mpc_color_codes = 5;
static MPC_STATION *mpc_color_codes = NULL;

/* Show the text for a given observation... which will be in
   'default_color' unless it's excluded.  If it is,  the residuals
   are shown in COLOR_EXCLUDED_OBS.  The three-character MPC code
   is also shown in a separate color. */

int mpc_column, resid_column;

static void show_residual_text( char *buff, const int line_no,
           const int column, const int default_color,
           const int is_included)
{
   put_colored_text( buff, line_no, column, strlen( buff), default_color);
   if( !is_included)
      {
      const int residual_field_size = 13;
      char tbuff[40];

      memcpy( tbuff, buff + resid_column - 2, residual_field_size);
      tbuff[0] = '(';                       /* put ()s around excluded obs */
      tbuff[residual_field_size - 1] = ')';
      tbuff[residual_field_size] = '\0';
      put_colored_text( tbuff, line_no, column + resid_column - 2,
               residual_field_size,
               (default_color == COLOR_SELECTED_OBS ?
                     COLOR_EXCLUDED_AND_SELECTED : COLOR_EXCLUDED_OBS) + 4096);
      }
   if( mpc_column >= 0)
      {
      buff += mpc_column;
      if( *buff != ' ')       /* show MPC code in color: */
         put_colored_text( buff, line_no, column + mpc_column, 3,
                        512 + 16 + find_mpc_color( mpc_color_codes, buff));
      }
}

#define SORT_BY_SCORE 0
#define SORT_BY_NAME  1

static void sort_mpc_codes( const int n_to_sort, const int sort_flag)
{
   int i, do_swap;

   for( i = 0; i < n_to_sort - 1; i++)
      {
      if( sort_flag == SORT_BY_SCORE)
         do_swap = (mpc_color_codes[i].score < mpc_color_codes[i + 1].score);
      else
         do_swap = (strcmp( mpc_color_codes[i].code,
                            mpc_color_codes[i + 1].code) > 0);
      if( do_swap)
         {
         MPC_STATION temp = mpc_color_codes[i];

         mpc_color_codes[i] = mpc_color_codes[i + 1];
         mpc_color_codes[i + 1] = temp;
         if( i)
            i -= 2;
         }
      }
}

static void add_to_mpc_color( const char *code, const int score_increment)
{
   int i;

   for( i = 0; mpc_color_codes[i].code[0]; i++)
      if( !strcmp( mpc_color_codes[i].code, code))
         mpc_color_codes[i].score += score_increment;
}

static void show_right_hand_scroll_bar( const int line_start,
      const int lines_to_show, const int first_line,
      const int n_lines)
{
   if( lines_to_show < n_lines)
      {
      int i;
      const int scroll0 =        first_line            * (lines_to_show - 2) / n_lines + 1;
      const int scroll1 = (first_line + lines_to_show) * (lines_to_show - 2) / n_lines + 1;

      for( i = 0; i < lines_to_show; i++)
         {
         int color = COLOR_SCROLL_BAR;
         const char *text = "|";

         if( !i || i == lines_to_show - 1)
            {
            color = COLOR_OBS_INFO;
            text = (i ? "v" : "^");
            }
         else
            if( i >= scroll0 && i <= scroll1)
               color = COLOR_HIGHLIT_BUTTON;
         put_colored_text( text, i + line_start, getmaxx( stdscr) - 1, -1,
                              color);
         }
      }
}

int first_residual_shown, n_stations_shown;

void show_residuals( const OBSERVE FAR *obs, const int n_obs,
              const int residual_format, const int curr_obs,
              const int top_line_residual_area,
              const int list_codes)
{
   int i, line_no = top_line_residual_area;
   int n_obs_shown = getmaxy( stdscr) - line_no;
   const int n_mpc_codes = find_mpc_color( mpc_color_codes, NULL);
   const int base_format = (residual_format & 3);
   char buff[120];

   n_stations_shown = n_obs_shown / 3;
   if( n_obs_shown < 0)
      return;
   for( i = 0; i < n_obs_shown; i++)         /* clear out space */
      put_colored_text( "", i + line_no, 0, -1, COLOR_BACKGROUND);
                  /* set 'scores' for each code to equal # of times */
                  /* that code was used: */
   for( i = 0; i < n_mpc_codes; i++)
      mpc_color_codes[i].score = 0;
   for( i = 0; i < n_obs; i++)
      add_to_mpc_color( obs[i].mpc_code, 1);

   if( n_stations_shown > n_mpc_codes)
      n_stations_shown = n_mpc_codes;
   if( !list_codes || n_stations_shown < 1)
      n_stations_shown = 1;
   n_obs_shown -= n_stations_shown;

   if( base_format != RESIDUAL_FORMAT_SHORT)    /* one residual/line */
      {
      int line_start = curr_obs - n_obs_shown / 2;

      if( line_start > n_obs - n_obs_shown)
         line_start = n_obs - n_obs_shown;
      if( line_start < 0)
         line_start = 0;

      for( i = 0; i < n_obs_shown; i++)
         if( line_start + i < n_obs)
            {
            int color = COLOR_BACKGROUND;

            add_to_mpc_color( obs[line_start + i].mpc_code,
                      (line_start + i == curr_obs ? n_obs * n_obs : n_obs));
            if( base_format == RESIDUAL_FORMAT_80_COL)
               {                         /* show in original 80-column MPC  */
               char resid_data[70];      /* format, w/added data if it fits */
               const int dropped_start = 12;     /* ...but omit designation */

               format_observation( obs + line_start + i, buff,
                         (residual_format & ~3) | RESIDUAL_FORMAT_HMS);
               strcpy( resid_data, buff + 42);
               recreate_observation_line( buff, obs + line_start + i);
               if( residual_format & RESIDUAL_FORMAT_HMS)
                  {
                  const long seconds = (long)( atof( buff + 25) * 86400.);
                  const char saved_byte = buff[32];

                  memmove( buff + 35, buff + 32, strlen( buff + 31));
                  sprintf( buff + 25, " %02ld:%02ld:%02ld ",
                       seconds / 3600L, (seconds / 60L) % 60L, seconds % 60L);
                  buff[35] = saved_byte;
                  }
               memmove( buff, buff + dropped_start, strlen( buff + dropped_start) + 1);
               strcpy( buff + strlen( buff), resid_data + 10);
               }
            else
               format_observation( obs + line_start + i, buff, residual_format);
            if( residual_format & RESIDUAL_FORMAT_SHOW_DELTAS)
               if( line_start + i != curr_obs)
                  {
                  double diff;
                  int column;

                  column = ((base_format == RESIDUAL_FORMAT_80_COL) ?
                                 15 : 0);
                  sprintf( buff + column, "%16.5lf",
                           obs[line_start + i].jd - obs[curr_obs].jd);
                  buff[column + 16] = ' ';

                  diff = obs[line_start + i].ra - obs[curr_obs].ra;
                  if( diff > PI)
                     diff -= PI + PI;
                  if( diff < -PI)
                     diff += PI + PI;
                  column = ((base_format == RESIDUAL_FORMAT_80_COL) ?
                                 32 : 24);
                  sprintf( buff + column, "%11.5lf", diff * 180. / PI);
                  buff[column + 11] = ' ';

                  column = ((base_format == RESIDUAL_FORMAT_80_COL) ?
                                 44 : 38);
                  diff = obs[line_start + i].dec - obs[curr_obs].dec;
                  sprintf( buff + column, "%11.5lf", diff * 180. / PI);
                  buff[column + 11] = ' ';
                  }


            if( line_start + i == curr_obs)
               color = COLOR_SELECTED_OBS;

            show_residual_text( buff, line_no++, 0, color,
                                             obs[line_start + i].is_included);
            }
      first_residual_shown = line_start;
      }
   else              /* put three observations/line */
      {
      const int n_cols = getmaxx( stdscr) / 25;
      const int n_rows = (n_obs - 1) / n_cols + 1;
      const int col_width = getmaxx( stdscr) / n_cols;
      int j;
      int line_start = curr_obs % n_rows - n_obs_shown / 2;

      if( line_start > n_rows - n_obs_shown)
         line_start = n_rows - n_obs_shown;
      if( line_start < 0)
         line_start = 0;
      first_residual_shown = line_start;

      for( i = 0; i < n_obs_shown && line_start < n_rows;
                                         i++, line_no++, line_start++)
         for( j = 0; j < n_cols; j++)
            {
            const int obs_number = j * n_rows + line_start;
            char buff[50];
            int color = COLOR_BACKGROUND, is_included = 1;

            if( obs_number < n_obs)
               {
               add_to_mpc_color( obs[obs_number].mpc_code,
                      (obs_number == curr_obs ? n_obs * n_obs : n_obs));
               buff[0] = ' ';
               format_observation( obs + obs_number, buff + 1, residual_format);
               if( obs_number == curr_obs)
                  {
                  buff[0] = '<';
                  color = COLOR_SELECTED_OBS;
                  }
               strcat( buff, (obs_number == curr_obs ? ">    " : "    "));
               is_included = obs[obs_number].is_included;
               }
            else
               memset( buff, ' ', col_width);
            buff[col_width] = '\0';
            show_residual_text( buff, line_no, j * col_width, color,
                     is_included);
            }
      }

               /* show "scroll bar" to right of observations: */
   show_right_hand_scroll_bar( top_line_residual_area,
                          n_obs_shown, first_residual_shown, n_obs);

   if( n_stations_shown)
      {
                /* OK,  sort obs by score,  highest first... */
      sort_mpc_codes( n_mpc_codes, SORT_BY_SCORE);
                /* ...then sort the ones we're going to display by name: */
      sort_mpc_codes( n_stations_shown, SORT_BY_NAME);

      for( i = 0; i < n_stations_shown; i++)
         {
         const int line_no = getmaxy( stdscr) - n_stations_shown + i;
         const int is_curr_code = !strcmp( mpc_color_codes[i].code,
                                        obs[curr_obs].mpc_code);

         sprintf( buff, "(%s)", mpc_color_codes[i].code);
         mpc_column = 1;
         show_residual_text( buff, line_no, 0, COLOR_BACKGROUND, 1);
         put_observer_data_in_text( mpc_color_codes[i].code, buff);
         mpc_column = -1;
         show_residual_text( buff, line_no, 6,
             (is_curr_code ? COLOR_FINAL_LINE : COLOR_BACKGROUND), 1);
         }
      }
}

void show_final_line( const OBSERVE FAR *obs, const int n_obs,
                                      const int curr_obs, const int color)
{
   char buff[90];
   int len;

   sprintf( buff, " %d/%d", curr_obs + 1, n_obs);
   len = strlen( buff);
   put_colored_text( buff, getmaxy( stdscr) - 1,
                              getmaxx( stdscr) - len, len, color);
}

static void show_residual_legend( const int line_no, const int residual_format)
{
   const int base_format = (residual_format & 3);
   static const char *legends[] = {
"YYYY MM DD.DDDDD    Obs   RA (J2000)     dec          Xres  Yres   delta  R",
NULL,       /* this is the 'with tabs' format used in Windows Find_Orb */
" YYMMDD Obs  Xres  Yres     ",
"   YYYY MM DD.DDDDD   RA (J2000)   dec               mag     ref Obs   Xres  Yres   delta  R" };
   char buff[290];

   strcpy( buff, legends[base_format]);
   if( residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
      {         /* residuals in time & cross-time, not RA and dec */
      char *text_loc;

      while( (text_loc = strstr( buff, "Xres")) != NULL)
         *text_loc = 'T';
      while( (text_loc = strstr( buff, "Yres")) != NULL)
         *text_loc = 'C';
      }
   if( residual_format & RESIDUAL_FORMAT_MAG_RESIDS)
      {
      if( base_format == RESIDUAL_FORMAT_FULL_NO_TABS)
         strcpy( buff + strlen( buff) - 2, "mresid");
      if( base_format == RESIDUAL_FORMAT_SHORT)
         memcpy( buff + 13, "Resi MRes", 9);
      }

   if( residual_format & RESIDUAL_FORMAT_HMS)
      text_search_and_replace( buff, ".DDDDD", " HH:MM:SS");

   if( line_no >= 0)
      {
      if( base_format == RESIDUAL_FORMAT_SHORT)
         {
         const int width = getmaxx( stdscr), n_cols = width / 25;
         const int cols = width / n_cols;
         int i;

         buff[cols] = '\0';
         for( i = 0; i < n_cols; i++)
            put_colored_text( buff, line_no, i * cols, -1,
                              COLOR_RESIDUAL_LEGEND);
         }
      else
         put_colored_text( buff, line_no, 0, -1,
                                 COLOR_RESIDUAL_LEGEND);
      }
   mpc_column = strstr( buff, "Obs") - buff;
   resid_column = strstr( buff, "res") - buff - 1;
}

static void show_a_file( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");
   char buff[160], err_text[100];
   int line_no = 0, keep_going = 1;
   int n_lines = 0, msg_num = 0;
   int n_lines_alloced = 0;
   int *index = NULL;

   if( !ifile)
      return;
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( n_lines == n_lines_alloced)
         {
         n_lines_alloced = 2 * n_lines_alloced + 50;
         index = (int *)realloc( index, (n_lines_alloced + 1) * sizeof( int *));
         }
      n_lines++;
      index[n_lines] = ftell( ifile);
      }
   index[0] = 0;
   *err_text = '\0';
   while( keep_going)
      {
      int i, c;
      const char *msgs[] = { "1Cursor keys to move",
                             "2Already at end of file!",
                             "2Already at top of file!" };
      extern const char *ephemeris_filename;
      const int is_ephem = !strcmp( filename, ephemeris_filename);
      const int top_possible_line = (is_ephem ? 3 : 0);
      const int n_lines_to_show = getmaxy( stdscr) - 1;
      int top_line;

      clear( );
      if( line_no < top_possible_line)
         line_no = top_possible_line;
      if( line_no > n_lines - 1)
         line_no = n_lines - 1;
      top_line = line_no - n_lines_to_show / 2;
      if( top_line >= n_lines - n_lines_to_show)
         top_line = n_lines - n_lines_to_show;
      if( top_line < 0)
         top_line = 0;
      if( is_ephem)
         {
         fseek( ifile, 0L, SEEK_SET);
         for( i = 0; i < 3; i++)
            {
            fgets_trimmed( buff, sizeof( buff), ifile);
            put_colored_text( buff, i, 0, -1, COLOR_BACKGROUND);
            }
         }
      fseek( ifile, index[top_line], SEEK_SET);
      for( i = 0; i < n_lines_to_show
                        && fgets_trimmed( buff, sizeof( buff), ifile); i++)
         {
         const int curr_line = top_line + i;

         if( i >= 3 || !is_ephem)
            put_colored_text( buff, i, 0, -1,
               (line_no == curr_line ? COLOR_ORBITAL_ELEMENTS : COLOR_BACKGROUND));
         }
               /* show "scroll bar" to right of text: */
      show_right_hand_scroll_bar( 0, n_lines_to_show, top_line, n_lines);
      sprintf( buff, "   Line %d of %d", line_no, n_lines);
      put_colored_text( buff, i, getmaxx( stdscr) - strlen( buff) - 1,
                                      strlen( buff), COLOR_FINAL_LINE);
      put_colored_text( "Quit", i, 25, 4, COLOR_FINAL_LINE);
      if( line_no < n_lines - 1)
         {
         put_colored_text( "pgDown", i, 30, 6, COLOR_FINAL_LINE);
         put_colored_text( "End", i, 37, 3, COLOR_FINAL_LINE);
         }
      if( line_no)
         {
         put_colored_text( "pgUp", i, 41, 4, COLOR_FINAL_LINE);
         put_colored_text( "Top", i, 46, 3, COLOR_FINAL_LINE);
         }

      strcpy( buff, msgs[msg_num] + 1);
      if( *err_text)
         strcpy( buff, err_text);
      put_colored_text( buff, i, 0, strlen( buff), msgs[msg_num][0] - '0');
      *err_text = '\0';
      msg_num = 0;
      refresh( );
      flushinp( );
      c = extended_getch( );
      if( c == KEY_MOUSE)
         {
         int x, y, z;
         unsigned long button;

         get_mouse_data( &x, &y, &z, &button);
#ifdef BUTTON5_PRESSED
         if( button & BUTTON4_PRESSED)   /* actually 'wheel up' */
            c = KEY_UP;
         else if( button & BUTTON5_PRESSED)   /* actually 'wheel down' */
            c = KEY_DOWN;
         else
#endif
          if( y == i)
            {
            if( x >= 25 && x <= 28)       /* "Quit" */
               c = 'q';
            else if( x >= 30 && x <= 35)  /* "pgDown" */
               c = 'd';
            else if( x >= 37 && x <= 39)  /* "End" */
               c = 'e';
            else if( x >= 41 && x <= 44)  /* "pgUp" */
               c = 'u';
            else if( x >= 46 && x <= 48)  /* "Top" */
               c = 't';
            }
         else if( x == getmaxx( stdscr) - 1)  /* clicked scroll bar */
            line_no = y * n_lines / n_lines_to_show;
         else
            line_no = y + top_line;
         }
      switch( c)
         {
         case KEY_C1:
         case KEY_END:
         case 'e': case 'E':
            line_no = n_lines - 1;
            break;
         case KEY_A1:
         case KEY_HOME:
         case 't': case 'T':
            line_no = top_possible_line;
            break;
         case KEY_UP:
#ifdef KEY_A2
         case KEY_A2:
#endif
            if( line_no > top_possible_line)
               line_no--;
            else
               msg_num = 2;
            break;
         case KEY_DOWN:
         case ' ':
#ifdef KEY_C2
         case KEY_C2:
#endif
            if( line_no >= n_lines - 1)
               msg_num = 1;
            else
               line_no++;
            break;
         case KEY_C3:         /* "PgDn" = lower right key in keypad */
         case KEY_NPAGE:
         case 'd': case 'D':
         case 13:
            if( line_no >= n_lines - 1)
               msg_num = 1;
            else
               line_no += n_lines_to_show - 1;
            break;
         case KEY_A3:         /* "PgUp" = upper right key in keypad */
         case KEY_PPAGE:
         case 'u': case 'U':
            if( line_no > top_possible_line)
               line_no -= n_lines_to_show - 1;
            else
               msg_num = 2;
            break;
#ifdef KEY_RESIZE
         case KEY_RESIZE:
            resize_term( 0, 0);
            break;
#endif
         case 'q': case 'Q': case 27:
#ifdef KEY_EXIT
         case KEY_EXIT:
#endif
            keep_going = 0;
            break;
         default:
            break;
         }
      }
   fclose( ifile);
   free( index);
   refresh( );
}

static int get_epoch_range_of_included_obs( const OBSERVE FAR *obs,
                  const int n_obs, double *start_jd, double *end_jd)
{
   int idx1, idx2, rval;

   rval = get_idx1_and_idx2( n_obs, obs, &idx1, &idx2);
   if( start_jd)
      *start_jd = obs[idx1].jd;
   if( end_jd)
      *end_jd   = obs[idx2].jd;
   return( rval);
}

#ifdef USE_MYCURSES
static void initialize_global_bmouse( void)
{
   memset( &global_bmouse, 0, sizeof( BMOUSE));
   global_bmouse.xmax = getmaxx( stdscr) - 1;
   global_bmouse.ymax = getmaxy( stdscr) - 1;
   global_bmouse.x = getmaxx( stdscr) / 2;
   global_bmouse.y = getmaxy( stdscr) / 2;
   global_bmouse.sensitivity = 0;
   init_mouse( &global_bmouse);
}
#endif

static void get_mouse_data( int *mouse_x, int *mouse_y,
                            int *mouse_z, unsigned long *button)
{
#ifdef USE_MYCURSES
            *mouse_x = global_bmouse.x / 8;
            *mouse_y = global_bmouse.y / 8;
            *mouse_z = 0;
            *button = global_bmouse.released;
#else       /* non-Mycurses case: */
            MEVENT mouse_event;

#ifdef __PDCURSES__
            nc_getmouse( &mouse_event);
#else
            getmouse( &mouse_event);
#endif
            *mouse_x = mouse_event.x;
            *mouse_y = mouse_event.y;
            *mouse_z = mouse_event.z;
            *button  = mouse_event.bstate;
#endif       /* end non-Mycurses case: */
}

static void put_colored_text( const char *text, const int line_no,
               const int column, const int n_bytes, const int color)
{
   attrset( COLOR_PAIR( color & 255));
   if( color & 256)
      attron( A_BLINK);
   if( color & 512)
      attron( A_BOLD);
#ifdef A_LEFTLINE
   if( color & 1024)
      attron( A_LEFTLINE);
   if( color & 2048)
      attron( A_RIGHTLINE);
   if( color & 4096)
      attron( A_ITALIC);
   if( color & 8192)
      attron( A_UNDERLINE);
#endif
#ifdef A_OVERLINE
   if( color & 16384)
      attron( A_OVERLINE);
#endif
   if( n_bytes > 0)
      {
      const int len = getmaxx( stdscr) - column;

      mvaddnstr( line_no, column, text, (n_bytes < len) ? n_bytes : len);
      }
   else              /* clear to end of line */
      {
      int len = strlen( text), remains = getmaxx( stdscr) - len - column;

      if( remains < 0)
         {
         len += remains;
         remains = 0;
         }
      mvaddnstr( line_no, column, text, len);
      if( remains > 0 && remains < 290)
         {
         char buff[300];

         memset( buff, ' ', remains);
         buff[remains] = '\0';
         mvaddnstr( line_no, column + len, buff, remains);
         }
      }
   if( color & 256)
      attroff( A_BLINK);
   if( color & 512)
      attroff( A_BOLD);
#ifdef A_LEFTLINE
   if( color & 1024)
      attroff( A_LEFTLINE);
   if( color & 2048)
      attroff( A_RIGHTLINE);
   if( color & 4096)
      attroff( A_ITALIC);
   if( color & 8192)
      attroff( A_UNDERLINE);
#endif
#ifdef A_OVERLINE
   if( color & 16384)
      attroff( A_OVERLINE);
#endif
}


static void cycle_residual_format( int *residual_format, char *message_to_user)
{
   switch( *residual_format & 3)
      {
      case RESIDUAL_FORMAT_FULL_NO_TABS:      /* 0 */
         (*residual_format) += 2;             /* cycle to short format */
         break;
      case RESIDUAL_FORMAT_SHORT:             /* 2 */
         (*residual_format)++;                /* cycles to MPC 80-col */
         break;
      case RESIDUAL_FORMAT_80_COL:            /* 3 */
         (*residual_format) -= 3;             /* cycles back to full-no-tab */
         break;
      }
   switch( *residual_format & 3)
      {
      case RESIDUAL_FORMAT_FULL_NO_TABS:        /* 0 */
         strcpy( message_to_user, "Standard observation data shown");
         *residual_format |= RESIDUAL_FORMAT_FOUR_DIGIT_YEARS;
         break;
      case RESIDUAL_FORMAT_SHORT:               /* 2 */
         strcpy( message_to_user, "Short MPC residual format selected");
         *residual_format &= ~RESIDUAL_FORMAT_FOUR_DIGIT_YEARS;
         break;
      case RESIDUAL_FORMAT_80_COL:              /* 3 */
         strcpy( message_to_user, "Display of original MPC reports selected");
         break;
      }

}

static OBSERVE FAR *load_object( FILE *ifile, OBJECT_INFO *id,
                       double *curr_epoch, double *epoch_shown, double *orbit)
{
   extern int n_obs_actually_loaded;
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

OBSERVE *add_observations( FILE *ifile, OBSERVE *obs,
                  const OBJECT_INFO *ids, int *n_obs)
{
   OBSERVE *obs2 = load_observations( ifile, ids->packed_desig, ids->n_obs);
   OBSERVE *temp_obs;
   extern int n_obs_actually_loaded;

   if( debug_level)
      printf( "Got %d new obs\n", n_obs_actually_loaded);
   fclose( ifile);
   temp_obs = (OBSERVE *)calloc( *n_obs + n_obs_actually_loaded,
                           sizeof( OBSERVE));
   memcpy( temp_obs, obs, *n_obs * sizeof( OBSERVE));
   memcpy( temp_obs + *n_obs, obs2, n_obs_actually_loaded * sizeof( OBSERVE));
   *n_obs += n_obs_actually_loaded;
   free( obs);
   free( obs2);
   obs = temp_obs;
   *n_obs = sort_obs_by_date_and_remove_duplicates( obs, *n_obs);
   return( obs);
}

int find_first_and_last_obs_idx( const OBSERVE *obs, const int n_obs,
         int *last)
{
   int i;

   if( last)
      {
      for( i = n_obs - 1; i > 0 && !obs[i].is_included; i--)
         ;
      *last = i;
      }
   for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
      ;
   return( i);
}

/* When doing (for example) a full six-parameter fit to an orbit,  it can
be helpful to use an epoch that is at the mid-point of the arc of
observations that is being fitted.  This improves stability.  You can't
do it if you're going to use the covariance matrix,  since you'd then
get uncertainties as of the mid-epoch;  but there are times (specifically,
when doing Monte Carlo orbits) when the covariance matrix isn't especially
meaningful.  */

static double mid_epoch_of_arc( const OBSERVE *obs, const int n_obs)

{
   int first, last;

   first = find_first_and_last_obs_idx( obs, n_obs, &last);
   return( (obs[first].jd + obs[last].jd) / 2.);
}


extern const char *elements_filename;

int sanity_test_observations( const char *filename);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);

/* main( ) begins by using the select_object_in_file( ) function (see above)
   to determine which object is going to be analyzed in the current 'run'.
   It loads up the observations for the object to be analyzed,
   and computes an initial guess for the orbit.

   It then goes into a loop,  setting up and displaying a screenful of
   data,  then asking the user to hit a key.  Based on which key is hit,
   an action such as taking an Herget step or a "full step",  or resetting
   the epoch,  or quitting,  is taken. */

#ifdef PURE_WINDOWS_VERSION
int dummy_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   char obj_name[80], tbuff[300], orbit_constraints[90];
   int c = 1, element_precision, get_new_object = 1, add_off_on = -1;
   int top_line_basic_info_perturbers;
   int top_line_orbital_elements, top_line_residual_legend;
   int top_line_residuals, is_monte_orbit = 0, list_codes = 1;
   OBSERVE FAR *obs;
   int n_obs, i, curr_obs, observation_display = 15, quit = 0;
   double epoch_shown, curr_epoch, orbit[6];
   double r1 = 1., r2 = 1.;
   char message_to_user[80];
   int update_element_display = 1, gauss_soln = 0;
   int residual_format = RESIDUAL_FORMAT_80_COL, bad_elements = 0;
   int element_format = 0, debug_mouse_messages = 0;
#ifdef USE_MYCURSES
   int curr_text_mode = 0;
#endif
   int auto_repeat_full_improvement = 0, n_ids, planet_orbiting = 0;
   OBJECT_INFO *ids;
   double noise_in_arcseconds = 1.;
   double monte_data[MONTE_DATA_SIZE];
   extern int monte_carlo_object_count;
   extern char default_comet_magnitude_type;
   extern double max_monte_rms;
#ifdef __PDCURSES__
   int original_xmax, original_ymax;
#else
   const int blinking_freq = CLOCKS_PER_SEC / 2;
#endif
   double max_residual_for_filtering = 2.5;
   extern int precise_residuals;    /* ephem0.cpp */

   for( i = 1; i < argc; i++)       /* check to see if we're debugging: */
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'd':
               debug_level = atoi( argv[i] + 2);
               if( !debug_level)
                  debug_level = 1;
               debug_printf( "findorb: debug_level = %d; %s %s\n",
                           debug_level, __DATE__, __TIME__);
               break;
            case 'c':
               {
               extern int combine_all_observations;

               combine_all_observations = 1;
               }
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
               break;
            case 'n':
               max_mpc_color_codes = atoi( argv[i] + 2);
               break;
            case 's':
               sanity_test_observations( argv[1]);
               printf( "Sanity check complete\n");
               exit( 0);
               break;
            case 'S':
               {
               extern int sanity_check_observations;

               sanity_check_observations = 0;
               }
               break;
            default:
               printf( "Unknown command-line option '%s'\n", argv[i]);
               return( -1);
            }

   get_defaults( &ephemeris_output_options);                       /* elem_out.c */
   sscanf( get_environment_ptr( "CONSOLE_OPTS"), "%s %d %d %d",
               mpc_code, &observation_display, &residual_format, &list_codes);
   sscanf( get_environment_ptr( "SETTINGS"), "%c,%d,%d,%d,%lf,%lf,%d",
               &default_comet_magnitude_type,
               &element_format, &element_precision,
               &ephemeris_output_options,
               &max_residual_for_filtering, &noise_in_arcseconds, &precise_residuals);

   strcpy( ephemeris_start, get_environment_ptr( "EPHEM_START"));
   sscanf( get_environment_ptr( "EPHEM_STEPS"), "%d %s",
               &n_ephemeris_steps, ephemeris_step_size);
   if( debug_level)
      debug_printf( "Options read\n");
   load_up_weight_records( "weight.txt");
   if( debug_level)
      debug_printf( "Default weighting table read\n");

   if( argc < 2)
      {
      printf( "'findorb' needs the name of an input file of MPC-formatted\n");
      printf( "astrometry as a command-line argument.\n");
      exit( 0);
      }

   *message_to_user = '\0';

#ifdef _WIN32
   if( !strcmp( argv[1], "c") || !strcmp( argv[1], "c+"))
      {
      const char *temp_clipboard_filename = "obs_temp.txt";

      clipboard_to_file( temp_clipboard_filename, argv[1][1] == '+');
      argv[1] = temp_clipboard_filename;
      }
#endif

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

   if( debug_level > 2)
      debug_printf( "%d objects in file\n", n_ids);
   if( debug_level > 3)
      for( i = 0; i < n_ids; i++)
         {
         debug_printf( "   Object %d: '%s', '%s'\n", i,
                        ids[i].packed_desig, ids[i].obj_name);
         object_comment_text( tbuff, ids + i);
         debug_printf( "   %s\n", tbuff);
         }
   if( n_ids > 0)
      set_solutions_found( ids, n_ids);
   if( debug_level > 2)
      debug_printf( "solutions set\n");

   if( debug_level > 2)
      debug_printf( "Initializing curses...");
#ifdef XCURSES
   Xinitscr( argc, (char **)argv);
#else
   initscr( );
#endif
   if( debug_level > 2)
      debug_printf( "Curses initialised, ");
   cbreak( );
   noecho( );
   clear( );
   refresh( );
   if( debug_level > 2)
      debug_printf( "(2), ");
   start_color( );
   init_pair( COLOR_BACKGROUND, COLOR_WHITE, COLOR_BLACK);
   init_pair( COLOR_ORBITAL_ELEMENTS, COLOR_BLACK, COLOR_YELLOW);
   init_pair( COLOR_FINAL_LINE, COLOR_WHITE, COLOR_BLUE);
   init_pair( COLOR_SELECTED_OBS, COLOR_WHITE, COLOR_MAGENTA);
   init_pair( COLOR_HIGHLIT_BUTTON, COLOR_BLACK, COLOR_GREEN);
   init_pair( COLOR_EXCLUDED_OBS, COLOR_RED, COLOR_GREEN);
   init_pair( COLOR_OBS_INFO, COLOR_WHITE, COLOR_RED);
   init_pair( COLOR_MESSAGE_TO_USER, COLOR_BLACK, COLOR_WHITE);
   init_pair( COLOR_RESIDUAL_LEGEND, COLOR_BLACK, COLOR_CYAN);
   init_pair( COLOR_MENU, COLOR_GREEN, COLOR_MAGENTA);

                  /* MPC color-coded station colors: */
   init_pair( 16, COLOR_YELLOW, COLOR_BLACK);
   init_pair( 17, COLOR_BLUE, COLOR_BLACK);
   init_pair( 18, COLOR_MAGENTA, COLOR_BLACK);
   init_pair( 19, COLOR_RED, COLOR_BLACK);
   init_pair( 20, COLOR_GREEN, COLOR_BLACK);

#ifndef USE_MYCURSES
   init_color( COLOR_GRAY, 500, 500, 500);
#endif
   init_pair( COLOR_SCROLL_BAR, COLOR_GREEN, COLOR_GRAY);

   if( debug_level > 2)
      debug_printf( "(3)\n");
   keypad( stdscr, 1);
// curs_set( 0);
#ifdef __PDCURSES__
   original_xmax = getmaxx( stdscr);
   original_ymax = getmaxy( stdscr);
   PDC_set_blink( TRUE);
   PDC_set_title( "Find_Orb -- Orbit Determination");
   if( strstr( longname( ), "SDL"))
      resize_term( 40, 110);
   if( !strcmp( longname( ), "Win32"))
      resize_term( 52, 126);
#endif
#if !defined (_WIN32) && !defined( __WATCOMC__)
   mousemask( ALL_MOUSE_EVENTS, NULL);
#endif
#ifdef _WIN32
   mouse_set( ALL_MOUSE_EVENTS);
#endif
#ifdef USE_MYCURSES
   initialize_global_bmouse( );
#endif
   while( !quit)
      {
      int obs_per_line;
      int base_format = (residual_format & 3);
      int line_no = 0, total_obs_lines;
      extern double solar_pressure[];
      extern int n_extra_params;

      if( debug_level > 3)
         debug_printf( "get_new_object = %d\n", get_new_object);
      if( get_new_object)
         {
         int id_number = 0;

         if( n_ids > 1)
            id_number = select_object_in_file( ids, n_ids);
         if( debug_level > 3 && id_number >= 0)
            debug_printf( "id_number = %d; '%s'\n", id_number,
                                    ids[id_number].obj_name);
         get_new_object = 0;
         *orbit_constraints = '\0';
         if( id_number < 0)
            {
            attrset( COLOR_PAIR( COLOR_BACKGROUND));
#ifdef __PDCURSES__
            resize_term( original_ymax, original_xmax);
#endif
            endwin( );
#ifdef XCURSES
            XCursesExit();
#endif
            printf( "User exited program\n");
            return( -1);
            }
         else
            {
            FILE *ifile;
            long file_offset;

            strcpy( obj_name, ids[id_number].obj_name);
            sprintf( tbuff, "Loading '%s'...", obj_name);
            put_colored_text( tbuff, getmaxy( stdscr) - 3,
                                 0, -1, COLOR_FINAL_LINE);
            if( debug_level)
               debug_printf( "%s: ", tbuff);
            refresh( );
            monte_carlo_object_count = 0;

            ifile = fopen( argv[1], "rb");
                /* Start a bit ahead of the actual data,  just in case */
                /* there's a #Weight: or similar command in there: */
            file_offset = ids[id_number].file_offset - 40L;
            if( file_offset < 0L)
               file_offset = 0L;
            fseek( ifile, file_offset, SEEK_SET);

            obs = load_object( ifile, ids + id_number, &curr_epoch,
                                                  &epoch_shown, orbit);
            fclose( ifile);
            n_obs = ids[id_number].n_obs;
            if( debug_level)
               debug_printf( "got obs; ");
            if( mpc_color_codes)
               free( mpc_color_codes);
            if( max_mpc_color_codes)
               mpc_color_codes = find_mpc_color_codes( n_obs, obs,
                          max_mpc_color_codes);
            if( debug_level)
               debug_printf( "got color codes; ");
            get_r1_and_r2( n_obs, obs, &r1, &r2);    /* orb_func.cpp */
            if( debug_level)
               debug_printf( "R1 = %lf; R2 = %lf\n", r1, r2);
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
               ;
            curr_obs = i;
            update_element_display = 1;
            clear( );
            }
         }

      if( curr_obs > n_obs - 1)
         curr_obs = n_obs - 1;
      if( curr_obs < 0)
         curr_obs = 0;
      if( debug_level > 2)
         debug_printf( "update_element_display = %d\n", update_element_display);
      if( update_element_display)
         bad_elements = write_out_elements_to_file( orbit, curr_epoch, epoch_shown,
             obs, n_obs, orbit_constraints, element_precision,
             is_monte_orbit, element_format);
      is_monte_orbit = 0;
      if( debug_level > 2)
         debug_printf( "elements written\n");
      update_element_display = 0;
      top_line_basic_info_perturbers = line_no;
      if( observation_display & 1)
         {
         if( c != KEY_MOUSE_MOVE)
            {
            show_basic_info( obs, n_obs);
            show_perturbers( );
            if( debug_level)
               refresh( );
            }
         line_no = 2;
         }
      if( observation_display & 2)
         for( i = 0; i < 4; i++)
            {
            char buff[190];

            generate_observation_text( obs, n_obs, curr_obs, i, buff);
            if( i == 3 && strlen( buff) < 72)
               {
               const time_t t0 = time( NULL);

               memset( buff + strlen( buff), ' ', 72 - strlen( buff));
               strcpy( buff + 72, ctime( &t0) + 11);
               buff[81] = '\0';
               }
            if( *buff)
               put_colored_text( buff, line_no++, 0, -1, COLOR_OBS_INFO);
            if( debug_level)
               refresh( );
            }
      top_line_orbital_elements = line_no;
      if( observation_display & 4)
         {
         FILE *ifile;

         ifile = fopen( elements_filename, "rb");
         if( ifile)
            {
            int iline = 0;
#ifndef __PDCURSES__
            const int is_blinking = (clock( ) / blinking_freq) % 2;
            int elem_color = ((bad_elements && is_blinking) ?
                                256 + COLOR_OBS_INFO : COLOR_ORBITAL_ELEMENTS);
#else
            int elem_color = (bad_elements ?
                                256 + COLOR_OBS_INFO : COLOR_ORBITAL_ELEMENTS);
#endif
            char covar_text[200];
            FILE *covar_file = NULL;
            extern double uncertainty_parameter;

            *covar_text = '\0';
            if( uncertainty_parameter < 98.)
               {
               covar_file = fopen( "covar.txt", "rb");
               while( memcmp( covar_text, "Sigmas:", 7) &&
                        fgets_trimmed( covar_text, sizeof( covar_text), covar_file))
                  ;
               }

            while( fgets_trimmed( tbuff, sizeof( tbuff), ifile))
               if( iline < 20 && *tbuff != '#')
                  {
                  char *tptr = strstr( tbuff, "Earth MOID:");

                  if( !memcmp( tbuff, "IMPACT", 6))
#ifndef __PDCURSES__
                     elem_color = (is_blinking ? COLOR_ATTENTION : COLOR_ORBITAL_ELEMENTS);
#else
                     elem_color = COLOR_ATTENTION + 256;
#endif
                  put_colored_text( tbuff, line_no + iline, 0, -1, elem_color);
                  if( covar_file)
                     {
                     int replace_pq;   /* replace PQ w/covar data */

                     if( element_format & ELEM_OUT_ALTERNATIVE_FORMAT)
                        replace_pq = (iline >= 4 && iline <= 6);
                     else
                        replace_pq = (iline >= 5 && iline <= 7);
                     if( replace_pq && fgets_trimmed(
                             covar_text, sizeof( covar_text), covar_file))
                        put_colored_text( covar_text, line_no + iline, 40, -1, elem_color);
                     if( fgets_trimmed(
                             covar_text, sizeof( covar_text), covar_file))
                        put_colored_text( covar_text, line_no + iline, 75, -1, elem_color);
                     }

                  if( tptr)         /* low Earth MOID:  show in flashing text to draw attn */
                     if( atof( tptr + 11) < .01)
                        put_colored_text( tptr, line_no + iline, tptr - tbuff,
                                       20, COLOR_ATTENTION + 256);
                  iline++;
                  }
               else if( !memcmp( tbuff, "# Tisserand", 11))
                  {
                  tbuff[30] = tolower( tbuff[24]);
                  tbuff[29] = 'T';
                  put_colored_text( tbuff + 29, line_no + (tbuff[30] == 'n'),
                                          60, 15, elem_color);
                  }
            line_no += iline;
            fclose( ifile);
            if( covar_file)
               fclose( covar_file);
            }
         if( debug_level)
            refresh( );
         }
      top_line_residual_legend = line_no;
      if( observation_display & 8)
         {
         if( c != KEY_MOUSE_MOVE && c != KEY_TIMER)
            show_residual_legend( line_no, residual_format);
         line_no++;
         }
      else
         if( c != KEY_MOUSE_MOVE && c != KEY_TIMER)
            show_residual_legend( -1, residual_format);
      if( debug_level)
         refresh( );
      if( debug_level > 2)
         debug_printf( "resid legend shown\n");

      top_line_residuals = line_no;
      if( c != KEY_MOUSE_MOVE && c != KEY_TIMER)
         show_residuals( obs, n_obs, residual_format, curr_obs, line_no,
                                       list_codes);
      if( debug_level > 2)
         debug_printf( "resids shown\n");
      if( debug_level)
         refresh( );
      if( c != KEY_MOUSE_MOVE && c != KEY_TIMER)
         show_final_line( obs, n_obs, curr_obs, COLOR_FINAL_LINE);
      if( debug_level)
         refresh( );
      if( *message_to_user)
         {
         int xloc;

         if( add_off_on >= 0)
            strcat( message_to_user, (add_off_on ? " on" : " off"));
         xloc = getmaxx( stdscr) - strlen( message_to_user) - 1;
         put_colored_text( message_to_user, getmaxy( stdscr) - 1,
                           (xloc < 0 ? 0 : xloc), -1, COLOR_MESSAGE_TO_USER);
         }
      add_off_on = -1;
      refresh( );
      move( getmaxy( stdscr) - 1, 0);
      if( c == AUTO_REPEATING)
         if( curses_kbhit( ) != ERR)
            {
            extended_getch( );
            c = 0;
            }
      if( c != AUTO_REPEATING)
         {
         int x1, y1, z1;
         unsigned long button;
#ifndef __PDCURSES__
         const int blink_state = (int)( clock( ) / blinking_freq);
#endif

//       flushinp( );
         c = 0;
         get_mouse_data( &x1, &y1, &z1, &button);
         while( !c && curses_kbhit( ) == ERR)
            {
            int x, y, z;

            get_mouse_data( &x, &y, &z, &button);
            if( x != x1 || y != y1)
               c = KEY_MOUSE_MOVE;
#ifndef __PDCURSES__
            else if( blink_state != (int)( clock( ) / blinking_freq))
               c = KEY_TIMER;
#else
            napms( 200);
#endif
            }
         if( !c)
            c = extended_getch( );
         auto_repeat_full_improvement = 0;
         }
      if( c != KEY_MOUSE_MOVE && c != KEY_TIMER)
         *message_to_user = '\0';

      if( base_format == RESIDUAL_FORMAT_SHORT)  /* multi-obs per line */
         {
         obs_per_line = getmaxx( stdscr) / 25;
         total_obs_lines = (n_obs - 1) / obs_per_line + 1;
         }
      else
         {
         obs_per_line = 1;
         total_obs_lines = n_obs;
         }

      if( c == KEY_MOUSE)
         {
         int x, y, z;
         const int station_start_line = getmaxy( stdscr) - n_stations_shown;
         unsigned long button;

         get_mouse_data( &x, &y, &z, &button);
#ifndef USE_MYCURSES
         if( debug_mouse_messages)
            sprintf( message_to_user, "x=%d y=%d z=%d button=%lx",
                              x, y, z, button);
#endif
#ifdef BUTTON5_PRESSED
         if( button & BUTTON4_PRESSED)   /* actually 'wheel up' */
            c = KEY_UP;
         else if( button & BUTTON5_PRESSED)   /* actually 'wheel down' */
            c = KEY_DOWN;
         else
#endif
            if( y >= station_start_line)
               {
               const char *search_code =
                       mpc_color_codes[y - station_start_line].code;
#ifdef USE_MYCURSES
               const int dir = 1;
#else
               const int button1_events = BUTTON1_CLICKED
                         | BUTTON1_DOUBLE_CLICKED
                         | BUTTON1_TRIPLE_CLICKED;
               const int dir =
                   (( button & button1_events) ? 1 : n_obs - 1);
#endif

               curr_obs = (curr_obs + dir) % n_obs;
               while( FSTRCMP( obs[curr_obs].mpc_code, search_code))
                  curr_obs = (curr_obs + dir) % n_obs;
               }
         else if( y >= top_line_residuals)
            {
            const int max_x = getmaxx( stdscr);

            if( x == max_x - 1)
               curr_obs = (y - top_line_residuals) * (n_obs - 1)
                           / (station_start_line - top_line_residuals - 1);
            else
               {
               curr_obs = first_residual_shown +
                      (y - top_line_residuals);
               if( base_format == RESIDUAL_FORMAT_SHORT)  /* multi-obs per line */
                  curr_obs += total_obs_lines * (x * obs_per_line / max_x);
#ifndef USE_MYCURSES
               if( button & BUTTON1_DOUBLE_CLICKED)
                  obs[curr_obs].is_included ^= 1;
#endif
               }
            }
         else if( y >= top_line_residual_legend &&
                                  y < top_line_residuals)
            c = 'k';       /* cycle the residual format */
         else if( (observation_display & 1) &&
                          y == top_line_basic_info_perturbers + 1)
            {                      /* clicked on a perturber 'radio button' */
            if( x / 7 == 9)
               c = '0';
            else if( x / 7 == 10)
               c = 'a';
            else
               c = '1' + (x / 7);
            }
         else if( (observation_display & 1) &&
                          y == top_line_basic_info_perturbers)
            {
            if( x < 24)             /* clicked on R1/R2 */
               c = 'r';
            else if( x - 24 < (int)strlen( command_remap))
               c = command_remap[x - 24];
            }                    /* clicked on top-line command list */
         else if( (observation_display & 4) &&
                          y == top_line_orbital_elements + 3 && x < 40)
            c = 'e';             /* Reset epoch of elements */
         }

      if( c >= '1' && c <= '9')
         perturbers ^= (1 << (c - '0'));
#ifdef ALT_0
      else if( c >= ALT_0 && c <= ALT_9)
         curr_obs = (n_obs - 1) * (c - ALT_0) / 10;
#endif
      else switch( c)
         {
         case '0':
            perturbers ^= 1024;
            break;
         case KEY_C1:
         case KEY_END:
            curr_obs = n_obs - 1;
            break;
         case KEY_A1:
         case KEY_HOME:
            curr_obs = 0;
            break;
         case KEY_UP:
#ifdef KEY_A2
         case KEY_A2:
#endif
            curr_obs--;
            break;
         case KEY_DOWN:
#ifdef KEY_C2
         case KEY_C2:
#endif
            curr_obs++;
            break;
         case KEY_LEFT:
#ifdef KEY_B1
         case KEY_B1:
#endif
            if( base_format == RESIDUAL_FORMAT_SHORT)  /* multi-obs per line */
               curr_obs -= total_obs_lines;
            else
               curr_obs--;
            break;
         case KEY_RIGHT:
#ifdef KEY_B3
         case KEY_B3:
#endif
            if( base_format == RESIDUAL_FORMAT_SHORT)  /* multi-obs per line */
               curr_obs += total_obs_lines;
            else
               curr_obs++;
            break;
         case KEY_C3:         /* "PgDn" = lower right key in keypad */
         case KEY_NPAGE:
            curr_obs += getmaxy( stdscr) - line_no - n_stations_shown;
            break;
         case KEY_A3:         /* "PgUp" = upper right key in keypad */
         case KEY_PPAGE:
            curr_obs -= getmaxy( stdscr) - line_no - n_stations_shown;
            break;
         case KEY_F(1):      /* turn on/off all obs prior to curr one */
         case KEY_F(13):     /* Shift-F1:  plain ol' F1 has trouble in XWin */
            obs[curr_obs].is_included ^= 1;
            for( i = 0; i < curr_obs; i++)
               obs[i].is_included = obs[curr_obs].is_included;
            strcpy( message_to_user, "All preceding observations toggled");
            add_off_on = obs[curr_obs].is_included;
            break;
         case KEY_F(16):     /* Shift-F4:  select another object to add in */
            if( (i = select_object_in_file( ids, n_ids)) >= 0)
               {
               FILE *ifile = fopen( argv[1], "rb");

               obs = add_observations( ifile, obs, ids + i, &n_obs);
               fclose( ifile);
               if( debug_level)
                  printf( "Now have %d obs\n", n_obs);
               if( mpc_color_codes)
                  free( mpc_color_codes);
               if( max_mpc_color_codes)
                  mpc_color_codes = find_mpc_color_codes( n_obs, obs,
                             max_mpc_color_codes);
               if( debug_level)
                  debug_printf( "got color codes; ");
               for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                  ;
               curr_obs = i;
               set_locs( orbit, curr_epoch, obs, n_obs);
               update_element_display = 1;
               clear( );
               }
            break;
         case KEY_F(2):          /* turn on/off all obs after curr one */
            obs[curr_obs].is_included ^= 1;
            for( i = curr_obs; i < n_obs; i++)
               obs[i].is_included = obs[curr_obs].is_included;
            strcpy( message_to_user, "All subsequent observations toggled");
            add_off_on = obs[curr_obs].is_included;
            break;
         case KEY_F(3):          /* turn on/off all obs w/same observatory ID */
            obs[curr_obs].is_included ^= 1;
            for( i = 0; i < n_obs; i++)
               if( !FSTRCMP( obs[i].mpc_code, obs[curr_obs].mpc_code))
                  obs[i].is_included = obs[curr_obs].is_included;
            strcpy( message_to_user, "All observations from xxx toggled");
            FMEMCPY( message_to_user + 22, obs[curr_obs].mpc_code, 3);
            add_off_on = obs[curr_obs].is_included;
            break;
         case KEY_F(4):          /* find prev obs from this observatory */
         case KEY_F(5):          /* find next obs from this observatory */
            {
            const int dir = (c == KEY_F(4) ? n_obs - 1 : 1);
            int new_obs = (curr_obs + dir) % n_obs;

            while( new_obs != curr_obs &&
                        FSTRCMP( obs[new_obs].mpc_code, obs[curr_obs].mpc_code))
               new_obs = (new_obs + dir) % n_obs;
            curr_obs = new_obs;
            }
            break;
         case KEY_F(6):          /* find prev excluded obs */
         case KEY_F(7):          /* find next excluded obs */
            {
            const int dir = (c == KEY_F(6) ? n_obs - 1 : 1);
            int new_obs = (curr_obs + dir) % n_obs;

            while( new_obs != curr_obs && obs[new_obs].is_included)
               new_obs = (new_obs + dir) % n_obs;
            curr_obs = new_obs;
            }
            break;
         case '*':         /* toggle use of solar radiation pressure */
            n_extra_params ^= 1;
            strcpy( message_to_user, "Solar radiation pressure is now");
            solar_pressure[0] = solar_pressure[1] = 0.;
            add_off_on = n_extra_params;
            break;
         case '^':
            if( n_extra_params == 2)
               {
               n_extra_params = 3;
               strcpy( message_to_user, "Three-parameter comet non-gravs");
               }
            else if( n_extra_params == 3)
               {
               n_extra_params = 0;
               strcpy( message_to_user, "Comet non-gravs off");
               }
            else
               {
               n_extra_params = 2;
               strcpy( message_to_user, "Two-parameter comet non-gravs");
               }
            solar_pressure[0] = solar_pressure[1] = solar_pressure[2] = 0.;
            break;
#ifdef __WATCOMC__
         case KEY_F(8):     /* show original screens */
            endwin( );
            extended_getch( );
            refresh( );
            break;
#endif
         case 'a': case 'A':
            perturbers ^= (7 << 20);
            strcpy( message_to_user, "Asteroids toggled");
            add_off_on = (perturbers >> 20) & 1;
            break;
         case 'b': case 'B':
            residual_format ^= RESIDUAL_FORMAT_HMS;
            strcpy( message_to_user,
                 (residual_format & RESIDUAL_FORMAT_HMS) ?
                 "Showing observation times as HH:MM:SS" :
                 "Showing observation times as decimal days");
            break;
         case 'c': case 'C':
#ifdef USE_MYCURSES
            if( c == 'c')
               curr_text_mode++;
            else
               curr_text_mode += N_TEXT_MODES - 1;
            curr_text_mode %= N_TEXT_MODES;
            set_text_mode( curr_text_mode);
            initialize_global_bmouse( );
#else
            {
            int new_xsize, new_ysize;

            inquire( "New screen size: ", tbuff, sizeof( tbuff),
                              COLOR_DEFAULT_INQUIRY);
            if( sscanf( tbuff, "%d %d", &new_xsize, &new_ysize) == 2)
               resize_term( new_ysize, new_xsize);
            }
#endif
            sprintf( message_to_user, "%d x %d text mode selected",
                        getmaxx( stdscr), getmaxy( stdscr));
            break;
         case '!':
            perturbers = ((perturbers == 0x3fe) ? 0 : 0x3fe);
            break;
         case KEY_F(9):           /* find start of included arc */
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
               ;
            curr_obs = i;
            break;
         case KEY_F(10):          /* find end of included arc */
            for( i = n_obs - 1; i > 0 && !obs[i].is_included; i--)
               ;
            curr_obs = i;
            break;
         case KEY_F(11):
            auto_repeat_full_improvement ^= 1;
            strcpy( message_to_user, "Automatic full improvement repeat is");
            add_off_on = auto_repeat_full_improvement;
            break;
         case 'd': case 'D':
            if( base_format == RESIDUAL_FORMAT_80_COL)
               {
               create_obs_file( obs, n_obs, 0);
               show_a_file( "observe.txt");
               }
            else
               {
               extern const char *residual_filename;

               write_residuals_to_file( residual_filename, argv[1], n_obs, obs,
                           residual_format);
               show_a_file( residual_filename);
               }
            break;
         case 'e': case'E':
            {
            double new_jd = epoch_shown;
            const double min_jd = 10000.;  /* = year -4685,  i.e., */
                                           /* obviously wrong */

            inquire( "Enter new epoch,  as YYYY MM DD, or JD,  or 'now':",
                             tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            if( extract_date( tbuff, &new_jd) >= 0)
               if( new_jd > min_jd)
                  epoch_shown = floor( new_jd * 100. + .5) / 100.;
            }
            update_element_display = 1;
            break;
         case 'f': case 'F':        /* do a "full improvement" step */
         case '|':                  /* Monte Carlo step */
         case AUTO_REPEATING:
            {
            double *stored_ra_decs = NULL;
            int err = 0;
            extern int using_sr;
            const clock_t t0 = clock( );

            if( c == '|')
               {
               const double rms = compute_rms( obs, n_obs);

               set_statistical_ranging( 0);
               inquire( "Gaussian noise level (arcsec): ",
                             tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
               noise_in_arcseconds = atof( tbuff);
               max_monte_rms =
                           sqrt( noise_in_arcseconds * noise_in_arcseconds
                                        + rms * rms);
               c = AUTO_REPEATING;
               }
            if( c == AUTO_REPEATING)
               stored_ra_decs =
                   add_gaussian_noise_to_obs( n_obs, obs, noise_in_arcseconds);
            push_orbit( curr_epoch, orbit);
            if( !strcmp( orbit_constraints, "e=1"))
               improve_parabolic( obs, n_obs, orbit, curr_epoch);
            else if( c == AUTO_REPEATING && using_sr)
               {
               debug_printf( "Getting %d SR\n", monte_carlo_object_count);
               find_nth_sr_orbit( orbit, obs, n_obs, monte_carlo_object_count);
               adjust_herget_results( obs, n_obs, orbit);
                        /* epoch is that of first included observation: */
               get_epoch_range_of_included_obs( obs, n_obs, &curr_epoch, NULL);
               }
            else
               {
               double saved_orbit[6];
               const double mid_epoch = mid_epoch_of_arc( obs, n_obs);
               const double epoch_to_use = (c != AUTO_REPEATING ? epoch_shown :
                                                      mid_epoch);

               memcpy( saved_orbit, orbit, 6 * sizeof( double));
               integrate_orbit( orbit, curr_epoch, epoch_to_use);
                        /* Only request sigmas for i=1 (last pass... */
                        /* _only_ pass for a real 'full improvement') */
               for( i = (c == AUTO_REPEATING ? 5 : 1); i && !err; i--)
                  {
                  int sigma_type = NO_ORBIT_SIGMAS_REQUESTED;

                  if( i == 1 && c != AUTO_REPEATING)
                     sigma_type = ((element_format & ELEM_OUT_HELIOCENTRIC_ONLY) ?
                           HELIOCENTRIC_SIGMAS_ONLY : ORBIT_SIGMAS_REQUESTED);
                  err = full_improvement( obs, n_obs, orbit, epoch_to_use,
                              orbit_constraints, sigma_type);
                  }
               if( err)    /* Full Monte Carlo isn't working.  Let's try SR, */
                  {        /* & recover from error by using the saved orbit */
                  set_statistical_ranging( 1);
                  memcpy( orbit, saved_orbit, 6 * sizeof( double));
                  }
               else
                  {
                  curr_epoch = mid_epoch;
                  integrate_orbit( orbit, epoch_to_use, curr_epoch);
                  }
               }
            get_r1_and_r2( n_obs, obs, &r1, &r2);
            if( c == AUTO_REPEATING)
               remove_gaussian_noise_from_obs( n_obs, obs, stored_ra_decs);
            update_element_display = (err ? 0 : 1);
            if( c == AUTO_REPEATING && !err)
               {
               extern double planet_mass[];
               ELEMENTS elem;
               double rel_orbit[6], orbit2[6];
               int curr_planet_orbiting;

               is_monte_orbit = 1;
               memcpy( orbit2, orbit, 6 * sizeof( double));
               integrate_orbit( orbit2, curr_epoch, epoch_shown);
               sprintf( message_to_user, "Monte Carlo %d",
                                    monte_carlo_object_count);
               curr_planet_orbiting = find_best_fit_planet( epoch_shown,
                                  orbit2, rel_orbit);
               if( !monte_carlo_object_count)
                  planet_orbiting = curr_planet_orbiting;

               if( planet_orbiting == curr_planet_orbiting)
                  {
                  calc_classical_elements( &elem, rel_orbit, epoch_shown, 1,
                               SOLAR_GM * planet_mass[planet_orbiting]);
                  add_monte_orbit( monte_data, &elem, monte_carlo_object_count);
                  }
               if( monte_carlo_object_count > 3)
                  {
                  double sigmas[MONTE_N_ENTRIES];
                  FILE *monte_file = fopen( "monte.txt", "wb");

                  fprintf( monte_file,
                          "Computed from %d orbits around object %d\n",
                          monte_carlo_object_count, planet_orbiting);
                  compute_monte_sigmas( sigmas, monte_data,
                                             monte_carlo_object_count);
                  dump_monte_data_to_file( monte_file, sigmas,
                        elem.major_axis, elem.ecc, planet_orbiting);
                  fclose( monte_file);
                  }
               }
            else
               {
               strcpy( message_to_user,
                               (err ? "Full step FAILED" : "Full step taken"));
               sprintf( message_to_user + strlen( message_to_user), "(%.3lf s)",
                      (double)( clock( ) - t0) / (double)CLOCKS_PER_SEC);
               }
            }
            break;
         case 'g': case 'G':        /* do a method of Gauss soln */
            {
            double new_epoch;

            perturbers = 0;
            push_orbit( curr_epoch, orbit);
            new_epoch = convenient_gauss( obs, n_obs, orbit, 1., gauss_soln);
            if( !new_epoch && gauss_soln)
               {
               gauss_soln = 0;
               new_epoch = convenient_gauss( obs, n_obs, orbit, 1., gauss_soln);
               }
            gauss_soln++;
            if( new_epoch)
               {
               curr_epoch = new_epoch;
               set_locs( orbit, curr_epoch, obs, n_obs);
               get_r1_and_r2( n_obs, obs, &r1, &r2);
               update_element_display = 1;
               strcpy( message_to_user, "Gauss solution found");
               }
            else
               strcpy( message_to_user, "Gauss method failed!");
            }
            break;
         case '#':
         case '/':
            push_orbit( curr_epoch, orbit);
                     /* epoch is that of first valid observation: */
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                ;
            if( c == '#')
               {
               simplex_method( obs + i, n_obs - i, orbit, r1, r2);
               strcpy( message_to_user, "Simplex method used");
               }
            else
               {
               integrate_orbit( orbit, curr_epoch, obs[i].jd);
               superplex_method( obs + i, n_obs - i, orbit);
               strcpy( message_to_user, "Superplex method used");
               }
//          integrate_orbit( orbit, obs[i].jd, curr_epoch);
            curr_epoch = obs[i].jd;
            get_r1_and_r2( n_obs, obs, &r1, &r2);    /* orb_func.cpp */
            update_element_display = 1;
            break;
         case 'h':               /* do a method of Herget step */
         case 'H':               /* same, but don't linearize */
         case ':':               /* just linearize */
            {
            int err = 0;

            push_orbit( curr_epoch, orbit);
            if( c != ':')
               {
               double d_r1, d_r2;

               err = herget_method( obs, n_obs, r1, r2, orbit, &d_r1, &d_r2,
                                          orbit_constraints);
               if( !err)
                  {
                  r1 += d_r1;
                  r2 += d_r2;
                  herget_method( obs, n_obs, r1, r2, orbit, NULL, NULL, NULL);
                  }
               }
            if( c != 'H')
               if( !err || (c == ':' && n_obs == 2))
                  err = adjust_herget_results( obs, n_obs, orbit);
                     /* epoch is that of first valid observation: */
            if( !err)
               {
               for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                  ;
               curr_epoch = obs[i].jd;
               update_element_display = 1;
               strcpy( message_to_user, (c == ':') ? "Orbit linearized" :
                                                  "Herget step taken");
               }
            else
               strcpy( message_to_user, (c == ':') ? "Linearizing FAILED" :
                                                  "Herget step FAILED");
            }
            break;
         case '<':
            push_orbit( curr_epoch, orbit);
            herget_method( obs, n_obs, r1, r2, orbit, NULL, NULL, NULL);
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
               ;
            curr_epoch = obs[i].jd;
            update_element_display = 1;
            sprintf( message_to_user, "Radii set: %lf %lf", r1, r2);
            break;
         case 'i': case 'I':
            observation_display ^= 2;
            strcpy( message_to_user, "Observation details toggled");
            add_off_on = (observation_display & 2);
            clear( );
            break;
         case 'j': case 'J':
            observation_display ^= 1;
            strcpy( message_to_user,
                     "Perturber/R1 & R2/step size data display toggled");
            add_off_on = (observation_display & 1);
            clear( );
            break;
         case 'k': case 'K':
            cycle_residual_format( &residual_format, message_to_user);
            break;
         case 'l': case 'L':
            if( !*orbit_constraints)
               inquire(
"Enter limits on the orbit (e.g.,  'e=0' or 'q=2.3' or 'q=.7,P=1.4').\n\
Constraints can be placed on e, q, Q, P, a, n, or i.",
                     orbit_constraints, sizeof( orbit_constraints),
                     COLOR_DEFAULT_INQUIRY);
            else
               {
               *orbit_constraints = '\0';
               strcpy( message_to_user, "Orbit is now unconstrained");
               }
            break;
         case 'm': case 'M':
            {
            extern const char *residual_filename;

            create_obs_file( obs, n_obs, 0);
            create_ephemeris( orbit, curr_epoch, obs, n_obs);
            write_residuals_to_file( residual_filename,
                             argv[1], n_obs, obs, RESIDUAL_FORMAT_SHORT);
            make_pseudo_mpec( "mpec.htm", obj_name);      /* ephem0.cpp */
            }
            break;
         case 'n': case 'N':   /* select a new object from the input file */
            get_new_object = 1;
            update_element_display = 1;
            pop_all_orbits( );
            break;
         case 'o': case 'O':
            observation_display ^= 4;
            strcpy( message_to_user, "Display of orbital elements toggled");
            add_off_on = (observation_display & 4);
            clear( );
            break;
         case 'p':
            if( element_precision < 15)
               {
               element_precision++;
               update_element_display = 1;
               }
            break;
         case 'P':
            if( element_precision > 1)
               {
               element_precision--;
               update_element_display = 1;
               }
            break;
         case CTRL( 'F'):
            {
            extern double one_sigma_eigenvect[];
            double sigma;
            static double total_sigma = 0.;

            inquire( "Enter sigma: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY);
            if( *tbuff == '0')
               total_sigma = 0.;
            else
               {
               sigma = atof( tbuff);
               integrate_orbit( orbit, curr_epoch, epoch_shown);
               for( i = 0; i < 6; i++)
                  orbit[i] += sigma * one_sigma_eigenvect[i];
               for( i = 0; i < n_extra_params; i++)
                  solar_pressure[i] += sigma * one_sigma_eigenvect[i + 6];
               integrate_orbit( orbit, epoch_shown, curr_epoch);
               set_locs( orbit, curr_epoch, obs, n_obs);
               total_sigma += sigma;
               update_element_display = 1;
               sprintf( message_to_user,  "Epoch = %lf (%lf); sigma %lf",
                        curr_epoch, epoch_shown, total_sigma);
               }
            }
            break;
         case CTRL( 'R'):
            {
            extern double sr_min_r, sr_max_r;
            extern int using_sr;

            inquire( "Enter SR R1, R2: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY);
            if( sscanf( tbuff, "%lf %lf %lf", &sr_min_r, &sr_max_r, &max_monte_rms) == 3)
               using_sr = 1;
            }
            break;
         case 'r': case 'R':
            inquire( "Enter new R1, R2: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY);
            if( sscanf( tbuff, "%lf%n", &r1, &i) == 1)
               {
               int j;

               while( tbuff[i] == ' ')
                  i++;
               if( tolower( tbuff[i]) == 'k')
                  {
                  r1 /= AU_IN_KM;
                  i++;
                  if( tolower( tbuff[i]) == 'm')
                     i++;
                  }
               if( !tbuff[i])    /* only one distance entered */
                  r2 = r1;
               else if( sscanf( tbuff + i, "%lf%n", &r2, &j) == 1)
                  {
                  i += j;
                  while( tbuff[i] == ' ')
                     i++;
                  if( tolower( tbuff[i]) == 'k')
                     r2 /= AU_IN_KM;
                  }
               sprintf( message_to_user, "R1 = %lf; R2 = %lf", r1, r2);
               }
            else if( *tbuff == 'g')
               {
               extern double general_relativity_factor;

               general_relativity_factor = atof( tbuff + 1);
               }
            update_element_display = 1;
            break;
         case 's': case 'S':     /* save orbital elements to a file */
            {
            char filename[80];
            FILE *ofile, *ifile;;
            double orbit2[6];

            inquire( "Enter filename for saving elements: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            sscanf( tbuff, "%s", filename);
            ofile = fopen( filename, "wb");
            if( ofile)
               {
               ifile = fopen( elements_filename, "rb");
               while( fgets( tbuff, sizeof( tbuff), ifile))
                  fputs( tbuff, ofile);
               fclose( ifile);
               fclose( ofile);
               }

            memcpy( orbit2, orbit, 6 * sizeof( double));
            integrate_orbit( orbit2, curr_epoch, epoch_shown);
            store_solution( obs, n_obs, orbit2, epoch_shown,
                                          perturbers);
            }
            break;
         case 't': case 'T':
            residual_format ^= RESIDUAL_FORMAT_TIME_RESIDS;
            if( residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
               strcpy( message_to_user, "Showing time/cross-track residuals");
            else
               strcpy( message_to_user, "Showing RA/dec residuals");
            break;
         case 'u': case 'U':
            if( !pop_orbit( &curr_epoch, orbit))
               {
               strcpy( message_to_user, "Last orbit operation undone");
               update_element_display = 1;
               set_locs( orbit, curr_epoch, obs, n_obs);
               }
            else
               strcpy( message_to_user, "No more orbits to undo!");
            break;
         case 'v':         /* apply Vaisala method */
         case 'V':         /* apply Vaisala method, without linearizing */
            {
            double vaisala_dist, angle_param;
            int n_fields, success = 0;

            inquire( "Enter peri/apohelion distance: ", tbuff, sizeof( tbuff),
                                    COLOR_DEFAULT_INQUIRY);
            n_fields = sscanf( tbuff, "%lf,%lf", &vaisala_dist, &angle_param);
            push_orbit( curr_epoch, orbit);
            if( n_fields == 1)      /* simple Vaisala */
               {
               herget_method( obs, n_obs, -vaisala_dist, 0., orbit, NULL, NULL,
                                                   NULL);
               if( c != 'V')
                  adjust_herget_results( obs, n_obs, orbit);
               success = 1;
               }
            else if( n_fields == 2)       /* "trial orbit" method */
               if( !find_trial_orbit( orbit, obs, n_obs,
                                           vaisala_dist, angle_param))
                  {
                  if( c != 'V')
                     adjust_herget_results( obs, n_obs, orbit);
                  success = 1;
                  }
            if( success)
               {
                        /* epoch is that of first included observation: */
               get_epoch_range_of_included_obs( obs, n_obs, &curr_epoch, NULL);
               get_r1_and_r2( n_obs, obs, &r1, &r2);
               update_element_display = 1;
               }
            }
            break;
         case 'w': case 'W':
            i = find_worst_observation( obs, n_obs);
            if( i > -1)
               curr_obs = i;
            break;
         case 'x': case 'X':
            if( obs[curr_obs].is_included)
               obs[curr_obs].is_included = 0;
            else
               obs[curr_obs].is_included = 1;
            strcpy( message_to_user, "Inclusion of observation toggled");
            add_off_on = obs[curr_obs].is_included;
            break;
         case 'y': case 'Y':
            show_a_file( "gauss.out");
            break;
         case 'z': case 'Z':
            {
            double state2[6], delta_squared = 0;

            inquire( "Time span: ", tbuff, sizeof( tbuff),
                              COLOR_DEFAULT_INQUIRY);
            memcpy( state2, orbit, 6 * sizeof( double));
            integrate_orbit( state2, curr_epoch, curr_epoch + atof( tbuff));
            integrate_orbit( state2, curr_epoch + atof( tbuff), curr_epoch);
            for( i = 0; i < 3; i++)
               {
               state2[i] -= orbit[i];
               delta_squared += state2[i] * state2[i];
               }
            sprintf( message_to_user, "Change = %.3e AU = %.3e km",
                              sqrt( delta_squared),
                              sqrt( delta_squared) * AU_IN_KM);
            }
            break;
         case ALT_D:
            inquire( "Debug level: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            debug_level = atoi( tbuff);
            break;
#ifdef __WATCOMC__
         case ALT_T:
#endif
         case '$':
            inquire( "Tolerance: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            if( atof( tbuff))
               {
               extern double integration_tolerance;

               integration_tolerance = atof( tbuff);
               }
            break;
         case '%':
            inquire( "Weight of this observation: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            if( atof( tbuff))
               {
               obs[curr_obs].weight = atof( tbuff);
               sprintf( message_to_user, "Weight reset to %.3e",
                                    obs[curr_obs].weight);
               }
            break;
         case '"':
            debug_mouse_messages ^= 1;
            strcpy( message_to_user, "Mouse debugging");
            add_off_on = debug_mouse_messages;
            break;
         case '@':
            {
            extern int setting_outside_of_arc;

            setting_outside_of_arc ^= 1;
            strcpy( message_to_user, "Setting outside of arc turned");
            add_off_on = setting_outside_of_arc;
            }
            break;
         case '(':
            set_locs( orbit, curr_epoch, obs, n_obs);
            break;
         case ')':
            inquire( "Enter name of file to be displayed: ",
                               tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            show_a_file( tbuff);
            break;
         case '`':
            {
            default_comet_magnitude_type =
                        'N' + 'T' - default_comet_magnitude_type;
            if( default_comet_magnitude_type == 'N')
               strcpy( message_to_user, "Using nuclear mags for comets");
            else
               strcpy( message_to_user, "Using total mags for comets");
            }
         case KEY_MOUSE:   /* already handled above */
            break;
         case 27:
         case 'q': case 'Q':
#ifdef KEY_EXIT
         case KEY_EXIT:
#endif
            quit = 1;
            break;
         case '=':
            residual_format ^= RESIDUAL_FORMAT_MAG_RESIDS;
            strcpy( message_to_user, "Magnitude residual display turned");
            add_off_on = (residual_format & RESIDUAL_FORMAT_MAG_RESIDS);
            break;
         case '+':
            element_format ^= ELEM_OUT_HELIOCENTRIC_ONLY;
            strcpy( message_to_user, (element_format & ELEM_OUT_HELIOCENTRIC_ONLY ?
                     "Heliocentric elements only" :
                     "Planetocentric elements allowed"));
            update_element_display = 1;
            break;
         case '[':
            show_a_file( "covar.txt");
            break;
         case ']':
            residual_format ^= RESIDUAL_FORMAT_SHOW_DELTAS;
            strcpy( message_to_user, "Delta display turned");
            add_off_on = (residual_format & RESIDUAL_FORMAT_SHOW_DELTAS);
            break;
         case '-':
            list_codes ^= 1;
            strcpy( message_to_user, "Codes listing turned");
            add_off_on = list_codes;
            break;
         case '\\':
            inquire( "Enter observatory code and time offset: ",
                               tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY);
            for( i = 0; i < n_obs; i++)
               if( !memcmp( tbuff, obs[i].mpc_code, 3))
                  obs[i].jd += atof( tbuff + 3) / 86400.;
                        /* Make sure obs remain sorted by time: */
            for( i = 0; i < n_obs - 1; i++)
               if( obs[i].jd > obs[i + 1].jd)
                  {
                  OBSERVE temp = obs[i];

                  obs[i] = obs[i + 1];
                  obs[i + 1] = temp;
                  if( i)
                     i -= 2;
                  }
            update_element_display = 1;
            break;
         case ',':
            show_a_file( "debug.txt");
            break;
         case '.':
            strcpy( message_to_user, longname( ));
            break;
         case '~':         /* remove duplicate obs */
            for( i = 0; i < n_obs; i++)
               if( !memcmp( obs + i, obs + i + 1, sizeof( OBSERVE)))
                  {
                  n_obs--;
                  memmove( obs + i, obs + i + 1,
                                 (n_obs - i) * sizeof( OBSERVE));
                  i--;
                  }
            update_element_display = 1;
            break;
         case KEY_MOUSE_MOVE:
         case KEY_TIMER:
            break;
#ifdef __PDCURSES__
         case KEY_RESIZE:
            resize_term( 0, 0);
            sprintf( message_to_user, "KEY_RESIZE: %d x %d",
                        getmaxx( stdscr), getmaxy( stdscr));
            break;
#endif
         case '\'':
            {
            extern int alt_mpcorb;

            alt_mpcorb ^= 1;
            strcpy( message_to_user, "Alternate MPCORB display");
            add_off_on = alt_mpcorb;
            update_element_display = 1;
            }
            break;
         case '{':
            {
            int n_excluded = 0;
            double reject_limit, rms;

            for( i = 0; i < n_obs; i++)
               if( !obs[i].is_included)
                  n_excluded++;
            reject_limit =
                    peirce_rayleigh_func( n_obs, n_excluded + 1, 6);
//          rms = compute_rms( obs, n_obs);
            rms = compute_rms( obs, n_obs) * 1.414213;
            if( filter_obs( obs, n_obs, rms * reject_limit)
                     != FILTERING_CHANGES_MADE)
               strcpy( message_to_user, "No filtering done");
            else
               sprintf( message_to_user, "Rejections at %.3lf sigmas (%.3lf\")",
                           reject_limit, rms * reject_limit);
            }
            break;
         case '&':
            if( !find_precovery_plates( "plates.txt", orbit, curr_epoch))
               show_a_file( "plates.txt");
            break;
         case ';':
            observation_display ^= 8;
            strcpy( message_to_user,"Residual legend");
            add_off_on = (observation_display & 8);
            clear( );
            break;
         case '}':
            {
            precise_residuals ^= 1;
            strcpy( message_to_user, "Precise residuals");
            add_off_on = precise_residuals;
            }
            break;
         case '_':
            link_arcs( obs, n_obs, r1, r2);
            show_a_file( "gauss.out");
            break;
         case '>':
            perturbers = 1;
            break;
         case KEY_F(14):      /* Shift-F2 */
         case KEY_F(15):      /* Shift-F3 */
            for( i = (c == KEY_F(15) ? 0 : curr_obs);
                 i <= (c == KEY_F(15) ? n_obs - 1 : curr_obs); i++)
               {
               MOTION_DETAILS m;

               compute_observation_motion_details( obs + i, &m);
               debug_printf( "Time resid %lf\n", m.time_residual);
               obs[i].jd += m.time_residual / 86400.;
               set_up_observation( obs + i);         /* mpc_obs.cpp */
               }
            break;
         case KEY_F(12):
            {
            int last, first;
            OBSERVE tobs1, tobs2;
            static int desired_soln = 0;

            first = find_first_and_last_obs_idx( obs, n_obs, &last);
            tobs1 = obs[first];
            tobs2 = obs[last];
            if( find_circular_orbits( &tobs1, &tobs2, orbit, desired_soln++))
               strcpy( message_to_user, "No circular orbits found");
            else
               {
               curr_epoch = tobs1.jd - tobs1.r / AU_PER_DAY;
               set_locs( orbit, curr_epoch, obs, n_obs);
               get_r1_and_r2( n_obs, obs, &r1, &r2);
               update_element_display = 1;
               }
            }
            break;
         default:
            debug_printf( "Key %d hit\n", c);
            show_a_file( "dos_help.txt");
            break;
         }
      }
   attrset( COLOR_PAIR( COLOR_BACKGROUND));
   show_final_line( obs, n_obs, curr_obs, COLOR_BACKGROUND);
#ifdef __PDCURSES__
   if( !strstr( longname( ), "Win32a"))
      resize_term( original_ymax, original_xmax);
#endif
   endwin( );
   free_weight_recs( );
   create_obs_file( obs, n_obs, 0);

   sprintf( tbuff, "%s %d %d %d", mpc_code,
             observation_display, residual_format, list_codes);
   set_environment_ptr( "CONSOLE_OPTS", tbuff);
   sprintf( tbuff, "%c,%d,%d,%d,%lf,%lf,%d",
               default_comet_magnitude_type,
               element_format, element_precision,
               ephemeris_output_options,
               max_residual_for_filtering, noise_in_arcseconds,
               precise_residuals);
   set_environment_ptr( "SETTINGS", tbuff);
   store_defaults( ephemeris_output_options);             /* elem_out.c */
   set_environment_ptr( "EPHEM_START", ephemeris_start);
   sprintf( tbuff, "%d %s", n_ephemeris_steps, ephemeris_step_size);
   set_environment_ptr( "EPHEM_STEPS", tbuff);
   free( ids);
#ifdef XCURSES
   XCursesExit();
#endif
   return( 0);
}

#ifdef PURE_WINDOWS_VERSION

int APIENTRY WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance,
                    LPSTR lpszCmdLine, int nCmdShow)
{
   char *argv[30];
   int i, argc = 1;

   argv[0] = "Find_Orb";
   for( i = 0; lpszCmdLine[i]; i++)
       if( lpszCmdLine[i] != ' ' && (!i || lpszCmdLine[i - 1] == ' '))
          argv[argc++] = lpszCmdLine + i;

   for( i = 0; lpszCmdLine[i]; i++)
       if( lpszCmdLine[i] == ' ')
          lpszCmdLine[i] = '\0';

   return dummy_main( argc, (const char **)argv);
}
#endif
