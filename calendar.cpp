#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "date.h"

#define X0 20
#define XSIZE 105
#define XEND (X0 + 7 * XSIZE)
#define TEXT_XOFFSET 3

#define Y0 20
#define YSIZE 104
#define YEND (Y0 + 5 * YSIZE)
#define TEXT_YOFFSET 15

#define DAY_OF_WEEK_HT 20
#define TOP_OF_DAYOFWK (Y0 + 5 * YSIZE + DAY_OF_WEEK_HT)

static const char *months[12] =  { "January", "February", "March",
               "April", "May", "June", "July", "August", "September",
               "October", "November", "December" };
static int show_jd_values = 0, dollhouse = 0, single_page = 0;

#ifdef RESTORE_EASTER
void easter_date( const long year, int *month, int *day);
#endif

static void show_month_text( const int month, const int year)
{
   char buff[80];

   sprintf( buff, "%s %d", months[month - 1], year);
   if( show_jd_values)
      sprintf( buff + strlen( buff), " (JD %ld.5)",
                  dmy_to_day( 0, month, (long)year, CALENDAR_JULIAN_GREGORIAN) - 1);
   printf( "%d %d moveto (%s) show\n",
            X0 + 7 * XSIZE / 2 - (int)strlen( buff) * 20 / 3,
            TOP_OF_DAYOFWK + 5, buff);
}

static void show_small_month( int month, int year, const int x0, int y0)
{
   long jd1, jd2;
   int n_days, starting_loc, i, j, day;
   char buff[100];

   while( month <= 0)
      {
      month += 12;
      year--;
      }
   while( month >= 13)
      {
      month -= 12;
      year++;
      }
   jd1 = dmy_to_day( 0, month, (long)year, CALENDAR_JULIAN_GREGORIAN);
   jd2 = dmy_to_day( 0, month + 1, (long)year, CALENDAR_JULIAN_GREGORIAN);
   n_days = (int)( jd2 - jd1);
   starting_loc = (int)( (jd1 + 2L) % 7L);
   strcpy( buff, months[month - 1]);
   sprintf( buff + strlen( buff), " %d", year);
   y0 -= 12;
   printf( "%d %d moveto (%s) show\n",
            x0 + (28 - (int)strlen( buff)) * TEXT_XOFFSET / 2, y0, buff);
   for( i = 0; i < 6; i++)
      {
      y0 -= 12;
      memset( buff, ' ', 20);
      for( j = 0; j < 7; j++)
         {
         day = i * 7 + j - starting_loc + 1;
         if( day >= 1 && day <= n_days)
            {
            buff[j * 3] = (char)( '0' + day / 10);
            buff[j * 3 + 1] = (char)( '0' + day % 10);
            if( day < 10)
               buff[j * 3] = ' ';
            }
         }
      for( j = 20; j && buff[j - 1] == ' '; j--)
         ;
      buff[j] = buff[20] ='\0';
      if( j)
         printf( "%d %d moveto (%s) show\n", x0, y0, buff);
      }
   if( show_jd_values)
      printf( "%d %d moveto ((JD %ld.5)) show\n", x0, y0 - 12, jd1 - 1);
}

int get_phase_data( double t, double *t_phases)
{
   const char *phase_file_name = "z:\\vsop\\phases.dat";
   FILE *phase_file;
   int rval = -1;

   phase_file = fopen( phase_file_name, "rb");
   if( phase_file)
      {
      long dt[12];
      long offset;
      double k;
      int i;
      const long phase0 = -55683L;

      k = floor( (t - 2451550.09765) / 29.530588853 - .5);
      offset = (4L * ((long)k - phase0) + 1L) * (long)sizeof( long);
      if( offset > 0L && offset < 1738116L - 48L)
         {
         fseek( phase_file, offset, SEEK_SET);
         if( !fread( dt, 12, sizeof( long), phase_file))
            {
            printf( "Read error on phase file\n");
            exit( -1);
            }
         rval = 0;

         for( i = 0; i < 12; i++, k += .25)
            {
            double t_phase1 = ( 2451550.09765 + 29.530588853 * k);

            t_phases[i] = t_phase1 + (double)dt[i] / 86400.;
            }
         }
      fclose( phase_file);
      }
   return( rval);
}

char **grab_dates( const int month, const int year)
{
   FILE *dates_file = fopen( "dates.txt", "rb");
   char **rval = (char **)calloc( 100, sizeof( char *));
   int n_found = 0, day, i, dst_day = 0;
   char buff[80];
   long jd;
   double t_phases[12];

   while( fgets( buff, sizeof( buff), dates_file))
      if( !memcmp( months[month - 1], buff, 3))
         {
         char *tptr = strstr( buff, " WARNING");

         day = atoi( buff + 3);
         if( tptr)
            *tptr = '\0';
         for( i = 0; buff[i] >= ' '; i++)
            ;
         buff[i] = '\0';
         rval[n_found] = (char *)malloc( i - 4);
         rval[n_found][0] = (char)day;
         strcpy( rval[n_found] + 1, buff + 7);
         n_found++;
         }
   fclose( dates_file);

#ifdef RESTORE_EASTER
   for( i = 0; i < 3; i++)
#else
   for( i = 0; i < 2; i++)
#endif
      {
      static const char *names[3] = { "Rosh Hashanah", "Yom Kippur", "Easter" };
      int tmonth;

#ifdef RESTORE_EASTER
      if( i == 2)        /* easter */
         easter_date( (long)year, &tmonth, &day);
      else            /* Rosh Hashanah/Yom Kippur */
#endif
         {
         long tyear;

         jd = dmy_to_day( (i ? 10 : 1), 1,
                             (long)year + 3761L, CALENDAR_HEBREW);
         day_to_dmy( jd, &day, &tmonth, &tyear, CALENDAR_JULIAN_GREGORIAN);
         }
      if( tmonth == month)
         {
         rval[n_found] = (char *)malloc( strlen( names[i]) + 2);
         rval[n_found][0] = (char)day;
         strcpy( rval[n_found] + 1, names[i]);
         n_found++;
         }
      }
            /* 2006 and before:  DST change times are first Sunday in */
            /* April,  last Sun in October.  2007 and later:  second  */
            /* Sunday in March,  first Sunday in November.            */
   if( year < 2007 && (month == 4 || month == 10))
      {
      const int day = year + year / 4 - year / 100 + year / 400;

      if( month == 4)
         dst_day = 7 - (day + 5) % 7;
      else
         dst_day = 31 - (day + 2) % 7;
      }
   if( year > 2006 && (month == 3 || month == 11))  /* incl start/end of D"S"T */
      {
      const int day = year + year / 4 - year / 100 + year / 400;

      if( month == 3)
         dst_day = 14 - (day + 2) % 7;
      else
         dst_day = 7 - (day + 2) % 7;
      }
   if( dst_day)
      {
      rval[n_found] = (char *)malloc( 20);
      rval[n_found][0] = (char)dst_day;
      strcpy( rval[n_found] + 1, (month < 6 ? "D'S'T begins" : "D'S'T ends"));
      n_found++;
      }

   jd = dmy_to_day( 0, month, year, CALENDAR_JULIAN_GREGORIAN);
   if( !get_phase_data( jd - 10, t_phases))
      for( i = 0; i < 12; i++)
         {
         day = (int)( (long)floor( t_phases[i] + .5) - jd);
         if( day > 0 && day < 32)
            {
            const char *phases[4] = { "*New moon", "*First quarter",
                       "*Full moon", "*Last quarter" };
            const char *phasestr = phases[i & 3];

            rval[n_found] = (char *)malloc( strlen( phasestr) + 2);
            rval[n_found][0] = (char)day;
            strcpy( rval[n_found] + 1, phasestr);
            n_found++;
            }
         }
   return( rval);
}

static const char *page_orientation_landscape = "%%PageOrientation: Landscape";
static const char *page_orientation_portrait = "%%PageOrientation: Portrait";

static const char *trailer_data[51] = {
                "%%Page: 1 1",
                page_orientation_landscape,
                "gsave calendar grestore",
                "",
                "%%Page: 2 2",
                page_orientation_portrait,
                "  gsave",
                "  /doublescale 1.34 def",
                "  -90 rotate -812 0 doublescale mul translate",
                "  doublescale doublescale scale",
                "  calendar grestore",
                "",
                "%%Page: 3 3",
                page_orientation_portrait,
                "  gsave",
                "  /doublescale 1.34 def",
                "  -90 rotate -812 -350 doublescale mul translate",
                "  doublescale doublescale scale",
                "  calendar grestore",
                "",
                "%%Page: 4 4",
                page_orientation_landscape,
                "  gsave",
                "  /quadscale 2 def",
                "  quadscale quadscale scale",
                "  calendar grestore",
                "",
                "%%Page: 5 5",
                page_orientation_landscape,
                "  gsave",
                "  /quadscale 2 def",
                "  0 -370 quadscale mul translate",
                "  quadscale quadscale scale",
                "  calendar grestore",
                "",
                "%%Page: 6 6",
                page_orientation_landscape,
                "  gsave",
                "  /quadscale 2 def",
                "  -290 quadscale mul -370 quadscale mul translate",
                "  quadscale quadscale scale",
                "  calendar grestore",
                "",
                "%%Page: 7 7",
                page_orientation_landscape,
                "  gsave",
                "  /quadscale 2 def",
                "  -290 quadscale mul 0 quadscale mul translate",
                "  quadscale quadscale scale",
                "  calendar grestore",
                NULL };


static const char *dollhouse_trailer_data[51] = {
                "%%Page: 1 1",
                page_orientation_portrait,
                "  gsave",
                "  /shrinkscale .2 def",
                "  -90 rotate -812 0 translate",
                "  shrinkscale shrinkscale scale",
                "  calendar grestore",
                NULL };

#define FONT_UNSET -1
#define FONT_ITALIC 1
#define FONT_PLAIN  2

int calendar( const int month, const int year)
{
   int i, j;
   const long jd1 = dmy_to_day( 0, month, (long)year, CALENDAR_JULIAN_GREGORIAN);
   const long jd2 = dmy_to_day( 0, month + 1, (long)year, CALENDAR_JULIAN_GREGORIAN);
   const int n_days = (int)( jd2 - jd1), starting_loc = (int)( (jd1 + 1L) % 7L);
   char buff[100];
   char lines_used[35], phases_shown[35];
   char **dates = grab_dates( month, year);


   printf( "%%!PS-Adobe-2.0\n");
   printf( "%%%%Pages: %d\n", (dollhouse || single_page ? 1 : 7));
   printf( "%%%%PageOrder: Ascend\n");
// printf( "%%%%BoundingBox: 48 69 556 746\n");
   printf( "%%%%Creator: calendar.cpp\n");
   printf( "%%%%Copyright: none\n");
   printf( "%%%%Title: Calendar for %s %d\n", months[month - 1], year);
   printf( "%%%%Version: none\n");
   printf( "%%%%DocumentData: Clean7Bit\n");
   printf( "%%%%EndComments\n");
   printf( "%%%%BeginDefaults\n");
   printf( "%%%%PageResources: font Times-Roman\n");
   printf( "%%%%PageResources: font Times-Italic\n");
   printf( "%%%%PageResources: font Courier-Bold\n");
   printf( "%%%%EndDefaults\n");
   printf( "/blood {\n");
   printf( "/aa 2 def /bb 1.7 def\n");
   printf( "aa 0 rmoveto bb 0 rlineto 0 aa rlineto aa 0 rlineto 0 bb rlineto\n");
   printf( "0 aa sub 0 rlineto 0 aa rlineto 0 bb sub 0 rlineto\n");
   printf( "0 0 aa sub rlineto 0 aa sub 0 rlineto\n");
   printf( "0 0 bb sub rlineto aa 0 rlineto 0 0 aa sub rlineto aa 3 mul 0 rmoveto\n");
   printf( "} def\n\n");
   printf( "/calendar {\n");
   printf( "612 0 translate 90 rotate\n");
   for( i = 0; i <= 5; i++)             /* horizontal lines separating weeks */
      printf( "%d %d moveto %d %d lineto\n", X0, Y0 + i * YSIZE,
                                           XEND, Y0 + i * YSIZE);
   printf( "%d %d moveto %d %d lineto\n", X0, TOP_OF_DAYOFWK,
                                        XEND, TOP_OF_DAYOFWK);

   for( i = 0; i <= 7; i++)       /* vertical lines */
      printf( "%d %d moveto %d %d lineto\n", X0 + i * XSIZE, Y0,
                                             X0 + i * XSIZE, YEND);
   printf( "/defaultfontsize  { 12 scalefont } def\n");
   printf( "/Times-Roman findfont defaultfontsize setfont\n");
   for( i = 0; i < 35; i++)
      lines_used[i] = phases_shown[i] = 0;

   for( i = 1; i <= n_days; i++)
      {
      int xcell = (i + starting_loc) % 7, curr_font = FONT_UNSET;
      int ycell = 4 - (i + starting_loc) / 7, double_cell = 0, x0, y0;
      int cell_idx;

      if( starting_loc == 6)        /* first of month is on a sunday */
         ycell++;
      if( !ycell && i + 7 <= n_days)     /* two days in one cell */
         {
         sprintf( buff, "%d/%d", i, i + 7);
         double_cell = 1;
         }
      else
         sprintf( buff, "%d", i);
      x0 = X0 + xcell * XSIZE;
      y0 = Y0 + (ycell + 1) * YSIZE;
      if( ycell >= 0)
         {
         int xboxsize = TEXT_XOFFSET + 8 * strlen( buff);

         printf( "%d %d moveto (%s) show\n",
                     x0 + TEXT_XOFFSET, y0 - TEXT_YOFFSET, buff);
         printf( "%d %d moveto %d %d rlineto %d %d rlineto\n",
               x0, y0 - (TEXT_YOFFSET + 2),
               xboxsize, 0, 0, TEXT_YOFFSET + 2);
         }
      else
         {
         double_cell = 1;
         ycell = 0;
         y0 = Y0 + YSIZE;
         }
      cell_idx = (4 - ycell) * 7 + xcell;
      for( j = 0; dates[j]; j++)
         if( dates[j][0] == i)
            {
            char *text_to_show = dates[j] + 1;
            const int is_lunar_phase = (*text_to_show == '*');
            const int font_to_use = (is_lunar_phase ? FONT_ITALIC : FONT_PLAIN);
            int xloc = X0 + xcell * XSIZE + TEXT_XOFFSET +
                           ((is_lunar_phase && !double_cell) ? 24 : 0);
            int yloc = (is_lunar_phase ?
                         y0 - TEXT_YOFFSET * (1 + double_cell)
                                        - phases_shown[cell_idx] * 11 :
                         y0 - YSIZE + 3 + lines_used[cell_idx] * 11);

            if( curr_font != font_to_use)
               {
               curr_font = font_to_use;
               if( curr_font == FONT_ITALIC)
                  printf( "/Times-Italic findfont 9 scalefont setfont\n");
               if( curr_font == FONT_PLAIN)
                  printf( "/Times-Roman findfont 9 scalefont setfont\n");
               }
            if( is_lunar_phase)     /* lunar phase */
               {
               text_to_show++;
               }

            printf( "%d %d moveto ", xloc, yloc);
            if( double_cell)
               printf( "(\050%d\051 ) show ", i);
            if( *text_to_show == 'B' && text_to_show[1] == ' ')
               {
               printf( "blood ");
               text_to_show += 2;
               }
            printf( "(%s) show\n", text_to_show);
            if( is_lunar_phase)
               phases_shown[cell_idx]++;
            else
               lines_used[cell_idx]++;
            }
      if( !lines_used[cell_idx])         /* mark it as being a "used" cell */
         lines_used[cell_idx] = 1;
      if( curr_font != FONT_UNSET)
         printf( "/Times-Roman findfont defaultfontsize setfont\n");
      }

   for( i = 0; i < 7; i++)
      {
      const char *day_of_week_text[7] = { "Sunday", "Monday", "Tuesday",
                           "Wednesday", "Thursday", "Friday", "Saturday" };

      printf( "%d %d moveto (%s) show\n",
                              X0 + i * XSIZE + XSIZE / 2 -
                              (int)strlen( day_of_week_text[i]) * 4,
                              Y0 + 5 * YSIZE + 5,
                              day_of_week_text[i]);
      }

   printf( "/Times-Roman findfont 24 scalefont setfont\n");

   show_month_text( month, year);

   printf( "/Courier-Bold findfont 8 scalefont setfont\n");

               /* show two preceding & two (or more) following months: */
   for( i = j = 0; i < 35; i++)
      if( !lines_used[i])
         {
         show_small_month( month + j - 2, year,
                     X0 + (i % 7) * XSIZE + 1, Y0 + YSIZE * (5 - i / 7));
         j++;
         if( j == 2)    /* avoid duplicating the 'main' month */
            j = 3;
         }
   printf( "stroke showpage } def\n");
   if( dollhouse)
      for( i = 0; dollhouse_trailer_data[i]; i++)
         printf( "%s\n", dollhouse_trailer_data[i]);
   else
      for( i = 0; (single_page ? i < 3 : trailer_data[i] != NULL); i++)
         printf( "%s\n", trailer_data[i]);
   return( 0);
}

int main( const int argc, const char **argv)
{
   int month = atoi( argv[1]), year = atoi( argv[2]);
   int i;

   if( month > year)       /* allow for possible reversed entry */
      {
      month = year;
      year = atoi( argv[1]);
      }
   for( i = 3; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'j':
               show_jd_values = 1;
               break;
            case 'd':
               dollhouse = 1;
               break;
            case 's':
               single_page = 1;
               break;
            }
   return( calendar( month, year));
}
