#define OBSERVE struct observe

OBSERVE
   {
   double jd, obs_posn[3], obs_vel[3], vect[3], ra, dec, obs_mag;
   double r,  obj_posn[3], obj_vel[3], solar_r, computed_ra, computed_dec;
   double weight, mag_weight, computed_mag;
   int flags, is_included;
   int time_precision, ra_precision, dec_precision, mag_precision;
   char mpc_code[4], packed_id[13], reference[6];
   char mag_band, mag_band2, discovery_asterisk, note1, note2, satellite_obs;
   char *second_line;
   };

#define OBJECT_INFO struct object_info

OBJECT_INFO
   {
   char packed_desig[13], obj_name[80];
   char solution_exists;
   int n_obs;
   double jd_start, jd_end;
   unsigned long file_offset;
   };

#define MOTION_DETAILS struct motion_details

MOTION_DETAILS
   {
   double xresid, yresid;
   double ra_motion, dec_motion, total_motion;     /* in arcmin/hour */
   double position_angle_of_motion;                /* in degrees,  0-360 */
   double radial_vel;                              /* in km/s */
   double time_residual;                           /* in seconds */
   double cross_residual;                          /* in arcseconds */
   };

/* Bitfield options for ephemeris_in_a_file( ): */
/* Bottom three bits define an ephemeris type.  "Observables" are the */
/* usual RA/dec,  radial velocity,  etc. type output.  "State vector  */
/* output" results in the time (as a JD) and position (in AU,  relative */
/* to the observer,  in Cartesian coordinates) being output.  "Position */
/* output" is the same thing,  minus the velocity components.  "MPCORB  */
/* output" means that the orbital elements will be written out for each */
/* ephemeris step,  on a single line.  "8-line output" is almost the    */
/* same,  except that the elements are written out in the MPC's usual   */
/* eight-line form.  "Close approaches" will result in the range minima */
/* (times and distances) being output.                                  */

#define OPTION_OBSERVABLES             0
#define OPTION_STATE_VECTOR_OUTPUT     1
#define OPTION_POSITION_OUTPUT         2
#define OPTION_MPCORB_OUTPUT           3
#define OPTION_8_LINE_OUTPUT           4
#define OPTION_CLOSE_APPROACHES        5

#define OPTION_ALT_AZ_OUTPUT           8
#define OPTION_RADIAL_VEL_OUTPUT     0x010
#define OPTION_MOTION_OUTPUT         0x020
#define OPTION_PHASE_ANGLE_OUTPUT    0x040
#define OPTION_GROUND_TRACK          0x100
#define OPTION_SEPARATE_MOTIONS      0x200

#define OPTION_ROUND_TO_NEAREST_STEP 0x400
#define OPTION_PHASE_ANGLE_BISECTOR  0x800
#define OPTION_HELIO_ECLIPTIC        0x1000
#define OPTION_TOPO_ECLIPTIC         0x2000

#define HELIOCENTRIC_SIGMAS_ONLY       0
#define ORBIT_SIGMAS_REQUESTED         1
#define NO_ORBIT_SIGMAS_REQUESTED    (-1)

/* Bitfields for 'flags' parameter of OBSERVE */
#define OBS_IS_COMET       1

#ifdef SEEK_CUR
OBSERVE FAR *load_observations( FILE *ifile, const char *packed_desig,
                        const int n_obs);
#endif
OBJECT_INFO *find_objects_in_file( const char *filename,
                                         int *n_found, const char *station);
int get_object_name( char *obuff, const char *packed_desig);
int get_observer_data( const char FAR *mpc_code, char *buff,
           double *lon_in_radians, double *rho_cos_phi, double *rho_sin_phi);
void recreate_observation_line( char *obuff,
                                   const OBSERVE FAR *obs);    /* ephem0.cpp */
void put_observer_data_in_text( const char FAR *mpc_code, char *buff);

void create_obs_file( const OBSERVE FAR *obs, int n_obs, const int append);
int find_worst_observation( const OBSERVE FAR *obs, const int n_obs);
double calc_absolute_magnitude( OBSERVE FAR *obs, int n_obs);
double compute_rms( const OBSERVE FAR *obs, const int n_obs);
int herget_method( OBSERVE FAR *obs, int n_obs, double r1, double r2,
         double *orbit, double *d_r1, double *d_r2, const char *limited_orbit);
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
void improve_parabolic( OBSERVE FAR *obs, int n_obs, double *orbit, double epoch);
int full_improvement( OBSERVE FAR *obs, int n_obs, double *orbit,
            const double epoch, const char *limited_orbit,
                 int sigmas_requested);
int set_locs( const double *orbit, const double t0, OBSERVE FAR *obs,
                                   const int n_obs);
void make_date_range_text( char *obuff, const double jd1, const double jd2);
                                                        /* orb_func.cpp */
int get_r1_and_r2( const int n_obs, const OBSERVE FAR *obs,
                             double *r1, double *r2);    /* elem_out.cpp */
int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                  int *idx1, int *idx2);  /* elem_out.cpp */
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
double get_step_size( const char *stepsize, char *step_units,
                                 int *step_digits);          /* ephem0.cpp */
int ephemeris_in_a_file( const char *filename, const double *orbit,
         OBSERVE *obs, const int n_obs,
         const int planet_no,
         const double epoch_jd, const double jd_start, const char *stepsize,
         const double lon,
         const double rho_cos_phi, const double rho_sin_phi,
         const int n_steps, const char *note_text,
         const int options);
int obj_desig_to_perturber( const char *packed_desig);     /* runge.cpp */
int find_best_fit_planet( const double jd, const double *ivect,
                     double *rel_vect);     /* runge.cpp */
int integrate_orbit( double *orbit, const double t0, const double t1);
int generate_observation_text( const OBSERVE FAR *obs, const int n_obs,
                      const int obs_idx, const int line_no, char *buff);
double convenient_gauss( const OBSERVE FAR *obs, int n_obs, double *orbit,
                  const double mu, const int desired_soln); /* gauss.cpp */
void set_solutions_found( OBJECT_INFO *ids, const int n_ids); /* orb_func.c */
int fetch_previous_solution( OBSERVE FAR *obs, const int n_obs, double *orbit,
               double *orbit_epoch, int *perturbers);
int store_solution( const OBSERVE FAR *obs, const int n_obs, const double *orbit,
       const double orbit_epoch, const int perturbers);
int compute_observation_motion_details( const OBSERVE FAR *obs,
               MOTION_DETAILS *m);                    /* mpc_obs.cpp */
int get_findorb_text( char *buff, const int ival);    /* ephem.cpp */
int create_mpc_packed_desig( char *packed_desig, const char *obj_name);
int write_out_elements_to_file( const double *orbit,
            const double curr_epoch,
            const double epoch_shown,
            OBSERVE FAR *obs, const int n_obs, const char *constraints,
            const int precision, const int monte_carlo,
            const int options);    /* elem_out.cpp */

#define ELEM_OUT_ALTERNATIVE_FORMAT    4
#define ELEM_OUT_NO_COMMENT_DATA       2
#define ELEM_OUT_HELIOCENTRIC_ONLY     1

/*
   Lowest two bits of the residual_format field:
      resid_format = 0 -> full-line output without tabs;
      resid_format = 1 -> full-line output with tabs;
      resid_format = 2 -> short MPC-like output,  CYY res res form
      resid_format = 3 -> standard 80-column format

   if( resid_format & 4),  four-digit year
   if( resid_format & 8),  residuals are expressed in time and cross-track
         instead of the default residuals in RA and dec
   if( resid_format & 16),  date is expressed as HH:MM:SS (instead of as a
         decimal fraction of a day)
   if( resid_format & 32),  magnitude residuals are shown instead of posn
*/

#define RESIDUAL_FORMAT_FULL_NO_TABS          0
#define RESIDUAL_FORMAT_FULL_WITH_TABS        1
#define RESIDUAL_FORMAT_SHORT                 2
#define RESIDUAL_FORMAT_80_COL                3
#define RESIDUAL_FORMAT_FOUR_DIGIT_YEARS      4
#define RESIDUAL_FORMAT_TIME_RESIDS           8
#define RESIDUAL_FORMAT_HMS                  16
#define RESIDUAL_FORMAT_MAG_RESIDS           32

int write_residuals_to_file( const char *filename, const char *ast_filename,
        const int n_obs, const OBSERVE FAR *obs_data, const int resid_format);
void format_observation( const OBSERVE FAR *obs, char *text,
                                   const int resid_format);   /* ephem0.cpp */

#define MPC_STATION struct mpc_station

MPC_STATION
   {
   char code[4];
   int color;
   int score;
   };

int find_mpc_color( const MPC_STATION *sdata, const char *mpc_code);
MPC_STATION *find_mpc_color_codes( const int n_obs, const OBSERVE FAR *obs,
                   const int max_n_colors);           /* elem_out.cpp */

#define FILTERING_CHANGES_MADE            1
#define FILTERING_NO_CHANGES_MADE         2
#define FILTERING_FAILED                  3

int filter_obs( OBSERVE FAR *obs, const int n_obs,           /* orb_fun2.cpp */
                  const double max_residual_in_arcseconds);

    /* Functions used to store and restore orbits for the 'undo' function. */
    /* Orbits are stored on a stack and can be retrieved from it.          */
void push_orbit( const double epoch, const double *orbit); /* orb_fun2.cpp */
int pop_orbit( double *epoch, double *orbit);              /* orb_fun2.cpp */
void pop_all_orbits( void);                                /* orb_fun2.cpp */
