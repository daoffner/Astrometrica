#include <stdio.h>
#include "watdefs.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define DegreesToRadians( x)   ((x) * PI / 180.)
#define RadiansToDegrees( x)   ((x) * 180. / PI)

int main( const int argc, const char **argv)
{
//---------------------
    // POLARIS @ J2000
    const double RAH = 2;
    const double RAM = 31;
    const double RAS = 49.09;

    const double DECD = 89;
    const double DECM = 15;
    const double DECS =  50.8;

    const double RA = RAH + (RAM / 60.0) + (RAS / 3600.0);
    const double DEC = DECD + (DECM / 60.0) + (DECS / 3600.0);
    const double fJD = 2456019.08790; // UTC julian date

    double fHA;
    DPT alt_az, loc_epoch, ra_dec, latlon;

    const double DEC_in_RAD = DegreesToRadians(DEC);
    const double LAT_in_RAD = DegreesToRadians(51.60138);
    const double LON_in_RAD = DegreesToRadians(-2.86391);
    const double RA_as_DEG_in_RAD = DegreesToRadians(RA * 15.0); // *15.0 put 24h into 360deg

    ra_dec.x = -RA_as_DEG_in_RAD;   /* NOTE THE NEGATION OF THE RA! */
    ra_dec.y = DEC_in_RAD;

    loc_epoch.x = 0;
    loc_epoch.y = 0;

    latlon.x = LON_in_RAD;
    latlon.y = LAT_in_RAD;


    full_ra_dec_to_alt_az(&ra_dec, &alt_az, &loc_epoch, &latlon, fJD, &fHA);
    printf( "az: %lf  alt: %lf\n",  RadiansToDegrees( alt_az.x),
                                    RadiansToDegrees( alt_az.y));
    printf( "Hour angle: %lf\n", RadiansToDegrees( fHA));
    return( 0);
}
