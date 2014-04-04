/****************************************************************************/
/*                                                                          */
/*     NGDC's Geomagnetic Field Modeling software for the IGRF and WMM      */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Disclaimer: This program has undergone limited testing. It is        */
/*     being distributed unoffically. The National Geophysical Data         */
/*     Center does not guarantee it's correctness.                          */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 6.0:                                                         */
/*     Bug fixes for the interpolation between models. Also added warnings  */
/*     for declination at low H and corrected behaviour at geogr. poles.    */
/*     Placed print-out commands into separate routines to facilitate       */
/*     fine-tuning of the tables                                            */
/*                                          Stefan Maus 8-24-2004           */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      This program calculates the geomagnetic field values from           */
/*      a spherical harmonic model.  Inputs required by the user are:       */
/*      a spherical harmonic model data file, coordinate preference,        */
/*      altitude, date/range-step, latitude, and longitude.                 */
/*                                                                          */
/*         Spherical Harmonic                                               */
/*         Model Data File       :  Name of the data file containing the    */
/*                                  spherical harmonic coefficients of      */
/*                                  the chosen model.  The model and path   */
/*                                  must be less than PATH chars.           */
/*                                                                          */
/*         Coordinate Preference :  Geodetic (WGS84 latitude and altitude   */
/*                                  above mean sea level (WGS84),           */
/*                                  or geocentric (spherical, altitude      */
/*                                  measured from the center of the Earth). */
/*                                                                          */
/*         Altitude              :  Altitude above sea level.               */
/*                                  if geocentric coordinate preference is  */
/*                                  used then the altitude must be in the   */
/*                                  range of 6370.20 km - 6971.20 km as     */
/*                                  measured from the center of the earth.  */
/*                                  Enter altitude in kilometers, meters,   */
/*                                  or feet                                 */
/*                                                                          */
/*         Date                  :  Date, in decimal years, for which to    */
/*                                  calculate the values of the magnetic    */
/*                                  field.  The date must be within the     */
/*                                  limits of the model chosen.             */
/*                                                                          */
/*         Latitude              :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for northern    */
/*                                  hemisphere, negative for the southern   */
/*                                  hemisphere.                             */
/*                                                                          */
/*         Longitude             :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for eastern     */
/*                                  hemisphere, negative for the western    */
/*                                  hemisphere.                             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      Subroutines called :  degrees_to_decimal,julday,getshc,interpsh,    */
/*                            extrapsh,shval3,dihf,safegets                 */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>            
#include <string.h>
#include <ctype.h>

/* The following include file must define a function 'isnan' */
/* This function, which returns '1' if the number is NaN and 0*/
/* otherwise, could be hand-written if not available. */
/* Comment out one of the two following lines, as applicable */
#include <math.h>               /* for gcc */
#include <mathimf.h>            /* for Intel icc */

#define NaN log(-1.0)

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#define IEXT 0
#define FALSE 0
#define TRUE 1                  /* constants */
#define RECL 81

#define MAXINBUFF RECL+14

/** Max size of in buffer **/

#define MAXREAD MAXINBUFF-2
/** Max to read 2 less than total size (just to be safe) **/

#define MAXMOD 30
/** Max number of models in a file **/

#define PATH MAXREAD
/** Max path and filename length **/

#define EXT_COEFF1 (float)0
#define EXT_COEFF2 (float)0
#define EXT_COEFF3 (float)0

#define MAXDEG 13
#define MAXCOEFF (MAXDEG*(MAXDEG+2)+1) /* index starts with 1!, (from old Fortran?) */
float gh1[MAXCOEFF];
float gh2[MAXCOEFF];
float gha[MAXCOEFF];              /* Geomag global variables */
float ghb[MAXCOEFF];
float d=0,f=0,h=0,i=0;
float dtemp,ftemp,htemp,itemp;
float x=0,y=0,z=0;
float xtemp,ytemp,ztemp;

FILE *stream = NULL;                /* Pointer to specified model data file */


/****************************************************************************/
/*                                                                          */
/*                             Program Geomag                               */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      This program, originally written in FORTRAN, was developed using    */
/*      subroutines written by                                              */
/*      A. Zunde                                                            */
/*      USGS, MS 964, Box 25046 Federal Center, Denver, Co.  80225          */
/*      and                                                                 */
/*      S.R.C. Malin & D.R. Barraclough                                     */
/*      Institute of Geological Sciences, United Kingdom.                   */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      Translated                                                          */
/*      into C by    : Craig H. Shaffer                                     */
/*                     Lockheed Missiles and Space Company                  */
/*                     Sunnyvale, California                                */
/*                     (408) 756 - 5210                                     */
/*                     29 July, 1988                                        */
/*                                                                          */
/*      Contact      : Greg Fountain                                        */
/*                     Lockheed Missiles and Space Company                  */
/*                     Dept. MO-30, Bldg. 581                               */
/*                     1111 Lockheed Way                                    */
/*                     P.O. Box 3504                                        */
/*                     Sunnyvale, Calif.  94088-3504                        */
/*                                                                          */
/*      Rewritten by : David Owens                                          */
/*                     dio@ngdc.noaa.gov                                    */
/*                     For Susan McClean                                    */
/*                                                                          */
/*      Contact      : Susan McLean                                         */
/*                     sjm@ngdc.noaa.gov                                    */
/*                     National Geophysical Data Center                     */
/*                     World Data Center-A for Solid Earth Geophysics       */
/*                     NOAA, E/GC1, 325 Broadway,                           */
/*                     Boulder, CO  80303                                   */
/*                                                                          */
/*      Original                                                            */
/*      FORTRAN                                                             */
/*      Programmer   : National Geophysical Data Center                     */
/*                     World Data Center-A for Solid Earth Geophysics       */
/*                     NOAA, E/GC1, 325 Broadway,                           */
/*                     Boulder, CO  80303                                   */
/*                                                                          */
/*      Contact      : National Geophysical Data Center                     */
/*                     World Data Center-A for Solid Earth Geophysics       */
/*                     NOAA, E/GC1, 325 Broadway,                           */
/*                     Boulder, CO  80303                                   */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      dio modifications include overhauling the interactive interface to  */
/*  support platform independence and improve fatal error dectection and    */
/*  prevention.  A command line interface was added and the interactive     */
/*  interface was streamlined.  The function safegets() was added and the   */
/*  function getshc's i/0 was modified.  A option to specify a date range   */
/*  range and step was also added.                                          */
/*                                                                          */
/*                                             dio 6-6-96                   */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      Some variables used in this program                                 */
/*                                                                          */
/*    Name         Type                    Usage                            */
/* ------------------------------------------------------------------------ */
/*                                                                          */
/*   a2,b2      Scalar Float           Squares of semi-major and semi-minor */
/*                                     axes of the reference spheroid used  */
/*                                     for transforming between geodetic or */
/*                                     geocentric coordinates.              */
/*                                                                          */
/*   minalt     Float array of MAXMOD  Minimum height of model.             */
/*                                                                          */
/*   altmin     Float                  Minimum height of selected model.    */
/*                                                                          */
/*   altmax     Float array of MAXMOD  Maximum height of model.             */
/*                                                                          */
/*   maxalt     Float                  Maximum height of selected model.    */
/*                                                                          */
/*   d          Scalar Float           Declination of the field from the    */
/*                                     geographic north (deg).              */
/*                                                                          */
/*   sdate  Scalar Float           start date inputted                      */
/*                                                                          */
/*   ddot       Scalar Float           annual rate of change of decl.       */
/*                                     (deg/yr)                             */
/*                                                                          */
/*   alt        Scalar Float           altitude above WGS84 Ellipsoid       */
/*                                                                          */
/*   epoch      Float array of MAXMOD  epoch of model.                      */
/*                                                                          */
/*   ext        Scalar Float           Three 1st-degree external coeff.     */
/*                                                                          */
/*   latitude   Scalar Float           Latitude.                            */
/*                                                                          */
/*   longitude  Scalar Float           Longitude.                           */
/*                                                                          */
/*   gh1        Float array            Schmidt quasi-normal internal        */
/*                                     spherical harmonic coeff.            */
/*                                                                          */
/*   gh2        Float array            Schmidt quasi-normal internal        */
/*                                     spherical harmonic coeff.            */
/*                                                                          */
/*   gha        Float array            Coefficients of resulting model.     */
/*                                                                          */
/*   ghb        Float array            Coefficients of rate of change model.*/
/*                                                                          */
/*   i          Scalar Float           Inclination (deg).                   */
/*                                                                          */
/*   idot       Scalar Float           Rate of change of i (deg/yr).        */
/*                                                                          */
/*   igdgc      Integer                Flag for geodetic or geocentric      */
/*                                     coordinate choice.                   */
/*                                                                          */
/*   inbuff     Char a of MAXINBUF     Input buffer.                        */
/*                                                                          */
/*   irec_pos   Integer array of MAXMOD Record counter for header           */
/*                                                                          */
/*   stream  Integer                   File handles for an opened file.     */
/*                                                                          */
/*   fileline   Integer                Current line in file (for errors)    */
/*                                                                          */
/*   max1       Integer array of MAXMOD Main field coefficient.             */
/*                                                                          */
/*   max2       Integer array of MAXMOD Secular variation coefficient.      */
/*                                                                          */
/*   max3       Integer array of MAXMOD Acceleration coefficient.           */
/*                                                                          */
/*   mdfile     Character array of PATH  Model file name.                   */
/*                                                                          */
/*   minyr      Float                  Min year of all models               */
/*                                                                          */
/*   maxyr      Float                  Max year of all models               */
/*                                                                          */
/*   yrmax      Float array of MAXMOD  Max year of model.                   */
/*                                                                          */
/*   yrmin      Float array of MAXMOD  Min year of model.                   */
/*                                                                          */
/****************************************************************************/


int main(int argv, char**argc)
{
#ifdef MAC
  ccommand(argv, argc);
#endif
  /*  Variable declaration  */
  
  /* Control variables */
  int   again = 1;
  int   decyears = 3;
  int   units = 4;
  int   decdeg = 3;
  int   range = -1;
  int   counter = 0;
  int   warn_H, warn_H_strong, warn_P;
  
  
  int   modelI;             /* Which model (Index) */
  int   nmodel;             /* Number of models in file */
  int   max1[MAXMOD];
  int   max2[MAXMOD];
  int   max3[MAXMOD];
  int   nmax;
  int   igdgc=3;
  int   isyear=-1;
  int   ismonth=-1;
  int   isday=-1;
  int   ieyear=-1;
  int   iemonth=-1;
  int   ieday=-1;
  int   ilat_deg=200;
  int   ilat_min=200;
  int   ilat_sec=200;
  int   ilon_deg=200;
  int   ilon_min=200;
  int   ilon_sec=200;
  int   fileline;
  int   ddeg, ideg;
  long  irec_pos[MAXMOD];
  
  char  mdfile[PATH];
  char  inbuff[MAXINBUFF];
  char  model[MAXMOD][9];
  char *begin;
  char *rest;
  
  float epoch[MAXMOD];
  float yrmin[MAXMOD];
  float yrmax[MAXMOD];
  float minyr;
  float maxyr;
  float altmin[MAXMOD];
  float altmax[MAXMOD];
  float minalt;
  float maxalt;
  float alt=-999999;
  float sdate=-1;
  float step=-1;
  float syr;
  float edate=-1;
  float latitude=200;
  float longitude=200;
  float ddot;
  float fdot;
  float hdot;
  float idot;
  float xdot;
  float ydot;
  float zdot;
  float dmin, imin;
  float warn_H_val, warn_H_strong_val;
  
  /*  Subroutines used  */
  
  void print_dashed_line();
  void print_long_dashed_line(void);
  void print_header();
  void print_result(double date, double d, double i, double h, double x, double y, double z, double f);
  void print_header_sv();
  void print_result_sv(double date, double ddot, double idot, double hdot, double xdot, double ydot, double zdot, double fdot);
  float degrees_to_decimal();
  float julday();
  int   interpsh();
  int   extrapsh();
  int   shval3();
  int   dihf();
  int   safegets(char *buffer,int n);
  int getshc();




  /* Initializations. */
  
  inbuff[MAXREAD+1]='\0';  /* Just to protect mem. */
  inbuff[MAXINBUFF-1]='\0';  /* Just to protect mem. */
  
  //printing out version number and header
  printf("\n\n Geomag v6.0 - Aug 24, 2005 ");
  
  if((argv==2)&&(*(argc[1])=='h')){
    printf("\n\nUSAGE: geomag model_file date coord alt lat lon\n\n");
    printf("Formats: \n");
    printf("   Date: xxxx.xxx for decimal  (1947.32)\n");
    printf("         YYYY,MM,DD for year, month, day  (1947.10.13)\n");
    printf("         or start_date-end_date (1943.21-1966.11)\n");
    printf("         or start_date-end_date-step (1943.21-1966.11-1.2)\n");
    printf("   Coord: D - Geodetic (WGS84 latitude and altitude AMSL)\n");
    printf("          C - Geocentric (spherical, altitude relative to Earth's center)\n");
    printf("   Altitude: Kxxxxxx.xxx for kilometers  (K1000.13)\n");
    printf("             Mxxxxxx.xxx for meters  (m1389.24)\n");
    printf("             Fxxxxxx.xxx for feet  (F192133.73)\n");
    printf("   Lat/Lon: xxx.xxx for decimal  (-76.53)\n");
    printf("            ddd,mm,ss for degrees, minutes, seconds (-60.23.22)\n");
    printf("            (Lat and Lon must be specified in the same format.)\n");
    printf("\nRange: \n");
    printf("   Date and altitude must fit model.\n");
    printf("   Lat: -90 to 90 (Use - to denote Southern latitude.)\n");
    printf("   Lon: -180 to 180 (Use - to denote Westen longitude.)\n");
    printf("   Minutes and Seconds: -59 to 59.\n");
    printf("\nNote: \n");
    printf("   All arguments are optional but must preserve order.\n");
    printf("   Invalid arguments are assumed to be 0.\n");
    exit(2);
  }
  
  if((argv>1)&&(argv<7)){
    printf("\n\nUSAGE (cammand line): geomag model_file date coord alt lat lon\n");
    printf("       or geomag h for help.\n");
    printf("USAGE (interactive): geomag\n\n");
  }
  
  /* Switch on how many arguments are supplied. */
  switch(argv){
  case 7 : strncpy(inbuff, argc[6], MAXREAD);
    if((rest=strchr(inbuff, ','))){     /* If it contains a comma */
      decdeg=2;                        /* Then not decimal degrees */
      begin=inbuff;
      rest[0]='\0';                    /* Chop into sub string */
      rest++;                          /* Move to next substring */
      ilon_deg=atoi(begin);
      begin=rest;
      if((rest=strchr(begin, ','))){
        rest[0]='\0';
        rest++;
        ilon_min=atoi(begin);
        ilon_sec=atoi(rest);
      } else {
        ilon_min=0;
        ilon_sec=0;
      }
    } else {
      decdeg=1;                        /* Else it's decimal */
      longitude=atof(argc[6]);
    }
    
  case 6 : strncpy(inbuff, argc[5], MAXREAD);
    if((rest=strchr(inbuff, ','))){
      decdeg=2;
      begin=inbuff;
      rest[0]='\0';
      rest++;
      ilat_deg=atoi(begin);
      begin=rest;
      if((rest=strchr(begin, ','))){
        rest[0]='\0';
        rest++;
        ilat_min=atoi(begin);
        ilat_sec=atoi(rest);
      } else {
        ilat_min=0;
        ilat_sec=0;
      }
    } else {
      decdeg=1;
      latitude=atof(argc[5]);
    }
    
  case 5 : strncpy(inbuff, argc[4], MAXREAD);
    inbuff[0]=toupper(inbuff[0]);
    if(inbuff[0]=='K') units=1;
    else if(inbuff[0]=='M') units=2;
    else if(inbuff[0]=='F') units=3;
    if(strlen(inbuff)>1){
      inbuff[0]='\0';
      begin=inbuff+1;
      alt=atof(begin);
    }
    
  case 4 : strncpy(inbuff, argc[3], MAXREAD);
    inbuff[0]=toupper(inbuff[0]);
    if(inbuff[0]=='D') igdgc=1;
    else if (inbuff[0]=='C') igdgc=2;
    
  case 3 : strncpy(inbuff, argc[2], MAXREAD);
    if((rest=strchr(inbuff, '-'))){   /* If it contains a dash */
      range = 2;                     /* They want a range */
      rest[0]='\0';                  /* Sep dates */
      rest++;
      begin=rest;
      if((rest=strchr(begin, '-'))){    /* If it contains 2 dashs */
        rest[0]='\0';                  /* Sep step */
        rest++;
        step=atof(rest);               /* Get step size */
      }
      if((rest=strchr(begin, ','))){    /* If it contains a comma */
        decyears=2;                    /* It's not decimal years */
        rest[0]='\0';
        rest++;
        ieyear=atoi(begin);
        begin=rest;
        if((rest=strchr(begin, ','))){
          rest[0]='\0';
          rest++;
          iemonth=atoi(begin);
          ieday=atoi(rest);
        } else {
          iemonth=0;
          ieday=0;
        }
        if((rest=strchr(inbuff, ','))){
          begin=inbuff;
          rest[0]='\0';
          rest++;
          isyear=atoi(begin);
          begin=rest;
          if((rest=strchr(begin, ','))){
            rest[0]='\0';
            rest++;
            ismonth=atoi(begin);
            isday=atoi(rest);
          } else {
            ismonth=0;
            isday=0;
          }
        } else {
          sdate=atof(inbuff);
        }
      } else {
        decyears=1;                    /* Else it's decimal years */
        sdate=atof(inbuff);
        edate=atof(begin);
      }
    } else {
      range = 1;
      if((rest=strchr(inbuff, ','))){   /* If it contains a comma */
        decyears=2;                    /* It's not decimal years */
        begin=inbuff;
        rest[0]='\0';
        rest++;
        isyear=atoi(begin);
        begin=rest;
        if((rest=strchr(begin, ','))){
          rest[0]='\0';
          rest++;
          ismonth=atoi(begin);
          isday=atoi(rest);
        } else {
          ismonth=0;
          isday=0;
        }
      } else {
        decyears=1;                    /* Else it's decimal years */
        sdate=atof(argc[2]);
      }
    }
    if(sdate==0){                        /* If date not valid */
      decyears=-1;
      range=-1;
    }
    
  case 2 : strncpy(mdfile,argc[1],MAXREAD);
    stream=fopen(mdfile, "rt");
    break;
  }
  
  
  /*  Obtain the desired model file and read the data  */
  
  printf("\n");
  while(again==1){
    again++;
    warn_H = 0;
    warn_H_val = 99999.0;
    warn_H_strong = 0;
    warn_H_strong_val = 99999.0;
    warn_P = 0;

    printf("\n");
    
    while(stream == NULL){
      printf("\n\n");
      printf("What is the name of the model data file to be opened? ==> ");
      safegets(inbuff,MAXREAD);
      strcpy(mdfile, inbuff);
      if(!(stream = fopen(mdfile, "rt"))){
        printf("\nError opening file %s.", mdfile);
      }
    }
    rewind(stream);
    
    fileline = 0;                            /* First line will be 1 */
    modelI = -1;                             /* First model will be 0 */
    while(fgets(inbuff,MAXREAD,stream)){     /* While not end of file
                                              * read to end of line or buffer */
      fileline++;                           /* On new line */
      
      
      if(strlen(inbuff) != RECL){       /* IF incorrect record size */
        printf("Corrupt record in file %s on line %d.\n", mdfile, fileline);
        fclose(stream);
        exit(5);
      }
      
      /* old statement Dec 1999 */
      /*       if(!strncmp(inbuff,"    ",4)){         /* If 1st 4 chars are spaces */
      /* New statement Dec 1999 changed by wmd  required by year 2000 models */
      if(!strncmp(inbuff,"   ",3)){         /* If 1st 3 chars are spaces */
        modelI++;                           /* New model */
        
        if(modelI > MAXMOD){                /* If too many headers */
          printf("Too many models in file %s on line %d.", mdfile, fileline);
          fclose(stream);
          exit(6);
        }
        
        irec_pos[modelI]=ftell(stream);
        /* Get fields from buffer into individual vars.  */
        sscanf(inbuff, "%s%f%d%d%d%f%f%f%f", model[modelI], &epoch[modelI],
               &max1[modelI], &max2[modelI], &max3[modelI], &yrmin[modelI],
               &yrmax[modelI], &altmin[modelI], &altmax[modelI]);
        
        /* Compute date range for all models */
        if(modelI == 0){                    /*If first model */
          minyr=yrmin[0];
          maxyr=yrmax[0];
        } else {
          if(yrmin[modelI]<minyr){
            minyr=yrmin[modelI];
          }
          if(yrmax[modelI]>maxyr){
            maxyr=yrmax[modelI];
          }
        }
        
      }
      
    }
    nmodel = modelI + 1;
    fclose(stream);
    
    /* if date specified in command line then warn if past end of validity */
    
    if((sdate>maxyr)&&(sdate<maxyr+1)){
      printf("\nWarning: The date %4.2f is out of range,\n", sdate);
      printf("         but still within one year of model expiration date.\n");
      printf("         An updated model file is available before 1.1.%4.0f\n",maxyr);
    }
    
    
    /*  Take in field data  */
    
    /* Get date */
    
    while(again==2){
      again++;
      while((decyears!=1)&&(decyears!=2)){
        printf("\nHow would you like to enter the date?\n");
        printf("       1) In decimal years.\n");
        printf("       2) In year, month, and day.\n");
        printf("\n                            ==> ");
        safegets(inbuff, MAXREAD);
        decyears=atoi(inbuff);
      }
      
      while((range!=1)&&(range!=2)){
        printf("\nWould you like output for a single date or for a range of dates?\n");
        printf("       1) A single date.\n");
        printf("       2) A range of dates.\n");
        printf("\n                            ==> ");
        safegets(inbuff, MAXREAD);
        range=atoi(inbuff);
      }
      
      
      if(range == 1){
        
        while((sdate<minyr)||(sdate>maxyr+1)){
          
          if(decyears==1){
            printf("\nEnter the decimal date (%4.2f to %4.0f): ",minyr, maxyr);
            safegets(inbuff, MAXREAD);
            sdate=atof(inbuff);
          } else {
            
            while((isyear>(int)maxyr+1)||(isyear<(int)minyr)){
              printf("\nEnter the date (%4.2f to %4.2f)\n ", minyr, maxyr);
              printf("\n   Year (%d to %d): ",(int)minyr,(int)maxyr);
              safegets(inbuff, MAXREAD);
              isyear=atoi(inbuff);
            }
            
            while((ismonth>12)||(ismonth<1)){
              printf("\n   Month (1-12): ");
              safegets(inbuff, MAXREAD);
              ismonth=atoi(inbuff);
            }
            
            while((isday>31)||(isday<1)){
              printf("\n   Day (1-31): ");
              safegets(inbuff, MAXREAD);
              isday=atoi(inbuff);
            }
            
            sdate = julday(ismonth,isday,isyear);
          }                  
          if((sdate<minyr)||(sdate>=maxyr+1)){
            ismonth=isday=isyear=0;
            printf("\nError: The date %4.2f is out of range.\n", sdate);
          }
          
          if((sdate>maxyr)&&(sdate<maxyr+1)){
            printf("\nWarning: The date %4.2f is out of range,\n", sdate);
            printf("         but still within one year of model expiration date.\n");
            printf("         An updated model file is available before 1.1.%4.0f\n",maxyr);
          }
          
        } /* if single date */
      } else {
        
        while((sdate<minyr)||(sdate>maxyr)){
          
          if(decyears==1){
            printf("\nEnter the decimal start date (%4.2f to %4.0f): ",minyr, maxyr);
            safegets(inbuff, MAXREAD);
            sdate=atof(inbuff);
          } else {
            
            while((isyear>(int)maxyr)||(isyear<(int)minyr)){
              ismonth=isday=isyear=0;
              printf("\nEnter the start date (%4.2f to %4.2f)\n ", minyr, maxyr);
              printf("\n   Year (%d to %d): ",(int)minyr,(int)(maxyr));
              safegets(inbuff, MAXREAD);
              isyear=atoi(inbuff);
            }
            
            while((ismonth>12)||(ismonth<1)){
              printf("\n   Month (1-12): ");
              safegets(inbuff, MAXREAD);
              ismonth=atoi(inbuff);
            }
            
            while((isday>31)||(isday<1)){
              printf("\n   Day (1-31): ");
              safegets(inbuff, MAXREAD);
              isday=atoi(inbuff);
            }
            
            sdate = julday(ismonth,isday,isyear);
            
            if((sdate<minyr)||(sdate>maxyr)){
              printf("\nThe start date %4.2f is out of range.\n", sdate);
            }
          }
        }
        
        while((edate<=sdate)||(edate>maxyr+1)){
          
          if(decyears==1){
            printf("\nEnter the decimal end date (%4.2f to %4.0f): ",sdate, maxyr);
            safegets(inbuff, MAXREAD);
            edate=atof(inbuff);
          } else {
            
            while((ieyear>(int)maxyr)||(ieyear<(int)sdate)){
              iemonth=ieday=ieyear=0;
              printf("\nEnter the end date (%4.2f to %4.0f)\n ", sdate, maxyr);
              printf("\n   Year (%d to %d): ",(int)sdate,(int)(maxyr));
              safegets(inbuff, MAXREAD);
              ieyear=atoi(inbuff);
            }
            
            while((iemonth>12)||(iemonth<1)){
              printf("\n   Month (1-12): ");
              safegets(inbuff, MAXREAD);
              iemonth=atoi(inbuff);
            }
            
            while((ieday>31)||(ieday<1)){
              printf("\n   Day (1-31): ");
              safegets(inbuff, MAXREAD);
              ieday=atoi(inbuff);
            }
            
            edate = julday(iemonth,ieday,ieyear);
            
          }
          
          if((edate<sdate)||(edate>maxyr+1)){
            printf("\nThe date %4.2f is out of range.\n", edate);
          }
          
          if(edate>maxyr && edate<=maxyr+1){
            printf("\nWarning: The end date %4.2f is out of range,\n", edate);
            printf("         but still within one year of model expiration date.\n");
            printf("         An updated model file is available before 1.1.%4.0f\n",maxyr);
          }
        }
        
        while((step<=0)||(step>(edate-sdate))){
          printf("\nEnter the step size in years. (0 to %4.2f): ",edate-sdate);
          safegets(inbuff, MAXREAD);
          step=atof(inbuff);
        }
        
      }
      
      
      /* Pick model */
      for (modelI=0; modelI<nmodel; modelI++)
        if (sdate<yrmax[modelI]) break;
      if (modelI == nmodel) modelI--;           /* if beyond end of last model use last model */
      
      /* Get altitude min and max for selected model. */
      minalt=altmin[modelI];
      maxalt=altmax[modelI];
      
      /* Get Coordinate prefs */
      while((igdgc!=1)&&(igdgc!=2)){
        printf("\n\nEnter Coordinate Preferences:");
        printf("\n    1) Geodetic (WGS84 latitude and altitude AMSL)");
        printf("\n    2) Geocentric (spherical, altitude relative to Earth's center)\n");
        printf("\n                            ==> ");
        safegets(inbuff, MAXREAD);
        igdgc=atoi(inbuff);
      }
      
      /* If needed modify ranges to reflect coords. */
      if(igdgc==2){
        minalt+=6371.2;  /* Add radius to ranges. */
        maxalt+=6371.2;
      }
      
      /* Get unit prefs */
      if(igdgc==1){
        while((units>3)||(units<1)){
          printf("\n\nEnter Unit Preferences:");
          printf("\n       1) Kilometers");
          printf("\n       2) Meters");
          printf("\n       3) Feet\n");
          printf("\n                            ==> ");
          safegets(inbuff, MAXREAD);
          units=atoi(inbuff);
        }
      }
      else units = 1; /* geocentric always in km */
      
      /* Do unit conversions if neccessary */
      if(units==2){
        minalt*=1000.0;
        maxalt*=1000.0;
      } else if(units==3){
        minalt*=3280.0839895;
        maxalt*=3280.0839895;
      }
      
      /* Get altitude */
      while((alt<minalt)||(alt>maxalt)){
        if (igdgc==2) printf("\n\nEnter geocentric altitude in km (%.2f to %.2f): ", minalt, maxalt);
        if (igdgc==1 && units==1) printf("\n\nEnter geodetic altitude AMSL (WGS84) in km (%.2f to %.2f): ", minalt, maxalt);
        if (igdgc==1 && units==2) printf("\n\nEnter geodetic altitude AMSL (WGS84) in meters (%.2f to %.2f): ", minalt, maxalt);
        if (igdgc==1 && units==3) printf("\n\nEnter geodetic altitude AMSL (WGS84) in feet (%.2f to %.2f): ", minalt, maxalt);
        safegets(inbuff, MAXREAD);
        alt=atof(inbuff);
      }
      
      /* Convert altitude to km */
      if(units==2){
        alt *= 0.001;
      } else if(units==3){
        alt /= 3280.0839895;
      }
      
      /* Get lat/long prefs */
      while((decdeg!=1)&&(decdeg!=2)){
        printf("\n\nHow would you like to enter the latitude and longitude?:");
        printf("\n       1) In decimal degrees.");
        printf("\n       2) In degrees, minutes, and seconds.\n");
        printf("\n                            ==> ");
        safegets(inbuff, MAXREAD);
        decdeg=atoi(inbuff);
      }
      
      
      /* Get lat/lon */
      
      if(decdeg==1){
        while((latitude<-90)||(latitude>90)){
          printf("\n\nEnter the decimal latitude (-90 to 90) (- for Southern hemisphere).\n");
          safegets(inbuff, MAXREAD);
          latitude=atof(inbuff);
        }
        while((longitude<-180)||(longitude>180)){
          printf("\n\nEnter the decimal longitude (-180 to 180) (- for Western hemisphere).\n");
          safegets(inbuff, MAXREAD);
          longitude=atof(inbuff);
        }
      } else {
        latitude=degrees_to_decimal(ilat_deg,ilat_min,ilat_sec);
        longitude=degrees_to_decimal(ilon_deg,ilon_min,ilon_sec);
        while((latitude<-90)||(latitude>90)){
          ilat_deg=ilat_min=ilat_sec=200;
          printf("\n\nEnter the decimal latitude (-90 to 90) (- for Southern hemisphere).\n");
          while((ilat_deg<-90)||(ilat_deg>90)){
            printf("\nDegrees (-90 to 90): ");
            safegets(inbuff, MAXREAD);
            ilat_deg=atoi(inbuff);
          }
          while((ilat_min<-59)||(ilat_min>59)){
            printf("\nMinutes (-59 to 59): ");
            safegets(inbuff, MAXREAD);
            ilat_min=atoi(inbuff);
          }
          while((ilat_sec<-59)||(ilat_sec>59)){
            printf("\nSeconds (-59 to 59): ");
            safegets(inbuff, MAXREAD);
            ilat_sec=atoi(inbuff);
          }
          
          latitude=degrees_to_decimal(ilat_deg,ilat_min,ilat_sec);
          
          if((latitude<-90)||(latitude>90)){
            printf("\nThe latitude %3.2f is out of range", latitude);
          }
          
        }
        while((longitude<-180)||(longitude>180)){
          ilon_deg=ilon_min=ilon_sec=200;
          printf("\n\nEnter the decimal longitude (-180 to 180) (- for Western hemisphere).\n");
          while((ilon_deg<-180)||(ilon_deg>180)){
            printf("\nDegrees (-180 to 180): ");
            safegets(inbuff, MAXREAD);
            ilon_deg=atoi(inbuff);
          }
          while((ilon_min<-59)||(ilon_min>59)){
            printf("\nMinutes (0 to 59): ");
            safegets(inbuff, MAXREAD);
            ilon_min=atoi(inbuff);
          }
          while((ilon_sec<-59)||(ilon_sec>59)){
            printf("\nSeconds (0 to 59): ");
            safegets(inbuff, MAXREAD);
            ilon_sec=atoi(inbuff);
          }
          
          longitude=degrees_to_decimal(ilon_deg,ilon_min,ilon_sec);
          
          if((longitude<-180)||(longitude>180)){
            printf("\nThe longitude %3.2f is out of range", longitude);
          }
          
        }
      }
      
      
      
      /** This will compute everything needed for 1 point in time. **/


      if(max2[modelI] == 0) {
        getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1);
        getshc(mdfile, 1, irec_pos[modelI+1], max1[modelI+1], 2);
        nmax = interpsh(sdate, yrmin[modelI], max1[modelI],
                        yrmin[modelI+1], max1[modelI+1], 3);
        nmax = interpsh(sdate+1, yrmin[modelI] , max1[modelI],
                        yrmin[modelI+1], max1[modelI+1],4);
      } else {
        getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1);
        getshc(mdfile, 0, irec_pos[modelI], max2[modelI], 2);
        nmax = extrapsh(sdate, epoch[modelI], max1[modelI], max2[modelI], 3);
        nmax = extrapsh(sdate+1, epoch[modelI], max1[modelI], max2[modelI], 4);
      }

      
      /* Do the first calculations */
      shval3(igdgc, latitude, longitude, alt, nmax, 3,
             IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
      dihf(3);
      shval3(igdgc, latitude, longitude, alt, nmax, 4,
             IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
      dihf(4);
      

      ddot = ((dtemp - d)*57.29578);
      if (ddot > 180.0) ddot -= 360.0;
      if (ddot <= -180.0) ddot += 360.0;
      ddot *= 60.0;
      
      idot = ((itemp - i)*57.29578)*60;
      d = d*(57.29578);   i = i*(57.29578);
      hdot = htemp - h;   xdot = xtemp - x;
      ydot = ytemp - y;   zdot = ztemp - z;
      fdot = ftemp - f;

      /* deal with geographic and magnetic poles */
      
      if (h < 100.0) /* at magnetic poles */
        {
          d = NaN;
          ddot = NaN;
          /* while rest is ok */
        }

      if (h < 1000.0) 
        {
          warn_H = 0;
          warn_H_strong = 1;
          if (h<warn_H_strong_val) warn_H_strong_val = h;
        }
      else if (h < 5000.0 && !warn_H_strong) 
        {
          warn_H = 1;
          if (h<warn_H_val) warn_H_val = h;
        }
                  
      if (90.0-fabs(latitude) <= 0.001) /* at geographic poles */
        {
          x = NaN;
          y = NaN;
          d = NaN;
          xdot = NaN;
          ydot = NaN;
          ddot = NaN;
          warn_P = 1;
          warn_H = 0;
          warn_H_strong = 0;
          /* while rest is ok */
        }
      
      /** Above will compute everything for 1 point in time.  **/
      
      
      /*  Output the final results. */
      
      printf("\n\n\n  Model: %s \n", model[modelI]);
      if(decdeg==1)
        {
          printf("  Latitude: %4.2f deg\n", latitude);
          printf("  Longitude: %4.2f deg\n", longitude);
        } 
      else 
        {
          printf("  Latitude: %d deg, %d min, %d sec\n",
                 ilat_deg,ilat_min, ilat_sec);
          printf("  Longitude: %d deg,  %d min, %d sec\n",
                 ilon_deg, ilon_min, ilon_sec);
        }
      printf("  Altitude: ");
      if(units==1)
          printf("%.2f km\n", alt);
      else if(units==2)
        printf("%.2f meters\n", alt*1000.0);
      else 
        printf("%.2f ft\n", (alt*3280.0839895));

      if(range==1)
        {
          printf("  Date of Interest: ");
          if(decyears==1)
            printf(" %4.2f\n\n", sdate);
          else 
            printf("%d/%d/%d\n\n", ismonth, isday, isyear);
          
          print_header();
          print_result(sdate,d, i, h, x, y, z, f);
          print_long_dashed_line();
          print_header_sv();
          print_result_sv(sdate,ddot,idot,hdot,xdot,ydot,zdot,fdot);
          print_dashed_line();

        } /* if range == 1 */
      else 
        {
          printf("  Range of Interest: ");
          if(decyears==1)
            printf("%4.2f - %4.2f, step %4.2f\n\n",sdate, edate, step);
          else
            printf("%d/%d/%d - %d/%d/%d, step %4.2f\n\n",ismonth, isday, isyear, iemonth, ieday, ieyear, step);

          print_header();
          print_result(sdate,d, i, h, x, y, z, f);

          for(syr=sdate+step;syr<(edate+step);syr+=step)
            {
              if((syr>edate)&&(edate!=(syr-step)))
                {
                  syr=edate;
                  print_long_dashed_line();
                }

              /* Do the calculations */

              for (counter=0;counter<step;counter++)
                {
                  if(max2[modelI] == 0){       /*If not last element in array */
                    if(syr>yrmin[modelI+1]){  /* And past model boundary */
                      modelI++;              /* Get next model */
                    }
                  }
                } /* for counter */

              if(max2[modelI] == 0)       /*If still not last element */
                {
                  getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1);
                  getshc(mdfile, 1, irec_pos[modelI+1], max1[modelI+1], 2);
                  nmax = interpsh(syr, yrmin[modelI], max1[modelI],
                                  yrmin[modelI+1], max1[modelI+1], 3);
                  nmax = interpsh(syr+1, yrmin[modelI] , max1[modelI],
                                  yrmin[modelI+1], max1[modelI+1],4);
                } 
              else 
                {
                  getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1);
                  getshc(mdfile, 0, irec_pos[modelI], max2[modelI], 2);
                  nmax = extrapsh(syr, epoch[modelI], max1[modelI],
                                  max2[modelI], 3);
                  nmax = extrapsh(syr+1, epoch[modelI], max1[modelI],
                                  max2[modelI], 4);
                }
              shval3(igdgc, latitude, longitude, alt, nmax, 3,
                     IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
              dihf(3);
              shval3(igdgc, latitude, longitude, alt, nmax, 4,
                     IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
              dihf(4);
          
              ddot = ((dtemp - d)*57.29578);
              if (ddot > 180.0) ddot -= 360.0;
              if (ddot <= -180.0) ddot += 360.0;
              ddot *= 60.0;
              
              idot = ((itemp - i)*57.29578)*60.;
              d = d*(57.29578);   i = i*(57.29578);
              hdot = htemp - h;   xdot = xtemp - x;
              ydot = ytemp - y;   zdot = ztemp - z;
              fdot = ftemp - f;

              /* deal with geographic and magnetic poles */
                    
              if (h < 100.0) /* at magnetic poles */
                {
                  d = NaN;
                  ddot = NaN;
                  /* while rest is ok */
                }

              if (90.0-fabs(latitude) <= 0.001) /* at geographic poles */
                {
                  x = NaN;
                  y = NaN;
                  d = NaN;
                  xdot = NaN;
                  ydot = NaN;
                  ddot = NaN;
                  warn_P = 1;
                  warn_H = 0;
                  warn_H_strong = 0;
                  /* while rest is ok */
                }

              print_result(syr, d, i, h, x, y, z, f);
            } /* for syr */
          
          print_long_dashed_line();
          print_header_sv();
          print_result_sv(edate,ddot,idot,hdot,xdot,ydot,zdot,fdot);
          print_dashed_line();
        } /* if range > 1 */

      if (warn_H)
        {
          printf("\nWarning: The horizontal field strength at this location is only %6.1f nT\n",warn_H_val);
          printf("         Compass readings have large uncertainties in areas where H is\n");
          printf("         smaller than 5000 nT\n\n");
        } 
      if (warn_H_strong)
        {
          printf("\nWarning: The horizontal field strength at this location is only %6.1f nT\n",warn_H_strong_val);
          printf("         Compass readings have VERY LARGE uncertainties in areas where H is\n");
          printf("         smaller than 1000 nT\n\n");
        }
      if (warn_P)
        {
          printf("\nWarning: Location is at geographic pole where X, Y, and declination are undefined\n\n");
        } 


      while((again!=0)&&(again!=1)&&(again!=2)&&(argv<7)){
        printf("\n  Enter");
        printf("\n     0) to quit.");
        printf("\n     1) to select a new model input file.");
        printf("\n     2) to compute for a new point using same data file.");
        printf("\n\n   ==> ");
        safegets(inbuff, MAXREAD);
        again=atoi(inbuff);
      }
      if(again!=0){
        /* Reset defaults to catch on all while loops */
        igdgc=decyears=units=decdeg=-1;
        ismonth=isday=isyear=sdate=edate=range=step=-1;
        latitude=ilat_deg=ilat_min=ilat_sec=200;
        longitude=ilon_deg=ilon_min=ilon_sec=200;
        alt=-9999999;
      }
    }
    /* Must specify new stream */
    stream=NULL;
  }
  if(argv<7){
    printf("\nThe same result could have been generated with the following command:\n");
    printf("%s %s ", argc[0], mdfile);
    if(range == 1){
      if(decyears==1){
        printf("%4.2f ", sdate);
      } else {
        printf("%d,%d,%d ", isyear,ismonth,isday);
      }
    } else {
      if(decyears==1){
        printf("%4.2f-%4.2f-%4.2f ", sdate,edate,step);
      } else {
        printf("%d,%d,%d-%d,%d,%d-%4.2f ", isyear,ismonth,isday,ieyear,iemonth,ieday,step);
      }
    }
    if(igdgc==1){
      printf("D ");
    } else {
      printf("C ");
    }
    if(units==1){
      printf("K%.2f ", alt);
    } else if(units==2){
      printf("M%.2f ", alt*1000.0);
    } else {
      printf("F%.2f ", (alt*3280.0839895));
    }
    if(decdeg==1){
      printf("%3.2f %4.2f\n\n", latitude, longitude);
  } else {
    printf("%d,%d,%d %d,%d,%d\n\n", ilat_deg, ilat_min, ilat_sec,
           ilon_deg, ilon_min, ilon_sec);
  }
    
  }
  return 0;
}

void print_dashed_line(void)
{
  printf(" -----------------------------------------------------------------------------\n");
  return;
}


void print_long_dashed_line(void)
{
  printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  return;
}
      
void print_header(void)
{ 
  print_dashed_line();
  printf("   Date         D          I           H        X        Y        Z        F\n");
  printf("   (yr)     (deg min)  (deg min)     (nT)     (nT)     (nT)     (nT)     (nT)\n");
  return;
}

void print_result(double date, double d, double i, double h, double x, double y, double z, double f)
{
  int   ddeg,ideg;
  float dmin,imin;

      /* Change d and i to deg and min */
      
  ddeg=(int)d;
  dmin=(d-(float)ddeg)*60;
  if(ddeg!=0) dmin=fabs(dmin);
  ideg=(int)i;
  imin=(i-(float)ideg)*60;
  if(ideg!=0) imin=fabs(imin);


  if (isnan(d))
    {
      if (isnan(x))
        printf("  %4.2f       NaN   %4dd %2.0fm  %8.1f      NaN      NaN %8.1f %8.1f\n",date,ideg,imin,h,z,f);
      else
        printf("  %4.2f       NaN   %4dd %2.0fm  %8.1f %8.1f %8.1f %8.1f %8.1f\n",date,ideg,imin,h,x,y,z,f);
    }
  else 
    printf("  %4.2f  %4dd %2.0fm  %4dd %2.0fm  %8.1f %8.1f %8.1f %8.1f %8.1f\n",date,ddeg,dmin,ideg,imin,h,x,y,z,f);
  return;
} /* print_result */

void print_header_sv(void)
{
  printf("   Date        dD         dI           dH       dX       dY       dZ       dF\n");
  printf("   (yr)     (min/yr)   (min/yr)    (nT/yr)  (nT/yr)  (nT/yr)  (nT/yr)  (nT/yr)\n");
} /* print_header_sv */

void print_result_sv(double date, double ddot, double idot, double hdot, double xdot, double ydot, double zdot, double fdot)
{
  if (isnan(ddot))
    {
      if (isnan(xdot))
        printf("  %4.2f       NaN  %7.1f     %8.1f      NaN      NaN %8.1f %8.1f\n",date,idot,hdot,zdot,fdot);
      else
        printf("  %4.2f       NaN  %7.1f     %8.1f %8.1f %8.1f %8.1f %8.1f\n",date,idot,hdot,xdot,ydot,zdot,fdot);
    }
  else 
    printf("  %4.2f  %7.1f   %7.1f     %8.1f %8.1f %8.1f %8.1f %8.1f\n",date,ddot,idot,hdot,xdot,ydot,zdot,fdot);
  return;
} /* print_result_sv */



/****************************************************************************/
/*                                                                          */
/*                       Subroutine safegets                                */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*  Gets characters from stdin untill it has reached n characters or \n,    */
/*     whichever comes first.  \n is converted to \0.                       */
/*                                                                          */
/*  Input: n - Integer number of chars                                      */
/*         *buffer - Character array ptr which can contain n+1 characters   */
/*                                                                          */
/*  Output: size - integer size of sting in buffer                          */
/*                                                                          */
/*  Note: All strings will be null terminated.                              */
/*                                                                          */
/*  By: David Owens                                                         */
/*      dio@ngdc.noaa.gov                                                   */
/****************************************************************************/

int safegets(char *buffer,int n){
   char *ptr;                    /** ptr used for finding '\n' **/

   fgets(buffer,n,stdin);        /** Get n chars **/
   buffer[n+1]='\0';             /** Set last char to null **/
   ptr=strchr(buffer,'\n');      /** If string contains '\n' **/
   if(ptr!=NULL){                /** If string contains '\n' **/
      ptr[0]='\0';               /** Change char to '\0' **/
   }

   return strlen(buffer);        /** Return the length **/
}


/****************************************************************************/
/*                                                                          */
/*                       Subroutine degrees_to_decimal                      */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Converts degrees,minutes, seconds to decimal degrees.                */
/*                                                                          */
/*     Input:                                                               */
/*            degrees - Integer degrees                                     */
/*            minutes - Integer minutes                                     */
/*            seconds - Integer seconds                                     */
/*                                                                          */
/*     Output:                                                              */
/*            decimal - degrees in decimal degrees                          */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 12, 1988                                                */
/*                                                                          */
/****************************************************************************/

float degrees_to_decimal(int degrees,int minutes,int seconds)
{
   float deg;
   float min;
   float sec;
   float decimal;

   deg = degrees;
   min = minutes/60.0;
   sec = seconds/3600.0;

   decimal = fabs(sec) + fabs(min) + fabs(deg);

   if (deg < 0) {
      decimal = -decimal;
   } else if(deg == 0){
      if(min < 0){
         decimal = -decimal;
      } else if(min == 0){
         if(sec<0){
            decimal = -decimal;
         }
      }
   }

   return(decimal);
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine julday                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the decimal day of year from month, day, year.              */
/*     Leap years accounted for 1900 and 2000 are not leap years.           */
/*                                                                          */
/*     Input:                                                               */
/*           year - Integer year of interest                                */
/*           month - Integer month of interest                              */
/*           day - Integer day of interest                                  */
/*                                                                          */
/*     Output:                                                              */
/*           date - Julian date to thousandth of year                       */
/*                                                                          */
/*     FORTRAN                                                              */
/*           S. McLean                                                      */
/*           NGDC, NOAA egc1, 325 Broadway, Boulder CO.  80301              */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 12, 1988                                                */
/*                                                                          */
/*     Julday Bug Fix                                                       */
/*           Thanks to Rob Raper                                            */
/****************************************************************************/


float julday(i_month, i_day, i_year)
int i_month;
int i_day;
int i_year;
{
   int   aggregate_first_day_of_month[13];
   int   leap_year = 0;
   int   truncated_dividend;
   float year;
   float day;
   float decimal_date;
   float remainder = 0.0;
   float divisor = 4.0;
   float dividend;
   float left_over;

   aggregate_first_day_of_month[1] = 1;
   aggregate_first_day_of_month[2] = 32;
   aggregate_first_day_of_month[3] = 60;
   aggregate_first_day_of_month[4] = 91;
   aggregate_first_day_of_month[5] = 121;
   aggregate_first_day_of_month[6] = 152;
   aggregate_first_day_of_month[7] = 182;
   aggregate_first_day_of_month[8] = 213;
   aggregate_first_day_of_month[9] = 244;
   aggregate_first_day_of_month[10] = 274;
   aggregate_first_day_of_month[11] = 305;
   aggregate_first_day_of_month[12] = 335;

   /* Test for leap year.  If true add one to day. */

   year = i_year;                                 /*    Century Years not   */
   if ((i_year != 1900) && (i_year != 2100))      /*  divisible by 400 are  */
   {                                              /*      NOT leap years    */
      dividend = year/divisor;
      truncated_dividend = dividend;
      left_over = dividend - truncated_dividend;
      remainder = left_over*divisor;
      if ((remainder > 0.0) && (i_month > 2))
      {
	 leap_year = 1;
      }
      else
      {
	 leap_year = 0;
      }
   }
   day = aggregate_first_day_of_month[i_month] + i_day - 1 + leap_year;
   if (leap_year)
   {
      decimal_date = year + (day/366.0);  /*In version 3.0 this was incorrect*/
   }
   else
   {
      decimal_date = year + (day/365.0);  /*In version 3.0 this was incorrect*/
   }
   return(decimal_date);
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine getshc                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Reads spherical harmonic coefficients from the specified             */
/*     model into an array.                                                 */
/*                                                                          */
/*     Input:                                                               */
/*           stream     - Logical unit number                               */
/*           iflag      - Flag for SV equal to ) or not equal to 0          */
/*                        for designated read statements                    */
/*           strec      - Starting record number to read from model         */
/*           nmax_of_gh - Maximum degree and order of model                 */
/*                                                                          */
/*     Output:                                                              */
/*           gh1 or 2   - Schmidt quasi-normal internal spherical           */
/*                        harmonic coefficients                             */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Bill Flanagan                                                  */
/*           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 15, 1988                                                */
/*                                                                          */
/****************************************************************************/


int getshc(file, iflag, strec, nmax_of_gh, gh)
char file[PATH];
int       iflag;
long int  strec;
int       nmax_of_gh;
int       gh;
{
   char  inbuff[MAXINBUFF];
   char irat[9];
   int ii,m,n,mm,nn;
   int ios;
   int line_num;
   float g,hh;
   float trash;

   stream = fopen(file, "rt");
   if (stream == NULL)
   {
      printf("\nError on opening file %s", file);
   }
   else
   {
      ii = 0;
      ios = 0;
      fseek(stream,strec,SEEK_SET);
      for ( nn = 1; nn <= nmax_of_gh; ++nn)
      {
	 for (mm = 0; mm <= nn; ++mm)
	 {
	    if (iflag == 1)
	    {
               fgets(inbuff, 3, stream);
               inbuff[3]='\0';
               sscanf(inbuff, "%d", &m);
               fgets(inbuff, 3, stream);
               inbuff[3]='\0';
               sscanf(inbuff, "%d", &n);
               fgets(inbuff, MAXREAD-4, stream);
               sscanf(inbuff, "%f%f%f%f%s%d",
                      &g, &hh, &trash, &trash, irat, &line_num);
	    }
	    else
	    {
               fgets(inbuff, 3, stream);
               inbuff[3]='\0';
               sscanf(inbuff, "%d", &m);
               fgets(inbuff, 3, stream);
               inbuff[3]='\0';
               sscanf(inbuff, "%d", &n);
               fgets(inbuff, MAXREAD-4, stream);
               sscanf(inbuff, "%f%f%f%f%s%d",
                      &trash, &trash, &g, &hh, irat, &line_num);
	    }
	    if ((nn != n) || (mm != m))
	    {
	       ios = -2;
	       fclose(stream);
	       return(ios);
	    }
	    ii = ii + 1;
	    switch(gh)
	    {
		case 1:  gh1[ii] = g;
			 break;
		case 2:  gh2[ii] = g;
			 break;
		default: printf("\nError in subroutine getshc");
			 break;
	    }
	    if (m != 0)
	    {
	       ii = ii+ 1;
	       switch(gh)
	       {
		   case 1:  gh1[ii] = hh;
			    break;
		   case 2:  gh2[ii] = hh;
			    break;
		   default: printf("\nError in subroutine getshc");
			    break;
	       }
	    }
	 }
      }
   }
   fclose(stream);
   return(ios);
}


/****************************************************************************/
/*                                                                          */
/*                           Subroutine extrapsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Extrapolates linearly a spherical harmonic model with a              */
/*     rate-of-change model.                                                */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of base model                                  */
/*           nmax1    - maximum degree and order of base model              */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of base model                 */
/*           nmax2    - maximum degree and order of rate-of-change model    */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of rate-of-change model       */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - Schmidt quasi-normal internal spherical             */
/*                    harmonic coefficients                                 */
/*           nmax   - maximum degree and order of resulting model           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 16, 1988                                                */
/*                                                                          */
/****************************************************************************/


int extrapsh(date, dte1, nmax1, nmax2, gh)
float date;
float dte1;
int   nmax1;
int   nmax2;
int   gh;
{
   int   nmax;
   int   k, l;
   int   ii;
   float factor;

   factor = date - dte1;
   if (nmax1 == nmax2)
   {
      k =  nmax1 * (nmax1 + 2);
      nmax = nmax1;
   }
   else
   {
      if (nmax1 > nmax2)
      {
	 k = nmax2 * (nmax2 + 2);
	 l = nmax1 * (nmax1 + 2);
	 switch(gh)
	 {
	    case 3:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			gha[ii] = gh1[ii];
		     }
		     break;
	    case 4:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			ghb[ii] = gh1[ii];
		     }
		     break;
	    default: printf("\nError in subroutine extrapsh");
		     break;
	 }
	 nmax = nmax1;
      }
      else
      {
	 k = nmax1 * (nmax1 + 2);
	 l = nmax2 * (nmax2 + 2);
	 switch(gh)
	 {
	    case 3:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			gha[ii] = factor * gh2[ii];
		     }
		     break;
	    case 4:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			ghb[ii] = factor * gh2[ii];
		     }
		     break;
	    default: printf("\nError in subroutine extrapsh");
		     break;
	 }
	 nmax = nmax2;
      }
   }
   switch(gh)
   {
      case 3:  for ( ii = 1; ii <= k; ++ii)
	       {
		  gha[ii] = gh1[ii] + factor * gh2[ii];
	       }
	       break;
      case 4:  for ( ii = 1; ii <= k; ++ii)
	       {
		  ghb[ii] = gh1[ii] + factor * gh2[ii];
	       }
	       break;
      default: printf("\nError in subroutine extrapsh");
	       break;
   }
   return(nmax);
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine interpsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Interpolates linearly, in time, between two spherical harmonic       */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of earlier model                               */
/*           nmax1    - maximum degree and order of earlier model           */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of earlier model              */
/*           dte2     - date of later model                                 */
/*           nmax2    - maximum degree and order of later model             */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of internal model             */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - coefficients of resulting model                     */
/*           nmax     - maximum degree and order of resulting model         */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/


int interpsh(date, dte1, nmax1, dte2, nmax2, gh)
float date;
float dte1;
int   nmax1;
float dte2;
int   nmax2;
int   gh;
{
   int   nmax;
   int   k, l;
   int   ii;
   float factor;

   factor = (date - dte1) / (dte2 - dte1);
   if (nmax1 == nmax2)
   {
      k =  nmax1 * (nmax1 + 2);
      nmax = nmax1;
   }
   else
   {
      if (nmax1 > nmax2)
      {
	 k = nmax2 * (nmax2 + 2);
	 l = nmax1 * (nmax1 + 2);
	 switch(gh)
	 {
	    case 3:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			gha[ii] = gh1[ii] + factor * (-gh1[ii]);
		     }
		     break;
	    case 4:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			ghb[ii] = gh1[ii] + factor * (-gh1[ii]);
		     }
		     break;
	    default: printf("\nError in subroutine extrapsh");
		     break;
	 }
	 nmax = nmax1;
      }
      else
      {
	 k = nmax1 * (nmax1 + 2);
	 l = nmax2 * (nmax2 + 2);
	 switch(gh)
	 {
	    case 3:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			gha[ii] = factor * gh2[ii];
		     }
		     break;
	    case 4:  for ( ii = k + 1; ii <= l; ++ii)
		     {
			ghb[ii] = factor * gh2[ii];
		     }
		     break;
	    default: printf("\nError in subroutine extrapsh");
		     break;
	 }
	 nmax = nmax2;
      }
   }
   switch(gh)
   {
      case 3:  for ( ii = 1; ii <= k; ++ii)
	       {
		  gha[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
	       }
	       break;
      case 4:  for ( ii = 1; ii <= k; ++ii)
	       {
		  ghb[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
	       }
	       break;
      default: printf("\nError in subroutine extrapsh");
	       break;
   }
   return(nmax);
}





/****************************************************************************/
/*                                                                          */
/*                           Subroutine shval3                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Calculates field components from spherical harmonic (sh)             */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           igdgc     - indicates coordinate system used; set equal        */
/*                       to 1 if geodetic, 2 if geocentric                  */
/*           latitude  - north latitude, in degrees                         */
/*           longitude - east longitude, in degrees                         */
/*           elev      - WGS84 altitude above mean sea level (igdgc=1), or  */
/*                       radial distance from earth's center (igdgc=2)      */
/*           a2,b2     - squares of semi-major and semi-minor axes of       */
/*                       the reference spheroid used for transforming       */
/*                       between geodetic and geocentric coordinates        */
/*                       or components                                      */
/*           nmax      - maximum degree and order of coefficients           */
/*           iext      - external coefficients flag (=0 if none)            */
/*           ext1,2,3  - the three 1st-degree external coefficients         */
/*                       (not used if iext = 0)                             */
/*                                                                          */
/*     Output:                                                              */
/*           x         - northward component                                */
/*           y         - eastward component                                 */
/*           z         - vertically-downward component                      */
/*                                                                          */
/*     based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,  */
/*     report no. 71/1, institute of geological sciences, U.K.              */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Norman W. Peddie                                               */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/


int shval3(igdgc, flat, flon, elev, nmax, gh, iext, ext1, ext2, ext3)
int   igdgc;
float flat;
float flon;
float elev;
int   nmax;
int   gh;
int   iext;
float ext1;
float ext2;
float ext3;
{
   float earths_radius = 6371.2;
   float dtr = 0.01745329;
   float slat;
   float clat;
   float ratio;
   float aa, bb, cc, dd;
   float sd;
   float cd;
   float r;
   float a2;
   float b2;
   float rr;
   float fm,fn;
   float sl[14];
   float cl[14];
   float p[119];
   float q[119];
   int ii,j,k,l,m,n;
   int npq;
   int ios;
   double arguement;
   double power;
   a2 = 40680631.59;            /* WGS84 */
   b2 = 40408299.98;            /* WGS84 */
   ios = 0;
   r = elev;
   arguement = flat * dtr;
   slat = sin( arguement );
   if ((90.0 - flat) < 0.001)
   {
      aa = 89.999;            /*  300 ft. from North pole  */
   }
   else
   {
      if ((90.0 + flat) < 0.001)
      {
	 aa = -89.999;        /*  300 ft. from South pole  */
      }
      else
      {
	 aa = flat;
      }
   }
   arguement = aa * dtr;
   clat = cos( arguement );
   arguement = flon * dtr;
   sl[1] = sin( arguement );
   cl[1] = cos( arguement );
   switch(gh)
   {
      case 3:  x = 0;
	       y = 0;
	       z = 0;
	       break;
      case 4:  xtemp = 0;
	       ytemp = 0;
	       ztemp = 0;
	       break;
      default: printf("\nError in subroutine shval3");
	       break;
   }
   sd = 0.0;
   cd = 1.0;
   l = 1;
   n = 0;
   m = 1;
   npq = (nmax * (nmax + 3)) / 2;
   if (igdgc == 1)
   {
      aa = a2 * clat * clat;
      bb = b2 * slat * slat;
      cc = aa + bb;
      arguement = cc;
      dd = sqrt( arguement );
      arguement = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
      r = sqrt( arguement );
      cd = (elev + dd) / r;
      sd = (a2 - b2) / dd * slat * clat / r;
      aa = slat;
      slat = slat * cd - clat * sd;
      clat = clat * cd + aa * sd;
   }
   ratio = earths_radius / r;
   arguement = 3.0;
   aa = sqrt( arguement );
   p[1] = 2.0 * slat;
   p[2] = 2.0 * clat;
   p[3] = 4.5 * slat * slat - 1.5;
   p[4] = 3.0 * aa * clat * slat;
   q[1] = -clat;
   q[2] = slat;
   q[3] = -3.0 * clat * slat;
   q[4] = aa * (slat * slat - clat * clat);
   for ( k = 1; k <= npq; ++k)
   {
      if (n < m)
      {
	 m = 0;
	 n = n + 1;
	 arguement = ratio;
	 power =  n + 2;
	 rr = pow(arguement,power);
	 fn = n;
      }
      fm = m;
      if (k >= 5)
      {
	 if (m == n)
	 {
	    arguement = (1.0 - 0.5/fm);
	    aa = sqrt( arguement );
	    j = k - n - 1;
	    p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
	    q[k] = aa * (clat * q[j] + slat/fm * p[j]);
	    sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
	    cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
	 }
	 else
	 {
	    arguement = fn*fn - fm*fm;
	    aa = sqrt( arguement );
	    arguement = ((fn - 1.0)*(fn-1.0)) - (fm * fm);
	    bb = sqrt( arguement )/aa;
	    cc = (2.0 * fn - 1.0)/aa;
	    ii = k - n;
	    j = k - 2 * n + 1;
	    p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
	    q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
	 }
      }
      switch(gh)
      {
	 case 3:  aa = rr * gha[l];
		  break;
	 case 4:  aa = rr * ghb[l];
		  break;
	 default: printf("\nError in subroutine shval3");
		  break;
      }
      if (m == 0)
      {
	 switch(gh)
	 {
	    case 3:  x = x + aa * q[k];
		     z = z - aa * p[k];
		     break;
	    case 4:  xtemp = xtemp + aa * q[k];
		     ztemp = ztemp - aa * p[k];
		     break;
	    default: printf("\nError in subroutine shval3");
		     break;
	 }
	 l = l + 1;
      }
      else
      {
	 switch(gh)
	 {
	    case 3:  bb = rr * gha[l+1];
		     cc = aa * cl[m] + bb * sl[m];
		     x = x + cc * q[k];
		     z = z - cc * p[k];
		     if (clat > 0)
		     {
			y = y + (aa * sl[m] - bb * cl[m]) *
			    fm * p[k]/((fn + 1.0) * clat);
		     }
		     else
		     {
			y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
		     }
		     l = l + 2;
		     break;
	    case 4:  bb = rr * ghb[l+1];
		     cc = aa * cl[m] + bb * sl[m];
		     xtemp = xtemp + cc * q[k];
		     ztemp = ztemp - cc * p[k];
		     if (clat > 0)
		     {
			ytemp = ytemp + (aa * sl[m] - bb * cl[m]) *
			    fm * p[k]/((fn + 1.0) * clat);
		     }
		     else
		     {
			ytemp = ytemp + (aa * sl[m] - bb * cl[m]) *
				q[k] * slat;
		     }
		     l = l + 2;
		     break;
	    default: printf("\nError in subroutine shval3");
		     break;
	 }
      }
      m = m + 1;
   }
   if (iext != 0)
   {
      aa = ext2 * cl[1] + ext3 * sl[1];
      switch(gh)
      {
	 case 3:   x = x - ext1 * clat + aa * slat;
		   y = y + ext2 * sl[1] - ext3 * cl[1];
		   z = z + ext1 * slat + aa * clat;
		   break;
	 case 4:   xtemp = xtemp - ext1 * clat + aa * slat;
		   ytemp = ytemp + ext2 * sl[1] - ext3 * cl[1];
		   ztemp = ztemp + ext1 * slat + aa * clat;
		   break;
	 default:  printf("\nError in subroutine shval3");
		   break;
      }
   }
   switch(gh)
   {
      case 3:   aa = x;
		x = x * cd + z * sd;
		z = z * cd - aa * sd;
		break;
      case 4:   aa = xtemp;
		xtemp = xtemp * cd + ztemp * sd;
		ztemp = ztemp * cd - aa * sd;
		break;
      default:  printf("\nError in subroutine shval3");
		break;
   }
   return(ios);
}


/****************************************************************************/
/*                                                                          */
/*                           Subroutine dihf                                */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the geomagnetic d, i, h, and f from x, y, and z.            */
/*                                                                          */
/*     Input:                                                               */
/*           x  - northward component                                       */
/*           y  - eastward component                                        */
/*           z  - vertically-downward component                             */
/*                                                                          */
/*     Output:                                                              */
/*           d  - declination                                               */
/*           i  - inclination                                               */
/*           h  - horizontal intensity                                      */
/*           f  - total intensity                                           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 22, 1988                                                */
/*                                                                          */
/****************************************************************************/

int dihf (gh)
int gh;
{
   int ios;
   int j;
   float sn;
   float h2;
   float hpx;
   double arguement, arguement2;
   double rad, pi;

   ios = gh;
   sn = 0.0001;
	 rad = 57.29577951;
	 pi = 3.141592654;
   switch(gh)
   {
      case 3:   for (j = 1; j <= 1; ++j)
		{
		   h2 = x*x + y*y;
		   arguement = h2;
		   h = sqrt(arguement);       /* calculate horizontal intensity */
		   arguement = h2 + z*z;
		   f = sqrt(arguement);      /* calculate total intensity */
		   if (f < sn)
		   {
		      d = NaN;        /* If d and i cannot be determined, */
		      i = NaN;        /*       set equal to NaN         */
		   }
		   else
		   {
		      arguement = z;
		      arguement2 = h;
		      i = atan2(arguement,arguement2);
		      if (h < sn)
		      {
              d = NaN;
		      }
		      else
		      {
			 hpx = h + x;
			 if (hpx < sn)
			 {
			    d = pi;
			 }
			 else
			 {
			    arguement = y;
			    arguement2 = hpx;
			    d = 2.0 * atan2(arguement,arguement2);
			 }
		      }
		   }
		}
		break;
      case 4:   for (j = 1; j <= 1; ++j)
		{
		   h2 = xtemp*xtemp + ytemp*ytemp;
		   arguement = h2;
		   htemp = sqrt(arguement);
		   arguement = h2 + ztemp*ztemp;
		   ftemp = sqrt(arguement);
		   if (ftemp < sn)
		   {
		      dtemp = NaN;    /* If d and i cannot be determined, */
		      itemp = NaN;    /*       set equal to 999.0         */
		   }
		   else
		   {
		      arguement = ztemp;
		      arguement2 = htemp;
		      itemp = atan2(arguement,arguement2);
		      if (htemp < sn)
		      {
              dtemp = NaN;
		      }
		      else
		      {
              hpx = htemp + xtemp;
              if (hpx < sn)
                {
                  dtemp = pi;
                }
              else
                {
			    arguement = ytemp;
			    arguement2 = hpx;
			    dtemp = 2.0 * atan2(arguement,arguement2);
			 }
		      }
		   }
		}
		break;
      default:  printf("\nError in subroutine dihf");
		break;
   }
return(ios);
}
