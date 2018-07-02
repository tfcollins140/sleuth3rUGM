/*******************************************************************************

  MODULE:                   stats_obj.c

  CONTAINING SYSTEM:        SLEUTH-3r (based on SLEUTH Model 3.0 Beta)
                            (Slope, Land-cover, Exclusion, Urbanization, 
                            Transportation, and Hillshade)
                            also known as UGM 3.0 Beta (for Urban Growth Model)

  VERSION:                  SLEUTH-3r [Includes Version D features]

  REVISION DATE:            August 31, 2006a
                            [Annotations added August 19, 2009]

  PURPOSE:

     This module provides functions which compute, accumulate, and provide
     statistics on model processes.

  NOTES:

  MODIFICATIONS:

  TO DO:

**************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "ugm_defines.h"
#include "pgrid_obj.h"
#include "proc_obj.h"
#include "igrid_obj.h"
#include "memory_obj.h"
#include "scenario_obj.h"
#include "ugm_macros.h"
#include "stats_obj.h"
#include "coeff_obj.h"
#include "utilities.h"

  /*VerD*/
  extern FILE *fpVerD2;
  extern FILE *fpVerD3;
  extern FILE *fpVerD4;
  
  extern BOOLEAN WriteSlopeFileFlag;
  extern BOOLEAN WriteRatioFileFlag;
  extern BOOLEAN WriteXypointsFileFlag;
  /*VerD*/

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                               SCCS ID                                     **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
char stats_obj_c_sccs_id[] = "@(#)stats_obj.c	1.72	12/4/00";

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                                 MACROS                                    **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
#define MAX_LINE_LEN 256
#define SIZE_CIR_Q 6000   /*VerD*/

#define Q_STORE(R,C)                                                     \
  if((sidx+1==ridx)||((sidx+1==SIZE_CIR_Q)&&!ridx)){                     \
    printf("Error Circular Queue Full\n");                               \
    printf("Increase SIZE_CIR_Q and recompile\n");                       \
    printf("sidx=%d ridx=%d SIZE_CIR_Q=%d\n",sidx,ridx,SIZE_CIR_Q);      \
    EXIT(1);}                                                            \
  cir_q[sidx].row = R;                                                   \
  cir_q[sidx].col = C;                                                   \
  sidx++;                                                                \
  depth++;                                                               \
  sidx %= SIZE_CIR_Q
#define Q_RETREIVE(R,C)                                                  \
  ridx = ridx%SIZE_CIR_Q;                                                \
  R = cir_q[ridx].row;                                                   \
  C = cir_q[ridx].col;                                                   \
  ridx++;                                                                \
  depth--


/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                      STATIC MEMORY FOR THIS OBJECT                        **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
static char *stats_val_t_names[] = {
  "sng",
  "sdg",
  "sdc",
  "og",
  "rt",
  "pop",
  "area",
  "edges",
  "clusters",
  "xmean",
  "ymean",
  "rad",
  "slope",
  "cl_size",
  "diffus",
  "spread",
  "breed",
  "slp_res",
  "rd_grav",
  "%urban",
  "%road",
  "grw_rate",
  "leesalee",
  "grw_pix"
};
static stats_info stats_actual[MAX_URBAN_YEARS];
static stats_info regression;
static stats_val_t average[MAX_URBAN_YEARS];
static stats_val_t std_dev[MAX_URBAN_YEARS];
static stats_val_t running_total[MAX_URBAN_YEARS];
static struct
{
  int run;
  int monte_carlo;
  int year;
  stats_val_t this_year;
}
record;

static struct
{
  double fmatch;
  double actual;
  double simulated;
  double compare;
  double leesalee;
  double product;
}
aggregate;

static struct
{
  long successes;
  long z_failure;
  long delta_failure;
  long slope_failure;
  long excluded_failure;
}
urbanization_attempt;

static int sidx;
static int ridx;

/* link element for Cluster routine */
typedef struct ugm_link
{
  int row;
  int col;
}
ugm_link;

static struct ugm_link cir_q[SIZE_CIR_Q];

/** D. Donato 8/21/2006 Added to deal with cumulative growth                 **/
  short *zgrwth_row;
  short *zgrwth_col;
  int    zgrwth_count;
/** D. Donato 8/21/2006 Added to deal with cumulative growth                 **/

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                        STATIC FUNCTION PROTOTYPES                         **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
static void stats_Save (char *filename);
static void stats_LogThisYearStats (FILE * fp);
static void stats_CalGrowthRate ();
static void stats_CalPercentUrban (int, int, int);
static void stats_CalAverages (int index);
static void stats_WriteControlStats (char *filename);
static void stats_WriteStatsValLine (char *filename, int run,
                           int year, stats_val_t * stats_ptr, int index);
static void stats_LogStatInfoHdr (FILE * fp);
static void stats_LogStatInfo (int run, int year, int index,
                               stats_info * stats_ptr, FILE * fp);
static void stats_LogStatVal (int run, int year, int index,
                              stats_val_t * stats_ptr, FILE * fp);
static void stats_LogStatValHdr (FILE * fp);
static void stats_ComputeThisYearStats ();
static void stats_SetNumGrowthPixels (int val);
static void stats_CalLeesalee ();
static void stats_ProcessGrowLog (int run, int year);
static void stats_DoAggregate (double fmatch);
static void stats_DoRegressions ();
static double stats_linefit (double *dependent,
                             double *independent,
                             int number_of_observations);
static void stats_LogControlStats (FILE * fp);
static void stats_LogControlStatsHdr (FILE * fp);
static void
    stats_compute_stats (GRID_P Z,                           /* IN     */
                         GRID_P slp,                         /* IN     */
                         double *stats_area,                 /* OUT    */
                         double *stats_edges,                /* OUT    */
                         double *stats_clusters,             /* OUT    */
                         double *stats_pop,                  /* OUT    */
                         double *stats_xmean,                /* OUT    */
                         double *stats_ymean,                /* OUT    */
                         double *stats_average_slope,        /* OUT    */
                         double *stats_rad,                  /* OUT    */
                         double *stats_mean_cluster_size,    /* OUT    */
                         GRID_P scratch_gif1,                /* MOD    */
                         GRID_P scratch_gif2);             /* MOD    */
static void
    stats_edge (GRID_P Z,                                    /* IN     */
                double *stats_area,                          /* OUT    */
                double *stats_edges);                      /* OUT    */
static void
    stats_circle (GRID_P Z,                                  /* IN     */
                  GRID_P slp,                                /* IN     */
                  int stats_area,                            /* IN     */
                  double *stats_xmean,                       /* OUT    */
                  double *stats_ymean,                       /* OUT    */
                  double *stats_average_slope,               /* OUT    */
                  double *stats_rad);                      /* OUT    */
static void
    stats_cluster (GRID_P Z,                                 /* IN     */
                   double *stats_clusters,                   /* OUT    */
                   double *stats_pop,                        /* OUT    */
                   double *stats_mean_cluster_size,          /* OUT    */
                   GRID_P scratch_gif1,                      /* MOD    */
                   GRID_P scratch_gif2);                   /* MOD    */
static void stats_ClearStatsValArrays ();
static void stats_ComputeBaseStats ();
static void stats_CalStdDev (int index);
static void
    stats_compute_leesalee (GRID_P Z,                        /* IN     */
                            GRID_P urban,                    /* IN     */
                            double *leesalee);             /* OUT    */




/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ConcatenateControlFiles
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
#if 1
void
  stats_ConcatenateControlFiles ()
{
  char func[] = "stats_ConcatenateControlFiles";
  char line[MAX_LINE_LEN];
  char source_file[MAX_FILENAME_LEN];
  char destination_file[MAX_FILENAME_LEN];
  char command[2 * MAX_FILENAME_LEN + 20];
  FILE *fp;
  FILE *source_fp;
  int line_count;
  int i;

  sprintf (destination_file, "%scontrol_stats.log", scen_GetOutputDir ());
  sprintf (source_file, "%scontrol_stats_pe_%u.log", scen_GetOutputDir (), 0);
  sprintf (command, "mv %s %s", source_file, destination_file);
  system (command);

  FILE_OPEN (fp, destination_file, "a");
  for (i = 1; i < glb_npes; i++)
  {
    sprintf (source_file, "%scontrol_stats_pe_%u.log", scen_GetOutputDir (), i);

    FILE_OPEN (source_fp, source_file, "r");

    line_count = 0;
    while (fgets (line, MAX_LINE_LEN, source_fp) != NULL)
    {
      line_count++;
      if (line_count <= 2)
        continue;
      fputs (line, fp);
    }
    fclose (source_fp);
    sprintf (command, "rm %s", source_file);
    printf ("%s %u command: %s\n", __FILE__, __LINE__, command);
    system (command);

  }
  fclose (fp);
}
#else
void
  stats_ConcatenateControlFiles (int current_run)
{
  char func[] = "stats_ConcatenateControlFiles";
  char line[MAX_LINE_LEN];
  char source_file[MAX_FILENAME_LEN];
  char destination_file[MAX_FILENAME_LEN];
  char command[2 * MAX_FILENAME_LEN + 20];
  FILE *fp;
  FILE *source_fp;
  int line_count;

  sprintf (destination_file, "%scontrol_stats.log", scen_GetOutputDir ());
  sprintf (source_file, "%scontrol_stats_pe_%u.log", scen_GetOutputDir (), 0);

  FILE_OPEN (fp, destination_file, "a");

  sprintf (source_file, "%scontrol_stats_pe_%u.log",
           scen_GetOutputDir (), current_run);

  FILE_OPEN (source_fp, source_file, "r");

  line_count = 0;
  while (fgets (line, MAX_LINE_LEN, source_fp) != NULL)
  {
    line_count++;
    if (line_count <= 2)
      continue;
    fputs (line, fp);
  }
  fclose (source_fp);
  fclose (fp);
  sprintf (command, "rm %s", source_file);
  system (command);
}
#endif
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ConcatenateStdDevFiles () 
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
#if 1
void
  stats_ConcatenateStdDevFiles ()
{
  char func[] = "stats_ConcatenateStdDevFiles";
  char line[MAX_LINE_LEN];
  char source_file[MAX_FILENAME_LEN];
  char destination_file[MAX_FILENAME_LEN];
  char command[2 * MAX_FILENAME_LEN + 20];
  FILE *fp;
  FILE *source_fp;
  int line_count;
  int i;

  sprintf (destination_file, "%sstd_dev.log", scen_GetOutputDir ());
  sprintf (source_file, "%sstd_dev_pe_%u.log", scen_GetOutputDir (), 0);
  sprintf (command, "mv %s %s", source_file, destination_file);
  system (command);

  FILE_OPEN (fp, destination_file, "a");
  for (i = 1; i < glb_npes; i++)
  {
    sprintf (source_file, "%sstd_dev_pe_%u.log", scen_GetOutputDir (), i);

    FILE_OPEN (source_fp, source_file, "r");

    line_count = 0;
    while (fgets (line, MAX_LINE_LEN, source_fp) != NULL)
    {
      line_count++;
      if (line_count <= 1)
        continue;
      fputs (line, fp);
    }
    fclose (source_fp);
    sprintf (command, "rm %s", source_file);
    printf ("%s %u command: %s\n", __FILE__, __LINE__, command);
    system (command);

  }
  fclose (fp);
#else
void
  stats_ConcatenateStdDevFiles (int current_run)
{
  char func[] = "stats_ConcatenateStdDevFiles";
  char line[MAX_LINE_LEN];
  char source_file[MAX_FILENAME_LEN];
  char destination_file[MAX_FILENAME_LEN];
  char command[2 * MAX_FILENAME_LEN + 20];
  FILE *fp;
  FILE *source_fp;
  int line_count;


  sprintf (destination_file, "%sstd_dev.log", scen_GetOutputDir ());
  sprintf (source_file, "%sstd_dev_pe_%u.log", scen_GetOutputDir (), 0);

  FILE_OPEN (fp, destination_file, "a");

  sprintf (source_file, "%sstd_dev_pe_%u.log",
           scen_GetOutputDir (), current_run);

  FILE_OPEN (source_fp, source_file, "r");

  line_count = 0;
  while (fgets (line, MAX_LINE_LEN, source_fp) != NULL)
  {
    line_count++;
    if (line_count <= 1)
      continue;
    fputs (line, fp);
  }
  fclose (source_fp);
  fclose (fp);
  sprintf (command, "rm %s", source_file);
  system (command);
#endif
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ConcatenateAvgFiles
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
#if 1
void
  stats_ConcatenateAvgFiles ()
{

  char func[] = "stats_ConcatenateAvgFiles";
  char line[MAX_LINE_LEN];
  char source_file[MAX_FILENAME_LEN];
  char destination_file[MAX_FILENAME_LEN];
  char command[2 * MAX_FILENAME_LEN + 20];
  FILE *fp;
  FILE *source_fp;
  int line_count;
  int i;

  sprintf (destination_file, "%savg.log", scen_GetOutputDir ());
  sprintf (source_file, "%savg_pe_%u.log", scen_GetOutputDir (), 0);
  sprintf (command, "mv %s %s", source_file, destination_file);
  system (command);

  FILE_OPEN (fp, destination_file, "a");
  for (i = 1; i < glb_npes; i++)
  {
    sprintf (source_file, "%savg_pe_%u.log", scen_GetOutputDir (), i);

    FILE_OPEN (source_fp, source_file, "r");

    line_count = 0;
    while (fgets (line, MAX_LINE_LEN, source_fp) != NULL)
    {
      line_count++;
      if (line_count <= 1)
        continue;
      fputs (line, fp);
    }
    fclose (source_fp);
    sprintf (command, "rm %s", source_file);
    printf ("%s %u command: %s\n", __FILE__, __LINE__, command);
    system (command);

  }
  fclose (fp);
#else

void
  stats_ConcatenateAvgFiles (int current_run)
{
  char func[] = "stats_ConcatenateAvgFiles";
  char line[MAX_LINE_LEN];
  char source_file[MAX_FILENAME_LEN];
  char destination_file[MAX_FILENAME_LEN];
  char command[2 * MAX_FILENAME_LEN + 20];
  FILE *fp;
  FILE *source_fp;
  int line_count;

  sprintf (destination_file, "%savg.log", scen_GetOutputDir ());
  sprintf (source_file, "%savg_pe_%u.log", scen_GetOutputDir (), 0);

  FILE_OPEN (fp, destination_file, "a");

  sprintf (source_file, "%savg_pe_%u.log", scen_GetOutputDir (), current_run);

  FILE_OPEN (source_fp, source_file, "r");

  line_count = 0;
  while (fgets (line, MAX_LINE_LEN, source_fp) != NULL)
  {
    line_count++;
    if (line_count <= 1)
      continue;
    fputs (line, fp);
  }
  fclose (source_fp);
  fclose (fp);
  sprintf (command, "rm %s", source_file);
#if 1
  system (command);
#endif
#endif
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_Update
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/



void
  stats_Update (int num_growth_pix)
{
  char func[] = "stats_Update";
  char filename[MAX_FILENAME_LEN];
  int total_pixels;
  int road_pixel_count;
  int excluded_pixel_count;

  total_pixels = mem_GetTotalPixels ();
  road_pixel_count = igrid_GetIGridRoadPixelCount (proc_GetCurrentYear ());
  excluded_pixel_count = igrid_GetIGridExcludedPixelCount ();

  stats_ComputeThisYearStats ();
  stats_SetNumGrowthPixels (num_growth_pix);
  stats_CalGrowthRate ();
  stats_CalPercentUrban (total_pixels, road_pixel_count, excluded_pixel_count);

  if (igrid_TestForUrbanYear (proc_GetCurrentYear ()))
  {
    stats_CalLeesalee ();
    sprintf (filename, "%sgrow_%u_%u.log",
    scen_GetOutputDir (), proc_GetCurrentRun (), proc_GetCurrentYear ());

    /*VerD*/
    if (proc_GetProcessingType () != PREDICTING)
	{
	    if (WriteXypointsFileFlag == 1)
		{
			fprintf (fpVerD2, "%6d %6d %6.0f %6.0f %6.0f %6.0f %6.0f %6d %6d\n",
				proc_GetCurrentRun(),
				proc_GetCurrentMonteCarlo(),

				coeff_GetCurrentDiffusion(),
				coeff_GetCurrentBreed(),
				coeff_GetCurrentSpread(),
				coeff_GetCurrentSlopeResist(),
				coeff_GetCurrentRoadGravity(),

				proc_GetCurrentYear (),

				stats_GetArea()      );
		}
	}
    /*VerD*/


    stats_Save (filename);
  }
  if (proc_GetProcessingType () == PREDICTING)
  {
    sprintf (filename, "%sgrow_%u_%u.log", scen_GetOutputDir (),
             proc_GetCurrentRun (), proc_GetCurrentYear ());
    stats_Save (filename);
  }
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_Init
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_Init ()
{
  static BOOLEAN first_call = TRUE;

  stats_ClearStatsValArrays ();
  if (first_call)
  {
    stats_ComputeBaseStats ();
    first_call = FALSE;
  }
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CalStdDev
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_CalStdDev (int index)
{
#define SD(val) pow(((val)*(val)/total_monte_carlo),0.5)
  int total_monte_carlo;

  total_monte_carlo = scen_GetMonteCarloIterations ();

  std_dev[index].sng = SD (record.this_year.sng - average[index].sng);
  std_dev[index].sdg = SD (record.this_year.sdg - average[index].sdg);
  std_dev[index].sdc = SD (record.this_year.sdc - average[index].sdc);
  std_dev[index].og = SD (record.this_year.og - average[index].og);
  std_dev[index].rt = SD (record.this_year.rt - average[index].rt);
  std_dev[index].pop = SD (record.this_year.pop - average[index].pop);
  std_dev[index].area = SD (record.this_year.area - average[index].area);
  std_dev[index].edges = SD (record.this_year.edges - average[index].edges);
  std_dev[index].clusters =
    SD (record.this_year.clusters - average[index].clusters);
  std_dev[index].xmean = SD (record.this_year.xmean - average[index].xmean);
  std_dev[index].ymean = SD (record.this_year.ymean - average[index].ymean);
  std_dev[index].rad = SD (record.this_year.rad - average[index].rad);
  std_dev[index].slope = SD (record.this_year.slope - average[index].slope);
  std_dev[index].mean_cluster_size =
    SD (record.this_year.mean_cluster_size - average[index].mean_cluster_size);
  std_dev[index].diffusion =
    SD (record.this_year.diffusion - average[index].diffusion);
  std_dev[index].spread = SD (record.this_year.spread - average[index].spread);
  std_dev[index].breed = SD (record.this_year.breed - average[index].breed);
  std_dev[index].slope_resistance =
    SD (record.this_year.slope_resistance - average[index].slope_resistance);
  std_dev[index].road_gravity =
    SD (record.this_year.road_gravity - average[index].road_gravity);
  std_dev[index].percent_urban =
    SD (record.this_year.percent_urban - average[index].percent_urban);
  std_dev[index].percent_road =
    SD (record.this_year.percent_road - average[index].percent_road);
  std_dev[index].growth_rate =
    SD (record.this_year.growth_rate - average[index].growth_rate);
  std_dev[index].leesalee =
    SD (record.this_year.leesalee - average[index].leesalee);
  std_dev[index].num_growth_pix =
    SD (record.this_year.num_growth_pix - average[index].num_growth_pix);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CalAverages
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_CalAverages (int index)
{
  int total_monte_carlo;

  total_monte_carlo = scen_GetMonteCarloIterations ();

  average[index].sng = running_total[index].sng / total_monte_carlo;
  average[index].sdg = running_total[index].sdg / total_monte_carlo;
  average[index].sdc = running_total[index].sdc / total_monte_carlo;
  average[index].og = running_total[index].og / total_monte_carlo;
  average[index].rt = running_total[index].rt / total_monte_carlo;
  average[index].pop = running_total[index].pop / total_monte_carlo;
  average[index].area = running_total[index].area / total_monte_carlo;
  average[index].edges = running_total[index].edges / total_monte_carlo;
  average[index].clusters = running_total[index].clusters / total_monte_carlo;
  average[index].xmean = running_total[index].xmean / total_monte_carlo;
  average[index].ymean = running_total[index].ymean / total_monte_carlo;
  average[index].rad = running_total[index].rad / total_monte_carlo;
  average[index].slope = running_total[index].slope / total_monte_carlo;
  average[index].mean_cluster_size =
    running_total[index].mean_cluster_size / total_monte_carlo;
  average[index].diffusion =
    running_total[index].diffusion / total_monte_carlo;
  average[index].spread = running_total[index].spread / total_monte_carlo;
  average[index].breed = running_total[index].breed / total_monte_carlo;
  average[index].slope_resistance =
    running_total[index].slope_resistance / total_monte_carlo;
  average[index].road_gravity =
    running_total[index].road_gravity / total_monte_carlo;
  average[index].percent_urban =
    running_total[index].percent_urban / total_monte_carlo;
  average[index].percent_road =
    running_total[index].percent_road / total_monte_carlo;
  average[index].growth_rate =
    running_total[index].growth_rate / total_monte_carlo;
  average[index].leesalee = running_total[index].leesalee / total_monte_carlo;
  average[index].num_growth_pix =
    running_total[index].num_growth_pix / total_monte_carlo;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_UpdateRunningTotal
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

static void
  stats_UpdateRunningTotal (int index)
{
#ifndef lint
  char func[] = "stats_UpdateRunningTotal";
#endif

  running_total[index].sng += record.this_year.sng;
  running_total[index].sdg += record.this_year.sdg;
  running_total[index].sdc += record.this_year.sdc;
  running_total[index].og += record.this_year.og;
  running_total[index].rt += record.this_year.rt;
  running_total[index].pop += record.this_year.pop;
  running_total[index].area += record.this_year.area;
  running_total[index].edges += record.this_year.edges;
  running_total[index].clusters += record.this_year.clusters;
  running_total[index].xmean += record.this_year.xmean;
  running_total[index].ymean += record.this_year.ymean;
  running_total[index].rad += record.this_year.rad;
  running_total[index].slope += record.this_year.slope;
  running_total[index].mean_cluster_size += record.this_year.mean_cluster_size;
  running_total[index].diffusion += record.this_year.diffusion;
  running_total[index].spread += record.this_year.spread;
  running_total[index].breed += record.this_year.breed;
  running_total[index].slope_resistance += record.this_year.slope_resistance;
  running_total[index].road_gravity += record.this_year.road_gravity;
  running_total[index].percent_urban += record.this_year.percent_urban;
  running_total[index].percent_road += record.this_year.percent_road;
  running_total[index].growth_rate += record.this_year.growth_rate;
  running_total[index].leesalee += record.this_year.leesalee;
  running_total[index].num_growth_pix += record.this_year.num_growth_pix;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ClearStatsValArrays
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

static void
  stats_ClearStatsValArrays ()
{
  char func[] = "stats_ClearStatsValArrays";
  int i;

  for (i = 0; i < MAX_URBAN_YEARS; i++)
  {
    memset ((void *) (&running_total[i]), 0, sizeof (stats_val_t));
    memset ((void *) (&average[i]), 0, sizeof (stats_val_t));
    memset ((void *) (&std_dev[i]), 0, sizeof (stats_val_t));
  }
  memset ((void *) (&regression), 0, sizeof (stats_info));
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetLeesalee
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetLeesalee ()
{
  return record.this_year.leesalee;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CalLeesalee
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_CalLeesalee ()
{
  char func[] = "stats_CalLeesalee";
  GRID_P z_ptr;
  GRID_P urban_ptr;

  z_ptr = pgrid_GetZPtr ();
  urban_ptr = igrid_GetUrbanGridPtrByYear (__FILE__, func,
                                       __LINE__, proc_GetCurrentYear ());
  record.this_year.leesalee = 1.0;
  if (proc_GetProcessingType () != PREDICTING)
  {
    stats_compute_leesalee (z_ptr,                           /* IN     */
                            urban_ptr,                       /* IN     */
                            &record.this_year.leesalee);   /* OUT    */
  }
  urban_ptr = igrid_GridRelease (__FILE__, func, __LINE__, urban_ptr);

}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetNumGrowthPixels
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_SetNumGrowthPixels (int val)
{
  record.this_year.num_growth_pix = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetNumGrowthPixels
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetNumGrowthPixels ()
{
  return record.this_year.num_growth_pix;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetPercentUrban
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

void
  stats_SetPercentUrban (int val)
{
  record.this_year.percent_urban = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CalPercentUrban
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_CalPercentUrban (int total_pixels, int road_pixels, int excld_pixels)
{
  record.this_year.percent_urban =
    (double) (100.0 * (record.this_year.pop + road_pixels) /
              (total_pixels - road_pixels - excld_pixels));
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetPercentUrban
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetPercentUrban ()
{
  return record.this_year.percent_urban;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CalGrowthRate
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

static void
  stats_CalGrowthRate ()
{
  record.this_year.growth_rate =
    record.this_year.num_growth_pix / record.this_year.pop * 100.0;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetGrowthRate
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetGrowthRate ()
{
  return record.this_year.growth_rate;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetSNG
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/


void
  stats_SetSNG (int val)
{
  record.this_year.sng = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetSDG
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetSDG (int val)
{
  record.this_year.sdg = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetOG
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetOG (int val)
{
  record.this_year.og = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetRT
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetRT (int val)
{
  record.this_year.rt = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetPOP
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetPOP (int val)
{
  record.this_year.pop = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetSNG
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

int
  stats_GetSNG ()
{
  return record.this_year.sng;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetSDG
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetSDG ()
{
  return record.this_year.sdg;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetOG
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetOG ()
{
  return record.this_year.og;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetRT
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetRT ()
{
  return record.this_year.rt;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetPOP
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetPOP ()
{
  return record.this_year.pop;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetArea
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

void
  stats_SetArea (int val)
{
  record.this_year.area = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetEdges
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetEdges (int val)
{
  record.this_year.edges = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetClusters
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetClusters (int val)
{
  record.this_year.clusters = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetPop
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetPop (int val)
{
  record.this_year.pop = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetXmean
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetXmean (double val)
{
  record.this_year.xmean = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetYmean
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetYmean (double val)
{
  record.this_year.ymean = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetRad
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetRad (double val)
{
  record.this_year.rad = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetAvgSlope
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetAvgSlope (double val)
{
  record.this_year.slope = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_SetMeanClusterSize
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_SetMeanClusterSize (double val)
{
  record.this_year.mean_cluster_size = val;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetArea
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

int
  stats_GetArea ()
{
  return record.this_year.area;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetEdges 
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetEdges ()
{
  return record.this_year.edges;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetClusters
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetClusters ()
{
  return record.this_year.clusters;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetPop
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  stats_GetPop ()
{
  return record.this_year.pop;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetXmean
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetXmean ()
{
  return record.this_year.xmean;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetYmean
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetYmean ()
{
  return record.this_year.ymean;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetRad
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetRad ()
{
  return record.this_year.rad;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetAvgSlope
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetAvgSlope ()
{
  return record.this_year.slope;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_GetMeanClusterSize
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
double
  stats_GetMeanClusterSize ()
{
  return record.this_year.mean_cluster_size;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ComputeThisYearStats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

static void
  stats_ComputeThisYearStats ()
{
  char func[] = "stats_ComputeThisYearStats";
  int total_pixels;
  GRID_P z_ptr;
  GRID_P slope_ptr;
  GRID_P stats_workspace1;
  GRID_P stats_workspace2;
  double area;
  double edges;
  double clusters;
  double pop;
  double xmean;
  double ymean;
  double slope;
  double rad;
  double mean_cluster_size;

  total_pixels = mem_GetTotalPixels ();
  assert (total_pixels > 0);
  z_ptr = pgrid_GetZPtr ();
  assert (z_ptr != NULL);

  slope_ptr = igrid_GetSlopeGridPtr (__FILE__, func, __LINE__);
  stats_workspace1 = mem_GetWGridPtr (__FILE__, func, __LINE__);
  stats_workspace2 = mem_GetWGridPtr (__FILE__, func, __LINE__);

  stats_compute_stats (z_ptr,                                /* IN     */
                       slope_ptr,                            /* IN     */
                       &area,                                /* OUT    */
                       &edges,                               /* OUT    */
                       &clusters,                            /* OUT    */
                       &pop,                                 /* OUT    */
                       &xmean,                               /* OUT    */
                       &ymean,                               /* OUT    */
                       &slope,                               /* OUT    */
                       &rad,                                 /* OUT    */
                       &mean_cluster_size,                   /* OUT    */
                       stats_workspace1,                     /* MOD    */
                       stats_workspace2);                  /* MOD    */
  record.this_year.area = area;
  record.this_year.edges = edges;
  record.this_year.clusters = clusters;
  record.this_year.pop = pop;
  record.this_year.xmean = xmean;
  record.this_year.ymean = ymean;
  record.this_year.slope = slope;
  record.this_year.rad = rad;
  record.this_year.mean_cluster_size = mean_cluster_size;
  record.this_year.diffusion = coeff_GetCurrentDiffusion ();
  record.this_year.spread = coeff_GetCurrentSpread ();
  record.this_year.breed = coeff_GetCurrentBreed ();
  record.this_year.slope_resistance = coeff_GetCurrentSlopeResist ();
  record.this_year.road_gravity = coeff_GetCurrentRoadGravity ();

  slope_ptr = igrid_GridRelease (__FILE__, func, __LINE__, slope_ptr);
  stats_workspace1 = mem_GetWGridFree (__FILE__, func, __LINE__,
                                       stats_workspace1);
  stats_workspace2 = mem_GetWGridFree (__FILE__, func, __LINE__,
                                       stats_workspace2);


}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CreateControlFile
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_CreateControlFile (char *filename)
{
  char func[] = "stats_CreateControlFile";
  FILE *fp;

  FILE_OPEN (fp, filename, "w");

  stats_LogControlStatsHdr (fp);
  fclose (fp);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_CreateStatsValFile
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_CreateStatsValFile (char *filename)
{
  char func[] = "stats_CreateStatsValFile";
  FILE *fp;

  FILE_OPEN (fp, filename, "w");

  stats_LogStatValHdr (fp);
  fclose (fp);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_WriteStatsValLine
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_WriteStatsValLine (char *filename, int run, int year,
                           stats_val_t * stats_ptr, int index)
{
  char func[] = "stats_WriteStatsValLine";
  FILE *fp;

  FILE_OPEN (fp, filename, "a");

  stats_LogStatVal (run, year, index, &(stats_ptr[index]), fp);
  fclose (fp);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogStatInfoHdr
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogStatInfoHdr (FILE * fp)
{
  fprintf (fp, "  run year index  area    edges clusters      pop    xmean");
  fprintf (fp, "    ymean      rad    slope cluster_size  %%urban\n");
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogStatValHdr
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogStatValHdr (FILE * fp)
{
#if 1
  int i;
  int num_elements;
  num_elements = sizeof (stats_val_t) / sizeof (double);
  fprintf (fp, "  run year index");
  for (i = 0; i < num_elements; i++)
  {
    fprintf (fp, "%8s ", stats_val_t_names[i]);
  }
  fprintf (fp, "\n");
#else
  fprintf (fp, "\n");
  fprintf (fp, "  run year index  area    edges clusters      pop    xmean");
  fprintf (fp, "    ymean      rad    slope cluster_size  sng      sdg");
  fprintf (fp, "      sdc       og       rt      pop\n");
#endif
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogStatInfo
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogStatInfo (int run, int year, int index,
                     stats_info * stats_ptr, FILE * fp)
{

  fprintf (fp, "%5u %4u %2u ", run, year, index);
  fprintf (fp, "%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
           stats_ptr->area,
           stats_ptr->edges,
           stats_ptr->clusters,
           stats_ptr->pop,
           stats_ptr->xmean,
           stats_ptr->ymean,
           stats_ptr->rad,
           stats_ptr->average_slope,
           stats_ptr->mean_cluster_size,
           stats_ptr->percent_urban
    );
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogStatVal
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogStatVal (int run, int year, int index,
                    stats_val_t * stats_ptr, FILE * fp)
{
  int i;
  int num_elements;
  double *ptr;

  /*num_elements = sizeof(struct stats_val_t)/sizeof(double); */
  num_elements = sizeof (stats_val_t) / sizeof (double);
  ptr = (double *) stats_ptr;

  fprintf (fp, "%5u %4u %2u   ", run, year, index);
#if 1
  for (i = 0; i < num_elements; i++)
  {
    fprintf (fp, "%8.2f ", *ptr);
    ptr++;
  }
  fprintf (fp, "\n");
#else
  fprintf (fp, "%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f ",
           stats_ptr->area,
           stats_ptr->edges,
           stats_ptr->clusters,
           stats_ptr->pop,
           stats_ptr->xmean,
           stats_ptr->ymean,
           stats_ptr->rad,
           stats_ptr->slope
    );
  fprintf (fp, "%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n",
           stats_ptr->mean_cluster_size,
           stats_ptr->sng,
           stats_ptr->sdg,
           stats_ptr->sdc,
           stats_ptr->og,
           stats_ptr->rt,
           stats_ptr->pop
    );
#endif
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogAverages
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_LogAverages (int index, FILE * fp)
{
  LOG_INT (fp, index);
  LOG_FLOAT (fp, average[index].area);
  LOG_FLOAT (fp, average[index].edges);
  LOG_FLOAT (fp, average[index].clusters);
  LOG_FLOAT (fp, average[index].pop);
  LOG_FLOAT (fp, average[index].xmean);
  LOG_FLOAT (fp, average[index].ymean);
  LOG_FLOAT (fp, average[index].rad);
  LOG_FLOAT (fp, average[index].slope);
  LOG_FLOAT (fp, average[index].mean_cluster_size);
  LOG_FLOAT (fp, average[index].sng);
  LOG_FLOAT (fp, average[index].sdg);
  LOG_FLOAT (fp, average[index].sdc);
  LOG_FLOAT (fp, average[index].og);
  LOG_FLOAT (fp, average[index].rt);
  LOG_FLOAT (fp, average[index].pop);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogThisYearStats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogThisYearStats (FILE * fp)
{
  LOG_FLOAT (fp, record.this_year.area);
  LOG_FLOAT (fp, record.this_year.edges);
  LOG_FLOAT (fp, record.this_year.clusters);
  LOG_FLOAT (fp, record.this_year.pop);
  LOG_FLOAT (fp, record.this_year.xmean);
  LOG_FLOAT (fp, record.this_year.ymean);
  LOG_FLOAT (fp, record.this_year.rad);
  LOG_FLOAT (fp, record.this_year.slope);
  LOG_FLOAT (fp, record.this_year.mean_cluster_size);
  LOG_FLOAT (fp, record.this_year.sng);
  LOG_FLOAT (fp, record.this_year.sdg);
  LOG_FLOAT (fp, record.this_year.sdc);
  LOG_FLOAT (fp, record.this_year.og);
  LOG_FLOAT (fp, record.this_year.rt);
  LOG_FLOAT (fp, record.this_year.pop);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogRecord
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_LogRecord (FILE * fp)
{
  LOG_INT (fp, record.run);
  LOG_INT (fp, record.monte_carlo);
  LOG_INT (fp, record.year);
  stats_LogThisYearStats (fp);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_compute_leesalee
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_compute_leesalee (GRID_P Z,                          /* IN     */
                          GRID_P urban,                      /* IN     */
                          double *leesalee)                /* OUT    */
{
  char func[] = "stats_compute_leesalee";
  int i;
  int the_union;
  int intersection;

  FUNC_INIT;
  assert (Z != NULL);
  assert (urban != NULL);
  assert (leesalee != NULL);

  the_union = 0;
  intersection = 0;
  for (i = 0; i < mem_GetTotalPixels (); i++)
  {
    if ((Z[i] != 0) || (urban[i] != 0))
    {
      the_union++;
    }

    if ((Z[i] != 0) && (urban[i] != 0))
    {
      intersection++;
    }
  }

  *leesalee = (double) intersection / the_union;
  FUNC_END;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_Analysis
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_Analysis (double fmatch)
{
#if 1
  char func[] = "stats_Analysis";
#endif
  char std_filename[MAX_FILENAME_LEN];
  char avg_filename[MAX_FILENAME_LEN];
  char cntrl_filename[MAX_FILENAME_LEN];
  char *output_dir;
  int yr;
  int i;
  int run;
  static int avg_log_created = 0;
  static int std_dev_log_created = 0;
  static int control_stats_log_created = 0;

  output_dir = scen_GetOutputDir ();
  run = proc_GetCurrentRun ();

  if (scen_GetWriteAvgFileFlag ())
  {
    sprintf (avg_filename, "%savg_pe_%u.log", output_dir, glb_mype);
    if (!avg_log_created)
    {
      stats_CreateStatsValFile (avg_filename);
      avg_log_created = 1;
    }
  }

  if (scen_GetWriteStdDevFileFlag ())
  {
    sprintf (std_filename, "%sstd_dev_pe_%u.log", output_dir, glb_mype);
    if (!std_dev_log_created)
    {
      stats_CreateStatsValFile (std_filename);
      std_dev_log_created = 1;
    }
  }

  if (proc_GetProcessingType () != PREDICTING)
  {
    sprintf (cntrl_filename, "%scontrol_stats_pe_%u.log", output_dir, glb_mype);
    if (!control_stats_log_created)
    {
      stats_CreateControlFile (cntrl_filename);
      control_stats_log_created = 1;
    }
  }

  if (proc_GetProcessingType () != PREDICTING)
  {
    /*
     *
     * start at i = 1, i = 0 is the initial seed
     *
     */
    for (i = 1; i < igrid_GetUrbanCount (); i++)
    {
      yr = igrid_GetUrbanYear (i);
      stats_CalAverages (i);
      stats_ProcessGrowLog (run, yr);

      if (scen_GetWriteAvgFileFlag ())
      {
        stats_WriteStatsValLine (avg_filename, run, yr, average, i);
      }
      if (scen_GetWriteStdDevFileFlag ())
      {
        stats_WriteStatsValLine (std_filename, run, yr, std_dev, i);
      }
    }
    stats_DoRegressions ();
    stats_DoAggregate (fmatch);
    stats_WriteControlStats (cntrl_filename);
  }
  if (proc_GetProcessingType () == PREDICTING)
  {
    for (yr = scen_GetPredictionStartDate () + 1;
         yr <= proc_GetStopYear (); yr++)
    {
#if 1
      stats_ClearStatsValArrays ();
#endif
      stats_ProcessGrowLog (run, yr);
      if (scen_GetWriteAvgFileFlag ())
      {
        stats_WriteStatsValLine (avg_filename, run, yr, average, 0);
      }
      if (scen_GetWriteStdDevFileFlag ())
      {
        stats_WriteStatsValLine (std_filename, run, yr, std_dev, 0);
      }
#if 1
      stats_ClearStatsValArrays ();
#endif
    }
  }
  stats_ClearStatsValArrays ();
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_Dump
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_Dump (char *file, int line)
{
  int i;
  int yr;

  fprintf (stdout, "%s %u stats_Dump\n", file, line);
  stats_LogStatValHdr (stdout);
  fprintf (stdout, "this_year:\n");
  stats_LogStatVal (proc_GetCurrentRun (), proc_GetCurrentYear (),
                    0, &record.this_year, stdout);
  fprintf (stdout, "running_total:\n");
  for (i = 0; i < MAX_URBAN_YEARS; i++)
  {
    yr = igrid_GetUrbanYear (i);
    if (i == 0)
    {
      yr = 0;
    }
    stats_LogStatVal (proc_GetCurrentRun (), yr, i,
                      &running_total[i], stdout);
  }
  fprintf (stdout, "average:\n");
  for (i = 0; i < MAX_URBAN_YEARS; i++)
  {
    yr = igrid_GetUrbanYear (i);
    if (i == 0)
    {
      yr = 0;
    }
    stats_LogStatVal (proc_GetCurrentRun (), yr, i, &average[i], stdout);
  }
  fprintf (stdout, "std_dev:\n");
  for (i = 0; i < MAX_URBAN_YEARS; i++)
  {
    yr = igrid_GetUrbanYear (i);
    if (i == 0)
    {
      yr = 0;
    }
    stats_LogStatVal (proc_GetCurrentRun (), yr, i, &std_dev[i], stdout);
  }
  stats_LogStatInfoHdr (stdout);
  fprintf (stdout, "stats_actual:\n");
  for (i = 0; i < MAX_URBAN_YEARS; i++)
  {
    yr = igrid_GetUrbanYear (i);
    if (i == 0)
    {
      yr = 0;
    }
    stats_LogStatInfo (proc_GetCurrentRun (), yr, i,
                       &stats_actual[i], stdout);
  }
  fprintf (stdout, "regression:\n");
  stats_LogStatInfo (proc_GetCurrentRun (), 0, 0, &regression, stdout);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogControlStatsHdr
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogControlStatsHdr (FILE * fp)
{
  fprintf (fp, "                                               Cluster\n");
  fprintf (fp, "  Run  Product Compare     Pop   Edges Clusters   ");
  fprintf (fp, "Size Leesalee  Slope ");
  fprintf (fp, " %%Urban   Xmean   Ymean     Rad  Fmatch ");
  fprintf (fp, "Diff  Brd Sprd  Slp   RG\n");
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogControlStats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_LogControlStats (FILE * fp)
{
  fprintf (fp, "%5u %8.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f ",
           proc_GetCurrentRun (),
           aggregate.product,
           aggregate.compare,
           regression.pop,
           regression.edges,
           regression.clusters,
           regression.mean_cluster_size,
           aggregate.leesalee,
           regression.average_slope,
           regression.percent_urban);
  fprintf (fp, "%7.5f %7.5f %7.5f %7.5f %4.0f %4.0f %4.0f %4.0f %4.0f\n",
           regression.xmean,
           regression.ymean,
           regression.rad,
           aggregate.fmatch,
           coeff_GetSavedDiffusion (),
           coeff_GetSavedBreed (),
           coeff_GetSavedSpread (),
           coeff_GetSavedSlopeResist (),
           coeff_GetSavedRoadGravity ());
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_DoAggregate
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_DoAggregate (double fmatch)
{
  char func[] = "stats_DoAggregate";
  int last_index;
  int i;
  double fmatch_tmp = 1.0;
  double numerator;
  double denominator;

  last_index = igrid_GetUrbanCount () - 1;
  aggregate.fmatch = fmatch;
  aggregate.actual = stats_actual[last_index].pop;
  aggregate.simulated = average[last_index].pop;
  aggregate.leesalee = 0.0;
  for (i = 1; i < igrid_GetUrbanCount (); i++)
  {
    aggregate.leesalee += average[i].leesalee;
  }
  aggregate.leesalee /= (igrid_GetUrbanCount () - 1);
  if (aggregate.actual > aggregate.simulated)
  {
    if (aggregate.actual != 0.0)
    {
      denominator = aggregate.actual;
      numerator = aggregate.simulated;
      aggregate.compare = numerator / denominator;
    }
    else
    {
      sprintf (msg_buf, "aggregate.actual = 0.0");
      LOG_ERROR (msg_buf);
      EXIT (1);
    }
  }
  else
  {
    if (aggregate.simulated != 0.0)
    {
      denominator = aggregate.simulated;
      numerator = aggregate.actual;
      aggregate.compare = numerator / denominator;
    }
    else
    {
      sprintf (msg_buf, "aggregate.simulated = 0.0");
      LOG_ERROR (msg_buf);
      EXIT (1);
    }
  }
  if (scen_GetDoingLanduseFlag ())
  {
    fmatch_tmp = fmatch;
  }
  aggregate.product =
    aggregate.compare *
    aggregate.leesalee *
    regression.edges *
    regression.clusters *
    regression.pop *
    regression.xmean *
    regression.ymean *
    regression.rad *
    regression.average_slope *
    regression.mean_cluster_size *
    regression.percent_urban *
    fmatch_tmp;

}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_DoRegressions
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_DoRegressions ()
{
  double dependent[MAX_URBAN_YEARS];
  double independent[MAX_URBAN_YEARS];
  int nobs;
  int i;

  /*VerD*/
  double DDD,RRR;

  if (proc_GetProcessingType () != PREDICTING)
  {

	if (WriteSlopeFileFlag == 1)
	{
  		if (proc_GetCurrentRun()  == 0 )
		{
			fprintf(fpVerD3,	"run ");
			fprintf(fpVerD3,	"diffusion ");
			fprintf(fpVerD3,	"breed ");
			fprintf(fpVerD3,	"spread ");
			fprintf(fpVerD3,	"slp_resst ");
			fprintf(fpVerD3,	"road_grav ");
			fprintf(fpVerD3,	"area area area ");
			fprintf(fpVerD3,	"edges edges edges ");
			fprintf(fpVerD3,	"cluster cluster cluster ");
			fprintf(fpVerD3,	"pop pop pop ");
			fprintf(fpVerD3,	"xmean xmean xmean ");
			fprintf(fpVerD3,	"ymean ymean ymean ");
			fprintf(fpVerD3,	"radius radius radius ");
			fprintf(fpVerD3,	"avg_slope avg_slope avg_slope ");
			fprintf(fpVerD3,	"mn_cl_sz mn_cl_sz mn_cl_sz ");
			fprintf(fpVerD3,	"pct_urban pct_urban pct_urban ");
			fprintf(fpVerD3,	"\n");

			fprintf(fpVerD3,	"number ");
			fprintf(fpVerD3,	"coeff ");
			fprintf(fpVerD3,	"coeff ");
			fprintf(fpVerD3,	"coeff ");
			fprintf(fpVerD3,	"coeff ");
			fprintf(fpVerD3,	"coeff ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"r slope intrcpt ");
			fprintf(fpVerD3,	"\n");
		}
		fprintf (fpVerD3,"%d %10.4f %10.4f %10.4f %10.4f %10.4f ",
					proc_GetCurrentRun (),
					coeff_GetSavedDiffusion(),
					coeff_GetSavedBreed(),
					 coeff_GetSavedSpread(),
					coeff_GetSavedSlopeResist(),
					coeff_GetSavedRoadGravity() );
	}
  }
  /*VerD*/

  nobs = igrid_GetUrbanCount () - 1;
  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].area;
    independent[i - 1] = average[i].area;
  }
  regression.area = stats_linefit (dependent, independent, nobs);

  nobs = igrid_GetUrbanCount () - 1;
  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].edges;
    independent[i - 1] = average[i].edges;
  }
  regression.edges = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].clusters;
    independent[i - 1] = average[i].clusters;
  }
  regression.clusters = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].pop;
    independent[i - 1] = average[i].pop;
  }
  regression.pop = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].xmean;
    independent[i - 1] = average[i].xmean;
  }
  regression.xmean = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].ymean;
    independent[i - 1] = average[i].ymean;
  }
  regression.ymean = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].rad;
    independent[i - 1] = average[i].rad;
  }
  regression.rad = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].average_slope;
    independent[i - 1] = average[i].slope;
  }
  regression.average_slope = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].mean_cluster_size;
    independent[i - 1] = average[i].mean_cluster_size;
  }
  regression.mean_cluster_size = stats_linefit (dependent, independent, nobs);

  for (i = 1; i <= nobs; i++)
  {
    dependent[i - 1] = stats_actual[i].percent_urban;
    independent[i - 1] = average[i].percent_urban;
  }
  regression.percent_urban = stats_linefit (dependent, independent, nobs);

  /*VerD*/

  if (WriteSlopeFileFlag == 1)
  {
	if (proc_GetProcessingType () != PREDICTING)   fprintf(fpVerD3,	"\n");
  }

	if (WriteRatioFileFlag == 1)
	{
		if (proc_GetProcessingType () != PREDICTING)
		{
			if (proc_GetCurrentRun()  == 0 )
			{
				fprintf(fpVerD4,	"run ");
				fprintf(fpVerD4,	"diffusion ");
				fprintf(fpVerD4,	"breed ");
				fprintf(fpVerD4,	"spread ");
				fprintf(fpVerD4,	"slp_resst ");
				fprintf(fpVerD4,	"road_grav ");
				fprintf(fpVerD4,	"control ");
				fprintf(fpVerD4,	"area area area ");
				fprintf(fpVerD4,	"edges edges edges ");
				fprintf(fpVerD4,	"cluster cluster cluster ");
				fprintf(fpVerD4,	"pop pop pop ");
				fprintf(fpVerD4,	"xmean xmean xmean ");
				fprintf(fpVerD4,	"ymean ymean ymean ");
				fprintf(fpVerD4,	"radius radius radius ");
				fprintf(fpVerD4,	"avg_slope avg_slope avg_slope ");
				fprintf(fpVerD4,	"mn_cl_sz mn_cl_sz mn_cl_sz ");
				fprintf(fpVerD4,	"pct_urban pct_urban pct_urban ");
				fprintf(fpVerD4,	"\n");

				fprintf(fpVerD4,	"number ");
				fprintf(fpVerD4,	"coeff ");
				fprintf(fpVerD4,	"coeff ");
				fprintf(fpVerD4,	"coeff ");
				fprintf(fpVerD4,	"coeff ");
				fprintf(fpVerD4,	"coeff ");
				fprintf(fpVerD4,	"year ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"diff ratio fract ");
				fprintf(fpVerD4,	"\n");
			}

			for (i = 1; i <= nobs; i++)
			{
				fprintf (fpVerD4,"%d %10.4f %10.4f %10.4f %10.4f %10.4f ",
							proc_GetCurrentRun (),
							coeff_GetSavedDiffusion(),
							coeff_GetSavedBreed(),
							coeff_GetSavedSpread(),
							coeff_GetSavedSlopeResist(),
							coeff_GetSavedRoadGravity() );

			    /*NOTE:  In the code below, RRR, and RRR-1 are really redundant,
					so maybe we don't need to print both.*/


				DDD = average[i].area - stats_actual[i].area;
				if  (stats_actual[i].area == 0)   RRR = 0;	
				else   RRR = average[i].area / stats_actual[i].area;
				fprintf(fpVerD4,"%d %f %f %f   ",igrid_GetUrbanYear(i), DDD, RRR, RRR-1 );
    
				DDD = average[i].edges - stats_actual[i].edges;	
				if  (stats_actual[i].edges == 0)   RRR = 0;
				else   RRR = average[i].edges / stats_actual[i].edges;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].clusters - stats_actual[i].clusters;
				if  (stats_actual[i].clusters == 0)   RRR = 0;
				else   RRR = average[i].clusters / stats_actual[i].clusters;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].pop - stats_actual[i].pop;
				if  (stats_actual[i].pop == 0)   RRR = 0;
				else   RRR = average[i].pop / stats_actual[i].pop;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].xmean - stats_actual[i].xmean;
				if  (stats_actual[i].xmean == 0)   RRR = 0;
				else   RRR = average[i].xmean / stats_actual[i].xmean;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].ymean - stats_actual[i].ymean;
				if  (stats_actual[i].ymean == 0)   RRR = 0;
				else   RRR = average[i].ymean / stats_actual[i].ymean;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].rad - stats_actual[i].rad;
				if  (stats_actual[i].rad == 0)   RRR = 0;
				else   RRR = average[i].rad / stats_actual[i].rad;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].slope - stats_actual[i].average_slope;
				if  (stats_actual[i].average_slope == 0)   RRR = 0;
				else   RRR = average[i].slope / stats_actual[i].average_slope;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].mean_cluster_size - stats_actual[i].mean_cluster_size;
				if  (stats_actual[i].mean_cluster_size == 0)   RRR = 0;
				else   RRR = average[i].mean_cluster_size / stats_actual[i].mean_cluster_size;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );
    
				DDD = average[i].percent_urban - stats_actual[i].percent_urban;
				if  (stats_actual[i].percent_urban == 0)   RRR = 0;
				else   RRR = average[i].percent_urban / stats_actual[i].percent_urban;
				fprintf(fpVerD4,"%f %f %f   ", DDD, RRR, RRR-1 );

				fprintf(fpVerD4,	"\n");
			}
		}

	/*VerD*/

	}

}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_Save
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_Save (char *filename)
{
  char func[] = "stats_Save";
  int num_written;
  int fseek_loc;
  int index;
  int i;
  FILE *fp;
  record.run = proc_GetCurrentRun ();
  record.monte_carlo = proc_GetCurrentMonteCarlo ();
  record.year = proc_GetCurrentYear ();
  index = 0;
  if (proc_GetProcessingType () != PREDICTING)
  {
    index = igrid_UrbanYear2Index (record.year);
  }

  stats_UpdateRunningTotal (index);

  if (record.monte_carlo == 0)
  {
    FILE_OPEN (fp, filename, "wb");
    for (i = 0; i < scen_GetMonteCarloIterations (); i++)
    {
      num_written = fwrite (&record, sizeof (record), 1, fp);
      if (num_written != 1)
      {
        printf ("%s %u ERROR\n", __FILE__, __LINE__);
      }
    }

  }
  else
  {
    FILE_OPEN (fp, filename, "r+b");
    rewind (fp);
    fseek_loc = fseek (fp, sizeof (record) * record.monte_carlo, SEEK_SET);
    num_written = fwrite (&record, sizeof (record), 1, fp);
    if (num_written != 1)
    {
      printf ("%s %u ERROR\n", __FILE__, __LINE__);
    }

  }
  fclose (fp);

}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ProcessGrowLog
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_ProcessGrowLog (int run, int year)
{
  char func[] = "stats_ProcessGrowLog";
  FILE *fp;
  int index;
  int fseek_loc;
  int i;
  int mc_count = 0;
  char filename[MAX_FILENAME_LEN];
  char command[MAX_FILENAME_LEN + 3];

  sprintf (filename, "%sgrow_%u_%u.log", scen_GetOutputDir (), run, year);
  sprintf (command, "rm %s", filename);

  FILE_OPEN (fp, filename, "rb");

  if (proc_GetProcessingType () != PREDICTING)
  {
    while (fread (&record, sizeof (record), 1, fp))
    {
      if (mc_count >= scen_GetMonteCarloIterations ())
      {
        sprintf (msg_buf, "mc_count >= scen_GetMonteCarloIterations ()");
        LOG_ERROR (msg_buf);
        EXIT (1);
      }
      if (feof (fp) || ferror (fp))
      {
        sprintf (msg_buf, "feof (fp) || ferror (fp)");
        LOG_ERROR (msg_buf);
        EXIT (1);
      }
      index = igrid_UrbanYear2Index (year);
      stats_CalStdDev (index);
      mc_count++;
    }
  }
  else
  {
    while (fread (&record, sizeof (record), 1, fp))
    {
      if (mc_count >= scen_GetMonteCarloIterations ())
      {
        sprintf (msg_buf, "mc_count >= scen_GetMonteCarloIterations ()");
        LOG_ERROR (msg_buf);
        EXIT (1);
      }
      if (feof (fp) || ferror (fp))
      {
        sprintf (msg_buf, "feof (fp) || ferror (fp)");
        LOG_ERROR (msg_buf);
        EXIT (1);
      }
      stats_UpdateRunningTotal (0);
    }
    stats_CalAverages (0);
    rewind (fp);
    mc_count = 0;
    while (fread (&record, sizeof (record), 1, fp))
    {
      if (mc_count >= scen_GetMonteCarloIterations ())
      {
        sprintf (msg_buf, "mc_count >= scen_GetMonteCarloIterations ()");
        LOG_ERROR (msg_buf);
        sprintf (msg_buf, "mc_count= %u scen_GetMonteCarloIterations= %u",
                 mc_count, scen_GetMonteCarloIterations ());
        LOG_ERROR (msg_buf);
        EXIT (1);
      }
      if (feof (fp) || ferror (fp))
      {
        sprintf (msg_buf, "feof (fp) || ferror (fp)");
        LOG_ERROR (msg_buf);
        EXIT (1);
      }
      stats_CalStdDev (0);
      mc_count++;
    }
  }
  fclose (fp);
  system (command);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_linefit
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static double
  stats_linefit (double *dependent,
                 double *independent,
                 int number_of_observations)
{
  char func[] = "Linefit";
  double dependent_avg;
  double independent_avg;
  double cross;
  double sum_dependent;
  double sum_independent;
  double r;
  int n;
  
  /*VerD*/
  double yslope;
  double yintercept;
  /*VerD*/

  FUNC_INIT;
  assert (dependent != NULL);
  assert (independent != NULL);
  assert (number_of_observations > 0);

  dependent_avg = 0;
  independent_avg = 0;

  for (n = 0; n < number_of_observations; n++)

  {
    dependent_avg += dependent[n];
    independent_avg += independent[n];
  }

  if (number_of_observations > 0)
  {
    dependent_avg /= (double) number_of_observations;
    independent_avg /= (double) number_of_observations;
  }
  else
  {
    sprintf (msg_buf, "number_of_observations = %d", number_of_observations);
    LOG_ERROR (msg_buf);
    EXIT (1);
  }
  cross = 0;
  sum_dependent = 0;
  sum_independent = 0;

  for (n = 0; n < number_of_observations; n++)
  {
    cross += ((dependent[n] - dependent_avg) * (independent[n] -
                                                independent_avg));
    sum_dependent += ((dependent[n] - dependent_avg) * (dependent[n] -
                                                        dependent_avg));
    sum_independent += ((independent[n] - independent_avg) * (independent[n]
                                                     - independent_avg));
  }

  /*VerD*/
  /*NOTE:  The following calculations assume that the desired least-squares
  fit minimizes the sum of squares of the VERTICAL distances from the straight
  line  y = mx + b, where m is the "yslope" and b is the "yintercept".  In this sense,
  the "dependent" variable y is the simulated quantity (e.g. modeled area count) and the
  "independent" variable x is the corresponding actual quantity (e.g. real area count).
  Thus, the terms "dependent" and "independent" are the OPPOSITE from those used in the
  calling program!  (The quantity r in the program is still correct because it is
  symetric in x ane y.)*/
  /*VerD*/
           
  if (sum_dependent == 0.0)
  {
      yslope = 0.0;
      yintercept = 0.0;
  }
  else
  {
      yslope = cross/sum_dependent;
      yintercept = independent_avg - yslope*dependent_avg;
  }
  /*VerD*/

  r = 0;

  if (sum_dependent * sum_independent < 1e-11)
    r = 0;
  else
    r = cross / pow (sum_dependent * sum_independent, 0.5);

  /*VerD*/
  if (WriteSlopeFileFlag == 1)
  {
	if (proc_GetProcessingType () != PREDICTING)
      fprintf(fpVerD3,"%f %f %f ", r, yslope, yintercept );
  }
  /*VerD*/

  FUNC_END;
  return (r * r);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_ComputeBaseStats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

static void
  stats_ComputeBaseStats ()
{
  char func[] = "stats_ComputeBaseStats";
  int i;
  int total_pixels;
  GRID_P urban_ptr;
  GRID_P slope_ptr;
  GRID_P stats_workspace1;
  GRID_P stats_workspace2;
  int road_pixel_count;
  int excluded_pixel_count;

  total_pixels = mem_GetTotalPixels ();
  assert (total_pixels > 0);

  for (i = 0; i < igrid_GetUrbanCount (); i++)
  {
    urban_ptr = igrid_GetUrbanGridPtr (__FILE__, func, __LINE__, i);
    slope_ptr = igrid_GetSlopeGridPtr (__FILE__, func, __LINE__);
    stats_workspace1 = mem_GetWGridPtr (__FILE__, func, __LINE__);
    stats_workspace2 = mem_GetWGridPtr (__FILE__, func, __LINE__);

    stats_compute_stats (urban_ptr,                          /* IN     */
                         slope_ptr,                          /* IN     */
                         &stats_actual[i].area,              /* OUT    */
                         &stats_actual[i].edges,             /* OUT    */
                         &stats_actual[i].clusters,          /* OUT    */
                         &stats_actual[i].pop,               /* OUT    */
                         &stats_actual[i].xmean,             /* OUT    */
                         &stats_actual[i].ymean,             /* OUT    */
                         &stats_actual[i].average_slope,     /* OUT    */
                         &stats_actual[i].rad,               /* OUT    */
                         &stats_actual[i].mean_cluster_size,   /* OUT    */
                         stats_workspace1,                   /* MOD    */
                         stats_workspace2);                /* MOD    */

    road_pixel_count = igrid_GetIGridRoadPixelCount (proc_GetCurrentYear ());
    excluded_pixel_count = igrid_GetIGridExcludedPixelCount ();
	
    /*VerD*/
	/*Note:  An extra factor of 100.0 was removed from the formula below.
	   This extra factor did not adversely affect any of the output statistics.
	   It did however cause a scaling error in some output data in the new
	   ratio/slope files.*/
    /*VerD*/

	stats_actual[i].percent_urban = 100.0 *
			(stats_actual[i].pop + road_pixel_count) /
			(igrid_GetNumRows () * igrid_GetNumCols () - road_pixel_count -
			excluded_pixel_count);
    /*VerD*/


    urban_ptr = igrid_GridRelease (__FILE__, func, __LINE__, urban_ptr);
    slope_ptr = igrid_GridRelease (__FILE__, func, __LINE__, slope_ptr);
    stats_workspace1 = mem_GetWGridFree (__FILE__, func, __LINE__,
                                         stats_workspace1);
    stats_workspace2 = mem_GetWGridFree (__FILE__, func, __LINE__,
                                         stats_workspace2);

  }
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_compute_stats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_compute_stats (GRID_P Z,                             /* IN     */
                       GRID_P slp,                           /* IN     */
                       double *stats_area,                   /* OUT    */
                       double *stats_edges,                  /* OUT    */
                       double *stats_clusters,               /* OUT    */
                       double *stats_pop,                    /* OUT    */
                       double *stats_xmean,                  /* OUT    */
                       double *stats_ymean,                  /* OUT    */
                       double *stats_average_slope,          /* OUT    */
                       double *stats_rad,                    /* OUT    */
                       double *stats_mean_cluster_size,      /* OUT    */
                       GRID_P scratch_gif1,                  /* MOD    */
                       GRID_P scratch_gif2)                /* MOD    */

{
  char func[] = "stats_compute_stats";

  FUNC_INIT;
  assert (Z != NULL);
  assert (slp != NULL);
  assert (stats_area != NULL);
  assert (stats_edges != NULL);
  assert (stats_clusters != NULL);
  assert (stats_pop != NULL);
  assert (stats_xmean != NULL);
  assert (stats_ymean != NULL);
  assert (stats_average_slope != NULL);
  assert (stats_rad != NULL);
  assert (stats_mean_cluster_size != NULL);
  assert (scratch_gif1 != NULL);
  assert (scratch_gif2 != NULL);

/** D. Donato 8/21/2006 Added to deal with cumulative growth                 **/
  zgrwth_row = (short *)  mem_GetGRZrowptr();
  zgrwth_col = (short *)  mem_GetGRZcolptr();
  zgrwth_count =          mem_GetGRZcount();
/******************* 8/21/2006  ***********************************************/

  /*
   *
   * compute the number of edge pixels
   *
   */
  stats_edge (Z,                                             /* IN     */
              stats_area,                                    /* OUT    */
              stats_edges);                                /* OUT    */


  /*
   *
   * compute the number of clusters
   *
   */
  stats_cluster (Z,                                          /* IN     */
                 stats_clusters,                             /* OUT    */
                 stats_pop,                                  /* OUT    */
                 stats_mean_cluster_size,                    /* OUT    */
                 scratch_gif1,                               /* MOD    */
                 scratch_gif2);                            /* MOD    */

  /*
   *
   * compute means
   *
   */
  stats_circle (Z,                                           /* IN     */
                slp,                                         /* IN     */
                *stats_area,                                 /* IN     */
                stats_xmean,                                 /* OUT    */
                stats_ymean,                                 /* OUT    */
                stats_average_slope,                         /* OUT    */
                stats_rad);                                /* OUT    */

  FUNC_END;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_edge
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_edge (GRID_P Z,                                      /* IN     */
              double *stats_area,                            /* OUT    */
              double *stats_edges)                         /* OUT    */
{
  char func[] = "stats_edge";
/* D.D. 8/21/2006 Variable "k" added. */
  int i, k;
  int j;
  int edge;
  int edges;
  int area;
  int rowi[4] = {-1, 1, 0, 0};
  int colj[4] = {0, 0, -1, 1};
  int loop;
  int row;
  int col;
  int nrows;
  int ncols;

  FUNC_INIT;
  assert (stats_area != NULL);
  assert (stats_edges != NULL);
  assert (Z != NULL);
  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);

  edges = 0;
  area = 0;

/*** D.Donato 8/21/2006 Alternative loop provided for Z grid.  ***
     The original two-level nested loop over i and j is now
     supplemented by an alternative loop. The original loop
     will be executed during start-up computation of base
     statistics using each urban grid in sequence; an the
     alternative loop will be used during growth computation
     using the cumulative-growth Z grid. The alternative
     loop is faster than the original.

*********************     8/21/2006    **************************/

 if (Z == mem_GetGRZpointer())
 {
    for ( k = 0; k < zgrwth_count; k++)
    {
      i = zgrwth_row[k];
      j = zgrwth_col[k];

      edge = FALSE;

      if (Z[OFFSET (i, j)] != 0)
      {
        area++;

        /* this does a 4 neighbor search (N, S, E, W) */
        for (loop = 0; loop <= 3; loop++)
        {
          row = i + rowi[loop];
          col = j + colj[loop];

          if (IMAGE_PT (row, col))
          {
            if (Z[OFFSET (row, col)] == 0)
            {
              edge = TRUE;
            }
          }
        }

        if (edge)
        {
          edges++;
        }
      }
    }
 }

 else

 {
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      edge = FALSE;

      if (Z[OFFSET (i, j)] != 0)
      {
        area++;

        /* this does a 4 neighbor search (N, S, E, W) */
        for (loop = 0; loop <= 3; loop++)
        {
          row = i + rowi[loop];
          col = j + colj[loop];

          if (IMAGE_PT (row, col))
          {
            if (Z[OFFSET (row, col)] == 0)
            {
              edge = TRUE;
            }
          }
        }

        if (edge)
        {
          edges++;
        }
      }
    }
  }
 }

  *stats_area = area;
  *stats_edges = edges;
  FUNC_END;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_circle
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_circle (GRID_P Z,                                    /* IN     */
                GRID_P slp,                                  /* IN     */
                int stats_area,                              /* IN     */
                double *stats_xmean,                         /* OUT    */
                double *stats_ymean,                         /* OUT    */
                double *stats_average_slope,                 /* OUT    */
                double *stats_rad)                         /* OUT    */
{
  char func[] = "stats_circle";
/* D.D. 8/21/2006 Variable "k" added. */
  int i, k;
  int j;
  int number;
  double xmean;
  double ymean;
  double addslope;
  int nrows;
  int ncols;

  FUNC_INIT;
  assert (stats_xmean != NULL);
  assert (stats_ymean != NULL);
  assert (stats_average_slope != NULL);
  assert (stats_rad != NULL);
  assert (Z != NULL);
  assert (slp != NULL);

  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);
  addslope = 0.0;
  ymean = 0.0;
  xmean = 0.0;
  number = 0.0;

  /*
   *
   * first, compute the means
   *
   */


/*** D.Donato 8/21/2006 Alternative loop provided for Z grid.  ***
     The original two-level nested loop over i and j is now
     supplemented by an alternative loop. The original loop
     will be executed during start-up computation of base
     statistics using each urban grid in sequence; an the
     alternative loop will be used during growth computation
     using the cumulative-growth Z grid. The alternative
     loop is faster than the original.

*********************     8/21/2006    **************************/


 if (Z == mem_GetGRZpointer())
 {

    for ( k = 0; k < zgrwth_count; k++)
    {
      i = zgrwth_row[k];
      j = zgrwth_col[k];

      if (Z[OFFSET (i, j)] > 0)
      {
        addslope += slp[OFFSET (i, j)];
        xmean += (double) j;
        ymean += (double) i;
        number++;
      }
    }
 }

 else

 {
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      if (Z[OFFSET (i, j)] > 0)
      {
        addslope += slp[OFFSET (i, j)];
        xmean += (double) j;
        ymean += (double) i;
        number++;
      }
    }
  }
 }


  if (number <= 0)
  {
    sprintf (msg_buf, "number = %d", number);
    LOG_ERROR (msg_buf);
    EXIT (1);
  }
  xmean /= (double) number;
  ymean /= (double) number;
  *stats_xmean = xmean;
  *stats_ymean = ymean;
  *stats_average_slope = addslope / number;

  /*
   *
   * compute the radius of the circle with same area as number
   *
   */
  *stats_rad = pow ((stats_area / PI), 0.5);

  FUNC_END;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_MemoryLog
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_MemoryLog (FILE * fp)
{
  LOG_MEM (fp, &cir_q[0], sizeof (ugm_link), SIZE_CIR_Q);
  LOG_MEM (fp, &stats_actual[0], sizeof (stats_info), MAX_URBAN_YEARS);
  LOG_MEM (fp, &regression, sizeof (stats_info), 1);
  LOG_MEM (fp, &average[0], sizeof (stats_val_t), MAX_URBAN_YEARS);
  LOG_MEM (fp, &std_dev[0], sizeof (stats_val_t), MAX_URBAN_YEARS);
  LOG_MEM (fp, &running_total[0], sizeof (stats_val_t), MAX_URBAN_YEARS);
  LOG_MEM (fp, &urbanization_attempt, sizeof (urbanization_attempt), 1);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_cluster
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_cluster (GRID_P Z,                                   /* IN     */
                 double *stats_clusters,                     /* OUT    */
                 double *stats_pop,                          /* OUT    */
                 double *stats_mean_cluster_size,            /* OUT    */
                 GRID_P scratch_gif1,                        /* MOD    */
                 GRID_P scratch_gif2)                      /* MOD    */
{
  char func[] = "stats_cluster";
/* D.D. 8/21/2006 Variable "k" added. */
  int i, k;
  int j;
  int depth;
  int num_clusters;
  int sum;
  int row;
  int col;
  int rrow;
  int ccol;
  int rowi[4] = {1, 0, -1, 0};
  int colj[4] = {0, 1, 0, -1};
  int loop;
/**^^**^^***
  long *visited;
  long *clusters;
***^^**^^**/
  PIXEL *visited;
  PIXEL *clusters;
  int total_pixels;
  int nrows;
  int ncols;

  FUNC_INIT;
  assert (stats_clusters != NULL);
  assert (stats_pop != NULL);
  assert (stats_mean_cluster_size != NULL);
  assert (Z != NULL);
  assert (scratch_gif1 != NULL);
  assert (scratch_gif2 != NULL);
  total_pixels = mem_GetTotalPixels ();
  assert (total_pixels > 0);
  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);

  sum = 0;
  *stats_pop = 0;
  depth = 0;
  num_clusters = 0;

  visited = scratch_gif1;
  clusters = scratch_gif2;
  for (i = 0; i < total_pixels; i++)
  {
    visited[i] = 0;
/** D.D. 8/21/2006 The next three lines are unnecessary and inefficient. ***
  }
  for (i = 0; i < total_pixels; i++)
  {
*** D.D. 8/21/2006                        *********************************/
    if (Z[i] != 0)
    {
      clusters[i] = 1;
      (*stats_pop)++;
    }
    else
    {
      clusters[i] = 0;
    }
  }
  for (j = 0; j < ncols; j++)
  {
    clusters[OFFSET (0, j)] = 0;
    clusters[OFFSET (nrows - 1, j)] = 0;
  }
  for (i = 0; i < nrows; i++)
  {
    clusters[OFFSET (i, 0)] = 0;
    clusters[OFFSET (i, ncols - 1)] = 0;
  }


/*** D.Donato 8/21/2006 Alternative loop provided for Z grid.  ***
     The original two-level nested loop over i and j is now
     supplemented by an alternative loop. The original loop
     will be executed during start-up computation of base
     statistics using each urban grid in sequence; an the
     alternative loop will be used during growth computation
     using the cumulative-growth Z grid. The alternative
     loop is faster than the original.

*********************     8/21/2006    **************************/

 if (Z == mem_GetGRZpointer())
 {
    for ( k = 0; k < zgrwth_count; k++)
    {
      i = zgrwth_row[k];
      j = zgrwth_col[k];
      
      if (clusters[OFFSET (i, j)] == 1 && visited[OFFSET (i, j)] == 0)
      {
        sum++;
        rrow = i;
        ccol = j;
        visited[OFFSET (i, j)] = 1;
        Q_STORE (rrow, ccol);
        do
        {
          Q_RETREIVE (row, col);
          for (loop = 0; loop <= 3; loop++)
          {
            rrow = row + rowi[loop];
            ccol = col + colj[loop];

            if (IMAGE_PT (rrow, ccol))
            {
              if (clusters[OFFSET (rrow, ccol)] == 1 &&
                  !visited[OFFSET (rrow, ccol)])
              {
                visited[OFFSET (rrow, ccol)] = 1;
                Q_STORE (rrow, ccol);

                sum++;
              }
            }
          }
        }
        while (depth > 0);
        num_clusters++;
      }
    }
 }

 else

 {

  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      if (clusters[OFFSET (i, j)] == 1 && visited[OFFSET (i, j)] == 0)
      {
        sum++;
        rrow = i;
        ccol = j;
        visited[OFFSET (i, j)] = 1;
        Q_STORE (rrow, ccol);
        do
        {
          Q_RETREIVE (row, col);
          for (loop = 0; loop <= 3; loop++)
          {
            rrow = row + rowi[loop];
            ccol = col + colj[loop];

            if (IMAGE_PT (rrow, ccol))
            {
              if (clusters[OFFSET (rrow, ccol)] == 1 &&
                  !visited[OFFSET (rrow, ccol)])
              {
                visited[OFFSET (rrow, ccol)] = 1;
                Q_STORE (rrow, ccol);

                sum++;
              }
            }
          }
        }
        while (depth > 0);
        num_clusters++;
      }
    }
  }

 }

  *stats_clusters = num_clusters;
  if (num_clusters > 0)
  {
    *stats_mean_cluster_size = sum / num_clusters;
  }
  else
  {
    sprintf (msg_buf, "num_clusters=%d", num_clusters);
    LOG_ERROR (msg_buf);
    EXIT (1);
  }
  FUNC_END;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogBaseStats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/

void
  stats_LogBaseStats (FILE * fp)
{
  char func[] = "stats_LogBaseStats";
  int i;
  int count;

  FUNC_INIT;
  count = igrid_GetUrbanCount ();
  assert (count > 0);


  fprintf (fp, "\n\n");
  fprintf (fp, "************************LOG OF BASE STATISTICS");
  fprintf (fp, " FOR URBAN INPUT DATA********************\n");
  fprintf (fp, " Year       Area      Edges   Clusters         ");
  fprintf (fp, "Pop       Mean Center        Radius");
  fprintf (fp, "   Avg Slope  MeanClusterSize\n");
  for (i = 0; i < count; i++)
  {
    fprintf (fp, "%5d   %8.2f   %8.2f   %8.2f    %8.2f  (%8.2f,%8.2f)",
             igrid_GetUrbanYear (i),
             stats_actual[i].area,
             stats_actual[i].edges,
             stats_actual[i].clusters,
             stats_actual[i].pop,
             stats_actual[i].xmean,
             stats_actual[i].ymean);
    fprintf (fp, "   %8.2f  %10.2f      %6.3f\n",
             stats_actual[i].rad,
             stats_actual[i].average_slope,
             stats_actual[i].mean_cluster_size);

    /*VerD*/
	if (proc_GetProcessingType () != PREDICTING)
	{
	    if (WriteXypointsFileFlag == 1)
		{
			fprintf (fpVerD2, "%6d %6d %6d %6d %6d %6d %6d %6d %6.0f\n",
				proc_GetTotalRuns(),
				scen_GetMonteCarloIterations(),
				0, 0, 0, 0, 0,
				igrid_GetUrbanYear(i),
				stats_actual[i].area  );
		}
	}
    /*VerD*/


  }
  fprintf (fp, "\n\n");
  FUNC_END;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_WriteControlStats
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  stats_WriteControlStats (char *filename)
{
  char func[] = "stats_WriteControlStats";
  FILE *fp;

  FILE_OPEN (fp, filename, "a");

  stats_LogControlStats (fp);
  fclose (fp);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_InitUrbanizationAttempts
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_InitUrbanizationAttempts ()
{
  urbanization_attempt.successes = 0;
  urbanization_attempt.z_failure = 0;
  urbanization_attempt.delta_failure = 0;
  urbanization_attempt.slope_failure = 0;
  urbanization_attempt.excluded_failure = 0;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_LogUrbanizationAttempts
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_LogUrbanizationAttempts (FILE * fp)
{
  int total;

  total = urbanization_attempt.successes +
    urbanization_attempt.z_failure +
    urbanization_attempt.delta_failure +
    urbanization_attempt.slope_failure +
    urbanization_attempt.excluded_failure;

  fprintf (fp, "\nLOG OF URBANIZATION ATTEMPTS\n");
  fprintf (fp, "Num Success                = %u\n",
           urbanization_attempt.successes);
  fprintf (fp, "Num Z Type Failures        = %u\n",
           urbanization_attempt.z_failure);
  fprintf (fp, "Num Delta Type Failures    = %u\n",
           urbanization_attempt.delta_failure);
  fprintf (fp, "Num Slope Type Failures    = %u\n",
           urbanization_attempt.slope_failure);
  fprintf (fp, "Num Exlcuded Type Failures = %u\n",
           urbanization_attempt.excluded_failure);
  fprintf (fp, "Total Attempts             = %u\n", total);
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_IncrementUrbanSuccess
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_IncrementUrbanSuccess ()
{
  urbanization_attempt.successes++;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_IncrementZFailure
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_IncrementZFailure ()
{
  urbanization_attempt.z_failure++;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_IncrementDeltaFailure
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_IncrementDeltaFailure ()
{
  urbanization_attempt.delta_failure++;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_IncrementSlopeFailure
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_IncrementSlopeFailure ()
{
  urbanization_attempt.slope_failure++;
}
/******************************************************************************
*******************************************************************************
** FUNCTION NAME: stats_IncrementEcludedFailure
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  stats_IncrementEcludedFailure ()
{
  urbanization_attempt.excluded_failure++;
}
