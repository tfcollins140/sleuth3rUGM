/*******************************************************************************

  MODULE:                   spread.c

  CONTAINING SYSTEM:        SLEUTH-3r (based on SLEUTH Model 3.0 Beta)
                            (Slope, Land-cover, Exclusion, Urbanization, 
                            Transportation, and Hillshade)
                            also known as UGM 3.0 Beta (for Urban Growth Model)

  VERSION:                  SLEUTH-3r [Includes Version D features]

  REVISION DATE:            August 31, 2006a
                            [Annotations added August 19, 2009]

  PURPOSE:

     This module models the spread of urbanization for a particular
     year within a particular iteration (Monte Carlo run) for a
     particular combination of input coefficients.

  NOTES:

  MODIFICATIONS:

  TO DO:

**************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "igrid_obj.h"
#include "landclass_obj.h"
#include "globals.h"
#include "random.h"
#include "utilities.h"
#include "memory_obj.h"
#include "igrid_obj.h"
#include "ugm_macros.h"
#include "coeff_obj.h"
#include "timer_obj.h"
#include "proc_obj.h"
#include "scenario_obj.h"
#include "stats_obj.h"

/*VerD*/
extern float aux_diffusion_coeff;
extern float aux_breed_coeff;
extern float aux_diffusion_mult;

float road_growth_breed_coefficient;
float road_growth_diffusion_coefficient;
/*VerD*/

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                                 MACROS                                    **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
#define SPREAD_MODULE
#define SWGHT_TYPE float
#define SLOPE_WEIGHT_ARRAY_SZ 256

/***                          D.D. July 28, 2006               (Begin)     ***/
#define WCOL(rr,kk)  rpocol_ptr[rporow_ptrIdx[rr]+kk]
#define MINCOL(rr)   rpocol_ptr[rporow_ptrMin[rr]]
#define MAXCOL(rr)   rpocol_ptr[rporow_ptrIdx[rr]+rporow_ptrNum[rr]-1]
/*******************          D.D. July 28, 2006      (End)  ******************/

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                      STATIC MEMORY FOR THIS OBJECT                        **
**                                                                           **
*******************************************************************************
\*****************************************************************************/

  static int int_road_gravity;
/** D. Donato 08/18/2006  Variable for checking whether to initialize delta.***
    Commented out 8/29/2006
  static GRID_P delta_previous;
*** D. Donato 08/18/2006  Variable for checking whether to initialize delta.**/

/***                          D.D. July 28, 2006               (Begin)      **/
/*** "int *" changed to "short *" 8/10/2006                                 **/
/*** "short *" corrected back to "int *" for rporow_ptrIdx" 8/14/2006       **/
  static int rpoStatus[15];         /* Indicates whether RPO files have been created for this road grid. */
  static int rpoIndex;              /* Points to the RPO file-set in use. */
  static GRID_P rpoList[15];        /* Contains pointers to road grids. */
  static char rpoInitIndic={'n'};   /* Indicates whether the RPO files have been initialized */
  static short *rporow_ptrNum;
  static short *rporow_ptrMin;
  static short *rporow_ptrMax;
  static int   *rporow_ptrIdx;
  static short *rpocol_ptr;
  static int tfoundN, tfoundRow, tfoundCol;
  static int  foundN,  foundRow,  foundCol;
/*******************          D.D. July 28, 2006      (End)  ******************/

/***  D. Donato  Aug. 14, 2006  Moved to module level from spr_phase5       **/
  static int    growth_count;
  static short *growth_row;
  static short *growth_col;
/* The following two lines were replaced by D. Donato on 06/21/2006
****  D. Donato  Aug. 14, 2006                                              **/

/** D. Donato 8/17/2006 Added to deal with cumulative growth                 **/
  static short *zgrwth_row;
  static short *zgrwth_col;
  static int    zgrwth_count;
/** D. Donato 8/17/2006 Added to deal with cumulative growth                 **/

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                        STATIC FUNCTION PROTOTYPES                         **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
static void
    spr_LogSlopeWeights (FILE * fp, int array_size, SWGHT_TYPE * lut);

static void
    spr_phase1n3 (COEFF_TYPE diffusion_coefficient,          /* IN     */
                  COEFF_TYPE breed_coefficient,              /* IN     */
                  GRID_P z,                                  /* IN     */
                  GRID_P delta,                              /* IN/OUT */
                  GRID_P slp,                                /* IN     */
                  GRID_P excld,                              /* IN     */
                  SWGHT_TYPE * swght,                        /* IN     */
                  int *sng,                                  /* IN/OUT */
                  int *sdc);                                 /* IN/OUT */

static void
    spr_phase4 (COEFF_TYPE spread_coefficient,               /* IN     */
                GRID_P z,                                    /* IN     */
                GRID_P excld,                                /* IN     */
                GRID_P delta,                                /* IN/OUT */
                GRID_P slp,                                  /* IN     */
                SWGHT_TYPE * swght,                          /* IN     */
                int *og);                                    /* IN/OUT */


static void
    spr_phase5 (COEFF_TYPE road_gravity,                     /* IN     */
                COEFF_TYPE diffusion_coefficient,            /* IN     */
                COEFF_TYPE breed_coefficient,                /* IN     */
                GRID_P z,                                    /* IN     */
                GRID_P delta,                                /* IN/OUT */
                GRID_P slp,                                  /* IN     */
                GRID_P excld,                                /* IN     */
                GRID_P roads,                                /* IN     */
                SWGHT_TYPE * swght,                          /* IN     */
                int *rt);                                    /* IN/OUT */
  /* D.D. Changed July 24, 2006 - No longer passing workspace
                int *rt,                                     ** IN/OUT **
                GRID_P workspace);                           ** MOD    **
  */
static void
    spr_get_slp_weights (int array_size,                     /* IN     */
                         SWGHT_TYPE * lut);                  /* OUT    */
static BOOLEAN spr_road_search (short i_grwth_center,        /* IN     */
                                short j_grwth_center,        /* IN     */
                                int *i_road,                 /* OUT    */
                                int *j_road,                 /* OUT    */
                                int max_search_index,        /* IN     */
                                GRID_P roads);               /* IN     */
static
  BOOLEAN spr_road_walk (int i_road_start,                   /* IN     */
                         int j_road_start,                   /* IN     */
                         int *i_road_end,                    /* OUT    */
                         int *j_road_end,                    /* OUT    */
                         GRID_P roads,                       /* IN     */
                         double diffusion_coefficient);      /* IN     */
static
  BOOLEAN spr_urbanize_nghbr (int i,                         /* IN     */
                              int j,                         /* IN     */
                              int *i_nghbr,                  /* OUT    */
                              int *j_nghbr,                  /* OUT    */
                              GRID_P z,                      /* IN     */
                              GRID_P delta,                  /* IN     */
                              GRID_P slp,                    /* IN     */
                              GRID_P excld,                  /* IN     */
                              SWGHT_TYPE * swght,            /* IN     */
                              PIXEL pixel_value,             /* IN     */
                              int *stat);                    /* OUT    */
static
  void spr_get_neighbor (int i_in,                           /* IN     */
                         int j_in,                           /* IN     */
                         int *i_out,                         /* OUT    */
                         int *j_out);                        /* OUT    */

static BOOLEAN
    spr_urbanize (int row,                                   /* IN     */
                  int col,                                   /* IN     */
                  GRID_P z,                                  /* IN     */
                  GRID_P delta,                              /* IN     */
                  GRID_P slp,                                /* IN     */
                  GRID_P excld,                              /* IN     */
                  SWGHT_TYPE * swght,                        /* IN     */
                  PIXEL pixel_value,                         /* IN     */
                  int *stat);                                /* OUT    */

static COEFF_TYPE
    spr_GetDiffusionValue (COEFF_TYPE diffusion_coeff);      /* IN    */
static COEFF_TYPE
    spr_GetRoadGravValue (COEFF_TYPE rg_coeff);              /* IN    */

/***                          D.D. July 28, 2006               (Begin)      **/
  static void spr_rpoList_Init(void);
  static void spr_rpoPopulate(int);
  int max( int, int );
/*******************          D.D. July 28, 2006      (End)  ******************/

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                               SCCS ID                                     **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
char spread_c_sccs_id[] = "@(#)spread.c	1.427	12/4/00";

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_phase1n3
** PURPOSE:       perform phase 1 & 3 growth types
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  spr_phase1n3 (COEFF_TYPE diffusion_coefficient,            /* IN     */
                COEFF_TYPE breed_coefficient,                /* IN     */
                GRID_P z,                                    /* IN     */
                GRID_P delta,                                /* IN/OUT */
                GRID_P slp,                                  /* IN     */
                GRID_P excld,                                /* IN     */
                SWGHT_TYPE * swght,                          /* IN     */
                int *sng,                                    /* IN/OUT */
                int *sdc)                                    /* IN/OUT */
{
  char func[] = "spr_phase1n3";
  int i;
  int j;
  int i_out;
  int j_out;
  int k;
  int count;
  int tries;
  int max_tries;
  COEFF_TYPE diffusion_value;
  BOOLEAN urbanized;

  FUNC_INIT;
  assert (MIN_DIFFUSION_VALUE <= diffusion_coefficient);
  assert (diffusion_coefficient <= MAX_DIFFUSION_VALUE);
  assert (MIN_BREED_VALUE <= breed_coefficient);
  assert (breed_coefficient <= MAX_BREED_VALUE);
  assert (z != NULL);
  assert (delta != NULL);
  assert (slp != NULL);
  assert (excld != NULL);
  assert (swght != NULL);
  assert (sng != NULL);
  assert (sdc != NULL);

  diffusion_value = spr_GetDiffusionValue (diffusion_coefficient);
  for (k = 0; k < 1 + (int) diffusion_value; k++)
  {
    i = RANDOM_ROW;
    j = RANDOM_COL;

    if (INTERIOR_PT (i, j))
    {
      if (spr_urbanize (i,                                     /* IN     */
                        j,                                     /* IN     */
                        z,                                     /* IN     */
                        delta,                                 /* IN/OUT */
                        slp,                                   /* IN     */
                        excld,                                 /* IN     */
                        swght,                                 /* IN     */
                        PHASE1G,                               /* IN     */
                        sng))                                  /* IN/OUT */
      {
        if (RANDOM_INT (101) < (int) breed_coefficient)
        {
          count = 0;
          max_tries = 8;
          for (tries = 0; tries < max_tries; tries++)
          {
            urbanized = FALSE;
            urbanized =
              spr_urbanize_nghbr (i,                         /* IN     */
                                  j,                         /* IN     */
                                  &i_out,                    /* OUT    */
                                  &j_out,                    /* OUT    */
                                  z,                         /* IN     */
                                  delta,                     /* IN/OUT */
                                  slp,                       /* IN     */
                                  excld,                     /* IN     */
                                  swght,                     /* IN     */
                                  PHASE3G,                   /* IN     */
                                  sdc);                      /* IN/OUT */
            if (urbanized)
            {
              count++;
              if (count == MIN_NGHBR_TO_SPREAD)
              {
                break;
              }
            }
          }
        }
      }
    }
  }

  FUNC_END;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_phase4
** PURPOSE:       perform phase 4 growth
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  spr_phase4 (COEFF_TYPE spread_coefficient,                 /* IN     */
              GRID_P z,                                      /* IN     */
              GRID_P excld,                                  /* IN     */
              GRID_P delta,                                  /* IN/OUT */
              GRID_P slp,                                    /* IN     */
              SWGHT_TYPE * swght,                            /* IN     */
              int *og)                                       /* IN/OUT */
{
  char func[] = "spr_phase4";
/* D.D. 08/18/2006 */
  int i;
/* D.D. 08/18/2006 */
  int row;
  int col;
  int row_nghbr;
  int col_nghbr;
  int pixel;
  int walkabout_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  int walkabout_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
  int urb_count;
  int nrows;
  int ncols;

  FUNC_INIT;
  assert (z != NULL);
  assert (excld != NULL);
  assert (delta != NULL);
  assert (slp != NULL);
  assert (swght != NULL);
  assert (og != NULL);

  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);

  /*
   *
   * LOOP OVER THE INTERIOR PIXELS LOOKING FOR URBAN FROM WHICH
   * TO PERFORM ORGANIC GROWTH
   *
   */
/** D. Donato 08/18/2006 Use the zgrwth arrays to loop over interior points **
***                      more efficiently.                                  **
  for (row = 1; row < nrows - 1; row++)
  {
    for (col = 1; col < ncols - 1; col++)
    {
***                                                                         */
  for (i=0; i< zgrwth_count; i++)
    {
      row = zgrwth_row[i];
      col = zgrwth_col[i];
      /*
       *
       * D.D. 8/18/2006 -- IS THIS AN INTERIOR PIXEL?
       *
       */
      if (row < 1 || row >= nrows - 1) continue;
      if (col < 1 || col >= ncols - 1) continue;
      /*
       *
       * IS THIS AN URBAN PIXEL AND DO WE PASS THE RANDOM 
       * SPREAD COEFFICIENT TEST
       *
       */
      if ((z[OFFSET (row, col)] > 0) &&
          (RANDOM_INT (101) < spread_coefficient))
      {
        /*
         * EXAMINE THE EIGHT CELL NEIGHBORS
         * SPREAD AT RANDOM IF AT LEAST TWO ARE URBAN
         * PIXEL ITSELF MUST BE URBAN (3)
         *
         */
        urb_count = util_count_neighbors (z, row, col, GT, 0);
        if ((urb_count >= 2) && (urb_count < 8))
        {
          pixel = RANDOM_INT (8);

          row_nghbr = row + walkabout_row[pixel];
          col_nghbr = col + walkabout_col[pixel];

          spr_urbanize (row_nghbr,                           /* IN     */
                        col_nghbr,                           /* IN     */
                        z,                                   /* IN     */
                        delta,                               /* IN/OUT */
                        slp,                                 /* IN     */
                        excld,                               /* IN     */
                        swght,                               /* IN     */
                        PHASE4G,                             /* IN     */
                        og);                                 /* IN/OUT */
        }
      }
    }
/* D.D. 08/18/2006 **
  }
** D.D. 08/18/2006 */
  FUNC_END;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_phase5
** PURPOSE:       perform phase 5 growth
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  spr_phase5 (COEFF_TYPE road_gravity,                       /* IN     */
              COEFF_TYPE diffusion_coefficient,              /* IN     */
              COEFF_TYPE breed_coefficient,                  /* IN     */
              GRID_P z,                                      /* IN     */
              GRID_P delta,                                  /* IN/OUT */
              GRID_P slp,                                    /* IN     */
              GRID_P excld,                                  /* IN     */
              GRID_P roads,                                  /* IN     */
              SWGHT_TYPE * swght,                            /* IN     */
              int *rt)                                       /* IN/OUT */
  /* D.D. Changed July 24, 2006 - No longer passing workspace
              int *rt,                                       ** IN/OUT **
              GRID_P workspace)                              ** MOD    **
  */

{
  char func[] = "spr_phase5";
  int iii;
/** D. Donato Aug. 14, 2006 Following line moved to module level        **
  int growth_count;
***     D. Donato  Aug. 14, 2006                                        */

/* D.D. The following two lines were restored July 24, 2006 */
/* On 8/10/2006 "int *" was changed to "short *".           */
/* On 8/14/2006 these two lines were moved to module level. **
  short *growth_row;
  short *growth_col;
** D. Donato 8/14/2006                                      */

/* The following two lines were replaced by D. Donato on 06/21/2006
  int *growth_row;
  int *growth_col;
*/
/* D.D. The following two lines were rescinded July 24, 2006.
  GRID_P growth_row;
  GRID_P growth_col;
*/
  int max_search_index;
  int growth_index;
  BOOLEAN road_found;
  int i_rd_start;
  int j_rd_start;
  int max_tries;
  BOOLEAN spread;
  BOOLEAN urbanized;
  int i_rd_end;
  int j_rd_end;
  int i_rd_end_nghbr;
  int j_rd_end_nghbr;
  int i_rd_end_nghbr_nghbr;
  int j_rd_end_nghbr_nghbr;
  int tries;
  int nrows;
  int ncols;
  int total_pixels;
float temp4;
int growth_count_fixed;

  FUNC_INIT;
  assert (road_gravity >= 0.0);
  assert (diffusion_coefficient >= 0.0);
  assert (breed_coefficient >= 0.0);
  assert (z != NULL);
  assert (delta != NULL);
  assert (slp != NULL);
  assert (excld != NULL);
  assert (roads != NULL);
  assert (swght != NULL);
  assert (rt != NULL);
/** D.D. Following line disabled.   **
  assert (workspace != NULL);
***                                **/
  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);

  total_pixels = mem_GetTotalPixels ();
  assert (total_pixels > 0);

  /*
   *
   * SET UP WORKSPACE
   *
   */
  
/* The following two lines were replaced by D. Donato 6/21/2006
  growth_row = (int *) workspace;
  growth_col = (int *) workspace + (nrows);
*/
/* D.D. The following two lines were replaced July 24, 2006
  growth_row = (GRID_P) workspace;
  growth_col = (GRID_P) workspace + (nrows);
  On 8/10/2006 "int *" was changed to "short *".
*/
/* the following two lines were moved to spr_spread 8/14/2006  **
  growth_row = (short *)  mem_GetGRCrowptr();
  growth_col = (short *)  mem_GetGRCcolptr();
** D. Donato 8/14/2006                                         */

  /*
   *
   * DETERMINE THE TOTAL GROWTH COUNT AND SAVE THE
   * ROW AND COL LOCATIONS OF THE NEW GROWTH
   *
   */
/** D. Donato Aug. 14, 2006 - Moved to spr_spread          ***
  growth_count = 0;
*** D. Donato Aug. 14, 2006                                **/

#ifdef CRAY_C90
#pragma _CRI ivdep
#endif

/** D. Donato Aug. 14, 2006 - The following loop was made unnecessary **
*** by adding statements to track growth each time a pixel in delta   **
*** is changed. This should save quite a bit of time since for large  **
*** images this loop requires tens of millions of iterations.         **
  for (iii = 0; iii < total_pixels; iii++)
  {
    if (delta[iii] > 0)
    {
      growth_row[growth_count] = iii / ncols;
      growth_col[growth_count] = iii % ncols;
      growth_count++;
    }
  }

*** D. Donato Aug. 14, 2006                                            */

  /*
   *
   * PHASE 5:  ROAD TRIPS
   * IF THERE IS NEW GROWTH, BEGIN PROCESSING ROAD TRIPS
   *
   */

   /*VerD*/
   /*
	fprintf(stdout,"total_pixels = %d\n",total_pixels);
    fprintf(stdout,"growth_count = %d\n",growth_count);
	*/
	/*VerD*/


  if (growth_count > 0)
  {
	  /*VerD*/
	  if (aux_breed_coeff >= 0)
	  {
		  road_growth_breed_coefficient =  aux_breed_coeff;
	  }
	  else
	  {
		  road_growth_breed_coefficient =  - aux_breed_coeff * breed_coefficient;
	  }

	/*
	    for (iii = 0; iii < 1 + (int) (breed_coefficient); iii++)
	*/

growth_count_fixed=growth_count;
	for (iii = 0; iii < 1 + (int) (road_growth_breed_coefficient); iii++)
	/*VerD*/

    {
      /*
       *
       * DETERMINE THE MAX INDEX INTO THE GLB_RD_SEARCH_INDICES ARRAY
       * for road_gravity of 1 we have  8 values
       * for road_gravity of 2 we have 16 values
       * for road_gravity of 3 we have 24 values
       *    and so on....
       *
       * if we need to cover N road_gravity values, then total number of 
       * indexed values would be
       * 8 + 16 + 24 + ... + 8*N = 8*(1+2+3+...+N) = 8*(N(1+N))/2
       *
       */
      int_road_gravity = spr_GetRoadGravValue (road_gravity);
      max_search_index = 4 * (int_road_gravity * (1 + int_road_gravity));
      max_search_index = MAX (max_search_index, nrows);
      max_search_index = MAX (max_search_index, ncols);

      /*
       *
       * RANDOMLY SELECT A GROWTH PIXEL TO START SEARCH
       * FOR ROAD
       *
       */
temp4=RANDOM_FLOAT;
      growth_index = (int) ((double) growth_count_fixed * temp4);

      /*
       *
       * SEARCH FOR ROAD ABOUT THIS GROWTH POINT
       *
       */


      road_found =
        spr_road_search (growth_row[growth_index],
                         growth_col[growth_index],
                         &i_rd_start,
                         &j_rd_start,
                         max_search_index,
                         roads);

      /*
       *
       * IF THERE'S A ROAD FOUND THEN WALK ALONG IT
       *
       */
      if (road_found)
      {
        
		/*VerD*/

		if (aux_diffusion_coeff  >= 0)
		{
			road_growth_diffusion_coefficient =  aux_diffusion_coeff;
		}
		else 
		{
			road_growth_diffusion_coefficient =
							- aux_diffusion_coeff * diffusion_coefficient;
		}

		/*VerD*/

        spread = spr_road_walk (i_rd_start,                  /* IN     */
                                j_rd_start,                  /* IN     */
                                &i_rd_end,                   /* OUT    */
                                &j_rd_end,                   /* OUT    */
                                roads,                       /* IN     */
                                road_growth_diffusion_coefficient);    /* IN     */

        if (spread == TRUE)
        {
          urbanized =
            spr_urbanize_nghbr (i_rd_end,                    /* IN     */
                                j_rd_end,                    /* IN     */
                                &i_rd_end_nghbr,             /* OUT    */
                                &j_rd_end_nghbr,             /* OUT    */
                                z,                           /* IN     */
                                delta,                       /* IN/OUT */
                                slp,                         /* IN     */
                                excld,                       /* IN     */
                                swght,                       /* IN     */
                                PHASE5G,                     /* IN     */
                                rt);                         /* IN/OUT */
          if (urbanized)
          {
            max_tries = 3;
            for (tries = 0; tries < max_tries; tries++)
            {
              urbanized =
                spr_urbanize_nghbr (i_rd_end_nghbr,          /* IN     */
                                    j_rd_end_nghbr,          /* IN     */
                                    &i_rd_end_nghbr_nghbr,   /* OUT    */
                                    &j_rd_end_nghbr_nghbr,   /* OUT    */
                                    z,                       /* IN     */
                                    delta,                   /* IN/OUT */
                                    slp,                     /* IN     */
                                    excld,                   /* IN     */
                                    swght,                   /* IN     */
                                    PHASE5G,                 /* IN     */
                                    rt);                     /* IN/OUT */

            }
          }
        }
      }
    }
  }
  FUNC_END;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_get_slp_weights
** PURPOSE:       calculate the slope weights
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  spr_get_slp_weights (int array_size,                       /* IN     */
                       SWGHT_TYPE * lut)                     /* OUT    */
{
  char func[] = "spr_get_slp_weights";
  float val;
  float exp;
  int i;

  FUNC_INIT;
  assert (lut != NULL);

  exp = coeff_GetCurrentSlopeResist () / (MAX_SLOPE_RESISTANCE_VALUE / 2.0);
  for (i = 0; i < array_size; i++)
  {
    if (i < scen_GetCriticalSlope ())
    {
      val = (scen_GetCriticalSlope () - (SWGHT_TYPE) i) / scen_GetCriticalSlope ();
      lut[i] = 1.0 - pow (val, exp);
    }
    else
    {
      lut[i] = 1.0;
    }
  }
  if (scen_GetLogFlag ())
  {
    if (scen_GetLogSlopeWeightsFlag ())
    {
      scen_Append2Log ();
      spr_LogSlopeWeights (scen_GetLogFP (), array_size, lut);
      scen_CloseLog ();
    }
  }
  FUNC_END;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_LogSlopeWeights
** PURPOSE:       log slope weights to FILE * fp
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  spr_LogSlopeWeights (FILE * fp, int array_size, SWGHT_TYPE * lut)
{
  int i;

  assert (fp != NULL);
  assert (array_size > 0);
  assert (lut != NULL);

  fprintf (fp, "\n%s %5u ***** LOG OF SLOPE WEIGHTS *****\n",
           __FILE__, __LINE__);
  fprintf (fp, "%s %5u CRITICAL_SLOPE= %f\n",
           __FILE__, __LINE__, scen_GetCriticalSlope ());
  fprintf (fp, "%s %5u coeff_GetCurrentSlopeResist= %f\n",
           __FILE__, __LINE__, coeff_GetCurrentSlopeResist ());
  fprintf (fp, "%s %5u MAX_SLOPE_RESISTANCE_VALUE= %f\n",
           __FILE__, __LINE__, MAX_SLOPE_RESISTANCE_VALUE);
  for (i = 0; i < array_size; i++)
  {
    if (i < scen_GetCriticalSlope ())
    {
      fprintf (fp, "%s %5u lut[%3u]= %f\n",
               __FILE__, __LINE__, i, lut[i]);
    }
  }
  fprintf (fp, "All values other values to lut[%u] = 1.000000\n", array_size);

}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_GetDiffusionValue
** PURPOSE:       calculate the diffusion value
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static COEFF_TYPE
  spr_GetDiffusionValue (COEFF_TYPE diffusion_coeff)
{

  COEFF_TYPE diffusion_value;
  double rows_sq;
  double cols_sq;

  rows_sq = igrid_GetNumRows () * igrid_GetNumRows ();
  cols_sq = igrid_GetNumCols () * igrid_GetNumCols ();

    /*VerD*/
  /*
   * diffusion_value's MAXIMUM (IF diffusion_coeff == 100)
   * WILL BE 5% OF THE IMAGE DIAGONAL.
   *  (unless overridden by AUX_DIFFUSION_MULT in scenario file)
   */

	if (aux_diffusion_mult <= 0)
	{
		diffusion_value = ((diffusion_coeff * 0.005) * sqrt (rows_sq + cols_sq));
	}
	else
	{
		diffusion_value = ((diffusion_coeff * aux_diffusion_mult) * sqrt (rows_sq + cols_sq));
	}

    /*VerD*/

  return diffusion_value;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_GetRoadGravValue
** PURPOSE:       calculate the road gravity value
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static COEFF_TYPE
  spr_GetRoadGravValue (COEFF_TYPE rg_coeff)
{

  int rg_value;
  int row;
  int col;

  row = igrid_GetNumRows ();
  col = igrid_GetNumCols ();

  /*
   * rg_value's MAXIMUM (IF rg_coeff == 100)
   * WILL BE 1/16 OF THE IMAGE DIMENSIONS. 
   */

  rg_value = (rg_coeff / MAX_ROAD_VALUE) * ((row + col) / 16.0);

  return rg_value;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_urbanize
** PURPOSE:       try to urbanize a pixel
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static BOOLEAN
  spr_urbanize (int row,                                     /* IN     */
                int col,                                     /* IN     */
                GRID_P z,                                    /* IN     */
                GRID_P delta,                                /* IN/OUT */
                GRID_P slp,                                  /* IN     */
                GRID_P excld,                                /* IN     */
                SWGHT_TYPE * swght,                          /* IN     */
                PIXEL pixel_value,                           /* IN     */
                int *stat)                                   /* IN/OUT */
{
  char func[] = "spr_urbanize";
  BOOLEAN val;
  int nrows;
  int ncols;

  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);

  FUNC_INIT;
  assert ((row >= 0) && (row < nrows));
  assert ((col >= 0) && (col < ncols));
  assert (z != NULL);
  assert (delta != NULL);
  assert (slp != NULL);
  assert (excld != NULL);
  assert (swght != NULL);
  assert (stat != NULL);


  val = FALSE;
  if (z[OFFSET ((row), (col))] == 0)
  {
    if (delta[OFFSET ((row), (col))] == 0)
    {
      if (RANDOM_FLOAT > swght[slp[OFFSET ((row), (col))]])
      {
        if (excld[OFFSET ((row), (col))] < RANDOM_INT (100))
        {
          val = TRUE;
          delta[OFFSET (row, col)] = pixel_value;
            if (pixel_value != 0) /** D. Donato 8/14/2006 - If statement added **/
              {
               growth_row[growth_count] = row;
               growth_col[growth_count] = col;
               growth_count++;
              }

          (*stat)++;
          stats_IncrementUrbanSuccess ();
        }
        else
        {
          stats_IncrementEcludedFailure ();
        }
      }
      else
      {
        stats_IncrementSlopeFailure ();
      }
    }
    else
    {
      stats_IncrementDeltaFailure ();
    }
  }
  else
  {
    stats_IncrementZFailure ();
  }

  FUNC_END;


  return val;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_get_neighbor
** PURPOSE:       find a neighboring pixel
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static
  void
  spr_get_neighbor (int i_in,                                /* IN     */
                    int j_in,                                /* IN     */
                    int *i_out,                              /* OUT    */
                    int *j_out)                              /* OUT    */
{
  char func[] = "spr_get_neighbor";
  int i;
  int j;
  int k;
  int nrows;
  int ncols;

  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();
  assert (nrows > 0);
  assert (ncols > 0);

  FUNC_INIT;
  assert (nrows > i_in);
  assert (ncols > j_in);
  assert (0 <= i_in);
  assert (0 <= j_in);
  assert (i_out != NULL);
  assert (j_out != NULL);

  util_get_next_neighbor (i_in, j_in, i_out, j_out, RANDOM_INT (8));
  for (k = 0; k < 8; k++)
  {
    i = (*i_out);
    j = (*j_out);
    if (IMAGE_PT (i, j))
    {
      break;
    }
    util_get_next_neighbor (i_in, j_in, i_out, j_out, -1);
  }
  FUNC_END;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_urbanize_nghbr
** PURPOSE:       try to urbanize a neighbor
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static
    BOOLEAN
  spr_urbanize_nghbr (int i,                                 /* IN     */
                      int j,                                 /* IN     */
                      int *i_nghbr,                          /* OUT    */
                      int *j_nghbr,                          /* OUT    */
                      GRID_P z,                              /* IN     */
                      GRID_P delta,                          /* IN/OUT */
                      GRID_P slp,                            /* IN     */
                      GRID_P excld,                          /* IN     */
                      SWGHT_TYPE * swght,                    /* IN     */
                      PIXEL pixel_value,                     /* IN     */
                      int *stat)                             /* IN/OUT */
{
  char func[] = "spr_urbanize_nghbr";
  BOOLEAN status = FALSE;

  FUNC_INIT;
  assert (i_nghbr != NULL);
  assert (j_nghbr != NULL);
  assert (z != NULL);
  assert (delta != NULL);
  assert (slp != NULL);
  assert (excld != NULL);
  assert (swght != NULL);
  assert (stat != NULL);

  if (IMAGE_PT (i, j))
  {
    spr_get_neighbor (i,                                     /* IN    */
                      j,                                     /* IN    */
                      i_nghbr,                               /* OUT   */
                      j_nghbr);                              /* OUT   */

    status = spr_urbanize ((*i_nghbr),                       /* IN     */
                           (*j_nghbr),                       /* IN     */
                           z,                                /* IN     */
                           delta,                            /* IN/OUT */
                           slp,                              /* IN     */
                           excld,                            /* IN     */
                           swght,                            /* IN     */
                           pixel_value,                      /* IN     */
                           stat);                            /* IN/OUT */
  }
  FUNC_END;

  return status;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_road_walk
** PURPOSE:       perform road walk
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static
    BOOLEAN
  spr_road_walk (int i_road_start,                           /* IN     */
                 int j_road_start,                           /* IN     */
                 int *i_road_end,                            /* OUT    */
                 int *j_road_end,                            /* OUT    */
                 GRID_P roads,                               /* IN     */
                 double diffusion_coefficient)               /* IN     */
{
  char func[] = "spr_road_walk";
  int i;
  int j;
  int i_nghbr;
  int j_nghbr;
  int k;
  BOOLEAN end_of_road;
  BOOLEAN spread = FALSE;
  int run_value;
  int run = 0;

  FUNC_INIT;
  assert (i_road_end != NULL);
  assert (j_road_end != NULL);
  assert (roads != NULL);

  i = i_road_start;
  j = j_road_start;
  end_of_road = FALSE;
  while (!end_of_road)
  {
    end_of_road = TRUE;
    util_get_next_neighbor (i, j, &i_nghbr, &j_nghbr, RANDOM_INT (8));
    for (k = 0; k < 8; k++)
    {
      if (IMAGE_PT (i_nghbr, j_nghbr))
      {
        if (roads[OFFSET (i_nghbr, j_nghbr)])
        {
          end_of_road = FALSE;
          run++;
          i = i_nghbr;
          j = j_nghbr;
          break;
        }
      }
      util_get_next_neighbor (i, j, &i_nghbr, &j_nghbr, -1);
    }
    run_value = (int) (roads[OFFSET (i, j)] / MAX_ROAD_VALUE *
                       diffusion_coefficient);
    if (run > run_value)
    {
      end_of_road = TRUE;
      spread = TRUE;
      (*i_road_end) = i;
      (*j_road_end) = j;
    }
  }
  FUNC_END;
  return spread;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_road_search
** PURPOSE:       perform road search
** AUTHOR:        David I. Donato
** PROGRAMMER:    David I. Donato, USGS Eastern Geographic Science Center
** CREATION DATE: 08/01/2006
** DESCRIPTION:
**
**
*/
static
    BOOLEAN
  spr_road_search (short i_grwth_center,                     /* IN     */
                   short j_grwth_center,                     /* IN     */
                   int *i_road,                              /* OUT    */
                   int *j_road,                              /* OUT    */
                   int max_search_index,                     /* IN     */
                   GRID_P roads)                             /* IN     */
{
  char func[] = "spr_road_search";
  int i;
  int j;
  int i_offset;
  int j_offset;
  BOOLEAN road_found = FALSE;
  int srch_index;
/***                          D.D. July 28, 2006               (Begin)       **/
/***  New variables for the new road-search algorithm                        **/
  BOOLEAN bn_found;
  int bn, tcount, total, N, k, kmax, kmin, ktest, n, r;
  int crow, ccol, trow, brow, lcol, rcol, nrows, ncols, srow;
  FILE *DD01DBG;
/*******************          D.D. July 28, 2006      (End)  ******************/

  FUNC_INIT;
  assert (i_road != NULL);
  assert (j_road != NULL);
  assert (max_search_index >= 0);

/***                          D.D. July 28, 2006               (Begin)       **/
/***  Set the variables to define the search area.                           **/
  nrows = igrid_GetNumRows();
  ncols = igrid_GetNumCols();
  crow  = i_grwth_center;
  ccol  = j_grwth_center;

  /*** Find the maximal search radius bn. Set N = bn. ***/
  N = MAX (max_search_index, nrows);
  N = MAX (max_search_index, ncols);
  n = ( (int) (sqrt((float)(N/4)) ) - 1);/** Choose an efficient starting   **/
  if (n < 1) {n = 1;}                    /** value for iteration to find bn.**/
  bn_found = FALSE;
  for (bn = n; bn < MAX (ncols, nrows); bn++)
  {
    total = 4 * ((1 + bn) * bn);
    if (total > N)
    {
      bn_found = TRUE;
      break;
    }
  }
  if (!bn_found)
  {
    sprintf (msg_buf, "Unable to find road search band bn in new RS algorithm.");
    LOG_ERROR (msg_buf);
    EXIT (1);
  }
  else
  {N = bn;}

  /** Set trow, brow, lcol, and rcol based on N. **/
  if (crow >= N)               {trow = crow - N;} else {trow = 0;}
  if (crow <= (nrows - 1 - N)) {brow = crow + N;} else {brow = nrows - 1;}
  if (ccol >= N)               {lcol = ccol - N;} else {lcol = 0;}
  if (ccol <= ncols - 1 - N)   {rcol = ccol + N;} else {rcol = ncols - 1;}


/*******************          D.D. July 28, 2006      (End)  ******************/

/***                          D.D. July 28, 2006               (Begin)       **/
/***  Search for road pixels                                                 **/

  foundN = -1;   /** Set the "road-pixel-found" indicator to "none found". **/

  for (srow=0; srow<=N; srow++)
       {
        tfoundN = -1;  /** Pre-set the temporary indicator to "not found". **/

        for (i=1; i<=2; i++ )  /** Process the rows srow above and srow below crow. **/
             {
              if (i == 1) {r = crow - srow;} else {r = crow + srow;}
              /** Process the first row containing the search center only once.     **/
              if (srow == 0 && i == 2) {break;}

              /** Don't waste time processing a row if it is above or below the
                  search area; if it  contains no road pixels; or if there is no
                  overlap between the section with road pixels and the search area. **/

              if (r < trow || r > brow)     {continue;}
              if (rporow_ptrNum[r] == 0)    {continue;}
              if (rporow_ptrMin[r] > rcol)  {continue;}
              if (rporow_ptrMax[r] < lcol)  {continue;}

              /** Perform a binary search to find the road pixel in this row
                  which is closest to the search center.                            **/

              kmin = 0; kmax = rporow_ptrNum[r] - 1;
              tcount = 0;
              while ((kmax - kmin) > 1)
                   {
                    ktest = (kmax + kmin)/2;
                    if (WCOL(r,ktest) <= ccol) {kmin = ktest;}
                    else                       {kmax = ktest;}
                    tcount++;
                   }

              if ( (WCOL(r,kmax) - ccol) < (ccol - WCOL(r,kmin)) )
                   {k = kmax;} else {k = kmin;}
              tfoundN = MAX( srow, abs(WCOL(r,k) - ccol));
              tfoundRow = r; tfoundCol = WCOL(r,k);

              /** Now save the best point found in this row if it is
                  closer to the center than the previous best point. ***/

              if ( (foundN < 0 && tfoundN > 0) || (foundN > 0 && tfoundN > 0 && tfoundN < foundN)) 
                   { foundN = tfoundN; foundRow = tfoundRow; foundCol = tfoundCol; }
             }

        /** If the point found is in the current search band (n = srow) then
        this point is the closest road pixel to the search center
        and there is no need to continue iterating through the remaining
        rows.                                                            ***/

        if (foundN == srow) {break;}
       }


  if ( foundN >= 0 && foundN <= N ) road_found = TRUE; 
  (*i_road) = foundRow;  (*j_road) = foundCol;

  FUNC_END;
  return road_found;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_spread
** PURPOSE:       main spread routine
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  spr_spread (
               float *average_slope,                         /* OUT    */
               int *num_growth_pix,                          /* OUT    */
               int *sng,
               int *sdc,
               int *og,
               int *rt,
               int *pop,
               GRID_P delta,     /* D.D. 8/29/2006  */
               GRID_P z                                      /* IN/OUT */
  ) 
{
  char func[] = "Spread";
/*GRID_P delta;   D.D. 8/29/2006 */
  int i;
  int total_pixels;
  int nrows;
  int ncols;

/* D.D. 8/18/2006 */
  int   row,col, colindex;
  short *ExcPixRow;
  short *ExcPixCol;
  short *UrbPixRow;
  short *UrbPixCol;
/* D.D. 8/18/2006 */

  double road_gravity;
  COEFF_TYPE diffusion_coefficient;
  COEFF_TYPE breed_coefficient;
  COEFF_TYPE spread_coefficient;
  GRID_P excld;
  GRID_P roads;
  GRID_P slp;
/*GRID_P scratch_gif1;  */
  /* D.D. Commented out July 24, 2006 - Now using int growth row and column
          arrays instead of a wgrid.
  GRID_P scratch_gif3;
  */
  SWGHT_TYPE swght[SLOPE_WEIGHT_ARRAY_SZ];

/** The following three lines were moved from spr_phase5                       **/
  growth_row = (short *)  mem_GetGRCrowptr();
  growth_col = (short *)  mem_GetGRCcolptr();
/*******************          D.D. Aug. 14, 2006      (End)  ******************/

  road_gravity = coeff_GetCurrentRoadGravity ();
  diffusion_coefficient = coeff_GetCurrentDiffusion ();
  breed_coefficient = coeff_GetCurrentBreed ();
  spread_coefficient = coeff_GetCurrentSpread ();

/** D.D. 8/29/2006                                                         ***
  scratch_gif1 = mem_GetWGridPtr (__FILE__, func, __LINE__);
*** D.D. 8/29/2006                                                         **/
  /* D.D. Commented out July 24, 2006 - Now using int growth row and column
          arrays instead of a wgrid.
  scratch_gif3 = mem_GetWGridPtr (__FILE__, func, __LINE__);
  */
  excld = igrid_GetExcludedGridPtr (__FILE__, func, __LINE__);
  roads = igrid_GetRoadGridPtrByYear (__FILE__, func,
                                      __LINE__, proc_GetCurrentYear ());
  slp = igrid_GetSlopeGridPtr (__FILE__, func, __LINE__);
  FUNC_INIT;
  assert (road_gravity >= 0.0);
  assert (diffusion_coefficient >= 0.0);
  assert (breed_coefficient >= 0.0);
  assert (spread_coefficient >= 0.0);
  assert (z != NULL);
  assert (excld != NULL);
  assert (roads != NULL);
  assert (slp != NULL);
/** D.D. 8/29/2006                     ***
  assert (scratch_gif1 != NULL);
*** D.D. 8/29/2006                     **/

  /* D.D. Commented out July 24, 2006 - Now using int growth row and column
          arrays instead of a wgrid.
  assert (scratch_gif3 != NULL);
  */
  total_pixels = mem_GetTotalPixels ();
  nrows = igrid_GetNumRows ();
  ncols = igrid_GetNumCols ();

  assert (total_pixels > 0);
  assert (nrows > 0);
  assert (ncols > 0);


  /*
   *
   * SET UP WORKSPACE
   *
   */
/**     D.D. 8/29/2006                                   ***
  delta = scratch_gif1;
***     D.D. 8/29/2006                                   **/

  /*
   *
   * ZERO THE GROWTH ARRAY FOR THIS TIME PERIOD
   *
   */


      for (i=0; i<growth_count; i++)
          {
           delta[OFFSET(growth_row[i], growth_col[i])] = 0;
          }

  growth_count = 0;

  /*
   *
   * GET SLOPE RATES
   *
   */
  spr_get_slp_weights (SLOPE_WEIGHT_ARRAY_SZ,                /* IN     */
                       swght);                               /* OUT    */

/***                          D.D. July 28, 2006               (Begin)       **/
/***  Call the routine to initialize the array of pointers to road grids.    **/

  spr_rpoList_Init();  /* Only runs on the first call */

/*******************          D.D. July 28, 2006      (End)  ******************/

/***                          D.D. July 28, 2006               (Begin)       **/
/*** Identify the correct files for road pixels. Create files if necessary.  **/
  for (rpoIndex=0; rpoIndex<scen_GetRoadDataFileCount(); rpoIndex++)
     {
      if (rpoList[rpoIndex] == roads)
           {
            if (rpoStatus[rpoIndex] == 1) break;
            spr_rpoPopulate(rpoIndex);
            break;
           }
     }
  if (rpoStatus[rpoIndex] != 1) 
/*******************          D.D. July 28, 2006      (End)  ******************/

/***                          D.D. July 28, 2006               (Begin)       **/
/*** Set the road-pixel-only row pointers for the current road grid.         **/
  rporow_ptrNum =  mem_GetRPOrowptrNum(rpoIndex);
  rporow_ptrMin =  mem_GetRPOrowptrMin(rpoIndex);
  rporow_ptrMax =  mem_GetRPOrowptrMax(rpoIndex);
  rporow_ptrIdx =  mem_GetRPOrowptrIdx(rpoIndex);
  rpocol_ptr    =  mem_GetRPOcolptr(rpoIndex);
/*******************          D.D. July 28, 2006      (End)  ******************/

/*******************          D.D. Aug. 14, 2006      (End)  ******************/

/** D. Donato 8/17/2006 Added to deal with cumulative growth                 **/
  zgrwth_row = (short *)  mem_GetGRZrowptr();
  zgrwth_col = (short *)  mem_GetGRZcolptr();
  zgrwth_count =          mem_GetGRZcount();
/** D. Donato 8/17/2006 Added to deal with cumulative growth                 **/

  /*
   *
   * PHASE 1N3 - SPONTANEOUS NEIGHBORHOOD GROWTH AND SPREADING
   *
   */

  timer_Start (SPR_PHASE1N3);
  spr_phase1n3 (diffusion_coefficient,                       /* IN     */
                breed_coefficient,                           /* IN     */
                z,                                           /* IN     */
                delta,                                       /* IN/OUT */
                slp,                                         /* IN     */
                excld,                                       /* IN     */
                swght,                                       /* IN     */
                sng,                                         /* IN/OUT */
                sdc);                                        /* IN/OUT */
  timer_Stop (SPR_PHASE1N3); 

  /*
   *
   * PHASE 4 - ORGANIC GROWTH
   *
   */

  timer_Start (SPR_PHASE4);
  spr_phase4 (spread_coefficient,                            /* IN     */
              z,                                             /* IN     */
              excld,                                         /* IN     */
              delta,                                         /* IN/OUT */
              slp,                                           /* IN     */
              swght,                                         /* IN     */
              og);                                           /* IN/OUT */
  timer_Stop (SPR_PHASE4);

  /*
   *
   * PHASE 5 - ROAD INFLUENCE GROWTH
   *
   */

  timer_Start (SPR_PHASE5);
  spr_phase5 (road_gravity,                                  /* IN     */
              diffusion_coefficient,                         /* IN     */
              breed_coefficient,                             /* IN     */
              z,                                             /* IN     */
              delta,                                         /* IN/OUT */
              slp,                                           /* IN     */
              excld,                                         /* IN     */
              roads,                                         /* IN     */
              swght,                                         /* IN     */
              rt);                                           /* IN/OUT */
   /* D.D. Changed July 24, 2006 - No longer passing 
           scratch_gif3. 
              rt,                                            ** IN/OUT **
              scratch_gif3);                                 ** MOD    **
  */
  timer_Stop (SPR_PHASE5);

/** D.D. 8/18/2006 Use the growth arrays to condition delta more efficiently. ***
  util_condition_gif (total_pixels,                          ** IN     **
                      delta,                                 ** IN     **
                      GT,                                    ** IN     **
                      PHASE5G,                               ** IN     **
                      delta,                                 ** IN/OUT **
                      0);                                    ** IN     */

      for (i=0; i<growth_count; i++)
          {
           if (delta[OFFSET(growth_row[i], growth_col[i])] > PHASE5G)
              {delta[OFFSET(growth_row[i], growth_col[i])] = 0;}
          }

/** D.D. 8/18/2006                                                            **/

/** D.D. 8/18/2006 Use the ExcPix  array to condition delta more efficiently. ***

  util_condition_gif (total_pixels,                          ** IN     **
                      excld,                                 ** IN     **
                      GE,                                    ** IN     **
                      100,                                   ** IN     **
                      delta,                                 ** IN/OUT **
                      0);                                    ** IN     */

  ExcPixRow = igrid_GetExcPixRowPtr();
  ExcPixCol = igrid_GetExcPixColPtr();

  colindex=0;
  for (row=0; row<nrows; row++)
     {
      if (ExcPixRow[row] > 0)
         {
          for (col=0; col<ExcPixRow[row]; col++)
              {
               if (excld[OFFSET(row, ExcPixCol[colindex])] >= 100)
                   delta[OFFSET(row, ExcPixCol[colindex])] = 0;
               colindex++;
              }
         }
     }

/** D.D. 8/18/2006                                                            **/

  /* now place growth array into current array */
  (*num_growth_pix) = 0;
  (*average_slope) = 0.0;

/*** D. Donato Aug. 14, 2006 The following loop has been replaced by a  **
***  more efficient loop. Rather than looping through what is sometimes **
***  tens of millions of 0 pixels (like the original loop), the new     **
***  loop only iterates through actual growth pixels.                   **
  for (i = 0; i < total_pixels; i++)
  {
    if ((z[i] == 0) && (delta[i] > 0))
    {
      ** new growth being placed into array **
      (*average_slope) += (float) slp[i];
      z[i] = delta[i];
      (*num_growth_pix)++;
    }
  }
***  D. Donato Aug. 14, 2006                                            */

/**  D. Donato Aug. 14, 2006             New loop   Corrected 8/17/2006 */
  for (i=0; i<growth_count; i++)
  {
    if (    (z[OFFSET(growth_row[i], growth_col[i])] == 0) &&
        (delta[OFFSET(growth_row[i], growth_col[i])]  > 0)   )
     {
      /* new growth being placed into array */
      (*average_slope) += (float) slp[OFFSET(growth_row[i], growth_col[i])];
          z[OFFSET(growth_row[i], growth_col[i])] =
      delta[OFFSET(growth_row[i], growth_col[i])];
         zgrwth_row[zgrwth_count] = growth_row[i];
         zgrwth_col[zgrwth_count] = growth_col[i];
         zgrwth_count++;
      (*num_growth_pix)++;
     }
  }
/** D. Donato Aug. 14, 2006          End of New Loop                    */

/** D.D. 8/18/2006 Use the zgrwth  arrays to count pixels in z more efficiently. ***
  *pop = util_count_pixels (total_pixels, z, GE, PHASE0G);
***                                                                              **/

  *pop = 0;
  for (i=0; i<zgrwth_count; i++)
  {
   if (z[OFFSET(zgrwth_row[i], zgrwth_col[i])] >= PHASE0G) {*pop++;}
  }

/** D.D. 8/18/2006 Use the zgrwth  arrays to count pixels in z more efficiently. **/

  if (*num_growth_pix == 0)
  {
    *average_slope = 0.0;
  }
  else
  {
    *average_slope /= (float) *num_growth_pix;
  }

/** D. Donato 8/17/2006                                                */
  mem_SetGRZcount(zgrwth_count);
/** D. Donato 8/17/2006                                                */

  roads = igrid_GridRelease (__FILE__, func, __LINE__, roads);
  excld = igrid_GridRelease (__FILE__, func, __LINE__, excld);
  slp = igrid_GridRelease (__FILE__, func, __LINE__, slp);
/** D.D. Following line commented out August 29, 2006                   ***
  scratch_gif1 = mem_GetWGridFree (__FILE__, func, __LINE__, scratch_gif1);
*** D.D.              August 29, 2006                                   **/
/** D.D. Following line commented out July 28, 2006                     ***
  scratch_gif3 = mem_GetWGridFree (__FILE__, func, __LINE__, scratch_gif3);
*** D.D.                  July 28, 2006                                 **/

  FUNC_END;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_rpoList_Init()
** PURPOSE:       initialize the list of RPO pointers to road grids
** AUTHOR:        David I. Donato
** PROGRAMMER:    David I. Donato, USGS Eastern Geographic Science Center
** CREATION DATE: 07/28/2006
** DESCRIPTION:
**
**
*/
static void
   spr_rpoList_Init()
{ 
/*** After execution of this function, rpoList[i] is a pointer to the ith  ***
**** road grid.  The initial status of each rpoList[i] entry is "0" for    ***
**** "not initialized". The road-search data structures still must be      ***
**** built for each road grid. They are built as needed.                   **/

  char func[]="spr_rpoList_Init";
  int i;
  FILE *DD01DBG;

  if (rpoInitIndic == 'Y') {return;} /* Not necessary to reinitialize */

  mem_AllocateRPOcol(); /* Allocate memory for RPO column arrays. */

  for (i=0; i<scen_GetRoadDataFileCount();i++)
       {
        rpoList[i] = igrid_GetRoadGridPtr(__FILE__, func, __LINE__, i);
        rpoStatus[i] = 0;
       }

  rpoInitIndic = 'Y';

  return;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: spr_rpoPopulate()
** PURPOSE:       populate the road-pixel-only (RPO) files for a road grid
** AUTHOR:        David I. Donato
** PROGRAMMER:    David I. Donato, USGS Eastern Geographic Science Center
** CREATION DATE: 07/28/2006
** DESCRIPTION:
**
**
*/
static void
   spr_rpoPopulate(int i)
{  
  char func[] = "spr_rpoPopulate";
  int row, col, rowmax, colmax, numpixels, mincol, maxcol, index;
  FILE *DD01DBG;

  index = 0;
  rowmax = igrid_GetNumRows();
  colmax = igrid_GetNumCols();
  rporow_ptrNum = mem_GetRPOrowptrNum (i);
  rporow_ptrMin = mem_GetRPOrowptrMin (i);
  rporow_ptrMax = mem_GetRPOrowptrMax (i);
  rporow_ptrIdx = mem_GetRPOrowptrIdx (i);
  rpocol_ptr    = mem_GetRPOcolptr(i);

  for (row=0; row < rowmax; row++)
       {
        mincol = colmax; maxcol = 0;
        rporow_ptrNum[row] = 0;
        rporow_ptrIdx[row] = index; 
        for (col=0; col < colmax; col++)
             {
              if (rpoList[i][OFFSET(row,col)] !=0)
                   {
                    rpocol_ptr[index++] = col;
                    if (col < mincol) {mincol = col;}
                    if (col > maxcol) {maxcol = col;}
                    rporow_ptrNum[row]++;
                   }
             }
        rporow_ptrMin[row]=mincol;
        rporow_ptrMax[row]=maxcol;
       }
  rpoStatus[i] = 1;

  return;
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: max
** PURPOSE:       find and return the maximum of two integers
** AUTHOR:        David I. Donato
** PROGRAMMER:    David I. Donato, USGS Eastern Geographic Science Center
** CREATION DATE: 07/28/2006
** DESCRIPTION:
**
**
*/
 int
   max(int i, int j)
{
  if (i > j) {return(i);}
  return(j);
}

