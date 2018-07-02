/*******************************************************************************
  INCLUDE-FILE NAME:        spread.h

  CONTAINING SYSTEM:        SLEUTH-3r (based on SLEUTH Model 3.0 Beta)
                            (Slope, Land-cover, Exclusion, Urbanization, 
                            Transportation, and Hillshade)
                            also known as UGM 3.0 Beta (for Urban Growth Model)

  VERSION:                  SLEUTH-3r [Includes Version D features]

  REVISION DATE:            August 31, 2006a
                            [Annotations added August 19, 2009]

*******************************************************************************/

#ifndef SPREAD_H
#define SPREAD_H

#ifdef SPREAD_MODULE
  /* stuff visable only to the growth module */
char spread_h_sccs_id[] = "@(#)spread.h	1.243	12/4/00";

  
#endif


/* #defines visable to any module including this header file*/

/*
 *
 * FUNCTION PROTOTYPES
 *
 */

void
  spr_spread  (
              float *average_slope,                          /* OUT    */
              int *num_growth_pix,                           /* OUT    */
              int* sng,
              int* sdc,
              int* og,
              int* rt,
              int* pop,
              GRID_P delta,   /** D.D. 8/29/2006  **/
              GRID_P z                                     /* IN/OUT */
              );                       /* MOD    */

#endif
