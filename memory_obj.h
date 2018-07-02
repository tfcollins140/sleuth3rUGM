/*******************************************************************************
  INCLUDE-FILE NAME:        memory_obj.h

  CONTAINING SYSTEM:        SLEUTH-3r (based on SLEUTH Model 3.0 Beta)
                            (Slope, Land-cover, Exclusion, Urbanization, 
                            Transportation, and Hillshade)
                            also known as UGM 3.0 Beta (for Urban Growth Model)

  VERSION:                  SLEUTH-3r [Includes Version D features]

  REVISION DATE:            August 31, 2006a
                            [Annotations added August 19, 2009]

*******************************************************************************/

#include "globals.h"
#include "ugm_typedefs.h"

void mem_Init();
void mem_MemoryLog(FILE* fp);

void mem_LogPartition(FILE* fp);
int mem_GetPackedBytesPerGrid();
GRID_P mem_GetIGridPtr( char* owner );

GRID_P mem_GetPGridPtr( char* owner );

GRID_P mem_GetWGridPtr( char* module, char* who, int line );

GRID_P mem_GetWGridFree( char* module, char* who, int line,GRID_P ptr );

int mem_GetTotalPixels();

void mem_CheckMemory(FILE* fp,char* module, char* function, int line);

void mem_ReinvalidateMemory();

int memGetBytesPerGridRound();

void mem_LogMinFreeWGrids(FILE* fp);
FILE* mem_GetLogFP();
void mem_CloseLog();

/* D.D. Added for growth Row and Column (GRC) arrays - July 28, 2006 */
/* D.D. Return type changed to short August 10, 2006                 */
short *mem_GetGRCrowptr ();
short *mem_GetGRCcolptr ();

/* D.D. Added for cumulative growth array - 8/17/2006                */
short *mem_GetGRZrowptr ();
short *mem_GetGRZcolptr ();
int    mem_GetGRZcount ();
void   mem_SetGRZcount (int);
GRID_P mem_GetGRZpointer();
void   mem_SetGRZpointer(GRID_P);
/* D.D. Added for cumulative growth array - 8/17/2006                */

short *mem_GetRPOrowptrNum (int i);
short *mem_GetRPOrowptrMin (int i);
short *mem_GetRPOrowptrMax (int i);
int   *mem_GetRPOrowptrIdx (int i);
short *mem_GetRPOcolptr (int i);
/**  D.D.  July 28, 2006                                   *******************/

/**  D.D.  Added for road-pixel-only (RPO) column arrays - Aug. 1, 2006    ***/
void mem_AllocateRPOcol();
/**  D.D.  Aug. 1, 2006                                                    ***/

