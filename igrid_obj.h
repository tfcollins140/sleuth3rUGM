/*******************************************************************************
  INCLUDE-FILE NAME:        igrid_obj.h

  CONTAINING SYSTEM:        SLEUTH-3r (based on SLEUTH Model 3.0 Beta)
                            (Slope, Land-cover, Exclusion, Urbanization, 
                            Transportation, and Hillshade)
                            also known as UGM 3.0 Beta (for Urban Growth Model)

  VERSION:                  SLEUTH-3r [Includes Version D features]

  REVISION DATE:            August 31, 2006a
                            [Annotations added August 19, 2009]

*******************************************************************************/

#ifndef IGRID_OBJ_H
#define IGRID_OBJ_H
#include "globals.h"
#include "grid_obj.h"
#include "utilities.h"
#include "ugm_typedefs.h"

typedef struct
{
   char       location[MAX_FILENAME_LEN];
   int        urban_count;
   int        road_count;
   int        landuse_count;
   int        excluded_count;
   int        slope_count;
   int        background_count;

   grid_info   urban[MAX_URBAN_YEARS];
   grid_info   road[MAX_ROAD_YEARS];
   grid_info   landuse[MAX_LANDUSE_YEARS];
   grid_info   excluded;
   grid_info   slope;
   grid_info   background;

}  igrid_info;

void igrid_MemoryLog(FILE* fp);
void igrid_ValidateGrids (FILE* fp);
void igrid_Debug (FILE * fp,char* caller, int location);
void igrid_ReadFiles();
int igrid_GetIGridCount();
igrid_info* igrid_GetStructPtr();
int igrid_GetNumRows();
int igrid_GetNumCols();
int igrid_GetNumTotalPixels();
char* igrid_GetLocation();
int igrid_GetLanduseYear(int i);
int igrid_GetUrbanYear(int i);
int igrid_GetUrbanCount();
void igrid_NormalizeRoads();
void igrid_LogIt(FILE*fp);
void igrid_VerifyInputs(FILE*fp);
int igrid_GetIGridRoadPixelCount(int year);
int igrid_GetIGridExcludedPixelCount();
road_percent_t igrid_GetIGridRoadPercentage(int year);
GRID_P igrid_GridRelease(char* file, char* fun, int line, GRID_P ptr);
GRID_P igrid_GetUrbanGridPtr(char* file, char* fun,int line, int index);
GRID_P igrid_GetUrbanGridPtrByYear(char* file, char* fun,int line, int year);
GRID_P igrid_GetRoadGridPtr(char* file, char* fun,int line, int index);
GRID_P igrid_GetRoadGridPtrByYear(char* file, char* fun,int line, int year);
GRID_P igrid_GetLanduseGridPtr(char* file, char* fun,int line, int index);
GRID_P igrid_GetSlopeGridPtr(char* file, char* fun,int line);
GRID_P igrid_GetExcludedGridPtr(char* file, char* fun,int line);
GRID_P igrid_GetBackgroundGridPtr(char* file, char* fun,int line);
void igrid_Dump(GRID_P ptr,FILE* fp);

BOOLEAN igrid_TestForUrbanYear(int year);
BOOLEAN igrid_TestForRoadYear(int year);
void igrid_init();
int igrid_UrbanYear2Index(int year);

/* D.D. Added July 24, 2006  */
int  igrid_GetIGridRoadPixelCountByIndex (int i);
/* D.D. Added July 24, 2006 */

/* D.D. Added August 16, 2006 */
int  igrid_BuildCompactExcPixFile(void);
int  igrid_BuildCompactUrbPixFile(void);
/* D.D. Added August 16, 2006 */

/* D.D. Added August 18, 2006 */
short* igrid_GetExcPixRowPtr();
short* igrid_GetExcPixColPtr();
short* igrid_GetUrbPixRowPtr();
short* igrid_GetUrbPixColPtr();
GRID_P igrid_GetUrbPixPointer();

/* D.D. Added August 18, 2006 */

#endif

