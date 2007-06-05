/* FILE: ADDAmain.c
 * AUTH: Maxim Yurkin
 * DESCR: Main. All the work moved to other modules.
 *
 *        Previous versions were developed by Alfons Hoekstra.
 *        Sequential version, Michel Grimminck Jan 1995
 *  
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <string.h>
#include "vars.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "io.h"

/* EXTERNAL FUNCTIONS */

/* calculator.c */
void Calculator(void);
/* make_particle.c */
void InitShape(void);
int MakeParticle(void);
/* param.c */
void InitVariables(void);
void ParseParameters(int argc,char **argv);
void VariablesInterconnect(void);
void DirectoryLog(int argc,char **argv);
void PrintInfo(void);
/* timing.c */
void StartTime(void);
void InitTiming(void);
void FinalStatistics(void);

/*============================================================*/

int main(int argc,char **argv)
{
  /* initialize error handling */
  logfile=NULL;
  /* start global time */
  StartTime();
  /* initialize communications */
  InitComm(&argc,&argv);
  /* welcome; initialize and printout version and copyright */
  PRINTZ("'Amsterdam DDA' v." ADDA_VERSION "\n"\
         "Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra\n\n");
  /* initialize and parse input parameters */
  InitVariables();
  ParseParameters(argc,argv);
  D("finished reading command line");
  VariablesInterconnect();
  /* initialize symmetries and box's; get number of dipoles; set some variables */
  InitShape();
  /* Create directory and start logfile (print command line) */
  DirectoryLog(argc,argv);
  /* initialize FFT grid and its subdivision over processors */
  ParSetup();
  /* MakeParticle; initialize dpl and nlocalRows */
  MakeParticle();
  /* print info to stdout and logfile */
  PrintInfo();
  /* initialize times and counters */
  /* Main calculation part */
  D("calculator started");
  Calculator();
  D("calculator finished");
  /* print timing and statistics; close logfile */
  FinalStatistics();
  /* finish execution normally */
  Stop(0);
  /* never actually reached */
  return 0;
}


