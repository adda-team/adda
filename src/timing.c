/* FILE: timing.c
 * AUTH: Maxim Yurkin
 * DESCR: Basic timing and statistics routines
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <time.h>
#include "vars.h"
#include "comm.h"
#include "const.h"

#define TO_SEC(p) ((p) / (double) CLOCKS_PER_SEC)

/* SEMI-GLOBAL VARIABLES */

/* used in CalculateE.c */
          /* time for Eplane calculation: total,calculation,communication */
clock_t Timing_EFieldPlane,Timing_calc_EField,Timing_comm_EField,
        Timing_IntField,Timing_IntFieldOne, /* for internal fields: total & one calculation */
        Timing_ScatQuan;        /* for integral scattering quantities */
unsigned long TotalEFieldPlane; /* total number of Efield planes calculations */
/* used in calculator.c */
clock_t Timing_Init,     /* for total initialization of the program (before CalculateE) */
        tstart_main;     /* starting time of the program */
unsigned long TotalEval; /* total number of orientation evaluations */
/* used in crosssec.c */
      /* time for all_dir and scat_grid calculations: total,calculation,communication */
clock_t Timing_EField_ad, Timing_calc_EField_ad, Timing_comm_EField_ad,
        Timing_EField_sg, Timing_calc_EField_sg, Timing_comm_EField_sg;
/* used in iterative.c */
clock_t Timing_OneIter,Timing_OneIterCalc, /* for one iteration: total & communication */
        Timing_InitIter;                   /* for initialization of iterative solver */
unsigned long TotalIter;                   /* total number of iterations performed */
/* used in fft.c */
clock_t Timing_FFT_Init,   /* for initialization of FFT routines */
        Timing_Dm_Init;    /* for building Dmatrix */
/* used in make_particle.c */
clock_t Timing_Particle;   /* for particle construction */

/*============================================================*/

void StartTime(void)
  /* start global time */
{
  time(&wt_start);
  last_chp_wt=wt_start;
  tstart_main = clock();
}

/*============================================================*/

void InitTiming(void)
  /* init timing variables and counters */
{
  TotalIter=TotalEval=TotalEFieldPlane=0;
  Timing_EField=Timing_FileIO=Timing_IntField=Timing_ScatQuan=0;
  Timing_Integration=0;
}

/*============================================================*/

void FinalStatistics(void)
  /* print final output and statistics */
{
  time_t wt_end;
  clock_t Timing_TotalTime;

  /* wait for all processes to show correct execution time */
  Synchronize();
  if (ringid==ROOT) {
    /* last time measurements */
    Timing_TotalTime = clock() - tstart_main;
    time(&wt_end);
    /* log statistics */
    fprintf(logfile,"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
                    "                Timing Results             \n"\
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    if (orient_avg) fprintf(logfile,"Total number of single particle evaluations: %d\n",
                            TotalEval);
    fprintf(logfile,
            "Total number of iterations: %d\n"\
            "Total planes of E field calculation (each %d points): %d\n\n"\
            "total time:          %4.4f\n"\
            "Wall time:           %.1f\n"\
            "Initialization time:   %4.4f\n"\
            "  init Dmatrix           %4.4f\n"\
            "  FFT setup:             %4.4f\n"\
            "  make particle:         %4.4f\n"\
            "Internal fields:       %4.4f\n"\
            "  one solution:          %4.4f\n"\
            "    init solver:           %4.4f\n"\
            "    one iteration:         %4.4f\n"\
            "      calculation:           %4.4f\n"\
            "      communication:         %4.4f\n"\
            "E field calculation:   %4.4f\n",
            TotalIter,nTheta,TotalEFieldPlane,
            TO_SEC(Timing_TotalTime),difftime(wt_end,wt_start),
            TO_SEC(Timing_Init),TO_SEC(Timing_Dm_Init),TO_SEC(Timing_FFT_Init),
            TO_SEC(Timing_Particle),
            TO_SEC(Timing_IntField),TO_SEC(Timing_IntFieldOne),TO_SEC(Timing_InitIter),
            TO_SEC(Timing_OneIter),TO_SEC(Timing_OneIterCalc),TO_SEC(Timing_OneIterComm),
            TO_SEC(Timing_EField));
    if (yzplane) fprintf(logfile,
            "  one plane:             %4.4f\n"\
            "    calculation:           %4.4f\n"\
            "    communication:         %4.4f\n",
            TO_SEC(Timing_EFieldPlane),TO_SEC(Timing_calc_EField),TO_SEC(Timing_comm_EField));
    if (all_dir) fprintf(logfile,
            "  one alldir:            %4.4f\n"\
            "    calculation:           %4.4f\n"\
            "    communication:         %4.4f\n",
            TO_SEC(Timing_EField_ad),TO_SEC(Timing_calc_EField_ad),TO_SEC(Timing_comm_EField_ad));
    if (scat_grid) fprintf(logfile,
            "  one scat_grid:            %4.4f\n"\
            "    calculation:           %4.4f\n"\
            "    communication:         %4.4f\n",
            TO_SEC(Timing_EField_sg),TO_SEC(Timing_calc_EField_sg),TO_SEC(Timing_comm_EField_sg));
    fprintf (logfile,
            "Other scat.quantities: %4.4f\n"\
            "file io:               %4.4f\n"\
            "Integration:           %4.4f\n",
            TO_SEC(Timing_ScatQuan),TO_SEC(Timing_FileIO),TO_SEC(Timing_Integration));
    /* close logfile */
    FCloseErr(logfile,F_LOG,ONE_POS);
  }
}

