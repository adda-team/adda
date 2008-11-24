/* FILE: timing.c
 * AUTH: Maxim Yurkin
 * DESCR: Basic timing and statistics routines
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <time.h>
#include "vars.h"
#include "comm.h"
#include "const.h"
#include "io.h"
#include "timing.h"

#ifdef MPI
#define TO_SEC(p) (p)
#else
#define TO_SEC(p) ((p) / (double) CLOCKS_PER_SEC)
#endif

/* SEMI-GLOBAL VARIABLES */

/* used in CalculateE.c */
TIME_TYPE Timing_EFieldPlane,Timing_comm_EField, /* for Eplane calculation: total and comm */
          Timing_IntField,Timing_IntFieldOne,    /* for internal fields: total & one calculation */
          Timing_ScatQuan;                       /* for integral scattering quantities */
unsigned long TotalEFieldPlane;    /* total number of Efield planes calculations */
/* used in calculator.c */
TIME_TYPE Timing_Init;     /* for total initialization of the program (before CalculateE) */
unsigned long TotalEval;   /* total number of orientation evaluations */
/* used in comm.c */
TIME_TYPE Timing_Dm_Init_comm;  /* communication time for initialization of D-matrix */
/* used in crosssec.c */
      /* total time for all_dir and scat_grid calculations */
TIME_TYPE Timing_EField_ad,Timing_comm_EField_ad,  /* time for all_dir: total & comm */
          Timing_EField_sg,Timing_comm_EField_sg,  /* time for scat_dir: total & comm */
          Timing_ScatQuan_comm;                    /* time for comm of scat.quantities */
/* used in iterative.c */
TIME_TYPE Timing_OneIter,       /* total for one iteration */
          Timing_InitIter,Timing_InitIter_comm;  /* for initialization of iterative solver:
                                                    total & comm */
unsigned long TotalIter;   /* total number of iterations performed */
/* used in fft.c */
TIME_TYPE Timing_FFT_Init,   /* for initialization of FFT routines */
          Timing_Dm_Init;    /* for building Dmatrix */
/* used in make_particle.c */
TIME_TYPE Timing_Particle,                  /* for particle construction */
          Timing_Granul,Timing_Granul_comm; /* for granule generation: total & comm */

/*============================================================*/

void StartTime(void)
  /* start global time */
{
  time(&wt_start);
  last_chp_wt=wt_start;
#ifndef MPI  /* otherwise this initialization is performed immediately after MPI_Init */
  tstart_main = GET_TIME();
#endif
}

/*============================================================*/

void InitTiming(void)
  /* init timing variables and counters */
{
  TotalIter=TotalEval=TotalEFieldPlane=0;
  Timing_EField=Timing_FileIO=Timing_IntField=Timing_ScatQuan=Timing_Integration=0;
  Timing_ScatQuan_comm=Timing_Dm_Init_comm=0;
}

/*============================================================*/

void FinalStatistics(void)
  /* print final output and statistics */
{
  time_t wt_end;
  TIME_TYPE Timing_TotalTime;

  /* wait for all processes to show correct execution time */
  Synchronize();
  if (ringid==ROOT) {
    /* last time measurements */
    Timing_TotalTime = GET_TIME() - tstart_main;
    time(&wt_end);
    /* log statistics */
    fprintf(logfile,
            "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
            "                Timing Results             \n"\
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    if (!prognose) {
      if (orient_avg) fprintf(logfile,
            "Total number of single particle evaluations: %lu\n",TotalEval);
      fprintf(logfile,
            "Total number of iterations: %lu\n"\
            "Total planes of E field calculation (each %d points): %lu\n\n",
            TotalIter,nTheta,TotalEFieldPlane);
    }
    fprintf(logfile,
            "Total wall time:     %.0f\n",difftime(wt_end,wt_start));
    fprintf(logfile,
#ifdef MPI
            "--Everything below is also wall times--\n"\
            "Time since MPI_Init: %.4f\n",
#else
            "--Everything below is processor times--\n"\
            "Total time:          %.4f\n",
#endif
            TO_SEC(Timing_TotalTime));
    fprintf(logfile,
            "  Initialization time: %.4f\n",TO_SEC(Timing_Init));
    if (!prognose) {
      fprintf(logfile,
            "    init Dmatrix         %.4f\n",TO_SEC(Timing_Dm_Init));
#ifdef PARALLEL
      fprintf(logfile,
            "      communication:       %.4f\n",TO_SEC(Timing_Dm_Init_comm));
#endif
      fprintf(logfile,
            "    FFT setup:           %.4f\n",TO_SEC(Timing_FFT_Init));
    }
    fprintf(logfile,
            "    make particle:       %.4f\n",TO_SEC(Timing_Particle));
    if (sh_granul) {
      fprintf(logfile,
            "      granule generator:   %.4f\n",TO_SEC(Timing_Granul));
#ifdef PARALLEL
      fprintf(logfile,
            "        communication:       %.4f\n",TO_SEC(Timing_Granul_comm));
#endif
      }
    if (!prognose) {
      fprintf(logfile,
            "  Internal fields:     %.4f\n"\
            "    one solution:        %.4f\n"\
            "      init solver:         %.4f\n",
            TO_SEC(Timing_IntField),TO_SEC(Timing_IntFieldOne),TO_SEC(Timing_InitIter));
#ifdef PARALLEL
      fprintf(logfile,
            "        communication:       %.4f\n",TO_SEC(Timing_InitIter_comm));
#endif
      fprintf(logfile,
            "      one iteration:       %.4f\n",TO_SEC(Timing_OneIter));
#ifdef PARALLEL
      fprintf(logfile,
            "        communication:       %.4f\n",TO_SEC(Timing_OneIterComm));
#endif
      fprintf(logfile,
            "  Scattered fields:    %.4f\n",TO_SEC(Timing_EField));
      if (yzplane) {
        fprintf(logfile,
            "    one plane:           %.4f\n",TO_SEC(Timing_EFieldPlane));
#ifdef PARALLEL
        fprintf(logfile,
            "      communication:       %.4f\n",TO_SEC(Timing_comm_EField));
#endif
      }
      if (all_dir) {
        fprintf(logfile,
            "    one alldir:          %.4f\n",TO_SEC(Timing_EField_ad));
#ifdef PARALLEL
        fprintf(logfile,
            "      communication:       %.4f\n",TO_SEC(Timing_comm_EField_ad));
#endif
      }
      if (scat_grid) {
        fprintf(logfile,
            "    one scat_grid:          %.4f\n",TO_SEC(Timing_EField_sg));
#ifdef PARALLEL
        fprintf(logfile,
            "      communication:       %.4f\n",TO_SEC(Timing_comm_EField_sg));
#endif
      }
      fprintf (logfile,
            "  Other sc.quantities: %.4f\n",TO_SEC(Timing_ScatQuan));
#ifdef PARALLEL
        fprintf(logfile,
            "    communication:       %.4f\n",TO_SEC(Timing_ScatQuan_comm));
#endif
      fprintf (logfile,
            "  File I/O:            %.4f\n"\
            "Integration:         %.4f\n",
            TO_SEC(Timing_FileIO),TO_SEC(Timing_Integration));
    }
    /* close logfile */
    FCloseErr(logfile,F_LOG,ONE_POS);
  }
}

