/* File: timing.c
 * $Date::                            $
 * Descr: basic timing and statistics routines
 *
 * Copyright (C) 2006,2008-2014 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include "const.h" // keep this first
#include "timing.h" // corresponding header
// project headers
#include "comm.h"
#include "io.h"
#include "vars.h"
// system headers
#include <math.h>
#include <time.h>
#include <stdio.h>

#ifdef ADDA_MPI
#	define TO_SEC(p) (p)
#else
#	define TO_SEC(p) ((p) / (double) CLOCKS_PER_SEC)
#endif

// SEMI-GLOBAL VARIABLES

// used in CalculateE.c
TIME_TYPE Timing_EPlane,Timing_EPlaneComm,    // for Eplane calculation: total and comm
          Timing_IntField,Timing_IntFieldOne, // for internal fields: total & one calculation
          Timing_IncBeam,                     // for generation (or reading) and saving (if needed) of incident beam
          Timing_ScatQuan;                    // for integral scattering quantities
size_t TotalEFieldPlane; // total number of planes for scattered field calculations
// used in calculator.c
TIME_TYPE Timing_Init, // for total initialization of the program (before CalculateE)
          Timing_Init_Int; // for initialization of interaction routines (including computing tables)
size_t TotalEval;      // total number of orientation evaluations
#ifdef OPENCL
TIME_TYPE Timing_OCL_Init; // for initialization of OpenCL (including building program)
#endif
// used in comm.c
TIME_TYPE Timing_InitDmComm; // communication time for initialization of D-matrix
// used in crosssec.c
	// total time for all_dir and scat_grid calculations
TIME_TYPE Timing_EFieldAD,Timing_EFieldADComm,  // time for all_dir: total & comm
          Timing_EFieldSG,Timing_EFieldSGComm,  // time for scat_dir: total & comm
          Timing_ScatQuanComm;                  // time for comm of scat.quantities
// used in fft.c
TIME_TYPE Timing_FFT_Init, // for initialization of FFT routines
          Timing_Dm_Init;  // for building Dmatrix
// used in iterative.c
time_t last_chp_wt; // wall time of the last checkpoint (1s precision is sufficient)
TIME_TYPE Timing_OneIter,Timing_OneIterComm,       // for one iteration: total & comm
          Timing_InitIter,Timing_InitIterComm,     // for initialization of iterations: total & comm
          Timing_IntFieldOneComm,                  // comm for one calculation of the internal fields
          Timing_MVP,Timing_MVPComm,               // total & comm time for MatVec during one run of iterative solver
          Timing_OneIterMVP,Timing_OneIterMVPComm; // total & comm time for MatVec during one iteration
size_t TotalIter;                               // total number of iterations performed
// used in make_particle.c
TIME_TYPE Timing_Particle,                 // for particle construction
          Timing_Granul,Timing_GranulComm; // for granule generation: total & comm
// used in matvec.c
size_t TotalMatVec; // total number of matrix-vector products

// LOCAL VARIABLES
SYSTEM_TIME wt_start; // starting wall time

#define FFORMT "%.4f" // format for timing results

//======================================================================================================================

double DiffSystemTime(const SYSTEM_TIME * restrict t1,const SYSTEM_TIME * restrict t2)
/* compute difference (in seconds) between two system times; not very fast (in contrast to functions in prec_time.c/h)
 * !!! order of arguments is inverse to that in standard difftime (for historical reasons)
 */
{
#ifdef WINDOWS
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	return (double)(t2->QuadPart - t1->QuadPart)/(double)(freq.QuadPart);
#elif defined(POSIX)
	return (double)(t2->tv_sec - t1->tv_sec) + MICRO*(double)(t2->tv_usec - t1->tv_usec);
#else // fallback for 1s-precision timer
	return difftime(*t2,*t1);
#endif
}

//======================================================================================================================

void StartTime(void)
// start global time
{
	GET_SYSTEM_TIME(&wt_start);
	time(&last_chp_wt);
#ifndef ADDA_MPI  // otherwise this initialization is performed immediately after MPI_Init
	tstart_main = GET_TIME();
#endif
}

//======================================================================================================================

void InitTiming(void)
// init timing variables and counters
{
	TotalIter=TotalMatVec=TotalEval=TotalEFieldPlane=0;
	Timing_EField=Timing_FileIO=Timing_IntField=Timing_ScatQuan=Timing_Integration=0;
	Timing_ScatQuanComm=Timing_InitDmComm=0;
#ifdef SPARSE
	Timing_Dm_Init=Timing_Granul=Timing_FFT_Init=Timing_GranulComm=0;
#endif	
}

//======================================================================================================================

void FinalStatistics(void)
// print final output and statistics
{
	SYSTEM_TIME wt_end;
	double totTime;
	TIME_TYPE Timing_TotalTime;

	// wait for all processes to show correct execution time
	Synchronize();
	if (IFROOT) {
		// last time measurements
		Timing_TotalTime = GET_TIME() - tstart_main;
		GET_SYSTEM_TIME(&wt_end);
		// log statistics
		fprintf(logfile,
			"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
			"                Timing Results             \n"
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		if (!prognosis) {
			if (orient_avg) fprintf(logfile,
				"Total number of single particle evaluations: %zu\n",TotalEval);
			fprintf(logfile,
				"Total number of iterations: %zu\n"
				"Total number of matrix-vector products: %zu\n"
				"Total planes of E field calculation (each %d points): %zu\n\n",
				TotalIter,TotalMatVec,nTheta,TotalEFieldPlane);
		}
		fprintf(logfile,
			"Total wall time:     "FFORMT"\n",totTime=DiffSystemTime(&wt_start,&wt_end));
#ifdef ADDA_MPI
		fprintf(logfile,
			"--Everything below is also wall times--\n"
			"Time since MPI_Init: "FFORMT"\n",TO_SEC(Timing_TotalTime));
#else // standard clock
		fprintf(logfile,
			"--Everything below is processor times--\n");
		/* Here we test for possible overflow of clock. If clock_t is only 4 bytes (e.g. long) and CLOCLS_PER_SEC = 10^6
		 * (not 1000), the overflow is expected whenever total time is larger than 72 min.
		 */
		if (CLOCKS_PER_SEC*totTime > pow(256,sizeof(clock_t))-1) fprintf(logfile,
			"--(some values are affected by timer overflow)--\n");
		fprintf(logfile,
			"Total time:          "FFORMT"\n",TO_SEC(Timing_TotalTime));
#endif
		fprintf(logfile,
			"  Initialization time: "FFORMT"\n",TO_SEC(Timing_Init));
		if (!prognosis) {
#ifdef OPENCL
			fprintf(logfile,
				"    init OpenCL          "FFORMT"\n",TO_SEC(Timing_OCL_Init));

#endif
			fprintf(logfile,
				"    init interaction     "FFORMT"\n",TO_SEC(Timing_Init_Int));
#ifndef SPARSE
			fprintf(logfile,
				"    init Dmatrix         "FFORMT"\n",TO_SEC(Timing_Dm_Init));
#	ifdef PARALLEL
			fprintf(logfile,
				"      communication:       "FFORMT"\n",TO_SEC(Timing_InitDmComm));
#	endif
			fprintf(logfile,
				"    FFT setup:           "FFORMT"\n",TO_SEC(Timing_FFT_Init));
#endif // !SPARSE
		}
		fprintf(logfile,
			"    make particle:       "FFORMT"\n",TO_SEC(Timing_Particle));
		if (sh_granul) {
			fprintf(logfile,
				"      granule generator:   "FFORMT"\n",TO_SEC(Timing_Granul));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_GranulComm));
#endif
		}
		if (!prognosis) {
			fprintf(logfile,
				"  Internal fields:     "FFORMT"\n"
				"    one solution:        "FFORMT"\n",
				TO_SEC(Timing_IntField),TO_SEC(Timing_IntFieldOne));
#ifdef PARALLEL
			fprintf(logfile,
				"      communication:       "FFORMT"\n",TO_SEC(Timing_IntFieldOneComm));
#endif
			fprintf(logfile,
				"      matvec products:     "FFORMT"\n",TO_SEC(Timing_MVP));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_MVPComm));
#endif
			fprintf(logfile,
				"      incident beam:       "FFORMT"\n",TO_SEC(Timing_IncBeam));
			fprintf(logfile,
				"      init solver:         "FFORMT"\n",TO_SEC(Timing_InitIter));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_InitIterComm));
#endif
			fprintf(logfile,
				"      one iteration:       "FFORMT"\n",TO_SEC(Timing_OneIter));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_OneIterComm));
#endif
			fprintf(logfile,
				"        matvec products:     "FFORMT"\n",TO_SEC(Timing_OneIterMVP));
#ifdef PARALLEL
			fprintf(logfile,
				"          communication:       "FFORMT"\n",TO_SEC(Timing_OneIterMVPComm));
#endif
			fprintf(logfile,
				"  Scattered fields:    "FFORMT"\n",TO_SEC(Timing_EField));
			if (yzplane || scat_plane) {
				fprintf(logfile,
					"    one plane:           "FFORMT"\n",TO_SEC(Timing_EPlane));
#ifdef PARALLEL
				fprintf(logfile,
					"      communication:       "FFORMT"\n",TO_SEC(Timing_EPlaneComm));
#endif
			}
			if (all_dir) {
				fprintf(logfile,
					"    one alldir:          "FFORMT"\n",TO_SEC(Timing_EFieldAD));
#ifdef PARALLEL
				fprintf(logfile,
					"      communication:       "FFORMT"\n",TO_SEC(Timing_EFieldADComm));
#endif
			}
			if (scat_grid) {
				fprintf(logfile,
					"    one scat_grid:          "FFORMT"\n",TO_SEC(Timing_EFieldSG));
#ifdef PARALLEL
				fprintf(logfile,
					"      communication:       "FFORMT"\n",TO_SEC(Timing_EFieldSGComm));
#endif
			}
			fprintf (logfile,
				"  Other sc.quantities: "FFORMT"\n",TO_SEC(Timing_ScatQuan));
#ifdef PARALLEL
			fprintf(logfile,
				"    communication:       "FFORMT"\n",TO_SEC(Timing_ScatQuanComm));
#endif
		}
		fprintf (logfile,
				"File I/O:            "FFORMT"\n",TO_SEC(Timing_FileIO));
		if (!prognosis) fprintf (logfile,
				"Integration:         "FFORMT"\n",TO_SEC(Timing_Integration));
		// close logfile
		FCloseErr(logfile,F_LOG,ONE_POS);
	}
}
