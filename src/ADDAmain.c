/* File: ADDAmain.c
 * $Date::                            $
 * Descr: main; all the work is done in other modules.
 *
 * Copyright (C) 2006-2013 ADDA contributors
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
// project headers
#include "comm.h"
#include "debug.h"
#include "io.h"
#include "vars.h"
// system headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* When linking both C and Fortran code under Windows, using gcc 4.6.1, EOL style in both stdout and stderr becomes that
 * of Unix. Probably this is due to automatic change of write mode for this streams to binary by gfortran (even when
 * Fortran code has nothing to do with i/o at all). The same may happen to stdin, but this is hard to test.
 *
 * This is not a big problem, but is unpleasant when doing diffs of produced stdout. So until the bug is fixed on
 * compiler side we use two different workarounds for mingw32 and mingw64.
 */
#ifdef __MINGW64_VERSION_STR
	// the simple solution of mingw32 (below) does not work here, so function calls are required
#	include <fcntl.h> // _O_TEXT
#elif defined(__MINGW32_VERSION)
#	include <fcntl.h> // _O_TEXT
	// this definition affects mingw32 program call sequence, which sets the mode of std streams
	unsigned int _CRT_fmode = _O_TEXT;
#endif

// EXTERNAL FUNCTIONS

// calculator.c
void Calculator(void);
// make_particle.c
void InitShape(void);
void MakeParticle(void);
// param.c
void InitVariables(void);
void ParseParameters(int argc,char **argv);
void VariablesInterconnect(void);
void FinalizeSymmetry(void);
void DirectoryLog(int argc,char **argv);
void PrintInfo(void);

//======================================================================================================================

int main(int argc,char **argv)
{
	/* Pointer argv can be declared restrict here and in all calling functions. However, that would be hard to verify,
	 * especially in newly-added functions for parsing command line option. Since the optimization gain is expected to
	 * be minor, if any, we stay conservative on this issue.
	 */
#ifdef __MINGW64_VERSION_STR
	_setmode(_fileno(stdin),_O_TEXT);
	_setmode(_fileno(stdout),_O_TEXT);
	_setmode(_fileno(stderr),_O_TEXT);
#endif
	// Initialize error handling and line wrapping
	logfile=NULL;
	term_width=DEF_TERM_WIDTH;
	// Start global time
	StartTime();
	// Initialize communications
	InitComm(&argc,&argv);
	// Initialize and parse input parameters
	InitVariables();
	ParseParameters(argc,argv);
	D("Reading command line finished");
	VariablesInterconnect(); // also initializes beam
	// Initialize box's; get number of dipoles; set some variables
	InitShape();
	D("Initialization of shape finished");
	FinalizeSymmetry(); // finalize symmetries and check for possible conflicts of symmetries with other options
	// !!! before this line errors should be printed in simple format, after - in advanced one
	// Create directory and start logfile (print command line)
	DirectoryLog(argc,argv);
	// Initialize FFT grid and its subdivision over processors
	ParSetup();
	// MakeParticle; initialize dpl and local_nRows
	MakeParticle();
	D("Make particle finished");
	// Print info to stdout and logfile
	PrintInfo();
	// Main calculation part
	D("Calculator started");
	Calculator();
	D("Calculator finished");
	// Print timing and statistics; close logfile
	FinalStatistics();
	// check error on stdout
	if (ferror(stdout))
		LogWarning(EC_WARN,ALL_POS,"Some errors occurred while writing to stdout during the execution of ADDA");
	// finish execution normally
	Stop(EXIT_SUCCESS);
	// never actually reached; just to make the compiler happy
	return 0;
}
