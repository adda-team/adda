/* File: ADDAmain.c
 * $Author$
 * $Date::                            $
 * Descr: main; all the work is done in other modules.
 *
 *        Previous versions were developed by Alfons Hoekstra.
 *        Sequential version, Michel Grimminck January 1995
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <string.h>
#include "vars.h"
#include "comm.h"
#include "debug.h"
#include "io.h"

// EXTERNAL FUNCTIONS

// calculator.c
void Calculator(void);
// make_particle.c
void InitShape(void);
int MakeParticle(void);
// param.c
void InitVariables(void);
void ParseParameters(int argc,char **argv);
void VariablesInterconnect(void);
void DirectoryLog(int argc,char **argv);
void PrintInfo(void);

//============================================================

int main(int argc,char **argv)
{
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
	// Initialize symmetries and box's; get number of dipoles; set some variables
	InitShape();
	// !!! before this line errors should be printed in simple format, after - in advanced one
	// Create directory and start logfile (print command line)
	DirectoryLog(argc,argv);
	// Initialize FFT grid and its subdivision over processors
	ParSetup();
	// MakeParticle; initialize dpl and nlocalRows
	MakeParticle();
	// Print info to stdout and logfile
	PrintInfo();
	// Main calculation part
	D("Calculator started");
	Calculator();
	D("Calculator finished");
	// Print timing and statistics; close logfile
	FinalStatistics();
	// check error on stdout
	if (ferror(stdout)) LogError(EC_WARN,ALL_POS,
		"Some errors occurred while writing to stdout during the execution of ADDA");
	// finish execution normally
	Stop(0);
	// never actually reached; just to make the compiler happy
	return 0;
}
