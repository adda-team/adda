/* File: io.c
 * $Author$
 * $Date::                            $
 * Descr: i/o routines
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * Copyright (C) 2009 Institute of Chemical Kinetics and Combustion & University of Amsterdam
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
// the following is for MkDirErr
#include "os.h"
#ifdef POSIX
#	include <sys/stat.h>
#	include <sys/types.h>
#endif

#include "io.h"
#include "comm.h"
#include "const.h"
#include "vars.h"
#include "memory.h"

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const char logname[];

// LOCAL VARIABLES

	// error buffer for warning message generated before logfile is opened
static char warn_buf[MAX_MESSAGE2]="";
//============================================================

void WrapLines(char *str)
/* wraps long lines in a string without breaking words; it replaces a number of spaces in string by
 * '\n' characters; line width is determined by variable term_width.
 */
{
	char *left,*right,*mid,*end;
	bool divided;

	end=str+strlen(str);
	left=str;
	while (left<end) {
		// left and right define beginning and end of current working line
		right=strchr(left,'\n');
		if (right==NULL) right=end;
		while ((right-left)>term_width) {
			divided=false;
			mid=left+term_width;
			// search backward for space
			while (mid>=left) {
				if(mid[0]==' ') {
					mid[0]='\n';
					left=mid+1;
					divided=true;
					break;
				}
				mid--;
			}
			// if backward search failed (too long word), search forward for space
			if (!divided) {
				mid=left+term_width+1;
				while (mid<right) {
					if(mid[0]==' ') {
						mid[0]='\n';
						left=mid+1;
						divided=true;
						break;
					}
					mid++;
				}
				// if no spaces are found at all, leave long line and proceed further
				if (!divided) break;
			}
		}
		if (right==end) left=end;
		else left=right+1;
	}
}

//============================================================

char *WrapLinesCopy(const char *str)
/* same as WrapLines, but creates a copy of the string, leaving the original intact; it is designed
 * to be run once during the program (or not many), since this memory is not freed afterwards
 */
{
	char *dup;

	dup=charVector(strlen(str)+1,ONE_POS,"string duplicate");
	strcpy(dup,str);
	WrapLines(dup);
	return dup;
}

//============================================================

void LogError(const int code,const int who,const char *fname,const int lineN,const char *fmt, ... )
/* performs output of error specified by 'code' at 'fname':'lineN'; 'fmt' & arguments ('...')
 * specify error message, 'who' specifies whether 1 (ringid=ROOT) or all processors should produce
 * output. If 'code' is EC_ERROR program aborts after output. We use sprintf a couple of times,
 * because we want each node to generate an atomic message, not a couple of messages after each
 * other, since other nodes may then interfere with our output. INFO is printed to stdout and
 * without showing the position in the source file, ERROR and WARN - to stderr and logfile.
 */
{
	va_list args;
	char line[MAX_MESSAGE2];
	char *pos;

	if (who==ALL || ringid==ROOT) { // controls whether output should be produced
		// first build output string
		va_start(args,fmt);
		if (code==EC_ERROR) strcpy(line,"ERROR: ");
		else if (code==EC_WARN) strcpy(line,"WARNING: ");
		else if (code==EC_INFO) strcpy(line,"INFO: ");
		else sprintf(line,"Error code=%d: ",code);
		pos=line+strlen(line);
#ifdef PARALLEL
		if (code!=EC_INFO) { // for EC_INFO position in source code is not saved
			pos+=sprintf(pos,"(%s:%d) ",fname,lineN);
			// rewrites last 2 chars
			if (who==ALL) {
				pos-=2;
				pos+=sprintf(pos," - ringID=%d) ",ringid);
			}
		}
		else if (who==ALL) pos+=sprintf(pos,"(ringID=%d) ",ringid);
#else
		if (code!=EC_INFO) pos+=sprintf(pos,"(%s:%d) ",fname,lineN);
#endif
		pos+=vsprintf(pos,fmt,args);
		strcpy(pos,"\n");
		va_end(args);
		// print line
		if (code==EC_INFO) {
			// put message to stdout, wrapping lines
			WrapLines(line);
			printf("%s",line);
			fflush(stdout);
		}
		else if (code==EC_ERROR || code==EC_WARN) {
			// first put error message in logfile
			if (logname[0]!=0) { // otherwise can't produce output at all
				if (ringid==ROOT) {
					/* logfile is initialized to NULL in the beginning of the program. Hence if
					 * logfile!=NULL, logfile is necessarily initialized (open or already closed).
					 * logfile==NULL (when logname!=0) means that error is in opening logfile
					 * itself
					 */
					if (logfile!=NULL) {
						if (fprintf(logfile,"%s",line)==EOF) {
							fclose(logfile); // in most cases this is redundant
							// try to reopen logfile and save message
							if ((logfile=fopen(logname,"a"))!=NULL) fprintf(logfile,"%s",line);
						}
						fflush(logfile); // needed for warnings to appear on time
					}
				} // other processors
				else if ((logfile=fopen(logname,"a"))!=NULL) {
					fprintf(logfile,"%s",line);
					fclose(logfile);
				}
			} // save line to buffer to save into logfile afterwards
			else if (code==EC_WARN) strcpy(warn_buf,line);
			// duplicate message to stderr, wrapping lines
			WrapLines(line);
			fprintf(stderr,"%s",line);
			fflush(stderr);
		}
	}
	if (code==EC_ERROR) {
		if (who==ONE && ringid!=ROOT) Synchronize();
		Stop(EXIT_FAILURE);
	}
}

//============================================================

void PrintError(const char *fmt, ... )
/* print anything to stderr (on root processor) and stop; much simpler than LogError (do not print
 * location of the error, and does not duplicate errors to files; assumes that all processors call
 * it.
 */
{
	va_list args;
	char line[MAX_MESSAGE];
	char *pos;

	if (ringid==ROOT) {
		va_start(args,fmt);
		strcpy(line,"ERROR: ");
		pos=line+strlen(line);
		pos+=vsprintf(pos,fmt,args);
		strcpy(pos,"\n");
		va_end(args);
		WrapLines(line);
		fprintf(stderr,"%s",line);
		fflush(stderr);
	}
	// wait for root to generate an error message
	Synchronize();
	Stop(EXIT_FAILURE);
}

//============================================================

void LogPending(void)
/* Logs pending warning messages (currently only one maximum). Should be called when logname is
 * created and logfile is opened on root.
 */
{
	if (warn_buf[0]!=0) {
		if (ringid==ROOT) fprintf(logfile,"%s",warn_buf);
		else if ((logfile=fopen(logname,"a"))!=NULL) {
			fprintf(logfile,"%s",warn_buf);
			fclose(logfile);
		}
		// empty buffer
		warn_buf[0]=0;
	}
}

//============================================================

void PrintBoth(FILE *file,const char *fmt, ... )
/* print anything both to file and to stdout; it is assumed that size of the message is limited to
 * MAX_PARAGRAPH (i.e. no filenames in the message)
 */
{
	va_list args;
	char line[MAX_PARAGRAPH];

	va_start(args,fmt);
	vsprintf(line,fmt,args);
	fprintf(file,"%s",line);
	printf("%s",line);
	va_end(args);
}

//============================================================

FILE *FOpenErr(const char *fname,const char *mode,const int who,const char *err_fname,
	const int lineN)
// open file and check for error
{
	FILE *file;

	if ((file=fopen(fname,mode))==NULL)
		LogError(EC_ERROR,who,err_fname,lineN,"Failed to open file '%s'",fname);
	return file;
}

//============================================================

void FCloseErr(FILE *file,const char *fname,const int who,const char *err_fname,const int lineN)
// close file and check for error
{
	if (fclose(file)) LogError(EC_WARN,who,err_fname,lineN,
		"Errors detected during work with file '%s'",fname);
}

//============================================================

void RemoveErr(const char *fname,const int who,const char *err_fname,const int lineN)
// remove file and check the result
{
	if(remove(fname) && errno!=ENOENT) LogError(EC_WARN,who,err_fname,lineN,
		"Failed to remove temporary file '%s' (%s). Remove it manually, if needed",
		fname,strerror(errno));
}

//============================================================

void MkDirErr(const char *dir,const int who,const char *err_fname,const int lineN)
// make directory and check for error.
{
#ifdef WINDOWS
	if (!CreateDirectory(dir,NULL)) {
		if (GetLastError()==ERROR_ALREADY_EXISTS)
			LogError(EC_WARN,who,err_fname,lineN,"Directory '%s' already exists",dir);
		else LogError(EC_ERROR,who,err_fname,lineN,"Failed to make directory '%s'",dir);
	}
#elif defined(POSIX)
	if (mkdir(dir,0755)==-1) {
		if (errno==EEXIST)
			LogError(EC_WARN,who,err_fname,lineN,"Directory '%s' already exists",dir);
		else LogError(EC_ERROR,who,err_fname,lineN,"Failed to make directory '%s'",dir);
	}
#else
	/* this should work for many cases (unknown OS), but e.g. system calls may fail when used in
	 * combination with MPI
	 */
	char sbuffer[MAX_DIRSYS];

	sprintf(sbuffer,"mkdir \"%s\"",dir);
	if (system(sbuffer)) LogError(EC_WARN,who,err_fname,lineN,"Failed to make directory '%s'",dir);
#endif
}
