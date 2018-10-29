/* File: io.c
 * $Date::                            $
 * Descr: i/o routines
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
#include "const.h" // keep this first
#include "io.h" // corresponding header
// project headers
#include "comm.h"
#include "memory.h"
#include "os.h"
#include "vars.h"
// system headers
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// the following is for MkDirErr
#ifdef POSIX
#	include <sys/stat.h>
#	include <sys/types.h>
#endif

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const char logfname[];

// LOCAL VARIABLES

// error buffer for warning message generated before logfile is opened
static char warn_buf[MAX_MESSAGE2]="";

//======================================================================================================================
/* The following two functions are based on implementations of vasprintf and asprintf from
 * http://stackoverflow.com/a/4900830, which is also similar to the implementation from gnulib. However, explicit error
 * checks has been added inside to avoid need to check the result, and pointer to a new string is returned instead of
 * passing it as an argument.
 */
static char *dyn_vsprintf(const char *format, va_list args)
/* same as vsprintf but allocates storage for the result; simple error handling is used, so it can be called from
 * LogError, etc.
 */
{
	va_list copy;
	va_copy(copy,args);
	char *buffer=NULL; // default value to have deterministic behavior
	int count=vsnprintf(NULL,0,format,args);
	if (count>=0) {
		buffer=(char*)malloc(((size_t)count+1)*sizeof(char));
		if (buffer==NULL) {
			fprintf(stderr,"ERROR: malloc failed in '%s'",__func__);
			Stop(EXIT_FAILURE);
		}
		count=vsnprintf(buffer,(size_t)count+1,format,copy);
	}
	va_end(copy);
	if (count<0) {
		fprintf(stderr,"ERROR: Code %d returned by vsnprintf in '%s'",count,__func__);
		Stop(count);
	}
	return buffer;
}

//======================================================================================================================

char *dyn_sprintf(const char *format, ...)
// same as sprintf, but allocates storage for the result
{
	va_list args;
	va_start(args,format);
	char *res=dyn_vsprintf(format,args);
	va_end(args);
	return res;
}

//======================================================================================================================
// The following two functions are simple modifications of the preceding two, with reallocation
char *rea_vsprintf(char *str,const char *format, va_list args)
/* same as vsprintf but result is added to string str, which is reallocated on the way; simple error handling is used,
 * so it can be called from LogError, etc.
 */{
	va_list copy;
	va_copy(copy,args);
	char *buffer=NULL; // default value to have deterministic behavior
	int count=vsnprintf(NULL,0,format,args);
	if (count>=0) {
		size_t len=strlen(str);
		buffer=(char*)realloc(str,((size_t)count+len+1)*sizeof(char));
		if (buffer==NULL) {
			fprintf(stderr,"ERROR: realloc failed in '%s'",__func__);
			Stop(EXIT_FAILURE);
		}
		count=vsnprintf(buffer+len,(size_t)count+1,format,copy);
	}
	va_end(copy);
	if (count<0) { // simple error handling, so it can be called from LogError, etc.
		fprintf(stderr,"ERROR: Code %d returned by vsnprintf in '%s'",count,__func__);
		Stop(count);
	}
	return buffer;
}

//======================================================================================================================

char *rea_sprintf(char *str,const char *format, ...)
// same as sprintf, but result is added to string str, which is reallocated on the way
{
	va_list args;
	va_start(args,format);
	char *res=rea_vsprintf(str,format,args);
	va_end(args);
	return res;
}

//======================================================================================================================

void WrapLines(char *restrict str)
/* wraps long lines in a string without breaking words; it replaces a number of spaces in string by '\n' characters;
 * line width is determined by variable term_width.
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

//======================================================================================================================

char *WrapLinesCopy(const char * restrict str)
/* same as WrapLines, but creates a copy of the string, leaving the original intact; it is designed to be run once
 * during the program (or not many), since this memory is not freed afterwards
 */
{
	char * restrict dup;

	dup=charVector(strlen(str)+1,ONE_POS,"string duplicate");
	strcpy(dup,str);
	WrapLines(dup);
	return dup;
}

//======================================================================================================================

static void ProcessError(const enum ec code,ERR_LOC_DECL,const char * restrict fmt,va_list args)
/* performs output of error specified by 'code' at 'file':'line'; 'fmt' & arguments ('...') specify error message, 'who'
 * specifies whether 1 (ringid=ADDA_ROOT) or all processors should produce output. We use sprintf-type functions a
 * couple of times, because we want each node to generate an atomic message, not a couple of messages after each other,
 * since other nodes may then interfere with our output. INFO is printed to stdout and without showing the position in
 * the source file, ERROR and WARN - to stderr and logfile.
 */
{
	char msg[MAX_MESSAGE2];
	int shift,tmp;

	if (who==ALL || IFROOT) { // controls whether output should be produced
		// first build output string
		switch (code) {
			case EC_ERROR: strcpy(msg,"ERROR: "); break;
			case EC_WARN: strcpy(msg,"WARNING: "); break;
			case EC_INFO: strcpy(msg,"INFO: "); break;
			case EC_OK: break; // redundant; for a full set of cases
		}
		shift=strlen(msg);
		if (code!=EC_INFO) { // for EC_INFO position in source code is not saved
			SNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE2,"(%s:%d) ",srcfile,srcline);
#ifdef PARALLEL
			if (who==ALL) {
				shift-=2; // rewrites last 2 chars
				SNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE2," - ringID=%d) ",ringid);
			}
		}
		else {
			if (who==ALL) SNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE2,"(ringID=%d) ",ringid);
#endif
		}
		VSNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE2,fmt,args);
		if (shift==MAX_MESSAGE2-1) shift--; // to avoid buffer overflows
		strcpy(msg+shift,"\n");
		// print message
		switch (code) {
			case EC_ERROR:
			case EC_WARN:
				// first put error message in logfile
				if (logfname[0]!=0) { // otherwise can't produce output at all
					if (IFROOT) {
						/* logfile is initialized to NULL in the beginning of the program. Hence if logfile!=NULL,
						 * logfile is necessarily initialized (open or already closed). logfile==NULL (when logname!=0)
						 * means that error is in opening logfile itself
						 */
						if (logfile!=NULL) {
							if (fprintf(logfile,"%s",msg)==EOF) {
								fclose(logfile); // in most cases this is redundant
								// try to reopen logfile and save message
								if ((logfile=fopen(logfname,"a"))!=NULL) fprintf(logfile,"%s",msg);
							}
							fflush(logfile); // needed for warnings to appear on time
						}
					} // other processors
					else if ((logfile=fopen(logfname,"a"))!=NULL) {
						fprintf(logfile,"%s",msg);
						fclose(logfile);
					}
				} // save warning message to buffer to save into logfile afterwards
				else if (code==EC_WARN) strcpy(warn_buf,msg);
				// write (duplicate) message to stderr, wrapping lines
				WrapLines(msg);
				fprintf(stderr,"%s",msg);
				fflush(stderr);
				break;
			case EC_INFO: // put message to stdout, wrapping lines
				WrapLines(msg);
				printf("%s",msg);
				fflush(stdout);
				break;
			case EC_OK: break; // redundant; for a full set of cases
		}
	}
}

//======================================================================================================================

void LogError(ERR_LOC_DECL,const char * restrict fmt, ... )
// A wrapper for ProcessError, which explicitly states fatal error and, hence, exits after the message is printed out.
{
	va_list args;

	va_start(args,fmt);
	ProcessError(EC_ERROR,ERR_LOC_CALL,fmt,args);
	va_end(args);
	/* There seems to be no reliable way to wait for all error messages to appear, since we can't be sure, which
	 * processors call this function. So the following condition is inverse to the one, that determines printing output
	 * by ProcessError function. This way we can be sure that those processors that trigger the final exit, will first
	 * produce some meaningful output.
	 */
	if (!(who==ALL || IFROOT)) Synchronize();
	Stop(EXIT_FAILURE);
}

//======================================================================================================================

void LogWarning(const enum ec code,ERR_LOC_DECL,const char * restrict fmt, ... )
/* Can process any error codes, but always return (do not have noreturn attribute). So it is misleading to call it with
 * EC_ERROR code.
 */
{
	va_list args;

	va_start(args,fmt);
	ProcessError(code,ERR_LOC_CALL,fmt,args);
	va_end(args);
}

//======================================================================================================================


size_t SnprintfErr(ERR_LOC_DECL,char * restrict str,const size_t size,const char * restrict fmt,...)
/* accepts the same arguments as snprintf + location arguments given by macro ERR_LOC_DECL; returns the result of
 * snprintf but first checks that everything went OK, otherwise error is produced. Return type is size_t because it is
 * always non-negative (otherwise error) and this type is more convenient for further use. Using this function adds
 * extra level of safety to building strings, additional to the automatic calculation of buffer size and controlling of
 * all input strings. It is recommended to use it, when potentially long strings (file and directory names) need to be
 * processed.
 */
{
	va_list args;
	int res;
	size_t out;

	va_start(args,fmt);
	res=vsnprintf(str,size,fmt,args);
	if (res<0) LogError(ERR_LOC_CALL,"Error code %d in call to snprintf",res);
	out=(size_t)res;
	if (out>=size) LogError(ERR_LOC_CALL,"Buffer overflow in call to snprintf");
	va_end(args);
	return out;
}

//======================================================================================================================

size_t SnprintfShiftErr(ERR_LOC_DECL,const size_t shift,char * restrict str,const size_t size,
	const char * restrict fmt,...)
/* Same as SnprintfErr, but accepts one extra argument - 'shift'. Writing of new data to 'str' starts after shift
 * characters (i.e. shift=0 is equivalent to SnprintfErr). 'size' specifies the memory available to the whole string
 * 'str', not to the remaining part. Returns the new (total) shift after write.
 */
{
	va_list args;
	int res;
	size_t out;

	va_start(args,fmt);
	res=vsnprintf(str+shift,size-shift,fmt,args);
	if (res<0) LogError(ERR_LOC_CALL,"Error code %d in call to snprintf",res);
	out=(size_t)res+shift;
	if (out>=size) LogError(ERR_LOC_CALL,"Buffer overflow in call to snprintf");
	va_end(args);
	return out;
}

//======================================================================================================================

size_t VsnprintfErr(ERR_LOC_DECL,char * restrict str,const size_t size,const char * restrict fmt,
	va_list args)
/* accepts the same arguments as vsnprintf + location arguments given by macro ERR_LOC_DECL; returns the result of
 * snprintf but first checks that everything went OK, otherwise error is produced. Return type is size_t because it is
 * always non-negative (otherwise error) and this type is more convenient for further use.
  */
{
	int res;
	size_t out;

	res=vsnprintf(str,size,fmt,args);
	if (res<0) LogError(ERR_LOC_CALL,"Error code %d in call to vsnprintf",res);
	out=(size_t)res;
	if (out>=size) LogError(ERR_LOC_CALL,"Buffer overflow in call to vsnprintf");
	return out;
}


//======================================================================================================================

void PrintError(const char * restrict fmt, ... )
/* print anything to stderr (on root processor) and stop; much simpler than LogError (do not print location of the
 * error, and does not duplicate errors to files; assumes that all processors call it.
 */
{
	va_list args;
	char msg[MAX_MESSAGE]="ERROR: ";
	int shift,tmp;

	if (IFROOT) {
		va_start(args,fmt);
		shift=strlen(msg);
		VSNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE,fmt,args);
		va_end(args);
		WrapLines(msg);
		fprintf(stderr,"%s\n",msg);
		fflush(stderr);
	}
	// wait for root to generate an error message
	Synchronize();
	Stop(EXIT_FAILURE);
}

//======================================================================================================================

void LogPending(void)
/* Logs pending warning messages (currently only one maximum). Should be called when logname is created and logfile is
 * opened on root.
 */
{
	if (warn_buf[0]!=0) {
		if (IFROOT) fprintf(logfile,"%s",warn_buf);
		else if ((logfile=fopen(logfname,"a"))!=NULL) {
			fprintf(logfile,"%s",warn_buf);
			fclose(logfile);
		}
		// empty buffer
		warn_buf[0]=0;
	}
}

//======================================================================================================================

void PrintBoth(FILE * restrict file,const char * restrict fmt, ... )
/* print anything both to file and to stdout; it is assumed that size of the message is limited to MAX_PARAGRAPH (i.e.
 * no filenames in the message). Makes little sense to call it by all processors.
 *
 * Not thread-safe! Should not be called in parallel from multiple threads (e.g. OpenMP)
 */
{
	va_list args;
	static char msg[MAX_PARAGRAPH]; // not to allocate at every call

	va_start(args,fmt);
	/* ALL_POS will not contain a lot of useful information here (since it is not related to the origin of function
	 * call). However, an error here is extremely improbable.
	 */
	VsnprintfErr(ALL_POS,msg,MAX_PARAGRAPH,fmt,args);
	fprintf(file,"%s",msg);
	printf("%s",msg);
	va_end(args);
}

//======================================================================================================================

FILE *FOpenErr(const char * restrict fname,const char * restrict mode,ERR_LOC_DECL)
// open file and check for error
{
	FILE * restrict file;

	if ((file=fopen(fname,mode))==NULL) LogError(ERR_LOC_CALL,"Failed to open file '%s'",fname);
	return file;
}

//======================================================================================================================

void FCloseErr(FILE * restrict file,const char * restrict fname,ERR_LOC_DECL)
// close file and check for error
{
	if (fclose(file))
		LogWarning(EC_WARN,ERR_LOC_CALL,"Errors detected during work with file '%s'",fname);
}

//======================================================================================================================

void RemoveErr(const char * restrict fname,ERR_LOC_DECL)
// remove file and check the result
{
	if(remove(fname) && errno!=ENOENT) LogWarning(EC_WARN,ERR_LOC_CALL,
		"Failed to remove temporary file '%s' (%s). Remove it manually, if needed",fname,strerror(errno));
}

//======================================================================================================================

void MkDirErr(const char * restrict dir,ERR_LOC_DECL)
// make directory and check for error.
{
#ifdef WINDOWS
	if (!CreateDirectory(dir,NULL)) {
		if (GetLastError()==ERROR_ALREADY_EXISTS) LogWarning(EC_WARN,ERR_LOC_CALL,"Directory '%s' already exists",dir);
		else LogError(ERR_LOC_CALL,"Failed to make directory '%s'",dir);
	}
#elif defined(POSIX)
	if (mkdir(dir,0755)==-1) {
		if (errno==EEXIST) LogWarning(EC_WARN,ERR_LOC_CALL,"Directory '%s' already exists",dir);
		else LogError(ERR_LOC_CALL,"Failed to make directory '%s'",dir);
	}
#else
	/* this should work for many cases (unknown OS), but e.g. system calls may fail when used in
	 * combination with MPI
	 */
	char sbuffer[MAX_DIRSYS];

	sprintf(sbuffer,"mkdir \"%s\"",dir);
	if (system(sbuffer)) LogWarning(EC_WARN,ERR_LOC_CALL,"Failed to make directory '%s'",dir);
#endif
}

//======================================================================================================================

static inline void SkipFullLine(FILE * restrict file,char * restrict buf,const int buf_size)
// skips full line in the file, starting from current position; uses buffer 'buf' with size 'buf_size'
{
	do fgets(buf,buf_size,file); while (strchr(buf,'\n')==NULL && !feof(file));
}

//======================================================================================================================

char *FGetsError(FILE * restrict file,const char * restrict fname,size_t *line,char * restrict buf,const int buf_size,
	ERR_LOC_DECL)
/* calls fgets, checks for errors and increments line number; s_fname and s_line are source fname and line number to be
 * shown in error message; uses buffer 'buf' with size 'buf_size'.
 */
{
	char *res;

	res=fgets(buf,buf_size,file);
	if (res!=NULL) {
		(*line)++;
		if (strchr(buf,'\n')==NULL && !feof(file)) LogError(ERR_LOC_CALL,
			"Buffer overflow while scanning lines in file '%s' (size of line %zu > %d)",fname,*line,BUF_LINE-1);
	}
	return res;
}

//======================================================================================================================

size_t SkipNLines(FILE * restrict file,const size_t n)
// skips n lines from the file starting from current position in a file; returns n
{
	char buf[BUF_LINE];
	size_t i;

	for (i=0;i<n;i++) SkipFullLine(file,buf,BUF_LINE);
	return n;
}

//======================================================================================================================

size_t SkipComments(FILE * restrict file)
// skips comments (#...), all lines, starting from current position in a file. returns number of lines skipped
{
	int ch;
	size_t lines=0;
	char buf[BUF_LINE];

	while ((ch=fgetc(file))=='#') {
		SkipFullLine(file,buf,BUF_LINE);
		lines++;
	}
	if (ch!=EOF) ungetc(ch,file);

	return lines;
}
