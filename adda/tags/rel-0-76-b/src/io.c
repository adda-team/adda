/* FILE: io.c
 * AUTH: Maxim Yurkin
 * DESCR: io routines
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
/* the following is for MkDirErr */
#include "os.h"
#ifdef POSIX
# include <sys/stat.h>
# include <sys/types.h>
#endif

#include "io.h"
#include "comm.h"
#include "const.h"
#include "vars.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in param.c */
extern char logname[];

/* LOCAL VARIABLES */

static char warn_buf[MAX_MESSAGE]=""; /* error buffer for warning message generated before
                                         logfile is opened */

/*============================================================*/

void LogError(const int code,const int who,const char *fname,
              const int lineN,const char *fmt, ... )
   /* performs output of error specified by code at fname:lineN
    * fmt + args (...) specify error message
    * who specifies whether 1 (ringid=ROOT) or all processors should produce output
    * if code is EC_ERROR program aborts after output.
    * We use sprintf a couple of times, because we want each node to
    * generate an atomic message, not a couple of messages after
    * each other, since other nodes may then interfere with our output
    * INFO is printed to stdout
    * ERROR and WARN - to stdout and logfile
    */
{
  va_list args;
  char line[MAX_MESSAGE];
  char id_str[MAX_LINE];

  if (who==ALL || ringid==ROOT) {  /* controls whether output should be produced */
    /* first build output string */
    va_start(args,fmt);
    id_str[0]=0;
#ifdef PARALLEL
    if (who==ALL) sprintf(id_str," - ringID=%d",ringid);
#endif
    if (code==EC_ERROR) strcpy(line,"ERROR:");
    else if (code==EC_WARN) strcpy(line,"WARNING:");
    else if (code==EC_INFO) strcpy(line,"INFO:");
    else sprintf(line,"Error code=%d:",code);
    sprintf(line+strlen(line)," (%s:%d%s) ",fname,lineN,id_str);
    vsprintf(line+strlen(line),fmt,args);
    strcat(line,"\n");
    va_end(args);
    /* print line */
    if (code==EC_INFO) {
      printf("%s",line);
      fflush(stdout);
    }
    else if (code==EC_ERROR || code==EC_WARN) {
      fprintf(stderr,"%s",line);
      fflush(stderr);
      /* duplicate error message in logfile */
      if (logname[0]!=0) {  /* otherwise can't produce output at all */
        if (ringid==ROOT) {
          /* logfile is initialized to NULL in the beginning of the program. Hence if logfile!=NULL
             then logfile is necessarily initialized (open or allready closed)
             logfile==NULL (when logname!=0) means that error is in opening logfile itself */
          if (logfile!=NULL) {
            if (fprintf(logfile,"%s",line)==EOF) {
              fclose(logfile); /* in most cases this is redundant */
              /* try to reopen logfile and save message */
              if ((logfile=fopen(logname,"a"))!=NULL) fprintf(logfile,"%s",line);
            }
            fflush(logfile);   /* needed for warnings to appear on time */
          }
        }   /* other processors */
        else if ((logfile=fopen(logname,"a"))!=NULL) {
          fprintf(logfile,"%s",line);
          fclose(logfile);
        }
      }    /* save line to buffer to save into logfile afterwards */
      else if (code==EC_WARN) strcpy(warn_buf,line);
    }
  }
  if (code==EC_ERROR) {
    if (who==ONE && ringid!=ROOT) Synchronize();
    Stop(1);
  }
}

/*============================================================*/

void PrintError(const char *fmt, ... )
   /* print anything to stderr (on root processor) and stop;
      much simpler than LogError (do not print location of the error,
      and does not duplicate errors to files;
      assumes that all processors call it;
      no internal buffer is used to be compatible with input of any length
      (can be argv[i] in message) */
{
  va_list args;

  if (ringid==ROOT) {
    va_start(args,fmt);
    fprintf(stderr,"ERROR: ");
    vfprintf(stderr,fmt,args);
    fprintf(stderr,"\n");
    va_end(args);
    fflush(stderr);
  }
  /* wait for root to generate an error message */
  Synchronize();
  Stop(1);
}

/*============================================================*/

void LogPending(void)
   /* Logs pending warning messages (currently only one maximum).
      Should be called when logname is created and logfile is opened
      on root */
{
  if (warn_buf[0]!=0) {
    if (ringid==ROOT) fprintf(logfile,"%s",warn_buf);
    else if ((logfile=fopen(logname,"a"))!=NULL) {
      fprintf(logfile,"%s",warn_buf);
      fclose(logfile);
    }
    /* empty buffer */
    warn_buf[0]=0;
  }
}

/*============================================================*/

void PrintBoth(FILE *file,const char *fmt, ... )
   /* print anything both to file and to stdout;
      assumed that size of the message is limited to MAX_PARAGRAPH
       (i.e. no filenames in the message) */
{
  va_list args;
  char line[MAX_PARAGRAPH];

  va_start(args,fmt);
  vsprintf(line,fmt,args);
  fprintf(file,"%s",line);
  printf("%s",line);
  va_end(args);
}

/*============================================================*/

FILE *FOpenErr(const char *fname,const char *mode,const int who,const char *err_fname,
               const int lineN)
   /* open file and check for error */
{
  FILE *file;

  if ((file=fopen(fname,mode))==NULL)
    LogError(EC_ERROR,who,err_fname,lineN,"Failed to open file '%s'",fname);

  return file;
}

/*============================================================*/

void FCloseErr(FILE *file,const char *fname,const int who,const char *err_fname,
                const int lineN)
  /* close file and check for error */
{
  if (fclose(file)) LogError(EC_WARN,who,err_fname,lineN,
                             "Errors detected during work with file '%s'",fname);
}

/*============================================================*/

void RemoveErr(const char *fname,const int who,const char *err_fname,const int lineN)
   /* remove file and check the result */
{

  if(remove(fname) && errno!=ENOENT) LogError(EC_WARN,who,err_fname,lineN,
    "Failed to remove temporary file '%s'. Remove it manually, if needed",fname);
}

/*============================================================*/

void MkDirErr(const char *dir,const int who,const char *err_fname,const int lineN)
   /* make directory and check for error. */
{
#ifdef WINDOWS
  if (!CreateDirectory(dir,NULL)) {
    if (GetLastError()==ERROR_ALREADY_EXISTS)
      LogError(EC_WARN,who,err_fname,lineN,"Directory '%s' allready exists",dir);
    else LogError(EC_ERROR,who,err_fname,lineN,"Failed to make directory '%s'",dir);
  }
#elif defined(POSIX)
  if (mkdir(dir,0755)==-1) {
    if (errno==EEXIST)
      LogError(EC_WARN,who,err_fname,lineN,"Directory '%s' allready exists",dir);
    else LogError(EC_ERROR,who,err_fname,lineN,"Failed to make directory '%s'",dir);
  }
#else
  /* this should work for many cases (unknown OS), but e.g. system calls may fail
     when used in combination with MPI */
  char sbuffer[MAX_DIRSYS];

  sprintf(sbuffer,"mkdir \"%s\"",dir);
  if (system(sbuffer)) LogError(EC_WARN,who,err_fname,lineN,
                                "Failed to make directory '%s'",dir);
#endif
}

