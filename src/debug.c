/* FILE: debug.h
 * AUTH: vesseur
 * VERS: 1.1
 * DATE: Mon Aug 16 11:45:30 MET DST 1993
 * HIST: 1.1 creation
 * BUGS: 
 * $Log$
 * Revision 1.1  2005/04/07 09:24:17  myurkin
 * Initial revision
 *
 * Revision 1.3  1993/10/13  13:14:17  vesseur
 * *** empty log message ***
 *
 * Revision 1.2  1993/08/19  10:12:54  vesseur
 * *** empty log message ***
 *
 * Revision 1.1  1993/08/18  14:48:22  vesseur
 * Initial revision
 *
 * Revision 1.1  1993/08/18  14:48:22  vesseur
 * Initial revision
 *
 */

#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

void
debug_printf (char *fname, int line, char *fmt, ... )
{ va_list args;
  char    pr_line[255];
  
  va_start(args, fmt);
  sprintf(pr_line, "%s:%d: ", fname, line);
  vsprintf(pr_line+strlen(pr_line), fmt, args);
  strcat(pr_line, "\n");
/*  assert(isasync(stderr)); */
  printf(pr_line);
  va_end(args);
}
