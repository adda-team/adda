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

#ifndef DEBUG
#  ifdef DEBUG_ALL
#    define DEBUG
#  endif
#endif

#ifdef DEBUG 
#  define D(p)        {extern int RingId;\
                       if(!RingId){debug_printf(__FILE__, __LINE__, p);}}
#  define D2(p,a)     {extern int RingId;\
                       if(!RingId){debug_printf(__FILE__, __LINE__, p,a);}}
#  define D3(p,a,b)   {extern int RingId;\
                       if(!RingId){debug_printf(__FILE__, __LINE__, p,a,b);}}
#  define D4(p,a,b,c) {extern int RingId;\
                       if(!RingId){debug_printf(__FILE__, __LINE__, p,a,b,c);}}
#else
#  define D(p)
#  define D2(p,a)
#  define D3(p,a,b)
#  define D4(p,a,b,c)
#  define S2(p,a)
#endif
