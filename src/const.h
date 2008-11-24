#ifndef __const_h
#define __const_h

#define false        0
#define true         1
#define PI           3.141592653589793
#define TWO_PI       6.283185307179586
#define UNDEF        -1

/* object types */
#define SPHERE        0
#define BOX           1
#define TRIANGLE      3
#define PRISMA        3
#define LINE          4
#define COATED        5
#define SPHEREBOX     6
#define DISK	      7
#define RBC	      8
#define AGGREGATE     10
#define ELLIPSOIDAL   11
#define RBC_ROT       12
#define STICK         13
#define SDISK_ROT     14
#define SPHEROID_ROT  15
#define LYMPHOCYTE1   16
#define LEUCOCYTE2    17
#define CYLINDER      18
#define LYMPHOCYTE2   19
#define READ          20

/* which way to calculate coupleconstant */
#define CM           0
#define RADCOR       1
#define LDR          2
#define CLDR	     3
#define SOrd	     4

/* how to calculate scattering quantities */
#define DRAINE       0
     /* SOrd defined above */

/* how to calculate interaction term */
#define POINT        0
     /* SOrd defined above */

/* ldr constants */
#define LDR_B1       1.8915316
#define LDR_B2      -0.1648469
#define LDR_B3       1.7700004

/* 2nd_order constants */
#define SO_B1        1.5867182
#define SO_B2        0.13488017
#define SO_B3        0.11895826

#define NORMAL       0
#define PAR_AND_PER  1

#define real         0
#define imag         1

/* path and size of tables */
#define TAB_PATH     "tables"
#define TAB_SIZE     142
#define TAB_RMAX     10

/* beam types */
#define PLANE        0
#define LMINUS       1
#define DAVIS1       2
#define DAVIS3       3
#define DAVIS5       4
#define BARTON1      5
#define BARTON3      6
#define BARTON5      7
#define WEIRD        8
#define BUGGY        9

/**GD************   Global Defines and Data structures (all for LogError) *****************/
 
#define POSIT __FILE__, __LINE__

#define ALL 0   
#define ONE 1

#define EC_MASK  0xF0000000
#define EC_FATAL 0xE0000000
#define EC_CRIT  0xD0000000
#define EC_ERROR 0xC0000000
#define EC_WARN  0xB0000000
#define EC_DEBUG 0xA0000000
#define EC_INFO  0x90000000
#define EC_MESS  0x80000000

#endif /*__const_h*/