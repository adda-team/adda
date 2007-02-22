/* Michel Grimminck */

#define false        0
#define true         1
#define PI           3.14159265359
#define UNDEF       -1.234125626E-30

/* object types */
#define SPHERE        0
#define BOX           1
#define ELLIPS        2
#define TRIANGLE      3
#define PRISMA        3
#define LINE          4
#define COATED        5
#define SPHEREBOX     6
#define DISK	      7
#define RBC	      8
#define SDISK	      9
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


/* which way to calculate coupleconstant */
#define CM           0
#define LDR          1
#define RADCOR       2

/* ldr constants */
#define LDR_B1      -1.891531
#define LDR_B2       0.1648469
#define LDR_B3      -1.7700004

#define NORMAL       0
#define PAR_AND_PER  1

#define real         0
#define imag         1

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

/**GD************   Global Defines and Data structures   *****************/
 
#define EC_MASK  0xF0000000
#define EC_FATAL 0xE0000000
#define EC_CRIT  0xD0000000
#define EC_ERROR 0xC0000000
#define EC_WARN  0xB0000000
#define EC_DEBUG 0xA0000000
#define EC_INFO  0x90000000
#define EC_MESS  0x80000000

