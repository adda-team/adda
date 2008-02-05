                                   ADDA 0.78a
                                  ***********
                                "Amsterdam DDA"

                    Maxim A. Yurkin(1,2) and Alfons G. Hoekstra(1)

             (1) Faculty of Science, Section Computational Science,
                       of the University of Amsterdam,
              Kruislaan 403, 1098 SJ, Amsterdam, The Netherlands,
                   tel: +31-20-525-7530, fax: +31-20-525-7490

             (2) Institute of Chemical Kinetics and Combustion,
               Siberian Branch of the Russian Academy of Sciences,
                  Institutskaya 3, Novosibirsk 630090 Russia,
                    tel: +7-383-333-3240, fax: +7-383-334-2350

                          email: adda@science.uva.nl

                          last revised: 5 February 2008

                   Copyright (C) 2006-2008 University of Amsterdam
        This software package is covered by the GNU General Public License.


                ##             ##### ##         ##### ##         ##
             /####          /#####  /##      /#####  /##      /####
            /  ###        //    /  / ###   //    /  / ###    /  ###
               /##       /     /  /   ### /     /  /   ###      /##
              /  ##           /  /     ###     /  /     ###    /  ##
              /  ##          ## ##      ##    ## ##      ##    /  ##
             /    ##         ## ##      ##    ## ##      ##   /    ##
             /    ##         ## ##      ##    ## ##      ##   /    ##
            /      ##        ## ##      ##    ## ##      ##  /      ##
            /########        ## ##      ##    ## ##      ##  /########
           /        ##       #  ##      ##    #  ##      ## /        ##
           #        ##          /       /        /       /  #        ##
          /####      ##    /###/       /    /###/       /  /####      ##
         /   ####    ## / /   ########/    /   ########/  /   ####    ## /
        /     ##      #/ /       ####     /       ####   /     ##      #/
        #                #                #              #
         ##               ##               ##             ##


                            WINDOWS 32 EXECUTABLES
                            **********************

This Win32 package contains executables of ADDA for sequential and parallel
execution (adda.exe and adda_mpi.exe) and DLL for FFTW 3.1.2 (libfftw3-3.dll).
ADDA was compiled with MinGW 5.1.3 (gcc 3.4.5), using the default Makefile.
Parallel version was compiled linking to MPICH2 1.0.6p1, and requires it to be
installed on the system. However, it may also work with other MPI
implementations (try it at your own risk). Please download the main ADDA package
from http://www.science.uva.nl/research/scs/Software/adda/ to complement this
one. All the information, including manual and sample input files, can be found
in the main package.
