//
// File:       fft_base_kernels.h
//
// Version:    <1.0>
//
// Disclaimer: IMPORTANT:  This Apple software is supplied to you by Apple Inc. ("Apple")
//             in consideration of your agreement to the following terms, and your use,
//             installation, modification or redistribution of this Apple software
//             constitutes acceptance of these terms.  If you do not agree with these
//             terms, please do not use, install, modify or redistribute this Apple
//             software.
//
//             In consideration of your agreement to abide by the following terms, and
//             subject to these terms, Apple grants you a personal, non - exclusive
//             license, under Apple's copyrights in this original Apple software ( the
//             "Apple Software" ), to use, reproduce, modify and redistribute the Apple
//             Software, with or without modifications, in source and / or binary forms;
//             provided that if you redistribute the Apple Software in its entirety and
//             without modifications, you must retain this notice and the following text
//             and disclaimers in all such redistributions of the Apple Software. Neither
//             the name, trademarks, service marks or logos of Apple Inc. may be used to
//             endorse or promote products derived from the Apple Software without specific
//             prior written permission from Apple.  Except as expressly stated in this
//             notice, no other rights or licenses, express or implied, are granted by
//             Apple herein, including but not limited to any patent rights that may be
//             infringed by your derivative works or by other works in which the Apple
//             Software may be incorporated.
//
//             The Apple Software is provided by Apple on an "AS IS" basis.  APPLE MAKES NO
//             WARRANTIES, EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED
//             WARRANTIES OF NON - INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A
//             PARTICULAR PURPOSE, REGARDING THE APPLE SOFTWARE OR ITS USE AND OPERATION
//             ALONE OR IN COMBINATION WITH YOUR PRODUCTS.
//
//             IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL OR
//             CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//             SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//             INTERRUPTION ) ARISING IN ANY WAY OUT OF THE USE, REPRODUCTION, MODIFICATION
//             AND / OR DISTRIBUTION OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER
//             UNDER THEORY OF CONTRACT, TORT ( INCLUDING NEGLIGENCE ), STRICT LIABILITY OR
//             OTHERWISE, EVEN IF APPLE HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright ( C ) 2008 Apple Inc. All Rights Reserved.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

// Note: some changes were made compared to the original Apple FFT
// - Changed floats to doubles
// - replaced native_* by * for double precision calculations.
//

#ifndef __CL_FFT_BASE_KERNELS_
#define __CL_FFT_BASE_KERNELS_

#include <string>

using namespace std;

static string baseKernels = string(
							
						"#ifdef AMD\n"
                          "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n"
						"#else\n"
                          "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
						"#endif\n"
                       //   "#ifndef M_PI\n"
					//	  "#define M_PI 0x1.921fb54442d18p+1\n"
					//	  "#endif\n"
						  "#define complexMul(a,b) ((double2)(mad(-(a).y, (b).y, (a).x * (b).x), mad((a).y, (b).x, (a).x * (b).y)))\n"
						  "#define conj(a) ((double2)((a).x, -(a).y))\n"
						  "#define conjTransp(a) ((double2)(-(a).y, (a).x))\n"		   
						  "\n"
						  "#define fftKernel2(a,dir) \\\n"
						  "{ \\\n"
						  "    double2 c = (a)[0];    \\\n"
						  "    (a)[0] = c + (a)[1];  \\\n"
						  "    (a)[1] = c - (a)[1];  \\\n"
						  "}\n"
						  "\n"						  
						  "#define fftKernel2S(d1,d2,dir) \\\n"
						  "{ \\\n"
						  "    double2 c = (d1);   \\\n"
						  "    (d1) = c + (d2);   \\\n"
						  "    (d2) = c - (d2);   \\\n"
						  "}\n"
						  "\n"						  
						  "#define fftKernel4(a,dir) \\\n"
						  "{ \\\n"
						  "    fftKernel2S((a)[0], (a)[2], dir); \\\n"
						  "    fftKernel2S((a)[1], (a)[3], dir); \\\n"
						  "    fftKernel2S((a)[0], (a)[1], dir); \\\n"
						  "    (a)[3] = (double2)(dir)*(conjTransp((a)[3])); \\\n"
						  "    fftKernel2S((a)[2], (a)[3], dir); \\\n"
						  "    double2 c = (a)[1]; \\\n"
						  "    (a)[1] = (a)[2]; \\\n"
						  "    (a)[2] = c; \\\n"
						  "}\n"
						  "\n"						  
						  "#define fftKernel4s(a0,a1,a2,a3,dir) \\\n"
						  "{ \\\n"
						  "    fftKernel2S((a0), (a2), dir); \\\n"
						  "    fftKernel2S((a1), (a3), dir); \\\n"
						  "    fftKernel2S((a0), (a1), dir); \\\n"
						  "    (a3) = (double2)(dir)*(conjTransp((a3))); \\\n"
						  "    fftKernel2S((a2), (a3), dir); \\\n"
						  "    double2 c = (a1); \\\n"
						  "    (a1) = (a2); \\\n"
						  "    (a2) = c; \\\n" 
						  "}\n"
						  "\n"						  
						  "#define bitreverse8(a) \\\n"
						  "{ \\\n"
						  "    double2 c; \\\n"
						  "    c = (a)[1]; \\\n"
						  "    (a)[1] = (a)[4]; \\\n"
						  "    (a)[4] = c; \\\n"
						  "    c = (a)[3]; \\\n"
						  "    (a)[3] = (a)[6]; \\\n"
						  "    (a)[6] = c; \\\n"
						  "}\n"
						  "\n"						  
						  "#define fftKernel8(a,dir) \\\n"
						  "{ \\\n"
                          //0x1.6a09e6p-1=1/sqrt(2)
						  //"	const double2 w1  = (double2)(0x1.6a09e6p-1f,  dir*0x1.6a09e6p-1f);  \\\n"
						  "	const double2 w1  = (double2)(0.7071067811865475244,  dir*0.7071067811865475244);  \\\n"
						  "	const double2 w3  = (double2)(-0.7071067811865475244, dir*0.7071067811865475244);  \\\n"
						  "	double2 c; \\\n"
						  "	fftKernel2S((a)[0], (a)[4], dir); \\\n"
						  "	fftKernel2S((a)[1], (a)[5], dir); \\\n"
						  "	fftKernel2S((a)[2], (a)[6], dir); \\\n"
						  "	fftKernel2S((a)[3], (a)[7], dir); \\\n"
						  "	(a)[5] = complexMul(w1, (a)[5]); \\\n"
						  "	(a)[6] = (double2)(dir)*(conjTransp((a)[6])); \\\n"
						  "	(a)[7] = complexMul(w3, (a)[7]); \\\n"
						  "	fftKernel2S((a)[0], (a)[2], dir); \\\n"
						  "	fftKernel2S((a)[1], (a)[3], dir); \\\n"
						  "	fftKernel2S((a)[4], (a)[6], dir); \\\n"
						  "	fftKernel2S((a)[5], (a)[7], dir); \\\n"
						  "	(a)[3] = (double2)(dir)*(conjTransp((a)[3])); \\\n"
						  "	(a)[7] = (double2)(dir)*(conjTransp((a)[7])); \\\n"
						  "	fftKernel2S((a)[0], (a)[1], dir); \\\n"
						  "	fftKernel2S((a)[2], (a)[3], dir); \\\n"
						  "	fftKernel2S((a)[4], (a)[5], dir); \\\n"
						  "	fftKernel2S((a)[6], (a)[7], dir); \\\n"
						  "	bitreverse8((a)); \\\n"
						  "}\n"
						  "\n"						  
						  "#define bitreverse4x4(a) \\\n"
						  "{ \\\n"
						  "	double2 c; \\\n"
						  "	c = (a)[1];  (a)[1]  = (a)[4];  (a)[4]  = c; \\\n"
						  "	c = (a)[2];  (a)[2]  = (a)[8];  (a)[8]  = c; \\\n"
						  "	c = (a)[3];  (a)[3]  = (a)[12]; (a)[12] = c; \\\n"
						  "	c = (a)[6];  (a)[6]  = (a)[9];  (a)[9]  = c; \\\n"
						  "	c = (a)[7];  (a)[7]  = (a)[13]; (a)[13] = c; \\\n"
						  "	c = (a)[11]; (a)[11] = (a)[14]; (a)[14] = c; \\\n"
						  "}\n"
						  "\n"						  
						  "#define fftKernel16(a,dir) \\\n"
						  "{ \\\n"
						  //"    const double w0 = 0x1.d906bcp-1f; \\\n"
						  "    const double w0 = 0.92387953251128675613; \\\n"//sin(67.5)
						  "    const double w1 = 0.38268343236508977173; \\\n"//sin(22.5)
						  "    const double w2 = 0.7071067811865475244; \\\n"//1/sqrt(2)
						  "    fftKernel4s((a)[0], (a)[4], (a)[8],  (a)[12], dir); \\\n"
						  "    fftKernel4s((a)[1], (a)[5], (a)[9],  (a)[13], dir); \\\n"
						  "    fftKernel4s((a)[2], (a)[6], (a)[10], (a)[14], dir); \\\n"
						  "    fftKernel4s((a)[3], (a)[7], (a)[11], (a)[15], dir); \\\n"
						  "    (a)[5]  = complexMul((a)[5], (double2)(w0, dir*w1)); \\\n"
						  "    (a)[6]  = complexMul((a)[6], (double2)(w2, dir*w2)); \\\n"
						  "    (a)[7]  = complexMul((a)[7], (double2)(w1, dir*w0)); \\\n"
						  "    (a)[9]  = complexMul((a)[9], (double2)(w2, dir*w2)); \\\n"
						  "    (a)[10] = (double2)(dir)*(conjTransp((a)[10])); \\\n"
						  "    (a)[11] = complexMul((a)[11], (double2)(-w2, dir*w2)); \\\n"
						  "    (a)[13] = complexMul((a)[13], (double2)(w1, dir*w0)); \\\n"
						  "    (a)[14] = complexMul((a)[14], (double2)(-w2, dir*w2)); \\\n"
						  "    (a)[15] = complexMul((a)[15], (double2)(-w0, dir*-w1)); \\\n"
						  "    fftKernel4((a), dir); \\\n"
						  "    fftKernel4((a) + 4, dir); \\\n"
						  "    fftKernel4((a) + 8, dir); \\\n"
						  "    fftKernel4((a) + 12, dir); \\\n"
						  "    bitreverse4x4((a)); \\\n"
						  "}\n"
						  "\n"						  
						  "#define bitreverse32(a) \\\n"
						  "{ \\\n"
						  "    double2 c1, c2; \\\n"
						  "    c1 = (a)[2];   (a)[2] = (a)[1];   c2 = (a)[4];   (a)[4] = c1;   c1 = (a)[8];   (a)[8] = c2;    c2 = (a)[16];  (a)[16] = c1;   (a)[1] = c2; \\\n"
						  "    c1 = (a)[6];   (a)[6] = (a)[3];   c2 = (a)[12];  (a)[12] = c1;  c1 = (a)[24];  (a)[24] = c2;   c2 = (a)[17];  (a)[17] = c1;   (a)[3] = c2; \\\n"
						  "    c1 = (a)[10];  (a)[10] = (a)[5];  c2 = (a)[20];  (a)[20] = c1;  c1 = (a)[9];   (a)[9] = c2;    c2 = (a)[18];  (a)[18] = c1;   (a)[5] = c2; \\\n"
						  "    c1 = (a)[14];  (a)[14] = (a)[7];  c2 = (a)[28];  (a)[28] = c1;  c1 = (a)[25];  (a)[25] = c2;   c2 = (a)[19];  (a)[19] = c1;   (a)[7] = c2; \\\n"
						  "    c1 = (a)[22];  (a)[22] = (a)[11]; c2 = (a)[13];  (a)[13] = c1;  c1 = (a)[26];  (a)[26] = c2;   c2 = (a)[21];  (a)[21] = c1;   (a)[11] = c2; \\\n"
						  "    c1 = (a)[30];  (a)[30] = (a)[15]; c2 = (a)[29];  (a)[29] = c1;  c1 = (a)[27];  (a)[27] = c2;   c2 = (a)[23];  (a)[23] = c1;   (a)[15] = c2; \\\n"
						  "}\n"
						  "\n"						  
						  "#define fftKernel32(a,dir) \\\n"
						  "{ \\\n"
						  "    fftKernel2S((a)[0],  (a)[16], dir); \\\n"
						  "    fftKernel2S((a)[1],  (a)[17], dir); \\\n"
						  "    fftKernel2S((a)[2],  (a)[18], dir); \\\n"
						  "    fftKernel2S((a)[3],  (a)[19], dir); \\\n"
						  "    fftKernel2S((a)[4],  (a)[20], dir); \\\n"
						  "    fftKernel2S((a)[5],  (a)[21], dir); \\\n"
						  "    fftKernel2S((a)[6],  (a)[22], dir); \\\n"
						  "    fftKernel2S((a)[7],  (a)[23], dir); \\\n"
						  "    fftKernel2S((a)[8],  (a)[24], dir); \\\n"
						  "    fftKernel2S((a)[9],  (a)[25], dir); \\\n"
						  "    fftKernel2S((a)[10], (a)[26], dir); \\\n"
						  "    fftKernel2S((a)[11], (a)[27], dir); \\\n"
						  "    fftKernel2S((a)[12], (a)[28], dir); \\\n"
						  "    fftKernel2S((a)[13], (a)[29], dir); \\\n"
						  "    fftKernel2S((a)[14], (a)[30], dir); \\\n"
						  "    fftKernel2S((a)[15], (a)[31], dir); \\\n"
						  "    (a)[17] = complexMul((a)[17], (double2)(0.19509032201612826785, dir*0.98078528040323044913)); \\\n"//sin(11.25),sin(78.75)
						  "    (a)[18] = complexMul((a)[18], (double2)(0.92387953251128675613, dir*0.38268343236508977173)); \\\n"//sin(67.5),sin(22.5)
						  "    (a)[19] = complexMul((a)[19], (double2)(0.83146961230254523708, dir*0.55557023301960222474)); \\\n"//sin(56.25),sin(33.75)
						  "    (a)[20] = complexMul((a)[20], (double2)(0.7071067811865475244, dir*0.7071067811865475244)); \\\n"//sin(45),sin(45)
						  "    (a)[21] = complexMul((a)[21], (double2)(0.55557023301960222474, dir*0.83146961230254523708)); \\\n"//sin(33.75),sin(56.25)
						  "    (a)[22] = complexMul((a)[22], (double2)(0.38268343236508977173, dir*0.92387953251128675613)); \\\n"//sin(22.5),sin(67.5)
						  "    (a)[23] = complexMul((a)[23], (double2)(0.98078528040323044913, dir*0.19509032201612826785)); \\\n"//sin(78.75),sin(11.25)
						  "    (a)[24] = complexMul((a)[24], (double2)(0.00000000000000000000, dir*1.00000000000000000000)); \\\n"//sin(0),sin(90)
						  "    (a)[25] = complexMul((a)[25], (double2)(-0.98078528040323044913, dir*0.19509032201612826785)); \\\n"
						  "    (a)[26] = complexMul((a)[26], (double2)(-0.38268343236508977173, dir*0.92387953251128675613)); \\\n"
						  "    (a)[27] = complexMul((a)[27], (double2)(-0.55557023301960222474, dir*0.83146961230254523708)); \\\n"
						  "    (a)[28] = complexMul((a)[28], (double2)(-0.7071067811865475244, dir*0.7071067811865475244)); \\\n"
						  "    (a)[29] = complexMul((a)[29], (double2)(-0.83146961230254523708, dir*0.55557023301960222474)); \\\n"
						  "    (a)[30] = complexMul((a)[30], (double2)(-0.92387953251128675613, dir*0.38268343236508977173)); \\\n"
						  "    (a)[31] = complexMul((a)[31], (double2)(-0.19509032201612826785, dir*0.98078528040323044913)); \\\n"
						  "    fftKernel16((a), dir); \\\n"
						  "    fftKernel16((a) + 16, dir); \\\n"
						  "    bitreverse32((a)); \\\n"
						  "}\n\n"
						  );

static string twistKernelInterleaved = string(
											  "__kernel void \\\n"
											  "clFFT_1DTwistInterleaved(__global double2 *in, unsigned int startRow, unsigned int numCols, unsigned int N, unsigned int numRowsToProcess, int dir) \\\n"
											  "{ \\\n"
											  "   double2 a, w; \\\n"
											  "   double ang; \\\n"
											  "   unsigned int j; \\\n"
											  "	unsigned int i = get_global_id(0); \\\n"
											  "	unsigned int startIndex = i; \\\n"
											  "	 \\\n"
											  "	if(i < numCols) \\\n"
											  "	{ \\\n"
											  "	    for(j = 0; j < numRowsToProcess; j++) \\\n"
											  "	    { \\\n"
											  "	        a = in[startIndex]; \\\n"
											  "	        ang = 2.0f * M_PI * dir * i * (startRow + j) / N; \\\n"
											  "	        w = (double2)(cos(ang), sin(ang)); \\\n"
											  "	        a = complexMul(a, w); \\\n"
											  "	        in[startIndex] = a; \\\n"
											  "	        startIndex += numCols; \\\n"
											  "	    } \\\n"
											  "	}	 \\\n"
											  "} \\\n"
											  );

static string twistKernelPlannar = string(
										  "__kernel void \\\n"
										  "clFFT_1DTwistSplit(__global double *in_real, __global double *in_imag , unsigned int startRow, unsigned int numCols, unsigned int N, unsigned int numRowsToProcess, int dir) \\\n"
										  "{ \\\n"
										  "    double2 a, w; \\\n"
										  "    double ang; \\\n"
										  "    unsigned int j; \\\n"
										  "	unsigned int i = get_global_id(0); \\\n"
										  "	unsigned int startIndex = i; \\\n"
										  "	 \\\n"
										  "	if(i < numCols) \\\n"
										  "	{ \\\n"
										  "	    for(j = 0; j < numRowsToProcess; j++) \\\n"
										  "	    { \\\n"
										  "	        a = (double2)(in_real[startIndex], in_imag[startIndex]); \\\n"
										  "	        ang = 2.0f * M_PI * dir * i * (startRow + j) / N; \\\n"
										  "	        w = (double2)(cos(ang), sin(ang)); \\\n"
										  "	        a = complexMul(a, w); \\\n"
										  "	        in_real[startIndex] = a.x; \\\n"
										  "	        in_imag[startIndex] = a.y; \\\n"
										  "	        startIndex += numCols; \\\n"
										  "	    } \\\n"
										  "	}	 \\\n"
										  "} \\\n"
										  );										  



#endif
