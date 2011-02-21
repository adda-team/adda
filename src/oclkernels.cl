/* File: oclkernels.cl
 * $Date::                            $
 * Descr: Kernel File for OpenCL kernels.
 *        includes all subfunctions of the Matvec routine as OpenCL kernels
 *
 * Copyright (C) 2010-2011 ADDA contributors
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

#ifdef AMD
#	pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#	pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

/*************************************
 *  functions used in kernels        *
 *                                   *
 *************************************/

void cMult(__global double2 *a,__global double2 *b,__global double2 *c)
// complex multiplication; c=ab; !!! c should be different from a and b !!!
{
	(*c).s0=(*a).s0*(*b).s0 - (*a).s1*(*b).s1;
	(*c).s1=(*a).s1*(*b).s0 + (*a).s0*(*b).s1;
}

//============================================================

void cMult2(__global double2 *a,__global double2 *b,__private double2 *c)
// complex multiplication; c=ab; !!! c should be different from a and b ! c is local here in cMult2
{
	(*c).s0=(*a).s0*(*b).s0 - (*a).s1*(*b).s1;
	(*c).s1=(*a).s1*(*b).s0 + (*a).s0*(*b).s1;
}

//============================================================

double cvNorm2(__global double2 *a)
// square of the norm of a complex vector[3]
{
	return ( a[0].s0*a[0].s0 + a[0].s1*a[0].s1 + a[1].s0*a[1].s0 + a[1].s1*a[1].s1
		+ a[2].s0*a[2].s0 + a[2].s1*a[2].s1 );
}

//============================================================

void cAdd(__global double2 *a,__private double2 *b,__global double2 *c)
// add two complex numbers; c=a+b
{
	(*c).s0 = (*a).s0 + (*b).s0;
	(*c).s1 = (*a).s1 + (*b).s1;
}

/*************************************
 *  Arith1 kernels                   *
 *                                   *
 *************************************/

__kernel void nConj(__global double2 *a)
{
	__private size_t id=get_global_id(0);

	a[id].s1= -a[id].s1;
}

//============================================================

__kernel void clzero(__global double2 *input)
{
	__private size_t id = get_global_id(0);

	input[id] = 0.0;
}

//============================================================

__kernel void arith1(__global unsigned char *material,__global unsigned short *position,
	__global double2 *cc_sqrt /* 15x3 */, __global double2 *argvec, __global double2 *Xmatrix,
	long local_Nsmall,long smallY,long gridX)
{
	__private size_t id=get_global_id(0);
	__private size_t j=3*id;
	__private char mat=material[id];
	__private size_t index;
	__private int xcomp;

	index = ((position[j+2]*smallY+position[j+1])*gridX+position[j]);
	for (xcomp=0;xcomp<3;xcomp++)
		cMult(&cc_sqrt[mat*15+xcomp],&argvec[j+xcomp],&Xmatrix[index+xcomp*local_Nsmall]);
}

/*************************************
 *  Arith2 kernel                    *
 *                                   *
 *************************************/

__kernel void arith2(__global double2 *Xmatrix,__global double2 *slices,long gridZ,long smallY,
	long gridX,long gridYZ,long local_Nsmall,long x)
{
	__private size_t z=get_global_id(0);
	__private size_t y=get_global_id(1);
	__private size_t i;
	__private size_t j;
	__private int xcomp;

	i = y*gridZ+z;
	j = (z*smallY+y)*gridX+x;
	for (xcomp=0;xcomp<3;xcomp++) slices[i+xcomp*gridYZ]=Xmatrix[j+xcomp*local_Nsmall];
}

/*************************************
 *  Arith3 kernels and functions     *
 *                                   *
 *************************************/

void cSymMatrVec(__private double2 *matr,__private double2 *vec,__private double2 *res)
// multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec
{
	res[0].s0 = matr[0].s0*vec[0].s0 - matr[0].s1*vec[0].s1
	          + matr[1].s0*vec[1].s0 - matr[1].s1*vec[1].s1
	          + matr[2].s0*vec[2].s0 - matr[2].s1*vec[2].s1;
	res[0].s1 = matr[0].s0*vec[0].s1 + matr[0].s1*vec[0].s0
	          + matr[1].s0*vec[1].s1 + matr[1].s1*vec[1].s0
	          + matr[2].s0*vec[2].s1 + matr[2].s1*vec[2].s0;

	res[1].s0 = matr[1].s0*vec[0].s0 - matr[1].s1*vec[0].s1
	          + matr[3].s0*vec[1].s0 - matr[3].s1*vec[1].s1
	          + matr[4].s0*vec[2].s0 - matr[4].s1*vec[2].s1;
	res[1].s1 = matr[1].s0*vec[0].s1 + matr[1].s1*vec[0].s0
	          + matr[3].s0*vec[1].s1 + matr[3].s1*vec[1].s0
	          + matr[4].s0*vec[2].s1 + matr[4].s1*vec[2].s0;

	res[2].s0 = matr[2].s0*vec[0].s0 - matr[2].s1*vec[0].s1
	          + matr[4].s0*vec[1].s0 - matr[4].s1*vec[1].s1
	          + matr[5].s0*vec[2].s0 - matr[5].s1*vec[2].s1;
	res[2].s1 = matr[2].s0*vec[0].s1 + matr[2].s1*vec[0].s0
	          + matr[4].s0*vec[1].s1 + matr[4].s1*vec[1].s0
	          + matr[5].s0*vec[2].s1 + matr[5].s1*vec[2].s0;
}

//============================================================

__kernel void arith3(__global double2 *slices_tr,__global double2 *Dmatrix,long local_x0,
	long smallY,long smallZ,long gridX,long DsizeY,long DsizeZ,long NDCOMP,char reduced_FFT,
	char transposed,long x)
{
	size_t z = get_global_id(0);
	size_t y = get_global_id(1);
	size_t gridZ = get_global_size(0);
	size_t gridY = get_global_size(1);
	__private double2 xv[3];
	__private double2 yv[3];
	__private double2 fmat[6];
	__private int Xcomp;
	__private size_t i=z *gridY + y;// indexSliceZY
	size_t j;
	size_t xa=x;
	size_t ya=y;
	size_t za=z;

	// works, because of the double2 vector type
	for (Xcomp=0;Xcomp<3;Xcomp++) xv[Xcomp]=slices_tr[i+Xcomp*gridY*gridZ];

	if (transposed==1) { // used only for G_SO
		if (x>0) xa=gridX-x;
		if (y>0) ya=gridY-y;
		if (z>0) za=gridZ-z;
	}
	else {
		if (y>=DsizeY) ya=gridY-y;
		if (z>=DsizeZ) za=gridZ-z;
	}
	j=NDCOMP*(xa*DsizeY*DsizeZ+za*DsizeY+ya);

	for (int f=0;f<6;f++) fmat[f]=Dmatrix[j+f]; // maybe could be done more efficient TODO

	if (reduced_FFT==1) {
		if (y>smallY) {
			fmat[1]*=-1;
			if (z>smallZ) fmat[2]*=-1;
			else fmat[4]*=-1;
		}
		else if (z>smallZ) {
			fmat[2]*=-1;
			fmat[4]*=-1;
		}
	}
	cSymMatrVec(fmat,xv,yv); // yv=fmat*xv
	for (Xcomp=0;Xcomp<3;Xcomp++) slices_tr[i+Xcomp*gridY*gridZ]=yv[Xcomp];
}

/*************************************
 *  Arith4 kernel                    *
 *                                   *
 *************************************/

__kernel void arith4(__global double2 *Xmatrix,__global double2 *slices,long gridZ,long smallY,
	long gridX,long gridYZ,long local_Nsmall,long x)
{
	__private size_t z =get_global_id(0);
	__private size_t y =get_global_id(1);
	__private size_t i;
	__private size_t j;
	__private int xcomp;

	i = y*gridZ+z;
	j = (z*smallY+y)*gridX+x;
	for (xcomp=0;xcomp<3;xcomp++) Xmatrix[j+xcomp*local_Nsmall]=slices[i+xcomp*gridYZ];
}

/*************************************
 *  Arith5 kernels                   *
 *                                   *
 *************************************/

__kernel void arith5(__global unsigned char *material,__global unsigned short *position,
	__global double2 *cc_sqrt,__global double2 *argvec,__global double2 *Xmatrix,long local_Nsmall,
	long smallY,long gridX,__global double2 *resultvec)
{
	__private size_t id = get_global_id(0);
	__private size_t j=3*id;
	__private char mat = material[id];
	__private size_t index;
	__private double2 temp;
	__private int xcomp;

	index = ((position[j+2]*smallY+position[j+1])*gridX+position[j]);
	for (xcomp=0;xcomp<3;xcomp++) {
		cMult2(&cc_sqrt[mat*15+xcomp],&Xmatrix[index+xcomp*local_Nsmall],&temp);
		resultvec[j+xcomp]=argvec[j+xcomp]+temp;
	}
}

//============================================================

__kernel void inpr(__global double *inprod, __global double2 *resultvec)
{
	__private size_t id = get_global_id(0);

	inprod[id]=cvNorm2(resultvec+(id*3));
}

/*************************************
 * transpose kernels                 *
 *                                   *
 *************************************/

__kernel void transpose(__global double2 *input,__global double2 *output,long width,long height)
{
	size_t idz = get_global_id(0);
	size_t idy = get_global_id(1);
	size_t wth = width*height;

	for (int k=0;k<3;k++) output[idy*height+idz+k*wth]=input[idz*width+idy+k*wth];
}

//optimised transpose kernel with cache and
//removed bank conflicts obtained from Nvidia SDK samples
#define BLOCK_DIM 16
__kernel void transposeo(
        __global double2 *idata,
        __global double2 *odata,
        long width,
        long height,
        __local double2 *block
    )
{
    // read tiles into local memory
    unsigned int xIndex = get_global_id(0);
    unsigned int yIndex = get_global_id(1);
    unsigned int zIndex = get_global_id(2);
    int htw = height*width;
    if((xIndex < width) && (yIndex < height))
    {
    unsigned int index_in = yIndex * width + xIndex + htw * zIndex;
    block[get_local_id(1)*(BLOCK_DIM+1)+get_local_id(0)] = idata[index_in];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    // write transposed tile back to global memory
    xIndex = get_group_id(1) * BLOCK_DIM + get_local_id(0);
    yIndex = get_group_id(0) * BLOCK_DIM + get_local_id(1);
    if((xIndex < height) && (yIndex < width))
    {
    unsigned int index_out = yIndex * height + xIndex + htw * zIndex;
    odata[index_out] = block[get_local_id(0)*(BLOCK_DIM+1)+get_local_id(1)];
    }
}
