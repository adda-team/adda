/* File: sparse_ops.h
 * Descr: low-level routines for the sparse mode ADDA; void in non-sparse mode
 *
 * Copyright (C) 2011-2013 ADDA contributors
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
#ifdef SPARSE

#ifndef __sparse_ops_h
#define __sparse_ops_h

#include "cmplx.h"
#include "interaction.h"
#include "vars.h"

#ifdef USE_SSE3

//=====================================================================================================================

static inline void CcMul(doublecomplex * restrict argvec_src,doublecomplex * restrict argvec_dest, const size_t j)
// Takes the j'th block in argvec_src, multiplies by cc_sqrt at that position and stores the result in argvec_dest.
{
	const size_t j3=j*3;
	*(__m128d *)&(argvec_dest[j3]) = cmul(*(__m128d *)&(argvec_src[j3]),
		*(__m128d *)&(cc_sqrt[material[j]][0]));
	*(__m128d *)&(argvec_dest[j3+1]) = cmul(*(__m128d *)&(argvec_src[j3+1]),
		*(__m128d *)&(cc_sqrt[material[j]][1]));
	*(__m128d *)&(argvec_dest[j3+2]) = cmul(*(__m128d *)&(argvec_src[j3+2]),
		*(__m128d *)&(cc_sqrt[material[j]][2]));
}

//=====================================================================================================================

static inline void AijProd(doublecomplex * restrict argvec,doublecomplex * restrict resultvec,const size_t i,
	const size_t j)
/* Handles the multiplication of the i'th block of argvec with the G_ij block of the G-matrix, and adds the result
 * to the j'th block of resultvec.
 */
{	
	doublecomplex iterm[6];
	const size_t i3 = 3*i, j3 = 3*j;

	(*CalcInterTerm)(position[i3]-position_full[j3], position[i3+1]-position_full[j3+1],
		position[i3+2]-position_full[j3+2], iterm);

	__m128d res, tmp;

	const __m128d argX = _mm_load_pd((double *)(argvec+j3));
	const __m128d argY = _mm_load_pd((double *)(argvec+j3+1));
	const __m128d argZ = _mm_load_pd((double *)(argvec+j3+2));

	res = cmul(argX, *(__m128d *)&(iterm[0]));
	tmp = cmul(argY, *(__m128d *)&(iterm[1]));
	res = cadd(tmp,res);
	tmp = cmul(argZ, *(__m128d *)&(iterm[2]));
	res = cadd(tmp,res);
	*(__m128d *)&(resultvec[i3]) = cadd(res, *(__m128d *)&(resultvec[i3]));

	res = cmul(argX, *(__m128d *)&(iterm[1]));
	tmp = cmul(argY, *(__m128d *)&(iterm[3]));
	res = cadd(tmp,res);
	tmp = cmul(argZ, *(__m128d *)&(iterm[4]));
	res = cadd(tmp,res);
	*(__m128d *)&(resultvec[i3+1]) = cadd(res, *(__m128d *)&(resultvec[i3+1]));

	res = cmul(argX, *(__m128d *)&(iterm[2]));
	tmp = cmul(argY, *(__m128d *)&(iterm[4]));
	res = cadd(tmp,res);
	tmp = cmul(argZ, *(__m128d *)&(iterm[5]));
	res = cadd(tmp,res);
	*(__m128d *)&(resultvec[i3+2]) = cadd(res, *(__m128d *)&(resultvec[i3+2]));
}

//=====================================================================================================================

static inline void DiagProd(doublecomplex * restrict argvec,doublecomplex * restrict resultvec,const size_t i)
/* Multiplies the result in the i'th block of resultvec by cc_sqrt at that block, subtracts the result from the i'th
 * block of argvec, and stores the result in the i'th block if resultvec.
 */
{
	const size_t i3 = i*3;

	const __m128d tmp1 = cmul(*(__m128d *)&(resultvec[i3]),*(__m128d *)&(cc_sqrt[material[i]][0]));
	const __m128d tmp2 = cmul(*(__m128d *)&(resultvec[i3+1]),*(__m128d *)&(cc_sqrt[material[i]][1]));
	const __m128d tmp3 = cmul(*(__m128d *)&(resultvec[i3+2]),*(__m128d *)&(cc_sqrt[material[i]][2]));

	*(__m128d *)&(resultvec[i3]) = _mm_sub_pd(*(__m128d *)&(argvec[i3]),tmp1);
	*(__m128d *)&(resultvec[i3+1]) = _mm_sub_pd(*(__m128d *)&(argvec[i3+1]),tmp2);
	*(__m128d *)&(resultvec[i3+2]) = _mm_sub_pd(*(__m128d *)&(argvec[i3+2]),tmp3);
}

#else //SSE3 not used

//=====================================================================================================================

static inline void CcMul(doublecomplex * restrict argvec_src,doublecomplex * restrict argvec_dest, const size_t j)
// Takes the j'th block in argvec_src, multiplies by cc_sqrt at that position and stores the result in argvec_dest.
{
	const size_t j3 = j*3;
	cMult(argvec_src[j3],cc_sqrt[material[j]][0],argvec_dest[j3]);
	cMult(argvec_src[j3+1],cc_sqrt[material[j]][1],argvec_dest[j3+1]);
	cMult(argvec_src[j3+2],cc_sqrt[material[j]][2],argvec_dest[j3+2]);
}

//=====================================================================================================================

static inline void AijProd(doublecomplex * restrict argvec,doublecomplex * restrict resultvec,const size_t i,
	const size_t j)
/* Handles the multiplication of the i'th block of argvec with the G_ij block of the G-matrix, and adds the result
 * to the j'th block of resultvec.
 */
{
	doublecomplex tmp1,resX,resY,resZ;
	doublecomplex iterm[6];
	const size_t i3=3*i,j3=3*j;

	//D("%d %d %d %d %d %d %d %d",i,j,position[3*i],position[3*i+1],position[3*i+2],
	//	position[3*j],position[3*j+1],position[3*j+2]);
	(*CalcInterTerm)(position[i3]-position_full[j3],position[i3+1]-position_full[j3+1],
		position[i3+2]-position_full[j3+2],iterm);

	cMult(argvec[j3],iterm[0],resX);
	cMult(argvec[j3+1],iterm[1],tmp1);
	cAdd(resX,tmp1,resX);
	cMult(argvec[j3+2],iterm[2],tmp1);
	cAdd(resX,tmp1,resX);

	cMult(argvec[j3],iterm[1],resY);
	cMult(argvec[j3+1],iterm[3],tmp1);
	cAdd(resY,tmp1,resY);
	cMult(argvec[j3+2],iterm[4],tmp1);
	cAdd(resY,tmp1,resY);

	cMult(argvec[j3],iterm[2],resZ);
	cMult(argvec[j3+1],iterm[4],tmp1);
	cAdd(resZ,tmp1,resZ);
	cMult(argvec[j3+2],iterm[5],tmp1);
	cAdd(resZ,tmp1,resZ);

	cAdd(resX,resultvec[i3],resultvec[i3]);
	cAdd(resY,resultvec[i3+1],resultvec[i3+1]);
	cAdd(resZ,resultvec[i3+2],resultvec[i3+2]);
}

//=====================================================================================================================

static inline void DiagProd(doublecomplex * restrict argvec,doublecomplex * restrict resultvec,const size_t i)
/* Multiplies the result in the i'th block of resultvec by cc_sqrt at that block, subtracts the result from the i'th
 * block of argvec, and stores the result in the i'th block if resultvec.
 */
{
	const size_t i3 = i*3;
	doublecomplex tmp1, tmp2, tmp3;

	cMult(resultvec[i3],cc_sqrt[material[i]][0],tmp1);
	cMult(resultvec[i3+1],cc_sqrt[material[i]][1],tmp2);
	cMult(resultvec[i3+2],cc_sqrt[material[i]][2],tmp3);
	cSubtr(argvec[i3], tmp1, resultvec[i3]);
	cSubtr(argvec[i3+1], tmp2, resultvec[i3+1]);
	cSubtr(argvec[i3+2], tmp3, resultvec[i3+2]);
}

#endif // !USE_SSE3

#endif // __sparse_ops_h

#endif // SPARSE
