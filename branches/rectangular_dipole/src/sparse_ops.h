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
/* Handles the multiplication of the j'th block of argvec with the G_ij block of the G-matrix, and adds the result
 * to the i'th block of resultvec.
 */
{	
	doublecomplex iterm[6];
	const size_t i3 = 3*i, j3 = 3*j;

	__m128d res, tmp;

	IGNORE_WARNING(-Wstrict-aliasing); // cast from doublecomplex* to double* is perfectly valid in C99
	const __m128d argX = _mm_load_pd((double *)(argvec+j3));
	const __m128d argY = _mm_load_pd((double *)(argvec+j3+1));
	const __m128d argZ = _mm_load_pd((double *)(argvec+j3+2));
	STOP_IGNORE;

	if (j!=local_nvoid_d0+i) { // main interaction is not computed for coinciding dipoles
		(*InterTerm_int)(position[i3]-position_full[j3], position[i3+1]-position_full[j3+1],
			position[i3+2]-position_full[j3+2], iterm);
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
	if (surface) { // surface interaction is computed always
		(*ReflTerm_int)(position[i3]-position_full[j3], position[i3+1]-position_full[j3+1],
			position[i3+2]+position_full[j3+2], iterm);
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

		res = cmul(argZ, *(__m128d *)&(iterm[5]));
		tmp = cmul(argX, *(__m128d *)&(iterm[2]));
		res=_mm_sub_pd(res,tmp);
		tmp = cmul(argY, *(__m128d *)&(iterm[4]));
		res=_mm_sub_pd(res,tmp);
		*(__m128d *)&(resultvec[i3+2]) = cadd(res, *(__m128d *)&(resultvec[i3+2]));
	}
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

	argvec_dest[j3]=argvec_src[j3]*cc_sqrt[material[j]][0];
	argvec_dest[j3+1]=argvec_src[j3+1]*cc_sqrt[material[j]][1];
	argvec_dest[j3+2]=argvec_src[j3+2]*cc_sqrt[material[j]][2];
}

//=====================================================================================================================

static inline void AijProd(doublecomplex * restrict argvec,doublecomplex * restrict resultvec,const size_t i,
	const size_t j)
/* Handles the multiplication of the j'th block of argvec with the G_ij block of the G-matrix, and adds the result
 * to the i'th block of resultvec.
 */
{
	doublecomplex res[3];
	doublecomplex iterm[6];
	const size_t i3=3*i,j3=3*j;

	if (j!=local_nvoid_d0+i) { // main interaction is not computed for coinciding dipoles
		(*InterTerm_int)(position[i3]-position_full[j3],position[i3+1]-position_full[j3+1],
			position[i3+2]-position_full[j3+2],iterm);
		cSymMatrVec(iterm,argvec+j3,res);
		cvAdd(res,resultvec+i3,resultvec+i3);
	}
	if (surface) { // surface interaction is computed always
		(*ReflTerm_int)(position[i3]-position_full[j3],position[i3+1]-position_full[j3+1],
			position[i3+2]+position_full[j3+2],iterm);
		cReflMatrVec(iterm,argvec+j3,res);
		cvAdd(res,resultvec+i3,resultvec+i3);
	}
}

//=====================================================================================================================

static inline void DiagProd(doublecomplex * restrict argvec,doublecomplex * restrict resultvec,const size_t i)
/* Multiplies the result in the i'th block of resultvec by cc_sqrt at that block, subtracts the result from the i'th
 * block of argvec, and stores the result in the i'th block if resultvec.
 */
{
	const size_t i3 = i*3;

	resultvec[i3]   = argvec[i3]   - resultvec[i3]*cc_sqrt[material[i]][0];
	resultvec[i3+1] = argvec[i3+1] - resultvec[i3+1]*cc_sqrt[material[i]][1];
	resultvec[i3+2] = argvec[i3+2] - resultvec[i3+2]*cc_sqrt[material[i]][2];
}

#endif // !USE_SSE3

#endif // __sparse_ops_h

#endif // SPARSE
