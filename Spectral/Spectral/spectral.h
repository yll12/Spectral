#ifndef SPECTRAL_H_
#define SPECTRAL_H_

#include "derivativeMatrix.h"

namespace Spectral
{
	//!  Spectral Methods class. 
	/*!
	Provides various numerical routines.
	*/
	class SpectralMethods
	{
	public:
		//! A static method that computes the differentiation
		//!matrices D1, D2, ..., DM on Chebyshev nodes. \n
		//!The code implements two strategies for enhanced
		//!accuracy suggested by W.Don and S.Solomonoff in
		//!SIAM J.Sci.Comp.Vol. 6, pp. 1253--1268 (1994).
		//!The two strategies are(a) the use of trigonometric
		//!identities to avoid the computation of differences
		//!x(k) - x(j) and(b) the use of the "flipping trick"
		//!which is necessary since sin t can be computed to high
		//!relative precision when t is small whereas sin(pi - t) cannot. \n
		//!
		//!J.A.C.Weideman, S.C.Reddy 1998.
		/*!
		\param n Size of differentiation matrix.        
		\param m Number of derivatives required (integer). \n
             Note: 0 < m <= n-1.
		\return dm[ell] contains ell-th derivative matrix, ell=1..m.
		*/
		static DerivativeMatrix chebdif(int n, int m);
	};

}

#endif 