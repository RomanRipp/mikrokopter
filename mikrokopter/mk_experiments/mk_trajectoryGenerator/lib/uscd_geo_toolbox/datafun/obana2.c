/* 
 *  OBANA2.C (obana2.mex*)
 *
 *  usage:
 *          znew = obana2(z,x,y,gridx,gridy,xrad,xcut,[yrad,ycut])
 *
 *  adapted from the original Matlab m-file obana.m by Martin Visbeck
 *  and additional extensions by Gerd Krahmann
 *  
 *  Matlab 4.2c mex-file 
 *  by Ulf Garternicht at IfM Kiel
 *  $Id: obana2.c,v 1.6 1997/08/25 09:22:40 garter Exp garter $
 */

#include <math.h>
#include "mex.h"

#define	ZOLD_IN		prhs[0]
#define	XOLD_IN		prhs[1]
#define	YOLD_IN		prhs[2]
#define	XNEW_IN		prhs[3]
#define	YNEW_IN		prhs[4]
#define	XRAD_IN		prhs[5]
#define	XCUT_IN		prhs[6]
#define	YRAD_IN		prhs[7]
#define	YCUT_IN		prhs[8]

#define	ZNEW_OUT	plhs[0]

/* ------------------------------------------------------------------------ * 
 *   OBANA2C (for complex ZOLD_IN)                                          *
 * ------------------------------------------------------------------------ */

double obana2c(
	double		znew_r[],
	double		znew_i[],
	double		zold_r[],
	double		zold_i[],
	double		xold[],
	double		yold[],
	unsigned int	mold,
	unsigned int	nold,
	double		xnew[],
	double		ynew[],
	unsigned int	mnew,
	unsigned int	nnew,
	double		xrad[],
	double		xcut[],
	double		yrad[],
	double		ycut[]
	)
{
	int		i1, i2, found_one, NaN_warning = 0;
	double		q2, q, qsum, qq, qdeg;
	const double	pifac = atan(1.0)/45.0; /* PI/180 */

	mexPrintf("Obana2 started, patience please ...\n");

	for (i1 = 0; i1 < (mnew*nnew); i1++) {

		/* print progress message */
		if ((i1%25) == 0)
			mexPrintf("%3.0f%c completed\r",(double)i1/(double)(mnew*nnew)*100.0,'%');

		/* reset values */
		found_one = 0;
		qsum = 0;
		znew_r[i1] = 0;
		znew_i[i1] = 0;

		for (i2 = 0; i2 < (mold*nold); i2++) {

			/* either for geographical coordinates (x=lon,y=lat) */
			if (yrad[0] == -1.0) {

				q2 = sin(yold[i2]*pifac)*sin(ynew[i1]*pifac)
					+cos(yold[i2]*pifac)*cos(ynew[i1]*pifac)
					*cos((xold[i2]-xnew[i1])*pifac); 
				qq = abs(q2);
				if (qq > 1.0) q2 = 1.0;
				q2 = atan(sqrt((1.0-q2)/(1.0+q2)))*222.24/pifac;					
				if (mexIsNaN(q2)) { /* zero distance ? */
					qdeg = 0.0;
					q2 = 0.0;
				} else {
					qdeg = q2 / 111.12;
					q2 = qdeg / xrad[i1];
				}

				if (qdeg < xcut[i1] ) {

					found_one = 1;

					q = exp(-q2*q2);

					/* sum it up */
					qsum += q;		
					znew_r[i1] += q*zold_r[i2];
					znew_i[i1] += q*zold_i[i2];
				}

			/* or for cartesian coordinates */
			} else {

				/* use values inside cut off radius only */
				if (sqrt((xold[i2]-xnew[i1])*(xold[i2]-xnew[i1])/(xcut[i1]*xcut[i1])
					+(yold[i2]-ynew[i1])*(yold[i2]-ynew[i1])/(ycut[i1]*ycut[i1])) < 1) {

					found_one = 1;

					/* compute Gaussian weighting factor */
					q = sqrt((xold[i2]-xnew[i1])*(xold[i2]-xnew[i1])/(xrad[i1]*xrad[i1])
						+(yold[i2]-ynew[i1])*(yold[i2]-ynew[i1])/(yrad[i1]*yrad[i1])); 
					q = exp(-q*q);

					/* sum it up */
					qsum += q;		
					znew_r[i1] += q*zold_r[i2];
					znew_i[i1] += q*zold_i[i2];
				}
			}
		}

		
		if (qsum > 0) {
			znew_r[i1] /= qsum;
			znew_i[i1] /= qsum;

		} else if (qsum == 0 && found_one) {
			/* in case of too small influence radius */
			znew_r[i1] = mexGetNaN();
			znew_i[i1] = mexGetNaN();
			NaN_warning = 1;   

		} else {
                        /* if no values found set to NaN */
			znew_r[i1] = mexGetNaN();
			znew_i[i1] = mexGetNaN();
		}
	}

	mexPrintf("100%c completed\n",'%');
	if (NaN_warning)
		mexErrMsgTxt("Influence radius is too small.\n");
	return;
}

/* ------------------------------------------------------------------------ * 
 *   OBANA2R (for real ZOLD_IN)                                             *
 * ------------------------------------------------------------------------ */

double obana2r(
	double		znew_r[],
	double		zold_r[],
	double		xold[],
	double		yold[],
	unsigned int	mold,
	unsigned int	nold,
	double		xnew[],
	double		ynew[],
	unsigned int	mnew,
	unsigned int	nnew,
	double		xrad[],
	double		xcut[],
	double		yrad[],
	double		ycut[]
	)
{
	int		i1, i2, found_one, NaN_warning = 0;
	double		q2, qdeg, qq, q, qsum;
	const double	pifac = atan(1.0)/45.0; /* PI/180 */

	mexPrintf("Obana2 started, patience please ...\n");

	for (i1 = 0; i1 < (mnew*nnew); i1++) {

		/* print progress message */
		if ((i1%25) == 0)
			mexPrintf("%3.0f%c completed\r",(double)i1/(double)(mnew*nnew)*100.0,'%');

		/* reset values */
		qsum = 0;
		znew_r[i1] = 0;
		found_one = 0;

		for (i2 = 0; i2 < (mold*nold); i2++) {

			/* either for geographical coordinates (x=lon,y=lat) */
			if (yrad[0] == -1.0) {
				q2 = sin(yold[i2]*pifac)*sin(ynew[i1]*pifac)
					+cos(yold[i2]*pifac)*cos(ynew[i1]*pifac)
					*cos((xold[i2]-xnew[i1])*pifac); 
				qq = abs(q2);
				if (qq > 1.0) q2 = 1.0;
				q2 = atan(sqrt((1.-q2)/(1.+q2)))*222.24/pifac;					
				if (mexIsNaN(q2)) { /* zero distance ? */
					qdeg = 0.0;
					q2 = 0.0;
				} else {
					qdeg = q2 / 111.12;
					q2 = qdeg / xrad[i1];
				}

				if (qdeg < xcut[i1] ) {

					found_one = 1;

					q = exp(-q2*q2);

					/* sum it up */
					qsum += q;		
					znew_r[i1] += q*zold_r[i2];
				}

			/* or for cartesian coordinates */
			} else {
				/* use values inside cut off radius only */
				if (sqrt((xold[i2]-xnew[i1])*(xold[i2]-xnew[i1])/(xcut[i1]*xcut[i1])
					+(yold[i2]-ynew[i1])*(yold[i2]-ynew[i1])/(ycut[i1]*ycut[i1])) < 1) {

					found_one = 1;

					/* compute Gaussian weighting factor */
					q = sqrt((xold[i2]-xnew[i1])*(xold[i2]-xnew[i1])/(xrad[i1]*xrad[i1])
						+(yold[i2]-ynew[i1])*(yold[i2]-ynew[i1])/(yrad[i1]*yrad[i1])); 
					q = exp(-q*q);

					/* sum it up */
					qsum += q;		
					znew_r[i1] += q*zold_r[i2];
				}
			}
		}

		
		if (qsum > 0) {
			znew_r[i1] /= qsum;
		} else if (qsum == 0 && found_one) {
			/* in case of too small influence radius */
			NaN_warning = 1;   
			znew_r[i1] = mexGetNaN();
		} else {
                        /* if no values found set to NaN */
			znew_r[i1] = mexGetNaN();
		}
	
	}

	mexPrintf("100%c completed\n",'%');
	if (NaN_warning)
		mexErrMsgTxt("Influence radius is too small.\n");
	return;
}

/* ------------------------------------------------------------------------ *
 *   MAIN                                                                   *
 * ------------------------------------------------------------------------ */

/* void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[]) */

#ifdef __STDC__ /* Linux */
void mexFunction(
	int	nlhs,
	Matrix	*plhs[],
	int	nrhs,
	Matrix	*prhs[]
	)
#else /* Alpha */
mexFunction(nlhs, plhs, nrhs, prhs)
int	nlhs, nrhs;
Matrix	*plhs[], *prhs[];
#endif

{
	double		*znew_r, *znew_i,
			*zold_r, *zold_i, *xold, *yold, *xnew, *ynew,
			*xrad, *xcut, *yrad, *ycut;
	unsigned int	mold, nold, mnew, nnew;

	/* Check for proper number of arguments */

	if (nrhs != 7 && nrhs != 9) {
		mexErrMsgTxt("USAGE: znew = obana2(z,x,y,gridx,gridy,xrad,xcut,yrad,ycut)");
	} else if (nlhs > 1) {
		mexErrMsgTxt("One output argument only.");
	}

	/* Check if the dimensions of ZOLD_IN, XOLD_IN, YOLD_IN, are the same */

	mold = mxGetM(ZOLD_IN);
	nold = mxGetN(ZOLD_IN);

	if (   mxGetM(XOLD_IN) != mold || mxGetM(YOLD_IN) != mold 
	    || mxGetN(XOLD_IN) != nold || mxGetN(YOLD_IN) != nold) {
		mexErrMsgTxt("Input parameters x,y,z have different size.");
	}

	/* Check if the dimensions of XNEW_IN, YNEW_IN, XRAD_IN, XCUT_IN, and
	optionally YRAD_IN, YCUT_IN are the same */

	mnew = mxGetM(XNEW_IN);
	nnew = mxGetN(XNEW_IN);

	if (   mxGetM(YNEW_IN) != mnew || mxGetN(YNEW_IN) != nnew
	    || mxGetM(XRAD_IN) != mnew || mxGetN(XRAD_IN) != nnew
	    || mxGetM(XCUT_IN) != mnew || mxGetN(XCUT_IN) != nnew) {
		mexErrMsgTxt("gridx,gridy,xrad,xcut must all have the same size.");
	}
	if (nrhs == 9 && (   mxGetM(YRAD_IN) != mnew || mxGetN(YRAD_IN) != nnew
	                  || mxGetM(YCUT_IN) != mnew || mxGetN(YCUT_IN) != nnew)) {
		mexErrMsgTxt("gridx,gridy,yrad,ycut must all have the same size.");
	}

	/* Assign pointers to the various input parameters */

	zold_r = mxGetPr(ZOLD_IN);
	zold_i = mxGetPi(ZOLD_IN);
	xold   = mxGetPr(XOLD_IN);
	yold   = mxGetPr(YOLD_IN);
	xnew   = mxGetPr(XNEW_IN);
	ynew   = mxGetPr(YNEW_IN);
	xrad   = mxGetPr(XRAD_IN);
	xcut   = mxGetPr(XCUT_IN);

	/* check if isotropic influence and cut-off radius is chosen */

	if (nrhs == 7) {
		yrad = xrad;
		ycut = xcut;
	} else {
		yrad = mxGetPr(YRAD_IN);
		ycut = mxGetPr(YCUT_IN);
	}

	/* 
	 *  1. Create a matrix for the return argument 
         *  2. Assign pointers to the output parameter   
	 *  3. Do the actual computations in a subroutine
         */

	if (mxIsComplex(ZOLD_IN)) {
		ZNEW_OUT = mxCreateFull(mnew, nnew, COMPLEX);

		znew_r = mxGetPr(ZNEW_OUT);
		znew_i = mxGetPi(ZNEW_OUT);

		obana2c(znew_r, znew_i, zold_r, zold_i, xold, yold, mold, nold, 
			xnew, ynew, mnew, nnew, xrad, xcut, yrad, ycut);

	} else {
		ZNEW_OUT = mxCreateFull(mnew, nnew, REAL);

		znew_r = mxGetPr(ZNEW_OUT);

		obana2r(znew_r, zold_r, xold, yold, mold, nold, 
			xnew, ynew, mnew, nnew, xrad, xcut, yrad, ycut);
	}
	return;
}
