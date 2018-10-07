/* In order to validate my matrix code, I use MATLAB to offer input&output */
#include "CMatrixOperation.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *inData;
	double *outData;
	int M, N;
	int i, j;
	inData = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
	outData = mxGetPr(plhs[0]);

	Mat inSrc(inData, M, N);
	Mat outSrc(M, N);
	outSrc = MatDotproduct(inSrc, 2.5);

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			outData[j*M + i] = outSrc.data[i][j];
		}
	}
}