#include "mex.h"
#include "matrix.h"
#include "Math.h"

/*This is a Matlab mexwarp to validate my Matrix Rank*/
/* Matlab code
mex MatrixLU.cpp
A = randi(6,[5 5])
[L,U]=MatrixLU(A)
L*U
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData;// Matrix Input M1*N1
	int M, N;
	int i, j, k;

	m_nInputData = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);

	double *m_nOutputData1;// L
	double *m_nOutputData2;// U
	plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
	m_nOutputData1 = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
	m_nOutputData2 = mxGetPr(plhs[1]);

	//B = MatrixLU(A)
	if (N != M)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input Matrix should be SQUARE.");
		// exit MEX file
	}
	else
	{
		int i, j, k;
		int n = M;
		for (k = 0; k <= n - 2; k++)
		{
			if (fabs(m_nInputData[k*n+k]) + 1 == 1)
			{
				mexErrMsgIdAndTxt("MyProg:ConvertString",
					"The Input Matrix can not be LU deposed.");
			}

			for (i = k + 1; i <= n - 1; i++)
			{
				m_nInputData[i*n+k] = m_nInputData[i*n+k] / m_nInputData[k*n+k];
			}

			for (i = k + 1; i <= n - 1; i++)
			{
				for (j = k + 1; j <= n - 1; j++)
				{
					m_nInputData[i*n+j] = m_nInputData[i*n + j] - m_nInputData[i*n + k] * m_nInputData[k*n + j];
				}
			}
		}

		for (i = 0; i <= n - 1; i++)
		{
			for (j = 0; j < i; j++)
			{
				m_nOutputData2[i*n+j] = m_nInputData[i*n+j];
				m_nOutputData1[i*n+j] = 0;
			}
			m_nOutputData2[i*n+i] = 1;
			m_nOutputData1[i*n+i] = m_nInputData[i*n + i];
			for (j = i + 1; j <= n - 1; j++)
			{
				m_nOutputData2[i*n + j] = 0;
				m_nOutputData1[i*n + j] = m_nInputData[i*n + j];
			}
		}
	}
}