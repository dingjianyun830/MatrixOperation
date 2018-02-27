#include "mex.h"
#include "matrix.h"
#include "Math.h"

/*This is a Matlab mexwarp to validate my Matrix Rank*/
/* Matlab code
mex RankMatrix.cpp
A = randi(6,[5 5])
RankMatrix(A)
rank(A)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData;// Matrix Input M1*N1
	double *m_nOutputData;
	int M, N;
	int i, j, k;

	m_nInputData = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	//B = RankMatrix(A)
	if (N ==1&& M==1)
	{
		m_nOutputData[0] = 1;
		// exit MEX file
	}
	else
	{
		int i, j, is, js, k, rank;
		double q, d;
		rank = M;
		if (M >= N)
		{
			rank = N;
		}
		for (k = 0; k <= rank - 1; k++)
		{
			q = 0;
			for (i = k; i <= M - 1; i++)
			{
				for (j = k; j <= N - 1; j++)
				{
					d = fabs(m_nInputData[i*M+j]);
					if (d > q)
					{
						q = d;
						is = i;
						js = j;
					}
				}
			}

			if (q + 1 == 1)
			{
				mexErrMsgIdAndTxt("MyProg:ConvertString",
					"The Input Matrix dont have rank.");
			}

			m_nOutputData[0] = m_nOutputData[0] + 1;

			if (is != k)
			{
				for (j = k; j <= N - 1; j++)
				{
					d = m_nInputData[k*M + j];
					m_nInputData[k*M+j] = m_nInputData[is*M+j];
					m_nInputData[is*M+j] = d;
				}
			}

			if (js != k)
			{
				for (i = k; i <= M - 1; i++)
				{
					d = m_nInputData[i*M + js];
					m_nInputData[i*M + js] = m_nInputData[i*M + k];
					m_nInputData[i*M + k] = d;
				}
			}

			for (i = k + 1; i <= M - 1; i++)
			{
				d = m_nInputData[i*M+k] / m_nInputData[k*M+k];
				for (j = k + 1; j <= N - 1; j++)
				{
					m_nInputData[i*M + j] = m_nInputData[i*M + j] - d * m_nInputData[k*M + j];
				}
			}
		}
	}
}