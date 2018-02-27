#include "mex.h"
#include "matrix.h"
#include "Math.h"

/*This is a Matlab mexwarp to validate my Matrix Rank*/
/* Matlab code
mex MatrixQR.cpp
A = randi(6,[5 5])
[Q,R]=MatrixQR(A)
Q*R
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData;// Matrix Input M*N
	int M, N;
	int i, j, k;

	m_nInputData = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);

	double *m_nOutputData1;// Q
	double *m_nOutputData2;// R
	plhs[0] = mxCreateDoubleMatrix(M, M, mxREAL);
	m_nOutputData1 = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
	m_nOutputData2 = mxGetPr(plhs[1]);

	//B = MatrixLU(A)
	if (M<N)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input Matrix should be M>N.");
		// exit MEX file
	}
	else
	{
		int i, j, k, nn, jj;
		double u, alpha, w, t;

		for (i = 0; i <= M - 1; i++)
		{
			for (j = 0; j <= M - 1; j++)
			{
				m_nOutputData1[j*M + i] = 0;
				if (i == j)
				{
					m_nOutputData1[j*M + i] = 1;
				}
			}
		}
		nn = N;
		if (M == N)
		{
			nn = M - 1;
		}

		for (k = 0; k <= nn - 1; k++)
		{
			u = 0;
			for (i = k; i <= M - 1; i++)
			{
				w = fabs(m_nInputData[k*M + i]);
				if (w > u)
				{
					u = w;
				}
			}
			alpha = 0;
			for (i = k; i <= M - 1; i++)
			{
				t = m_nInputData[k*M + i] / u;
				alpha = alpha + t*t;
			}
			if (m_nInputData[k*M + k] > 0)
			{
				u = -u;
			}
			alpha = u*sqrt(alpha);
			if (fabs(alpha) + 1 == 1)
			{
				mexErrMsgIdAndTxt("MyProg:ConvertString",
					"The Input Matrix can not be QR.");
			}
			u = sqrt(2 * alpha*(alpha - m_nInputData[k*M + k]));
			if ((u + 1) != 1)
			{
				m_nInputData[k*M + k] = (m_nInputData[k*M + k] - alpha) / u;
				for (i = k + 1; i <= M - 1; i++)
				{
					m_nInputData[k*M + i] = m_nInputData[k*M + i] / u;
				}

				for (j = 0; j <= M - 1; j++)
				{
					t = 0;
					for (jj = k; jj <= M - 1; jj++)
					{
						t = t + m_nInputData[k*M + jj] * m_nOutputData1[j*M + jj];
					}
					for (i = k; i <= M - 1; i++)
					{
						m_nOutputData1[j*M + i] = m_nOutputData1[j*M + i] - 2 * t * m_nInputData[k*M + i];
					}
				}

				for (j = k + 1; j <= N - 1; j++)
				{
					t = 0;
					for (jj = k; jj <= M - 1; jj++)
					{
						t = t + m_nInputData[k*M + jj] * m_nInputData[j*M + jj];
					}
					for (i = k; i <= M - 1; i++)
					{
						m_nInputData[j*M + i] = m_nInputData[j*M + i] - 2 * t*m_nInputData[k*M + i];
					}
				}

				m_nInputData[k*M + k] = alpha;

				for (i = k + 1; i <= M - 1; i++)
				{
					m_nInputData[k*M + i] = 0;
				}
			}
		}

		for (i = 0; i <= M - 2; i++)
		{
			for (j = i + 1; j <= M - 1; j++)
			{
				t = m_nOutputData1[j*M + i];
				m_nOutputData1[j*M + i] = m_nOutputData1[i*M + j];
				m_nOutputData1[i*M + j] = t;
			}
		}

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				m_nOutputData2[j*M + i] = m_nInputData[j*M + i];
			}
		}
	}
}