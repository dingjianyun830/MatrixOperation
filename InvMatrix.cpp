#include "mex.h"
#include "matrix.h"
#include "Math.h"

/*This is a Matlab mexwarp to validate my Matrix Inverse*/
/*Matlab code
mex InvMatrix.cpp
A = randi(6,[5 5])
InvMatrix(A)
inv(A)
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
	plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	//B = InvMatrix(A)
	if (N != M)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input Matrix be SQUARE.");
		// exit MEX file
	}
	else
	{
		int n = M;
		int *is, *js;	
		double d, p;
		is = new int[n];
		js = new int[n];
		for (k = 0; k <= n - 1; k++)
		{
			d = 0;
			for (i = k; i <= n - 1; i++)
			{
				for (j = k; j <= n - 1; j++)
				{
					p = fabs(m_nInputData[j*n+i]);
					if (p > d)
					{
						d = p;
						is[k] = i;
						js[k] = j;
					}
				}
			}
			if (d + 1 == 1)
			{
				delete[] is, js;
				mexErrMsgIdAndTxt("MyProg:ConvertString",
					"The Input Matrix can not Inverse.");
			}

			if (is[k] != k)
			{
				for (j = 0; j <= n - 1; j++)
				{
					p = m_nInputData[k*n + j];
					m_nInputData[k*n + j] = m_nInputData[is[k] * n + j];
					m_nInputData[is[k] * n + j] = p;
				}
			}

			if (js[k] != k)
			{
				for (i = 0; i <= n - 1; i++)
				{
					p = m_nInputData[i*n + k];
					m_nInputData[i*n + k] = m_nInputData[i*n + js[k]];
					m_nInputData[i*n + js[k]] = p;
				}
			}

			m_nInputData[k*n + k] = 1 / m_nInputData[k*n + k];

			for (j = 0; j <= n - 1; j++)
			{
				if (j != k)
				{
					m_nInputData[k*n + j] = m_nInputData[k*n + j] * m_nInputData[k*n + k];
				}
			}

			for (i = 0; i <= n - 1; i++)
			{
				if (i != k)
				{
					for (j = 0; j <= n - 1; j++)
					{
						if (j != k)
						{
							m_nInputData[i*n + j] = m_nInputData[i*n + j] - m_nInputData[i*n + k] * m_nInputData[k*n + j];
						}
					}
				}
			}

			for (i = 0; i <= n - 1; i++)
			{
				if (i != k)
				{
					m_nInputData[i*n + k] = -m_nInputData[i*n + k] * m_nInputData[k*n + k];
				}
			}
		}

		for (k = n - 1; k >= 0; k--)
		{
			if (js[k] != k)
			{
				for (j = 0; j <= n - 1; j++)
				{
					p = m_nInputData[k*n + j];
					m_nInputData[k*n + j] = m_nInputData[js[k] * n + j];
					m_nInputData[js[k] * n + j] = p;
				}
			}

			if (is[k] != k)
			{
				for (i = 0; i <= n - 1; i++)
				{
					p = m_nInputData[i*n + k];
					m_nInputData[i*n + k] = m_nInputData[i*n + is[k]];
					m_nInputData[i*n + is[k]] = p;
				}
			}
		}
		delete[] is, js;
	}

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			m_nOutputData[j*M + i] = m_nInputData[j*M + i];
		}
	}
}