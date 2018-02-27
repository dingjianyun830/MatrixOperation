#include "mex.h"
#include "matrix.h"
#include "Math.h"

/*This is a Matlab mexwarp to validate my Matrix Det*/
/* Matlab code
mex DetMatrix.cpp
A = randi(6,[5 5])
DetMatrix(A)
det(A)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData;// Matrix Input M1*N1
	double *m_nOutputData;
	int M, N;

	m_nInputData = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);
	
	//c = det(A)
	if (N != M)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input Matrix must be SQUARE.");
		// exit MEX file
	}
	else
	{
		int i, j, k;
		int is, js;
		double  f, q, d;
		f = 1;
		double det = 1;
		int n = M;
		for (k = 0; k <= n - 2; k++)
		{
			q = 0;
			for (i = k; i <= n - 1; i++)
			{
				for (j = k; k <= n - 1; j++)
				{
					d = fabs(m_nInputData[j*n + i]);
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
				det = 0;
				m_nOutputData[0] = det;
			}
			else
			{
				if (is != k)
				{
					f = -f;
					for (j = k; j <= n - 1; j++)
					{
						d = m_nInputData[k*n + j];
						m_nInputData[k*n + j] = m_nInputData[is*n + j];
						m_nInputData[is*n + j] = d;
					}
				}

				if (js != k)
				{
					f = -f;
					for (i = k; i <= n - 1; i++)
					{
						d = m_nInputData[i*n + js];
						m_nInputData[i*n + js] = m_nInputData[i*n + k];
						m_nInputData[i*n + k] = d;
					}
				}

				det = det*m_nInputData[k*n + k];

				for (i = k + 1; i <= n - 1; i++)
				{
					d = m_nInputData[i*n + k] / m_nInputData[k*n + k];
					for (j = k + 1; j <= n - 1; j++)
					{
						m_nInputData[i*n + j] = m_nInputData[i*n + j] - d * m_nInputData[k*n + j];
					}
				}
			}
		}
		det = f*det*m_nInputData[(n-1)*n + n-1];
		m_nOutputData[0] = det;
	}	
}