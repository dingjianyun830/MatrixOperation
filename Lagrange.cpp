#include "mex.h"
#include "matrix.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex Lagrange.cpp
N = 20;
x = [1:N]
y = randi(6,[1 N])
t = randi(6, 1)
Lagrange(x,y,t)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData1;// Matrix Input M1*N1
	double *m_nInputData2;// Matrix Input M2*N2
	double *m_nInputData3;
	double *m_nOutputData;
	int M1, N1, M2, N2;

	m_nInputData1 = mxGetPr(prhs[0]);
	m_nInputData2 = mxGetPr(prhs[1]);
	m_nInputData3 = mxGetPr(prhs[2]);
	M1 = mxGetM(prhs[0]);
	N1 = mxGetN(prhs[0]);
	M2 = mxGetM(prhs[1]);
	N2 = mxGetN(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	int i, j, k, m, n;
	double s;
	double z = 0;
	double t = m_nInputData3[0];

	n = N1;
	if (N1 != N2)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input x y diminsion should be match.");
		// exit MEX file
	}

	if (n < 1)
	{
		z = 0;
		m_nOutputData[0] = z;
	}
	if (n == 1)
	{
		z = m_nInputData2[0];
		m_nOutputData[0] = z;
	}
	if (n == 2)
	{
		z = (m_nInputData2[0] * (t - m_nInputData1[1]) - m_nInputData2[1] * (t - m_nInputData1[0])) / (m_nInputData1[0] - m_nInputData1[1]);
		m_nOutputData[0] = z;
	}

	i = 0;
	while ((m_nInputData1[i] < t) && (i < n))
	{
		i = i + 1;
	}
	k = i - 4;
	if (k < 0)
	{
		k = 0;
	}
	m = i + 3;
	if (m > n - 1)
	{
		m = n - 1;
	}
	for (i = k; i <= m; i++)
	{
		s = 1;
		for (j = k; j <= m; j++)
		{
			if (j != i)
			{
				s = s*(t - m_nInputData1[j]) / (m_nInputData1[i] - m_nInputData1[j]);
			}
		}
		z = z + s*m_nInputData2[i];
	}
	m_nOutputData[0] = z;
}