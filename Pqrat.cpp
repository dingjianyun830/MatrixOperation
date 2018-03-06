#include "mex.h"
#include "matrix.h"
#include "math.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex Pqrat.cpp
N = 8;
x = [1:N]
y = randi(6,[1 N])
t = randi(6, 1)
Pqrat(x,y,t)
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

	double t = m_nInputData3[0];
	int n = N1;

	if (N1 != N2)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input x y diminsion should be match.");
		// exit MEX file
	}

	int i, j, k, m, l;
	double h, b[8];
	double z = 0;
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
	if (n <= 8)
	{
		k = 0;
		m = 8;
	}
	else if (t < m_nInputData1[4])
	{
		k = 0;
		m = 8;
	}
	else if (t > m_nInputData1[n - 5])
	{
		k = n - 8;
		m = 8;
	}
	else
	{
		k = 1;
		j = n;
		while (j - k != 1)
		{
			i = (k + j) / 2;
			if (t < m_nInputData1[i - 1])
			{
				j = i;
			}
			else
			{
				k = i;
			}
		}
		k = k - 4;
		m = 8;
	}
	b[0] = m_nInputData2[k];
	for (i = 2; i <= m; i++)
	{
		h = m_nInputData2[i + k - 1];
		l = 0;
		j = 1;
		while ((l == 0) && (j <= i - 1))
		{
			if (fabs(h - b[j - 1]) + 1 == 1)
			{
				l = 1;
			}
			else
			{
				h = (m_nInputData1[i + k - 1] - m_nInputData1[j + k - 1]) / (h - b[j - 1]);
			}
			j = j + 1;
		}
		b[i - 1] = h;
		if (l != 0)
		{
			b[i - 1] = 1.0e+35;
		}
	}
	z = b[m - 1];
	for (i = m - 1; i >= 1; i--)
	{
		z = b[i - 1] + (t - m_nInputData1[i + k - 1]) / z;
	}
	m_nOutputData[0] = z;
}