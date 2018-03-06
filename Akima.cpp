#include "mex.h"
#include "matrix.h"
#include "math.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex Akima.cpp
N = 8;
x = [1:N]
y = randi(6,[1 N])
t = 1.3;
Akima(x,y,u,v)
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
			"The Input m_nInputData1 m_nInputData2 diminsion should be match.");
		// exit MEX file
	}

	int m, kk;
	int k;
	double u[5];
	double p, q;
	double z = 0;
	double s[4];
	if (n < 1)
	{
		k = 0;
		m_nOutputData[0] = z;
	}
	if (n == 1)
	{
		k = 0;
		s[0] = m_nInputData2[0];
		z = m_nInputData2[0];
		m_nOutputData[0] = z;
	}
	if (n == 2)
	{
		k = 0;
		s[0] = m_nInputData2[0];
		s[1] = (m_nInputData2[1] - m_nInputData2[0]) / (m_nInputData1[1] - m_nInputData1[0]);
		z = (m_nInputData2[0] * (t - m_nInputData1[1]) - m_nInputData2[1] * (t - m_nInputData1[0])) / (m_nInputData1[0] - m_nInputData1[1]);
		m_nOutputData[0] = z;
	}
	if (t <= m_nInputData1[1])
	{
		k = 0;
	}
	else if (t >= m_nInputData1[n - 1])
	{
		k = n - 2;
	}
	else
	{
		k = 1;
		m = n;
		while (((k - m) != 1) && ((k - m) != -1))
		{
			kk = (k + m) / 2;
			if (t < m_nInputData1[kk - 1])
			{
				m = kk;
			}
			else
			{
				k = kk;
			}
		}
		k = k - 1;
	}
	u[2] = (m_nInputData2[k + 1] - m_nInputData2[k]) / (m_nInputData1[k + 1] - m_nInputData1[k]);

	if (n == 3)
	{
		if (k == 0)
		{
			u[3] = (m_nInputData2[2] - m_nInputData2[1]) / (m_nInputData1[2] - m_nInputData1[1]);
			u[4] = 2 * u[3] - u[2];
			u[1] = 2 * u[2] - u[3];
			u[0] = 2 * u[1] - u[2];
		}
		else
		{
			u[1] = (m_nInputData2[1] - m_nInputData2[0]) / (m_nInputData1[1] - m_nInputData1[0]);
			u[0] = 2 * u[1] - u[2];
			u[3] = 2 * u[2] - u[1];
			u[4] = 2 * u[3] - u[2];
		}
	}
	else
	{
		if (k <= 1)
		{
			u[3] = (m_nInputData2[k + 2] - m_nInputData2[k + 1]) / (m_nInputData1[k + 2] - m_nInputData1[k + 1]);
			if (k == 1)
			{
				u[1] = (m_nInputData2[1] - m_nInputData2[0]) / (m_nInputData1[1] - m_nInputData1[0]);
				u[0] = 2 * u[1] - u[2];
				if (n == 4)
				{
					u[4] = 2 * u[3] - u[2];
				}
				else
				{
					u[4] = (m_nInputData2[4] - m_nInputData2[3]) / (m_nInputData1[4] - m_nInputData1[3]);
				}
			}
			else if (k >= (n - 3))
			{
				u[1] = (m_nInputData2[k] - m_nInputData2[k - 1]) / (m_nInputData1[k] - m_nInputData1[k - 1]);
				if (k == (n - 3))
				{
					u[3] = (m_nInputData2[n - 1] - m_nInputData2[n - 2]) / (m_nInputData1[n - 1] - m_nInputData1[n - 2]);
					u[4] = 2 * u[3] - u[2];
					if (n == 4)
					{
						u[0] = 2 * u[1] - u[2];
					}
					else
					{
						u[0] = (m_nInputData2[k - 1] - m_nInputData2[k - 2]) / (m_nInputData1[k - 1] - m_nInputData1[k - 2]);
					}
				}
				else
				{
					u[3] = 2 * u[2] - u[1];
					u[4] = 2 * u[3] - u[2];
					u[0] = (m_nInputData2[k - 1] - m_nInputData2[k - 2]) / (m_nInputData1[k - 1] - m_nInputData1[k - 2]);
				}
			}
			else
			{
				u[1] = (m_nInputData2[k] - m_nInputData2[k - 1]) / (m_nInputData1[k] - m_nInputData1[k - 1]);
				u[0] = (m_nInputData2[k - 1] - m_nInputData2[k - 2]) / (m_nInputData1[k - 1] - m_nInputData1[k - 2]);
				u[3] = (m_nInputData2[k + 2] - m_nInputData2[k + 1]) / (m_nInputData1[k + 2] - m_nInputData1[k + 1]);
				u[4] = (m_nInputData2[k + 3] - m_nInputData2[k + 2]) / (m_nInputData1[k + 3] - m_nInputData1[k + 2]);
			}
		}
	}
	s[0] = fabs(u[3] - u[2]);
	s[1] = fabs(u[0] - u[1]);
	if ((s[0] + 1 == 1) && (s[1] + 1 == 1))
	{
		p = (u[1] + u[2]) / 2;
	}
	else
	{
		p = (s[0] * u[1] + s[1] * u[2]) / (s[0] + s[1]);
	}
	s[0] = fabs(u[3] - u[4]);
	s[1] = fabs(u[0] - u[1]);
	if ((s[0] + 1 == 1) && (s[1] + 1 == 1))
	{
		q = (u[2] + u[3]) / 2;
	}
	else
	{
		q = (s[0] * u[2] + s[1] * u[3]) / (s[0] + s[1]);
	}
	s[0] = m_nInputData2[k];
	s[1] = p;
	s[3] = m_nInputData1[k + 1] - m_nInputData1[k];
	s[2] = (3 * u[2] - 2 * p - q) / s[3];
	s[3] = (q + p - 2 * u[2]) / (s[3] * s[3]);
	p = t - m_nInputData1[k];
	z = s[0] + s[1] * p + s[2] * p*p + s[3] * p*p*p;

	m_nOutputData[0] = z;
}