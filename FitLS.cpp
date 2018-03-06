#include "mex.h"
#include "matrix.h"
#include "math.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex FitLS.cpp
N=8;
x = [1:N]
y = randi(6,[1,N])
m = 2;
[a,dt]=FitLS(x,y,m)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData1;// Matrix Input M1*N1---x
	double *m_nInputData2;// Matrix Input M2*N2---y
	double *m_nInputData3;// m
	m_nInputData1 = mxGetPr(prhs[0]);
	m_nInputData2 = mxGetPr(prhs[1]);
	m_nInputData3 = mxGetPr(prhs[2]);
	int M1, N1, M2, N2;
	M1 = mxGetM(prhs[0]);
	N1 = mxGetN(prhs[0]);
	M2 = mxGetM(prhs[1]);
	N2 = mxGetN(prhs[1]);
	int m = m_nInputData3[0];

	double *m_nOutputData1;//Matrix Output ---a
	double *m_nOutputData2;//Matrix Output ---dt
	plhs[0] = mxCreateDoubleMatrix(1, m+1, mxREAL);
	m_nOutputData1 = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL);
	m_nOutputData2 = mxGetPr(plhs[1]);


	int n = N1;

	int i, j, k;
	double z, p, c, g, q, d1, d2;
	double s[20], t[20], b[20];

	for (i = 0; i <= m; i++)
	{
		m_nOutputData1[i] = 0;
	}

	if (m + 1 > n)
	{
		m = n - 1;
	}

	if (m > 19)
	{
		m = 19;
	}

	z = 0.0;

	for (i = 0; i <= n - 1; i++)
	{
		z = z + m_nInputData1[i] / (1.0 * n);
	}

	b[0] = 1.0;
	d1 = 1.0*n;
	p = 0;
	c = 0;

	for (i = 0; i <= n - 1; i++)
	{
		p = p + (m_nInputData1[i] - z);
		c = c + m_nInputData2[i];
	}
	c = c / d1;
	p = p / d1;
	m_nOutputData1[0] = c*b[0];
	if (m > 0)
	{
		t[1] = 1;
		t[0] = -p;
		d2 = 0;
		c = 0;
		g = 0;
		for (i = 0; i <= n - 1; i++)
		{
			q = m_nInputData1[i] - z - p;
			d2 = d2 + q*q;
			c = c + m_nInputData2[i] * q;
			g = g + (m_nInputData1[i] - z)*q*q;
		}
		c = c / d2;
		p = g / d2;
		q = d2 / d1;
		d1 = d2;
		m_nOutputData1[1] = c*t[1];
		m_nOutputData1[0] = c*t[0] + m_nOutputData1[0];
	}

	for (j = 2; j <= m; j++)
	{
		s[j] = t[j - 1];
		s[j - 1] = -p*t[j - 1] + t[j - 2];
		if (j >= 3)
		{
			for (k = j - 2; k >= 1; k--)
			{
				s[k] = -p*t[k] + t[k - 1] - q*b[k];
			}
		}
		s[0] = -p*t[0] - q*b[0];
		d2 = 0;
		c = 0;
		g = 0;
		for (i = 0; i <= n - 1; i++)
		{
			q = s[j];
			for (k = j - 1; k >= 0; k--)
			{
				q = q*(m_nInputData1[i] - z) + s[k];
			}
			d2 = d2 + q*q;
			c = c + m_nInputData2[i] * q;
			g = g + (m_nInputData1[i] - z)*q*q;
		}

		c = c / d2;
		p = g / d2;
		q = d2 / d1;
		d1 = d2;
		m_nOutputData1[j] = c*s[j];
		t[j] = s[j];
		for (k = j - 1; k >= 0; k--)
		{
			m_nOutputData1[k] = c*s[k] + m_nOutputData1[k];
			b[k] = t[k];
			t[k] = s[k];
		}
	}

	m_nOutputData2[0] = 0;
	m_nOutputData2[1] = 0;
	m_nOutputData2[2] = 0;
	for (i = 0; i <= n - 1; i++)
	{
		q = m_nOutputData1[m];
		for (k = m - 1; k >= 0; k--)
		{
			q = m_nOutputData1[k] + q*(m_nInputData1[i] - z);
		}
		p = q - m_nInputData2[i];
		if (fabs(p) > m_nOutputData2[2])
		{
			m_nOutputData2[2] = fabs(p);
		}
		m_nOutputData2[0] = m_nOutputData2[0] + p*p;
		m_nOutputData2[1] = m_nOutputData2[1] + fabs(p);
	}
}