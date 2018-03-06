#include "mex.h"
#include "matrix.h"
#include "math.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex Lagrange2.cpp
N=8;
M=10;
x = [1:N]
y = [1:M]
z = randi(6,[N,M])
u = 1.3;
v = 2.3;
Lagrange2(x,y,z,u,v)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData1;// Matrix Input M1*N1
	double *m_nInputData2;// Matrix Input M2*N2
	double *m_nInputData3;// Matrix Input M3*N3
	double *m_nInputData4;
	double *m_nInputData5;
	double *m_nOutputData;
	int M1, N1, M2, N2, M3, N3;

	m_nInputData1 = mxGetPr(prhs[0]);
	m_nInputData2 = mxGetPr(prhs[1]);
	m_nInputData3 = mxGetPr(prhs[2]);
	M1 = mxGetM(prhs[0]);
	N1 = mxGetN(prhs[0]);
	M2 = mxGetM(prhs[1]);
	N2 = mxGetN(prhs[1]);
	M3 = mxGetM(prhs[2]);
	N3 = mxGetN(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	m_nInputData4 = mxGetPr(prhs[3]);
	m_nInputData5 = mxGetPr(prhs[3]);
	double u = m_nInputData4[0];
	double v = m_nInputData5[0];
	int n = M3;
	int m = N3;

	int i, j, k, kk;
	int ip, ipp, iq, iqq;
	double h, b[10];
	double f;

	if (u <= m_nInputData1[0])
	{
		ip = 1;
		ipp = 4;
	}
	else if (u >= m_nInputData1[n - 1])
	{
		ip = n - 3;
		ipp = n;
	}
	else
	{
		i = 1;
		j = n;
		while (((i - j) != 1) && ((i - j) != -1))
		{
			kk = (i + j) / 2;
			if (u < m_nInputData1[kk - 1])
			{
				j = kk;
			}
			else
			{
				i = kk;
			}
		}
		ip = i - 3;
		ipp = i + 4;
	}

	if (ip < 1)
	{
		ip = 1;
	}

	if (ipp > n)
	{
		ipp = n;
	}

	if (v <= m_nInputData2[0])
	{
		iq = 1;
		iqq = 4;
	}
	else if (v >= m_nInputData2[m - 1])
	{
		iq = m - 3;
		iqq = m;
	}
	else
	{
		i = 1;
		j = m;
		while (((i - j) != 1) && ((i - j) != -1))
		{
			kk = (i + j) / 2;
			if (v < m_nInputData2[kk - 1])
			{
				j = kk;
			}
			else
			{
				i = kk;
			}
		}
		iq = i - 3;
		iqq = i + 4;
	}

	if (iq < 1)
	{
		iq = 1;
	}

	if (iqq > m)
	{
		iqq = m;
	}

	for (i = ip - 1; i <= ipp - 1; i++)
	{
		b[i - ip + 1] = 0;
		for (j = iq - 1; j <= iqq - 1; j++)
		{
			h = m_nInputData3[j*M3+i];
			for (k = iq - 1; k <= iqq - 1; k++)
			{
				if (k != j)
				{
					h = h*(v - m_nInputData2[k]) / (m_nInputData2[j] - m_nInputData2[k]);
				}
			}
			b[i - ip + 1] = b[i - ip + 1] + h;
		}
	}

	f = 0;
	for (i = ip - 1; i <= ipp - 1; i++)
	{
		h = b[i - ip + 1];
		for (j = ip - 1; j <= ipp - 1; j++)
		{
			if (j != i)
			{
				h = h*(u - m_nInputData1[j]) / (m_nInputData1[i] - m_nInputData1[j]);
			}
		}
		f = f + h;
	}
	
	m_nOutputData[0] = f;
}