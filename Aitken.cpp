#include "mex.h"
#include "matrix.h"
#include "math.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex Aitken.cpp
N = 8;
x = [1:N]
y = randi(6,[1 N])
t = 1.3;
Aitken(x,y,t)
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
	double xx[10], yy[10];
	double z = 0;
	double eps = 0.000001;
	if (n < 1)
	{
		m_nOutputData[0] = z;
	}
	if (n == 1)
	{
		z = m_nInputData2[0];
		m_nOutputData[0] = z;
	}
	m = 10;
	if (m > n)
	{
		m = n;
	}
	if (t <= m_nInputData1[0])
	{
		k = 1;
	}
	else if (t >= m_nInputData1[n - 1])
	{
		k = n;
	}
	else
	{
		k = 1;
		j = n;
		while ((k - j != 1) && (k - j != -1))
		{
			l = (k + j) / 2;
			if (t < m_nInputData1[l - 1])
			{
				j = l;
			}
			else
			{
				k = l;
			}
		}
		if (fabs(t - m_nInputData1[l - 1])>fabs(t - m_nInputData1[j - 1]))
		{
			k = j;
		}
	}
	j = 1;
	l = 0;
	for (i = 1; i <= m; i++)
	{
		k = k + j*l;
		if ((k<1) || (k>n))
		{
			l = l + 1;
			j = -j;
			k = k + j*l;
		}
		xx[i - 1] = m_nInputData1[k - 1];
		yy[i - 1] = m_nInputData2[k - 1];
		l = l + 1;
		j = -j;
	}
	i = 0;
	do
	{
		i = i + 1;
		z = yy[i];
		for (j = 0; j <= i - 1; j++)
		{
			z = yy[j] + (t - xx[j])*(yy[j] - z) / (xx[j] - xx[i]);
		}
		yy[i] = z;
	} while ((i != m - 1) && (fabs(yy[i] - yy[i]) > eps));
	m_nOutputData[0] = z;
}