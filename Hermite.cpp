#include "mex.h"
#include "matrix.h"
#include "math.h"

/*This is a Matlab mexwarp to validate my Interpolation Algorithm*/
/*Matlab code
mex Hermite.cpp
N = 8;
x = [1:N]
y = randi(6,[1 N])
t = 1.3;
Hermite(x,y,dy,t)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData1;// Matrix Input M1*N1
	double *m_nInputData2;// Matrix Input M2*N2
	double *m_nInputData3;// Matrix Input M3*N3
	double *m_nInputData4;
	double *m_nOutputData;
	int M1, N1, M2, N2;

	m_nInputData1 = mxGetPr(prhs[0]);
	m_nInputData2 = mxGetPr(prhs[1]);
	m_nInputData3 = mxGetPr(prhs[2]);
	m_nInputData4 = mxGetPr(prhs[3]);
	M1 = mxGetM(prhs[0]);
	N1 = mxGetN(prhs[0]);
	M2 = mxGetM(prhs[1]);
	N2 = mxGetN(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	double t = m_nInputData4[0];
	int n = N1;

	if (N1 != N2)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input x y diminsion should be match.");
		// exit MEX file
	}

	int i, j;
	double p, q, s;
	double z = 0;
	for (i = 1; i <= n; i++)
	{
		s = 1.0;
		for (j = 1; j <= n; j++)
		{
			if (j != i)
			{
				s = s*(t - m_nInputData1[j - 1]) / (m_nInputData1[i - 1] - m_nInputData1[j - 1]);
			}
		}
		s = s*s;
		p = 0;
		for (j = 1; j <= n; j++)
		{
			if (j != i)
			{
				p = p + 1 / (m_nInputData1[i - 1] - m_nInputData1[j - 1]);
			}
		}
		q = m_nInputData2[i - 1] + (t - m_nInputData1[i - 1])*(m_nInputData3[i - 1] - 2.0*m_nInputData2[i - 1] * p);
		z = z + q*s;
	}
	m_nOutputData[0] = z;
}