#include "mex.h"
#include "matrix.h"

/*This is a Matlab mexwarp to validate my Matrix dotProduct*/
/*Matlab code
mex Prod1.cpp
A = randi(6, [5 5])
B = randi(6, [5 5])
Prod1(A,B)
A.*B
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData1;// Matrix Input M1*N1
	double *m_nInputData2;// Matrix Input M2*N2
	double *m_nOutputData;
	int M1, N1, M2, N2;
	int i, j;

	m_nInputData1 = mxGetPr(prhs[0]);
	m_nInputData2 = mxGetPr(prhs[1]);
	M1 = mxGetM(prhs[0]);
	N1 = mxGetN(prhs[0]);
	M2 = mxGetM(prhs[1]);
	N2 = mxGetN(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(M1, N1, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	//C = Prod(A,B)
	if (M2 == 1 && N2 == 1)
	{
		for (i = 0; i < M1; i++)
		{
			for (j = 0; j < N1; j++)
			{
				m_nOutputData[j*M1 + i] = m_nInputData1[j*M1 + i] * m_nInputData2[0];
			}
		}
	}
	else if (M1 == M2&&N1 == N2)
	{
		for (i = 0; i < M1; i++)
		{
			for (j = 0; j < N1; j++)
			{
				m_nOutputData[j*M1 + i] = m_nInputData1[j*M1 + i] * m_nInputData2[j*M1 + i];
			}
		}
	}
	else
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input diminsion should be match.");
		// exit MEX file
	}
}