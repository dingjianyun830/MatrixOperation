#include "mex.h"
#include "matrix.h"

/*This is a Matlab mexwarp to validate my Matrix Product*/
/*Matlab code
mex Prod2.cpp
A = randi(6,[5 5])
Prod2(A,2.5)
A*2.5
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData1;// Matrix Input M1*N1
	double *m_nInputData2;// Matrix Input M2*N2
	double *m_nOutputData;
	int M1, N1, M2, N2;

	m_nInputData1 = mxGetPr(prhs[0]);
	m_nInputData2 = mxGetPr(prhs[1]);
	M1 = mxGetM(prhs[0]);
	N1 = mxGetN(prhs[0]);
	M2 = mxGetM(prhs[1]);
	N2 = mxGetN(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(M1, N2, mxREAL);
	m_nOutputData = mxGetPr(plhs[0]);

	//C = Prod(A,B)
	if (N1 != M2)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input diminsion should be match.");
		// exit MEX file
	}
	else
	{
		for (int i = 0; i < M1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				for (int k = 0; k < N1; k++)
				{
					m_nOutputData[j*M1 + i] += m_nInputData1[k*M1 + i] * m_nInputData2[j*M2 + k];
				}
			}
		}
	}
}