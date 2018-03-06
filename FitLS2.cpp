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
	plhs[0] = mxCreateDoubleMatrix(1, m + 1, mxREAL);
	m_nOutputData1 = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL);
	m_nOutputData2 = mxGetPr(plhs[1]);


	int n = N1;

	int i, j, k, l, kk;
	double apx[20], apy[20], bx[20], by[20], u[20][20];
	double t[20], t1[20], t2[20];
	double xx, yy, d1, d2, g, g1, g2;
	double x2, dd, y1, x1;
	double **v;

	for (i = 0; i <= p; i++)
	{
		for (j = 0; j <= q; j++)
		{
			a[i][j] = 0;
		}
	}

	if (p > n - 1)
	{
		p = n - 1;
	}

	if (p > 19)
	{
		p = 19;
	}

	if (q > m - 1)
	{
		q = m - 1;
	}

	if (q > 19)
	{
		q = 10;
	}

	xx = 0;
	for (i = 0; i <= n - 1; i++)
	{
		xx = xx + x[i] / (1.0*n);
	}

	yy = 0;
	for (i = 0; i <= m - 1; i++)
	{
		yy = yy + y[i] / (1.0*m);
	}

	d1 = 1 * n;
	apx[0] = 0;
	for (i = 0; i <= n - 1; i++)
	{
		apx[0] = apx[0] + x[i] - xx;
	}

	apx[0] = apx[0] / d1;
	for (j = 0; j <= m - 1; j++)
	{
		v[0][j] = 0;
		for (i = 0; i <= n - 1; i++)
		{
			v[0][j] = v[0][j] + z[i][j];
		}
		v[0][j] = v[0][j] / d1;
	}

	if (p > 0)
	{
		d2 = 0;
		apx[1] = 0;
		for (i = 0; i <= n - 1; i++)
		{
			g = x[i] - xx - apx[0];
			d2 = d2 + g*g;
			apx[1] = apx[1] + (x[i] - xx)*g*g;
		}

		apx[1] = apx[1] / d2;
		bx[1] = d2 / d1;
		for (j = 0; j <= m - 1; j++)
		{
			v[i][j] = 0;
			for (i = 0; i <= n - 1; i++)
			{
				g = x[i] - xx - apx[0];
				v[1][j] = v[1][j] + z[i][j] * g;
			}
			v[l][j] = v[1][j] / d2;
		}
		d1 = d2;
	}

	for (k = 2; k <= p; k++)
	{
		d2 = 0;
		apx[k] = 0;
		for (j = 0; j <= m - 1; j++)
		{
			v[k][j] = 0;
		}
		for (i = 0; i <= n - 1; i++)
		{
			g1 = 1;
			g2 = x[i] - xx - apx[0];
			for (j = 2; j <= k; j++)
			{
				g = (x[i] - xx - apx[j - 1])*g2 - bx[j - 1] * g1;
				g1 = g2;
				g2 = g;
			}
			d2 = d2 + g*g;
			apx[k] = apx[k] + (x[i] - xx)*g*g;
			for (j = 0; j <= m - 1; j++)
			{
				v[k][j] = v[k][j] + z[i][j] * g;
			}
		}

		for (j = 0; j <= m - 1; j++)
		{
			v[k][j] = v[k][j] / d2;
		}
		apx[k] = apx[k] / d2;
		bx[k] = d2 / d1;
		d1 = d2;
	}

	d1 = m;
	apy[0] = 0;
	for (i = 0; i <= m - 1; i++)
	{
		apy[0] = apy[0] + y[i] - yy;
	}
	apy[0] = apy[0] / d1;
	for (j = 0; j <= p; j++)
	{
		u[j][0] = 0;
		for (i = 0; i <= m - 1; i++)
		{
			u[j][0] = u[j][0] + v[j][i];
		}
		u[j][0] = u[j][0] / d1;
	}

	if (q > 1)
	{
		d2 = 0;
		apy[1] = 0;
		for (i = 0; i <= m - 1; i++)
		{
			g = y[i] - yy - apy[0];
			d2 = d2 + g*g;
			apy[1] = apy[1] + (y[i] - yy)*g*g;
		}
		apy[1] = apy[1] / d2;
		by[1] = d2 / d1;
		for (j = 0; j <= p; j++)
		{
			u[j][1] = 0;
			for (i = 0; i <= m - 1; i++)
			{
				g = y[i] - yy - apy[0];
				u[j][1] = u[j][1] + v[j][i] * g;
			}
			u[j][1] = u[j][1] / d2;
		}
		d1 = d2;
	}

	for (k = 2; k <= q; k++)
	{
		d2 = 0;
		apy[k] = 0;
		for (j = 0; j <= p; j++)
		{
			u[j][k] = 0;
		}

		for (i = 0; i <= m - 1; i++)
		{
			g1 = 1;
			g2 = y[i] - yy - apy[0];
			for (j = 2; j <= k; j++)
			{
				g = (y[i] - yy - apy[j - 1])*g2 - by[j - 1] * g1;
				g1 = g2;
				g2 = g;
			}
			d2 = d2 + g*g;
			apy[k] = apy[k] + (y[i] - yy)*g*g;
			for (j = 0; j <= p; j++)
			{
				u[j][k] = u[j][k] + v[j][i] * g;
			}
		}

		for (j = 0; j <= p; j++)
		{
			u[j][k] = u[j][k] / d2;
		}
		apy[k] = apy[k] / d2;
		by[k] = d2 / d1;
		d1 = d2;
	}
	v[0][0] = 1;
	v[1][0] = -apy[0];
	v[1][1] = 1;

	for (i = 0; i <= p; i++)
	{
		for (j = 0; j <= q; j++)
		{
			a[i][j] = 0;
		}
	}

	for (i = 2; i <= q; i++)
	{
		v[i][i] = v[i - 1][i - 1];
		v[i][i - 1] = -apy[i - 1] * v[i - 1][i - 1] + v[i - 1][i - 2];
		if (i >= 3)
		{
			for (k = i - 2; k >= 1; k--)
			{
				v[i][k] = -apy[i - 1] * v[i - 1][k] + v[i - 1][k - 1] - by[i - 1] * v[i - 2][k];
			}
		}

		v[i][0] = -apy[i - 1] * v[i - 1][0] - by[i - 1] * v[i - 2][0];
	}

	for (i = 0; i <= p; i++)
	{
		if (i == 0)
		{
			t[0] = 1;
			t1[0] = 1;
		}
		else
		{
			if (i == 1)
			{
				t[0] = -apx[0];
				t[1] = 1;
				t2[0] = t[0];
				t2[1] = t[1];
			}
			else
			{
				t[i] = t2[i - 1];
				t[i - 1] = -apx[i - 1] * t2[i - 1] + t2[i - 2];
				if (i >= 3)
				{
					for (k = i - 2; k >= 1; k--)
					{
						t[k] = -apx[i - 1] * t2[k] + t2[k - 1] - bx[i - 1] * t1[k];
					}
				}
				t[0] = -apx[i - 1] * t2[0] - bx[i - 1] * t1[0];
				t2[i] = t[i];
				for (k = i - 1; k >= 0; k--)
				{
					t1[k] = t2[k];
					t2[k] = t[k];
				}
			}
		}

		for (j = 0; j <= q; j++)
		{
			for (k = i; k >= 0; k--)
			{
				for (l = j; l >= 0; l--)
				{
					a[k][l] = a[k][l] + u[i][j] * t[k] * v[j][l];
				}
			}
		}
	}

	dt[0] = 0;
	dt[1] = 0;
	dt[2] = 0;
	for (i = 0; i <= n - 1; i++)
	{
		x1 = x[i] - xx;

		for (j = 0; j <= m - 1; j++)
		{
			y1 = y[i] - yy;
			x2 = 1;
			dd = 0;
			for (k = 0; k <= p; k++)
			{
				g = a[k][q];
				for (kk = q - 1; kk >= 0; kk--)
				{
					g = g*y1 + a[k][kk];
				}
				g = g*x2;
				dd = dd + g;
				x2 = x2*x1;
			}
			dd = dd - z[i][j];
			if (fabs(dd) > dt[2])
			{
				dt[2] = fabs(dd);
			}
			dt[0] = dt[0] + dd*dd;
			dt[1] = dt[1] + fabs(dd);
		}
	}

	return 1;
}