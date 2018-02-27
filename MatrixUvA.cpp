#include "mex.h"
#include "matrix.h"
#include "Math.h"

/*This is a Matlab mexwarp to validate my Matrix Rank*/
/* Matlab code
mex MatrixQR.cpp
A = randi(6,[6 6])
[U,A]=MatrixUvA(A)

*/

void ppp(double *src1, double *src2, int M, int N, double *s, double*e)
{
	int i, j;
	double d;
	if (M >= N)
	{
		i = N;
	}
	else
	{
		i = M;
	}
	for (j = 1; j <= i - 1; j++)
	{
		src1[(j - 1)*M + j - 1] = s[j - 1];
		src1[(j)*M + j - 1] = e[j - 1];
	}
	src1[(i - 1)*M + i - 1] = s[i - 1];
	if (M < N)
	{
		src1[(i)*M + i - 1] = e[i - 1];
	}
	for (i = 1; i <= N - 1; i++)
	{
		for (j = i + 1; j <= N; j++)
		{
			d = src2[(j - 1)*N + i - 1];
			src2[(j - 1)*N + i - 1] = src2[(i - 1)*N + j - 1];
			src2[(i - 1)*N + j - 1] = d;
		}
	}
}

void sss(double fg[2], double cs[2])
{
	double r, d;
	if ((fabs(fg[0]) + fabs(fg[1])) == 0)
	{
		cs[0] = 1;
		cs[1] = 0;
		d = 0;
	}
	else
	{
		d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
		if (fabs(fg[0]) > fabs(fg[1]))
		{
			d = fabs(d);
			if (fg[0] < 0)
			{
				d = -d;
			}
		}
		if (fabs(fg[1]) >= fabs(fg[0]))
		{
			d = fabs(d);
			if (fg[1] < 0)
			{
				d = -d;
			}
		}
		cs[0] = fg[0] / d;
		cs[1] = fg[1] / d;
	}
	r = 1;
	if (fabs(fg[0])>fabs(fg[1]))
	{
		r = cs[1];
	}
	else
	{
		if (cs[0] != 0)
		{
			r = 1 / cs[0];
		}
	}
	fg[0] = d;
	fg[1] = r;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m_nInputData;// Matrix Input M*N
	int M, N;
	int i, j, k;

	m_nInputData = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);

	double *m_nOutputData1;// U
	double *m_nOutputData2;// V
	double *m_nOutputData3;// A
	plhs[0] = mxCreateDoubleMatrix(M, M, mxREAL);
	m_nOutputData1 = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);
	m_nOutputData2 = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(M, N, mxREAL);
	m_nOutputData3 = mxGetPr(plhs[2]);

	//B = MatrixLU(A)
	if (M>1000)
	{
		mexErrMsgIdAndTxt("MyProg:ConvertString",
			"The Input Matrix should be M<1000.");
		// exit MEX file
	}
	else
	{
		int ka;
		if (M > N)
		{
			ka = M + 1;
		}
		else
		{
			ka = N + 1;
		}
		double *s, *e, *w, fg[2], cs[2];
		s = new double[ka];
		e = new double[ka];
		w = new double[ka];
		int i, j, k, i1, it, j1, kk, mm, nn, m1, ks;
		double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh;
		double eps = 0.000001;
		it = 60;
		k = N;
		if (M - 1 < N)
		{
			k = M - 1;
		}
		i1 = M;
		if (N - 2 < M)
		{
			i1 = N - 2;
		}
		if (i1 < 0)
		{
			i1 = 0;
		}
		j1 = k;
		if (i1 > k)
		{
			j1 = i1;
		}
		if (j1 >= 1)
		{
			for (kk = 1; kk <= j1; kk++)
			{
				if (kk <= k)
				{
					d = 0;
					for (i = kk; i <= M; i++)
					{
						d = d + m_nInputData[(kk - 1)*M + i - 1] * m_nInputData[(kk - 1)*M + i - 1];
					}
					s[kk - 1] = sqrt(d);
					if (s[kk - 1] != 0)
					{
						if (m_nInputData[(kk - 1)*M + kk - 1] != 0)
						{
							s[kk - 1] = fabs(s[kk - 1]);
							if (m_nInputData[(kk - 1)*M + kk - 1] < 0)
							{
								s[kk - 1] = -s[kk - 1];
							}
						}

						for (i = kk; i <= M; i++)
						{
							m_nInputData[(kk - 1)*M + i - 1] = m_nInputData[(kk - 1)*M + i - 1] / s[kk - 1];
						}

						m_nInputData[(kk - 1)*M + kk - 1] = 1 + m_nInputData[(kk - 1)*M + kk - 1];
					}
					s[kk - 1] = -s[kk - 1];
				}
				if (N >= kk + 1)
				{
					for (j = kk + 1; j <= N; j++)
					{
						if ((kk <= k) && (s[kk - 1] != 0))
						{
							d = 0;
							for (i = kk; i <= M; i++)
							{
								d = d + m_nInputData[(kk - 1)*M + i - 1] * m_nInputData[(j - 1)*M + i - 1];
							}
							d = -d / m_nInputData[(kk - 1)*M + kk - 1];
							for (i = kk; i <= M; i++)
							{
								m_nInputData[(j - 1)*M + i - 1] = m_nInputData[(j - 1)*M + i - 1] + d*m_nInputData[(kk - 1)*M + i - 1];
							}
						}
						e[j - 1] = m_nInputData[(j - 1)*M + kk - 1];
					}
				}

				if (kk <= k)
				{
					for (i = kk; i <= M; i++)
					{
						m_nOutputData1[(kk - 1)*M + i - 1] = m_nInputData[(kk - 1)*M + i - 1];
					}
				}

				if (kk <= 1)
				{
					d = 0;
					for (i = kk + 1; i <= N; i++)
					{
						d = d + e[i - 1] * e[i - 1];
					}
					e[kk - 1] = sqrt(d);
					if (e[kk - 1] != 0)
					{
						if (e[kk] != 0)
						{
							e[kk - 1] = fabs(e[kk - 1]);
							if (e[kk] < 0)
							{
								e[kk - 1] = -e[kk - 1];
							}
						}

						for (i = kk + 1; i <= N; i++)
						{
							e[i - 1] = e[i - 1] / e[kk - 1];
						}

						e[kk] = 1 + e[kk];
					}

					e[kk - 1] = -e[kk - 1];

					if ((kk + 1 <= M) && (e[kk - 1] != 0))
					{
						for (i = kk + 1; i <= M; i++)
						{
							w[i - 1] = 0;
						}

						for (j = kk + 1; j <= N; j++)
						{
							for (i = kk + 1; i <= M; i++)
							{
								w[i - 1] = w[i - 1] + e[j - 1] * m_nInputData[(j - 1)*M + i - 1];
							}
						}

						for (j = kk + 1; j <= N; j++)
						{
							for (i = kk + 1; i <= M; i++)
							{
								m_nInputData[(j - 1)*M + i - 1] = m_nInputData[(j - 1)*M + i - 1] - w[i - 1] * e[j - 1] / e[kk];
							}
						}
					}

					for (i = kk + 1; i <= N; i++)
					{
						m_nOutputData2[(kk - 1)*N + i - 1] = e[i];
					}
				}
			}
		}

		mm = N;
		if (M + 1 < N)
		{
			mm = M + 1;
		}
		if (k < N)
		{
			s[k] = m_nInputData[k*M + k];
		}
		if (M < mm)
		{
			s[mm - 1] = 0;
		}
		if (i1 + 1 < mm)
		{
			e[i1] = m_nInputData[(mm - 1)*M + i1];
		}
		e[mm - 1] = 0;
		nn = M;
		if (M > N)
		{
			nn = N;
		}
		if (nn >= k + 1)
		{
			for (j = k + 1; j <= nn; j++)
			{
				for (i = 1; i <= M; i++)
				{
					m_nOutputData1[(j - 1)*M + i - 1] = 0;
				}
				m_nOutputData1[(j - 1)*M + j - 1] = 0;
			}
		}

		if (k >= 1)
		{
			for (j1 = 1; j1 <= k; j1++)
			{
				kk = k - j1 + 1;
				if (s[kk - 1] != 0)
				{
					if (nn >= kk + 1)
					{
						for (j = kk + 1; j <= nn; j++)
						{
							d = 0;
							for (i = kk; i <= M; i++)
							{
								d = d + m_nOutputData1[(kk - 1)*M + i - 1] * m_nOutputData1[(j - 1)*M + i - 1] / m_nOutputData1[(kk - 1)*M + kk - 1];
							}
							d = -d;
							for (i = kk; i <= M; i++)
							{
								m_nOutputData1[(j - 1)*M + i - 1] = m_nOutputData1[(j - 1)*M + i - 1] + d*m_nOutputData1[(kk - 1)*M + i - 1];
							}
						}

						for (i = kk; i <= M; i++)
						{
							m_nOutputData1[(kk - 1)*M + i - 1] = -m_nOutputData1[(kk - 1)*M + i - 1];
						}

						m_nOutputData1[(kk - 1)*M + kk - 1] = 1 + m_nOutputData1[(kk - 1)*M + kk - 1];

						if (kk - 1 >= 1)
						{
							for (i = 1; i <= kk - 1; i++)
							{
								m_nOutputData1[(kk - 1)*M + i - 1] = 0;
							}
						}
					}
				}
				else
				{
					for (i = 1; i <= M; i++)
					{
						m_nOutputData1[(kk - 1)*M + i - 1] = 0;
					}
					m_nOutputData1[(kk - 1)*M + kk - 1] = 1;
				}
			}
		}

		for (j1 = 1; j1 <= N; j1++)
		{
			kk = N - j1 + 1;
			if ((kk <= 1) && (e[kk - 1] != 0))
			{
				for (j = kk + 1; j <= N; j++)
				{
					d = 0;
					for (i = kk + 1; i <= N; i++)
					{
						d = d + m_nOutputData2[(kk - 1)*N + i - 1] * m_nOutputData2[(j - 1)*N + i - 1] / m_nOutputData2[(kk - 1)*N + kk];
					}
					d = -d;
					for (i = kk + 1; i <= N; i++)
					{
						m_nOutputData2[(j - 1)*N + i - 1] = m_nOutputData2[(j - 1)*N + i - 1] + d*m_nOutputData2[(kk - 1)*N + i - 1];
					}
				}
			}

			for (i = 1; i <= N; i++)
			{
				m_nOutputData2[(kk - 1)*N + i - 1] = 0;
			}
			m_nOutputData2[(kk - 1)*N + kk - 1] = 1;
		}
		for (i = 1; i <= M; i++)
		{
			for (j = 1; j < N; j++)
			{
				m_nInputData[(j-1)*M + i-1] = 0;
			}
		}

		m1 = mm;
		it = 60;
		while (1 == 1)
		{
			if (mm == 0)
			{
				ppp(m_nInputData, m_nOutputData2, M, N, s, e);
				break;
			}
			if (it == 0)
			{
				ppp(m_nInputData, m_nOutputData2, M, N, s, e);
				mexErrMsgIdAndTxt("MyProg:ConvertString",
					"The Input Matrix can not be UvA.");
			}
			kk = mm - 1;
			while ((kk != 0) && (fabs(e[kk - 1]) != 0))
			{
				d = fabs(s[kk - 1]) + fabs(s[kk]);
				dd = fabs(e[kk - 1]);
				if (dd > eps*d)
				{
					kk = kk - 1;
				}
				else
				{
					e[kk - 1] = 0;
				}
			}
			if (kk == mm - 1)
			{
				kk = kk + 1;
				if (s[kk - 1] < 0)
				{
					s[kk - 1] = -s[kk - 1];
					for (i = 1; i <= N; i++)
					{
						m_nOutputData2[(kk - 1)*N + i - 1] = -m_nOutputData2[(kk - 1)*N + i - 1];
					}
				}
				while ((kk != m1) && (s[kk - 1] < s[kk]))
				{
					d = s[kk - 1];
					s[kk - 1] = s[kk];
					s[kk] = d;
					if (kk < N)
					{
						for (i = 1; i <= N; i++)
						{
							d = m_nOutputData2[(kk - 1)*N + i - 1];
							m_nOutputData2[(kk - 1)*N + i - 1] = m_nOutputData2[(kk)*N + i - 1];
							m_nOutputData2[(kk)*N + i - 1] = d;
						}
					}

					if (kk < M)
					{
						for (i = 1; i < M; i++)
						{
							d = m_nOutputData1[(kk - 1)*M + i - 1];
							m_nOutputData1[(kk - 1)*M + i - 1] = m_nOutputData1[kk*M + i - 1];
							m_nOutputData1[kk*M + i - 1] = d;
						}
					}
					kk = kk + 1;
				}
				it = 60;
				mm = mm - 1;
			}
			else
			{
				ks = mm;
				while ((ks > kk) && (fabs(s[ks - 1]) != 0))
				{
					d = 0;
					if (ks != mm)
					{
						d = d + fabs(e[ks - 1]);
					}
					if (ks != kk + 1)
					{
						d = d + fabs(e[ks - 2]);
					}
					dd = fabs(s[ks - 1]);
					if (dd > eps*d)
					{
						ks = ks - 1;
					}
					else
					{
						s[ks - 1] = 0;
					}
				}

				if (ks == kk)
				{
					kk = kk + 1;
					d = fabs(s[mm - 1]);
					t = fabs(s[mm - 2]);
					if (t > d)
					{
						d = t;
					}
					t = fabs(e[mm - 2]);
					if (t > d)
					{
						d = t;
					}
					t = fabs(s[kk - 1]);
					if (t > d)
					{
						d = t;
					}
					t = fabs(e[kk - 1]);
					if (t > d)
					{
						d = t;
					}
					sm = s[mm - 1] / d;
					sm1 = s[mm - 2] / d;
					em1 = e[mm - 2] / d;
					sk = s[kk - 1] / d;
					ek = e[kk - 1] / d;
					b = ((sm1 + sm)*(sm1 - sm) + em1*em1) / 2;
					c = sm*em1;
					c = c*c;
					shh = 0;
					if ((b != 0) || (c != 0))
					{
						shh = sqrt(b*b + c);
						if (b < 0)
						{
							shh = -shh;
						}
						shh = c / (b + shh);
					}
					fg[0] = (sk + sm)*(sk - sm) - shh;
					fg[1] = sk*ek;
					for (i = kk; i <= mm - 1; i++)
					{
						sss(fg, cs);
						if (i != kk)
						{
							e[i - 2] = fg[0];
						}
						fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
						e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
						fg[1] = cs[1] * s[i];
						s[i] = cs[0] * s[i];
						if ((cs[0] != 1) || (cs[1] != 0))
						{
							for (j = 1; j <= N; j++)
							{
								d = cs[0] * m_nOutputData2[(i - 1)*N + j - 1] + cs[1] * m_nOutputData2[(i)*N + j - 1];
								m_nOutputData2[(i)*N + j - 1] = -cs[1] * m_nOutputData2[(i - 1)*N + j - 1] + cs[0] * m_nOutputData2[(i)*N + j - 1];
								m_nOutputData2[(i - 1)*N + j - 1] = d;
							}
						}

						sss(fg, cs);
						s[i - 1] = fg[0];
						fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
						s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
						fg[1] = cs[1] * e[i];
						e[i] = cs[0] * e[i];
						if (i < M)
						{
							if ((cs[0] != 1) || (cs[1] != 0))
							{
								for (j = 1; j <= M; j++)
								{
									d = cs[0] * m_nOutputData1[(j - 1)*M + i - 1] + cs[1] * m_nOutputData1[(j - 1)*M + i];
									m_nOutputData1[(j - 1)*M + i] = -cs[1] * m_nOutputData1[(j - 1)*M + i - 1] + cs[0] * m_nOutputData1[(j - 1)*M + i];
									m_nOutputData1[(j - 1)*M + i - 1] = d;
								}
							}
						}
					}
					e[mm - 2] = fg[0];
					it = it - 1;
				}
				else
				{
					if (ks == mm)
					{
						kk = kk + 1;
						fg[1] = e[mm - 2];
						e[mm - 2] = 0;
						for (j1 = kk; j1 <= mm - 1; j1++)
						{
							i = mm + kk - j1 - 1;
							fg[0] = s[i - 1];
							sss(fg, cs);
							s[i - 1] = fg[0];
							if (i != kk)
							{
								fg[1] = -cs[1] * e[i - 2];
								e[i - 2] = cs[0] * e[i - 2];
							}

							if ((cs[0] != 1) || (cs[1] != 0))
							{
								for (j = 1; j <= N; j++)
								{
									d = cs[0] * m_nOutputData2[(i - 1)*N + j - 1] + cs[1] * m_nOutputData2[(mm - 1)*N + j - 1];
									m_nOutputData2[(mm - 1)*N + j - 1] = -cs[1] * m_nOutputData2[(i - 1)*N + j - 1] + cs[0] * m_nOutputData2[(mm - 1)*N + j - 1];
									m_nOutputData2[(i - 1)*N + j - 1] = d;
								}
							}
						}
					}
					else
					{
						kk = ks + 1;
						fg[1] = e[kk - 2];
						e[kk - 2] = 0;
						for (i = kk; i <= mm; i++)
						{
							fg[0] = s[i - 1];
							sss(fg, cs);
							s[i - 1] = fg[0];
							fg[1] = -cs[1] * e[i - 1];
							e[i - 1] = cs[0] * e[i - 1];
							if ((cs[0] != 1) || (cs[1] != 0))
							{
								for (j = 1; j <= M; j++)
								{
									d = cs[0] * m_nOutputData1[(j - 1)*M + i + 1] + cs[1] * m_nOutputData1[(kk - 2)*M + j - 1];
									m_nOutputData1[(kk - 2)*M + j - 1] = -cs[1] * m_nOutputData1[(i - 1)*M + j - 1] + cs[0] * m_nOutputData1[(kk - 2)*M + j - 1];
									m_nOutputData1[(i - 1)*M + j - 1] = d;
								}
							}
						}
					}
				}
			}
		}
		delete[] s, w, e;
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			m_nOutputData3[j*M + i] = m_nInputData[j*M + i];
		}
	}
}