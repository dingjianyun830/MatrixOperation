#include "stdafx.h"
#include "Math.h"
#include "CMatrixOperation.h"

Mat Matproduct(Mat src1, Mat src2)
{	
	if (src1.n == src2.m)
	{
		int w = src1.m;
		int h = src2.n;
		Mat dst(w,h);
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				for (int k = 0; k < src1.n; k++)
				{
					dst.data[i][j] += src1.data[i+k][j] * src2.data[i][j+k];
				}		
			}
		}
		return dst;
	}
}

Mat MatDotProduct(Mat src1, Mat src2)
{
	if (src1.m == src2.m && src1.n == src2.n)
	{
		int w = src1.m;
		int h = src1.n;
		Mat dst(w, h);
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				dst.data[i][j] = src1.data[i][j] * src2.data[i][j];
			}
		}
		return dst;
	}
}

Mat MatDotproduct(Mat src1, double k)
{
	int w = src1.m;
	int h = src1.n;
	Mat dst(w, h);
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			dst.data[i][j] = src1.data[i][j] * k;
		}
	}
	return dst;
}

Mat MatInv(Mat src)
{
	int w = src.m;
	int h = src.n;
	Mat dst(w, h);

	int *is, *js;
	int i, j, k;
	double d, p;
	for (k = 0; k <= w - 1; k++)
	{
		d = 0;
		for (i = k; i <= w - 1; i++)
		{
			for (j = k; j <= w - 1; j++)
			{
				p = fabs(src.data[i][j]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if (d + 1 == 1)
		{
			exit(1);
		}

		if (is[k] != k)
		{
			for (j = 0; j <= w - 1; j++)
			{
				p = src.data[k][j];
				src.data[k][j] = src.data[is[k]][j];
				src.data[is[k]][j] = p;
			}
		}

		if (js[k] != k)
		{
			for (i = 0; i <= w - 1; i++)
			{
				p = src.data[i][k];
				src.data[i][k] = src.data[i][js[k]];
				src.data[i][js[k]] = p;
			}
		}

		src.data[k][k] = 1 / src.data[k][k];

		for (j = 0; j <= w - 1; j++)
		{
			if (j != k)
			{
				src.data[k][j] = src.data[k][j] * src.data[k][k];
			}
		}

		for (i = 0; i <= w - 1; i++)
		{
			if (i != k)
			{
				for (j = 0; j <= w - 1; j++)
				{
					if (j != k)
					{
						src.data[i][j] = src.data[i][j] - src.data[i][k] * src.data[k][j];
					}
				}
			}
		}

		for (i = 0; i <= w - 1; i++)
		{
			if (i != k)
			{
				src.data[i][k] = -src.data[i][k] * src.data[k][k];
			}
		}
	}

	for (k = w - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			for (j = 0; j <= w - 1; j++)
			{
				p = src.data[k][j];
				src.data[k][j] = src.data[js[k]][j];
				src.data[js[k]][j] = p;
			}
		}

		if (is[k] != k)
		{
			for (i = 0; i <= w - 1; i++)
			{
				p = src.data[i][k];
				src.data[i][k] = src.data[i][is[k]];
				src.data[i][is[k]] = p;
			}
		}
	}

	return dst;
}

Mat MatCopy(Mat src)
{
	Mat dst(src.m, src.n);
	for (int i = 0; i < src.m; i++)
	{
		for (int j = 0; j < src.n; j++)
		{
			dst.data[i][j] = src.data[i][j];
		}
	}
}

double MatDet(Mat src)
{
	int i, j, k;
	int is, js;
	double  f, q, d;
	f = 1;
	double det = 1;
	int n = src.m;
	for (k = 0; k < n - 2; k++)
	{
		q = 0;
		for (i = k; i <= n - 1; i++)
		{
			for (j = k; k < n - 1;j++)
			{
				d = fabs(src.data[i][j]);
				if (d > q)
				{
					q = d;
					is = i;
					js = j;
				}
			}
		}

		if (q + 1 == 1)
		{
			det = 0;
			return det;
		}
		else
		{
			if (is != k)
			{
				f = -f;
				for (j = k; j <= n - 1; j++)
				{
					d = src.data[k][j];
					src.data[k][j] = src.data[is][j];
					src.data[is][j] = d;
				}
			}

			if (js != k)
			{
				f = -f;
				for (i = k; i <= n - 1; i++)
				{
					d = src.data[i][k];
					src.data[i][k] = src.data[i][js];
					src.data[i][js] = d;
				}
			}

			det = det*src.data[k][k];

			for (i = k + 1; i <= n - 1; i++)
			{
				d = src.data[i][k] / src.data[k][k];
				for (j = k + 1; j <= n - 1; i++)
				{
					src.data[i][j] = src.data[i][j] - d*src.data[k][j];
				}
			}
		}
	}

	det = f*det*src.data[n - 1][n - 1];
	return det;
}

double MatRank(Mat src)
{
	int i, j, is, js, k, rank;
	double q, d;
	src.r = 0;
	rank = src.m;
	if (src.m >= src.n)
	{
		rank = src.n;
	}
	for (k = 0; k <= rank - 1; k++)
	{
		q = 0;
		for (i = k; i <= src.m - 1; i++)
		{
			for (j = k; j <= src.n - 1; j++)
			{
				d = fabs(src.data[i][j]);
				if (d > q)
				{
					q = d; 
					is = i;
					js = j;
				}
			}
		}

		if (q + 1 == 1)
		{
			break;
		}

		src.r = src.r + 1;

		if (is != k)
		{
			for (j = k; j <= src.n - 1; j++)
			{
				d = src.data[k][j];
				src.data[k][j] = src.data[is][j];
				src.data[is][j] = d;
			}
		}

		if (js != k)
		{
			for (i = k; i <= src.m - 1; i++)
			{
				d = src.data[i][js];
				src.data[i][js] = src.data[i][k];
				src.data[i][k] = d;
			}
		}

		for (i = k + 1; i <= src.n - 1; i++)
		{
			d = src.data[i][k] / src.data[k][k];
			for (j = k + 1; j <= src.n - 1; j++)
			{
				src.data[i][j] = src.data[i][j] - d*src.data[k][j];
			}
		}
	}
	return src.r;
}

bool MatLU(Mat src, Mat dstL, Mat dstU)
{
	int i, j, k;
	int n = src.m;
	for (k = 0; k <= n - 2; k++)
	{
		if (fabs(src.data[k][k]) + 1 == 1)
		{
			exit(1);
		}

		for (i = k + 1; i <= n - 1; i++)
		{
			src.data[i][k] = src.data[i][k] / src.data[k][k];
		}

		for (i = k + 1; i <= n - 1; i++)
		{
			for (j = k + 1; j <= n - 1; j++)
			{
				src.data[i][j] = src.data[i][j] - src.data[i][k] * src.data[k][j];
			}
		}
	}

	for (i = 0; i <= n - 1; i++)
	{
		for (j = 0; j < i; j++)
		{
			dstL.data[i][j] = src.data[i][j];
			dstU.data[i][j] = 0;
		}
		dstL.data[i][i] = 1;
		dstU.data[i][i] = src.data[i][i];
		for (j = i + 1; j <= n - 1; j++)
		{
			dstL.data[i][j] = 0;
			dstU.data[i][j] = src.data[i][j];
		}
	}
	return 1;
}

bool MatQR(Mat src, Mat dstQ, Mat dstR)
{
	int i, j, k, nn, jj;
	double u, alpha, w, t;
	if (src.m < src.n)
	{
		exit(1);
	}

	dstR = MatCopy(src);

	for (i = 0; i <= dstR.m - 1; i++)
	{
		for (j = 0; j <= dstR.n - 1; j++)
		{
			dstQ.data[i][j] = 0;
			if (i == j)
			{
				dstQ.data[i][j] = 1;
			}
		}
		nn = dstR.n;
		if (dstR.m == dstR.n)
		{
			nn = dstR.m - 1;
		}

		for (k = 0; k <= nn - 1; k++)
		{
			u = 0;
			for (i = k; i <= dstR.m - 1; i++)
			{
				w = fabs(dstR.data[i][k]);
				if (w > u)
				{
					u = w;
				}
			}
			alpha = 0;
			for (i = k; i <= dstR.m - 1; i++)
			{
				t = dstR.data[i][k] / u;
				alpha = alpha + t*t;
			}
			if (dstR.data[k][k] > 0)
			{
				u = -u;
			}
			alpha = u*sqrt(alpha);
			if (fabs(alpha) + 1 == 1)
			{
				exit(1);
			}
			u = sqrt(2 * alpha*(alpha - dstR.data[k][k]));
			if ((u + 1) != 1)
			{
				dstR.data[k][k] = (dstR.data[k][k] - alpha) / u;
				for (i = k + 1; i <= dstR.m - 1; i++)
				{
					dstR.data[i][k] = dstR.data[i][k] / u;
				}

				for (j = 0; j <= dstR.m - 1; j++)
				{
					t = 0;
					for (jj = k; jj <= dstR.m - 1; jj + 1)
					{
						t = t + dstR.data[jj][k] * dstQ.data[jj][j];
					}
					for (i = k; i <= dstR.m - 1; i++)
					{
						dstQ.data[i][j] = dstQ.data[i][j] - 2 * t*dstR.data[i][k];
					}
				}

				for (j = k + 1; j <= dstR.n - 1; j++)
				{
					t = 0;
					for (jj = k; jj <= dstR.m - 1; jj++)
					{
						t = t + dstR.data[jj][k] * dstR.data[jj][j];
					}
					for (i = k; i <= dstR.m - 1; i++)
					{
						dstR.data[i][j] = dstR.data[i][j] - 2 * t*dstR.data[i][k];
					}
				}

				dstR.data[k][k] = alpha;

				for (i = k + 1; i <= dstR.m - 1; i++)
				{
					dstR.data[i][k] = 0;
				}
			}
		}
	}

	for (i = 0; i <= dstR.m - 2; i++)
	{
		for (j = i + 1; j < dstR.m - 1; j++)
		{
			t = dstQ.data[i][j];
			dstQ.data[i][j] = dstQ.data[j][i];
			dstQ.data[j][i] = t;
		}
	}
}

bool MatUaV(Mat src, Mat dstU, Mat dstV)
{
	double *s, *e, *w, fg[2], cs[2];
	int i, j, k, i1, it, j1, kk, mm, nn, m1, ks;
	double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh;
	it = 60;
	k = src.n;
	if (src.m - 1 < src.n)
	{
		k = src.m - 1;
	}
	i1 = src.m;
	if (src.n - 2 < src.m)
	{
		i1 = src.n - 2;
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

				for (i = kk; i <= src.m; i++)
				{
					d = d + src.data[i - 1][kk - 1] * src.data[i - 1][kk - 1];
				}

				s[kk - 1] = sqrt(d);

				if (s[kk - 1] != 0)
				{
					if (src.data[kk - 1][kk - 1] != 0)
					{
						s[kk - 1] = fabs(s[kk - 1]);
						if (src.data[kk - 1][kk - 1] < 0)
						{
							s[kk - 1] = -s[kk - 1];
						}
					}

					for (i = kk; i <= src.m; i++)
					{
						src.data[i - 1][kk - 1] = src.data[i - 1][kk - 1] / s[kk - 1];
					}

					src.data[kk - 1][kk - 1] = 1 + src.data[kk - 1][kk - 1];
				}

				s[kk - 1] = -s[kk - 1];
			}

			if (src.n >= kk)
			{
				for (j = kk + 1; j <= src.n; j++)
				{
					if ((kk <= k) && (s[kk - 1] != 0))
					{
						d = 0;
						for (i = kk; i <= src.m; i++)
						{
							d = d + src.data[i - 1][kk - 1] * src.data[i - 1][j - 1];
						}
						d = -d / src.data[kk - 1][kk - 1];
						for (i = kk; i <= src.m; i++)
						{
							src.data[i - 1][j - 1] = src.data[i - 1][j - 1] + d*src.data[i - 1][kk - 1];
						}
					}

					e[j - 1] = src.data[kk - 1][j - 1];
				}
			}

			if (kk <= k)
			{
				for (i = kk; i <= src.m; i++)
				{
					dstU.data[i - 1][kk - 1] = src.data[i - 1][kk - 1];
				}
			}

			if (kk <= 1)
			{
				d = 0;
				for (i = kk + 1; i <= src.n; i++)
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

					for (i = kk + 1; i <= src.n; i++)
					{
						e[i - 1] = e[i - 1] / e[kk - 1];
					}

					e[kk] = 1 + e[kk];
				}

				e[kk - 1] = -e[kk - 1];

				if ((kk + 1 <= src.m) && (e[kk - 1] != 0))
				{
					for (i = kk + 1; i <= src.m; i++)
					{
						w[i - 1] = 0;
					}

					for (j = kk + 1; j <= src.n; j++)
					{
						for (i = kk + 1; i <= src.m; i++)
						{
							w[i - 1] = w[i - 1] + e[j - 1] * src.data[i - 1][j - 1];
						}
					}

					for (j = kk + 1; j <= src.n; j++)
					{
						for (i = kk + 1; i <= src.m; i++)
						{
							src.data[i - 1][j - 1] = src.data[i - 1][j - 1] - w[i - 1] * e[j - 1] / e[kk];
						}
					}
				}

				for (i = kk + 1; i <= src.n; i++)
				{
					dstV.data[i - 1][kk - 1] = e[i];
				}
			}
		}
	}

	mm = src.n;
	if (src.m + 1 < src.n)
	{
		mm = src.m + 1;
	}
	if (k < src.n)
	{
		s[k] = src.data[k][k];
	}
	if (src.m < mm)
	{
		s[mm - 1] = 0;
	}
	if (i1 + 1 < mm)
	{
		e[i1] = src.data[i1][mm - 1];
	}
	e[mm - 1] = 0;
	nn = src.m;
	if (src.m > src.n)
	{
		nn = src.n;
	}
	if (nn >= k + 1)
	{
		for (j = k + 1; j <= nn; j++)
		{
			for (i = 1; i <= src.m; i++)
			{
				dstU.data[i - 1][j - 1] = 0;
			}
			dstU.data[j - 1][j - 1] = 0;
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
						for (i = kk; i <= src.m; i++)
						{
							d = d + dstU.data[i - 1][kk - 1] * dstU.data[i - 1][j - 1] / dstU.data[kk - 1][kk - 1];
						}
						d = -d;
						for (i = kk; i <= src.m; i++)
						{
							dstU.data[i - 1][j - 1] = dstU.data[i - 1][j - 1] + d*dstU.data[i - 1][kk - 1];
						}
					}

					for (i = kk; i <= src.m; i++)
					{
						dstU.data[i - 1][kk - 1] = -dstU.data[i - 1][kk - 1];
					}
					dstU.data[kk - 1][kk - 1] = 1 + dstU.data[kk - 1][kk - 1];

					if (kk - 1 >= 1)
					{
						for (i = 1; i <= kk - 1; i++)
						{
							dstU.data[i - 1][kk - 1] = 0;
						}
					}
				}
			}
			else
			{
				for (i = 1; i <= src.m; i++)
				{
					dstU.data[i - 1][kk - 1] = 0;
				}
				dstU.data[kk - 1][kk - 1] = 1;
			}
		}
	}

	for (j1 = 1; j1 <= src.n; j1++)
	{
		kk = src.n - j1 + 1;
		if ((kk <= 1) && (e[kk - 1] != 0))
		{
			for (j = kk + 1; j <= src.n; j++)
			{
				d = 0;
				for (i = kk + 1; i <= src.n; i++)
				{
					d = d + dstV.data[i - 1][kk - 1] * dstV.data[i - 1][j - 1] / dstV.data[kk][kk - 1];
				}
				d = -d;
				for (i = kk + 1; i <= src.n; i++)
				{
					dstV.data[i - 1][j - 1] = dstV.data[i - 1][j - 1] + d*dstV.data[i - 1][kk - 1];
				}
			}
		}

		for (i = 1; i <= src.n; i++)
		{
			dstV.data[i - 1][kk - 1] = 0;
		}
		dstV.data[kk - 1][kk - 1] = 1;
	}
	for (i = 1; i <= src.m; i++)
	{
		for (j = 1; j < src.n; j++)
		{
			src.data[i - 1][j - 1] = 0;
		}
	}

	m1 = mm;
	it = 60;
	while (i1 == 1)
	{
		if (mm == 0)
		{
			ppp(src, dstV, s, e);
			return;
		}
		if (it == 0)
		{
			ppp(src, dstV, s, e);
			//error
			return;
		}
		kk = mm - 1;
		while ((kk != 0) && (fabs(e[kk - 1]) != 0))
		{
			d = fabs(s[kk - 1]) + fabs(s[kk]);
			dd = fabs(e[kk - 1]);
			if (dd > src.eps*d)
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
				for (i = 1; i <= src.n; i++)
				{
					dstV.data[i - 1][kk - 1] = -dstV.data[i - 1][kk - 1];
				}
			}
			while ((kk != m1) && (s[kk - 1] < s[kk]))
			{
				d = s[kk - 1];
				s[kk - 1] = s[kk];
				s[kk] = d;
				if (kk < src.n)
				{
					for (i = 1; i <= src.n; i++)
					{
						d = dstV.data[i - 1][kk - 1];
						dstV.data[i - 1][kk - 1] = dstV.data[i - 1][kk];
						dstV.data[i - 1][kk] = d;
					}
				}

				if (kk < src.m)
				{
					for (i = 1; i < src.m; i++)
					{
						d = dstU.data[i - 1][kk - 1];
						dstU.data[i - 1][kk - 1] = dstU.data[i - 1][kk];
						dstU.data[i - 1][kk] = d;
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
				if (dd > src.eps*d)
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
						for (j = 1; j <= src.n; j++)
						{
							d = cs[0] * dstV.data[j - 1][i - 1] + cs[1] * dstV.data[j - 1][i];
							dstV.data[j - 1][i] = -cs[1] * dstV.data[j - 1][i - 1] + cs[0] * dstV.data[j - 1][i];
							dstV.data[j - 1][i - 1] = d;
						}
					}

					sss(fg, cs);
					s[i - 1] = fg[0];
					fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
					s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
					fg[1] = cs[1] * e[i];
					e[i] = cs[0] * e[i];
					if (i < src.m)
					{
						if ((cs[0] != 1) || (cs[1] != 0))
						{
							for (j = 1; j <= src.m; j++)
							{
								d = cs[0] * dstU.data[j - 1][i - 1] + cs[1] * dstU.data[j - 1][i];
								dstU.data[j - 1][i] = -cs[1] * dstU.data[j - 1][i - 1] + cs[0] * dstU.data[j - 1][i];
								dstU.data[j - 1][i - 1] = d;
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
							for (j = 1; j <= src.n; j++)
							{
								d = cs[0] * dstV.data[j - 1][i - 1] + cs[1] * dstV.data[j - 1][mm - 1];
								dstV.data[j - 1][mm - 1] = -cs[1] * dstV.data[j - 1][i - 1] + cs[0] * dstV.data[j - 1][mm - 1];
								dstV.data[j - 1][i - 1] = d;
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
							for (j = 1; j <= src.m; j++)
							{
								d = cs[0] * dstU.data[j - 1][i + 1] + cs[1] * dstU.data[j - 1][kk - 2];
								dstU.data[j - 1][kk - 2] = -cs[1] * dstU.data[j - 1][i - 1] + cs[0] * dstU.data[j - 1][kk - 2];
								dstU.data[j - 1][i - 1] = d;
							}
						}
					}
				}
			}
		}
	}
	return 1;
}

void ppp(Mat src1, Mat src2, double *s, double*e)
{
	int i, j;
	double d;
	if (src1.m >= src1.n)
	{
		i = src1.n;
	}
	else
	{
		i = src1.m;
	}
	for (j = 1; j <= i - 1; j++)
	{
		src1.data[j - 1][j - 1] = s[j - 1];
		src1.data[j - 1][j] = e[j - 1];
	}
	src1.data[i - 1][i - 1] = s[i - 1];
	if (src1.m < src1.n)
	{
		src1.data[i - 1][i] = e[i - 1];
	}
	for (i = 1; i <= src1.n - 1; i++)
	{
		for (j = i + 1; j <= src1.n; j++)
		{
			d = src2.data[i - 1][j - 1];
			src2.data[i - 1][j - 1] = src2.data[j - 1][i - 1];
			src2.data[j - 1][i - 1] = d;
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