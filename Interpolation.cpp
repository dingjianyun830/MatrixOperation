#include "Interpolation.h"

double Lagrange(double *x, double *y, int n, double t)
{
	int i, j, k, m;
	double s;
	double z = 0;

	if (n < 1)
	{
		z = 0;
		return z;
	}
	if (n == 1)
	{
		z = y[0];
		return z;
	}
	if (n == 2)
	{
		z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return z;
	}

	i = 0;
	while ((x[i] < t) && (i < n))
	{
		i = i + 1;
	}
	k = i - 4;
	if (k < 0)
	{
		k = 0;
	}
	m = i + 3;
	if (m > n - 1)
	{
		m = n - 1;
	}
	for (i = k; i <= m; i++)
	{
		s = 1;
		for (j = k; j <= m; j++)
		{
			if (j != i)
			{
				s = s*(t - x[j]) / (x[i] - x[j]);
			}
		}
		z = z + s*y[i];
	}
	return z;
}

double Pqrat(double *x, double *y, int n, double t)
{
	int i, j, k, m, l;
	double h, b[8];
	double z = 0;
	if (n < 1)
	{
		z = 0; 
		return z;
	}
	if (n == 1)
	{
		z = y[0];
		return z;
	}
	if (n <= 8)
	{
		k = 0;
		m = 8;
	}
	else if (t < x[4])
	{
		k = 0;
		m = 8;
	}
	else if (t > x[n - 5])
	{
		k = n - 8;
		m = 8;
	}
	else
	{
		k = 1;
		j = n;
		while (j - k != 1)
		{
			i = (k + j) / 2;
			if (t < x[i - 1])
			{
				j = i;
			}
			else
			{
				k = i;
			}
		}
		k = k - 4;
		m = 8;
	}
	b[0] = y[k];
	for (i = 2; i <= m; i++)
	{
		h = y[i + k - 1];
		l = 0;
		j = 1;
		while ((l == 0) && (j <= i - 1))
		{
			if (fabs(h - b[j - 1]) + 1 == 1)
			{
				l = 1;
			}
			else
			{
				h = (x[i + k - 1] - x[j + k - 1]) / (h - b[j - 1]);
			}
			j = j + 1;
		}
		b[i - 1] = h;
		if (l != 0)
		{
			b[i - 1] = 1.0e+35;
		}
	}
	z = b[m - 1];
	for (i = m - 1; i >= 1; i--)
	{
		z = b[i - 1] + (t - x[i + k - 1]) / z;
	}
	return z;
}

double Hermite(double *x, double *y, double *dy, int n, double t)
{
	int i, j;
	double p, q, s;
	double z = 0;
	for (i = 1; i <= n; i++)
	{
		s = 1;
		for (j = 1; j <= n; j++)
		{
			if (j != i)
			{
				s = s*(t - x[j - 1]) / (x[i - 1] - x[j - 1]);
			}
			s = s*s;
			p = 0;
			for (j = 1; j <= n; j++)
			{
				if (j != i)
				{
					p = p + 1 / (x[i - 1] - x[j - 1]);
				}
			}
			q = y[i - 1] + (t - x[i - 1])*(dy[i - 1] - 2.0*y[i - 1] * p);
			z = z + q*s;
		}
	}
	return z;
}

double Aitken(double *x, double*y, int n, double t) 
{
	int i, j, k, m, l;
	double xx[10], yy[10];
	double z = 0;
	if (n < 1)
	{
		return z;
	}
	if (n == 1)
	{
		z = y[0];
		return z;
	}
	m = 10;
	if (m > n)
	{
		m = n;
	}
	if (t <= x[0])
	{
		k = 1;
	}
	else if (t >= x[n - 1])
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
			if (t < x[l - 1])
			{
				j = l;
			}
			else
			{
				k = l;
			}
		}
		if (fabs(t - x[l - 1])>fabs(t - x[j - 1]))
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
		xx[i - 1] = x[k - 1];
		yy[i - 1] = y[k - 1];
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
	return z;
}

double Akima(double *x, double *y, int n, double t, int k)
{
	int m, kk;
	double u[5];
	double p, q;
	double z = 0;
	double s[4];
	if (n < 1)
	{
		k = 0;
		return z;
	}
	if (n == 1)
	{
		k = 0;
		s[0] = y[0];
		z = y[0];
		return z;
	}
	if (n == 2)
	{
		k = 0;
		s[0] = y[0];
		s[1] = (y[1] - y[0]) / (x[1] - x[0]);
		z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return z;
	}
	if (t <= x[1])
	{
		k = 0;
	}
	else if (t >= x[n - 1])
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
			if (t < x[kk - 1])
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
	u[2] = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);

	if (n == 3)
	{
		if (k == 0)
		{
			u[3] = (y[2] - y[1]) / (x[2] - x[1]);
			u[4] = 2 * u[3] - u[2];
			u[1] = 2 * u[2] - u[3];
			u[0] = 2 * u[1] - u[2];
		}
		else
		{
			u[1] = (y[1] - y[0]) / (x[1] - x[0]);
			u[0] = 2 * u[1] - u[2];
			u[3] = 2 * u[2] - u[1];
			u[4] = 2 * u[3] - u[2];
		}
	}
	else
	{
		if (k <= 1)
		{
			u[3] = (y[k + 2] - y[k + 1]) / (x[k + 2] - x[k + 1]);
			if (k == 1)
			{
				u[1] = (y[1] - y[0]) / (x[1] - x[0]);
				u[0] = 2 * u[1] - u[2];
				if (n == 4)
				{
					u[4] = 2 * u[3] - u[2];
				}
				else
				{
					u[4] = (y[4] - y[3]) / (x[4] - x[3]);
				}
			}
			else if (k >= (n - 3))
			{
				u[1] = (y[k] - y[k - 1]) / (x[k] - x[k - 1]);
				if (k == (n - 3))
				{
					u[3] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
					u[4] = 2 * u[3] - u[2];
					if (n == 4)
					{
						u[0] = 2 * u[1] - u[2];
					}
					else
					{
						u[0] = (y[k - 1] - y[k - 2]) / (x[k - 1] - x[k - 2]);
					}
				}
				else
				{
					u[3] = 2 * u[2] - u[1];
					u[4] = 2 * u[3] - u[2];
					u[0] = (y[k - 1] - y[k - 2]) / (x[k - 1] - x[k - 2]);
				}
			}
			else
			{
				u[1] = (y[k] - y[k - 1]) / (x[k] - x[k - 1]);
				u[0] = (y[k - 1] - y[k - 2]) / (x[k - 1] - x[k - 2]);
				u[3] = (y[k + 2] - y[k + 1]) / (x[k + 2] - x[k + 1]);
				u[4] = (y[k + 3] - y[k + 2]) / (x[k + 3] - x[k + 2]);
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
	s[0] = y[k];
	s[1] = p;
	s[3] = x[k + 1] - x[k];
	s[2] = (3 * u[2] - 2 * p - q) / s[3];
	s[3] = (q + p - 2 * u[2]) / (s[3] * s[3]);
	p = t - x[k];
	z = s[0] + s[1] * p + s[2] * p*p + s[3] * p*p*p;

	return z;
}

double Lagrange2(double *x, double *y, double **z, int n, int m, double u, double v)
{
	int i, j, k, kk;
	int ip, ipp, iq, iqq;
	double h, b[10];
	double f;

	if (u <= x[0])
	{
		ip = 1;
		ipp = 4;
	}
	else if (u >= x[n - 1])
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
			if (u < x[kk - 1])
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

	if (v <= y[0])
	{
		iq = 1;
		iqq = 4;
	}
	else if (v >= y[m - 1])
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
			if (v < y[kk - 1])
			{
				j = kk;
			}
			else
			{
				i = kk;
			}
		}
		iq = 3;
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
			h = z[i][j];
			for (k = iq - 1; k <= iqq - 1; k++)
			{
				if (k != j)
				{
					h = h*(v - y[k]) / (y[j] - y[k]);
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
				h = h*(u - x[j]) / (x[i] - x[j]);
			}
		}
		f = f + h;
	}
	return f;
}

bool FitLS(double *x, double *y, int n, double *a, int m, double *dt)
{
	int i, j, k;
	double z, p, c, g, q, d1, d2;
	double s[20], t[20], b[20];

	for (i = 0; i <= m; i++)
	{
		a[i] = 0;
	}

	if (m + 1 > n)
	{
		m = n - 1;
	}

	if (m > 19)
	{
		m = 19;
	}

	z = 0;

	for (i = 0; i <= n - 1; i++)
	{
		z = z + x[i] / (1.0 * n);
	}

	b[0] = 1.0;
	d1 = 1.0*n;
	p = 0;
	c = 0;

	for (i = 0; i < n - 1; i++)
	{
		p = p + (x[i] - z);
		c = c + y[i];
	}
	c = c / d1;
	p = p / d1;
	a[0] = c*b[0];
	if (m > 0)
	{
		t[1] = 1;
		t[0] = -p;
		d2 = 0;
		c = 0; 
		q = 0;
		for (i = 0; i <= n - 1; i++)
		{
			q = x[i] - z - p;
			d2 = d2 + q*q;
			c = c + y[i] * q;
			g = g + (x[i] - z)*q*q;
		}
		c = c / d2;
		p = g / d2;
		q = d2 / d1;
		d1 = d2;
		a[1] = c*t[1];
		a[0] = c*t[0] + a[0];
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
				q = q*(x[i] - z) + s[k];
			}
			d2 = d2 + q*q;
			c = c + y[i] * q;
			g = g + (x[i] - z)*q*q;
		}

		c = c / d2;
		p = g / d2;
		q = d2 / d1;
		d1 = d2;
		for (k = j - 1; k >= 0; k--)
		{
			a[k] = c*s[k] + a[k];
			b[k] = t[k];
			t[k] = s[k];
		}
	}

	dt[0] = 0;
	dt[1] = 0;
	dt[2] = 0;
	for (i = 0; i <= n - 1; i++)
	{
		q = a[m];
		for (k = m - 1; k >= 0; k--)
		{
			q = a[k] + q*(x[i] - z);
		}
		p = q - y[i];
		if (fabs(p) > dt[2])
		{
			dt[2] = fabs(p);
		}
		dt[0] = dt[0] + p*p;
		dt[1] = dt[1] + fabs(p);
	}

	return 1;
}

bool FitLS2(double *x, double *y, double **z, int n, int m, double **a, int p, int q, double *dt)
{
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