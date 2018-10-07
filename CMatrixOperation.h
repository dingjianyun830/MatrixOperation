#pragma once

#include <iostream>
using namespace std;

class Mat
{
public:
	int m;
	int n;
	int r;
	double **data;
	double eps = 0.000001;
	Mat(int w, int h)
	{
		m = w;
		n = h;
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				data[i][j] = 0;
			}
		}
	}
	//TransDataFromMex
	Mat(double *src, int M, int N)
	{
		m = N;
		n = M;
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				data[i][j] = src[j*M+i];
			}
		}
	}
	void MatPrint();
	// transform Mat to Matlab format
	bool TransDataToMex(Mat A, double *src);
};

Mat Matproduct(Mat src1, Mat src2);
Mat MatDotproduct(Mat src1, double k);
Mat MatDotProduct(Mat src1, Mat src2);
Mat MatInv(Mat src);
Mat MatCopy(Mat src);
double MatDet(Mat src);
double MatRank(Mat src);
bool MatLU(Mat src, Mat dstL, Mat dstU);
bool MatQR(Mat src, Mat dstQ, Mat dstR);
bool MatUaV(Mat src, Mat dstU, Mat dstV);
void ppp(Mat src1, Mat src2, double *s, double*e);
void sss(double fg[2], double cs[2]);

bool MatGauss(Mat A, Mat b, int n);
bool MatGausJordan(Mat A, Mat B, int n, int m);
