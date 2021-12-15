#include <iostream>
#include <set>
#include <algorithm>
#include<stdio.h>
#include<stdlib.h>
#include <windows.h>
#include<time.h>
#include "sheep.h"


double GPu(double x, double a, double k, double m) {
	if (x > a)
		return k * pow((x - a), m);
	else if (x < -a) {
		return k * pow((-x - a), m);
	}
	else
		return 0;
}
double GPy(double x) {
	return 1 + 0.25 * (x + 1);
}
void mul_matrix(double z[DIM], double M[DIM][DIM]) {
	double z2[DIM];
	int i = 0, j = 0;
	for (int i = 0; i < DIM; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < DIM; j++)
	{
		for (i = 0; i < DIM; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}
}
void mul_matrix01(double z[1], double M[1][1]) {
	double z2[1];
	int i = 0, j = 0;
	for (int i = 0; i < 1; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 1; j++)
	{
		for (i = 0; i < 1; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}
}
void mul_matrix02(double z[2], double M[2][2]) {
	double z2[2];
	int i = 0, j = 0;
	for (int i = 0; i < 2; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 2; j++)
	{
		for (i = 0; i < 2; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[3] M[3][3]
void mul_matrix03(double z[3], double M[3][3]) {
	double z2[3];
	int i = 0, j = 0;
	for (int i = 0; i < 3; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 3; j++)
	{
		for (i = 0; i < 3; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
void mul_matrix04(double z[4], double M[4][4]) {
	double z2[4];
	int i = 0, j = 0;
	for (int i = 0; i < 4; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 4; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}
}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[5] M[5][5]
void mul_matrix05(double z[5], double M[5][5]) {
	double z2[5];
	int i = 0, j = 0;
	for (int i = 0; i < 5; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 5; j++)
	{
		for (i = 0; i < 5; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}
}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[6] M[6][6]
void mul_matrix06(double z[6], double M[6][6]) {
	double z2[6];
	int i = 0, j = 0;
	for (int i = 0; i < 6; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 6; j++)
	{
		for (i = 0; i < 6; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
void mul_matrix08(double z[8], double M[8][8]) {
	double z2[8];
	int i = 0, j = 0;
	for (int i = 0; i < 8; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 8; j++)
	{
		for (i = 0; i < 8; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[9] M[9][9]
void mul_matrix09(double z[9], double M[9][9]) {
	double z2[9];
	int i = 0, j = 0;
	for (int i = 0; i < 9; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 9; j++)
	{
		for (i = 0; i < 9; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
void mul_matrix10(double z[10], double M[10][10]) {
	double z2[10];
	int i = 0, j = 0;
	for (int i = 0; i < 10; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 10; j++)
	{
		for (i = 0; i < 10; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[12] M[12][12]
void mul_matrix12(double z[12], double M[12][12]) {
	double z2[12];
	int i = 0, j = 0;
	for (int i = 0; i < 12; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 12; j++)
	{
		for (i = 0; i < 12; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}

void mul_matrix15(double z[15], double M[15][15]) {
	double z2[15];
	int i = 0, j = 0;
	for (int i = 0; i < 15; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 15; j++)
	{
		for (i = 0; i < 15; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[10] M[10][10]
void mul_matrix1(double z[10], double M[10][10]) {
	double z2[10];
	int i = 0, j = 0;
	for (int i = 0; i < 10; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 10; j++)
	{
		for (i = 0; i < 10; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
//MÐý×ª¾ØÕó£¬ÊµÏÖz=M(z) z[20] M[20][20]
void mul_matrix2(double z[20], double M[20][20]) {
	double z2[20];
	int i = 0, j = 0;
	for (int i = 0; i < 20; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 20; j++)
	{
		for (i = 0; i < 20; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
void mul_matrix20(double z[20], double M[20][20]) {
	double z2[20];
	int i = 0, j = 0;
	for (int i = 0; i < 20; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 20; j++)
	{
		for (i = 0; i < 20; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
void mul_matrix30(double z[30], double M[30][30]) {
	double z2[30];
	int i = 0, j = 0;
	for (int i = 0; i < 30; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 30; j++)
	{
		for (i = 0; i < 30; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}
void mul_matrix40(double z[40], double M[40][40]) {
	double z2[40];
	int i = 0, j = 0;
	for (int i = 0; i < 40; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0; j < 40; j++)
	{
		for (i = 0; i < 40; i++)
		{
			z[j] += z2[i] * M[i][j];
		}
	}

}





double Computafitness(double X[]) {

	//CEC2017µÄº¯Êý
	//Unimodal Function
	//shifted and rotated Bent Cigar
#if shifted_rotated_BentCigar
	extern double BC_O[100];
	extern double BC_M10[10][10];
	extern double BC_M20[20][20];
	extern double BC_M30[30][30];
	extern double BC_M50[50][50];
	extern double BC_M100[100][100];

	double z[DIM];
	for (int i = 0; i < DIM; i++)
		z[i] = X[i] - BC_O[i];

#if DIM==10
	mul_matrix(z, BC_M10);
#endif

#if DIM==20
	mul_matrix(z, BC_M20);
#endif

#if DIM==30
	mul_matrix(z, BC_M30);
#endif

#if DIM==50
	mul_matrix(z, BC_M50);
#endif

#if DIM==100
	mul_matrix(z, BC_M100);
#endif

	double f = 0;
	double y1 = z[0] * z[0];
	double y2 = 0.0;
	for (int i = 1; i < DIM; i++)
		y2 += 1E6*z[i] * z[i];

	f = y1 + y2;
	return f + 100;
#endif

#if shifted_rotated_sum_of_Different_Power
	extern double SD_O[100];
	extern double SD_M10[10][10];
	extern double SD_M20[20][20];
	extern double SD_M30[30][30];
	extern double SD_M50[50][50];
	extern double SD_M100[100][100];
	double z[DIM];
	for (int i = 0; i < DIM; i++)
		z[i] = X[i] - SD_O[i];

#if DIM==10
	mul_matrix(z, SD_M10);
#endif

#if DIM==20
	mul_matrix(z, SD_M20);
#endif

#if DIM==30
	mul_matrix(z, SD_M30);
#endif

#if DIM==50
	mul_matrix(z, SD_M50);
#endif

#if DIM==100
	mul_matrix(z, SD_M100);
#endif
	double f = 0.0;
	for (int i = 0; i < DIM; i++)
		f += pow(fabs(z[i]), i + 2);
	return f + 200;
#endif

#if shifted_rotated_Zakharov
	extern double Za_O[100];
	extern double Za_M10[10][10];
	extern double Za_M20[20][20];
	extern double Za_M30[30][30];
	extern double Za_M50[50][50];
	extern double Za_M100[100][100];


	double z[DIM];
	for (int i = 0; i < DIM; i++)
		z[i] = X[i] - Za_O[i];

#if DIM==10
	mul_matrix(z, Za_M10);
#endif

#if DIM==20
	mul_matrix(z, Za_M20);
#endif

#if DIM==30
	mul_matrix(z, Za_M30);
#endif

#if DIM==50
	mul_matrix(z, Za_M50);
#endif

#if DIM==100
	mul_matrix(z, Za_M100);
#endif
	double y1 = 0.0;
	double y2 = 0.0;
	double f = 0;
	for (int i = 0; i < DIM; i++)
	{
		y1 += z[i] * z[i];
		y2 += 0.5*z[i];
	}
	f = y1 + y2 * y2 + y2 * y2*y2*y2;
	return f + 300;
#endif
	//Simple Multimodel Functions
#if shifted_rotated_Rosenbrock
	extern double Ro_O[100];
	extern double Ro_M10[10][10];
	extern double Ro_M20[20][20];
	extern double Ro_M30[30][30];
	extern double Ro_M50[50][50];
	extern double Ro_M100[100][100];
	double z[DIM];
	for (int i = 0; i < DIM; i++)
		z[i] = 2.048*(X[i] - Ro_O[i]) / 100;

#if DIM==10
	mul_matrix(z, Ro_M10);
#endif

#if DIM==20
	mul_matrix(z, Ro_M20);
#endif

#if DIM==30
	mul_matrix(z, Ro_M30);
#endif

#if DIM==50
	mul_matrix(z, Ro_M50);
#endif

#if DIM==100
	mul_matrix(z, Ro_M100);
#endif

	for (int i = 0; i < DIM; i++)
		z[i] = z[i] + 1;
	double result;
	double f1 = 0;
	for (int i = 0; i < DIM - 1; i++)
	{
		f1 += 100 * (z[i] * z[i] - z[i + 1]) * (z[i] * z[i] - z[i + 1]) + (z[i] - 1) * (z[i] - 1);
	}
	result = f1;
	return result + 400;
#endif

#if shifted_rotated_Rastrigin
	double f = 0;
	extern double Ra1_O[100];
	extern double Ra1_M10[10][10];
	extern double Ra1_M20[20][20];
	extern double Ra1_M30[30][30];
	extern double Ra1_M50[50][50];
	extern double Ra1_M100[100][100];

	double z[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z[i] = X[i] - Ra1_O[i];
	}

#if DIM == 10
	mul_matrix(z, Ra1_M10);
#endif
#if DIM == 20
	mul_matrix(z, Ra1_M20);
#endif
#if DIM == 30
	mul_matrix(z, Ra1_M30);
#endif
#if DIM == 50
	mul_matrix(z, Ra1_M50);
#endif
#if DIM == 100
	mul_matrix(z, Ra1_M100);
#endif
	for (int i = 0; i < DIM; i++) {
		f = f + (z[i] * z[i] - 10 * cos(2 * pai*z[i]) + 10);
	}
	return f + 500;
#endif

#if shifted_Rotated_Scaffer7
	extern double Sc_O[100];
	extern double Sc_M10[10][10];
	extern double Sc_M20[20][20];
	extern double Sc_M30[30][30];
	extern double Sc_M50[50][50];
	extern double Sc_M100[100][100];

	double z[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z[i] = 0.005* (X[i] - Sc_O[i]);
	}
#if DIM == 10
	mul_matrix(z, Sc_M10);
#endif
#if DIM == 20
	mul_matrix(z, Sc_M20);
#endif
#if DIM == 30
	mul_matrix(z, Sc_M30);
#endif
#if DIM == 50
	mul_matrix(z, Sc_M50);
#endif
#if DIM == 100
	mul_matrix(z, Sc_M100);
#endif
	double f = 0.0;
	double y = 0;
	for (int i = 0; i < DIM - 1; i++)
	{
		y = z[i] * z[i] + z[i + 1] * z[i + 1];
		y = pow(y, 0.5);
		f += pow(y, 0.5)*(sin(50 * pow(y, 0.2) + 1));
	}
	f = f / (DIM - 1);
	f = pow(f, 2);
	return f + 600;
#endif


#if shifted_Rotated_Lunacek_Bi_Rastrigin

	double f = 0;
	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(DIM + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[DIM];
	double F1 = 0;
	double F2 = 0;
	double temp = 0;
	extern double BR_O[100];
	extern double BR_M10[10][10];
	extern double BR_M20[20][20];
	extern double BR_M30[30][30];
	extern double BR_M50[50][50];
	extern double BR_M100[100][100];
	double v[DIM][DIM];
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (DIM - 1)));
	}
	double X1[DIM];
	for (int i = 0; i < DIM; i++)
		X1[i] = 6 * (X[i] - BR_O[i]);
#if  DIM==10	
	mul_matrix(X1, BR_M10);
#endif

#if DIM==20
	mul_matrix(X1, BR_M20);
#endif

#if  DIM==30	
	mul_matrix(X1, BR_M30);
#endif

#if DIM==50
	mul_matrix(X1, BR_M50);
#endif

#if DIM==100
	mul_matrix(X1, BR_M100);
#endif


	double z[DIM];
	double y[DIM];
	double xz[DIM];//x»¡

	for (int i = 0; i < DIM; i++)
		y[i] = 0.1*(X1[i] - BR_O[i]);

	for (int i = 0; i < DIM; i++)
	{
		if (X1[i] < 0)
			sign[i] = -1;
		else if (X1[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * y[i] + u0;

	}

	for (int i = 0; i < DIM; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < DIM; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * DIM;
	if (F1 < F2)
		f = F1;
	else f = F2;

	for (int i = 0; i < DIM; i++)
		temp = temp + cos(2 * pai*z[i]);
	f = f + 10 * (DIM - temp);


	return f + 700;
#endif

#if shifted_Rotated_Non_Continuous_Rastrigin
	double f = 0;
	double c1 = 0;
	double c2 = 0;
	double sign[DIM];

	double xnew[DIM];//xÃ°
	double xnew1[DIM];//x¼â
	double y[DIM];

	double mid[DIM][DIM];
	double mid1[DIM][DIM];

	double v[DIM][DIM];//¶Ô½Ç¾ØÕó
	double Tasy[DIM];
	double Tosz[DIM];
	extern double NR_O[100];
	extern double NR_M10[10][10];
	extern double NR_M20[20][20];
	extern double NR_M30[30][30];
	extern double NR_M50[50][50];
	extern double NR_M100[100][100];

	double X1[DIM];
	for (int i = 0; i < DIM; i++)
	{
		X1[i] = 0.0512*(X[i] - NR_O[i]);
	}
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
			v[i][j] = 0;
		v[i][i] = pow(10, i / (2 * (DIM - 1)));
	}
	for (int i = 0; i < DIM; i++)
		xnew[i] = 0.0512 * (X1[i] - NR_O[i]);

#if DIM==10
	mul_matrix(xnew, NR_M10);
#endif

#if DIM==20
	mul_matrix(xnew, NR_M20);
#endif

#if DIM==30
	mul_matrix(xnew, NR_M30);
#endif

#if DIM==50
	mul_matrix(xnew, NR_M50);
#endif

#if DIM==100
	mul_matrix(xnew, NR_M100);
#endif

	for (int i = 0; i < DIM; i++)
	{
		if (fabs(xnew[i]) < 0.5)
			y[i] = xnew[i];
		else
			y[i] = round(2 * xnew[i]) / 2;
	}
	for (int i = 0; i < DIM; i++)
	{
		if (y[i] != 0)
		{
			if (y[i] < 0)
			{
				c1 = 5.5;
				c2 = 3.1;
				sign[i] = -1;
			}
			else
			{
				c1 = 10;
				c2 = 7.9;
				sign[i] = 1;
			}
			xnew1[i] = log(fabs(y[i]));
		}
		else
		{
			xnew1[i] = 0;
			sign[i] = 0;
			c1 = 5.5;
			c2 = 3.1;
		}
		Tosz[i] = sign[i] * exp(xnew1[i] + 0.049*(sin(c1*xnew1[i]) + sin(c2*xnew1[i])));
		if (Tosz[i] > 0)
			Tasy[i] = pow(Tosz[i], 1 + 0.2*i / (DIM - 1)*pow(Tosz[i], 0.5));
		else Tasy[i] = Tosz[i];

	}
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
			mid[i][j] = NR_M30[i][j] * v[j][j];

	}
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			double temp = 0;
			for (int k = 0; k < DIM; k++)
				temp += mid[i][k] * NR_M30[k][j];
			mid1[i][j] = temp;
		}
	}
	mul_matrix(Tasy, mid1);
	for (int i = 0; i < DIM; i++)
		f += Tasy[i] * Tasy[i] - 10 * cos(2 * pai*Tasy[i]) + 10;

	return f + 800;
#endif	


#if Shifted_Rotated_Levy
	extern double Le_O[100];
	extern double Le_M10[10][10];
	extern double Le_M20[20][20];
	extern double Le_M30[30][30];
	extern double Le_M50[50][50];
	extern double Le_M100[100][100];
	double z[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z[i] = 5.12*(X[i] - Le_O[i]) / 100;
	}
#if DIM == 10
	mul_matrix(z, Le_M10);
#endif

#if DIM == 20
	mul_matrix(z, Le_M20);
#endif

#if DIM == 30
	mul_matrix(z, Le_M30);
#endif

#if DIM == 50
	mul_matrix(z, Le_M50);
#endif

#if DIM == 100
	mul_matrix(z, Le_M100);
#endif

	double W[DIM];
	double f = 0.0;
	for (int i = 0; i < DIM; i++)
	{
		W[i] = 1 + (z[i] - 1) / 4;
	}

	double f1 = sin(pai*W[0])*sin(pai*W[0]);
	for (int j = 0; j < DIM - 1; j++)
		f += (W[j] - 1)*(W[j] - 1)*(1 + 10 * sin(pai*W[j] + 1)* sin(pai*W[j] + 1));
	double f2 = (W[DIM - 1] - 1) *(W[DIM - 1] - 1) * (1 + sin(2 * pai*W[DIM - 1])*sin(2 * pai*W[DIM - 1]));
	f = f + f1 + f2;
	return f + 900;

#endif


	//Shifted_Rotated_Schwefel
#if Shifted_Rotated_Schwefel
	extern double MSc_O[100];
	double f = 0;
	double f1 = 0;
	double z[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z[i] = 10 * (X[i] - MSc_O[i]) + 4.209687462275036E+2;
		double mod = fmod(fabs(z[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z[i]) <= 500)
			f += z[i] * sin(pow(fabs(z[i]), 0.5));
		if (z[i] > 500)
			f += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z[i] - 500)*(z[i] - 500) / (10000 * DIM);
		if (z[i] < -500)
			f += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z[i] + 500)*(z[i] + 500) / (10000 * DIM);
	}
	f = 418.9829*DIM - f;
	return f + 1000;
#endif


#if hybrid_function_1
	//F(X)=Zakharov+Rosenbrock+Rastrigin	
	double p[3] = { 0.2,0.4,0.4 };
	extern double Hy1_0[100];
	extern double Hy1_M10[10][10];
	extern double Hy1_M20[20][20];
	extern double Hy1_M30[30][30];
	extern double Hy1_M50[50][50];
	extern double Hy1_M100[100][100];
	int Dim1;
	int Dim2;
	int Dim3;

	Dim1 = DIM / 10 * 2;
	Dim2 = DIM / 10 * 4;
	Dim3 = DIM / 10 * 4;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f = 0;

#if DIM==10
	double Hy11_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy11_M2[i][j] = Hy1_M10[i][j];
	}
	double z1[2];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy1_0[i];
	mul_matrix02(z1, Hy11_M2);
	double y1 = 0.0;
	double y2 = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		y1 += z1[i] * z1[i];
		y2 += 0.5*z1[i];
	}
	f1 = y1 + y2 * y2 + y2 * y2 * y2 * y2;


	double Hy12_M4[4][4];
	for (int i = 2; i < 6; i++)
	{
		for (int j = 2; j < 6; j++)
			Hy12_M4[i - 2][j - 2] = Hy1_M10[i][j];
	}

	double z2[4];
	for (int i = 0; i < 4; i++)
		z2[i] = X[i + 2] - Hy1_0[i + 2];
	mul_matrix04(z2, Hy12_M4);
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);


	double Hy13_M4[4][4];
	for (int i = 6; i < 10; i++)
	{
		for (int j = 6; j < 10; j++)
			Hy13_M4[i - 6][j - 6] = Hy1_M10[i][j];
	}

	double z3[4];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 6] - Hy1_0[i + 6];
	mul_matrix04(z3, Hy13_M4);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	f = f1 + f2 + f3;
#endif

#if DIM==20
	double Hy11_M4[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			Hy11_M4[i][j] = Hy1_M20[i][j];
	}
	double z1[4];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy1_0[i];
	mul_matrix04(z1, Hy11_M4);
	double y1 = 0.0;
	double y2 = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		y1 += z1[i] * z1[i];
		y2 += 0.5*z1[i];
	}
	f1 = y1 + y2 * y2 + y2 * y2 * y2 * y2;


	double Hy12_M8[8][8];
	for (int i = 4; i < 12; i++)
	{
		for (int j = 4; j < 12; j++)
			Hy12_M8[i - 4][j - 4] = Hy1_M20[i][j];
	}
	//extern double Ro_M20[20][20];
	double z2[8];
	for (int i = 0; i < 8; i++)
		z2[i] = X[i + 4] - Hy1_0[i + 4];
	mul_matrix08(z2, Hy12_M8);
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);


	double Hy13_M8[8][8];
	for (int i = 12; i < 20; i++)
	{
		for (int j = 12; j < 20; j++)
			Hy13_M8[i - 12][j - 12] = Hy1_M20[i][j];
	}

	double z3[8];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 12] - Hy1_0[i + 12];
	mul_matrix08(z3, Hy13_M8);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	f = f1 + f2 + f3;
#endif

#if DIM==30
	double Hy11_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy11_M6[i][j] = Hy1_M30[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy1_0[i];
	mul_matrix06(z1, Hy11_M6);
	double y1 = 0.0;
	double y2 = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		y1 += z1[i] * z1[i];
		y2 += 0.5*z1[i];
	}
	f1 = y1 + y2 * y2 + y2 * y2 * y2 * y2;


	double Hy12_M12[12][12];
	for (int i = 6; i < 18; i++)
	{
		for (int j = 6; j < 18; j++)
			Hy12_M12[i - 6][j - 6] = Hy1_M30[i][j];
	}

	double z2[12];
	for (int i = 0; i < 12; i++)
		z2[i] = X[i + 6] - Hy1_0[i + 6];
	mul_matrix12(z2, Hy12_M12);
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);


	double Hy13_M12[12][12];
	for (int i = 18; i < 30; i++)
	{
		for (int j = 18; j < 30; j++)
			Hy13_M12[i - 18][j - 18] = Hy1_M30[i][j];
	}

	double z3[12];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 18] - Hy1_0[i + 18];
	mul_matrix12(z3, Hy13_M12);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	f = f1 + f2 + f3;
#endif
#if DIM==50
	double Hy11_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy11_M10[i][j] = Hy1_M50[i][j];
	}
	double z1[10];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy1_0[i];
	mul_matrix10(z1, Hy11_M10);
	double y1 = 0.0;
	double y2 = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		y1 += z1[i] * z1[i];
		y2 += 0.5*z1[i];
	}
	f1 = y1 + y2 * y2 + y2 * y2 * y2 * y2;


	double Hy12_M20[20][20];
	for (int i = 10; i < 30; i++)
	{
		for (int j = 10; j < 30; j++)
			Hy12_M20[i - 10][j - 10] = Hy1_M50[i][j];
	}
	//extern double Ro_M20[20][20];
	double z2[20];
	for (int i = 0; i < 20; i++)
		z2[i] = X[i + 10] - Hy1_0[i + 10];
	mul_matrix20(z2, Hy12_M20);
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);


	double Hy13_M20[20][20];
	for (int i = 30; i < 50; i++)
	{
		for (int j = 30; j < 50; j++)
			Hy13_M20[i - 30][j - 30] = Hy1_M50[i][j];
	}

	double z3[20];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 30] - Hy1_0[i + 30];
	mul_matrix20(z3, Hy13_M20);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	f = f1 + f2 + f3;
#endif

#if DIM==100
	double Hy11_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy11_M20[i][j] = Hy1_M100[i][j];
	}
	double z1[20];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy1_0[i];
	mul_matrix20(z1, Hy11_M20);
	double y1 = 0.0;
	double y2 = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		y1 += z1[i] * z1[i];
		y2 += 0.5*z1[i];
	}
	f1 = y1 + y2 * y2 + y2 * y2 * y2 * y2;


	double Hy12_M40[40][40];
	for (int i = 20; i < 60; i++)
	{
		for (int j = 20; j < 60; j++)
			Hy12_M40[i - 20][j - 20] = Hy1_M100[i][j];
	}
	//extern double Ro_M20[20][20];
	double z2[40];
	for (int i = 0; i < 40; i++)
		z2[i] = X[i + 20] - Hy1_0[i + 20];
	mul_matrix40(z2, Hy12_M40);
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);


	double Hy13_M40[40][40];
	for (int i = 60; i < 100; i++)
	{
		for (int j = 60; j < 100; j++)
			Hy13_M40[i - 60][j - 60] = Hy1_M100[i][j];
	}

	double z3[40];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 60] - Hy1_0[i + 60];
	mul_matrix40(z3, Hy13_M40);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	f = f1 + f2 + f3;
#endif
	return f + 1100;
#endif

#if hybrid_function_2
	//F(X)=High Conditioned Ellipdtic+Modified Schwefel's F10+BentCigar
	double p[3] = { 0.3,0.3,0.4 };
	int Dim1;
	int Dim2;
	int Dim3;
	Dim1 = DIM / 10 * 3;
	Dim2 = DIM / 10 * 3;
	Dim3 = DIM / 10 * 4;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f = 0;

	extern double Hy2_O[100];
	extern double Hy2_M10[10][10];
	extern double Hy2_M20[20][20];
	extern double Hy2_M30[30][30];
	extern double Hy2_M50[50][50];
	extern double Hy2_M100[100][100];

#if DIM==10
	//High Conditioned Elliptic
	double Hy21_M3[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			Hy21_M3[i][j] = Hy2_M10[i][j];
	}
	double z1[3];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy2_O[i];
	mul_matrix03(z1, Hy21_M3);
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	//Modified Schwefel's F10
	double Hy22_M3[3][3];
	for (int i = 3; i < 6; i++)
	{
		for (int j = 3; j < 6; j++)
			Hy22_M3[i - 3][j - 3] = Hy2_M10[i][j];
	}

	double z2[3];
	for (int i = 0; i < 3; i++)
		z2[i] = X[i + 3] - Hy2_O[i + 3];
	mul_matrix03(z2, Hy22_M3);

	for (int i = 0; i < 3; i++)
	{
		z2[i] = z2[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z2[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z2[i]) <= 500)
			f2 += z2[i] * sin(pow(fabs(z2[i]), 0.5));
		if (z2[i] > 500)
			f2 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] - 500)*(z2[i] - 500) / (10000 * Dim2);
		if (z2[i] < -500)
			f2 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] + 500)*(z2[i] + 500) / (10000 * Dim2);
	}
	f2 = 418.9829*Dim2 - f2;


	//Bent Cigar 
	double Hy23_M4[4][4];
	for (int i = 6; i < 10; i++)
	{
		for (int j = 6; j < 10; j++)
			Hy23_M4[i - 6][j - 6] = Hy2_M10[i][j];
	}

	double z3[4];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 6] - Hy2_O[i + 6];
	mul_matrix04(z3, Hy23_M4);
	double y1 = z3[0] * z3[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim3; i++)
		y2 += 1E6*z3[i] * z3[i];
	f3 = y1 + y2;

	f = f1 + f2 + f3;
#endif

#if DIM==20
	//High Conditioned Elliptic
	double Hy21_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy21_M6[i][j] = Hy2_M20[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy2_O[i];
	mul_matrix06(z1, Hy21_M6);
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	//Modified Schwefel's F10
	double Hy22_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy22_M6[i - 6][j - 6] = Hy2_M20[i][j];
	}

	double z2[6];
	for (int i = 0; i < 6; i++)
		z2[i] = X[i + 6] - Hy2_O[i + 6];
	mul_matrix06(z2, Hy22_M6);
	for (int i = 0; i < 6; i++)
	{
		z2[i] = z2[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z2[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z2[i]) <= 500)
			f2 += z2[i] * sin(pow(fabs(z2[i]), 0.5));
		if (z2[i] > 500)
			f2 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] - 500)*(z2[i] - 500) / (10000 * Dim2);
		if (z2[i] < -500)
			f2 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] + 500)*(z2[i] + 500) / (10000 * Dim2);
	}
	f2 = 418.9829*Dim2 - f2;


	//Bent Cigar 
	double Hy23_M8[8][8];
	for (int i = 12; i < 20; i++)
	{
		for (int j = 12; j < 20; j++)
			Hy23_M8[i - 12][j - 12] = Hy2_M20[i][j];
	}

	double z3[8];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 12] - Hy2_O[i + 12];
	mul_matrix08(z3, Hy23_M8);
	double y1 = z3[0] * z3[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim3; i++)
		y2 += 1E6*z3[i] * z3[i];
	f3 = y1 + y2;

	f = f1 + f2 + f3;
#endif

#if DIM==30
	//High Conditioned Elliptic
	double Hy21_M9[9][9];
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
			Hy21_M9[i][j] = Hy2_M30[i][j];
	}
	double z1[9];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy2_O[i];
	mul_matrix09(z1, Hy21_M9);
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	//Modified Schwefel's F10
	double Hy22_M9[9][9];
	for (int i = 9; i < 18; i++)
	{
		for (int j = 9; j < 18; j++)
			Hy22_M9[i - 9][j - 9] = Hy2_M30[i][j];
	}

	double z2[9];
	for (int i = 0; i < 9; i++)
		z2[i] = X[i + 9] - Hy2_O[i + 9];
	mul_matrix09(z2, Hy22_M9);
	for (int i = 0; i < 9; i++)
	{
		z2[i] = z2[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z2[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z2[i]) <= 500)
			f2 += z2[i] * sin(pow(fabs(z2[i]), 0.5));
		if (z2[i] > 500)
			f2 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] - 500)*(z2[i] - 500) / (10000 * Dim2);
		if (z2[i] < -500)
			f2 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] + 500)*(z2[i] + 500) / (10000 * Dim2);
	}
	f2 = 418.9829*Dim2 - f2;


	//Bent Cigar 
	double Hy23_M12[12][12];
	for (int i = 18; i < 30; i++)
	{
		for (int j = 18; j < 30; j++)
			Hy23_M12[i - 18][j - 18] = Hy2_M30[i][j];
	}

	double z3[12];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 18] - Hy2_O[i + 18];
	mul_matrix12(z3, Hy23_M12);
	double y1 = z3[0] * z3[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim3; i++)
		y2 += 1E6*z3[i] * z3[i];
	f3 = y1 + y2;

	f = f1 + f2 + f3;
#endif

#if DIM==50
	//High Conditioned Elliptic
	double Hy21_M15[15][15];
	for (int i = 0; i < 15; i++)
	{
		for (int j = 0; j < 15; j++)
			Hy21_M15[i][j] = Hy2_M50[i][j];
	}
	double z1[15];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy2_O[i];
	mul_matrix15(z1, Hy21_M15);
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	//Modified Schwefel's F10
	double Hy22_M15[15][15];
	for (int i = 15; i < 30; i++)
	{
		for (int j = 15; j < 30; j++)
			Hy22_M15[i - 15][j - 15] = Hy2_M50[i][j];
	}

	double z2[15];
	for (int i = 0; i < 15; i++)
		z2[i] = X[i + 15] - Hy2_O[i + 15];
	mul_matrix15(z2, Hy22_M15);

	for (int i = 0; i < 15; i++)
	{
		z2[i] = z2[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z2[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z2[i]) <= 500)
			f2 += z2[i] * sin(pow(fabs(z2[i]), 0.5));
		if (z2[i] > 500)
			f2 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] - 500)*(z2[i] - 500) / (10000 * Dim2);
		if (z2[i] < -500)
			f2 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] + 500)*(z2[i] + 500) / (10000 * Dim2);
	}
	f2 = 418.9829*Dim2 - f2;


	//Bent Cigar 
	double Hy23_M20[20][20];
	for (int i = 30; i < 50; i++)
	{
		for (int j = 30; j < 50; j++)
			Hy23_M20[i - 30][j - 30] = Hy2_M50[i][j];
	}

	double z3[20];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 30] - Hy2_O[i + 30];
	mul_matrix20(z3, Hy23_M20);
	double y1 = z3[0] * z3[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim3; i++)
		y2 += 1E6*z3[i] * z3[i];
	f3 = y1 + y2;

	f = f1 + f2 + f3;
#endif

#if DIM==100
	//High Conditioned Elliptic
	double Hy21_M30[30][30];
	for (int i = 0; i < 30; i++)
	{
		for (int j = 0; j < 30; j++)
			Hy21_M30[i][j] = Hy2_M100[i][j];
	}
	double z1[30];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy2_O[i];
	mul_matrix30(z1, Hy21_M30);
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	//Modified Schwefel's F10
	double Hy22_M30[30][30];
	for (int i = 30; i < 60; i++)
	{
		for (int j = 30; j < 60; j++)
			Hy22_M30[i - 30][j - 30] = Hy2_M100[i][j];
	}

	double z2[30];
	for (int i = 0; i < 30; i++)
		z2[i] = X[i + 30] - Hy2_O[i + 30];
	mul_matrix30(z2, Hy22_M30);
	for (int i = 0; i < 15; i++)
	{
		z2[i] = z2[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z2[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z2[i]) <= 500)
			f2 += z2[i] * sin(pow(fabs(z2[i]), 0.5));
		if (z2[i] > 500)
			f2 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] - 500)*(z2[i] - 500) / (10000 * Dim2);
		if (z2[i] < -500)
			f2 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z2[i] + 500)*(z2[i] + 500) / (10000 * Dim2);
	}
	f2 = 418.9829*Dim2 - f2;


	//Bent Cigar 
	double Hy23_M40[40][40];
	for (int i = 60; i < 100; i++)
	{
		for (int j = 60; j < 100; j++)
			Hy23_M40[i - 60][j - 60] = Hy2_M100[i][j];
	}

	double z3[40];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 40] - Hy2_O[i + 60];
	mul_matrix40(z3, Hy23_M40);

	double y1 = z3[0] * z3[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim3; i++)
		y2 += 1E6*z3[i] * z3[i];
	f3 = y1 + y2;

	f = f1 + f2 + f3;
#endif
	return f + 1200;
#endif

#if hybrid_function_3
	//F(X)=BentCigar+Rosenbrock+Lunache Bi-Rastrigin
	double p[3] = { 0.3,0.3,0.4 };
	int Dim1;
	int Dim2;
	int Dim3;
	Dim1 = DIM / 10 * 3;
	Dim2 = DIM / 10 * 3;
	Dim3 = DIM / 10 * 4;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f = 0;

	extern double Hy3_O[100];
	extern double Hy3_M10[10][10];
	extern double Hy3_M20[20][20];
	extern double Hy3_M30[30][30];
	extern double Hy3_M50[50][50];
	extern double Hy3_M100[100][100];

#if DIM==10
	double Hy31_M3[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			Hy31_M3[i][j] = Hy3_M10[i][j];
	}
	double z1[3];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy3_O[i];
	mul_matrix03(z1, Hy31_M3);
	//bentcigar
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	double Hy32_M3[3][3];
	for (int i = 3; i < 6; i++)
	{
		for (int j = 3; j < 6; j++)
			Hy32_M3[i - 3][j - 3] = Hy3_M10[i][j];
	}
	double z2[3];
	for (int i = 0; i < 3; i++)
		z2[i] = X[i + 3] - Hy3_O[i + 3];
	mul_matrix03(z2, Hy32_M3);
	//rosenbrock
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);



	double Hy33_M4[4][4];
	for (int i = 6; i < 10; i++)
	{
		for (int j = 6; j < 10; j++)
			Hy33_M4[i - 6][j - 6] = Hy3_M10[i][j];
	}
	double z3[4];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 6] - Hy3_O[i + 6];
	mul_matrix04(z3, Hy33_M4);
	//Lunache Bi-Rastrigin function
	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0 * u0 - d) / s, 0.5);
	double sign[4];
	double F1 = 0;
	double F2 = 0;
	double temp = 0;
	double v[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[4];
	double y[4];
	double xz[4];//x»¡
	for (int i = 0; i < Dim3; i++)
		y[i] = 0.1*(z3[i] - Hy3_O[i + 6]);

	for (int i = 0; i < Dim3; i++)
	{
		if (z3[i] < 0)
			sign[i] = -1;
		else if (z3[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * y[i] + u0;
	}
	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		temp = temp + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - temp);


	f = f1 + f2 + f3;
#endif

#if DIM==20
	double Hy31_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy31_M6[i][j] = Hy3_M20[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy3_O[i];
	mul_matrix06(z1, Hy31_M6);
	//bentcigar
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	double Hy32_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy32_M6[i - 6][j - 6] = Hy3_M20[i][j];
	}
	double z2[6];
	for (int i = 0; i < 6; i++)
		z2[i] = X[i + 6] - Hy3_O[i + 6];
	mul_matrix06(z2, Hy32_M6);
	//rosenbrock
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);



	double Hy33_M8[8][8];
	for (int i = 12; i < 20; i++)
	{
		for (int j = 12; j < 20; j++)
			Hy33_M8[i - 12][j - 12] = Hy3_M20[i][j];
	}
	double z3[8];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 12] - Hy3_O[i + 12];
	mul_matrix08(z3, Hy33_M8);
	//Lunache Bi-Rastrigin function
	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0 * u0 - d) / s, 0.5);
	double sign[8];
	double F1 = 0;
	double F2 = 0;
	double temp = 0;
	double v[8][8];
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[8];
	double y[8];
	double xz[8];//x»¡
	for (int i = 0; i < Dim3; i++)
		y[i] = 0.1*(z3[i] - Hy3_O[i + 12]);

	for (int i = 0; i < Dim3; i++)
	{
		if (z3[i] < 0)
			sign[i] = -1;
		else if (z3[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * y[i] + u0;
	}
	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		temp = temp + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - temp);


	f = f1 + f2 + f3;
#endif

#if DIM==30
	double Hy31_M9[9][9];
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
			Hy31_M9[i][j] = Hy3_M30[i][j];
	}
	double z1[9];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy3_O[i];
	mul_matrix09(z1, Hy31_M9);
	//bentcigar
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	double Hy32_M9[9][9];
	for (int i = 9; i < 18; i++)
	{
		for (int j = 9; j < 18; j++)
			Hy32_M9[i - 9][j - 9] = Hy3_M30[i][j];
	}
	double z2[9];
	for (int i = 0; i < 9; i++)
		z2[i] = X[i + 9] - Hy3_O[i + 9];
	mul_matrix09(z2, Hy32_M9);
	//rosenbrock
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);



	double Hy33_M12[12][12];
	for (int i = 18; i < 30; i++)
	{
		for (int j = 18; j < 30; j++)
			Hy33_M12[i - 18][j - 18] = Hy3_M30[i][j];
	}
	double z3[12];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 18] - Hy3_O[i + 18];
	mul_matrix12(z3, Hy33_M12);
	//Lunache Bi-Rastrigin function
	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0 * u0 - d) / s, 0.5);
	double sign[12];
	double F1 = 0;
	double F2 = 0;
	double temp = 0;
	double v[12][12];
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 12; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[12];
	double y[12];
	double xz[12];//x»¡
	for (int i = 0; i < Dim3; i++)
		y[i] = 0.1*(z3[i] - Hy3_O[i + 18]);

	for (int i = 0; i < Dim3; i++)
	{
		if (z3[i] < 0)
			sign[i] = -1;
		else if (z3[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * y[i] + u0;
	}
	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		temp = temp + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - temp);


	f = f1 + f2 + f3;
#endif


#if DIM==50
	double Hy31_M15[15][15];
	for (int i = 0; i < 15; i++)
	{
		for (int j = 0; j < 15; j++)
			Hy31_M15[i][j] = Hy3_M50[i][j];
	}
	double z1[15];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy3_O[i];
	mul_matrix15(z1, Hy31_M15);
	//bentcigar
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	double Hy32_M15[15][15];
	for (int i = 15; i < 30; i++)
	{
		for (int j = 15; j < 30; j++)
			Hy32_M15[i - 15][j - 15] = Hy3_M50[i][j];
	}
	double z2[15];
	for (int i = 0; i < 15; i++)
		z2[i] = X[i + 15] - Hy3_O[i + 15];
	mul_matrix15(z2, Hy32_M15);
	//rosenbrock
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);



	double Hy33_M20[20][20];
	for (int i = 30; i < 50; i++)
	{
		for (int j = 30; j < 50; j++)
			Hy33_M20[i - 30][j - 30] = Hy3_M50[i][j];
	}
	double z3[20];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 30] - Hy3_O[i + 30];
	mul_matrix20(z3, Hy33_M20);
	//Lunache Bi-Rastrigin function
	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0 * u0 - d) / s, 0.5);
	double sign[20];
	double F1 = 0;
	double F2 = 0;
	double temp = 0;
	double v[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[20];
	double y[20];
	double xz[20];//x»¡
	for (int i = 0; i < Dim3; i++)
		y[i] = 0.1*(z3[i] - Hy3_O[i + 30]);

	for (int i = 0; i < Dim3; i++)
	{
		if (z3[i] < 0)
			sign[i] = -1;
		else if (z3[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * y[i] + u0;
	}
	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		temp = temp + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - temp);


	f = f1 + f2 + f3;
#endif
#if DIM==100
	double Hy31_M30[30][30];
	for (int i = 0; i < 30; i++)
	{
		for (int j = 0; j < 30; j++)
			Hy31_M30[i][j] = Hy3_M100[i][j];
	}
	double z1[30];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy3_O[i];
	mul_matrix30(z1, Hy31_M30);
	//bentcigar
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	double Hy32_M30[30][30];
	for (int i = 30; i < 60; i++)
	{
		for (int j = 30; j < 60; j++)
			Hy32_M30[i - 30][j - 30] = Hy3_M100[i][j];
	}
	double z2[30];
	for (int i = 0; i < 30; i++)
		z2[i] = X[i + 30] - Hy3_O[i + 30];
	mul_matrix30(z2, Hy32_M30);
	//rosenbrock
	for (int i = 0; i < Dim2 - 1; i++)
		f2 += 100 * (z2[i] * z2[i] - z2[i + 1]) * (z2[i] * z2[i] - z2[i + 1]) + (z2[i] - 1) * (z2[i] - 1);



	double Hy33_M40[40][40];
	for (int i = 60; i < 100; i++)
	{
		for (int j = 60; j < 100; j++)
			Hy33_M40[i - 60][j - 60] = Hy3_M100[i][j];
	}
	double z3[40];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 60] - Hy3_O[i + 60];
	mul_matrix40(z3, Hy33_M40);
	//Lunache Bi-Rastrigin function
	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0 * u0 - d) / s, 0.5);
	double sign[40];
	double F1 = 0;
	double F2 = 0;
	double temp = 0;
	double v[40][40];
	for (int i = 0; i < 40; i++)
	{
		for (int j = 0; j < 40; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[40];
	double y[40];
	double xz[40];//x»¡
	for (int i = 0; i < Dim3; i++)
		y[i] = 0.1*(z3[i] - Hy3_O[i + 60]);

	for (int i = 0; i < Dim3; i++)
	{
		if (z3[i] < 0)
			sign[i] = -1;
		else if (z3[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * y[i] + u0;
	}
	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		temp = temp + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - temp);


	f = f1 + f2 + f3;
#endif
	return f + 1300;
#endif

#if hybrid_function_4
	//F(X)=High_Condition_Elliptic+Ackley+Schaffer7+Rastrigin
	double p[4] = { 0.2,0.2,0.2,0.4 };
	int Dim1;
	int Dim2;
	int Dim3;
	int Dim4;
	Dim1 = DIM / 10 * 2;
	Dim2 = DIM / 10 * 2;
	Dim3 = DIM / 10 * 2;
	Dim4 = DIM / 10 * 4;

	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f = 0;

	extern double Hy4_O[100];
	extern double Hy4_M10[10][10];
	extern double Hy4_M20[20][20];
	extern double Hy4_M30[30][30];
	extern double Hy4_M50[50][50];
	extern double Hy4_M100[100][100];

#if DIM==10
	double Hy41_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy41_M2[i][j] = Hy4_M10[i][j];
	}
	double z1[2];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy4_O[i];
	mul_matrix02(z1, Hy41_M2);
	//High Conditioned Elliptic
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	double Hy42_M2[2][2];
	for (int i = 2; i < 4; i++)
	{
		for (int j = 2; j < 4; j++)
			Hy42_M2[i - 2][j - 2] = Hy4_M10[i][j];
	}
	double z2[2];
	for (int i = 0; i < 2; i++)
		z2[i] = X[i + 2] - Hy4_O[i + 2];
	mul_matrix02(z2, Hy42_M2);
	//Ackley
	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / Dim2) * ff1)) - exp((1.0 / Dim2) * ff2) + exp(1) + 20;


	double Hy43_M2[2][2];
	for (int i = 4; i < 6; i++)
	{
		for (int j = 4; j < 6; j++)
			Hy43_M2[i - 4][j - 4] = Hy4_M10[i][j];
	}
	double z3[2];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 4] - Hy4_O[i + 4];
	mul_matrix02(z3, Hy43_M2);


	//Schaffer7
	double fw = 0.0;
	double Si = 0;
	for (int i = 0; i < Dim3 - 1; i++)
	{
		Si = pow(z3[i] * z3[i] + z3[i + 1] * z3[i + 1], 0.5);
		fw += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	fw = fw / Dim3;
	f3 = fw * fw;

	double Hy44_M4[4][4];
	for (int i = 6; i < 10; i++)
	{
		for (int j = 6; j < 10; j++)
			Hy44_M4[i - 6][j - 6] = Hy4_M10[i][j];
	}
	double z4[4];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 6] - Hy4_O[i + 6];
	mul_matrix04(z4, Hy44_M4);
	//Rastrigin
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	f = f1 + f2 + f3 + f4;
#endif


#if DIM==20
	double Hy41_M4[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			Hy41_M4[i][j] = Hy4_M20[i][j];
	}
	double z1[4];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy4_O[i];
	mul_matrix04(z1, Hy41_M4);
	//High Conditioned Elliptic
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	double Hy42_M4[4][4];
	for (int i = 4; i < 8; i++)
	{
		for (int j = 4; j < 8; j++)
			Hy42_M4[i - 4][j - 4] = Hy4_M20[i][j];
	}
	double z2[4];
	for (int i = 0; i < 4; i++)
		z2[i] = X[i + 4] - Hy4_O[i + 4];
	mul_matrix04(z2, Hy42_M4);
	//Ackley
	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / Dim2) * ff1)) - exp((1.0 / Dim2) * ff2) + exp(1) + 20;


	double Hy43_M4[4][4];
	for (int i = 8; i < 12; i++)
	{
		for (int j = 8; j < 12; j++)
			Hy43_M4[i - 8][j - 8] = Hy4_M20[i][j];
	}
	double z3[4];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 8] - Hy4_O[i + 8];
	mul_matrix04(z3, Hy43_M4);


	//Schaffer7
	double fw = 0.0;
	double Si = 0;
	for (int i = 0; i < Dim3 - 1; i++)
	{
		Si = pow(z3[i] * z3[i] + z3[i + 1] * z3[i + 1], 0.5);
		fw += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	fw = fw / Dim3;
	f3 = fw * fw;

	double Hy44_M8[8][8];
	for (int i = 12; i < 20; i++)
	{
		for (int j = 12; j < 20; j++)
			Hy44_M8[i - 12][j - 12] = Hy4_M20[i][j];
	}
	double z4[8];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 12] - Hy4_O[i + 12];
	mul_matrix08(z4, Hy44_M8);
	//Rastrigin
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	f = f1 + f2 + f3 + f4;
#endif

#if DIM==30
	double Hy41_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy41_M6[i][j] = Hy4_M30[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy4_O[i];
	mul_matrix06(z1, Hy41_M6);
	//High Conditioned Elliptic
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	double Hy42_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy42_M6[i - 6][j - 6] = Hy4_M30[i][j];
	}
	double z2[6];
	for (int i = 0; i < 6; i++)
		z2[i] = X[i + 6] - Hy4_O[i + 6];
	mul_matrix06(z2, Hy42_M6);
	//Ackley
	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / Dim2) * ff1)) - exp((1.0 / Dim2) * ff2) + exp(1) + 20;


	double Hy43_M6[6][6];
	for (int i = 12; i < 18; i++)
	{
		for (int j = 12; j < 18; j++)
			Hy43_M6[i - 12][j - 12] = Hy4_M30[i][j];
	}
	double z3[6];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 12] - Hy4_O[i + 12];
	mul_matrix06(z3, Hy43_M6);


	//Schaffer7
	double fw = 0.0;
	double Si = 0;
	for (int i = 0; i < Dim3 - 1; i++)
	{
		Si = pow(z3[i] * z3[i] + z3[i + 1] * z3[i + 1], 0.5);
		fw += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	fw = fw / Dim3;
	f3 = fw * fw;

	double Hy44_M12[12][12];
	for (int i = 18; i < 30; i++)
	{
		for (int j = 18; j < 30; j++)
			Hy44_M12[i - 18][j - 18] = Hy4_M30[i][j];
	}
	double z4[12];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 18] - Hy4_O[i + 18];
	mul_matrix12(z4, Hy44_M12);
	//Rastrigin
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	f = f1 + f2 + f3 + f4;
#endif

#if DIM==50
	double Hy41_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy41_M10[i][j] = Hy4_M50[i][j];
	}
	double z1[10];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy4_O[i];
	mul_matrix10(z1, Hy41_M10);
	//High Conditioned Elliptic
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	double Hy42_M10[10][10];
	for (int i = 10; i < 20; i++)
	{
		for (int j = 10; j < 20; j++)
			Hy42_M10[i - 10][j - 10] = Hy4_M50[i][j];
	}
	double z2[10];
	for (int i = 0; i < 10; i++)
		z2[i] = X[i + 10] - Hy4_O[i + 10];
	mul_matrix10(z2, Hy42_M10);
	//Ackley
	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / Dim2) * ff1)) - exp((1.0 / Dim2) * ff2) + exp(1) + 20;


	double Hy43_M10[10][10];
	for (int i = 20; i < 30; i++)
	{
		for (int j = 20; j < 30; j++)
			Hy43_M10[i - 20][j - 20] = Hy4_M20[i][j];
	}
	double z3[10];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 20] - Hy4_O[i + 20];
	mul_matrix10(z3, Hy43_M10);


	//Schaffer7
	double fw = 0.0;
	double Si = 0;
	for (int i = 0; i < Dim3 - 1; i++)
	{
		Si = pow(z3[i] * z3[i] + z3[i + 1] * z3[i + 1], 0.5);
		fw += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	fw = fw / Dim3;
	f3 = fw * fw;

	double Hy44_M20[20][20];
	for (int i = 30; i < 50; i++)
	{
		for (int j = 30; j < 50; j++)
			Hy44_M20[i - 30][j - 30] = Hy4_M50[i][j];
	}
	double z4[20];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 30] - Hy4_O[i + 30];
	mul_matrix20(z4, Hy44_M20);
	//Rastrigin
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	f = f1 + f2 + f3 + f4;
#endif

#if DIM==100
	double Hy41_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy41_M20[i][j] = Hy4_M100[i][j];
	}
	double z1[20];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy4_O[i];
	mul_matrix20(z1, Hy41_M20);
	//High Conditioned Elliptic
	for (int i = 0; i < Dim1; i++)
		f1 += pow(pow(10, 6), i / (Dim1 - 1))*z1[i] * z1[i];


	double Hy42_M20[20][20];
	for (int i = 20; i < 40; i++)
	{
		for (int j = 20; j < 40; j++)
			Hy42_M20[i - 20][j - 20] = Hy4_M100[i][j];
	}
	double z2[20];
	for (int i = 0; i < 20; i++)
		z2[i] = X[i + 20] - Hy4_O[i + 20];
	mul_matrix20(z2, Hy42_M20);
	//Ackley
	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / Dim2) * ff1)) - exp((1.0 / Dim2) * ff2) + exp(1) + 20;


	double Hy43_M20[20][20];
	for (int i = 40; i < 60; i++)
	{
		for (int j = 40; j < 60; j++)
			Hy43_M20[i - 40][j - 40] = Hy4_M20[i][j];
	}
	double z3[20];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 40] - Hy4_O[i + 40];
	mul_matrix20(z3, Hy43_M20);


	//Schaffer7
	double fw = 0.0;
	double Si = 0;
	for (int i = 0; i < Dim3 - 1; i++)
	{
		Si = pow(z3[i] * z3[i] + z3[i + 1] * z3[i + 1], 0.5);
		fw += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	fw = fw / Dim3;
	f3 = fw * fw;

	double Hy44_M40[40][40];
	for (int i = 60; i < 100; i++)
	{
		for (int j = 60; j < 100; j++)
			Hy44_M40[i - 60][j - 60] = Hy4_M100[i][j];
	}
	double z4[40];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 60] - Hy4_O[i + 60];
	mul_matrix40(z4, Hy44_M40);
	//Rastrigin
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	f = f1 + f2 + f3 + f4;
#endif
	return f + 1400;
#endif

#if hybrid_function_5
	//F(X)=BentCigar+HGBat+Rastrigin+Rosenbrock
	double p[4] = { 0.2,0.2,0.3,0.3 };
	int Dim1;
	int Dim2;
	int Dim3;
	int Dim4;
	Dim1 = DIM / 10 * 2;
	Dim2 = DIM / 10 * 2;
	Dim3 = DIM / 10 * 3;
	Dim4 = DIM / 10 * 3;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f = 0;

	extern double Hy5_O[100];
	extern double Hy5_M10[10][10];
	extern double Hy5_M20[20][20];
	extern double Hy5_M30[30][30];
	extern double Hy5_M50[50][50];
	extern double Hy5_M100[100][100];

#if DIM==10
	//Bent Cigar function
	double Hy51_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy51_M2[i][j] = Hy5_M10[i][j];
	}
	double z1[2];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy5_O[i];
	mul_matrix02(z1, Hy51_M2);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	//HGBat function
	double Hy52_M2[2][2];
	for (int i = 2; i < 4; i++)
	{
		for (int j = 2; j < 4; j++)
			Hy52_M2[i - 2][j - 2] = Hy5_M10[i][j];
	}

	double z2[2];
	for (int i = 0; i < 2; i++)
		z2[i] = X[i + 2] - Hy5_O[i + 2];
	mul_matrix02(z2, Hy52_M2);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	f2 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim2 + 0.5;


	//Rastrigin function 
	double Hy53_M3[3][3];
	for (int i = 4; i < 7; i++)
	{
		for (int j = 4; j < 7; j++)
			Hy53_M3[i - 4][j - 4] = Hy5_M10[i][j];
	}

	double z3[3];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 4] - Hy5_O[i + 4];
	mul_matrix03(z3, Hy53_M3);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M3[3][3];
	for (int i = 7; i < 10; i++)
	{
		for (int j = 7; j < 10; j++)
			Hy54_M3[i - 7][j - 7] = Hy5_M10[i][j];
	}

	double z4[3];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 7] - Hy5_O[i + 7];
	mul_matrix03(z4, Hy54_M3);

	for (int i = 0; i < Dim4 - 1; i++)
		f4 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);


	f = f1 + f2 + f3 + f4;
#endif

#if DIM==20
	//Bent Cigar function
	double Hy51_M4[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			Hy51_M4[i][j] = Hy5_M20[i][j];
	}
	double z1[4];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy5_O[i];
	mul_matrix04(z1, Hy51_M4);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	//HGBat function
	double Hy52_M4[4][4];
	for (int i = 4; i < 8; i++)
	{
		for (int j = 4; j < 8; j++)
			Hy52_M4[i - 4][j - 4] = Hy5_M20[i][j];
	}

	double z2[4];
	for (int i = 0; i < 4; i++)
		z2[i] = X[i + 4] - Hy5_O[i + 4];
	mul_matrix04(z2, Hy52_M4);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	f2 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim2 + 0.5;


	//Rastrigin function 
	double Hy53_M6[6][6];
	for (int i = 8; i < 14; i++)
	{
		for (int j = 8; j < 14; j++)
			Hy53_M6[i - 8][j - 8] = Hy5_M20[i][j];
	}

	double z3[6];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 8] - Hy5_O[i + 8];
	mul_matrix06(z3, Hy53_M6);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M6[6][6];
	for (int i = 14; i < 20; i++)
	{
		for (int j = 14; j < 20; j++)
			Hy54_M6[i - 14][j - 14] = Hy5_M20[i][j];
	}

	double z4[6];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 14] - Hy5_O[i + 14];
	mul_matrix06(z4, Hy54_M6);

	for (int i = 0; i < Dim4 - 1; i++)
		f4 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);


	f = f1 + f2 + f3 + f4;
#endif

#if DIM==30
	//Bent Cigar function
	double Hy51_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy51_M6[i][j] = Hy5_M30[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy5_O[i];
	mul_matrix06(z1, Hy51_M6);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	//HGBat function
	double Hy52_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy52_M6[i - 6][j - 6] = Hy5_M30[i][j];
	}

	double z2[6];
	for (int i = 0; i < 6; i++)
		z2[i] = X[i + 6] - Hy5_O[i + 6];
	mul_matrix06(z2, Hy52_M6);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	f2 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim2 + 0.5;


	//Rastrigin function 
	double Hy53_M9[9][9];
	for (int i = 12; i < 21; i++)
	{
		for (int j = 12; j < 21; j++)
			Hy53_M9[i - 12][j - 12] = Hy5_M30[i][j];
	}

	double z3[9];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 12] - Hy5_O[i + 12];
	mul_matrix09(z3, Hy53_M9);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M9[9][9];
	for (int i = 21; i < 30; i++)
	{
		for (int j = 21; j < 30; j++)
			Hy54_M9[i - 21][j - 21] = Hy5_M30[i][j];
	}

	double z4[9];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 21] - Hy5_O[i + 21];
	mul_matrix09(z4, Hy54_M9);

	for (int i = 0; i < Dim4 - 1; i++)
		f4 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);


	f = f1 + f2 + f3 + f4;
#endif

#if DIM==50
	//Bent Cigar function
	double Hy51_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy51_M10[i][j] = Hy5_M100[i][j];
	}
	double z1[10];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy5_O[i];
	mul_matrix10(z1, Hy51_M10);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	//HGBat function
	double Hy52_M10[10][10];
	for (int i = 10; i < 20; i++)
	{
		for (int j = 10; j < 20; j++)
			Hy52_M10[i - 10][j - 10] = Hy5_M100[i][j];
	}

	double z2[10];
	for (int i = 0; i < 10; i++)
		z2[i] = X[i + 10] - Hy5_O[i + 10];
	mul_matrix10(z2, Hy52_M10);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	f2 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim2 + 0.5;


	//Rastrigin function 
	double Hy53_M15[15][15];
	for (int i = 20; i < 35; i++)
	{
		for (int j = 20; j < 35; j++)
			Hy53_M15[i - 20][j - 20] = Hy5_M50[i][j];
	}

	double z3[15];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 20] - Hy5_O[i + 20];
	mul_matrix15(z3, Hy53_M15);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M15[15][15];
	for (int i = 35; i < 50; i++)
	{
		for (int j = 35; j < 50; j++)
			Hy54_M15[i - 35][j - 35] = Hy5_M50[i][j];
	}

	double z4[15];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 35] - Hy5_O[i + 35];
	mul_matrix15(z4, Hy54_M15);

	for (int i = 0; i < Dim4 - 1; i++)
		f4 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);


	f = f1 + f2 + f3 + f4;
#endif

#if DIM==100
	//Bent Cigar function
	double Hy51_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy51_M20[i][j] = Hy5_M100[i][j];
	}
	double z1[20];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy5_O[i];
	mul_matrix20(z1, Hy51_M20);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	//HGBat function
	double Hy52_M20[20][20];
	for (int i = 20; i < 40; i++)
	{
		for (int j = 20; j < 40; j++)
			Hy52_M20[i - 20][j - 20] = Hy5_M100[i][j];
	}

	double z2[20];
	for (int i = 0; i < 20; i++)
		z2[i] = X[i + 20] - Hy5_O[i + 20];
	mul_matrix20(z2, Hy52_M20);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	f2 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim2 + 0.5;


	//Rastrigin function 
	double Hy53_M30[30][30];
	for (int i = 40; i < 70; i++)
	{
		for (int j = 40; j < 70; j++)
			Hy53_M30[i - 40][j - 40] = Hy5_M100[i][j];
	}

	double z3[30];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 40] - Hy5_O[i + 40];
	mul_matrix30(z3, Hy53_M30);
	for (int i = 0; i < Dim3; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M30[30][30];
	for (int i = 70; i < 100; i++)
	{
		for (int j = 70; j < 100; j++)
			Hy54_M30[i - 70][j - 70] = Hy5_M100[i][j];
	}

	double z4[30];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 70] - Hy5_O[i + 70];
	mul_matrix30(z4, Hy54_M30);

	for (int i = 0; i < Dim4 - 1; i++)
		f4 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);


	f = f1 + f2 + f3 + f4;
#endif
	return f + 1500;

#endif

#if hybrid_function_6
	//F(X)=Schaffer f6+HGBat+Rosenbrock+Modified Schwefel's f10
	double p[4] = { 0.2,0.2,0.3,0.3 };
	int Dim1;
	int Dim2;
	int Dim3;
	int Dim4;

	Dim1 = DIM / 10 * 2;
	Dim2 = DIM / 10 * 2;
	Dim3 = DIM / 10 * 3;
	Dim4 = DIM / 10 * 3;

	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f = 0;

	extern double Hy6_O[100];
	extern double Hy6_M10[10][10];
	extern double Hy6_M20[20][20];
	extern double Hy6_M30[30][30];
	extern double Hy6_M50[50][50];
	extern double Hy6_M100[100][100];

#if DIM==10
	//Expanded Schaffer F6
	double Hy61_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy61_M2[i][j] = Hy6_M10[i][j];
	}
	double z1[2];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy6_O[i];
	mul_matrix02(z1, Hy61_M2);
	double g = 0.0;
	for (int i = 0; i < Dim1 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z1[i] * z1[i] + z1[i + 1] * z1[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[i] * z1[i] + z1[i + 1] * z1[i + 1])), 2);
		f1 += g;
	}
	f1 += 0.5 + (pow(sin(pow(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0])), 2);



	//HGBat
	double Hy62_M2[2][2];
	for (int i = 2; i < 4; i++)
	{
		for (int j = 2; j < 4; j++)
			Hy62_M2[i - 2][j - 2] = Hy6_M10[i][j];
	}

	double z2[2];
	for (int i = 0; i < 2; i++)
		z2[i] = X[i + 2] - Hy6_O[i + 2];
	mul_matrix02(z2, Hy62_M2);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += z2[i];
	}
	f2 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim2 + 0.5;



	//Rosenbrock
	double Hy63_M3[3][3];
	for (int i = 4; i < 7; i++)
	{
		for (int j = 4; j < 7; j++)
			Hy63_M3[i - 4][j - 4] = Hy6_M10[i][j];
	}

	double z3[3];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 4] - Hy6_O[i + 4];
	mul_matrix03(z3, Hy63_M3);
	for (int i = 0; i < Dim3 - 1; i++)
		f3 += 100 * (z3[i] * z3[i] - z3[i + 1]) * (z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1) * (z3[i] - 1);


	//Modified Schwefel F10
	double Hy64_M3[3][3];
	for (int i = 7; i < 10; i++)
	{
		for (int j = 7; j < 10; j++)
			Hy64_M3[i - 7][j - 7] = Hy6_M10[i][j];
	}

	double z4[3];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 7] - Hy6_O[i + 7];
	mul_matrix03(z4, Hy64_M3);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;

	f = f1 + f2 + f3 + f4;
#endif

#if DIM==20
	//Expanded Schaffer F6
	double Hy61_M4[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			Hy61_M4[i][j] = Hy6_M20[i][j];
	}
	double z1[4];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy6_O[i];
	mul_matrix04(z1, Hy61_M4);
	double g = 0.0;
	for (int i = 0; i < Dim1 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z1[i] * z1[i] + z1[i + 1] * z1[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[i] * z1[i] + z1[i + 1] * z1[i + 1])), 2);
		f1 += g;
	}
	f1 += 0.5 + (pow(sin(pow(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0])), 2);



	//HGBat
	double Hy62_M4[4][4];
	for (int i = 4; i < 8; i++)
	{
		for (int j = 4; j < 8; j++)
			Hy62_M4[i - 4][j - 4] = Hy6_M20[i][j];
	}

	double z2[4];
	for (int i = 0; i < 4; i++)
		z2[i] = X[i + 4] - Hy6_O[i + 4];
	mul_matrix04(z2, Hy62_M4);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += z2[i];
	}
	f2 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim2 + 0.5;



	//Rosenbrock
	double Hy63_M6[6][6];
	for (int i = 8; i < 14; i++)
	{
		for (int j = 8; j < 14; j++)
			Hy63_M6[i - 8][j - 8] = Hy6_M20[i][j];
	}

	double z3[6];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 8] - Hy6_O[i + 8];
	mul_matrix06(z3, Hy63_M6);
	for (int i = 0; i < Dim3 - 1; i++)
		f3 += 100 * (z3[i] * z3[i] - z3[i + 1]) * (z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1) * (z3[i] - 1);


	//Modified Schwefel F10
	double Hy64_M6[6][6];
	for (int i = 14; i < 20; i++)
	{
		for (int j = 14; j < 20; j++)
			Hy64_M6[i - 14][j - 14] = Hy6_M20[i][j];
	}

	double z4[6];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 14] - Hy6_O[i + 14];
	mul_matrix06(z4, Hy64_M6);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;

	f = f1 + f2 + f3 + f4;
#endif

#if DIM==30
	//Expanded Schaffer F6
	double Hy61_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy61_M6[i][j] = Hy6_M30[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy6_O[i];
	mul_matrix06(z1, Hy61_M6);
	double g = 0.0;
	for (int i = 0; i < Dim1 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z1[i] * z1[i] + z1[i + 1] * z1[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[i] * z1[i] + z1[i + 1] * z1[i + 1])), 2);
		f1 += g;
	}
	f1 += 0.5 + (pow(sin(pow(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0])), 2);



	//HGBat
	double Hy62_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy62_M6[i - 6][j - 6] = Hy6_M30[i][j];
	}

	double z2[6];
	for (int i = 0; i < 6; i++)
		z2[i] = X[i + 6] - Hy6_O[i + 6];
	mul_matrix06(z2, Hy62_M6);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += z2[i];
	}
	f2 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim2 + 0.5;



	//Rosenbrock
	double Hy63_M9[9][9];
	for (int i = 12; i < 21; i++)
	{
		for (int j = 12; j < 21; j++)
			Hy63_M9[i - 12][j - 12] = Hy6_M30[i][j];
	}

	double z3[9];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 12] - Hy6_O[i + 12];
	mul_matrix09(z3, Hy63_M9);
	for (int i = 0; i < Dim3 - 1; i++)
		f3 += 100 * (z3[i] * z3[i] - z3[i + 1]) * (z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1) * (z3[i] - 1);


	//Modified Schwefel F10
	double Hy64_M9[9][9];
	for (int i = 21; i < 30; i++)
	{
		for (int j = 21; j < 30; j++)
			Hy64_M9[i - 21][j - 21] = Hy6_M30[i][j];
	}

	double z4[9];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 21] - Hy6_O[i + 21];
	mul_matrix09(z4, Hy64_M9);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;

	f = f1 + f2 + f3 + f4;
#endif


#if DIM==50
	//Expanded Schaffer F6
	double Hy61_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy61_M10[i][j] = Hy6_M50[i][j];
	}
	double z1[10];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy6_O[i];
	mul_matrix10(z1, Hy61_M10);
	double g = 0.0;
	for (int i = 0; i < Dim1 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z1[i] * z1[i] + z1[i + 1] * z1[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[i] * z1[i] + z1[i + 1] * z1[i + 1])), 2);
		f1 += g;
	}
	f1 += 0.5 + (pow(sin(pow(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0])), 2);



	//HGBat
	double Hy62_M10[10][10];
	for (int i = 10; i < 20; i++)
	{
		for (int j = 10; j < 20; j++)
			Hy62_M10[i - 10][j - 10] = Hy6_M50[i][j];
	}

	double z2[10];
	for (int i = 0; i < 10; i++)
		z2[i] = X[i + 10] - Hy6_O[i + 10];
	mul_matrix10(z2, Hy62_M10);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += z2[i];
	}
	f2 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim2 + 0.5;



	//Rosenbrock
	double Hy63_M15[15][15];
	for (int i = 20; i < 35; i++)
	{
		for (int j = 20; j < 35; j++)
			Hy63_M15[i - 20][j - 20] = Hy6_M50[i][j];
	}

	double z3[15];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 20] - Hy6_O[i + 20];
	mul_matrix15(z3, Hy63_M15);
	for (int i = 0; i < Dim3 - 1; i++)
		f3 += 100 * (z3[i] * z3[i] - z3[i + 1]) * (z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1) * (z3[i] - 1);


	//Modified Schwefel F10
	double Hy64_M15[15][15];
	for (int i = 35; i < 50; i++)
	{
		for (int j = 35; j < 50; j++)
			Hy64_M15[i - 35][j - 35] = Hy6_M50[i][j];
	}

	double z4[15];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 35] - Hy6_O[i + 35];
	mul_matrix15(z4, Hy64_M15);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;

	f = f1 + f2 + f3 + f4;
#endif	

#if DIM==100
	//Expanded Schaffer F6
	double Hy61_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy61_M20[i][j] = Hy6_M100[i][j];
	}
	double z1[20];
	for (int i = 0; i < Dim1; i++)
		z1[i] = X[i] - Hy6_O[i];
	mul_matrix20(z1, Hy61_M20);
	double g = 0.0;
	for (int i = 0; i < Dim1 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z1[i] * z1[i] + z1[i + 1] * z1[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[i] * z1[i] + z1[i + 1] * z1[i + 1])), 2);
		f1 += g;
	}
	f1 += 0.5 + (pow(sin(pow(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z1[Dim1 - 1] * z1[Dim1 - 1] + z1[0] * z1[0])), 2);



	//HGBat
	double Hy62_M20[20][20];
	for (int i = 20; i < 40; i++)
	{
		for (int j = 20; j < 40; j++)
			Hy62_M20[i - 20][j - 20] = Hy6_M100[i][j];
	}

	double z2[20];
	for (int i = 0; i < 20; i++)
		z2[i] = X[i + 20] - Hy6_O[i + 20];
	mul_matrix20(z2, Hy62_M20);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim2; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += z2[i];
	}
	f2 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim2 + 0.5;



	//Rosenbrock
	double Hy63_M30[30][30];
	for (int i = 40; i < 70; i++)
	{
		for (int j = 40; j < 70; j++)
			Hy63_M30[i - 40][j - 40] = Hy6_M100[i][j];
	}

	double z3[30];
	for (int i = 0; i < Dim3; i++)
		z3[i] = X[i + 40] - Hy6_O[i + 40];
	mul_matrix30(z3, Hy63_M30);
	for (int i = 0; i < Dim3 - 1; i++)
		f3 += 100 * (z3[i] * z3[i] - z3[i + 1]) * (z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1) * (z3[i] - 1);


	//Modified Schwefel F10
	double Hy64_M30[30][30];
	for (int i = 70; i < 100; i++)
	{
		for (int j = 70; j < 100; j++)
			Hy64_M30[i - 70][j - 70] = Hy6_M100[i][j];
	}

	double z4[30];
	for (int i = 0; i < Dim4; i++)
		z4[i] = X[i + 70] - Hy6_O[i + 70];
	mul_matrix30(z4, Hy64_M30);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;

	f = f1 + f2 + f3 + f4;
#endif

	return f + 1600;

#endif

#if hybrid_function_7
	//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
	double p[5] = { 0.1,0.2,0.2,0.2,0.3 };
	int Dim1 = DIM / 10 * 1;
	int Dim2 = DIM / 10 * 2;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 3;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;
	double f = 0;

	extern double Hy7_O[100];
	extern double Hy7_M10[10][10];
	extern double Hy7_M20[20][20];
	extern double Hy7_M30[30][30];
	extern double Hy7_M50[50][50];
	extern double Hy7_M100[100][100];


#if  DIM==10
	double Hy71_M1[1][1];
	for (int i = 0; i < 1; i++)
		for (int j = 0; j < 1; j++)
			Hy71_M1[i][j] = Hy7_M30[i][j];
	double Hy72_M2[2][2];
	for (int i = 1; i < 3; i++)
		for (int j = 1; j < 3; j++)
			Hy72_M2[i - 1][j - 1] = Hy7_M10[i][j];
	double Hy73_M2[2][2];
	for (int i = 3; i < 5; i++)
		for (int j = 3; j < 5; j++)
			Hy73_M2[i - 3][j - 3] = Hy7_M10[i][j];
	double Hy74_M2[2][2];
	for (int i = 5; i < 7; i++)
		for (int j = 5; j < 7; j++)
			Hy74_M2[i - 5][j - 5] = Hy7_M10[i][j];
	double Hy75_M3[3][3];
	for (int i = 7; i < 10; i++)
		for (int j = 7; j < 10; j++)
			Hy75_M3[i - 7][j - 7] = Hy7_M10[i][j];

	double z1[1];
	for (int i = 0; i < 1; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix01(z1, Hy71_M1);

	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim1*Dim1);
	y3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f1 = y1 * y - y1;


	double z2[2];
	for (int i = 1; i < 3; i++)
		z2[i - 1] = X[i] - Hy7_O[i];
	mul_matrix02(z2, Hy72_M2);

	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt(ff1 / Dim2)) - exp(ff2 / Dim2) + exp(1) + 20;

	double z3[2];
	for (int i = 3; i < 5; i++)
		z3[i - 3] = X[i] - Hy7_O[i];
	mul_matrix02(z3, Hy73_M2);
	//Expanded Griewank's plus Rosenbrock
	double temp[2];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[2];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[2];
	double yy[2];
	double xz[2];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 3]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[2];
	for (int i = 5; i < 7; i++)
		z4[i - 5] = X[i] - Hy7_O[i];
	mul_matrix02(z4, Hy74_M2);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;


	double z5[3];
	for (int i = 7; i < 10; i++)
		z5[i - 7] = X[i] - Hy7_O[i];
	mul_matrix03(z5, Hy75_M3);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f5 = f5 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f = f1 + f2 + f3 + f4 + f5;
#endif
#if  DIM==20
	double Hy71_M2[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Hy71_M2[i][j] = Hy7_M20[i][j];
	double Hy72_M4[4][4];
	for (int i = 2; i < 6; i++)
		for (int j = 2; j < 6; j++)
			Hy72_M4[i - 2][j - 2] = Hy7_M20[i][j];
	double Hy73_M4[4][4];
	for (int i = 6; i < 10; i++)
		for (int j = 6; j < 10; j++)
			Hy73_M4[i - 6][j - 6] = Hy7_M20[i][j];
	double Hy74_M4[4][4];
	for (int i = 10; i < 14; i++)
		for (int j = 10; j < 14; j++)
			Hy74_M4[i - 10][j - 10] = Hy7_M20[i][j];
	double Hy75_M6[6][6];
	for (int i = 14; i < 20; i++)
		for (int j = 14; j < 20; j++)
			Hy75_M6[i - 14][j - 14] = Hy7_M20[i][j];

	double z1[2];
	for (int i = 0; i < 2; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix02(z1, Hy71_M2);

	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim1*Dim1);
	y3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f1 = y1 * y - y1;


	double z2[4];
	for (int i = 2; i < 6; i++)
		z2[i - 2] = X[i] - Hy7_O[i];
	mul_matrix04(z2, Hy72_M4);

	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt(ff1 / Dim2)) - exp(ff2 / Dim2) + exp(1) + 20;

	double z3[4];
	for (int i = 6; i < 10; i++)
		z3[i - 6] = X[i] - Hy7_O[i];
	mul_matrix04(z3, Hy73_M4);
	//Expanded Griewank's plus Rosenbrock
	double temp[4];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[4];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[4];
	double yy[4];
	double xz[4];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 6]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[4];
	for (int i = 10; i < 14; i++)
		z4[i - 10] = X[i] - Hy7_O[i];
	mul_matrix04(z4, Hy74_M4);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;


	double z5[6];
	for (int i = 14; i < 20; i++)
		z5[i - 14] = X[i] - Hy7_O[i];
	mul_matrix06(z5, Hy75_M6);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f5 = f5 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f = f1 + f2 + f3 + f4 + f5;
#endif


#if  DIM==30
	double Hy71_M3[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Hy71_M3[i][j] = Hy7_M30[i][j];
	double Hy72_M6[6][6];
	for (int i = 3; i < 9; i++)
		for (int j = 3; j < 9; j++)
			Hy72_M6[i - 3][j - 3] = Hy7_M30[i][j];
	double Hy73_M6[6][6];
	for (int i = 9; i < 15; i++)
		for (int j = 9; j < 15; j++)
			Hy73_M6[i - 9][j - 9] = Hy7_M30[i][j];
	double Hy74_M6[6][6];
	for (int i = 15; i < 21; i++)
		for (int j = 15; j < 21; j++)
			Hy74_M6[i - 15][j - 15] = Hy7_M30[i][j];
	double Hy75_M9[9][9];
	for (int i = 21; i < 30; i++)
		for (int j = 21; j < 30; j++)
			Hy75_M9[i - 21][j - 21] = Hy7_M30[i][j];

	double z1[3];
	for (int i = 0; i < 3; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix03(z1, Hy71_M3);

	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim1*Dim1);
	y3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f1 = y1 * y - y1;


	double z2[6];
	for (int i = 3; i < 9; i++)
		z2[i - 3] = X[i] - Hy7_O[i];
	mul_matrix06(z2, Hy72_M6);

	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt(ff1 / Dim2)) - exp(ff2 / Dim2) + exp(1) + 20;

	double z3[6];
	for (int i = 9; i < 15; i++)
		z3[i - 9] = X[i] - Hy7_O[i];
	mul_matrix06(z3, Hy73_M6);
	//Expanded Griewank's plus Rosenbrock
	double temp[6];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[6];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[6];
	double yy[6];
	double xz[6];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 9]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[6];
	for (int i = 15; i < 21; i++)
		z4[i - 15] = X[i] - Hy7_O[i];
	mul_matrix06(z4, Hy74_M6);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;


	double z5[9];
	for (int i = 21; i < 30; i++)
		z5[i - 21] = X[i] - Hy7_O[i];
	mul_matrix09(z5, Hy75_M9);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f5 = f5 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f = f1 + f2 + f3 + f4 + f5;
#endif
#if  DIM==50
	double Hy71_M5[5][5];
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			Hy71_M5[i][j] = Hy7_M50[i][j];
	double Hy72_M10[10][10];
	for (int i = 5; i < 15; i++)
		for (int j = 5; j < 15; j++)
			Hy72_M10[i - 5][j - 5] = Hy7_M50[i][j];
	double Hy73_M10[10][10];
	for (int i = 15; i < 25; i++)
		for (int j = 15; j < 25; j++)
			Hy73_M10[i - 15][j - 15] = Hy7_M50[i][j];
	double Hy74_M10[10][10];
	for (int i = 25; i < 35; i++)
		for (int j = 25; j < 35; j++)
			Hy74_M10[i - 25][j - 25] = Hy7_M50[i][j];
	double Hy75_M15[15][15];
	for (int i = 35; i < 50; i++)
		for (int j = 35; j < 50; j++)
			Hy75_M15[i - 35][j - 35] = Hy7_M50[i][j];

	double z1[5];
	for (int i = 0; i < 5; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix05(z1, Hy71_M5);

	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim1*Dim1);
	y3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f1 = y1 * y - y1;


	double z2[10];
	for (int i = 5; i < 15; i++)
		z2[i - 5] = X[i] - Hy7_O[i];
	mul_matrix10(z2, Hy72_M10);

	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt(ff1 / Dim2)) - exp(ff2 / Dim2) + exp(1) + 20;

	double z3[10];
	for (int i = 15; i < 25; i++)
		z3[i - 15] = X[i] - Hy7_O[i];
	mul_matrix10(z3, Hy73_M10);
	//Expanded Griewank's plus Rosenbrock
	double temp[10];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[10];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[10];
	double yy[10];
	double xz[10];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 15]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[10];
	for (int i = 25; i < 35; i++)
		z4[i - 25] = X[i] - Hy7_O[i];
	mul_matrix10(z4, Hy74_M10);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;


	double z5[15];
	for (int i = 35; i < 50; i++)
		z5[i - 35] = X[i] - Hy7_O[i];
	mul_matrix15(z5, Hy75_M15);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f5 = f5 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f = f1 + f2 + f3 + f4 + f5;
#endif
#if  DIM==100
	double Hy71_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy71_M10[i][j] = Hy7_M100[i][j];
	double Hy72_M20[20][20];
	for (int i = 10; i < 30; i++)
		for (int j = 10; j < 30; j++)
			Hy72_M20[i - 10][j - 10] = Hy7_M100[i][j];
	double Hy73_M20[20][20];
	for (int i = 30; i < 50; i++)
		for (int j = 30; j < 50; j++)
			Hy73_M20[i - 30][j - 30] = Hy7_M100[i][j];
	double Hy74_M20[20][20];
	for (int i = 50; i < 70; i++)
		for (int j = 50; j < 70; j++)
			Hy74_M20[i - 50][j - 50] = Hy7_M100[i][j];
	double Hy75_M30[30][30];
	for (int i = 70; i < 100; i++)
		for (int j = 70; j < 100; j++)
			Hy75_M30[i - 70][j - 70] = Hy7_M100[i][j];

	double z1[10];
	for (int i = 0; i < 10; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix10(z1, Hy71_M10);

	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim1*Dim1);
	y3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f1 = y1 * y - y1;


	double z2[20];
	for (int i = 10; i < 30; i++)
		z2[i - 10] = X[i] - Hy7_O[i];
	mul_matrix20(z2, Hy72_M20);

	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt(ff1 / Dim2)) - exp(ff2 / Dim2) + exp(1) + 20;

	double z3[20];
	for (int i = 30; i < 50; i++)
		z3[i - 30] = X[i] - Hy7_O[i];
	mul_matrix20(z3, Hy73_M20);
	//Expanded Griewank's plus Rosenbrock
	double temp[20];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[20];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[20];
	double yy[20];
	double xz[20];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 30]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[20];
	for (int i = 50; i < 70; i++)
		z4[i - 50] = X[i] - Hy7_O[i];
	mul_matrix20(z4, Hy74_M20);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f4 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f4 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f4 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f4 = 418.9829*Dim4 - f4;


	double z5[30];
	for (int i = 70; i < 100; i++)
		z5[i - 70] = X[i] - Hy7_O[i];
	mul_matrix30(z5, Hy75_M30);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f5 = f5 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f = f1 + f2 + f3 + f4 + f5;
#endif

	return f + 1700;
#endif

#if hybrid_function_8
	//F(X)=High_Conditioned_Elliptic + Ackley + Rastrigin + HGBat +Discus	
	double p[5] = { 0.2,0.2,0.2,0.2,0.2 };

	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;
	double f = 0;

	extern double Hy8_O[100];
	extern double Hy8_M10[10][10];
	extern double Hy8_M20[20][20];
	extern double Hy8_M30[30][30];
	extern double Hy8_M50[50][50];
	extern double Hy8_M100[100][100];

#if DIM==10
	double Hy81_M2[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Hy81_M2[i][j] = Hy8_M10[i][j];
	double Hy82_M2[2][2];
	for (int i = 2; i < 4; i++)
		for (int j = 2; j < 4; j++)
			Hy82_M2[i - 2][j - 2] = Hy8_M10[i][j];
	double Hy83_M2[2][2];
	for (int i = 4; i < 6; i++)
		for (int j = 4; j < 6; j++)
			Hy83_M2[i - 4][j - 4] = Hy8_M10[i][j];
	double Hy84_M2[2][2];
	for (int i = 6; i < 8; i++)
		for (int j = 6; j < 8; j++)
			Hy84_M2[i - 6][j - 6] = Hy8_M10[i][j];
	double Hy85_M2[2][2];
	for (int i = 8; i < 10; i++)
		for (int j = 8; j < 10; j++)
			Hy85_M2[i - 8][j - 8] = Hy8_M10[i][j];

	double z1[2];
	for (int i = 0; i < 2; i++)
		z1[i] = X[i] - Hy8_O[i];
	mul_matrix02(z1, Hy81_M2);
	for (int i = 0; i < 2; i++)
		f1 += pow(pow(10, 6), i / (2 - 1))*z1[i] * z1[i];

	double z2[2];
	for (int i = 2; i < 4; i++)
		z2[i - 2] = X[i] - Hy8_O[i];
	mul_matrix02(z2, Hy82_M2);
	double ff1 = 0;
	double ff2 = 0;
	for (int i = 0; i < 2; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / 2) * ff1)) - exp((1.0 / 2) * ff2) + exp(1) + 20;

	double z3[2];
	for (int i = 4; i < 6; i++)
		z3[i - 4] = X[i] - Hy8_O[i];
	mul_matrix02(z3, Hy83_M2);
	for (int i = 0; i < 2; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	double z4[2];
	for (int i = 6; i < 8; i++)
		z4[i - 6] = X[i] - Hy8_O[i];
	mul_matrix02(z4, Hy84_M2);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < 2; i++)
	{
		fd1 += z4[i] * z4[i];
		fd2 += z4[i];
	}
	f4 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / 2 + 0.5;

	//Discus
	double z5[2];
	for (int i = 8; i < 10; i++)
		z5[i - 8] = X[i] - Hy8_O[i];
	mul_matrix02(z5, Hy85_M2);
	double fz1 = 0.0;
	for (int i = 1; i < 2; i++)
		fz1 += z5[i] * z5[i];
	f5 = 1E6 * z5[0] * z5[0] + fz1;

	f = f1 + f2 + f3 + f4 + f5;
#endif

#if DIM==20
	double Hy81_M4[4][4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			Hy81_M4[i][j] = Hy8_M20[i][j];
	double Hy82_M4[4][4];
	for (int i = 4; i < 8; i++)
		for (int j = 4; j < 8; j++)
			Hy82_M4[i - 4][j - 4] = Hy8_M20[i][j];
	double Hy83_M4[4][4];
	for (int i = 8; i < 12; i++)
		for (int j = 8; j < 12; j++)
			Hy83_M4[i - 8][j - 8] = Hy8_M20[i][j];
	double Hy84_M4[4][4];
	for (int i = 12; i < 16; i++)
		for (int j = 12; j < 16; j++)
			Hy84_M4[i - 12][j - 12] = Hy8_M20[i][j];
	double Hy85_M4[4][4];
	for (int i = 16; i < 20; i++)
		for (int j = 16; j < 20; j++)
			Hy85_M4[i - 16][j - 16] = Hy8_M20[i][j];

	double z1[4];
	for (int i = 0; i < 4; i++)
		z1[i] = X[i] - Hy8_O[i];
	mul_matrix04(z1, Hy81_M4);
	for (int i = 0; i < 4; i++)
		f1 += pow(pow(10, 6), i / (4 - 1))*z1[i] * z1[i];

	double z2[4];
	for (int i = 4; i < 8; i++)
		z2[i - 4] = X[i] - Hy8_O[i];
	mul_matrix04(z2, Hy82_M4);
	double ff1 = 0;
	double ff2 = 0;
	for (int i = 0; i < 4; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / 4) * ff1)) - exp((1.0 / 4) * ff2) + exp(1) + 20;

	double z3[4];
	for (int i = 8; i < 12; i++)
		z3[i - 8] = X[i] - Hy8_O[i];
	mul_matrix04(z3, Hy83_M4);
	for (int i = 0; i < 4; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	double z4[4];
	for (int i = 12; i < 16; i++)
		z4[i - 12] = X[i] - Hy8_O[i];
	mul_matrix04(z4, Hy84_M4);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < 4; i++)
	{
		fd1 += z4[i] * z4[i];
		fd2 += z4[i];
	}
	f4 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / 4 + 0.5;

	//Discus
	double z5[4];
	for (int i = 16; i < 20; i++)
		z5[i - 16] = X[i] - Hy8_O[i];
	mul_matrix04(z5, Hy85_M4);
	double fz1 = 0.0;
	for (int i = 1; i < 4; i++)
		fz1 += z5[i] * z5[i];
	f5 = 1E6 * z5[0] * z5[0] + fz1;

	f = f1 + f2 + f3 + f4 + f5;
#endif


#if DIM==30
	double Hy81_M6[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Hy81_M6[i][j] = Hy8_M30[i][j];
	double Hy82_M6[6][6];
	for (int i = 6; i < 12; i++)
		for (int j = 6; j < 12; j++)
			Hy82_M6[i - 6][j - 6] = Hy8_M30[i][j];
	double Hy83_M6[6][6];
	for (int i = 12; i < 18; i++)
		for (int j = 12; j < 18; j++)
			Hy83_M6[i - 12][j - 12] = Hy8_M30[i][j];
	double Hy84_M6[6][6];
	for (int i = 18; i < 24; i++)
		for (int j = 18; j < 24; j++)
			Hy84_M6[i - 18][j - 18] = Hy8_M30[i][j];
	double Hy85_M6[6][6];
	for (int i = 24; i < 30; i++)
		for (int j = 24; j < 30; j++)
			Hy85_M6[i - 24][j - 24] = Hy8_M30[i][j];

	double z1[6];
	for (int i = 0; i < 6; i++)
		z1[i] = X[i] - Hy8_O[i];
	mul_matrix06(z1, Hy81_M6);
	for (int i = 0; i < 6; i++)
		f1 += pow(pow(10, 6), i / (6 - 1))*z1[i] * z1[i];

	double z2[6];
	for (int i = 6; i < 12; i++)
		z2[i - 6] = X[i] - Hy8_O[i];
	mul_matrix06(z2, Hy82_M6);
	double ff1 = 0;
	double ff2 = 0;
	for (int i = 0; i < 6; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / 6) * ff1)) - exp((1.0 / 6) * ff2) + exp(1) + 20;

	double z3[6];
	for (int i = 12; i < 18; i++)
		z3[i - 12] = X[i] - Hy8_O[i];
	mul_matrix06(z3, Hy83_M6);
	for (int i = 0; i < 6; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	double z4[6];
	for (int i = 18; i < 24; i++)
		z4[i - 18] = X[i] - Hy8_O[i];
	mul_matrix06(z4, Hy84_M6);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < 6; i++)
	{
		fd1 += z4[i] * z4[i];
		fd2 += z4[i];
	}
	f4 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / 6 + 0.5;

	//Discus
	double z5[6];
	for (int i = 24; i < 30; i++)
		z5[i - 24] = X[i] - Hy8_O[i];
	mul_matrix06(z5, Hy85_M6);
	double fz1 = 0.0;
	for (int i = 1; i < 6; i++)
		fz1 += z5[i] * z5[i];
	f5 = 1E6 * z5[0] * z5[0] + fz1;

	f = f1 + f2 + f3 + f4 + f5;
#endif
#if DIM==50
	double Hy81_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy81_M10[i][j] = Hy8_M50[i][j];
	double Hy82_M10[10][10];
	for (int i = 10; i < 20; i++)
		for (int j = 10; j < 20; j++)
			Hy82_M10[i - 10][j - 10] = Hy8_M50[i][j];
	double Hy83_M10[10][10];
	for (int i = 20; i < 30; i++)
		for (int j = 20; j < 30; j++)
			Hy83_M10[i - 20][j - 20] = Hy8_M50[i][j];
	double Hy84_M10[10][10];
	for (int i = 30; i < 40; i++)
		for (int j = 30; j < 40; j++)
			Hy84_M10[i - 30][j - 30] = Hy8_M50[i][j];
	double Hy85_M10[10][10];
	for (int i = 40; i < 50; i++)
		for (int j = 40; j < 50; j++)
			Hy85_M10[i - 40][j - 40] = Hy8_M50[i][j];

	double z1[10];
	for (int i = 0; i < 10; i++)
		z1[i] = X[i] - Hy8_O[i];
	mul_matrix10(z1, Hy81_M10);
	for (int i = 0; i < 10; i++)
		f1 += pow(pow(10, 6), i / (10 - 1))*z1[i] * z1[i];

	double z2[10];
	for (int i = 10; i < 20; i++)
		z2[i - 10] = X[i] - Hy8_O[i];
	mul_matrix10(z2, Hy82_M10);
	double ff1 = 0;
	double ff2 = 0;
	for (int i = 0; i < 10; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / 10) * ff1)) - exp((1.0 / 10) * ff2) + exp(1) + 20;

	double z3[10];
	for (int i = 20; i < 30; i++)
		z3[i - 20] = X[i] - Hy8_O[i];
	mul_matrix10(z3, Hy83_M10);
	for (int i = 0; i < 10; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	double z4[10];
	for (int i = 30; i < 40; i++)
		z4[i - 30] = X[i] - Hy8_O[i];
	mul_matrix10(z4, Hy84_M10);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < 10; i++)
	{
		fd1 += z4[i] * z4[i];
		fd2 += z4[i];
	}
	f4 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / 10 + 0.5;

	//Discus
	double z5[10];
	for (int i = 40; i < 50; i++)
		z5[i - 40] = X[i] - Hy8_O[i];
	mul_matrix10(z5, Hy85_M10);
	double fz1 = 0.0;
	for (int i = 1; i < 10; i++)
		fz1 += z5[i] * z5[i];
	f5 = 1E6 * z5[0] * z5[0] + fz1;

	f = f1 + f2 + f3 + f4 + f5;
#endif

#if DIM==100
	double Hy81_M20[20][20];
	for (int i = 0; i < 20; i++)
		for (int j = 0; j < 20; j++)
			Hy81_M20[i][j] = Hy8_M100[i][j];
	double Hy82_M20[20][20];
	for (int i = 20; i < 40; i++)
		for (int j = 20; j < 40; j++)
			Hy82_M20[i - 20][j - 20] = Hy8_M100[i][j];
	double Hy83_M20[20][20];
	for (int i = 40; i < 60; i++)
		for (int j = 40; j < 60; j++)
			Hy83_M20[i - 40][j - 40] = Hy8_M100[i][j];
	double Hy84_M20[20][20];
	for (int i = 60; i < 80; i++)
		for (int j = 60; j < 80; j++)
			Hy84_M20[i - 60][j - 60] = Hy8_M100[i][j];
	double Hy85_M20[20][20];
	for (int i = 80; i < 100; i++)
		for (int j = 80; j < 100; j++)
			Hy85_M20[i - 80][j - 80] = Hy8_M100[i][j];

	double z1[20];
	for (int i = 0; i < 20; i++)
		z1[i] = X[i] - Hy8_O[i];
	mul_matrix20(z1, Hy81_M20);
	for (int i = 0; i < 20; i++)
		f1 += pow(pow(10, 6), i / (20 - 1))*z1[i] * z1[i];

	double z2[20];
	for (int i = 20; i < 40; i++)
		z2[i - 20] = X[i] - Hy8_O[i];
	mul_matrix20(z2, Hy82_M20);
	double ff1 = 0;
	double ff2 = 0;
	for (int i = 0; i < 20; i++)
	{
		ff1 += z2[i] * z2[i];
		ff2 += cos(2 * pai * z2[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / 20) * ff1)) - exp((1.0 / 20 * ff2) + exp(1) + 20;

	double z3[20];
	for (int i = 40; i < 60; i++)
		z3[i - 40] = X[i] - Hy8_O[i];
	mul_matrix20(z3, Hy83_M20);
	for (int i = 0; i < 20; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	double z4[20];
	for (int i = 60; i < 80; i++)
		z4[i - 60] = X[i] - Hy8_O[i];
	mul_matrix20(z4, Hy84_M20);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < 20; i++)
	{
		fd1 += z4[i] * z4[i];
		fd2 += z4[i];
	}
	f4 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / 20 + 0.5;

	//Discus
	double z5[20];
	for (int i = 80; i < 100; i++)
		z5[i - 80] = X[i] - Hy8_O[i];
	mul_matrix20(z5, Hy85_M20);
	double fz1 = 0.0;
	for (int i = 1; i < 20; i++)
		fz1 += z5[i] * z5[i];
	f5 = 1E6 * z5[0] * z5[0] + fz1;

	f = f1 + f2 + f3 + f4 + f5;
#endif

	return f + 1800;

#endif

#if hybrid_function_9
	//F(X)=BentCigar + Rastrigin + expanded Griewank's plus Rosenbrock + Weierstrass + Expanded Schaffer6
	double p[5] = { 0.2,0.2,0.2,0.2,0.2 };
	int Dim1 = DIM / 10 * 2;
	int Dim2 = DIM / 10 * 2;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 2;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;
	double f = 0;

	extern double Hy9_O[100];
	extern double Hy9_M10[10][10];
	extern double Hy9_M20[20][20];
	extern double Hy9_M30[30][30];
	extern double Hy9_M50[50][50];
	extern double Hy9_M100[100][100];


#if DIM==10
	double Hy91_M2[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Hy91_M2[i][j] = Hy9_M10[i][j];
	double Hy92_M2[2][2];
	for (int i = 2; i < 4; i++)
		for (int j = 2; j < 4; j++)
			Hy92_M2[i - 2][j - 2] = Hy9_M10[i][j];
	double Hy93_M2[2][2];
	for (int i = 4; i < 6; i++)
		for (int j = 4; j < 6; j++)
			Hy93_M2[i - 4][j - 4] = Hy9_M10[i][j];
	double Hy94_M2[2][2];
	for (int i = 6; i < 8; i++)
		for (int j = 6; j < 8; j++)
			Hy94_M2[i - 6][j - 6] = Hy9_M10[i][j];
	double Hy95_M2[2][2];
	for (int i = 8; i < 10; i++)
		for (int j = 8; j < 10; j++)
			Hy95_M2[i - 8][j - 8] = Hy9_M10[i][j];
	//Bent Cigar
	double z1[2];
	for (int i = 0; i < 2; i++)
		z1[i] = X[i] - Hy9_O[i];
	mul_matrix02(z1, Hy91_M2);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	//Rastrigin
	double z2[2];
	for (int i = 2; i < 4; i++)
		z2[i - 2] = X[i] - Hy9_O[i];
	mul_matrix02(z2, Hy92_M2);
	for (int i = 0; i < 2; i++)
		f2 = f2 + (z2[i] * z2[i] - 10 * cos(2 * pai*z2[i]) + 10);

	double z3[2];
	for (int i = 4; i < 6; i++)
		z3[i - 4] = X[i] - Hy9_O[i];
	mul_matrix02(z3, Hy93_M2);
	//Griewank's plus Rosenbrock
	double temp[2];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[2];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[2];
	double yy[2];
	double xz[2];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 4]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);

	//weierstrass
	double z4[2];
	for (int i = 6; i < 8; i++)
		z4[i - 6] = X[i] - Hy9_O[i];
	mul_matrix02(z4, Hy94_M2);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yy2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yy2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim4; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(z4[i] + 0.5));
		ff1 += yy1;
	}
	f4 = ff1 - Dim4 * yy2;



	double z5[2];
	for (int i = 8; i < 10; i++)
		z5[i - 8] = X[i] - Hy9_O[i];
	mul_matrix02(z5, Hy95_M2);
	double g = 0.0;
	for (int i = 0; i < Dim5 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z5[i] * z5[i] + z5[i + 1] * z5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[i] * z5[i] + z5[i + 1] * z5[i + 1])), 2);
		f5 += g;
	}
	f5 += 0.5 + (pow(sin(pow(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0])), 2);

	f = f1 + f2 + f3 + f4 + f5;
#endif

#if DIM==20
	double Hy91_M4[4][4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			Hy91_M4[i][j] = Hy9_M20[i][j];
	double Hy92_M4[4][4];
	for (int i = 4; i < 8; i++)
		for (int j = 4; j < 8; j++)
			Hy92_M4[i - 4][j - 4] = Hy9_M20[i][j];
	double Hy93_M4[4][4];
	for (int i = 8; i < 12; i++)
		for (int j = 8; j < 12; j++)
			Hy93_M4[i - 8][j - 8] = Hy9_M20[i][j];
	double Hy94_M4[4][4];
	for (int i = 12; i < 16; i++)
		for (int j = 12; j < 16; j++)
			Hy94_M4[i - 12][j - 12] = Hy9_M20[i][j];
	double Hy95_M4[4][4];
	for (int i = 16; i < 20; i++)
		for (int j = 16; j < 20; j++)
			Hy95_M4[i - 16][j - 16] = Hy9_M20[i][j];
	//Bent Cigar
	double z1[4];
	for (int i = 0; i < 4; i++)
		z1[i] = X[i] - Hy9_O[i];
	mul_matrix04(z1, Hy91_M4);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;

	//Rastrigin
	double z2[4];
	for (int i = 4; i < 8; i++)
		z2[i - 4] = X[i] - Hy9_O[i];
	mul_matrix04(z2, Hy92_M4);
	for (int i = 0; i < 4; i++)
		f2 = f2 + (z2[i] * z2[i] - 10 * cos(2 * pai*z2[i]) + 10);

	double z3[4];
	for (int i = 8; i < 12; i++)
		z3[i - 8] = X[i] - Hy9_O[i];
	mul_matrix04(z3, Hy93_M4);
	//Griewank's plus Rosenbrock
	double temp[4];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[4];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[4];
	double yy[4];
	double xz[4];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 8]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);

	//weierstrass
	double z4[4];
	for (int i = 12; i < 16; i++)
		z4[i - 12] = X[i] - Hy9_O[i];
	mul_matrix04(z4, Hy94_M4);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yy2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yy2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim4; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(z4[i] + 0.5));
		ff1 += yy1;
	}
	f4 = ff1 - Dim4 * yy2;



	double z5[4];
	for (int i = 16; i < 20; i++)
		z5[i - 16] = X[i] - Hy9_O[i];
	mul_matrix04(z5, Hy95_M4);
	double g = 0.0;
	for (int i = 0; i < Dim5 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z5[i] * z5[i] + z5[i + 1] * z5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[i] * z5[i] + z5[i + 1] * z5[i + 1])), 2);
		f5 += g;
	}
	f5 += 0.5 + (pow(sin(pow(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0])), 2);

	f = f1 + f2 + f3 + f4 + f5;
#endif

#if DIM==30
	double Hy91_M6[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Hy91_M6[i][j] = Hy9_M30[i][j];
	double Hy92_M6[6][6];
	for (int i = 6; i < 12; i++)
		for (int j = 6; j < 12; j++)
			Hy92_M6[i - 6][j - 6] = Hy9_M30[i][j];
	double Hy93_M6[6][6];
	for (int i = 12; i < 18; i++)
		for (int j = 12; j < 18; j++)
			Hy93_M6[i - 12][j - 12] = Hy9_M30[i][j];
	double Hy94_M6[6][6];
	for (int i = 18; i < 24; i++)
		for (int j = 18; j < 24; j++)
			Hy94_M6[i - 18][j - 18] = Hy9_M30[i][j];
	double Hy95_M6[6][6];
	for (int i = 24; i < 30; i++)
		for (int j = 24; j < 30; j++)
			Hy95_M6[i - 24][j - 24] = Hy9_M30[i][j];

	double z1[6];
	for (int i = 0; i < 6; i++)
		z1[i] = X[i] - Hy9_O[i];
	mul_matrix06(z1, Hy91_M6);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	double z2[6];
	for (int i = 6; i < 12; i++)
		z2[i - 6] = X[i] - Hy9_O[i];
	mul_matrix06(z2, Hy92_M6);
	for (int i = 0; i < 6; i++)
		f2 = f2 + (z2[i] * z2[i] - 10 * cos(2 * pai*z2[i]) + 10);

	double z3[6];
	for (int i = 12; i < 18; i++)
		z3[i - 12] = X[i] - Hy9_O[i];
	mul_matrix06(z3, Hy93_M6);
	//Griewank's plus Rosenbrock
	double temp[6];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[6];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[6];
	double yy[6];
	double xz[6];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 12]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;
	}
	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[6];
	for (int i = 18; i < 24; i++)
		z4[i - 18] = X[i] - Hy9_O[i];
	mul_matrix06(z4, Hy94_M6);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yy2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yy2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim4; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(z4[i] + 0.5));
		ff1 += yy1;
	}
	f4 = ff1 - Dim4 * yy2;



	double z5[6];
	for (int i = 24; i < 30; i++)
		z5[i - 24] = X[i] - Hy9_O[i];
	mul_matrix06(z5, Hy95_M6);
	double g = 0.0;
	for (int i = 0; i < Dim5 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z5[i] * z5[i] + z5[i + 1] * z5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[i] * z5[i] + z5[i + 1] * z5[i + 1])), 2);
		f5 += g;
	}
	f5 += 0.5 + (pow(sin(pow(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0])), 2);

	f = f1 + f2 + f3 + f4 + f5;
#endif

#if DIM==50
	double Hy91_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy91_M10[i][j] = Hy9_M50[i][j];
	double Hy92_M10[10][10];
	for (int i = 10; i < 20; i++)
		for (int j = 10; j < 20; j++)
			Hy92_M10[i - 10][j - 10] = Hy9_M50[i][j];
	double Hy93_M10[10][10];
	for (int i = 20; i < 30; i++)
		for (int j = 20; j < 30; j++)
			Hy93_M10[i - 20][j - 20] = Hy9_M50[i][j];
	double Hy94_M10[10][10];
	for (int i = 30; i < 40; i++)
		for (int j = 30; j < 40; j++)
			Hy94_M10[i - 30][j - 30] = Hy9_M50[i][j];
	double Hy95_M10[10][10];
	for (int i = 40; i < 50; i++)
		for (int j = 40; j < 50; j++)
			Hy95_M10[i - 40][j - 40] = Hy9_M50[i][j];

	double z1[10];
	for (int i = 0; i < 10; i++)
		z1[i] = X[i] - Hy9_O[i];
	mul_matrix10(z1, Hy91_M10);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	double z2[10];
	for (int i = 10; i < 20; i++)
		z2[i - 10] = X[i] - Hy9_O[i];
	mul_matrix10(z2, Hy92_M10);
	for (int i = 0; i < 10; i++)
		f2 = f2 + (z2[i] * z2[i] - 10 * cos(2 * pai*z2[i]) + 10);

	double z3[10];
	for (int i = 20; i < 30; i++)
		z3[i - 20] = X[i] - Hy9_O[i];
	mul_matrix10(z3, Hy93_M10);
	//Griewank's plus Rosenbrock
	double temp[10];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[10];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[10];
	double yy[10];
	double xz[10];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 20]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[10];
	for (int i = 30; i < 40; i++)
		z4[i - 30] = X[i] - Hy9_O[i];
	mul_matrix10(z4, Hy94_M10);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yy2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yy2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim4; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(z4[i] + 0.5));
		ff1 += yy1;
	}
	f4 = ff1 - Dim4 * yy2;



	double z5[10];
	for (int i = 40; i < 50; i++)
		z5[i - 40] = X[i] - Hy9_O[i];
	mul_matrix10(z5, Hy95_M10);
	double g = 0.0;
	for (int i = 0; i < Dim5 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z5[i] * z5[i] + z5[i + 1] * z5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[i] * z5[i] + z5[i + 1] * z5[i + 1])), 2);
		f5 += g;
	}
	f5 += 0.5 + (pow(sin(pow(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0])), 2);

	f = f1 + f2 + f3 + f4 + f5;
#endif

#if DIM==100
	double Hy91_M20[20][20];
	for (int i = 0; i < 20; i++)
		for (int j = 0; j < 20; j++)
			Hy91_M20[i][j] = Hy9_M100[i][j];
	double Hy92_M20[20][20];
	for (int i = 20; i < 40; i++)
		for (int j = 20; j < 40; j++)
			Hy92_M20[i - 20][j - 20] = Hy9_M100[i][j];
	double Hy93_M20[20][20];
	for (int i = 40; i < 60; i++)
		for (int j = 40; j < 60; j++)
			Hy93_M20[i - 40][j - 40] = Hy9_M100[i][j];
	double Hy94_M20[20][20];
	for (int i = 60; i < 80; i++)
		for (int j = 60; j < 80; j++)
			Hy94_M20[i - 60][j - 60] = Hy9_M100[i][j];
	double Hy95_M20[20][20];
	for (int i = 80; i < 100; i++)
		for (int j = 80; j < 100; j++)
			Hy95_M20[i - 80][j - 80] = Hy9_M100[i][j];

	double z1[20];
	for (int i = 0; i < 20; i++)
		z1[i] = X[i] - Hy9_O[i];
	mul_matrix20(z1, Hy91_M20);
	double y1 = z1[0] * z1[0];
	double y2 = 0.0;
	for (int i = 1; i < Dim1; i++)
		y2 += 1E6*z1[i] * z1[i];
	f1 = y1 + y2;


	double z2[20];
	for (int i = 20; i < 40; i++)
		z2[i - 20] = X[i] - Hy9_O[i];
	mul_matrix20(z2, Hy92_M20);
	for (int i = 0; i < 20; i++)
		f2 = f2 + (z2[i] * z2[i] - 10 * cos(2 * pai*z2[i]) + 10);

	double z3[20];
	for (int i = 40; i < 60; i++)
		z3[i - 40] = X[i] - Hy9_O[i];
	mul_matrix20(z3, Hy93_M20);
	//Griewank's plus Rosenbrock
	double temp[20];
	for (int i = 0; i < Dim3 - 2; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	temp[Dim3 - 2] = 100 * (z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1])*(z3[Dim3 - 2] * z3[Dim3 - 2] - z3[Dim3 - 1]) + (z3[Dim3 - 2] - 1)*(z3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[20];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[20];
	double yy[20];
	double xz[20];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 40]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f3 = F1;
	else f3 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f3 = f3 + 10 * (Dim3 - ff);


	double z4[20];
	for (int i = 60; i < 80; i++)
		z4[i - 60] = X[i] - Hy9_O[i];
	mul_matrix20(z4, Hy94_M20);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yy2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yy2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim4; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(z4[i] + 0.5));
		ff1 += yy1;
	}
	f4 = ff1 - Dim4 * yy2;



	double z5[20];
	for (int i = 80; i < 100; i++)
		z5[i - 80] = X[i] - Hy9_O[i];
	mul_matrix20(z5, Hy95_M20);
	double g = 0.0;
	for (int i = 0; i < Dim5 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z5[i] * z5[i] + z5[i + 1] * z5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[i] * z5[i] + z5[i + 1] * z5[i + 1])), 2);
		f5 += g;
	}
	f5 += 0.5 + (pow(sin(pow(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z5[Dim5 - 1] * z5[Dim5 - 1] + z5[0] * z5[0])), 2);

	f = f1 + f2 + f3 + f4 + f5;
#endif

	return f + 1900;
#endif

#if hybrid_function_10
	//F(X)=HappyCat+Katsuura+Ackley+Rastrigin+Modified Schwefel's f10+Schaffer's f7
	double p[6] = { 0.1,0.1,0.2,0.2,0.2,0.2 };

	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;
	double f6 = 0;
	double f = 0;
	int Dim1 = DIM / 10 * 1;
	int Dim2 = DIM / 10 * 1;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 2;
	int Dim6 = DIM / 10 * 2;

	extern double Hy10_0[100];
	extern double Hy10_M10[10][10];
	extern double Hy10_M20[20][20];
	extern double Hy10_M30[30][30];
	extern double Hy10_M50[50][50];
	extern double Hy10_M100[100][100];

#if DIM==10
	double Hy101_M1[1][1];
	for (int i = 0; i < 1; i++)
		for (int j = 0; j < 1; j++)
			Hy101_M1[i][j] = Hy10_M10[i][j];
	double Hy102_M1[1][1];
	for (int i = 1; i < 2; i++)
		for (int j = 1; j < 2; j++)
			Hy102_M1[i - 1][j - 1] = Hy10_M10[i][j];
	double Hy103_M2[2][2];
	for (int i = 2; i < 4; i++)
		for (int j = 2; j < 4; j++)
			Hy103_M2[i - 2][j - 2] = Hy10_M10[i][j];
	double Hy104_M2[2][2];
	for (int i = 4; i < 6; i++)
		for (int j = 4; j < 6; j++)
			Hy104_M2[i - 4][j - 4] = Hy10_M10[i][j];
	double Hy105_M2[2][2];
	for (int i = 6; i < 8; i++)
		for (int j = 6; j < 8; j++)
			Hy105_M2[i - 6][j - 6] = Hy10_M10[i][j];
	double Hy106_M2[2][2];
	for (int i = 8; i < 10; i++)
		for (int j = 8; j < 10; j++)
			Hy106_M2[i - 8][j - 8] = Hy10_M10[i][j];

	double z1[1];
	for (int i = 0; i < 1; i++)
		z1[i] = X[i] - Hy10_0[i];
	mul_matrix01(z1, Hy101_M1);
	double fd1 = 0.0;
	double fd2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		fd1 += z1[i] * z1[i];
		fd2 += z1[i];
	}
	k = fabs(fd1 - Dim1);
	f1 = pow(k, 0.25) + (0.5*fd1 + fd2) / Dim1 + 0.5;

	double z2[1];
	for (int i = 1; i < 2; i++)
		z2[i - 1] = X[i] - Hy10_0[i];
	mul_matrix01(z2, Hy102_M1);
	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim2 * Dim2);
	y3 = 10 / pow(Dim2, 1.2);

	for (int i = 0; i < Dim2; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z2[i] - ceil(pow(2, j)*z2[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f2 = y1 * y - y1;

	double z3[2];
	for (int i = 2; i < 4; i++)
		z3[i - 2] = X[i] - Hy10_0[i];
	mul_matrix02(z3, Hy103_M2);
	double fx1 = 0;
	double fx2 = 0;
	for (int i = 0; i < Dim3; i++)
	{
		fx1 += z3[i] * z3[i];
		fx2 += cos(2 * pai * z3[i]);
	}
	f3 = -20 * exp(-0.2*sqrt((1.0 / Dim3) * fx1)) - exp((1.0 / Dim3) * fx2) + exp(1) + 20;

	double z4[2];
	for (int i = 4; i < 6; i++)
		z4[i - 4] = X[i] - Hy10_0[i];
	mul_matrix02(z4, Hy104_M2);
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	double z5[2];
	for (int i = 6; i < 8; i++)
		z5[i - 6] = X[i] - Hy10_0[i];
	mul_matrix02(z5, Hy105_M2);
	for (int i = 0; i < Dim4; i++)
	{
		z5[i] = z5[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z5[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z5[i]) <= 500)
			f5 += z5[i] * sin(pow(fabs(z5[i]), 0.5));
		if (z5[i] > 500)
			f5 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] - 500)*(z5[i] - 500) / (10000 * Dim5);
		if (z5[i] < -500)
			f5 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] + 500)*(z5[i] + 500) / (10000 * Dim5);
	}
	f5 = 418.9829*Dim5 - f5;

	double z6[2];
	for (int i = 8; i < 10; i++)
		z6[i - 8] = X[i] - Hy10_0[i];
	mul_matrix02(z6, Hy106_M2);
	double Si = 0;
	for (int i = 0; i < Dim6 - 1; i++)
	{
		Si = pow(z6[i] * z6[i] + z6[i + 1] * z6[i + 1], 0.5);
		f6 += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	f6 = f6 / (Dim6 - 1);
	f6 = f6 * f6;

	f = f1 + f2 + f3 + f4 + f5 + f6;
#endif

#if DIM==20
	double Hy101_M2[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Hy101_M2[i][j] = Hy10_M20[i][j];
	double Hy102_M2[2][2];
	for (int i = 2; i < 4; i++)
		for (int j = 2; j < 4; j++)
			Hy102_M2[i - 2][j - 2] = Hy10_M20[i][j];
	double Hy103_M4[4][4];
	for (int i = 4; i < 8; i++)
		for (int j = 4; j < 8; j++)
			Hy103_M4[i - 4][j - 4] = Hy10_M20[i][j];
	double Hy104_M4[4][4];
	for (int i = 8; i < 12; i++)
		for (int j = 8; j < 12; j++)
			Hy104_M4[i - 8][j - 8] = Hy10_M20[i][j];
	double Hy105_M4[4][4];
	for (int i = 12; i < 16; i++)
		for (int j = 12; j < 16; j++)
			Hy105_M4[i - 12][j - 12] = Hy10_M20[i][j];
	double Hy106_M4[4][4];
	for (int i = 16; i < 20; i++)
		for (int j = 16; j < 20; j++)
			Hy106_M4[i - 16][j - 16] = Hy10_M20[i][j];

	double z1[2];
	for (int i = 0; i < 2; i++)
		z1[i] = X[i] - Hy10_0[i];
	mul_matrix02(z1, Hy101_M2);
	double fd1 = 0.0;
	double fd2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		fd1 += z1[i] * z1[i];
		fd2 += z1[i];
	}
	k = fabs(fd1 - Dim1);
	f1 = pow(k, 0.25) + (0.5*fd1 + fd2) / Dim1 + 0.5;

	double z2[2];
	for (int i = 2; i < 4; i++)
		z2[i - 2] = X[i] - Hy10_0[i];
	mul_matrix02(z2, Hy102_M2);
	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim2 * Dim2);
	y3 = 10 / pow(Dim2, 1.2);

	for (int i = 0; i < Dim2; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z2[i] - ceil(pow(2, j)*z2[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f2 = y1 * y - y1;

	double z3[4];
	for (int i = 4; i < 8; i++)
		z3[i - 4] = X[i] - Hy10_0[i];
	mul_matrix04(z3, Hy103_M4);
	double fx1 = 0;
	double fx2 = 0;
	for (int i = 0; i < Dim3; i++)
	{
		fx1 += z3[i] * z3[i];
		fx2 += cos(2 * pai * z3[i]);
	}
	f3 = -20 * exp(-0.2*sqrt((1.0 / Dim3) * fx1)) - exp((1.0 / Dim3) * fx2) + exp(1) + 20;

	double z4[4];
	for (int i = 8; i < 12; i++)
		z4[i - 8] = X[i] - Hy10_0[i];
	mul_matrix04(z4, Hy104_M4);
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	double z5[4];
	for (int i = 12; i < 16; i++)
		z5[i - 12] = X[i] - Hy10_0[i];
	mul_matrix04(z5, Hy105_M4);
	for (int i = 0; i < Dim4; i++)
	{
		z5[i] = z5[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z5[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z5[i]) <= 500)
			f5 += z5[i] * sin(pow(fabs(z5[i]), 0.5));
		if (z5[i] > 500)
			f5 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] - 500)*(z5[i] - 500) / (10000 * Dim5);
		if (z5[i] < -500)
			f5 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] + 500)*(z5[i] + 500) / (10000 * Dim5);
	}
	f5 = 418.9829*Dim5 - f5;

	double z6[4];
	for (int i = 16; i < 20; i++)
		z6[i - 16] = X[i] - Hy10_0[i];
	mul_matrix04(z6, Hy106_M4);
	double Si = 0;
	for (int i = 0; i < Dim6 - 1; i++)
	{
		Si = pow(z6[i] * z6[i] + z6[i + 1] * z6[i + 1], 0.5);
		f6 += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.5)) + 1);
	}
	f6 = f6 / (Dim6 - 1);
	f6 = f6 * f6;

	f = f1 + f2 + f3 + f4 + f5 + f6;
#endif

#if DIM==30
	double Hy101_M3[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Hy101_M3[i][j] = Hy10_M30[i][j];
	double Hy102_M3[3][3];
	for (int i = 3; i < 6; i++)
		for (int j = 3; j < 6; j++)
			Hy102_M3[i - 3][j - 3] = Hy10_M30[i][j];
	double Hy103_M6[6][6];
	for (int i = 6; i < 12; i++)
		for (int j = 6; j < 12; j++)
			Hy103_M6[i - 6][j - 6] = Hy10_M30[i][j];
	double Hy104_M6[6][6];
	for (int i = 12; i < 18; i++)
		for (int j = 12; j < 18; j++)
			Hy104_M6[i - 12][j - 12] = Hy10_M30[i][j];
	double Hy105_M6[6][6];
	for (int i = 18; i < 24; i++)
		for (int j = 18; j < 24; j++)
			Hy105_M6[i - 18][j - 18] = Hy10_M30[i][j];
	double Hy106_M6[6][6];
	for (int i = 24; i < 30; i++)
		for (int j = 24; j < 30; j++)
			Hy106_M6[i - 24][j - 24] = Hy10_M30[i][j];

	double z1[3];
	for (int i = 0; i < 3; i++)
		z1[i] = X[i] - Hy10_0[i];
	mul_matrix03(z1, Hy101_M3);
	double fd1 = 0.0;
	double fd2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		fd1 += z1[i] * z1[i];
		fd2 += z1[i];
	}
	k = fabs(fd1 - Dim1);
	f1 = pow(k, 0.25) + (0.5*fd1 + fd2) / Dim1 + 0.5;

	double z2[3];
	for (int i = 3; i < 6; i++)
		z2[i - 3] = X[i] - Hy10_0[i];
	mul_matrix03(z2, Hy102_M3);
	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim2 * Dim2);
	y3 = 10 / pow(Dim2, 1.2);

	for (int i = 0; i < Dim2; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z2[i] - ceil(pow(2, j)*z2[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f2 = y1 * y - y1;

	double z3[6];
	for (int i = 6; i < 12; i++)
		z3[i - 6] = X[i] - Hy10_0[i];
	mul_matrix06(z3, Hy103_M6);
	double fx1 = 0;
	double fx2 = 0;
	for (int i = 0; i < Dim3; i++)
	{
		fx1 += z3[i] * z3[i];
		fx2 += cos(2 * pai * z3[i]);
	}
	f3 = -20 * exp(-0.2*sqrt((1.0 / Dim3) * fx1)) - exp((1.0 / Dim3) * fx2) + exp(1) + 20;

	double z4[6];
	for (int i = 12; i < 18; i++)
		z4[i - 12] = X[i] - Hy10_0[i];
	mul_matrix06(z4, Hy104_M6);
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	double z5[6];
	for (int i = 18; i < 24; i++)
		z5[i - 18] = X[i] - Hy10_0[i];
	mul_matrix06(z5, Hy105_M6);
	for (int i = 0; i < Dim4; i++)
	{
		z5[i] = z5[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z5[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z5[i]) <= 500)
			f5 += z5[i] * sin(pow(fabs(z5[i]), 0.5));
		if (z5[i] > 500)
			f5 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] - 500)*(z5[i] - 500) / (10000 * Dim5);
		if (z5[i] < -500)
			f5 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] + 500)*(z5[i] + 500) / (10000 * Dim5);
	}
	f5 = 418.9829*Dim5 - f5;

	double z6[6];
	for (int i = 24; i < 30; i++)
		z6[i - 24] = X[i] - Hy10_0[i];
	mul_matrix06(z6, Hy106_M6);
	double Si = 0;
	for (int i = 0; i < Dim6 - 1; i++)
	{
		Si = pow(z6[i] * z6[i] + z6[i + 1] * z6[i + 1], 0.5);
		f6 += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	f6 = f6 / (Dim6 - 1);
	f6 = f6 * f6;

	f = f1 + f2 + f3 + f4 + f5 + f6;
#endif
#if DIM==50
	double Hy101_M5[5][5];
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			Hy101_M5[i][j] = Hy10_M50[i][j];
	double Hy102_M5[5][5];
	for (int i = 5; i < 10; i++)
		for (int j = 5; j < 10; j++)
			Hy102_M5[i - 5][j - 5] = Hy10_M50[i][j];
	double Hy103_M10[10][10];
	for (int i = 10; i < 20; i++)
		for (int j = 10; j < 20; j++)
			Hy103_M10[i - 10][j - 10] = Hy10_M50[i][j];
	double Hy104_M10[10][10];
	for (int i = 20; i < 30; i++)
		for (int j = 20; j < 30; j++)
			Hy104_M10[i - 20][j - 20] = Hy10_M50[i][j];
	double Hy105_M10[10][10];
	for (int i = 30; i < 40; i++)
		for (int j = 30; j < 40; j++)
			Hy105_M10[i - 30][j - 30] = Hy10_M50[i][j];
	double Hy106_M10[10][10];
	for (int i = 40; i < 50; i++)
		for (int j = 40; j < 50; j++)
			Hy106_M10[i - 40][j - 40] = Hy10_M50[i][j];

	double z1[5];
	for (int i = 0; i < 5; i++)
		z1[i] = X[i] - Hy10_0[i];
	mul_matrix05(z1, Hy101_M5);
	double fd1 = 0.0;
	double fd2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		fd1 += z1[i] * z1[i];
		fd2 += z1[i];
	}
	k = fabs(fd1 - Dim1);
	f1 = pow(k, 0.25) + (0.5*fd1 + fd2) / Dim1 + 0.5;

	double z2[5];
	for (int i = 5; i < 10; i++)
		z2[i - 5] = X[i] - Hy10_0[i];
	mul_matrix05(z2, Hy102_M5);
	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim2 * Dim2);
	y3 = 10 / pow(Dim2, 1.2);

	for (int i = 0; i < Dim2; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z2[i] - ceil(pow(2, j)*z2[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f2 = y1 * y - y1;

	double z3[10];
	for (int i = 10; i < 20; i++)
		z3[i - 10] = X[i] - Hy10_0[i];
	mul_matrix10(z3, Hy103_M10);
	double fx1 = 0;
	double fx2 = 0;
	for (int i = 0; i < Dim3; i++)
	{
		fx1 += z3[i] * z3[i];
		fx2 += cos(2 * pai * z3[i]);
	}
	f3 = -20 * exp(-0.2*sqrt((1.0 / Dim3) * fx1)) - exp((1.0 / Dim3) * fx2) + exp(1) + 20;

	double z4[10];
	for (int i = 20; i < 30; i++)
		z4[i - 20] = X[i] - Hy10_0[i];
	mul_matrix10(z4, Hy104_M10);
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	double z5[10];
	for (int i = 30; i < 40; i++)
		z5[i - 30] = X[i] - Hy10_0[i];
	mul_matrix10(z5, Hy105_M10);
	for (int i = 0; i < Dim4; i++)
	{
		z5[i] = z5[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z5[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z5[i]) <= 500)
			f5 += z5[i] * sin(pow(fabs(z5[i]), 0.5));
		if (z5[i] > 500)
			f5 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] - 500)*(z5[i] - 500) / (10000 * Dim5);
		if (z5[i] < -500)
			f5 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] + 500)*(z5[i] + 500) / (10000 * Dim5);
	}
	f5 = 418.9829*Dim5 - f5;

	double z6[10];
	for (int i = 40; i < 50; i++)
		z6[i - 40] = X[i] - Hy10_0[i];
	mul_matrix10(z6, Hy106_M10);
	double Si = 0;
	for (int i = 0; i < Dim6 - 1; i++)
	{
		Si = pow(z6[i] * z6[i] + z6[i + 1] * z6[i + 1], 0.5);
		f6 += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2)) + 1);
	}
	f6 = f6 / (Dim6 - 1);
	f6 = f6 * f6;

	f = f1 + f2 + f3 + f4 + f5 + f6;
#endif

#if DIM==100
	double Hy101_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy101_M10[i][j] = Hy10_M100[i][j];
	double Hy102_M10[10][10];
	for (int i = 10; i < 20; i++)
		for (int j = 10; j < 20; j++)
			Hy102_M10[i - 10][j - 10] = Hy10_M100[i][j];
	double Hy103_M20[20][20];
	for (int i = 20; i < 40; i++)
		for (int j = 20; j < 40; j++)
			Hy103_M20[i - 20][j - 20] = Hy10_M100[i][j];
	double Hy104_M20[20][20];
	for (int i = 40; i < 60; i++)
		for (int j = 40; j < 60; j++)
			Hy104_M20[i - 40][j - 40] = Hy10_M100[i][j];
	double Hy105_M20[20][20];
	for (int i = 60; i < 80; i++)
		for (int j = 60; j < 80; j++)
			Hy105_M20[i - 60][j - 60] = Hy10_M100[i][j];
	double Hy106_M20[20][20];
	for (int i = 80; i < 100; i++)
		for (int j = 80; j < 100; j++)
			Hy106_M20[i - 80][j - 80] = Hy10_M100[i][j];

	double z1[10];
	for (int i = 0; i < 10; i++)
		z1[i] = X[i] - Hy10_0[i];
	mul_matrix10(z1, Hy101_M10);
	double fd1 = 0.0;
	double fd2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < Dim1; i++)
	{
		fd1 += z1[i] * z1[i];
		fd2 += z1[i];
	}
	k = fabs(fd1 - Dim1);
	f1 = pow(k, 0.25) + (0.5*fd1 + fd2) / Dim1 + 0.5;

	double z2[10];
	for (int i = 10; i < 20; i++)
		z2[i - 10] = X[i] - Hy10_0[i];
	mul_matrix10(z2, Hy102_M10);
	double y = 1.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	y1 = 10 / (Dim2 * Dim2);
	y3 = 10 / pow(Dim2, 1.2);

	for (int i = 0; i < Dim2; i++)
	{
		for (int j = 1; j <= 32; j++)
			y2 += (fabs(pow(2, j)*z2[i] - ceil(pow(2, j)*z2[i]))) / pow(2, j);
		y *= pow((1 + (i + 1)*y2), y3);
	}
	f2 = y1 * y - y1;

	double z3[20];
	for (int i = 20; i < 40; i++)
		z3[i - 20] = X[i] - Hy10_0[i];
	mul_matrix20(z3, Hy103_M20);
	double fx1 = 0;
	double fx2 = 0;
	for (int i = 0; i < Dim3; i++)
	{
		fx1 += z3[i] * z3[i];
		fx2 += cos(2 * pai * z3[i]);
	}
	f3 = -20 * exp(-0.2*sqrt((1.0 / Dim3) * fx1)) - exp((1.0 / Dim3) * fx2) + exp(1) + 20;

	double z4[20];
	for (int i = 40; i < 60; i++)
		z4[i - 40] = X[i] - Hy10_0[i];
	mul_matrix20(z4, Hy104_M20);
	for (int i = 0; i < Dim4; i++)
		f4 = f4 + (z4[i] * z4[i] - 10 * cos(2 * pai*z4[i]) + 10);

	double z5[20];
	for (int i = 60; i < 80; i++)
		z5[i - 60] = X[i] - Hy10_0[i];
	mul_matrix20(z5, Hy105_M20);
	for (int i = 0; i < Dim4; i++)
	{
		z5[i] = z5[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z5[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z5[i]) <= 500)
			f5 += z5[i] * sin(pow(fabs(z5[i]), 0.5));
		if (z5[i] > 500)
			f5 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] - 500)*(z5[i] - 500) / (10000 * Dim5);
		if (z5[i] < -500)
			f5 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z5[i] + 500)*(z5[i] + 500) / (10000 * Dim5);
	}
	f5 = 418.9829*Dim5 - f5;

	double z6[20];
	for (int i = 80; i < 100; i++)
		z6[i - 80] = X[i] - Hy10_0[i];
	mul_matrix20(z6, Hy106_M20);
	double Si = 0;
	for (int i = 0; i < Dim6 - 1; i++)
	{
		Si = pow(z6[i] * z6[i] + z6[i + 1] * z6[i + 1], 0.5);
		f6 += pow(Si, 0.5)*(sin(50.0*pow(Si, 0.2) + 1);
	}
	f6 = f6 / (Dim6 - 1);
	f6 = f6 * f6;

	f = f1 + f2 + f3 + f4 + f5 + f6;
#endif
	return f + 2000;
#endif

#if Composition_function_1
	//F(X)=Rosenbrock+High_Conditioned_Elliptic+Rastrigin
	double R[3] = { 10,20,30 };
	double T[3] = { 1,1E-6,1 };
	double bias[3] = { 0,100,200 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	extern double C11_O[100];
	extern double C12_O[100];
	extern double C13_O[100];

	double z1[DIM];
	for (int i = 0; i < DIM; i++)
		z1[i] = X[i];
	for (int i = 0; i < DIM; i++)
		y1 += (z1[i] - C11_O[i]) * (z1[i] - C11_O[i]);


	double z2[DIM];
	for (int i = 0; i < DIM; i++)
		z2[i] = X[i];
	for (int i = 0; i < DIM; i++)
		y2 += (z2[i] - C12_O[i]) * (z2[i] - C12_O[i]);


	double z3[DIM];
	for (int i = 0; i < DIM; i++)
		z3[i] = X[i];

	for (int i = 0; i < DIM; i++)
		y3 += (z3[i] - C13_O[i]) * (z3[i] - C13_O[i]);

	for (int i = 0; i < DIM - 1; i++)
		f1 += 100 * (z1[i] * z1[i] - z1[i + 1]) * (z1[i] * z1[i] - z1[i + 1]) + (z1[i] - 1) * (z1[i] - 1);


	for (int i = 1; i < DIM; i++)
		f2 += pow(1E6, i / (DIM - 1))*z2[i] * z2[i];


	for (int i = 0; i < DIM; i++)
		f3 = f3 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w = w1 + w2 + w3;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);

	return f + 2100;

#endif


#if Composition_function_2
	//F(X)=Rastrigin + Griewank + Modifed Schwefel's f10
	double R[3] = { 10,20,30 };
	double T[3] = { 1,10,1 };
	double bias[3] = { 0,100,200 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;

	double y1 = 0;
	double y2 = 0;
	double y3 = 0;

	double w1 = 0;
	double w2 = 0;
	double w3 = 0;

	double w = 0;
	extern double C21_O[100];
	extern double C22_O[100];
	extern double C23_O[100];

	double z1[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z1[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y1 += (z1[i] - C21_O[i]) *(z1[i] - C21_O[i]);

	//for (int i = 0; i < DIM - 1; i++)
	for (int i = 0; i < DIM; i++)
		f1 = f1 + (z1[i] * z1[i] - 10 * cos(2 * pai*z1[i]) + 10);

	double z2[DIM];
	for (int i = 0; i < DIM; i++)
		z2[i] = X[i];

	for (int i = 0; i < DIM; i++)
		y2 += (z2[i] - C22_O[i]) *(z2[i] - C22_O[i]);

	double ff1 = 0, ff2 = 1;
	for (int i = 0; i < DIM; i++) {
		ff1 += z2[i] * z2[i] / 4000;
		ff2 *= cos(z2[i] / sqrt(i + 1));
	}
	f2 = ff1 - ff2 + 1;
	//cout << "----f2=" << f2 << endl;


	double z3[DIM];
	for (int i = 0; i < DIM; i++) {
		z3[i] = X[i];
	}
	for (int i = 0; i < DIM; i++)
		y3 += (z3[i] - C23_O[i]) *(z3[i] - C23_O[i]);

	for (int i = 0; i < DIM; i++)
	{
		z3[i] = z3[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z3[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z3[i]) <= 500)
			f3 += z3[i] * sin(pow(fabs(z3[i]), 0.5));
		if (z3[i] > 500)
			f3 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z3[i] - 500)*(z3[i] - 500) / (10000 * DIM);
		if (z3[i] < -500)
			f3 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z3[i] + 500)*(z3[i] + 500) / (10000 * DIM);
	}
	f3 = 418.9829*DIM - f3;




	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

	w = w1 + w2 + w3;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
	return f + 2200;

#endif

#if Composition_function_3
	//F(X)=Rosenbrock+Ackley+Modifed Schwefel+Rastrigin
	double R[4] = { 10,20,30,40 };
	double T[4] = { 1,10,1,1 };
	double bias[4] = { 0,100,200,300 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;


	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;


	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w4 = 0;

	double w = 0;
	extern double C31_O[100];
	extern double C32_O[100];
	extern double C33_O[100];
	extern double C34_O[100];


	for (int i = 0; i < DIM; i++)
		y1 += (X[i] - C31_O[i]) * (X[i] - C31_O[i]);

	for (int i = 0; i < DIM - 1; i++)
		f1 += 100 * (X[i] * X[i] - X[i + 1]) * (X[i] * X[i] - X[i + 1]) + (X[i] - 1) * (X[i] - 1);

	//cout<<"f1="<<f1 <<endl;

	for (int i = 0; i < DIM; i++)
		y2 += (X[i] - C32_O[i]) * (X[i] - C32_O[i]);

	double fff1 = 0, fff2 = 0;
	for (int i = 0; i < DIM; i++) {
		fff1 += X[i] * X[i];
		fff2 += cos(2 * pai * X[i]);
	}
	f2 = -20 * exp(-0.2*sqrt((1.0 / DIM) * fff1)) - exp((1.0 / DIM) * fff2) + exp(1) + 20;
	//	cout << "-----f2=" << endl;

	double X1[DIM];
	for (int i = 0; i < DIM; i++)
		y3 += (X[i] - C33_O[i]) *(X[i] - C33_O[i]);
	for (int i = 0; i < DIM; i++)
	{
		X1[i] = X[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(X1[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(X1[i]) <= 500)
			f3 += X1[i] * sin(pow(fabs(X1[i]), 0.5));
		if (X1[i] > 500)
			f3 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (X1[i] - 500)*(X1[i] - 500) / (10000 * DIM);
		if (X1[i] < -500)
			f3 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (X1[i] + 500)*(X1[i] + 500) / (10000 * DIM);
	}
	f3 = 418.9829*DIM - f3;
	//	cout << "---f3=" << f3 << endl;


	for (int i = 0; i < DIM; i++)
		y4 += (X[i] - C34_O[i]) * (X[i] - C34_O[i]);
	for (int i = 0; i < DIM; i++)
		f4 = f4 + (X[i] * X[i] - 10 * cos(2 * pai*X[i]) + 10);
	//	cout << "----f4=" << f4 << endl;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w4 = 1 / pow(y4, 0.5)*exp(-y4 / (2 * DIM*R[3] * R[3]));

	w = w1 + w2 + w3 + w4;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	w4 = w4 / w;


	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]) + w4 * (T[3] * f4 + bias[3]);
	return f + 2300;
#endif

#if Composition_function_4
	//F(X)=Ackley + High Conditioned Elliptic + Griewank + Rastrigin
	double R[4] = { 10,20,30,40 };
	double T[4] = { 10,1E-6,10,1 };
	double bias[4] = { 0,100,200,300 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;


	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;


	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w4 = 0;

	double w = 0;
	extern double C41_O[100];
	extern double C42_O[100];
	extern double C43_O[100];
	extern double C44_O[100];


	for (int i = 0; i < DIM; i++)
		y1 += (X[i] - C41_O[i]) * (X[i] - C41_O[i]);
	//Ackley
	double fg1 = 0, fg2 = 0;
	for (int i = 0; i < DIM; i++)
	{
		fg1 += X[i] * X[i];
		fg2 += cos(2 * pai * X[i]);
	}
	f1 = -20 * exp(-0.2*sqrt((1.0 / DIM) * fg1)) - exp((1.0 / DIM) * fg2) + exp(1) + 20;



	for (int i = 0; i < DIM; i++)
		y2 += (X[i] - C42_O[i]) * (X[i] - C42_O[i]);
	//High Conditioned Elliptic function
	for (int i = 0; i < DIM; i++)
		f2 += pow(pow(10, 6), i / (DIM - 1))*X[i] * X[i];


	for (int i = 0; i < DIM; i++)
		y3 += (X[i] - C43_O[i]) *(X[i] - C43_O[i]);

	double f111 = 0, f222 = 1;
	for (int i = 0; i < DIM; i++)
	{
		f111 += X[i] * X[i] / 4000;
		f222 *= cos(X[i] / sqrt(i + 1));
	}
	f3 = f111 - f222 + 1;



	for (int i = 0; i < DIM; i++)
		y4 += (X[i] - C44_O[i]) * (X[i] - C44_O[i]);

	for (int i = 0; i < DIM; i++)
		f4 = f4 + (X[i] * X[i] - 10 * cos(2 * pai*X[i]) + 10);


	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w4 = 1 / pow(y4, 0.5)*exp(-y4 / (2 * DIM*R[3] * R[3]));

	w = w1 + w2 + w3 + w4;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	w4 = w4 / w;


	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]) + w4 * (T[3] * f4 + bias[3]);
	return f + 2400;
#endif

#if Composition_function_5
	//F(X)=Rastrigin + HappyCat + Ackley + Discus + Rosenbrock
	double f = 0;
	double R[5] = { 10,20,30,40,50 };
	double T[5] = { 10,1,10,1E-6,1 };
	double bias[5] = { 0,100,200,300,400 };

	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;
	double y5 = 0;

	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;


	extern double C51_O[100];
	extern double C52_O[100];
	extern double C53_O[100];
	extern double C54_O[100];
	extern double C55_O[100];

	for (int i = 0; i < DIM; i++)
	{
		y1 += (X[i] - C51_O[i])*(X[i] - C51_O[i]);
		y2 += (X[i] - C52_O[i])*(X[i] - C52_O[i]);
		y3 += (X[i] - C53_O[i])*(X[i] - C53_O[i]);
		y4 += (X[i] - C54_O[i])*(X[i] - C54_O[i]);
		y5 += (X[i] - C55_O[i])*(X[i] - C55_O[i]);

	}


	//Rastrigin
	for (int i = 0; i < DIM; i++)
		f1 = f1 + (X[i] * X[i] - 10 * cos(2 * pai*X[i]) + 10);

	//HappyCat
	double fd1 = 0.0;
	double fd2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < DIM; i++)
	{
		fd1 += X[i] * X[i];
		fd2 += X[i];
	}
	k = fabs(fd1 - DIM);
	f2 = pow(k, 0.25) + (0.5*fd1 + fd2) / DIM + 0.5;


	//Ackley
	double fa1 = 0, fa2 = 0;
	for (int i = 0; i < DIM; i++) {
		fa1 += X[i] * X[i];
		fa2 += cos(2 * pai * X[i]);
	}
	f3 = -20 * exp(-0.2*sqrt((1.0 / DIM) * fa1)) - exp((1.0 / DIM) * fa2) + exp(1) + 20;

	//Discus
	double ff = 0.0;
	for (int i = 1; i < DIM; i++)
		ff += X[i] * X[i];

	f4 = 1E6*X[0] * X[0] + ff;


	//Rosenbrock
	for (int i = 0; i < DIM - 1; i++)
		f5 += 100 * (X[i] * X[i] - X[i + 1]) * (X[i] * X[i] - X[i + 1]) + (X[i] - 1) * (X[i] - 1);

	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w4 = 0;
	double w5 = 0;
	double w = 0;

	w1 = 1 / sqrt(y1)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / sqrt(y2)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / sqrt(y3)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w4 = 1 / sqrt(y4)*exp(-y4 / (2 * DIM*R[3] * R[3]));
	w5 = 1 / sqrt(y5)*exp(-y5 / (2 * DIM*R[4] * R[4]));

	w = w1 + w2 + w3 + w4 + w5;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	w4 = w4 / w;
	w5 = w5 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]) + w4 * (T[3] * f4 + bias[3]) + w5 * (T[4] * f5 + bias[4]);

	return f + 2500;

#endif

#if Composition_function_6
	//F(x)=Expanded Schaffer's F6 + Modified Schwefel's f10 + Griewank + Rosenbrock + Rastrigin
	double R[5] = { 10,20,20,30,40 };
	double T[5] = { 1E-26,10,1E-6,10,5E-4 };
	double bias[5] = { 0,100,200,300,400 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;

	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;
	double y5 = 0;

	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w4 = 0;
	double w5 = 0;
	double w = 0;

	extern double C61_O[100];
	extern double C62_O[100];
	extern double C63_O[100];
	extern double C64_O[100];
	extern double C65_O[100];


	for (int i = 0; i < DIM; i++)
	{

		y1 += (X[i] - C61_O[i]) * (X[i] - C61_O[i]);
		y2 += (X[i] - C62_O[i]) * (X[i] - C62_O[i]);
		y3 += (X[i] - C63_O[i]) *(X[i] - C63_O[i]);
		y4 += (X[i] - C64_O[i]) *(X[i] - C64_O[i]);
		y5 += (X[i] - C65_O[i]) *(X[i] - C65_O[i]);

	}


	double g = 0.0;
	for (int i = 0; i < DIM - 1; i++)
	{
		double yy = X[i] * X[i] + X[i + 1] * X[i + 1];
		g = 0.5 + (sin(pow(yy, 0.5))*sin(pow(yy, 0.5)) - 0.5) / pow((1 + 0.001*yy), 2);
		f1 += g;
	}
	f1 += 0.5 + (sin(pow(X[DIM - 1] * X[DIM - 1] + X[0] * X[0], 0.5))*sin(pow(X[DIM - 1] * X[DIM - 1] + X[0] * X[0], 0.5)) - 0.5) / pow((1 + 0.001*(X[DIM - 1] * X[DIM - 1] + X[0] * X[0])), 2);
	//	cout << "----f1=" << f1 << endl;

	double X1[DIM];
	double ff = 0;
	for (int i = 0; i < DIM; i++)
	{
		X1[i] = X[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(X1[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(X1[i]) <= 500)
			ff += X1[i] * sin(pow(fabs(X1[i]), 0.5));
		else if (X1[i] > 500)
			ff += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (X1[i] - 500)*(X1[i] - 500) / (10000 * DIM);
		else if (X1[i] < -500)
			ff += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (X1[i] + 500)*(X1[i] + 500) / (10000 * DIM);
	}
	f2 = 418.9829*DIM - ff;
	//	cout << "----f2=" << f2 << endl;

	double ff1 = 0, ff2 = 1;
	for (int i = 0; i < DIM; i++)
	{
		ff1 += X[i] * X[i] / 4000;
		ff2 *= cos(X[i] / sqrt(i + 1));
	}
	f3 = ff1 - ff2 + 1;
	//	cout << "----f3=" << f3 << endl;


	for (int i = 0; i < DIM - 1; i++) {
		f4 += 100 * (X[i] * X[i] - X[i + 1]) * (X[i] * X[i] - X[i + 1]) + (X[i] - 1) * (X[i] - 1);
	}
	//	cout << "----f4=" << f4 << endl;



	for (int i = 0; i < DIM; i++)
		f5 = f5 + (X[i] * X[i] - 10 * cos(2 * pai*X[i]) + 10);
	//cout << "----f5=" << f5 << endl;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w4 = 1 / pow(y4, 0.5)*exp(-y4 / (2 * DIM*R[3] * R[3]));
	w5 = 1 / pow(y5, 0.5)*exp(-y5 / (2 * DIM*R[4] * R[4]));

	w = w1 + w2 + w3 + w4 + w5;
	//cout << "---w=" << w << endl;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	w4 = w4 / w;
	w5 = w5 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]) + w4 * (T[3] * f4 + bias[3]) + w5 * (T[4] * f5 + bias[4]);

	return f + 2600;
#endif

#if Composition_function_7
	//F(X)=HGBat+Rastrigin+Modified Schwefel's f10+Bent Cigar+High Conditioned Elliptic+Scaffer's f6
	int Dim = DIM;
	double R[6] = { 10,20,30,40,50,60 };
	double T[6] = { 10,10,2.5,1E-26,1E-6,5E-4 };
	double bias[6] = { 0,100,200,300,400,500 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;
	double f6 = 0;

	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;
	double y5 = 0;
	double y6 = 0;

	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w4 = 0;
	double w5 = 0;
	double w6 = 0;

	double w = 0;
	extern double C71_O[100];
	extern double C72_O[100];
	extern double C73_O[100];
	extern double C74_O[100];
	extern double C75_O[100];
	extern double C76_O[100];


	double z1[DIM];
	for (int i = 0; i < DIM; i++) {
		z1[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y1 += (z1[i] - C71_O[i]) * (z1[i] - C71_O[i]);

	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < DIM; i++)
	{
		ff1 += z1[i] * z1[i];
		ff2 += z1[i];
	}
	f1 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / DIM + 0.5;



	double z2[DIM];
	for (int i = 0; i < DIM; i++) {
		z2[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y2 += (z2[i] - C72_O[i]) *(z2[i] - C72_O[i]);
	for (int i = 0; i < Dim; i++)
		f2 = f2 + (z2[i] * z2[i] - 10 * cos(2 * pai*z2[i]) + 10);


	double z3[DIM];
	for (int i = 0; i < DIM; i++) {
		z3[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y3 += (z3[i] - C73_O[i]) *(z3[i] - C73_O[i]);

	for (int i = 0; i < DIM; i++)
	{
		z3[i] = z3[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z3[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z3[i]) <= 500)
			f3 += z3[i] * sin(pow(fabs(z3[i]), 0.5));
		if (z3[i] > 500)
			f3 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z3[i] - 500)*(z3[i] - 500) / (10000 * Dim);
		if (z3[i] < -500)
			f3 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z3[i] + 500)*(z3[i] + 500) / (10000 * Dim);
	}
	f3 = 418.9829*Dim - f3;

	double z4[DIM];
	for (int i = 0; i < DIM; i++) {
		z4[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y4 += (z4[i] - C74_O[i]) *(z4[i] - C74_O[i]);

	double ye1 = z4[0] * z4[0];
	double ye2 = 0.0;
	for (int i = 1; i < DIM; i++)
		ye2 += 1E6*z4[i] * z4[i];
	f4 = ye1 + ye2;


	double z5[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z5[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y5 += (z5[i] - C75_O[i]) * (z5[i] - C75_O[i]);

	for (int i = 0; i < DIM; i++)
		f5 += pow(pow(10, 6), i / (Dim - 1))*z5[i] * z5[i];


	double z6[DIM];
	for (int i = 0; i < DIM; i++)
	{
		z6[i] = X[i];
	}

	for (int i = 0; i < DIM; i++)
		y6 += (z6[i] - C76_O[i]) * (z6[i] - C76_O[i]);
	double g = 0.0;
	for (int i = 0; i < DIM - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z6[i] * z6[i] + z6[i + 1] * z6[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z6[i] * z6[i] + z6[i + 1] * z6[i + 1])), 2);
		f6 += g;
	}
	f6 += 0.5 + (pow(sin(pow(z6[Dim - 1] * z6[Dim - 1] + z6[0] * z6[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z6[Dim - 1] * z6[Dim - 1] + z6[0] * z6[0])), 2);


	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * Dim*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * Dim*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * Dim*R[2] * R[2]));
	w4 = 1 / pow(y4, 0.5)*exp(-y4 / (2 * Dim*R[3] * R[3]));
	w5 = 1 / pow(y5, 0.5)*exp(-y5 / (2 * Dim*R[4] * R[4]));
	w6 = 1 / pow(y6, 0.5)*exp(-y6 / (2 * Dim*R[5] * R[5]));
	w = w1 + w2 + w3 + w4 + w5 + w6;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	w4 = w4 / w;
	w5 = w5 / w;
	w6 = w6 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]) + w4 * (T[3] * f4 + bias[3]) + w5 * (T[4] * f5 + bias[4]) + w6 * (T[5] * f6 + bias[5]);
	return f + 2700;
#endif

#if Composition_function_8
	//F(X)=Ackley+Griewank+Disus+Rosenbrock+HappyCat+Expanded Schaffer's F6
	double R[6] = { 10,20,30,40,50,60 };
	double T[6] = { 10,10,1E-6,1,1,5E-4 };
	double bias[6] = { 0,100,200,300,400,500 };

	double f = 0;
	double f1 = 0;
	double f2 = 0;
	double f3 = 0;
	double f4 = 0;
	double f5 = 0;
	double f6 = 0;

	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;
	double y5 = 0;
	double y6 = 0;


	extern double C81_O[100];
	extern double C82_O[100];
	extern double C83_O[100];
	extern double C84_O[100];
	extern double C85_O[100];
	extern double C86_O[100];
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w4 = 0;
	double w5 = 0;
	double w6 = 0;

	double w = 0;
	for (int i = 0; i < DIM; i++)
	{
		//cout << X[i] << ",";
		y1 += (X[i] - C81_O[i]) *(X[i] - C81_O[i]);
		y2 += (X[i] - C82_O[i]) * (X[i] - C82_O[i]);
		y3 += (X[i] - C83_O[i]) * (X[i] - C83_O[i]);
		y4 += (X[i] - C84_O[i]) *(X[i] - C84_O[i]);
		y5 += (X[i] - C85_O[i])* (X[i] - C85_O[i]);
		y6 += (X[i] - C86_O[i])* (X[i] - C86_O[i]);
	}


	//Ackley
	double ff1 = 0, ff2 = 0;
	for (int i = 0; i < DIM; i++)
	{
		ff1 += X[i] * X[i];
		ff2 += cos(2 * pai * X[i]);
	}
	f1 = -20 * exp(-0.2*sqrt((1.0 / DIM) * ff1)) - exp((1.0 / DIM) * ff2) + exp(1) + 20;
	//cout << "---f1=" << f1 << endl;


	//Griewannk
	double fr1 = 0, fr2 = 1;
	for (int i = 0; i < DIM; i++)
	{
		fr1 += X[i] * X[i] / 4000;
		fr2 *= cos(X[i] / sqrt(i + 1));
	}
	f2 = fr1 - fr2 + 1;
	//cout << "---f2="<<f2<<endl;


	//Discus
	double fd1 = 0.0;
	for (int i = 1; i < DIM; i++)
		fd1 += X[i] * X[i];
	f3 = 1E6*X[0] * X[0] + fd1;
	//cout << "-----f3=" << f3 << endl;


	//Rosenbrock
	for (int i = 0; i < DIM - 1; i++)
		f4 += 100 * (X[i] * X[i] - X[i + 1]) * (X[i] * X[i] - X[i + 1]) + (X[i] - 1) * (X[i] - 1);
	//cout << "---f4=" << f4 << endl;

	//HappyCat
	double fs1 = 0.0;
	double fs2 = 0.0;
	double k = 0.0;
	for (int i = 0; i < DIM; i++)
	{
		fs1 += X[i] * X[i];
		fs2 += X[i];
	}
	k = fabs(fs1 - DIM);
	if (k == 0)
		f5 = (0.5*fs1 + fs2) / DIM + 0.5;
	else
		f5 = pow(k, 0.25) + (0.5*fs1 + fs2) / DIM + 0.5;
	//cout << "----f5=" << f5 << endl;

	//Expanded Scaffer's F6
	double g = 0.0;
	for (int i = 0; i < DIM - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(X[i] * X[i] + X[i + 1] * X[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(X[i] * X[i] + X[i + 1] * X[i + 1])), 2);
		f6 += g;
	}
	f6 += 0.5 + (pow(sin(pow(X[DIM - 1] * X[DIM - 1] + X[0] * X[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(X[DIM - 1] * X[DIM - 1] + X[0] * X[0])), 2);
	//cout <<"----f6="<< f6 << endl;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w4 = 1 / pow(y4, 0.5)*exp(-y4 / (2 * DIM*R[3] * R[3]));
	w5 = 1 / pow(y5, 0.5)*exp(-y5 / (2 * DIM*R[4] * R[4]));
	w6 = 1 / pow(y6, 0.5)*exp(-y6 / (2 * DIM*R[5] * R[5]));
	w = w1 + w2 + w3 + w4 + w5 + w6;

	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	w4 = w4 / w;
	w5 = w5 / w;
	w6 = w6 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]) + w4 * (T[3] * f4 + bias[3]) + w5 * (T[4] * f5 + bias[4]) + w6 * (T[5] * f6 + bias[5]);
	//cout << f<<endl;
	return f + 2800;
#endif

#if Composition_function_9
	//F(X)=Hybrid function 5+hybrid function 8+hybrid function 9
	double R[3] = { 10,30,50 };
	double T[3] = { 1,1,1 };
	double bias[3] = { 0,100,200 };

	double result = 0;

	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	extern double C91_O[100];
	extern double C92_O[100];
	extern double C93_O[100];
	for (int i = 0; i < DIM; i++)
	{
		y1 += (X[i] - C91_O[i])*(X[i] - C91_O[i]);
		y2 += (X[i] - C92_O[i])*(X[i] - C92_O[i]);
		y3 += (X[i] - C93_O[i])*(X[i] - C93_O[i]);
	}


#if DIM==10
	//hybrid function 5
	extern double Hy5_O[100];
	extern double Hy5_M10[10][10];
	int Dim111;
	int Dim211;
	int Dim311;
	int Dim411;
	Dim111 = DIM / 10 * 2;
	Dim211 = DIM / 10 * 2;
	Dim311 = DIM / 10 * 3;
	Dim411 = DIM / 10 * 3;

	double fh51 = 0;
	double fh52 = 0;
	double fh53 = 0;
	double fh54 = 0;
	double fh5 = 0;

	//Bent Cigar function
	double Hy51_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy51_M2[i][j] = Hy5_M10[i][j];
	}
	double z1[2];
	for (int i = 0; i < Dim111; i++)
	{
		z1[i] = X[i] - Hy5_O[i];

	}
	mul_matrix02(z1, Hy51_M2);
	double yy1 = z1[0] * z1[0];
	double yy2 = 0.0;
	for (int i = 1; i < Dim111; i++)
		yy2 += 1E6*z1[i] * z1[i];
	fh51 = yy1 + yy2;


	//HGBat function
	double Hy52_M2[2][2];
	for (int i = 2; i < 4; i++)
	{
		for (int j = 2; j < 4; j++)
			Hy52_M2[i - 2][j - 2] = Hy5_M10[i][j];
	}

	double z2[2];
	for (int i = 0; i < 2; i++)
	{
		z2[i] = X[i + 2] - Hy5_O[i + 2];
		//	y1 += z2[i] * z2[i];
	}
	mul_matrix02(z2, Hy52_M2);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim211; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	fh52 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim211 + 0.5;


	//Rastrigin function 
	double Hy53_M3[3][3];
	for (int i = 4; i < 7; i++)
	{
		for (int j = 4; j < 7; j++)
			Hy53_M3[i - 4][j - 4] = Hy5_M10[i][j];
	}

	double z3[3];
	for (int i = 0; i < Dim311; i++)
	{
		z3[i] = X[i + 4] - Hy5_O[i + 4];
		//	y1 += z3[i] * z3[i];
	}
	mul_matrix03(z3, Hy53_M3);
	for (int i = 0; i < Dim311; i++)
		fh53 = fh53 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M3[3][3];
	for (int i = 7; i < 10; i++)
	{
		for (int j = 7; j < 10; j++)
			Hy54_M3[i - 7][j - 7] = Hy5_M10[i][j];
	}

	double z4[3];
	for (int i = 0; i < Dim411; i++)
	{
		z4[i] = X[i + 7] - Hy5_O[i + 7];
	}
	mul_matrix03(z4, Hy54_M3);

	for (int i = 0; i < Dim411 - 1; i++)
		fh54 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);
	fh5 = fh51 + fh52 + fh53 + fh54;
	//cout << "fh5=" << fh5 << endl;

	//f2=hybrid function 8
	extern double Hy8_O[100];
	extern double Hy8_M10[10][10];
	double fh81 = 0;
	double fh82 = 0;
	double fh83 = 0;
	double fh84 = 0;
	double fh85 = 0;
	double fh8 = 0;

	double Hy81_M2[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Hy81_M2[i][j] = Hy8_M10[i][j];
	double Hy82_M2[2][2];
	for (int i = 2; i < 4; i++)
		for (int j = 2; j < 4; j++)
			Hy82_M2[i - 2][j - 2] = Hy8_M10[i][j];
	double Hy83_M2[2][2];
	for (int i = 4; i < 6; i++)
		for (int j = 4; j < 6; j++)
			Hy83_M2[i - 4][j - 4] = Hy8_M10[i][j];
	double Hy84_M2[2][2];
	for (int i = 6; i < 8; i++)
		for (int j = 6; j < 8; j++)
			Hy84_M2[i - 6][j - 6] = Hy8_M10[i][j];
	double Hy85_M2[2][2];
	for (int i = 8; i < 10; i++)
		for (int j = 8; j < 10; j++)
			Hy85_M2[i - 8][j - 8] = Hy8_M10[i][j];

	double s1[2];
	for (int i = 0; i < 2; i++)
	{
		s1[i] = X[i] - Hy8_O[i];
		//	y2 += s1[i] * s1[i];
	}
	mul_matrix02(s1, Hy81_M2);

	for (int i = 0; i < 2; i++)
		fh81 += pow(pow(10, 6), i / (2 - 1))*s1[i] * s1[i];


	double s2[2];
	for (int i = 2; i < 4; i++)
	{
		s2[i - 2] = X[i] - Hy8_O[i];
	}
	mul_matrix02(s2, Hy82_M2);
	double ffx1 = 0;
	double ffx2 = 0;
	for (int i = 0; i < 2; i++)
	{
		ffx1 += s2[i] * s2[i];
		ffx2 += cos(2 * pai * s2[i]);
	}
	fh82 = -20 * exp(-0.2*sqrt((1.0 / 2) * ffx1)) - exp((1.0 / 2) * ffx2) + exp(1) + 20;


	double s3[2];
	for (int i = 4; i < 6; i++)
	{
		s3[i - 4] = X[i] - Hy8_O[i];
		//y2 += s3[i - 4] * s3[i - 4];
	}
	mul_matrix02(s3, Hy83_M2);
	for (int i = 0; i < 2; i++)
		fh83 = fh83 + (s3[i] * s3[i] - 10 * cos(2 * pai*s3[i]) + 10);

	double s4[2];
	for (int i = 6; i < 8; i++)
	{
		s4[i - 6] = X[i] - Hy8_O[i];
	}
	mul_matrix02(s4, Hy84_M2);
	double ffd1 = 0.0;
	double ffd2 = 0.0;
	for (int i = 0; i < 2; i++)
	{
		ffd1 += s4[i] * s4[i];
		ffd2 += s4[i];
	}
	fh84 = pow(fabs(ffd1*ffd1 - ffd2 * ffd2), 0.5) + (0.5*ffd1 + ffd2) / 2 + 0.5;

	//Discus
	double s5[2];
	for (int i = 8; i < 10; i++)
	{
		s5[i - 8] = X[i] - Hy8_O[i];
		//y2 += s5[i - 8] * s5[i - 8];
	}
	mul_matrix02(s5, Hy85_M2);
	double fz1 = 0.0;
	for (int i = 1; i < 2; i++)
		fz1 += s5[i] * s5[i];
	fh85 = 1E6 * s5[0] * s5[0] + fz1;

	fh8 = fh81 + fh82 + fh83 + fh84 + fh85;
	//	cout << "---fh8=" << fh8 << endl;

	extern double Hy9_O[100];
	extern double Hy9_M10[10][10];
	//f3=hybrid function 9
	int Dim11 = DIM / 10 * 2;
	int Dim21 = DIM / 10 * 2;
	int Dim31 = DIM / 10 * 2;
	int Dim41 = DIM / 10 * 2;
	int Dim51 = DIM / 10 * 2;
	double fh91 = 0;
	double fh92 = 0;
	double fh93 = 0;
	double fh94 = 0;
	double fh95 = 0;
	double fh9 = 0;

	double Hy91_M2[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Hy91_M2[i][j] = Hy9_M10[i][j];
	double Hy92_M2[2][2];
	for (int i = 2; i < 4; i++)
		for (int j = 2; j < 4; j++)
			Hy92_M2[i - 2][j - 2] = Hy9_M10[i][j];
	double Hy93_M2[2][2];
	for (int i = 4; i < 6; i++)
		for (int j = 4; j < 6; j++)
			Hy93_M2[i - 4][j - 4] = Hy9_M10[i][j];
	double Hy94_M2[2][2];
	for (int i = 6; i < 8; i++)
		for (int j = 6; j < 8; j++)
			Hy94_M2[i - 6][j - 6] = Hy9_M10[i][j];
	double Hy95_M2[2][2];
	for (int i = 8; i < 10; i++)
		for (int j = 8; j < 10; j++)
			Hy95_M2[i - 8][j - 8] = Hy9_M10[i][j];

	double p1[2];
	for (int i = 0; i < 2; i++)
	{
		p1[i] = X[i] - Hy9_O[i];
	}
	mul_matrix02(p1, Hy91_M2);
	double yu1 = p1[0] * p1[0];
	double yu2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yu2 += 1E6*p1[i] * p1[i];
	fh91 = yu1 + yu2;


	double p2[2];
	for (int i = 2; i < 4; i++)
	{
		p2[i - 2] = X[i] - Hy9_O[i];

	}
	mul_matrix02(p2, Hy92_M2);
	for (int i = 0; i < 2; i++)
		fh92 = fh92 + (p2[i] * p2[i] - 10 * cos(2 * pai*p2[i]) + 10);

	double p3[2];
	for (int i = 4; i < 6; i++)
	{
		p3[i - 4] = X[i] - Hy9_O[i];

	}
	mul_matrix02(p3, Hy93_M2);
	//Griewank's plus Rosenbrock
	double temp[2];
	for (int i = 0; i < Dim31 - 2; i++)
		temp[i] = 100 * (p3[i] * p3[i] - p3[i + 1])*(p3[i] * p3[i] - p3[i + 1]) + (p3[i] - 1)*(p3[i] - 1) + 100 * (p3[i + 1] * p3[i + 1] - p3[i + 2])*(p3[i + 1] * p3[i + 1] - p3[i + 2]) + (p3[i + 1] - 1)*(p3[i + 1] - 1);
	temp[Dim31 - 2] = 100 * (p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1])*(p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1]) + (p3[Dim31 - 2] - 1)*(p3[Dim31 - 2] - 1);
	temp[Dim31 - 1] = 100 * (p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0])*(p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0]) + (p3[Dim31 - 1] - 1)*(p3[Dim31 - 1] - 1);

	double u0 = 2.5;//u0
	double d = 1;//d
	double s = 1 - 1 / (2 * pow(Dim31 + 20, 0.5) - 8.2);//s
	double u1 = -pow((u0*u0 - d) / s, 0.5);//u1
	double sign[2];//signº¯Êý
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[2][2];//¶Ô½Ç¾ØÕó
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim31 - 1)));
	}
	double zn[2];//z
	double yy[2];//y  y=10(x-o)/100
	double xz[2];//x»¡
	for (int i = 0; i < Dim31; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 4]);

	for (int i = 0; i < Dim31; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim31; i++)
		zn[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim31; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim31;
	if (F1 < F2)
		fh93 = F1;
	else fh93 = F2;

	for (int i = 0; i < Dim31; i++)
		ff = ff + cos(2 * pai*zn[i]);
	fh93 = fh93 + 10 * (Dim31 - ff);


	double p4[2];
	for (int i = 6; i < 8; i++)
	{
		p4[i - 6] = X[i] - Hy9_O[i];
		//y3 += p4[i - 6] * p4[i - 6];
	}
	mul_matrix02(p4, Hy94_M2);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yyr2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yyr2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim41; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(p4[i] + 0.5));
		ff1 += yy1;
	}
	fh94 = ff1 - Dim41 * yyr2;



	double p5[2];
	for (int i = 8; i < 10; i++)
	{
		p5[i - 8] = X[i] - Hy9_O[i];
		//y3 += p5[i - 8] * p5[i - 8];
	}
	mul_matrix02(p5, Hy95_M2);
	double g = 0.0;
	for (int i = 0; i < Dim51 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(p5[i] * p5[i] + p5[i + 1] * p5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[i] * p5[i] + p5[i + 1] * p5[i + 1])), 2);
		fh95 += g;
	}
	fh95 += 0.5 + (pow(sin(pow(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0])), 2);

	fh9 = fh91 + fh92 + fh93 + fh94 + fh95;


	//	cout << "---fh9=" << fh9 << endl;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	result = w1 * (T[0] * fh5 + bias[0]) + w2 * (T[1] * fh8 + bias[1]) + w3 * (T[2] * fh9 + bias[2]);
#endif 



#if DIM==30
	//hybrid function 5
	extern double Hy5_O[100];
	extern double Hy5_M30[30][30];
	int Dim111;
	int Dim211;
	int Dim311;
	int Dim411;
	Dim111 = DIM / 10 * 2;
	Dim211 = DIM / 10 * 2;
	Dim311 = DIM / 10 * 3;
	Dim411 = DIM / 10 * 3;

	double fh51 = 0;
	double fh52 = 0;
	double fh53 = 0;
	double fh54 = 0;
	double fh5 = 0;

	//Bent Cigar function
	double Hy51_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy51_M6[i][j] = Hy5_M30[i][j];
	}
	double z1[6];
	for (int i = 0; i < Dim111; i++)
	{
		z1[i] = X[i] - Hy5_O[i];

	}
	mul_matrix06(z1, Hy51_M6);
	double yy1 = z1[0] * z1[0];
	double yy2 = 0.0;
	for (int i = 1; i < Dim111; i++)
		yy2 += 1E6*z1[i] * z1[i];
	fh51 = yy1 + yy2;


	//HGBat function
	double Hy52_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy52_M6[i - 6][j - 6] = Hy5_M30[i][j];
	}

	double z2[6];
	for (int i = 0; i < 6; i++)
	{
		z2[i] = X[i + 6] - Hy5_O[i + 6];
		//	y1 += z2[i] * z2[i];
	}
	mul_matrix06(z2, Hy52_M6);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim211; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	fh52 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim211 + 0.5;


	//Rastrigin function 
	double Hy53_M9[9][9];
	for (int i = 12; i < 21; i++)
	{
		for (int j = 12; j < 21; j++)
			Hy53_M9[i - 12][j - 12] = Hy5_M30[i][j];
	}

	double z3[9];
	for (int i = 0; i < Dim311; i++)
	{
		z3[i] = X[i + 12] - Hy5_O[i + 12];
		//	y1 += z3[i] * z3[i];
	}
	mul_matrix09(z3, Hy53_M9);
	for (int i = 0; i < Dim311; i++)
		fh53 = fh53 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M9[9][9];
	for (int i = 21; i < 30; i++)
	{
		for (int j = 21; j < 30; j++)
			Hy54_M9[i - 21][j - 21] = Hy5_M30[i][j];
	}

	double z4[9];
	for (int i = 0; i < Dim411; i++)
	{
		z4[i] = X[i + 21] - Hy5_O[i + 21];
	}
	mul_matrix09(z4, Hy54_M9);

	for (int i = 0; i < Dim411 - 1; i++)
		fh54 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);
	fh5 = fh51 + fh52 + fh53 + fh54;
	//cout << "fh5=" << fh5 << endl;

	//f2=hybrid function 8
	extern double Hy8_O[100];
	extern double Hy8_M30[30][30];
	double fh81 = 0;
	double fh82 = 0;
	double fh83 = 0;
	double fh84 = 0;
	double fh85 = 0;
	double fh8 = 0;

	double Hy81_M6[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Hy81_M6[i][j] = Hy8_M30[i][j];
	double Hy82_M6[6][6];
	for (int i = 6; i < 12; i++)
		for (int j = 6; j < 12; j++)
			Hy82_M6[i - 6][j - 6] = Hy8_M30[i][j];
	double Hy83_M6[6][6];
	for (int i = 12; i < 18; i++)
		for (int j = 12; j < 18; j++)
			Hy83_M6[i - 12][j - 12] = Hy8_M30[i][j];
	double Hy84_M6[6][6];
	for (int i = 18; i < 24; i++)
		for (int j = 18; j < 24; j++)
			Hy84_M6[i - 18][j - 18] = Hy8_M30[i][j];
	double Hy85_M6[6][6];
	for (int i = 24; i < 30; i++)
		for (int j = 24; j < 30; j++)
			Hy85_M6[i - 24][j - 24] = Hy8_M30[i][j];

	double s1[6];
	for (int i = 0; i < 6; i++)
	{
		s1[i] = X[i] - Hy8_O[i];
		//	y2 += s1[i] * s1[i];
	}
	mul_matrix06(s1, Hy81_M6);

	for (int i = 0; i < 6; i++)
		fh81 += pow(pow(10, 6), i / (6 - 1))*s1[i] * s1[i];


	double s2[6];
	for (int i = 6; i < 12; i++)
	{
		s2[i - 6] = X[i] - Hy8_O[i];
	}
	mul_matrix06(s2, Hy82_M6);
	double ffx1 = 0;
	double ffx2 = 0;
	for (int i = 0; i < 6; i++)
	{
		ffx1 += s2[i] * s2[i];
		ffx2 += cos(2 * pai * s2[i]);
	}
	fh82 = -20 * exp(-0.2*sqrt((1.0 / 6) * ffx1)) - exp((1.0 / 6) * ffx2) + exp(1) + 20;


	double s3[6];
	for (int i = 12; i < 18; i++)
	{
		s3[i - 12] = X[i] - Hy8_O[i];
		//y2 += s3[i - 12] * s3[i - 12];
	}
	mul_matrix06(s3, Hy83_M6);
	for (int i = 0; i < 6; i++)
		fh83 = fh83 + (s3[i] * s3[i] - 10 * cos(2 * pai*s3[i]) + 10);

	double s4[6];
	for (int i = 18; i < 24; i++)
	{
		s4[i - 18] = X[i] - Hy8_O[i];
	}
	mul_matrix06(s4, Hy84_M6);
	double ffd1 = 0.0;
	double ffd2 = 0.0;
	for (int i = 0; i < 6; i++)
	{
		ffd1 += s4[i] * s4[i];
		ffd2 += s4[i];
	}
	fh84 = pow(fabs(ffd1*ffd1 - ffd2 * ffd2), 0.5) + (0.5*ffd1 + ffd2) / 6 + 0.5;

	//Discus
	double s5[6];
	for (int i = 24; i < 30; i++)
	{
		s5[i - 24] = X[i] - Hy8_O[i];
		//y2 += s5[i - 24] * s5[i - 24];
	}
	mul_matrix06(s5, Hy85_M6);
	double fz1 = 0.0;
	for (int i = 1; i < 6; i++)
		fz1 += s5[i] * s5[i];
	fh85 = 1E6 * s5[0] * s5[0] + fz1;

	fh8 = fh81 + fh82 + fh83 + fh84 + fh85;
	//	cout << "---fh8=" << fh8 << endl;

	extern double Hy9_O[100];
	extern double Hy9_M30[30][30];
	//f3=hybrid function 9
	int Dim11 = DIM / 10 * 2;
	int Dim21 = DIM / 10 * 2;
	int Dim31 = DIM / 10 * 2;
	int Dim41 = DIM / 10 * 2;
	int Dim51 = DIM / 10 * 2;
	double fh91 = 0;
	double fh92 = 0;
	double fh93 = 0;
	double fh94 = 0;
	double fh95 = 0;
	double fh9 = 0;

	double Hy91_M6[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Hy91_M6[i][j] = Hy9_M30[i][j];
	double Hy92_M6[6][6];
	for (int i = 6; i < 12; i++)
		for (int j = 6; j < 12; j++)
			Hy92_M6[i - 6][j - 6] = Hy9_M30[i][j];
	double Hy93_M6[6][6];
	for (int i = 12; i < 18; i++)
		for (int j = 12; j < 18; j++)
			Hy93_M6[i - 12][j - 12] = Hy9_M30[i][j];
	double Hy94_M6[6][6];
	for (int i = 18; i < 24; i++)
		for (int j = 18; j < 24; j++)
			Hy94_M6[i - 18][j - 18] = Hy9_M30[i][j];
	double Hy95_M6[6][6];
	for (int i = 24; i < 30; i++)
		for (int j = 24; j < 30; j++)
			Hy95_M6[i - 24][j - 24] = Hy9_M30[i][j];

	double p1[6];
	for (int i = 0; i < 6; i++)
	{
		p1[i] = X[i] - Hy9_O[i];
	}
	mul_matrix06(p1, Hy91_M6);
	double yu1 = p1[0] * p1[0];
	double yu2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yu2 += 1E6*p1[i] * p1[i];
	fh91 = yu1 + yu2;


	double p2[6];
	for (int i = 6; i < 12; i++)
	{
		p2[i - 6] = X[i] - Hy9_O[i];

	}
	mul_matrix06(p2, Hy92_M6);
	for (int i = 0; i < 6; i++)
		fh92 = fh92 + (p2[i] * p2[i] - 10 * cos(2 * pai*p2[i]) + 10);

	double p3[6];
	for (int i = 12; i < 18; i++)
	{
		p3[i - 12] = X[i] - Hy9_O[i];

	}
	mul_matrix06(p3, Hy93_M6);
	//Griewank's plus Rosenbrock
	double temp[6];
	for (int i = 0; i < Dim31 - 2; i++)
		temp[i] = 100 * (p3[i] * p3[i] - p3[i + 1])*(p3[i] * p3[i] - p3[i + 1]) + (p3[i] - 1)*(p3[i] - 1) + 100 * (p3[i + 1] * p3[i + 1] - p3[i + 2])*(p3[i + 1] * p3[i + 1] - p3[i + 2]) + (p3[i + 1] - 1)*(p3[i + 1] - 1);
	temp[Dim31 - 2] = 100 * (p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1])*(p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1]) + (p3[Dim31 - 2] - 1)*(p3[Dim31 - 2] - 1);
	temp[Dim31 - 1] = 100 * (p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0])*(p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0]) + (p3[Dim31 - 1] - 1)*(p3[Dim31 - 1] - 1);

	double u0 = 2.5;//u0
	double d = 1;//d
	double s = 1 - 1 / (2 * pow(Dim31 + 20, 0.5) - 8.2);//s
	double u1 = -pow((u0*u0 - d) / s, 0.5);//u1
	double sign[6];//signº¯Êý
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[6][6];//¶Ô½Ç¾ØÕó
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim31 - 1)));
	}
	double zn[6];//z
	double yy[6];//y  y=10(x-o)/100
	double xz[6];//x»¡
	for (int i = 0; i < Dim31; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 12]);

	for (int i = 0; i < Dim31; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim31; i++)
		zn[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim31; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim31;
	if (F1 < F2)
		fh93 = F1;
	else fh93 = F2;

	for (int i = 0; i < Dim31; i++)
		ff = ff + cos(2 * pai*zn[i]);
	fh93 = fh93 + 10 * (Dim31 - ff);


	double p4[6];
	for (int i = 18; i < 24; i++)
	{
		p4[i - 18] = X[i] - Hy9_O[i];
		//y3 += p4[i - 18] * p4[i - 18];
	}
	mul_matrix06(p4, Hy94_M6);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yyr2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yyr2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim41; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(p4[i] + 0.5));
		ff1 += yy1;
	}
	fh94 = ff1 - Dim41 * yyr2;



	double p5[6];
	for (int i = 24; i < 30; i++)
	{
		p5[i - 24] = X[i] - Hy9_O[i];
		//y3 += p5[i - 24] * p5[i - 24];
	}
	mul_matrix06(p5, Hy95_M6);
	double g = 0.0;
	for (int i = 0; i < Dim51 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(p5[i] * p5[i] + p5[i + 1] * p5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[i] * p5[i] + p5[i + 1] * p5[i + 1])), 2);
		fh95 += g;
	}
	fh95 += 0.5 + (pow(sin(pow(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0])), 2);

	fh9 = fh91 + fh92 + fh93 + fh94 + fh95;


	//	cout << "---fh9=" << fh9 << endl;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	result = w1 * (T[0] * fh5 + bias[0]) + w2 * (T[1] * fh8 + bias[1]) + w3 * (T[2] * fh9 + bias[2]);
#endif 

#if DIM==50
	//hybrid function 5
	extern double Hy5_O[100];
	extern double Hy5_M50[50][50];
	int Dim111;
	int Dim211;
	int Dim311;
	int Dim411;
	Dim111 = DIM / 10 * 2;
	Dim211 = DIM / 10 * 2;
	Dim311 = DIM / 10 * 3;
	Dim411 = DIM / 10 * 3;

	double fh51 = 0;
	double fh52 = 0;
	double fh53 = 0;
	double fh54 = 0;
	double fh5 = 0;

	//Bent Cigar function
	double Hy51_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy51_M10[i][j] = Hy5_M50[i][j];
	}
	double z1[10];
	for (int i = 0; i < Dim111; i++)
	{
		z1[i] = X[i] - Hy5_O[i];

	}
	mul_matrix10(z1, Hy51_M10);
	double yy1 = z1[0] * z1[0];
	double yy2 = 0.0;
	for (int i = 1; i < Dim111; i++)
		yy2 += 1E6*z1[i] * z1[i];
	fh51 = yy1 + yy2;


	//HGBat function
	double Hy52_M10[10][10];
	for (int i = 10; i < 20; i++)
	{
		for (int j = 10; j < 20; j++)
			Hy52_M10[i - 10][j - 10] = Hy5_M50[i][j];
	}

	double z2[10];
	for (int i = 0; i < 10; i++)
	{
		z2[i] = X[i + 10] - Hy5_O[i + 10];
		//	y1 += z2[i] * z2[i];
	}
	mul_matrix10(z2, Hy52_M10);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim211; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	fh52 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim211 + 0.5;


	//Rastrigin function 
	double Hy53_M15[15][15];
	for (int i = 20; i < 35; i++)
	{
		for (int j = 20; j < 35; j++)
			Hy53_M15[i - 20][j - 20] = Hy5_M50[i][j];
	}

	double z3[15];
	for (int i = 0; i < Dim311; i++)
	{
		z3[i] = X[i + 20] - Hy5_O[i + 20];
		//	y1 += z3[i] * z3[i];
	}
	mul_matrix15(z3, Hy53_M15);
	for (int i = 0; i < Dim311; i++)
		fh53 = fh53 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M15[15][15];
	for (int i = 35; i < 50; i++)
	{
		for (int j = 35; j < 50; j++)
			Hy54_M15[i - 35][j - 35] = Hy5_M50[i][j];
	}

	double z4[15];
	for (int i = 0; i < Dim411; i++)
	{
		z4[i] = X[i + 35] - Hy5_O[i + 35];
	}
	mul_matrix15(z4, Hy54_M15);

	for (int i = 0; i < Dim411 - 1; i++)
		fh54 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);
	fh5 = fh51 + fh52 + fh53 + fh54;
	//cout << "fh5=" << fh5 << endl;

	//f2=hybrid function 8
	extern double Hy8_O[100];
	extern double Hy8_M50[50][50];
	double fh81 = 0;
	double fh82 = 0;
	double fh83 = 0;
	double fh84 = 0;
	double fh85 = 0;
	double fh8 = 0;

	double Hy81_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy81_M10[i][j] = Hy8_M50[i][j];
	double Hy82_M10[10][10];
	for (int i = 10; i < 20; i++)
		for (int j = 10; j < 20; j++)
			Hy82_M10[i - 10][j - 10] = Hy8_M50[i][j];
	double Hy83_M10[10][10];
	for (int i = 20; i < 30; i++)
		for (int j = 20; j < 30; j++)
			Hy83_M10[i - 20][j - 20] = Hy8_M50[i][j];
	double Hy84_M10[10][10];
	for (int i = 30; i < 40; i++)
		for (int j = 30; j < 40; j++)
			Hy84_M10[i - 30][j - 30] = Hy8_M50[i][j];
	double Hy85_M10[10][10];
	for (int i = 40; i < 50; i++)
		for (int j = 40; j < 50; j++)
			Hy85_M10[i - 40][j - 40] = Hy8_M50[i][j];

	double s1[10];
	for (int i = 0; i < 10; i++)
	{
		s1[i] = X[i] - Hy8_O[i];
		//	y2 += s1[i] * s1[i];
	}
	mul_matrix10(s1, Hy81_M10);

	for (int i = 0; i < 10; i++)
		fh81 += pow(pow(10, 6), i / (10 - 1))*s1[i] * s1[i];


	double s2[10];
	for (int i = 10; i < 20; i++)
	{
		s2[i - 10] = X[i] - Hy8_O[i];
	}
	mul_matrix10(s2, Hy82_M10);
	double ffx1 = 0;
	double ffx2 = 0;
	for (int i = 0; i < 10; i++)
	{
		ffx1 += s2[i] * s2[i];
		ffx2 += cos(2 * pai * s2[i]);
	}
	fh82 = -20 * exp(-0.2*sqrt((1.0 / 10) * ffx1)) - exp((1.0 / 10) * ffx2) + exp(1) + 20;


	double s3[10];
	for (int i = 20; i < 30; i++)
	{
		s3[i - 20] = X[i] - Hy8_O[i];
		//y2 += s3[i - 20] * s3[i - 20];
	}
	mul_matrix10(s3, Hy83_M10);
	for (int i = 0; i < 10; i++)
		fh83 = fh83 + (s3[i] * s3[i] - 10 * cos(2 * pai*s3[i]) + 10);

	double s4[10];
	for (int i = 30; i < 40; i++)
	{
		s4[i - 30] = X[i] - Hy8_O[i];
	}
	mul_matrix10(s4, Hy84_M10);
	double ffd1 = 0.0;
	double ffd2 = 0.0;
	for (int i = 0; i < 10; i++)
	{
		ffd1 += s4[i] * s4[i];
		ffd2 += s4[i];
	}
	fh84 = pow(fabs(ffd1*ffd1 - ffd2 * ffd2), 0.5) + (0.5*ffd1 + ffd2) / 2 + 0.5;

	//Discus
	double s5[10];
	for (int i = 40; i < 50; i++)
	{
		s5[i - 40] = X[i] - Hy8_O[i];
		//y2 += s5[i - 40] * s5[i - 40];
	}
	mul_matrix10(s5, Hy85_M10);
	double fz1 = 0.0;
	for (int i = 1; i < 10; i++)
		fz1 += s5[i] * s5[i];
	fh85 = 1E6 * s5[0] * s5[0] + fz1;

	fh8 = fh81 + fh82 + fh83 + fh84 + fh85;
	//	cout << "---fh8=" << fh8 << endl;

	extern double Hy9_O[100];
	extern double Hy9_M50[50][50];
	//f3=hybrid function 9
	int Dim11 = DIM / 10 * 2;
	int Dim21 = DIM / 10 * 2;
	int Dim31 = DIM / 10 * 2;
	int Dim41 = DIM / 10 * 2;
	int Dim51 = DIM / 10 * 2;
	double fh91 = 0;
	double fh92 = 0;
	double fh93 = 0;
	double fh94 = 0;
	double fh95 = 0;
	double fh9 = 0;

	double Hy91_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy91_M10[i][j] = Hy9_M50[i][j];
	double Hy92_M10[10][10];
	for (int i = 10; i < 20; i++)
		for (int j = 10; j < 20; j++)
			Hy92_M10[i - 10][j - 10] = Hy9_M50[i][j];
	double Hy93_M10[10][10];
	for (int i = 20; i < 30; i++)
		for (int j = 20; j < 30; j++)
			Hy93_M10[i - 20][j - 20] = Hy9_M50[i][j];
	double Hy94_M10[10][10];
	for (int i = 30; i < 40; i++)
		for (int j = 30; j < 40; j++)
			Hy94_M10[i - 30][j - 30] = Hy9_M50[i][j];
	double Hy95_M10[10][10];
	for (int i = 40; i < 50; i++)
		for (int j = 40; j < 50; j++)
			Hy95_M10[i - 40][j - 40] = Hy9_M50[i][j];

	double p1[10];
	for (int i = 0; i < 10; i++)
	{
		p1[i] = X[i] - Hy9_O[i];
	}
	mul_matrix10(p1, Hy91_M10);
	double yu1 = p1[0] * p1[0];
	double yu2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yu2 += 1E6*p1[i] * p1[i];
	fh91 = yu1 + yu2;


	double p2[10];
	for (int i = 10; i < 20; i++)
	{
		p2[i - 10] = X[i] - Hy9_O[i];

	}
	mul_matrix10(p2, Hy92_M10);
	for (int i = 0; i < 10; i++)
		fh92 = fh92 + (p2[i] * p2[i] - 10 * cos(2 * pai*p2[i]) + 10);

	double p3[10];
	for (int i = 20; i < 30; i++)
	{
		p3[i - 20] = X[i] - Hy9_O[i];

	}
	mul_matrix10(p3, Hy93_M10);
	//Griewank's plus Rosenbrock
	double temp[10];
	for (int i = 0; i < Dim31 - 2; i++)
		temp[i] = 100 * (p3[i] * p3[i] - p3[i + 1])*(p3[i] * p3[i] - p3[i + 1]) + (p3[i] - 1)*(p3[i] - 1) + 100 * (p3[i + 1] * p3[i + 1] - p3[i + 2])*(p3[i + 1] * p3[i + 1] - p3[i + 2]) + (p3[i + 1] - 1)*(p3[i + 1] - 1);
	temp[Dim31 - 2] = 100 * (p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1])*(p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1]) + (p3[Dim31 - 2] - 1)*(p3[Dim31 - 2] - 1);
	temp[Dim31 - 1] = 100 * (p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0])*(p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0]) + (p3[Dim31 - 1] - 1)*(p3[Dim31 - 1] - 1);

	double u0 = 2.5;//u0
	double d = 1;//d
	double s = 1 - 1 / (2 * pow(Dim31 + 20, 0.5) - 8.2);//s
	double u1 = -pow((u0*u0 - d) / s, 0.5);//u1
	double sign[10];//signº¯Êý
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[10][10];//¶Ô½Ç¾ØÕó
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim31 - 1)));
	}
	double zn[10];//z
	double yy[10];//y  y=10(x-o)/100
	double xz[10];//x»¡
	for (int i = 0; i < Dim31; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 20]);

	for (int i = 0; i < Dim31; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim31; i++)
		zn[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim31; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim31;
	if (F1 < F2)
		fh93 = F1;
	else fh93 = F2;

	for (int i = 0; i < Dim31; i++)
		ff = ff + cos(2 * pai*zn[i]);
	fh93 = fh93 + 10 * (Dim31 - ff);


	double p4[10];
	for (int i = 30; i < 40; i++)
	{
		p4[i - 30] = X[i] - Hy9_O[i];
		//y3 += p4[i - 30] * p4[i - 30];
	}
	mul_matrix10(p4, Hy94_M10);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yyr2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yyr2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim41; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(p4[i] + 0.5));
		ff1 += yy1;
	}
	fh94 = ff1 - Dim41 * yyr2;



	double p5[10];
	for (int i = 40; i < 50; i++)
	{
		p5[i - 40] = X[i] - Hy9_O[i];
		//y3 += p5[i - 40] * p5[i - 40];
	}
	mul_matrix10(p5, Hy95_M10);
	double g = 0.0;
	for (int i = 0; i < Dim51 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(p5[i] * p5[i] + p5[i + 1] * p5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[i] * p5[i] + p5[i + 1] * p5[i + 1])), 2);
		fh95 += g;
	}
	fh95 += 0.5 + (pow(sin(pow(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0])), 2);

	fh9 = fh91 + fh92 + fh93 + fh94 + fh95;


	//	cout << "---fh9=" << fh9 << endl;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	result = w1 * (T[0] * fh5 + bias[0]) + w2 * (T[1] * fh8 + bias[1]) + w3 * (T[2] * fh9 + bias[2]);
#endif 


#if DIM==100
	//hybrid function 5
	extern double Hy5_O[100];
	extern double Hy5_M100[100][100];
	int Dim111;
	int Dim211;
	int Dim311;
	int Dim411;
	Dim111 = DIM / 10 * 2;
	Dim211 = DIM / 10 * 2;
	Dim311 = DIM / 10 * 3;
	Dim411 = DIM / 10 * 3;

	double fh51 = 0;
	double fh52 = 0;
	double fh53 = 0;
	double fh54 = 0;
	double fh5 = 0;

	//Bent Cigar function
	double Hy51_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy51_M20[i][j] = Hy5_M100[i][j];
	}
	double z1[20];
	for (int i = 0; i < Dim111; i++)
	{
		z1[i] = X[i] - Hy5_O[i];

	}
	mul_matrix20(z1, Hy51_M20);
	double yy1 = z1[0] * z1[0];
	double yy2 = 0.0;
	for (int i = 1; i < Dim111; i++)
		yy2 += 1E6*z1[i] * z1[i];
	fh51 = yy1 + yy2;


	//HGBat function
	double Hy52_M20[20][20];
	for (int i = 20; i < 40; i++)
	{
		for (int j = 20; j < 40; j++)
			Hy52_M20[i - 20][j - 20] = Hy5_M100[i][j];
	}

	double z2[20];
	for (int i = 0; i < 20; i++)
	{
		z2[i] = X[i + 20] - Hy5_O[i + 20];
		//	y1 += z2[i] * z2[i];
	}
	mul_matrix20(z2, Hy52_M20);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim211; i++)
	{
		fd1 += z2[i] * z2[i];
		fd2 += z2[i];
	}
	fh52 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim211 + 0.5;


	//Rastrigin function 
	double Hy53_M30[30][30];
	for (int i = 40; i < 70; i++)
	{
		for (int j = 40; j < 70; j++)
			Hy53_M30[i - 40][j - 40] = Hy5_M100[i][j];
	}

	double z3[30];
	for (int i = 0; i < Dim311; i++)
	{
		z3[i] = X[i + 40] - Hy5_O[i + 40];
		//	y1 += z3[i] * z3[i];
	}
	mul_matrix30(z3, Hy53_M30);
	for (int i = 0; i < Dim311; i++)
		fh53 = fh53 + (z3[i] * z3[i] - 10 * cos(2 * pai*z3[i]) + 10);

	//Rosenbrock function 
	double Hy54_M30[30][30];
	for (int i = 70; i < 100; i++)
	{
		for (int j = 70; j < 100; j++)
			Hy54_M30[i - 70][j - 70] = Hy5_M100[i][j];
	}

	double z4[30];
	for (int i = 0; i < Dim411; i++)
	{
		z4[i] = X[i + 70] - Hy5_O[i + 70];
	}
	mul_matrix30(z4, Hy54_M30);

	for (int i = 0; i < Dim411 - 1; i++)
		fh54 += 100 * (z4[i] * z4[i] - z4[i + 1]) * (z4[i] * z4[i] - z4[i + 1]) + (z4[i] - 1) * (z4[i] - 1);
	fh5 = fh51 + fh52 + fh53 + fh54;
	//cout << "fh5=" << fh5 << endl;

	//f2=hybrid function 8
	extern double Hy8_O[100];
	extern double Hy8_M100[100][100];
	double fh81 = 0;
	double fh82 = 0;
	double fh83 = 0;
	double fh84 = 0;
	double fh85 = 0;
	double fh8 = 0;

	double Hy81_M20[20][20];
	for (int i = 0; i < 20; i++)
		for (int j = 0; j < 20; j++)
			Hy81_M20[i][j] = Hy8_M100[i][j];
	double Hy82_M20[20][20];
	for (int i = 20; i < 40; i++)
		for (int j = 20; j < 40; j++)
			Hy82_M20[i - 20][j - 20] = Hy8_M100[i][j];
	double Hy83_M20[20][20];
	for (int i = 40; i < 60; i++)
		for (int j = 40; j < 60; j++)
			Hy83_M20[i - 40][j - 40] = Hy8_M100[i][j];
	double Hy84_M20[20][20];
	for (int i = 60; i < 80; i++)
		for (int j = 60; j < 80; j++)
			Hy84_M20[i - 60][j - 60] = Hy8_M100[i][j];
	double Hy85_M20[20][20];
	for (int i = 80; i < 100; i++)
		for (int j = 80; j < 100; j++)
			Hy85_M20[i - 80][j - 80] = Hy8_M100[i][j];

	double s1[20];
	for (int i = 0; i < 20; i++)
	{
		s1[i] = X[i] - Hy8_O[i];
		//	y2 += s1[i] * s1[i];
	}
	mul_matrix20(s1, Hy81_M20);

	for (int i = 0; i < 20; i++)
		fh81 += pow(pow(10, 6), i / (20 - 1))*s1[i] * s1[i];


	double s2[20];
	for (int i = 20; i < 40; i++)
	{
		s2[i - 20] = X[i] - Hy8_O[i];
	}
	mul_matrix20(s2, Hy82_M20);
	double ffx1 = 0;
	double ffx2 = 0;
	for (int i = 0; i < 20; i++)
	{
		ffx1 += s2[i] * s2[i];
		ffx2 += cos(2 * pai * s2[i]);
	}
	fh82 = -20 * exp(-0.2*sqrt((1.0 / 20) * ffx1)) - exp((1.0 / 20) * ffx2) + exp(1) + 20;


	double s3[20];
	for (int i = 40; i < 60; i++)
	{
		s3[i - 40] = X[i] - Hy8_O[i];
		//y2 += s3[i - 40] * s3[i - 40];
	}
	mul_matrix20(s3, Hy83_M20);
	for (int i = 0; i < 20; i++)
		fh83 = fh83 + (s3[i] * s3[i] - 10 * cos(2 * pai*s3[i]) + 10);

	double s4[20];
	for (int i = 60; i < 80; i++)
	{
		s4[i - 60] = X[i] - Hy8_O[i];
	}
	mul_matrix20(s4, Hy84_M20);
	double ffd1 = 0.0;
	double ffd2 = 0.0;
	for (int i = 0; i < 20; i++)
	{
		ffd1 += s4[i] * s4[i];
		ffd2 += s4[i];
	}
	fh84 = pow(fabs(ffd1*ffd1 - ffd2 * ffd2), 0.5) + (0.5*ffd1 + ffd2) / 2 + 0.5;

	//Discus
	double s5[20];
	for (int i = 80; i < 100; i++)
	{
		s5[i - 80] = X[i] - Hy8_O[i];
		//y2 += s5[i - 80] * s5[i - 80];
	}
	mul_matrix20(s5, Hy85_M20);
	double fz1 = 0.0;
	for (int i = 1; i < 20; i++)
		fz1 += s5[i] * s5[i];
	fh85 = 1E6 * s5[0] * s5[0] + fz1;

	fh8 = fh81 + fh82 + fh83 + fh84 + fh85;
	//	cout << "---fh8=" << fh8 << endl;

	extern double Hy9_O[100];
	extern double Hy9_M100[100][100];
	//f3=hybrid function 9
	int Dim11 = DIM / 10 * 2;
	int Dim21 = DIM / 10 * 2;
	int Dim31 = DIM / 10 * 2;
	int Dim41 = DIM / 10 * 2;
	int Dim51 = DIM / 10 * 2;
	double fh91 = 0;
	double fh92 = 0;
	double fh93 = 0;
	double fh94 = 0;
	double fh95 = 0;
	double fh9 = 0;

	double Hy91_M20[20][20];
	for (int i = 0; i < 20; i++)
		for (int j = 0; j < 20; j++)
			Hy91_M20[i][j] = Hy9_M100[i][j];
	double Hy92_M20[20][20];
	for (int i = 20; i < 40; i++)
		for (int j = 20; j < 40; j++)
			Hy92_M20[i - 20][j - 20] = Hy9_M100[i][j];
	double Hy93_M20[20][20];
	for (int i = 40; i < 60; i++)
		for (int j = 40; j < 60; j++)
			Hy93_M20[i - 40][j - 40] = Hy9_M100[i][j];
	double Hy94_M20[20][20];
	for (int i = 60; i < 80; i++)
		for (int j = 60; j < 80; j++)
			Hy94_M20[i - 60][j - 60] = Hy9_M100[i][j];
	double Hy95_M20[20][20];
	for (int i = 80; i < 100; i++)
		for (int j = 80; j < 100; j++)
			Hy95_M20[i - 80][j - 80] = Hy9_M100[i][j];

	double p1[20];
	for (int i = 0; i < 20; i++)
	{
		p1[i] = X[i] - Hy9_O[i];
	}
	mul_matrix20(p1, Hy91_M20);
	double yu1 = p1[0] * p1[0];
	double yu2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yu2 += 1E6*p1[i] * p1[i];
	fh91 = yu1 + yu2;


	double p2[20];
	for (int i = 20; i < 40; i++)
	{
		p2[i - 20] = X[i] - Hy9_O[i];

	}
	mul_matrix20(p2, Hy92_M20);
	for (int i = 0; i < 20; i++)
		fh92 = fh92 + (p2[i] * p2[i] - 10 * cos(2 * pai*p2[i]) + 10);

	double p3[20];
	for (int i = 40; i < 60; i++)
	{
		p3[i - 40] = X[i] - Hy9_O[i];

	}
	mul_matrix20(p3, Hy93_M20);
	//Griewank's plus Rosenbrock
	double temp[20];
	for (int i = 0; i < Dim31 - 2; i++)
		temp[i] = 100 * (p3[i] * p3[i] - p3[i + 1])*(p3[i] * p3[i] - p3[i + 1]) + (p3[i] - 1)*(p3[i] - 1) + 100 * (p3[i + 1] * p3[i + 1] - p3[i + 2])*(p3[i + 1] * p3[i + 1] - p3[i + 2]) + (p3[i + 1] - 1)*(p3[i + 1] - 1);
	temp[Dim31 - 2] = 100 * (p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1])*(p3[Dim31 - 2] * p3[Dim31 - 2] - p3[Dim31 - 1]) + (p3[Dim31 - 2] - 1)*(p3[Dim31 - 2] - 1);
	temp[Dim31 - 1] = 100 * (p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0])*(p3[Dim31 - 1] * p3[Dim31 - 1] - p3[0]) + (p3[Dim31 - 1] - 1)*(p3[Dim31 - 1] - 1);

	double u0 = 2.5;//u0
	double d = 1;//d
	double s = 1 - 1 / (2 * pow(Dim31 + 20, 0.5) - 8.2);//s
	double u1 = -pow((u0*u0 - d) / s, 0.5);//u1
	double sign[20];//signº¯Êý
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[20][20];//¶Ô½Ç¾ØÕó
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim31 - 1)));
	}
	double zn[20];//z
	double yy[20];//y  y=10(x-o)/100
	double xz[20];//x»¡
	for (int i = 0; i < Dim31; i++)
		yy[i] = 0.1*(temp[i] - Hy9_O[i + 40]);

	for (int i = 0; i < Dim31; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim31; i++)
		zn[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim31; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim31;
	if (F1 < F2)
		fh93 = F1;
	else fh93 = F2;

	for (int i = 0; i < Dim31; i++)
		ff = ff + cos(2 * pai*zn[i]);
	fh93 = fh93 + 10 * (Dim31 - ff);


	double p4[20];
	for (int i = 60; i < 80; i++)
	{
		p4[i - 60] = X[i] - Hy9_O[i];
		//y3 += p4[i - 60] * p4[i - 60];
	}
	mul_matrix20(p4, Hy94_M20);
	double a = 0.5;
	int b = 3;
	int k_max = 20;
	double ff1 = 0.0;
	double yyr2 = 0.0;
	for (int j = 0; j <= k_max; j++)
	{
		yyr2 += pow(a, j)*cos(pai*pow(b, j));
	}
	for (int i = 0; i < Dim41; i++)
	{
		double yy1 = 0.0;
		for (int j = 0; j <= k_max; j++)
			yy1 += pow(a, j)*cos(2 * pai*pow(b, j)*(p4[i] + 0.5));
		ff1 += yy1;
	}
	fh94 = ff1 - Dim41 * yyr2;



	double p5[20];
	for (int i = 80; i < 100; i++)
	{
		p5[i - 80] = X[i] - Hy9_O[i];
		//y3 += p5[i - 80] * p5[i - 80];
	}
	mul_matrix20(p5, Hy95_M20);
	double g = 0.0;
	for (int i = 0; i < Dim51 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(p5[i] * p5[i] + p5[i + 1] * p5[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[i] * p5[i] + p5[i + 1] * p5[i + 1])), 2);
		fh95 += g;
	}
	fh95 += 0.5 + (pow(sin(pow(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(p5[Dim51 - 1] * p5[Dim51 - 1] + p5[0] * p5[0])), 2);

	fh9 = fh91 + fh92 + fh93 + fh94 + fh95;


	//	cout << "---fh9=" << fh9 << endl;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;

	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));
	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	result = w1 * (T[0] * fh5 + bias[0]) + w2 * (T[1] * fh8 + bias[1]) + w3 * (T[2] * fh9 + bias[2]);
#endif

	return result + 2900;
#endif
	/*
	#if Composition_function_10
		double f = 0;
		double y1 = 0;
		double y2 = 0;
		double y3 = 0;
		double R[3] = { 10,30,50 };
		double T[3] = { 1,1,1 };
		double bias[3] = { 0,100,200 };
		extern double C101_O[100];
		extern double C102_O[100];
		extern double C103_O[100];
		for (int i = 0; i < DIM; i++)
		{
			y1 += (X[i] - C101_O[i])*(X[i] - C101_O[i]);
			y2 += (X[i] - C102_O[i])*(X[i] - C102_O[i]);
			y3 += (X[i] - C103_O[i])*(X[i] - C103_O[i]);

		}


	#if DIM==10
		//hybrid function 5
			//BentCigar+HGBat+Rastrigin+Rosenbrock
		double p[4] = { 0.2,0.2,0.3,0.3 };
		int Dim11;
		int Dim12;
		int Dim13;
		int Dim14;
		Dim11 = DIM / 10 * 2;//6Î¬
		Dim12 = DIM / 10 * 2;//6Î¬
		Dim13 = DIM / 10 * 3;//12Î¬
		Dim14 = DIM / 10 * 3;//12Î¬
		double f11 = 0;
		double f12 = 0;
		double f13 = 0;
		double f14 = 0;
		double f1 = 0;

		extern double Hy5_O[100];
		extern double Hy5_M10[10][10];

		//Bent Cigar function
		double Hy51_M2[2][2];
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				Hy51_M2[i][j] = Hy5_M10[i][j];
		}
		double z11[2];
		for (int i = 0; i < Dim11; i++)
			z11[i] = X[i] - Hy5_O[i];
		mul_matrix02(z11, Hy51_M2);
		double yyw1 = z11[0] * z11[0];
		double yyw2 = 0.0;
		for (int i = 1; i < Dim11; i++)
			yyw2 += 1E6*z11[i] * z11[i];
		f11 = yyw1 + yyw2;


		//HGBat function
		double Hy52_M2[2][2];
		for (int i = 2; i < 4; i++)
		{
			for (int j = 2; j < 4; j++)
				Hy52_M2[i - 2][j - 4] = Hy5_M10[i][j];
		}

		double z12[2];
		for (int i = 0; i < 2; i++)
			z12[i] = X[i + 2] - Hy5_O[i + 2];
		mul_matrix02(z12, Hy52_M2);
		double fd1 = 0.0;
		double fd2 = 0.0;
		for (int i = 0; i < Dim12; i++)
		{
			fd1 += z12[i] * z12[i];
			fd2 += z12[i];
		}
		f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


		//Rastrigin function
		double Hy53_M3[3][3];
		for (int i = 4; i < 7; i++)
		{
			for (int j = 4; j < 7; j++)
				Hy53_M3[i - 4][j - 4] = Hy5_M10[i][j];
		}

		double z13[3];
		for (int i = 0; i < Dim13; i++)
			z13[i] = X[i + 4] - Hy5_O[i + 4];
		mul_matrix03(z13, Hy53_M3);
		for (int i = 0; i < Dim13; i++)
			f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

		//Rosenbrock function
		double Hy54_M3[3][3];
		for (int i = 7; i < 10; i++)
		{
			for (int j = 7; j < 10; j++)
				Hy54_M3[i - 7][j - 7] = Hy5_M10[i][j];
		}

		double z14[3];
		for (int i = 0; i < Dim14; i++)
			z14[i] = X[i + 7] - Hy5_O[i + 7];
		mul_matrix03(z14, Hy54_M3);

		for (int i = 0; i < Dim14 - 1; i++)
			f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


		f1 = f11 + f12 + f13 + f14;


		//hybrid function 6
		//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
		double p1[4] = { 0.2,0.2,0.3,0.3 };
		int Dim21;
		int Dim22;
		int Dim23;
		int Dim24;

		Dim21 = DIM / 10 * 2;//6Î¬
		Dim22 = DIM / 10 * 2;//6Î¬
		Dim23 = DIM / 10 * 3;//9Î¬
		Dim24 = DIM / 10 * 3;//9Î¬

		double f21 = 0;
		double f22 = 0;
		double f23 = 0;
		double f24 = 0;
		double f2 = 0;

		extern double Hy6_O[100];
		extern double Hy6_M10[10][10];

		//Expanded Schaffer F6
		double Hy61_M2[2][2];
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				Hy61_M2[i][j] = Hy2_M10[i][j];
		}
		double z21[2];
		for (int i = 0; i < Dim21; i++)
			z21[i] = X[i] - Hy6_O[i];
		mul_matrix02(z21, Hy61_M2);
		double g = 0.0;
		for (int i = 0; i < Dim21 - 1; i++)
		{
			g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
			f21 += g;
		}
		f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



		//HGBat
		double Hy62_M2[2][2];
		for (int i = 2; i < 4; i++)
		{
			for (int j = 2; j < 4; j++)
				Hy62_M6[i - 2][j - 2] = Hy6_M10[i][j];
		}

		double z22[2];
		for (int i = 0; i < 2; i++)
			z22[i] = X[i + 2] - Hy6_O[i + 2];
		mul_matrix02(z22, Hy62_M2);
		double ff1 = 0.0;
		double ff2 = 0.0;
		for (int i = 0; i < Dim22; i++)
		{
			ff1 += z22[i] * z22[i];
			ff2 += z22[i];
		}
		f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



		//Rosenbrock
		double Hy63_M3[3][3];
		for (int i = 4; i < 7; i++)
		{
			for (int j = 4; j < 7; j++)
				Hy63_M3[i - 4][j - 4] = Hy6_M10[i][j];
		}

		double z23[3];
		for (int i = 0; i < Dim23; i++)
			z23[i] = X[i + 4] - Hy6_O[i + 4];
		mul_matrix03(z23, Hy63_M3);
		for (int i = 0; i < Dim23 - 1; i++)
			f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


		//Modified Schwefel F10
		double Hy64_M3[3][3];
		for (int i = 7; i < 10; i++)
		{
			for (int j = 7; j < 10; j++)
				Hy64_M9[i - 7][j - 7] = Hy6_M10[i][j];
		}

		double z24[3];
		for (int i = 0; i < Dim24; i++)
			z24[i] = X[i + 7] - Hy6_O[i + 7];
		mul_matrix03(z24, Hy64_M3);

		for (int i = 0; i < Dim24; i++)
		{
			z24[i] = z24[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z24[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z24[i]) <= 500)
				f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
			if (z24[i] > 500)
				f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
			if (z24[i] < -500)
				f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
		}
		f24 = 418.9829*Dim24 - f24;

		f2 = f21 + f22 + f23 + f24;


		// hybrid_function_7
		//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
		double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
		int Dim1 = DIM / 10 * 1;
		int Dim2 = DIM / 10 * 2;
		int Dim3 = DIM / 10 * 2;
		int Dim4 = DIM / 10 * 2;
		int Dim5 = DIM / 10 * 3;
		double f31 = 0;
		double f32 = 0;
		double f33 = 0;
		double f34 = 0;
		double f35 = 0;
		double f3 = 0;

		extern double Hy7_O[100];
		extern double Hy7_M10[10][10];
		double Hy71_M1[1][1];
		for (int i = 0; i < 1; i++)
			for (int j = 0; j < 1; j++)
				Hy71_M1[i][j] = Hy7_M10[i][j];
		double Hy72_M2[2][2];
		for (int i = 1; i < 3; i++)
			for (int j = 1; j < 3; j++)
				Hy72_M2[i - 1][j - 1] = Hy7_M10[i][j];
		double Hy73_M2[2][2];
		for (int i = 3; i < 5; i++)
			for (int j = 3; j < 5; j++)
				Hy73_M2[i - 3][j - 3] = Hy7_M10[i][j];
		double Hy74_M2[2][2];
		for (int i = 5; i < 7; i++)
			for (int j = 5; j < 7; j++)
				Hy74_M2[i - 5][j - 5] = Hy7_M10[i][j];
		double Hy75_M3[3][3];
		for (int i = 7; i < 10; i++)
			for (int j = 7; j < 10; j++)
				Hy75_M3[i - 7][j - 7] = Hy7_M10[i][j];

		double z1[1];
		for (int i = 0; i < 1; i++)
			z1[i] = X[i] - Hy7_O[i];
		mul_matrix01(z1, Hy71_M1);

		double yyy = 1.0;
		double yy1 = 0.0;
		double yy2 = 0.0;
		double yy3 = 0.0;
		yy1 = 10 / (Dim1*Dim1);
		yy3 = 10 / pow(Dim1, 1.2);
		for (int i = 0; i < Dim1; i++)
		{
			for (int j = 1; j <= 32; j++)
				yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
			yyy *= pow((1 + (i + 1)*yy2), yy3);
		}
		f31 = yy1 * yyy - yy1;


		double z2[2];
		for (int i = 1; i < 3; i++)
			z2[i - 1] = X[i] - Hy7_O[i];
		mul_matrix01(z2, Hy72_M1);

		double ffs1 = 0, ffs2 = 0;
		for (int i = 0; i < Dim2; i++) {
			ffs1 += z2[i] * z2[i];
			ffs2 += cos(2 * pai * z2[i]);
		}
		f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

		double z3[2];
		for (int i = 3; i < 5; i++)
			z3[i - 3] = X[i] - Hy7_O[i];
		mul_matrix02(z3, Hy73_M2);
		//Expanded Griewank's plus Rosenbrock
		double temp[2];
		for (int i = 0; i < Dim3 - 1; i++)
			temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
		temp[Dim3 - 1] = 100 * (z3[0] * z3[0] - z3[1])*(z3[0] * z3[0] - z3[1]) + (z3[0] - 1)*(z3[0] - 1);

		double u0 = 2.5;
		double d = 1;
		double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
		double u1 = -pow((u0*u0 - d) / s, 0.5);
		double sign[2];
		double F1 = 0;
		double F2 = 0;
		double ff = 0;
		double v[2][2];
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				v[i][j] = 0;
			v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
		}
		double z[2];
		double yy[2];
		double xz[2];//x»¡
		for (int i = 0; i < Dim3; i++)
			yy[i] = 0.1*(temp[i] - Hy7_O[i + 3]);

		for (int i = 0; i < Dim3; i++)
		{
			if (temp[i] < 0)
				sign[i] = -1;
			else if (temp[i] == 0)
				sign[i] = 0;
			else sign[i] = 1;
			xz[i] = 2 * sign[i] * yy[i] + u0;

		}

		for (int i = 0; i < Dim3; i++)
			z[i] = v[i][i] * (xz[i] - u0);
		for (int i = 0; i < Dim3; i++)
		{
			F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
			F2 = s * (xz[i] - u1)*(xz[i] - u1);
		}
		F2 = F2 + d * Dim3;
		if (F1 < F2)
			f33 = F1;
		else f33 = F2;

		for (int i = 0; i < Dim3; i++)
			ff = ff + cos(2 * pai*z[i]);
		f33 = f33 + 10 * (Dim3 - ff);


		double z4[2];
		for (int i = 5; i < 7; i++)
			z4[i - 5] = X[i] - Hy7_O[i];
		mul_matrix02(z4, Hy74_M2);

		for (int i = 0; i < Dim4; i++)
		{
			z4[i] = z4[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z4[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z4[i]) <= 500)
				f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
			if (z4[i] > 500)
				f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
			if (z4[i] < -500)
				f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
		}
		f34 = 418.9829*Dim4 - f34;


		double z5[3];
		for (int i = 7; i < 10; i++)
			z5[i - 7] = X[i] - Hy7_O[i];
		mul_matrix03(z5, Hy75_M3);
		//Rastrigin
		for (int i = 0; i < Dim5; i++)
			f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

		f3 = f31 + f32 + f33 + f34 + f35;
		double w1 = 0;
		double w2 = 0;
		double w3 = 0;
		double w = 0;
		w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
		w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
		w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

		w = w1 + w2 + w3;
		w1 = w1 / w;
		w2 = w2 / w;
		w3 = w3 / w;

		f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
	#endif


	#if DIM==30
		//hybrid function 5
			//BentCigar+HGBat+Rastrigin+Rosenbrock
		double p[4] = { 0.2,0.2,0.3,0.3 };
		int Dim11;
		int Dim12;
		int Dim13;
		int Dim14;
		Dim11 = DIM / 10 * 2;//6Î¬
		Dim12 = DIM / 10 * 2;//6Î¬
		Dim13 = DIM / 10 * 3;//12Î¬
		Dim14 = DIM / 10 * 3;//12Î¬
		double f11 = 0;
		double f12 = 0;
		double f13 = 0;
		double f14 = 0;
		double f1 = 0;

		extern double Hy5_O[100];
		extern double Hy5_M30[30][30];

		//Bent Cigar function
		double Hy51_M6[6][6];
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
				Hy51_M6[i][j] = Hy5_M30[i][j];
		}
		double z11[6];
		for (int i = 0; i < Dim11; i++)
			z11[i] = X[i] - Hy5_O[i];
		mul_matrix06(z11, Hy51_M6);
		double yyw1 = z11[0] * z11[0];
		double yyw2 = 0.0;
		for (int i = 1; i < Dim11; i++)
			yyw2 += 1E6*z11[i] * z11[i];
		f11 = yyw1 + yyw2;


		//HGBat function
		double Hy52_M6[6][6];
		for (int i = 6; i < 12; i++)
		{
			for (int j = 6; j < 12; j++)
				Hy52_M6[i - 6][j - 6] = Hy5_M30[i][j];
		}

		double z12[6];
		for (int i = 0; i < 6; i++)
			z12[i] = X[i + 6] - Hy5_O[i + 6];
		mul_matrix06(z12, Hy52_M6);
		double fd1 = 0.0;
		double fd2 = 0.0;
		for (int i = 0; i < Dim12; i++)
		{
			fd1 += z12[i] * z12[i];
			fd2 += z12[i];
		}
		f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


		//Rastrigin function
		double Hy53_M9[9][9];
		for (int i = 12; i < 21; i++)
		{
			for (int j = 12; j < 21; j++)
				Hy53_M9[i - 12][j - 12] = Hy5_M30[i][j];
		}

		double z13[9];
		for (int i = 0; i < Dim13; i++)
			z13[i] = X[i + 12] - Hy5_O[i + 12];
		mul_matrix09(z13, Hy53_M9);
		for (int i = 0; i < Dim13; i++)
			f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

		//Rosenbrock function
		double Hy54_M9[9][9];
		for (int i = 21; i < 30; i++)
		{
			for (int j = 21; j < 30; j++)
				Hy54_M9[i - 21][j - 21] = Hy5_M30[i][j];
		}

		double z14[9];
		for (int i = 0; i < Dim14; i++)
			z14[i] = X[i + 21] - Hy5_O[i + 21];
		mul_matrix09(z14, Hy54_M9);

		for (int i = 0; i < Dim14 - 1; i++)
			f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


		f1 = f11 + f12 + f13 + f14;


		//hybrid function 6
		//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
		double p1[4] = { 0.2,0.2,0.3,0.3 };
		int Dim21;
		int Dim22;
		int Dim23;
		int Dim24;

		Dim21 = DIM / 10 * 2;//6Î¬
		Dim22 = DIM / 10 * 2;//6Î¬
		Dim23 = DIM / 10 * 3;//9Î¬
		Dim24 = DIM / 10 * 3;//9Î¬

		double f21 = 0;
		double f22 = 0;
		double f23 = 0;
		double f24 = 0;
		double f2 = 0;

		extern double Hy6_O[100];
		extern double Hy6_M30[30][30];

		//Expanded Schaffer F6
		double Hy61_M6[6][6];
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
				Hy61_M6[i][j] = Hy6_M30[i][j];
		}
		double z21[6];
		for (int i = 0; i < Dim21; i++)
			z21[i] = X[i] - Hy6_O[i];
		mul_matrix06(z21, Hy61_M6);
		double g = 0.0;
		for (int i = 0; i < Dim21 - 1; i++)
		{
			g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
			f21 += g;
		}
		f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



		//HGBat
		double Hy62_M6[6][6];
		for (int i = 6; i < 12; i++)
		{
			for (int j = 6; j < 12; j++)
				Hy62_M6[i - 6][j - 6] = Hy6_M30[i][j];
		}

		double z22[6];
		for (int i = 0; i < 6; i++)
			z22[i] = X[i + 6] - Hy6_O[i + 6];
		mul_matrix06(z22, Hy62_M6);
		double ff1 = 0.0;
		double ff2 = 0.0;
		for (int i = 0; i < Dim22; i++)
		{
			ff1 += z22[i] * z22[i];
			ff2 += z22[i];
		}
		f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



		//Rosenbrock
		double Hy63_M9[9][9];
		for (int i = 12; i < 21; i++)
		{
			for (int j = 12; j < 21; j++)
				Hy63_M9[i - 12][j - 12] = Hy6_M30[i][j];
		}

		double z23[9];
		for (int i = 0; i < Dim23; i++)
			z23[i] = X[i + 12] - Hy6_O[i + 12];
		mul_matrix09(z23, Hy63_M9);
		for (int i = 0; i < Dim23 - 1; i++)
			f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


		//Modified Schwefel F10
		double Hy64_M9[9][9];
		for (int i = 21; i < 30; i++)
		{
			for (int j = 21; j < 30; j++)
				Hy64_M9[i - 21][j - 21] = Hy6_M30[i][j];
		}

		double z24[9];
		for (int i = 0; i < Dim24; i++)
			z24[i] = X[i + 21] - Hy6_O[i + 21];
		mul_matrix09(z24, Hy64_M9);

		for (int i = 0; i < Dim24; i++)
		{
			z24[i] = z24[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z24[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z24[i]) <= 500)
				f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
			if (z24[i] > 500)
				f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
			if (z24[i] < -500)
				f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
		}
		f24 = 418.9829*Dim24 - f24;

		f2 = f21 + f22 + f23 + f24;


		// hybrid_function_7
		//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
		//double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
		int Dim1 = DIM / 10 * 1;
		int Dim2 = DIM / 10 * 2;
		int Dim3 = DIM / 10 * 2;
		int Dim4 = DIM / 10 * 2;
		int Dim5 = DIM / 10 * 3;
		double f31 = 0;
		double f32 = 0;
		double f33 = 0;
		double f34 = 0;
		double f35 = 0;
		double f3 = 0;

		extern double Hy7_O[100];
		extern double Hy7_M30[30][30];
		double Hy71_M3[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				Hy71_M3[i][j] = Hy7_M30[i][j];
		double Hy72_M6[6][6];
		for (int i = 3; i < 9; i++)
			for (int j = 3; j < 9; j++)
				Hy72_M6[i - 3][j - 3] = Hy7_M30[i][j];
		double Hy73_M6[6][6];
		for (int i = 9; i < 15; i++)
			for (int j = 9; j < 15; j++)
				Hy73_M6[i - 9][j - 9] = Hy7_M30[i][j];
		double Hy74_M6[6][6];
		for (int i = 15; i < 21; i++)
			for (int j = 15; j < 21; j++)
				Hy74_M6[i - 15][j - 15] = Hy7_M30[i][j];
		double Hy75_M9[9][9];
		for (int i = 21; i < 30; i++)
			for (int j = 21; j < 30; j++)
				Hy75_M9[i - 21][j - 21] = Hy7_M30[i][j];

		double z1[3];
		for (int i = 0; i < 3; i++)
			z1[i] = X[i] - Hy7_O[i];
		mul_matrix03(z1, Hy71_M3);

		double yyy = 1.0;
		double yy1 = 0.0;
		double yy2 = 0.0;
		double yy3 = 0.0;
		yy1 = 10 / (Dim1*Dim1);
		yy3 = 10 / pow(Dim1, 1.2);
		for (int i = 0; i < Dim1; i++)
		{
			for (int j = 1; j <= 32; j++)
				yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
			yyy *= pow((1 + (i + 1)*yy2), yy3);
		}
		f31 = yy1 * yyy - yy1;


		double z2[6];
		for (int i = 3; i < 9; i++)
			z2[i - 3] = X[i] - Hy7_O[i];
		mul_matrix06(z2, Hy72_M6);

		double ffs1 = 0, ffs2 = 0;
		for (int i = 0; i < Dim2; i++) {
			ffs1 += z2[i] * z2[i];
			ffs2 += cos(2 * pai * z2[i]);
		}
		f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

		double z3[6];
		for (int i = 9; i < 15; i++)
			z3[i - 9] = X[i] - Hy7_O[i];
		mul_matrix06(z3, Hy73_M6);
		//Expanded Griewank's plus Rosenbrock
		double temp[6];
		for (int i = 0; i < Dim3 - 1; i++)
			temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
		temp[Dim3 - 1] = 100 * (z3[0] * z3[0] - z3[1])*(z3[0] * z3[0] - z3[1]) + (z3[0] - 1)*(z3[0] - 1);

		double u0 = 2.5;
		double d = 1;
		double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
		double u1 = -pow((u0*u0 - d) / s, 0.5);
		double sign[6];
		double F1 = 0;
		double F2 = 0;
		double ff = 0;
		double v[6][6];
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
				v[i][j] = 0;
			v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
		}
		double z[6];
		double yy[6];
		double xz[6];//x»¡
		for (int i = 0; i < Dim3; i++)
			yy[i] = 0.1*(temp[i] - Hy7_O[i + 9]);

		for (int i = 0; i < Dim3; i++)
		{
			if (temp[i] < 0)
				sign[i] = -1;
			else if (temp[i] == 0)
				sign[i] = 0;
			else sign[i] = 1;
			xz[i] = 2 * sign[i] * yy[i] + u0;

		}

		for (int i = 0; i < Dim3; i++)
			z[i] = v[i][i] * (xz[i] - u0);
		for (int i = 0; i < Dim3; i++)
		{
			F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
			F2 = s * (xz[i] - u1)*(xz[i] - u1);
		}
		F2 = F2 + d * Dim3;
		if (F1 < F2)
			f33 = F1;
		else f33 = F2;

		for (int i = 0; i < Dim3; i++)
			ff = ff + cos(2 * pai*z[i]);
		f33 = f33 + 10 * (Dim3 - ff);


		double z4[6];
		for (int i = 15; i < 21; i++)
			z4[i - 15] = X[i] - Hy7_O[i];
		mul_matrix06(z4, Hy74_M6);

		for (int i = 0; i < Dim4; i++)
		{
			z4[i] = z4[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z4[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z4[i]) <= 500)
				f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
			if (z4[i] > 500)
				f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
			if (z4[i] < -500)
				f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
		}
		f34 = 418.9829*Dim4 - f34;


		double z5[9];
		for (int i = 21; i < 30; i++)
			z5[i - 21] = X[i] - Hy7_O[i];
		mul_matrix09(z5, Hy75_M9);
		//Rastrigin
		for (int i = 0; i < Dim5; i++)
			f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

		f3 = f31 + f32 + f33 + f34 + f35;
		double w1 = 0;
		double w2 = 0;
		double w3 = 0;
		double w = 0;
		w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
		w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
		w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

		w = w1 + w2 + w3;
		w1 = w1 / w;
		w2 = w2 / w;
		w3 = w3 / w;

		f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
	#endif

	#if DIM==50
		//hybrid function 5
			//BentCigar+HGBat+Rastrigin+Rosenbrock
		double p[4] = { 0.2,0.2,0.3,0.3 };
		int Dim11;
		int Dim12;
		int Dim13;
		int Dim14;
		Dim11 = DIM / 10 * 2;//6Î¬
		Dim12 = DIM / 10 * 2;//6Î¬
		Dim13 = DIM / 10 * 3;//12Î¬
		Dim14 = DIM / 10 * 3;//12Î¬
		double f11 = 0;
		double f12 = 0;
		double f13 = 0;
		double f14 = 0;
		double f1 = 0;

		extern double Hy5_O[100];
		extern double Hy5_M50[10][10];

		//Bent Cigar function
		double Hy51_M10[10][10];
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
				Hy51_M10[i][j] = Hy5_M50[i][j];
		}
		double z11[10];
		for (int i = 0; i < Dim11; i++)
			z11[i] = X[i] - Hy5_O[i];
		mul_matrix10(z11, Hy51_M10);
		double yyw1 = z11[0] * z11[0];
		double yyw2 = 0.0;
		for (int i = 1; i < Dim11; i++)
			yyw2 += 1E6*z11[i] * z11[i];
		f11 = yyw1 + yyw2;


		//HGBat function
		double Hy52_M10[10][10];
		for (int i = 10; i < 20; i++)
		{
			for (int j = 10; j < 20; j++)
				Hy52_M2[i - 10][j - 10] = Hy5_M50[i][j];
		}

		double z12[10];
		for (int i = 0; i < 10; i++)
			z12[i] = X[i + 10] - Hy5_O[i + 10];
		mul_matrix10(z12, Hy52_M2);
		double fd1 = 0.0;
		double fd2 = 0.0;
		for (int i = 0; i < Dim12; i++)
		{
			fd1 += z12[i] * z12[i];
			fd2 += z12[i];
		}
		f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


		//Rastrigin function
		double Hy53_M15[15][15];
		for (int i = 20; i < 35; i++)
		{
			for (int j = 20; j < 35; j++)
				Hy53_M15[i - 20][j - 20] = Hy5_M50[i][j];
		}

		double z13[15];
		for (int i = 0; i < Dim13; i++)
			z13[i] = X[i + 20] - Hy5_O[i + 20];
		mul_matrix03(z13, Hy53_M3);
		for (int i = 0; i < Dim13; i++)
			f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

		//Rosenbrock function
		double Hy54_M15[15][15];
		for (int i = 35; i < 50; i++)
		{
			for (int j = 35; j < 50; j++)
				Hy54_M3[i - 35][j - 35] = Hy5_M50[i][j];
		}

		double z14[15];
		for (int i = 0; i < Dim14; i++)
			z14[i] = X[i + 35] - Hy5_O[i + 35];
		mul_matrix03(z14, Hy54_M3);

		for (int i = 0; i < Dim14 - 1; i++)
			f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


		f1 = f11 + f12 + f13 + f14;


		//hybrid function 6
		//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
		double p1[4] = { 0.2,0.2,0.3,0.3 };
		int Dim21;
		int Dim22;
		int Dim23;
		int Dim24;

		Dim21 = DIM / 10 * 2;//6Î¬
		Dim22 = DIM / 10 * 2;//6Î¬
		Dim23 = DIM / 10 * 3;//9Î¬
		Dim24 = DIM / 10 * 3;//9Î¬

		double f21 = 0;
		double f22 = 0;
		double f23 = 0;
		double f24 = 0;
		double f2 = 0;

		extern double Hy6_O[100];
		extern double Hy6_M50[50][50];

		//Expanded Schaffer F6
		double Hy61_M15[15][15];
		for (int i = 0; i < 15; i++)
		{
			for (int j = 0; j < 15; j++)
				Hy61_M15[i][j] = Hy2_M50[i][j];
		}
		double z21[15];
		for (int i = 0; i < Dim21; i++)
			z21[i] = X[i] - Hy6_O[i];
		mul_matrix15(z21, Hy61_M15);
		double g = 0.0;
		for (int i = 0; i < Dim21-1; i++)
		{
			g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
			f21 += g;
		}
		f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



		//HGBat
		double Hy62_M10[10][10];
		for (int i = 10; i < 20; i++)
		{
			for (int j = 10; j < 20; j++)
				Hy62_M6[i - 10][j - 10] = Hy6_M50[i][j];
		}

		double z22[10];
		for (int i = 0; i < 10; i++)
			z22[i] = X[i + 10] - Hy6_O[i + 10];
		mul_matrix10(z22, Hy62_M10);
		double ff1 = 0.0;
		double ff2 = 0.0;
		for (int i = 0; i < Dim22; i++)
		{
			ff1 += z22[i] * z22[i];
			ff2 += z22[i];
		}
		f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



		//Rosenbrock
		double Hy63_M15[15][15];
		for (int i = 20; i < 35; i++)
		{
			for (int j = 20; j < 35; j++)
				Hy63_M3[i - 20][j - 20] = Hy6_M50[i][j];
		}

		double z23[15];
		for (int i = 0; i < Dim23; i++)
			z23[i] = X[i + 20] - Hy6_O[i + 20];
		mul_matrix15(z23, Hy63_M15);
		for (int i = 0; i < Dim23 - 1; i++)
			f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


		//Modified Schwefel F10
		double Hy64_M15[15][15];
		for (int i = 35; i < 50; i++)
		{
			for (int j = 35; j < 50; j++)
				Hy64_M9[i - 35][j - 35] = Hy6_M50[i][j];
		}

		double z24[15];
		for (int i = 0; i < Dim24; i++)
			z24[i] = X[i + 35] - Hy6_O[i + 35];
		mul_matrix15(z24, Hy64_M15);

		for (int i = 0; i < Dim24; i++)
		{
			z24[i] = z24[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z24[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z24[i]) <= 500)
				f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
			if (z24[i] > 500)
				f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
			if (z24[i] < -500)
				f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
		}
		f24 = 418.9829*Dim24 - f24;

		f2 = f21 + f22 + f23 + f24;


		// hybrid_function_7
		//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
		double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
		int Dim1 = DIM / 10 * 1;
		int Dim2 = DIM / 10 * 2;
		int Dim3 = DIM / 10 * 2;
		int Dim4 = DIM / 10 * 2;
		int Dim5 = DIM / 10 * 3;
		double f31 = 0;
		double f32 = 0;
		double f33 = 0;
		double f34 = 0;
		double f35 = 0;
		double f3 = 0;

		extern double Hy7_O[100];
		extern double Hy7_M50[50][50];
		double Hy71_M1[5][5];
		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				Hy71_M5[i][j] = Hy7_M50[i][j];
		double Hy72_M10[10][10];
		for (int i = 5; i < 15; i++)
			for (int j = 5; j < 15; j++)
				Hy72_M10[i - 1][j - 1] = Hy7_M50[i][j];
		double Hy73_M10[10][10];
		for (int i = 15; i < 25; i++)
			for (int j = 15; j < 25; j++)
				Hy73_M2[i - 15][j - 15] = Hy7_M50[i][j];
		double Hy74_M10[10][10];
		for (int i = 25; i < 35; i++)
			for (int j = 25; j < 35; j++)
				Hy74_M2[i - 25][j - 35] = Hy7_M50[i][j];
		double Hy75_M15[15][15];
		for (int i = 35; i < 50; i++)
			for (int j = 35; j < 50; j++)
				Hy75_M15[i - 35][j - 35] = Hy7_M50[i][j];

		double z1[5];
		for (int i = 0; i < 5; i++)
			z1[i] = X[i] - Hy7_O[i];
		mul_matrix05(z1, Hy71_M5);

		double yyy = 1.0;
		double yy1 = 0.0;
		double yy2 = 0.0;
		double yy3 = 0.0;
		yy1 = 10 / (Dim1*Dim1);
		yy3 = 10 / pow(Dim1, 1.2);
		for (int i = 0; i < Dim1; i++)
		{
			for (int j = 1; j <= 32; j++)
				yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
			yyy *= pow((1 + (i + 1)*yy2), yy3);
		}
		f31 = yy1 * yyy - yy1;


		double z10[10];
		for (int i = 5; i < 15; i++)
			z2[i - 5] = X[i] - Hy7_O[i];
		mul_matrix05(z2, Hy72_M5);

		double ffs1 = 0, ffs2 = 0;
		for (int i = 0; i < Dim2; i++) {
			ffs1 += z2[i] * z2[i];
			ffs2 += cos(2 * pai * z2[i]);
		}
		f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

		double z3[10];
		for (int i = 15; i < 25; i++)
			z3[i - 15] = X[i] - Hy7_O[i];
		mul_matrix10(z3, Hy73_M10);
		//Expanded Griewank's plus Rosenbrock
		double temp[10];
		for (int i = 0; i < Dim3 - 1; i++)
			temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
		temp[Dim3 - 1] = 100 * (z3[0] * z3[0] - z3[1])*(z3[0] * z3[0] - z3[1]) + (z3[0] - 1)*(z3[0] - 1);

		double u0 = 2.5;
		double d = 1;
		double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
		double u1 = -pow((u0*u0 - d) / s, 0.5);
		double sign[2];
		double F1 = 0;
		double F2 = 0;
		double ff = 0;
		double v[10][10];
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
				v[i][j] = 0;
			v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
		}
		double z[10];
		double yy[10];
		double xz[10];//x»¡
		for (int i = 0; i < Dim3; i++)
			yy[i] = 0.1*(temp[i] - Hy7_O[i + 15]);

		for (int i = 0; i < Dim3; i++)
		{
			if (temp[i] < 0)
				sign[i] = -1;
			else if (temp[i] == 0)
				sign[i] = 0;
			else sign[i] = 1;
			xz[i] = 2 * sign[i] * yy[i] + u0;

		}

		for (int i = 0; i < Dim3; i++)
			z[i] = v[i][i] * (xz[i] - u0);
		for (int i = 0; i < Dim3; i++)
		{
			F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
			F2 = s * (xz[i] - u1)*(xz[i] - u1);
		}
		F2 = F2 + d * Dim3;
		if (F1 < F2)
			f33 = F1;
		else f33 = F2;

		for (int i = 0; i < Dim3; i++)
			ff = ff + cos(2 * pai*z[i]);
		f33 = f33 + 10 * (Dim3 - ff);


		double z4[10];
		for (int i = 25; i < 35; i++)
			z4[i - 25] = X[i] - Hy7_O[i];
		mul_matrix10(z4, Hy74_M10);

		for (int i = 0; i < Dim4; i++)
		{
			z4[i] = z4[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z4[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z4[i]) <= 500)
				f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
			if (z4[i] > 500)
				f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
			if (z4[i] < -500)
				f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
		}
		f34 = 418.9829*Dim4 - f34;


		double z5[15];
		for (int i = 35; i < 50; i++)
			z5[i - 35] = X[i] - Hy7_O[i];
		mul_matrix15(z5, Hy75_M15);
		//Rastrigin
		for (int i = 0; i < Dim5; i++)
			f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

		f3 = f31 + f32 + f33 + f34 + f35;
		double w1 = 0;
		double w2 = 0;
		double w3 = 0;
		double w = 0;
		w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
		w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
		w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

		w = w1 + w2 + w3;
		w1 = w1 / w;
		w2 = w2 / w;
		w3 = w3 / w;

		f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
	#endif

	#if DIM==100
		//hybrid function 5
			//BentCigar+HGBat+Rastrigin+Rosenbrock
		double p[4] = { 0.2,0.2,0.3,0.3 };
		int Dim11;
		int Dim12;
		int Dim13;
		int Dim14;
		Dim11 = DIM / 10 * 2;//6Î¬
		Dim12 = DIM / 10 * 2;//6Î¬
		Dim13 = DIM / 10 * 3;//12Î¬
		Dim14 = DIM / 10 * 3;//12Î¬
		double f11 = 0;
		double f12 = 0;
		double f13 = 0;
		double f14 = 0;
		double f1 = 0;

		extern double Hy5_O[100];
		extern double Hy5_M100[100][100];

		//Bent Cigar function
		double Hy51_M20[20][20];
		for (int i = 0; i < 20; i++)
		{
			for (int j = 0; j < 20; j++)
				Hy51_M20[i][j] = Hy5_M100[i][j];
		}
		double z11[20];
		for (int i = 0; i < Dim11; i++)
			z11[i] = X[i] - Hy5_O[i];
		mul_matrix20(z11, Hy51_M20);
		double yyw1 = z11[0] * z11[0];
		double yyw2 = 0.0;
		for (int i = 1; i < Dim11; i++)
			yyw2 += 1E6*z11[i] * z11[i];
		f11 = yyw1 + yyw2;


		//HGBat function
		double Hy52_M20[20][20];
		for (int i = 20; i < 40; i++)
		{
			for (int j = 20; j < 40; j++)
				Hy52_M2[i - 20][j - 40] = Hy5_M100[i][j];
		}

		double z12[20];
		for (int i = 0; i < 20; i++)
			z12[i] = X[i + 20] - Hy5_O[i + 20];
		mul_matrix20(z12, Hy52_M20);
		double fd1 = 0.0;
		double fd2 = 0.0;
		for (int i = 0; i < Dim12; i++)
		{
			fd1 += z12[i] * z12[i];
			fd2 += z12[i];
		}
		f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


		//Rastrigin function
		double Hy53_M30[30][30];
		for (int i = 40; i < 70; i++)
		{
			for (int j = 40; j < 70; j++)
				Hy53_M30[i - 40][j - 40] = Hy5_M100[i][j];
		}

		double z13[30];
		for (int i = 0; i < Dim13; i++)
			z13[i] = X[i + 40] - Hy5_O[i + 40];
		mul_matrix30(z13, Hy53_M30);
		for (int i = 0; i < Dim13; i++)
			f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

		//Rosenbrock function
		double Hy54_M30[30][30];
		for (int i = 70; i < 100; i++)
		{
			for (int j = 70; j < 100; j++)
				Hy54_M3[i - 70][j - 70] = Hy5_M100[i][j];
		}

		double z14[30];
		for (int i = 0; i < Dim14; i++)
			z14[i] = X[i + 70] - Hy5_O[i + 70];
		mul_matrix30(z14, Hy54_M30);

		for (int i = 0; i < Dim14 - 1; i++)
			f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


		f1 = f11 + f12 + f13 + f14;


		//hybrid function 6
		//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
		double p1[4] = { 0.2,0.2,0.3,0.3 };
		int Dim21;
		int Dim22;
		int Dim23;
		int Dim24;

		Dim21 = DIM / 10 * 2;//6Î¬
		Dim22 = DIM / 10 * 2;//6Î¬
		Dim23 = DIM / 10 * 3;//9Î¬
		Dim24 = DIM / 10 * 3;//9Î¬

		double f21 = 0;
		double f22 = 0;
		double f23 = 0;
		double f24 = 0;
		double f2 = 0;

		extern double Hy6_O[100];
		extern double Hy6_M100[100][100];

		//Expanded Schaffer F6
		double Hy61_M20[20][20];
		for (int i = 0; i < 20; i++)
		{
			for (int j = 0; j < 20; j++)
				Hy61_M20[i][j] = Hy2_M100[i][j];
		}
		double z21[20];
		for (int i = 0; i < Dim21; i++)
			z21[i] = X[i] - Hy6_O[i];
		mul_matrix20(z21, Hy61_M20);
		double g = 0.0;
		for (int i = 0; i < Dim21; i++)
		{
			g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
			f21 += g;
		}
		f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



		//HGBat
		double Hy62_M20[20][20];
		for (int i = 20; i < 40; i++)
		{
			for (int j = 20; j < 40; j++)
				Hy62_M6[i - 20][j - 20] = Hy6_M10[i][j];
		}

		double z22[20];
		for (int i = 0; i < 20; i++)
			z22[i] = X[i + 20] - Hy6_O[i + 20];
		mul_matrix20(z22, Hy62_M20);
		double ff1 = 0.0;
		double ff2 = 0.0;
		for (int i = 0; i < Dim22; i++)
		{
			ff1 += z22[i] * z22[i];
			ff2 += z22[i];
		}
		f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



		//Rosenbrock
		double Hy63_M30[30][30];
		for (int i = 40; i < 70; i++)
		{
			for (int j = 40; j < 70; j++)
				Hy63_M3[i - 40][j - 40] = Hy6_M100[i][j];
		}

		double z23[30];
		for (int i = 0; i < Dim23; i++)
			z23[i] = X[i + 40] - Hy6_O[i + 40];
		mul_matrix30(z23, Hy63_M30);
		for (int i = 0; i < Dim23 - 1; i++)
			f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


		//Modified Schwefel F10
		double Hy64_M30[30][30];
		for (int i = 70; i < 100; i++)
		{
			for (int j = 70; j < 100; j++)
				Hy64_M9[i - 70][j - 70] = Hy6_M100[i][j];
		}

		double z24[30];
		for (int i = 0; i < Dim24; i++)
			z24[i] = X[i + 70] - Hy6_O[i + 70];
		mul_matrix30(z24, Hy64_M30);

		for (int i = 0; i < Dim24; i++)
		{
			z24[i] = z24[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z24[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z24[i]) <= 500)
				f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
			if (z24[i] > 500)
				f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
			if (z24[i] < -500)
				f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
		}
		f24 = 418.9829*Dim24 - f24;

		f2 = f21 + f22 + f23 + f24;


		// hybrid_function_7
		//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
		double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
		int Dim1 = DIM / 10 * 1;
		int Dim2 = DIM / 10 * 2;
		int Dim3 = DIM / 10 * 2;
		int Dim4 = DIM / 10 * 2;
		int Dim5 = DIM / 10 * 3;
		double f31 = 0;
		double f32 = 0;
		double f33 = 0;
		double f34 = 0;
		double f35 = 0;
		double f3 = 0;

		extern double Hy7_O[100];
		extern double Hy7_M100[100][100];
		double Hy71_M1[10][10];
		for (int i = 0; i < 10; i++)
			for (int j = 0; j < 10; j++)
				Hy71_M10[i][j] = Hy7_M100[i][j];
		double Hy72_M20[20][20];
		for (int i = 10; i < 30; i++)
			for (int j = 10; j < 30; j++)
				Hy72_M20[i - 10][j - 10] = Hy7_M100[i][j];
		double Hy73_M20[20][20];
		for (int i = 30; i < 50; i++)
			for (int j = 30; j < 50; j++)
				Hy73_M20[i - 30][j - 30] = Hy7_M100[i][j];
		double Hy74_M20[20][20];
		for (int i = 50; i < 70; i++)
			for (int j = 50; j < 70; j++)
				Hy74_M2[i - 50][j - 50] = Hy7_M100[i][j];
		double Hy75_M30[30][30];
		for (int i = 70; i < 100; i++)
			for (int j = 70; j < 100; j++)
				Hy75_M30[i - 70][j - 70] = Hy7_M100[i][j];

		double z1[10];
		for (int i = 0; i < 10; i++)
			z1[i] = X[i] - Hy7_O[i];
		mul_matrix10(z1, Hy71_M10);

		double yyy = 1.0;
		double yy1 = 0.0;
		double yy2 = 0.0;
		double yy3 = 0.0;
		yy1 = 10 / (Dim1*Dim1);
		yy3 = 10 / pow(Dim1, 1.2);
		for (int i = 0; i < Dim1; i++)
		{
			for (int j = 1; j <= 32; j++)
				yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
			yyy *= pow((1 + (i + 1)*yy2), yy3);
		}
		f31 = yy1 * yyy - yy1;


		double z2[20];
		for (int i = 10; i < 30; i++)
			z2[i - 10] = X[i] - Hy7_O[i];
		mul_matrix10(z2, Hy72_M10);

		double ffs1 = 0, ffs2 = 0;
		for (int i = 0; i < Dim2; i++) {
			ffs1 += z2[i] * z2[i];
			ffs2 += cos(2 * pai * z2[i]);
		}
		f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

		double z3[20];
		for (int i = 30; i < 50; i++)
			z3[i - 30] = X[i] - Hy7_O[i];
		mul_matrix20(z3, Hy73_M20);
		//Expanded Griewank's plus Rosenbrock
		double temp[20];
		for (int i = 0; i < Dim3 - 1; i++)
			temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
		temp[Dim3 - 1] = 100 * (z3[0] * z3[0] - z3[1])*(z3[0] * z3[0] - z3[1]) + (z3[0] - 1)*(z3[0] - 1);

		double u0 = 2.5;
		double d = 1;
		double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
		double u1 = -pow((u0*u0 - d) / s, 0.5);
		double sign[20];
		double F1 = 0;
		double F2 = 0;
		double ff = 0;
		double v[20][20];
		for (int i = 0; i < 20; i++)
		{
			for (int j = 0; j < 20; j++)
				v[i][j] = 0;
			v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
		}
		double z[20];
		double yy[20];
		double xz[20];//x»¡
		for (int i = 0; i < Dim3; i++)
			yy[i] = 0.1*(temp[i] - Hy7_O[i + 30]);

		for (int i = 0; i < Dim3; i++)
		{
			if (temp[i] < 0)
				sign[i] = -1;
			else if (temp[i] == 0)
				sign[i] = 0;
			else sign[i] = 1;
			xz[i] = 2 * sign[i] * yy[i] + u0;

		}

		for (int i = 0; i < Dim3; i++)
			z[i] = v[i][i] * (xz[i] - u0);
		for (int i = 0; i < Dim3; i++)
		{
			F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
			F2 = s * (xz[i] - u1)*(xz[i] - u1);
		}
		F2 = F2 + d * Dim3;
		if (F1 < F2)
			f33 = F1;
		else f33 = F2;

		for (int i = 0; i < Dim3; i++)
			ff = ff + cos(2 * pai*z[i]);
		f33 = f33 + 10 * (Dim3 - ff);


		double z4[20];
		for (int i = 50; i < 70; i++)
			z4[i - 50] = X[i] - Hy7_O[i];
		mul_matrix20(z4, Hy74_M20);

		for (int i = 0; i < Dim4; i++)
		{
			z4[i] = z4[i] + 4.209687462275036E+2;
			double mod = fmod(fabs(z4[i]), 500);
			double fabsmod = fabs(mod);
			if (fabs(z4[i]) <= 500)
				f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
			if (z4[i] > 500)
				f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
			if (z4[i] < -500)
				f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
		}
		f34 = 418.9829*Dim4 - f34;


		double z5[30];
		for (int i = 70; i < 100; i++)
			z5[i - 70] = X[i] - Hy7_O[i];
		mul_matrix30(z5, Hy75_M30);
		//Rastrigin
		for (int i = 0; i < Dim5; i++)
			f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

		f3 = f31 + f32 + f33 + f34 + f35;
		double w1 = 0;
		double w2 = 0;
		double w3 = 0;
		double w = 0;
		w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
		w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
		w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

		w = w1 + w2 + w3;
		w1 = w1 / w;
		w2 = w2 / w;
		w3 = w3 / w;

		f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
	#endif
		return f + 3000;

	#endif
	*/

#if Composition_function_10
	double f = 0;
	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double R[3] = { 10,30,50 };
	double T[3] = { 1,1,1 };
	double bias[3] = { 0,100,200 };
	extern double C101_O[100];
	extern double C102_O[100];
	extern double C103_O[100];
	for (int i = 0; i < DIM; i++)
	{
		y1 += (X[i] - C101_O[i])*(X[i] - C101_O[i]);
		y2 += (X[i] - C102_O[i])*(X[i] - C102_O[i]);
		y3 += (X[i] - C103_O[i])*(X[i] - C103_O[i]);

	}


#if DIM==10
	//hybrid function 5
		//BentCigar+HGBat+Rastrigin+Rosenbrock
	double p[4] = { 0.2,0.2,0.3,0.3 };
	int Dim11;
	int Dim12;
	int Dim13;
	int Dim14;
	Dim11 = DIM / 10 * 2;
	Dim12 = DIM / 10 * 2;
	Dim13 = DIM / 10 * 3;
	Dim14 = DIM / 10 * 3;
	double f11 = 0;
	double f12 = 0;
	double f13 = 0;
	double f14 = 0;
	double f1 = 0;

	extern double Hy5_O[100];
	extern double Hy5_M10[10][10];

	//Bent Cigar function
	double Hy51_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy51_M2[i][j] = Hy5_M10[i][j];
	}
	double z11[2];
	for (int i = 0; i < Dim11; i++)
		z11[i] = X[i] - Hy5_O[i];
	mul_matrix02(z11, Hy51_M2);
	double yyw1 = z11[0] * z11[0];
	double yyw2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yyw2 += 1E6*z11[i] * z11[i];
	f11 = yyw1 + yyw2;


	//HGBat function
	double Hy52_M2[2][2];
	for (int i = 2; i < 4; i++)
	{
		for (int j = 2; j < 4; j++)
			Hy52_M2[i - 2][j - 2] = Hy5_M10[i][j];
	}

	double z12[2];
	for (int i = 0; i < 2; i++)
		z12[i] = X[i + 2] - Hy5_O[i + 2];
	mul_matrix02(z12, Hy52_M2);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim12; i++)
	{
		fd1 += z12[i] * z12[i];
		fd2 += z12[i];
	}
	f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


	//Rastrigin function 
	double Hy53_M3[3][3];
	for (int i = 4; i < 7; i++)
	{
		for (int j = 4; j < 7; j++)
			Hy53_M3[i - 4][j - 4] = Hy5_M10[i][j];
	}

	double z13[3];
	for (int i = 0; i < Dim13; i++)
		z13[i] = X[i + 4] - Hy5_O[i + 4];
	mul_matrix03(z13, Hy53_M3);
	for (int i = 0; i < Dim13; i++)
		f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

	//Rosenbrock function 
	double Hy54_M3[3][3];
	for (int i = 7; i < 10; i++)
	{
		for (int j = 7; j < 10; j++)
			Hy54_M3[i - 7][j - 7] = Hy5_M10[i][j];
	}

	double z14[3];
	for (int i = 0; i < Dim14; i++)
		z14[i] = X[i + 7] - Hy5_O[i + 7];
	mul_matrix03(z14, Hy54_M3);

	for (int i = 0; i < Dim14 - 1; i++)
		f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


	f1 = f11 + f12 + f13 + f14;


	//hybrid function 6
	//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
	double p1[4] = { 0.2,0.2,0.3,0.3 };
	int Dim21;
	int Dim22;
	int Dim23;
	int Dim24;

	Dim21 = DIM / 10 * 2;
	Dim22 = DIM / 10 * 2;
	Dim23 = DIM / 10 * 3;
	Dim24 = DIM / 10 * 3;

	double f21 = 0;
	double f22 = 0;
	double f23 = 0;
	double f24 = 0;
	double f2 = 0;

	extern double Hy6_O[100];
	extern double Hy6_M10[10][10];

	//Expanded Schaffer F6
	double Hy61_M2[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			Hy61_M2[i][j] = Hy6_M10[i][j];
	}
	double z21[2];
	for (int i = 0; i < Dim21; i++)
		z21[i] = X[i] - Hy6_O[i];
	mul_matrix02(z21, Hy61_M2);
	double g = 0.0;
	for (int i = 0; i < Dim21 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
		f21 += g;
	}
	f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



	//HGBat
	double Hy62_M2[2][2];
	for (int i = 2; i < 4; i++)
	{
		for (int j = 2; j < 4; j++)
			Hy62_M2[i - 2][j - 2] = Hy6_M10[i][j];
	}

	double z22[2];
	for (int i = 0; i < 2; i++)
		z22[i] = X[i + 2] - Hy6_O[i + 2];
	mul_matrix02(z22, Hy62_M2);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim22; i++)
	{
		ff1 += z22[i] * z22[i];
		ff2 += z22[i];
	}
	f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



	//Rosenbrock
	double Hy63_M3[3][3];
	for (int i = 4; i < 7; i++)
	{
		for (int j = 4; j < 7; j++)
			Hy63_M3[i - 4][j - 4] = Hy6_M10[i][j];
	}

	double z23[3];
	for (int i = 0; i < Dim23; i++)
		z23[i] = X[i + 4] - Hy6_O[i + 4];
	mul_matrix03(z23, Hy63_M3);
	for (int i = 0; i < Dim23 - 1; i++)
		f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


	//Modified Schwefel F10
	double Hy64_M3[3][3];
	for (int i = 7; i < 10; i++)
	{
		for (int j = 7; j < 10; j++)
			Hy64_M3[i - 7][j - 7] = Hy6_M10[i][j];
	}

	double z24[3];
	for (int i = 0; i < Dim24; i++)
		z24[i] = X[i + 7] - Hy6_O[i + 7];
	mul_matrix03(z24, Hy64_M3);

	for (int i = 0; i < Dim24; i++)
	{
		z24[i] = z24[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z24[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z24[i]) <= 500)
			f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
		if (z24[i] > 500)
			f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
		if (z24[i] < -500)
			f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
	}
	f24 = 418.9829*Dim24 - f24;

	f2 = f21 + f22 + f23 + f24;


	// hybrid_function_7
	//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
	double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
	int Dim1 = DIM / 10 * 1;
	int Dim2 = DIM / 10 * 2;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 3;
	double f31 = 0;
	double f32 = 0;
	double f33 = 0;
	double f34 = 0;
	double f35 = 0;
	double f3 = 0;

	extern double Hy7_O[100];
	extern double Hy7_M10[10][10];
	double Hy71_M1[1][1];
	for (int i = 0; i < 1; i++)
		for (int j = 0; j < 1; j++)
			Hy71_M1[i][j] = Hy7_M10[i][j];
	double Hy72_M2[2][2];
	for (int i = 1; i < 3; i++)
		for (int j = 1; j < 3; j++)
			Hy72_M2[i - 1][j - 1] = Hy7_M10[i][j];
	double Hy73_M2[2][2];
	for (int i = 3; i < 5; i++)
		for (int j = 3; j < 5; j++)
			Hy73_M2[i - 3][j - 3] = Hy7_M10[i][j];
	double Hy74_M2[2][2];
	for (int i = 5; i < 7; i++)
		for (int j = 5; j < 7; j++)
			Hy74_M2[i - 5][j - 5] = Hy7_M10[i][j];
	double Hy75_M3[3][3];
	for (int i = 7; i < 10; i++)
		for (int j = 7; j < 10; j++)
			Hy75_M3[i - 7][j - 7] = Hy7_M10[i][j];

	double z1[1];
	for (int i = 0; i < 1; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix01(z1, Hy71_M1);

	double yyy = 1.0;
	double yy1 = 0.0;
	double yy2 = 0.0;
	double yy3 = 0.0;
	yy1 = 10 / (Dim1*Dim1);
	yy3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		yyy *= pow((1 + (i + 1)*yy2), yy3);
	}
	f31 = yy1 * yyy - yy1;


	double z2[2];
	for (int i = 1; i < 3; i++)
		z2[i - 1] = X[i] - Hy7_O[i];
	mul_matrix02(z2, Hy72_M2);

	double ffs1 = 0, ffs2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ffs1 += z2[i] * z2[i];
		ffs2 += cos(2 * pai * z2[i]);
	}
	f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

	double z3[2];
	for (int i = 3; i < 5; i++)
		z3[i - 3] = X[i] - Hy7_O[i];
	mul_matrix02(z3, Hy73_M2);
	//Expanded Griewank's plus Rosenbrock
	double temp[2];
	for (int i = 0; i < Dim3 - 1; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	//temp[Dim3 - 2] = 100 * (p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1])*(p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1]) + (p3[Dim3 - 2] - 1)*(p3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[2];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[2];
	double yy[2];
	double xz[2];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 3]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f33 = F1;
	else f33 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f33 = f33 + 10 * (Dim3 - ff);


	double z4[2];
	for (int i = 5; i < 7; i++)
		z4[i - 5] = X[i] - Hy7_O[i];
	mul_matrix02(z4, Hy74_M2);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f34 = 418.9829*Dim4 - f34;


	double z5[3];
	for (int i = 7; i < 10; i++)
		z5[i - 7] = X[i] - Hy7_O[i];
	mul_matrix03(z5, Hy75_M3);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f3 = f31 + f32 + f33 + f34 + f35;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;
	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
#endif


#if DIM==30
	//hybrid function 5
		//BentCigar+HGBat+Rastrigin+Rosenbrock
	double p[4] = { 0.2,0.2,0.3,0.3 };
	int Dim11;
	int Dim12;
	int Dim13;
	int Dim14;
	Dim11 = DIM / 10 * 2;//6Î¬
	Dim12 = DIM / 10 * 2;//6Î¬
	Dim13 = DIM / 10 * 3;//12Î¬
	Dim14 = DIM / 10 * 3;//12Î¬
	double f11 = 0;
	double f12 = 0;
	double f13 = 0;
	double f14 = 0;
	double f1 = 0;

	extern double Hy5_O[100];
	extern double Hy5_M30[30][30];

	//Bent Cigar function
	double Hy51_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy51_M6[i][j] = Hy5_M30[i][j];
	}
	double z11[6];
	for (int i = 0; i < Dim11; i++)
		z11[i] = X[i] - Hy5_O[i];
	mul_matrix06(z11, Hy51_M6);
	double yyw1 = z11[0] * z11[0];
	double yyw2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yyw2 += 1E6*z11[i] * z11[i];
	f11 = yyw1 + yyw2;


	//HGBat function
	double Hy52_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy52_M6[i - 6][j - 6] = Hy5_M30[i][j];
	}

	double z12[6];
	for (int i = 0; i < 6; i++)
		z12[i] = X[i + 6] - Hy5_O[i + 6];
	mul_matrix06(z12, Hy52_M6);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim12; i++)
	{
		fd1 += z12[i] * z12[i];
		fd2 += z12[i];
	}
	f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


	//Rastrigin function 
	double Hy53_M9[9][9];
	for (int i = 12; i < 21; i++)
	{
		for (int j = 12; j < 21; j++)
			Hy53_M9[i - 12][j - 12] = Hy5_M30[i][j];
	}

	double z13[9];
	for (int i = 0; i < Dim13; i++)
		z13[i] = X[i + 12] - Hy5_O[i + 12];
	mul_matrix09(z13, Hy53_M9);
	for (int i = 0; i < Dim13; i++)
		f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

	//Rosenbrock function 
	double Hy54_M9[9][9];
	for (int i = 21; i < 30; i++)
	{
		for (int j = 21; j < 30; j++)
			Hy54_M9[i - 21][j - 21] = Hy5_M30[i][j];
	}

	double z14[9];
	for (int i = 0; i < Dim14; i++)
		z14[i] = X[i + 21] - Hy5_O[i + 21];
	mul_matrix09(z14, Hy54_M9);

	for (int i = 0; i < Dim14 - 1; i++)
		f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


	f1 = f11 + f12 + f13 + f14;
	//cout << "f1=" << f1 << endl;

	//hybrid function 6
	//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
	double p1[4] = { 0.2,0.2,0.3,0.3 };
	int Dim21;
	int Dim22;
	int Dim23;
	int Dim24;

	Dim21 = DIM / 10 * 2;//6Î¬
	Dim22 = DIM / 10 * 2;//6Î¬
	Dim23 = DIM / 10 * 3;//9Î¬
	Dim24 = DIM / 10 * 3;//9Î¬

	double f21 = 0;
	double f22 = 0;
	double f23 = 0;
	double f24 = 0;
	double f2 = 0;

	extern double Hy6_O[100];
	extern double Hy6_M30[30][30];

	//Expanded Schaffer F6
	double Hy61_M6[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			Hy61_M6[i][j] = Hy6_M30[i][j];
	}
	double z21[6];
	for (int i = 0; i < Dim21; i++)
		z21[i] = X[i] - Hy6_O[i];
	mul_matrix06(z21, Hy61_M6);
	double g = 0.0;
	for (int i = 0; i < Dim21 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
		f21 += g;
	}
	f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



	//HGBat
	double Hy62_M6[6][6];
	for (int i = 6; i < 12; i++)
	{
		for (int j = 6; j < 12; j++)
			Hy62_M6[i - 6][j - 6] = Hy6_M30[i][j];
	}

	double z22[6];
	for (int i = 0; i < 6; i++)
		z22[i] = X[i + 6] - Hy6_O[i + 6];
	mul_matrix06(z22, Hy62_M6);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim22; i++)
	{
		ff1 += z22[i] * z22[i];
		ff2 += z22[i];
	}
	f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



	//Rosenbrock
	double Hy63_M9[9][9];
	for (int i = 12; i < 21; i++)
	{
		for (int j = 12; j < 21; j++)
			Hy63_M9[i - 12][j - 12] = Hy6_M30[i][j];
	}

	double z23[9];
	for (int i = 0; i < Dim23; i++)
		z23[i] = X[i + 12] - Hy6_O[i + 12];
	mul_matrix09(z23, Hy63_M9);
	for (int i = 0; i < Dim23 - 1; i++)
		f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


	//Modified Schwefel F10
	double Hy64_M9[9][9];
	for (int i = 21; i < 30; i++)
	{
		for (int j = 21; j < 30; j++)
			Hy64_M9[i - 21][j - 21] = Hy6_M30[i][j];
	}

	double z24[9];
	for (int i = 0; i < Dim24; i++)
		z24[i] = X[i + 21] - Hy6_O[i + 21];
	mul_matrix09(z24, Hy64_M9);

	for (int i = 0; i < Dim24; i++)
	{
		z24[i] = z24[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z24[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z24[i]) <= 500)
			f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
		if (z24[i] > 500)
			f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
		if (z24[i] < -500)
			f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
	}
	f24 = 418.9829*Dim24 - f24;

	f2 = f21 + f22 + f23 + f24;

	//cout << "f2=" << f2 << endl;
	// hybrid_function_7
	//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
	double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
	int Dim1 = DIM / 10 * 1;
	int Dim2 = DIM / 10 * 2;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 3;
	double f31 = 0;
	double f32 = 0;
	double f33 = 0;
	double f34 = 0;
	double f35 = 0;
	double f3 = 0;

	extern double Hy7_O[100];
	extern double Hy7_M30[30][30];
	double Hy71_M3[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Hy71_M3[i][j] = Hy7_M30[i][j];
	double Hy72_M6[6][6];
	for (int i = 3; i < 9; i++)
		for (int j = 3; j < 9; j++)
			Hy72_M6[i - 3][j - 3] = Hy7_M30[i][j];
	double Hy73_M6[6][6];
	for (int i = 9; i < 15; i++)
		for (int j = 9; j < 15; j++)
			Hy73_M6[i - 9][j - 9] = Hy7_M30[i][j];
	double Hy74_M6[6][6];
	for (int i = 15; i < 21; i++)
		for (int j = 15; j < 21; j++)
			Hy74_M6[i - 15][j - 15] = Hy7_M30[i][j];
	double Hy75_M9[9][9];
	for (int i = 21; i < 30; i++)
		for (int j = 21; j < 30; j++)
			Hy75_M9[i - 21][j - 21] = Hy7_M30[i][j];

	double z1[3];
	for (int i = 0; i < 3; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix03(z1, Hy71_M3);

	double yyy = 1.0;
	double yy1 = 0.0;
	double yy2 = 0.0;
	double yy3 = 0.0;
	yy1 = 10 / (Dim1*Dim1);
	yy3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		yyy *= pow((1 + (i + 1)*yy2), yy3);
	}
	f31 = yy1 * yyy - yy1;


	double z2[6];
	for (int i = 3; i < 9; i++)
		z2[i - 3] = X[i] - Hy7_O[i];
	mul_matrix06(z2, Hy72_M6);

	double ffs1 = 0, ffs2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ffs1 += z2[i] * z2[i];
		ffs2 += cos(2 * pai * z2[i]);
	}
	f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

	double z3[6];
	for (int i = 9; i < 15; i++)
		z3[i - 9] = X[i] - Hy7_O[i];
	mul_matrix06(z3, Hy73_M6);
	//Expanded Griewank's plus Rosenbrock
	double temp[6];
	for (int i = 0; i < Dim3 - 1; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	//temp[Dim3 - 2] = 100 * (p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1])*(p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1]) + (p3[Dim3 - 2] - 1)*(p3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[6];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[6];
	double yy[6];
	double xz[6];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 9]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f33 = F1;
	else f33 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f33 = f33 + 10 * (Dim3 - ff);


	double z4[6];
	for (int i = 15; i < 21; i++)
		z4[i - 15] = X[i] - Hy7_O[i];
	mul_matrix06(z4, Hy74_M6);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f34 = 418.9829*Dim4 - f34;


	double z5[9];
	for (int i = 21; i < 30; i++)
		z5[i - 21] = X[i] - Hy7_O[i];
	mul_matrix09(z5, Hy75_M9);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f3 = f31 + f32 + f33 + f34 + f35;
	//cout << "f3=" << f3 << endl;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;
	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;
	//cout << "w1=" << w1 << endl;
	//cout << "w2=" << w2 << endl;
	//cout << "w3=" << w3 << endl;
	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
#endif

#if DIM==50
	//hybrid function 5
		//BentCigar+HGBat+Rastrigin+Rosenbrock
	double p[4] = { 0.2,0.2,0.3,0.3 };
	int Dim11;
	int Dim12;
	int Dim13;
	int Dim14;
	Dim11 = DIM / 10 * 2;//6Î¬
	Dim12 = DIM / 10 * 2;//6Î¬
	Dim13 = DIM / 10 * 3;//12Î¬
	Dim14 = DIM / 10 * 3;//12Î¬
	double f11 = 0;
	double f12 = 0;
	double f13 = 0;
	double f14 = 0;
	double f1 = 0;

	extern double Hy5_O[100];
	extern double Hy5_M50[50][50];

	//Bent Cigar function
	double Hy51_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy51_M10[i][j] = Hy5_M50[i][j];
	}
	double z11[10];
	for (int i = 0; i < Dim11; i++)
		z11[i] = X[i] - Hy5_O[i];
	mul_matrix10(z11, Hy51_M10);
	double yyw1 = z11[0] * z11[0];
	double yyw2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yyw2 += 1E6*z11[i] * z11[i];
	f11 = yyw1 + yyw2;


	//HGBat function
	double Hy52_M10[10][10];
	for (int i = 10; i < 20; i++)
	{
		for (int j = 10; j < 20; j++)
			Hy52_M10[i - 10][j - 10] = Hy5_M50[i][j];
	}

	double z12[10];
	for (int i = 0; i < 10; i++)
		z12[i] = X[i + 10] - Hy5_O[i + 10];
	mul_matrix10(z12, Hy52_M10);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim12; i++)
	{
		fd1 += z12[i] * z12[i];
		fd2 += z12[i];
	}
	f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


	//Rastrigin function 
	double Hy53_M15[15][15];
	for (int i = 20; i < 35; i++)
	{
		for (int j = 20; j < 35; j++)
			Hy53_M15[i - 20][j - 20] = Hy5_M50[i][j];
	}

	double z13[15];
	for (int i = 0; i < Dim13; i++)
		z13[i] = X[i + 20] - Hy5_O[i + 20];
	mul_matrix15(z13, Hy53_M15);
	for (int i = 0; i < Dim13; i++)
		f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

	//Rosenbrock function 
	double Hy54_M15[15][15];
	for (int i = 35; i < 50; i++)
	{
		for (int j = 35; j < 50; j++)
			Hy54_M15[i - 35][j - 35] = Hy5_M50[i][j];
	}

	double z14[15];
	for (int i = 0; i < Dim14; i++)
		z14[i] = X[i + 35] - Hy5_O[i + 35];
	mul_matrix15(z14, Hy54_M15);

	for (int i = 0; i < Dim14 - 1; i++)
		f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


	f1 = f11 + f12 + f13 + f14;


	//hybrid function 6
	//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
	double p1[4] = { 0.2,0.2,0.3,0.3 };
	int Dim21;
	int Dim22;
	int Dim23;
	int Dim24;

	Dim21 = DIM / 10 * 2;//6Î¬
	Dim22 = DIM / 10 * 2;//6Î¬
	Dim23 = DIM / 10 * 3;//9Î¬
	Dim24 = DIM / 10 * 3;//9Î¬

	double f21 = 0;
	double f22 = 0;
	double f23 = 0;
	double f24 = 0;
	double f2 = 0;

	extern double Hy6_O[100];
	extern double Hy6_M50[50][50];

	//Expanded Schaffer F6
	double Hy61_M10[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			Hy61_M10[i][j] = Hy6_M50[i][j];
	}
	double z21[10];
	for (int i = 0; i < Dim21; i++)
		z21[i] = X[i] - Hy6_O[i];
	mul_matrix10(z21, Hy61_M10);
	double g = 0.0;
	for (int i = 0; i < Dim21 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
		f21 += g;
	}
	f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



	//HGBat
	double Hy62_M10[10][10];
	for (int i = 10; i < 20; i++)
	{
		for (int j = 10; j < 20; j++)
			Hy62_M10[i - 10][j - 10] = Hy6_M50[i][j];
	}

	double z22[10];
	for (int i = 0; i < 10; i++)
		z22[i] = X[i + 10] - Hy6_O[i + 10];
	mul_matrix10(z22, Hy62_M10);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim22; i++)
	{
		ff1 += z22[i] * z22[i];
		ff2 += z22[i];
	}
	f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



	//Rosenbrock
	double Hy63_M15[15][15];
	for (int i = 20; i < 35; i++)
	{
		for (int j = 20; j < 35; j++)
			Hy63_M15[i - 20][j - 20] = Hy6_M50[i][j];
	}

	double z23[15];
	for (int i = 0; i < Dim23; i++)
		z23[i] = X[i + 20] - Hy6_O[i + 20];
	mul_matrix15(z23, Hy63_M15);
	for (int i = 0; i < Dim23 - 1; i++)
		f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


	//Modified Schwefel F10
	double Hy64_M15[15][15];
	for (int i = 35; i < 50; i++)
	{
		for (int j = 35; j < 50; j++)
			Hy64_M15[i - 35][j - 35] = Hy6_M50[i][j];
	}

	double z24[15];
	for (int i = 0; i < Dim24; i++)
		z24[i] = X[i + 35] - Hy6_O[i + 35];
	mul_matrix15(z24, Hy64_M15);

	for (int i = 0; i < Dim24; i++)
	{
		z24[i] = z24[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z24[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z24[i]) <= 500)
			f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
		if (z24[i] > 500)
			f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
		if (z24[i] < -500)
			f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
	}
	f24 = 418.9829*Dim24 - f24;

	f2 = f21 + f22 + f23 + f24;


	// hybrid_function_7
	//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
	double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
	int Dim1 = DIM / 10 * 1;
	int Dim2 = DIM / 10 * 2;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 3;
	double f31 = 0;
	double f32 = 0;
	double f33 = 0;
	double f34 = 0;
	double f35 = 0;
	double f3 = 0;

	extern double Hy7_O[100];
	extern double Hy7_M50[50][50];
	double Hy71_M5[5][5];
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			Hy71_M5[i][j] = Hy7_M50[i][j];
	double Hy72_M10[10][10];
	for (int i = 5; i < 15; i++)
		for (int j = 5; j < 15; j++)
			Hy72_M10[i - 5][j - 5] = Hy7_M50[i][j];
	double Hy73_M10[10][10];
	for (int i = 15; i < 25; i++)
		for (int j = 15; j < 25; j++)
			Hy73_M10[i - 15][j - 15] = Hy7_M50[i][j];
	double Hy74_M10[10][10];
	for (int i = 25; i < 35; i++)
		for (int j = 25; j < 35; j++)
			Hy74_M10[i - 25][j - 25] = Hy7_M50[i][j];
	double Hy75_M15[15][15];
	for (int i = 35; i < 50; i++)
		for (int j = 35; j < 50; j++)
			Hy75_M15[i - 35][j - 35] = Hy7_M50[i][j];

	double z1[5];
	for (int i = 0; i < 5; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix05(z1, Hy71_M5);

	double yyy = 1.0;
	double yy1 = 0.0;
	double yy2 = 0.0;
	double yy3 = 0.0;
	yy1 = 10 / (Dim1*Dim1);
	yy3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		yyy *= pow((1 + (i + 1)*yy2), yy3);
	}
	f31 = yy1 * yyy - yy1;


	double z2[10];
	for (int i = 5; i < 15; i++)
		z2[i - 5] = X[i] - Hy7_O[i];
	mul_matrix10(z2, Hy72_M10);

	double ffs1 = 0, ffs2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ffs1 += z2[i] * z2[i];
		ffs2 += cos(2 * pai * z2[i]);
	}
	f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

	double z3[10];
	for (int i = 15; i < 25; i++)
		z3[i - 15] = X[i] - Hy7_O[i];
	mul_matrix10(z3, Hy73_M10);
	//Expanded Griewank's plus Rosenbrock
	double temp[10];
	for (int i = 0; i < Dim3 - 1; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	//temp[Dim3 - 2] = 100 * (p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1])*(p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1]) + (p3[Dim3 - 2] - 1)*(p3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[10];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[10][10];
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[10];
	double yy[10];
	double xz[10];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 15]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f33 = F1;
	else f33 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f33 = f33 + 10 * (Dim3 - ff);


	double z4[10];
	for (int i = 25; i < 35; i++)
		z4[i - 25] = X[i] - Hy7_O[i];
	mul_matrix10(z4, Hy74_M10);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f34 = 418.9829*Dim4 - f34;


	double z5[15];
	for (int i = 35; i < 50; i++)
		z5[i - 35] = X[i] - Hy7_O[i];
	mul_matrix15(z5, Hy75_M15);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f3 = f31 + f32 + f33 + f34 + f35;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;
	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
#endif

#if DIM==100
	//hybrid function 5
		//BentCigar+HGBat+Rastrigin+Rosenbrock
	double p[4] = { 0.2,0.2,0.3,0.3 };
	int Dim11;
	int Dim12;
	int Dim13;
	int Dim14;
	Dim11 = DIM / 10 * 2;
	Dim12 = DIM / 10 * 2;
	Dim13 = DIM / 10 * 3;
	Dim14 = DIM / 10 * 3;
	double f11 = 0;
	double f12 = 0;
	double f13 = 0;
	double f14 = 0;
	double f1 = 0;

	extern double Hy5_O[100];
	extern double Hy5_M100[100][100];

	//Bent Cigar function
	double Hy51_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy51_M20[i][j] = Hy5_M100[i][j];
	}
	double z11[20];
	for (int i = 0; i < Dim11; i++)
		z11[i] = X[i] - Hy5_O[i];
	mul_matrix20(z11, Hy51_M20);
	double yyw1 = z11[0] * z11[0];
	double yyw2 = 0.0;
	for (int i = 1; i < Dim11; i++)
		yyw2 += 1E6*z11[i] * z11[i];
	f11 = yyw1 + yyw2;


	//HGBat function
	double Hy52_M20[20][20];
	for (int i = 20; i < 40; i++)
	{
		for (int j = 20; j < 40; j++)
			Hy52_M20[i - 20][j - 20] = Hy5_M100[i][j];
	}

	double z12[20];
	for (int i = 0; i < 20; i++)
		z12[i] = X[i + 20] - Hy5_O[i + 20];
	mul_matrix20(z12, Hy52_M20);
	double fd1 = 0.0;
	double fd2 = 0.0;
	for (int i = 0; i < Dim12; i++)
	{
		fd1 += z12[i] * z12[i];
		fd2 += z12[i];
	}
	f12 = pow(fabs(fd1*fd1 - fd2 * fd2), 0.5) + (0.5*fd1 + fd2) / Dim12 + 0.5;


	//Rastrigin function 
	double Hy53_M30[30][30];
	for (int i = 40; i < 70; i++)
	{
		for (int j = 40; j < 70; j++)
			Hy53_M30[i - 40][j - 40] = Hy5_M100[i][j];
	}

	double z13[30];
	for (int i = 0; i < Dim13; i++)
		z13[i] = X[i + 40] - Hy5_O[i + 40];
	mul_matrix30(z13, Hy53_M30);
	for (int i = 0; i < Dim13; i++)
		f13 = f13 + (z13[i] * z13[i] - 10 * cos(2 * pai*z13[i]) + 10);

	//Rosenbrock function 
	double Hy54_M30[30][30];
	for (int i = 70; i < 100; i++)
	{
		for (int j = 70; j < 100; j++)
			Hy54_M30[i - 70][j - 70] = Hy5_M100[i][j];
	}

	double z14[30];
	for (int i = 0; i < Dim14; i++)
		z14[i] = X[i + 70] - Hy5_O[i + 70];
	mul_matrix30(z14, Hy54_M30);

	for (int i = 0; i < Dim14 - 1; i++)
		f14 += 100 * (z14[i] * z14[i] - z14[i + 1]) * (z14[i] * z14[i] - z14[i + 1]) + (z14[i] - 1) * (z14[i] - 1);


	f1 = f11 + f12 + f13 + f14;


	//hybrid function 6
	//F(X)=Schaff F6+HGBat+Rosenbrock+Modified Schawefel's f10
	double p1[4] = { 0.2,0.2,0.3,0.3 };
	int Dim21;
	int Dim22;
	int Dim23;
	int Dim24;

	Dim21 = DIM / 10 * 2;
	Dim22 = DIM / 10 * 2;
	Dim23 = DIM / 10 * 3;
	Dim24 = DIM / 10 * 3;

	double f21 = 0;
	double f22 = 0;
	double f23 = 0;
	double f24 = 0;
	double f2 = 0;

	extern double Hy6_O[100];
	extern double Hy6_M100[100][100];

	//Expanded Schaffer F6
	double Hy61_M20[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			Hy61_M20[i][j] = Hy6_M100[i][j];
	}
	double z21[20];
	for (int i = 0; i < Dim21; i++)
		z21[i] = X[i] - Hy6_O[i];
	mul_matrix20(z21, Hy61_M20);
	double g = 0.0;
	for (int i = 0; i < Dim21 - 1; i++)
	{
		g = 0.5 + (pow(sin(pow(z21[i] * z21[i] + z21[i + 1] * z21[i + 1], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[i] * z21[i] + z21[i + 1] * z21[i + 1])), 2);
		f21 += g;
	}
	f21 += 0.5 + (pow(sin(pow(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0], 0.5)), 2) - 0.5) / pow((1 + 0.001*(z21[Dim21 - 1] * z21[Dim21 - 1] + z21[0] * z21[0])), 2);



	//HGBat
	double Hy62_M20[20][20];
	for (int i = 20; i < 40; i++)
	{
		for (int j = 20; j < 40; j++)
			Hy62_M20[i - 20][j - 20] = Hy6_M100[i][j];
	}

	double z22[20];
	for (int i = 0; i < 20; i++)
		z22[i] = X[i + 20] - Hy6_O[i + 20];
	mul_matrix20(z22, Hy62_M20);
	double ff1 = 0.0;
	double ff2 = 0.0;
	for (int i = 0; i < Dim22; i++)
	{
		ff1 += z22[i] * z22[i];
		ff2 += z22[i];
	}
	f22 = pow(fabs(ff1*ff1 - ff2 * ff2), 0.5) + (0.5*ff1 + ff2) / Dim22 + 0.5;



	//Rosenbrock
	double Hy63_M30[30][30];
	for (int i = 40; i < 70; i++)
	{
		for (int j = 40; j < 70; j++)
			Hy63_M30[i - 40][j - 40] = Hy6_M100[i][j];
	}

	double z23[30];
	for (int i = 0; i < Dim23; i++)
		z23[i] = X[i + 40] - Hy6_O[i + 40];
	mul_matrix30(z23, Hy63_M30);
	for (int i = 0; i < Dim23 - 1; i++)
		f23 += 100 * (z23[i] * z23[i] - z23[i + 1]) * (z23[i] * z23[i] - z23[i + 1]) + (z23[i] - 1) * (z23[i] - 1);


	//Modified Schwefel F10
	double Hy64_M30[30][30];
	for (int i = 70; i < 100; i++)
	{
		for (int j = 70; j < 100; j++)
			Hy64_M30[i - 70][j - 70] = Hy6_M100[i][j];
	}

	double z24[30];
	for (int i = 0; i < Dim24; i++)
		z24[i] = X[i + 70] - Hy6_O[i + 70];
	mul_matrix30(z24, Hy64_M30);

	for (int i = 0; i < Dim24; i++)
	{
		z24[i] = z24[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z24[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z24[i]) <= 500)
			f24 += z24[i] * sin(pow(fabs(z24[i]), 0.5));
		if (z24[i] > 500)
			f24 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] - 500)*(z24[i] - 500) / (10000 * Dim24);
		if (z24[i] < -500)
			f24 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z24[i] + 500)*(z24[i] + 500) / (10000 * Dim24);
	}
	f24 = 418.9829*Dim24 - f24;

	f2 = f21 + f22 + f23 + f24;


	// hybrid_function_7
	//F(X)=Katsuura + Ackley + Expanded Griewank's plus Rosenbrock + Modified Schwefel + Rastrigin
	double p3[5] = { 0.1,0.2,0.2,0.2,0.3 };
	int Dim1 = DIM / 10 * 1;
	int Dim2 = DIM / 10 * 2;
	int Dim3 = DIM / 10 * 2;
	int Dim4 = DIM / 10 * 2;
	int Dim5 = DIM / 10 * 3;
	double f31 = 0;
	double f32 = 0;
	double f33 = 0;
	double f34 = 0;
	double f35 = 0;
	double f3 = 0;

	extern double Hy7_O[100];
	extern double Hy7_M100[100][100];
	double Hy71_M10[10][10];
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Hy71_M10[i][j] = Hy7_M100[i][j];
	double Hy72_M20[20][20];
	for (int i = 10; i < 30; i++)
		for (int j = 10; j < 30; j++)
			Hy72_M20[i - 10][j - 10] = Hy7_M100[i][j];
	double Hy73_M20[20][20];
	for (int i = 30; i < 50; i++)
		for (int j = 30; j < 50; j++)
			Hy73_M20[i - 30][j - 30] = Hy7_M100[i][j];
	double Hy74_M20[20][20];
	for (int i = 50; i < 70; i++)
		for (int j = 50; j < 70; j++)
			Hy74_M20[i - 50][j - 50] = Hy7_M100[i][j];
	double Hy75_M30[30][30];
	for (int i = 70; i < 100; i++)
		for (int j = 70; j < 100; j++)
			Hy75_M30[i - 70][j - 70] = Hy7_M100[i][j];

	double z1[10];
	for (int i = 0; i < 10; i++)
		z1[i] = X[i] - Hy7_O[i];
	mul_matrix10(z1, Hy71_M10);

	double yyy = 1.0;
	double yy1 = 0.0;
	double yy2 = 0.0;
	double yy3 = 0.0;
	yy1 = 10 / (Dim1*Dim1);
	yy3 = 10 / pow(Dim1, 1.2);
	for (int i = 0; i < Dim1; i++)
	{
		for (int j = 1; j <= 32; j++)
			yy2 += (fabs(pow(2, j)*z1[i] - ceil(pow(2, j)*z1[i]))) / pow(2, j);
		yyy *= pow((1 + (i + 1)*yy2), yy3);
	}
	f31 = yy1 * yyy - yy1;


	double z2[20];
	for (int i = 10; i < 30; i++)
		z2[i - 10] = X[i] - Hy7_O[i];
	mul_matrix20(z2, Hy72_M20);

	double ffs1 = 0, ffs2 = 0;
	for (int i = 0; i < Dim2; i++) {
		ffs1 += z2[i] * z2[i];
		ffs2 += cos(2 * pai * z2[i]);
	}
	f32 = -20 * exp(-0.2*sqrt(ffs1 / Dim2)) - exp(ffs2 / Dim2) + exp(1) + 20;

	double z3[20];
	for (int i = 30; i < 50; i++)
		z3[i - 30] = X[i] - Hy7_O[i];
	mul_matrix20(z3, Hy73_M20);
	//Expanded Griewank's plus Rosenbrock
	double temp[20];
	for (int i = 0; i < Dim3 - 1; i++)
		temp[i] = 100 * (z3[i] * z3[i] - z3[i + 1])*(z3[i] * z3[i] - z3[i + 1]) + (z3[i] - 1)*(z3[i] - 1) + 100 * (z3[i + 1] * z3[i + 1] - z3[i + 2])*(z3[i + 1] * z3[i + 1] - z3[i + 2]) + (z3[i + 1] - 1)*(z3[i + 1] - 1);
	//temp[Dim3 - 2] = 100 * (p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1])*(p3[Dim3 - 2] * p3[Dim3 - 2] - p3[Dim3 - 1]) + (p3[Dim3 - 2] - 1)*(p3[Dim3 - 2] - 1);
	temp[Dim3 - 1] = 100 * (z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0])*(z3[Dim3 - 1] * z3[Dim3 - 1] - z3[0]) + (z3[Dim3 - 1] - 1)*(z3[Dim3 - 1] - 1);

	double u0 = 2.5;
	double d = 1;
	double s = 1 - 1 / (2 * pow(Dim3 + 20, 0.5) - 8.2);
	double u1 = -pow((u0*u0 - d) / s, 0.5);
	double sign[20];
	double F1 = 0;
	double F2 = 0;
	double ff = 0;
	double v[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
			v[i][j] = 0;
		v[i][i] = pow(100, i / (2 * (Dim3 - 1)));
	}
	double z[20];
	double yy[20];
	double xz[20];//x»¡
	for (int i = 0; i < Dim3; i++)
		yy[i] = 0.1*(temp[i] - Hy7_O[i + 30]);

	for (int i = 0; i < Dim3; i++)
	{
		if (temp[i] < 0)
			sign[i] = -1;
		else if (temp[i] == 0)
			sign[i] = 0;
		else sign[i] = 1;
		xz[i] = 2 * sign[i] * yy[i] + u0;

	}

	for (int i = 0; i < Dim3; i++)
		z[i] = v[i][i] * (xz[i] - u0);
	for (int i = 0; i < Dim3; i++)
	{
		F1 = F1 + (xz[i] - u0)*(xz[i] - u0);
		F2 = s * (xz[i] - u1)*(xz[i] - u1);
	}
	F2 = F2 + d * Dim3;
	if (F1 < F2)
		f33 = F1;
	else f33 = F2;

	for (int i = 0; i < Dim3; i++)
		ff = ff + cos(2 * pai*z[i]);
	f33 = f33 + 10 * (Dim3 - ff);


	double z4[20];
	for (int i = 50; i < 70; i++)
		z4[i - 50] = X[i] - Hy7_O[i];
	mul_matrix20(z4, Hy74_M20);

	for (int i = 0; i < Dim4; i++)
	{
		z4[i] = z4[i] + 4.209687462275036E+2;
		double mod = fmod(fabs(z4[i]), 500);
		double fabsmod = fabs(mod);
		if (fabs(z4[i]) <= 500)
			f34 += z4[i] * sin(pow(fabs(z4[i]), 0.5));
		if (z4[i] > 500)
			f34 += (500 - mod)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] - 500)*(z4[i] - 500) / (10000 * Dim4);
		if (z4[i] < -500)
			f34 += (mod - 500)*sin(pow((500 - fabsmod), 0.5)) - (z4[i] + 500)*(z4[i] + 500) / (10000 * Dim4);
	}
	f34 = 418.9829*Dim4 - f34;


	double z5[30];
	for (int i = 70; i < 100; i++)
		z5[i - 70] = X[i] - Hy7_O[i];
	mul_matrix30(z5, Hy75_M30);
	//Rastrigin
	for (int i = 0; i < Dim5; i++)
		f35 = f35 + (z5[i] * z5[i] - 10 * cos(2 * pai*z5[i]) + 10);

	f3 = f31 + f32 + f33 + f34 + f35;
	double w1 = 0;
	double w2 = 0;
	double w3 = 0;
	double w = 0;
	w1 = 1 / pow(y1, 0.5)*exp(-y1 / (2 * DIM*R[0] * R[0]));
	w2 = 1 / pow(y2, 0.5)*exp(-y2 / (2 * DIM*R[1] * R[1]));
	w3 = 1 / pow(y3, 0.5)*exp(-y3 / (2 * DIM*R[2] * R[2]));

	w = w1 + w2 + w3;
	w1 = w1 / w;
	w2 = w2 / w;
	w3 = w3 / w;

	f = w1 * (T[0] * f1 + bias[0]) + w2 * (T[1] * f2 + bias[1]) + w3 * (T[2] * f3 + bias[2]);
#endif
	return f + 3000;

#endif


}


double *ComputGradient(double X[])
{
	//CEC2017µÄº¯Êý
	//Unimodal Function
	//shifted and rotated Bent Cigar

#if shifted_rotated_BentCigar
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;
	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if shifted_rotated_sum_of_Different_Power
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if shifted_rotated_Zakharov
	double *gradient = new double[DIM];
	double  newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if shifted_rotated_Rosenbrock
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
		for (int i = 0; i < DIM; i++)
		{
			gradient[i] = newgradient;
		}
	return gradient;
#endif


#if shifted_rotated_Rastrigin
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if shifted_Rotated_Scaffer7
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if shifted_Rotated_Lunacek_Bi_Rastrigin
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if shifted_Rotated_Non_Continuous_Rastrigin
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if Shifted_Rotated_Levy
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if Shifted_Rotated_Schwefel
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if hybrid_function_1
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif


#if hybrid_function_2
	double *gradient = new double[DIM];
	double newgradient = 0.0;

	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
	return gradient;
#endif



#if hybrid_function_3
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
	return gradient;
#endif




#if hybrid_function_4
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if hybrid_function_5
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if hybrid_function_6
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif


#if hybrid_function_7
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if hybrid_function_8
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif

#if hybrid_function_9
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//	gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif

#if hybrid_function_10
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif




#if Composition_function_1
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif




#if Composition_function_2
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if Composition_function_3
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if Composition_function_4
	double *gradient = new double[DIM];
	double Xnew[DIM];
	double newgradient = 0.0;
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//	gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if Composition_function_5
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif



#if Composition_function_6
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//	gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif




#if Composition_function_7
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif







#if Composition_function_8
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;
		//cout<< gradient[i] <<",";
	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif





#if Composition_function_9
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif


#if Composition_function_10
	double *gradient = new double[DIM];
	double newgradient = 0.0;
	double Xnew[DIM];
	double fnew = 0;
	double fold = 0;
	for (int i = 0; i < DIM; i++)
		Xnew[i] = X[i];
	fold = Computafitness(X);

	for (int i = 0; i < DIM; i++)
	{
		Xnew[i] = step + X[i];
		fnew = Computafitness(Xnew);
		newgradient += (fnew - fold) / step;
		//gradient[i] = (fnew - fold) / step;

	}
	for (int i = 0; i < DIM; i++)
	{
		gradient[i] = newgradient;
	}
	return gradient;
#endif
}


//ÎÞ²ÎÊýµÄ¹¹Ôìº¯Êý£¬È·¶¨½â¿Õ¼äµÄ·¶Î§
GROUPSHEEP::GROUPSHEEP() {
	//Ã¿Î¬×ø±êµÄÉÏÏÞºÍÏÂÏÞ
	left = left_range;
	right = right_range;
}

void GROUPSHEEP::initofgroup() {					//³õÊ¼»¯ÖÖÈº
	generationTimes = 1;
	for (int i = 0; i < SNUM; i++) {				//³õÊ¼»¯ÑòÈº
		for (int j = 0; j < DIM; j++) {				//³õÊ¼»¯Ã¿Î¬×ø±êµÄÖµ
			sheep[i].coordinate[j] = rand() / (double)RAND_MAX*(right - left) + left;	//Ëæ»ú¸ø³öÃ¿Ò»Ö»ÑòµÄ×ø±ê
		}
		sheep[i].fitness = Computafitness(sheep[i].coordinate);							//¼ÆËãµÚÒ»´úÑòÈºÃ¿¸ö½âµÄÊÊ¶ÈÖµ
		sheep[i].number = i + 1;					//¸øÖÖÈºµÄÃ¿¸ö¸öÌå±àºÅ£¬±àºÅ±ÈÊý×éÏÂ±ê´ó1
	}
}

void GROUPSHEEP::leader() {
	//cout << "leader" << endl;
	GradientDescent();
	double loc[DIM], res;//´æ´¢ÏòÁìÍ·Ñò¿¿½üºóµÄÁÙÊ±×ø±êºÍÁÙÊ±ÊÊÓ¦¶ÈÖµ
	for (int i = 0; i < SNUM; i++) {
		sheep[i].scatter = 0;
		for (int j = 0; j < DIM; j++) {
			loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[bellwethernumber - 1].coordinate[j] - sheep[i].coordinate[j]);
			if (loc[j] < left)
				loc[j] = right - left + loc[j];
			else if (loc[j] > right)
				loc[j] = loc[j] - right + left;
		}
		res = Computafitness(loc);
		if (res < sheep[i].fitness) {					//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï£¬È¡×îÓÅÎªÐÂ
			for (int k = 0; k < DIM; k++) {
				sheep[i].coordinate[k] = loc[k];
			}
			sheep[i].fitness = res;
		}
	}
	generationTimes++;
}

void GROUPSHEEP::wander() {
	//cout << "wander" << endl;
	int i, j, random;
	double loc[DIM], res;
	for (i = 0; i < SNUM; i++) //Ã¿Ö»ÑòÏòÆäËûÑòÑ§Ï°Ò»´Î
	{
		for (j = 0; j < DIM; j++)
		{
			do {
				random = (int)(rand() / (double)RAND_MAX*SNUM);//Ëæ»úÉú³É±»Ñ§Ï°µÄÑòµÄ±àºÅ[0..39]
			} while (random == i || random == SNUM);	//±»Ñ§Ï°µÄÑòÊÇÈÎÒâ·Ç¼ºµÄ
			if (sheep[i].fitness < sheep[random].fitness)
			{
				loc[j] = sheep[random].coordinate[j];
				sheep[random].coordinate[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				if (sheep[random].coordinate[j] < left)					//´¦Àí³¬¹ý½â¿Õ¼ä
					sheep[random].coordinate[j] = right - left + sheep[random].coordinate[j];
				else if (sheep[random].coordinate[j] > right)
					sheep[random].coordinate[j] = sheep[random].coordinate[j] - right + left;
				res = Computafitness(sheep[random].coordinate);
				if (res < sheep[random].fitness) 		//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
					sheep[random].fitness = res;
				else
					sheep[random].coordinate[j] = loc[j];

				loc[j] = sheep[i].coordinate[j];
				sheep[i].coordinate[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				if (sheep[i].coordinate[j] < left)					//´¦Àí³¬¹ý½â¿Õ¼ä
					sheep[i].coordinate[j] = right - left + sheep[i].coordinate[j];
				else if (sheep[i].coordinate[j] > right)
					sheep[i].coordinate[j] = sheep[i].coordinate[j] - right + left;
				res = Computafitness(sheep[i].coordinate);
				if (res < sheep[i].fitness) 		//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
					sheep[i].fitness = res;
				else
					sheep[i].coordinate[j] = loc[j];
			}
			else
			{
				loc[j] = sheep[i].coordinate[j];
				sheep[i].coordinate[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
				if (sheep[i].coordinate[j] < left)					//´¦Àí³¬¹ý½â¿Õ¼ä
					sheep[i].coordinate[j] = right - left + sheep[i].coordinate[j];
				else if (sheep[i].coordinate[j] > right)
					sheep[i].coordinate[j] = sheep[i].coordinate[j] - right + left;
				res = Computafitness(sheep[i].coordinate);
				if (res < sheep[i].fitness) 		//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
					sheep[i].fitness = res;
				else
					sheep[i].coordinate[j] = loc[j];

				loc[j] = sheep[random].coordinate[j];
				sheep[random].coordinate[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				if (sheep[random].coordinate[j] < left)					//´¦Àí³¬¹ý½â¿Õ¼ä
					sheep[random].coordinate[j] = right - left + sheep[random].coordinate[j];
				else if (sheep[random].coordinate[j] > right)
					sheep[random].coordinate[j] = sheep[random].coordinate[j] - right + left;
				res = Computafitness(sheep[random].coordinate);
				if (res < sheep[random].fitness) 		//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
					sheep[random].fitness = res;
				else
					sheep[random].coordinate[j] = loc[j];
			}
		}
	}
}

void GROUPSHEEP::bellwether() {						//ÕÒ³öÏÂÒ»´Îµü´úµÄÁìÍ·Ñò
	double bellwetherfitness;
	double sum = 0;
	bellwetherfitness = sheep[0].fitness;
	bellwethernumber = 1;
	for (int k = 1; k < SNUM; k++) {
		if (sheep[k].fitness < bellwetherfitness) {
			bellwetherfitness = sheep[k].fitness;
			bellwethernumber = k + 1;				//´æÈëÁìÍ·ÑòµÄ±àºÅ
		}
	}
	for (int j = 0; j < DIM; j++)
		GBest[j] = sheep[bellwethernumber - 1].coordinate[j];
	//ÇóÆ½¾ùÊÊ¶ÈÖµ
	for (int i = 0; i < SNUM; i++)
		sum += sheep[i].fitness;
	meanfitness = sum / SNUM;
	//Çó×î²îÊÊ¶ÈÖµ
	worstfitness = sheep[0].fitness;
	for (int j = 1; j < SNUM; j++) {
		if (sheep[j].fitness > worstfitness) {
			worstfitness = sheep[j].fitness;
		}
	}
}


int GROUPSHEEP::huntaway()
{

	static int lastHunt = 0;//×îºóÒ»´ÎÄÁÑòÈ®½éÈë´úÊý
	if (fabs(sheep[bellwethernumber - 1].fitness - oldbellwetherfitness) >= U) return 0;
	if (generationTimes < lastHunt) lastHunt = 0;
	if (generationTimes - lastHunt <= TIMES_HUNTAWAY) return 0; //Èç¹û¾àÉÏ´ÎÄÁÑòÈ®½éÈëµü´ú´ÎÊýÐ¡ÓÚµÈÓÚÄÁÑòÈ®µü´ú¼ä¸ô£¬²»Ö´ÐÐÄÁÑòÈ®½éÈë

	int i, j, random, randomsheep;
	double loc[DIM], res;
	lastHunt = generationTimes;		//¸üÐÂ×îºóÒ»´ÎÄÁÑòÈ®½éÈë´úÊý
	int k_max = 20;             //ÁÚÓòÉèÖÃ×î´óÖµÎªSNUM/2
	//int number = 0;


	if (fabs(sheep[bellwethernumber - 1].fitness - oldbellwetherfitness) < U)
	{
		//cout << "huntaway" << endl;
		int k = 1;									//µÚ0ÁÚÓò£¬ÒýÈëÄÁÑòÈ®¼à¶½Ö®Ç°		
		while (k <= k_max) {							//µÚkÁÚÓò
			double bellwetherfitness1 = sheep[bellwethernumber - 1].fitness;		
			int isscatter = 0;               //³õÊ¼»¯ÑòÊÇ·ñ±»´òÉ¢£¬ÊÇ1·ñ0
			srand((int)time(0));

			for (i = 0; i < SNUM; i++) {
				if (i + 1 != bellwethernumber)			//±£³ÖÁìÍ·Ñò²»±ä
				{
					if (isscatter >= k)
						break;
					else {
						//double pp = k / SNUM;
						if (rand() / (double)RAND_MAX < 0.3)
						{
							if (sheep[i].scatter == 1)
								continue;

							else
							{
								sheep[i].scatter = 1;
								
								isscatter++;               //´òÉ¢µÄÑòµÄÊýÁ¿
								for (j = 0; j < DIM; j++)
								{
									sheep[i].coordinate[j] = rand() / (double)RAND_MAX*(right - left) + left;//Ëæ»ú¸ø³öÑòµÄ×ø±ê
									if (sheep[i].coordinate[j] < left)			//´¦ÀíÎ»ÖÃ³¬³ö½â¿Õ¼äµÄÇé¿ö
										sheep[i].coordinate[j] = right - left + sheep[i].coordinate[j];
									else if (sheep[i].coordinate[j] > right)
										sheep[i].coordinate[j] = sheep[i].coordinate[j] - right + left;
								}
								sheep[i].fitness = Computafitness(sheep[i].coordinate);  //µÃ³öÑòµÄÊÊÓ¦¶ÈÖµ
							}
							i += int(rand() / (double)RAND_MAX);
						}
					}
				}
			}
			int num = 0;
			if (k > isscatter) {
				for (i = SNUM; i >= 0; i--)
				{
					if (i + 1 != bellwethernumber)			//±£³ÖÁìÍ·Ñò²»±ä
					{
						if (num < k - isscatter)
						{
							if (sheep[i].scatter == 0)
							{
								sheep[i].scatter = 1;
								for (j = 0; j < DIM; j++)
								{
									sheep[i].coordinate[j] = rand() / (double)RAND_MAX*(right - left) + left;//Ëæ»ú¸ø³öÑòµÄ×ø±ê
									if (sheep[i].coordinate[j] < left)			//´¦ÀíÎ»ÖÃ³¬³ö½â¿Õ¼äµÄÇé¿ö
										sheep[i].coordinate[j] = right - left + sheep[i].coordinate[j];
									else if (sheep[i].coordinate[j] > right)
										sheep[i].coordinate[j] = sheep[i].coordinate[j] - right + left;
								}
								sheep[i].fitness = Computafitness(sheep[i].coordinate);  //µÃ³öÑòµÄÊÊÓ¦¶ÈÖµ
								num++;
							}
						}
					}
				}
			}
			isscatter = isscatter + num;

			for (i = 0; i < SNUM; i++)
			{
				//Ã»ÓÐ±»ÇýÉ¢µÄÑò
				if (isscatter == 0)
					return 0;
				if (sheep[i].scatter != 1)
				{
					//Ã»ÓÐ±»´òÉ¢µÄÑòiÑ¡ÔñÒ»Ö»±»´òÉ¢µÄÑòrandom½øÐÐ¿¿½ü
					while (1)
					{
						random = (int)(rand() / (double)RAND_MAX*(SNUM));
						if (sheep[random].scatter == 1)
							break;

					}
					if (sheep[i].fitness <= sheep[random].fitness)
					{
						
						//Èç¹ûÒªÑ§Ï°µÄÑò±È×Ô¼ºµÄÊÊ¶ÈÖµ²î»òÕßÒ»Ñù£¬ÄÇÃ´ÈÃËü¿¿½üÎÒ
						for (j = 0; j < DIM; j++)
						{
							loc[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
							//printf("%.6lf", rand() / (double)RAND_MAX);
							if (loc[j] < left)			//´¦ÀíÎ»ÖÃ³¬³ö½â¿Õ¼äµÄÇé¿ö
								loc[j] = right - left + loc[j];
							else if (loc[j] > right)
								loc[j] = loc[j] - right + left;
						}
						res = Computafitness(loc);		//¼ÆËãÊÊ¶ÈÖµ
						if (res < sheep[random].fitness)	//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
						{
							for (j = 0; j < DIM; j++)
							{
								sheep[random].coordinate[j] = loc[j];
							}

							sheep[random].fitness = res;
						}

						//ÎÒÔ¶ÀëËü
						for (j = 0; j < DIM; j++)
						{
							loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
							if (loc[j] < left)					//´¦Àí³¬¹ý½â¿Õ¼ä
								loc[j] = right - left + loc[j];
							else if (loc[j] > right)
								loc[j] = loc[j] - right + left;
						}
						res = Computafitness(loc);
						if (res < sheep[i].fitness)			//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
						{
							for (j = 0; j < DIM; j++)
							{
								sheep[i].coordinate[j] = loc[j];
							}
							sheep[i].fitness = res;
						}

					}
					else	//Èç¹ûÒªÑ§Ï°µÄÑò±È×Ô¼ºµÄÊÊ¶ÈÖµÒªºÃ£¬ÄÇÃ´¿¿½üËü
					{
						for (j = 0; j < DIM; j++)
						{
							loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
							//printf("%.6lf", rand() / (double)RAND_MAX);
							if (loc[j] < left)			//´¦ÀíÎ»ÖÃ³¬³ö½â¿Õ¼äµÄÇé¿ö
								loc[j] = right - left + loc[j];
							else if (loc[j] > right)
								loc[j] = loc[j] - right + left;
						}
						res = Computafitness(loc);	//¼ÆËãÊÊ¶ÈÖµ
						if (res < sheep[i].fitness)	//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
						{
							for (j = 0; j < DIM; j++)
							{
								sheep[i].coordinate[j] = loc[j];
							}
							sheep[i].fitness = res;
						}
						//ËüÔ¶ÀëÎÒ
						for (j = 0; j < DIM; j++)
						{
							loc[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
							if (loc[j] < left)					//´¦Àí³¬¹ý½â¿Õ¼ä
								loc[j] = right - left + loc[j];
							else
								if (loc[j] > right)
									loc[j] = loc[j] - right + left;
						}
						res = Computafitness(loc);
						if (res < sheep[random].fitness)			//ÐÂÎ»ÖÃÓë¾ÉÎ»ÖÃ½øÐÐ±È½Ï
						{
							for (j = 0; j < DIM; j++)
							{
								sheep[random].coordinate[j] = loc[j];
							}
							sheep[random].fitness = res;
						}

					}
				}
			}
			bellwether();

			double nowbellfitness = sheep[bellwethernumber - 1].fitness;
			if (nowbellfitness < bellwetherfitness1)
			{
				sheep[bellwethernumber - 1].fitness = nowbellfitness;
				bellwetherfitness1 = nowbellfitness;
				k = 1;     //ÔÚµÚkÁÚÓòµÃµ½¸üºÃµÄ½âºó·µ»ØµÚ1ÁÚÓòÖØÐÂ½øÐÐËÑË÷
				//cout << "back to k=1" << endl;

				for (int i = 0; i < SNUM; i++)
				{
					sheep[i].scatter = 0;
				}
			}
			else {
				k++;				//ÔÚµÚkÁÚÓòµÃµ½µÄ½â²»Èçµ±Ç°×îÓÅ½âºÃ£¬ÔòÌøÈëk+1ÁÚÓò½øÐÐËÑË÷
				for (int i = 0; i < SNUM; i++)
				{
					sheep[i].scatter = 0;
				}
			}

		}

	}
	return 1;
}

//ÌÝ¶ÈÏÂ½µ
void GROUPSHEEP::GradientDescent()
{
	double PreviousFitness, NowFitness;
	double Step;
	int count = 0;
	BeforeGD = Computafitness(GBest);
	while (true)
	{
		Step = 1.0 / (1.0 + exp(++count));
		//Step = 1E-4;
		//++count;
		//cout <<" count= "<< count << endl;
		PreviousFitness = Computafitness(GBest);
		double* Gradient = ComputGradient(GBest);
		for (int j = 0; j < DIM; j++)
		{
			GBest[j] -= Step * Gradient[j];
			//for (int j = 0; j < DIM; j++)
				//cout <<"GBest[j]="<< GBest[j] << endl;
			if (GBest[j] < left)			//´¦ÀíÎ»ÖÃ³¬³ö½â¿Õ¼äµÄÇé¿ö
				GBest[j] = rand() / (double)RAND_MAX * (right - left) + left;//new
			else if (GBest[j] > right)
				GBest[j] = rand() / (double)RAND_MAX * (right - left) + left;//new
		}
		delete[] Gradient;
		NowFitness = Computafitness(GBest);
		if (NowFitness < sheep[bellwethernumber - 1].fitness)
		{
			for (int j = 0; j < DIM; j++)
				sheep[bellwethernumber - 1].coordinate[j] = GBest[j];
			sheep[bellwethernumber - 1].fitness = NowFitness;
		}
		//µü´úÖÕÖ¹Ìõ¼þ
		if ((fabs(NowFitness - PreviousFitness) < Epsilon) || (NowFitness >= PreviousFitness))
		{
			break;
		}

	}
	AfterGD = NowFitness;

}
