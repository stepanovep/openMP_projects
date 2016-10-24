#include "stdafx.h"
#include <ctime>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace std;

const int N = 3000;

#define forn(i, n) for(int i = 0; i < n; ++i)
#define endl '\n'

class Vector {
public:
	int n;
	double *v;
	Vector (int size) {
		n = size;
		v = new double[size];
		fill(v, v+size, 0);
	}

	Vector (const Vector &a) {
		n = a.n;
		v = new double[n];
		forn(i, n)
			v[i] = a.v[i];
	}

	void init(double *a) {
		for (int i = 0; i < n; ++i)
			v[i] = a[i];
	}

	void print() {
		forn(i, n)
			cout << v[i] << " ";
		cout << endl;
	}
	
	inline double & operator [] (const int idx) {
		return v[idx];
	}
};

class Matrix {
public:
	int n;
	double *m;
	Matrix (int size)
		: n(size) {
		m = new double[size * size];
	}

	Matrix (Matrix &B) {
		n = B.n;
		m = new double[n*n];
		forn(i, n)
			forn(j, n)
				m[i*n + j] = B(i,j);
	}

	void init(double a[][N]) {
		forn(i, n) {
			forn(j, n)
				m[i*n + j] = a[i][j];
		}
	}

	inline double & operator () (int i, int j) { return m[i*n + j]; }
	inline double   operator () (int i, int j) const { return m[i*n + j]; }
};


ostream & operator << (ostream& out, const Vector &v) {
	int n = v.n;
	forn(i, n)
		out << v.v[i] << ' ';

	out << endl;
    return out;
}
ostream & operator << (ostream& out, const Matrix &M) {
	int n = M.n;
	forn(i, n) {
		forn(j, n)
			out << M.m[i*n+j] << ' ';
		out << endl;
	}
	out << endl;
	return out;
}

Vector operator + (const Vector &a, const Vector &b) {
	int n = a.n;
	Vector c(a);
	forn(i, n)
		c.v[i] += b.v[i];

	return c;
}
Vector operator * (Matrix &m, Vector &a) {
	int n = m.n;
	Vector res(n);
	forn(i, n) {
		forn(j, n) {
			//res.v[i] += m[i][j] * a[j];
			res.v[i] += m.m[i*n + j] * a[j];
		}
	}
	return res;
}
Matrix operator - (Matrix& A, Matrix& B) {
	Matrix C(A);
	int n = A.n;
	forn(i, n) {
		forn(j, n) {
			//C[i][j] -= B[i][j];
			C(i,j) -= B(i,j);
		}
	}
	return C;
}

double norm (Matrix & M) {
	double maxSum = 0;
	int n = M.n;
	forn(i, n) {
		double sumi = 0;
		forn(j, n) {
			sumi += fabs(M.m[i*n + j]);
		}
		maxSum = max( maxSum, sumi);
	}

	return maxSum;
}
double normParallel (Matrix & M) {
	double maxSum1 = 0, maxSum2 = 0;
	int n = M.n;

	#pragma omp parallel
		{
		#pragma omp sections
			{
			#pragma omp section
				{
					for (int i = 0; i < n / 2; ++i) {
						double sumi = 0;
						for (int j = 0; j < n; ++j)
							//sumi += fabs(M[i][j]);
							sumi += fabs(M(i,j));
						maxSum1 = max (sumi, maxSum1);
					}

				}
			#pragma omp section
				{
					for (int i = n/2; i < n; ++i) {
						double sumi = 0;
						for (int j = 0; j < n; ++j)
							//sumi += fabs(M[i][j]);
							sumi += fabs(M(i, j));
						maxSum2 = max(sumi, maxSum2);
					}
				}
			}
		}

	return max (maxSum1, maxSum2);
}


Vector MatrVecMult(const Matrix &M, const Vector &v) {
	int n = M.n;
	Vector res(n);
	forn(i, n) {
		forn(j, n) {
			res.v[i] += M.m[i*n+j] * v.v[j];
		}
	}
	return res;
}
Vector MatrVecMultParallel(const Matrix &M, const Vector & v) {
	int n = M.n;
	Vector res(n);
	omp_set_dynamic(1);

	/*
	omp_set_num_threads(2);
	#pragma omp parallel sections
	{

		#pragma omp section
		{
			for (int i = 0; i < n / 2; ++i)
				forn(j, n)
					res.v[i] += M.m[i*n + j] * v.v[j];
		}
		#pragma omp section
		{
			for (int i = n / 2; i < n; ++i)
				forn(j, n)
					res.v[i] += M.m[i*n + j] * v.v[j];
		}
	}*/

	
	omp_set_num_threads(4);
	#pragma omp parallel for default(none) shared(res, M, v)
	forn(i, n) {
		forn(j, n) {
			res.v[i] += M.m[i*n + j] * v.v[j];
		}
	}

	return res;
}

void init (Matrix &B, Vector &F, Vector &X) {

	// ----------- Matrix B -------------- //
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			if (i == j)
				B.m[i*N + j] = 0;
			else {
				if (j == i+1)
					B.m[i*N + j] = 0.25;
				else
					B.m[i*N + j] = B.m[i*N + j - 1] / 2;

				B.m[j*N + i] = B.m[i*N + j];
			}
		}
	}

	// ------- Vector F ---------- //
	for (int i = 0; i < N; ++i)
		F[i] = 1;

	// ------- Vector X ---------- //
	for (int i = 0; i < N; ++i)
		X[i] = 0;
}


int main () {


	ios::sync_with_stdio(false);
	cin.tie(nullptr); cout.tie(nullptr);
	setlocale(LC_ALL, "Russian");

	/*double tick = omp_get_wtick();
	printf("Точность таймера %lf\n", tick);*/
	/*Vector v(3), u(3);
	double a[] = {1,2,3};
	
	v.init(a);

	double mm[][3] = {{1,2,3},
					 {2,3,4},
					 {4,5,6}};
	Matrix m(N);
	m.init(mm);*/

	Matrix B(N);					// Матрица B
	Vector F(N), X(N);				// Вектора F,X,   X - начальное приближение.    X_k+1 = B * X_k + F
					
	double t1, t2;
	{
		init(B, F, X);				// Инициализация начальных данных
		//cout << B << F << X;
		double s_t = omp_get_wtime();
		for (int i = 0; i < 50; ++i) {
			X = MatrVecMult(B, X) + F;
			//cout << X << endl;
		}
		double e_t = omp_get_wtime();
		t1 = e_t - s_t;
	}


	{
		init(B, F, X);				// Инициализация начальных данных
		//cout << B << F << X;
		double s_t = omp_get_wtime();
		for (int i = 0; i < 50; ++i) {
			X = MatrVecMultParallel(B, X) + F;
			//cout << X << endl;
		}
		double e_t = omp_get_wtime();
		t2 = e_t - s_t;
	}

	cout << "Обычное вычисление  : " << t1 << endl << endl;
	cout << "Параллельное вычисление  : " << t2 << endl << endl;


	return 0;
}