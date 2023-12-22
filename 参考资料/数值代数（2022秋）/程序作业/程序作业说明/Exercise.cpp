#include "Exercise.h"
clock_t start, ende;


void exercise_1()
{
	int N = 1; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);

	//初始化A和b
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 6;
		A[i + 1][i] = 8;
		A[i][i + 1] = 1;
		b[i] = 15;
	}
	A[N - 1][N - 1] = 8;
	b[1] = 7;
	b[N - 1] = 14;


}

void exercise_2_1()
{
	
}

void exercise_2_2()
{
	
}

void exercise_3_1()
{
	
}

void exercise_3_2()
{
	
}