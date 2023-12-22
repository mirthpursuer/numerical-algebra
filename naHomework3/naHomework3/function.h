#pragma once
#include<iostream>
#include<vector>
using namespace std;

//用于调试的函数

void PrintMatrix(vector<vector<double>>& matrix);//打印矩阵

void PrintVector(std::vector<double>& vector);//打印向量

// 第一章的函数

void forward_subs(vector<vector<double>>& L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>>& U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>>& U, vector<double>& y);//对角元为1的回代法

void gauss_elim(vector<vector<double>>& A);//Gauss消去法

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//全主元Gauss消去法

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法

std::vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& b); //用列主元Gauss消去求解方程

void vector_pb(vector<int>& u, vector<double>& b);//计算向量P*b【可选】

void vector_qb(vector<int>& v, vector<double>& b);//计算向量Q*b【可选】

void cholesky_decomp(vector<vector<double>>& A);//对称正定阵标准Cholesky分解

void transposeMatrix(std::vector<std::vector<double>>& inputMatrix, std::vector<std::vector<double>>& transposedMatrix);
//实现矩阵的转置并输出转置后的新矩阵

void modified_cholesky_decomp(vector<vector<double>>& A);//改进的平方根法

void Hilbert_Matrix(vector<vector<double>>& A);

void matrix_DLT(vector<vector<double>>& A);//计算矩阵D*L^T【可选】

//第二章的函数

std::vector<double> sign(std::vector<double>& w); //生成w的符号向量

double InnerProduct(std::vector<double>& a, std::vector<double>& b);//计算向量内积

double VectorInfinityNorm(std::vector<double>& vec); //计算向量的无穷范数

double VectorOneNorm(std::vector<double>& vec); //计算向量的一范数

std::vector<double> UnitVectorGenerating(std::vector<double>& vec, int n); //生成无穷范数对应下标的单位向量

double MatrixOneNorm(int n, std::vector<std::vector<double>>& A);//计算矩阵的一范数

double MatrixInfinityNorm(std::vector<std::vector<double>>& matrix); //计算矩阵的无穷范数

std::vector<double> MatrixVectorMultiply(std::vector<std::vector<double>>& A, std::vector<double>& b);//计算矩阵A和向量b相乘得到的向量

vector<double> VectorSubtraction(vector<double>x, vector<double>y); //计算向量相减得到的向量

//第三章的函数

void house(std::vector<double>& x, std::vector<double>& v, double& beta); //计算Householder变换

void QRDecomposition(std::vector<std::vector<double>>& A, std::vector<double>& d); //计算QR分解

std::vector<std::vector<double>> HouseholderMatrix(std::vector<std::vector<double>>& A, std::vector<double>& d, int k); //生成Householder变换矩阵

std::vector<double> QR_equation_solving(vector<vector<double>>& A, vector<double>& b); //基于QR分解求解方程组

void equation_generating1(vector<vector<double>>& A, vector<double>& b); //生成第一章的第一个方程组

void equation_generating2(vector<vector<double>>& A, vector<double>& x, vector<double>& b); //生成第一章的第二个方程组

void equation_generating3(vector<vector<double>>& A, vector<double>& b); //生成第一章的第三个方程组

std::vector<std::vector<double>> MultiplyMatrices(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2); //计算矩阵相乘

void ReduceMatrix(std::vector<std::vector<double>>& Q, int n); //将矩阵Q只保留前n列

vector<double> LS_proplem_solving(vector<vector<double>>& A, vector<double>& b); //求解线性最小二乘问题

double VectorTwoNorm(std::vector<double>& x); //计算向量的二范数