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

vector<vector<double>> MultiplyMatrices(vector<vector<double>>& matrix1, vector<std::vector<double>>& matrix2); //计算矩阵相乘

void ReduceMatrix(std::vector<std::vector<double>>& Q, int n); //将矩阵Q只保留前n列

vector<double> LS_proplem_solving(vector<vector<double>>& A, vector<double>& b); //求解线性最小二乘问题

double VectorTwoNorm(std::vector<double>& x); //计算向量的二范数

//第四章的函数

vector<double> Jacobi_Iteration(vector<vector<double>>& A, vector<double>& b); // Jacobi迭代法

vector<double> GS_Iteration(vector<vector<double>>& A, vector<double>& b); // G-S迭代法

vector<double> SOR_Iteration(vector<vector<double>>& A, vector<double>& b, double omega); // SOR迭代法

void Iterations(double epsilon); // 用三种迭代法求解差分方程

int SOR_Performance(vector<vector<double>>& A, vector<double>& b, double omega); //返回对于某个omega，SOR迭代求解的迭代次数

double BisearchOmega(vector<vector<double>>& A, vector<double>& b); // 在(1,2)区间里寻找SOR迭代法的最佳omgea

vector<vector<double>> MatrixSubtraction(vector<vector<double>> x, vector<vector<double>> y); // 返回两个矩阵之差

void Jacobi_Iteration2(vector<vector<double>>& u); // Jacobi迭代法

void GS_Iteration2(vector<vector<double>>& u); // G-S迭代法

void SOR_Iteration2(vector<vector<double>>& u, double omega); // SOR迭代法

int SOR_Performance2(vector<vector<double>>& A, double omega); //返回对于某个omega，SOR迭代求解的迭代次数

double BisearchOmega2(vector<vector<double>>& A); // 在(1,2)区间里寻找SOR迭代法的最佳omgea

void Iterations2(int n); // 用三种迭代法求解二维差分方程

//第五章的函数

double dotProduct(vector<double>& v1, vector<double>& v2); //求两个向量的内积

void ConjugateGradient1(vector<vector<double>>& A, vector<double>& x, vector<double>& b); // 共轭梯度法求解方程组, 并把解向量整合进矩阵输出

void ConjugateGradient2(vector<vector<double>>& A, vector<double>& x, vector<double>& b); // 共轭梯度法求解方程组

std::vector<std::vector<double>> generateMatrixA(int n);

//第六章的函数

vector<vector<double>> equation_to_matrix(vector<double>& equation); // 把首一多项式转化为以其特征多项式的矩阵

double maxModulus(vector<double>& vec); // 找到向量的模最大分量

double powerMethod(vector<vector<double>>& matrix); // 用幂法求解矩阵的模最大特征值

vector<vector<double>> getSubMatrix(vector<vector<double>>& A, int startRow, int endRow, int startCol, int endCol); // 提取子矩阵

void setSubMatrix(vector<vector<double>>& A, vector<vector<double>>& submatrix, int startRow, int startCol); // 赋值子矩阵

void hessenberg(vector<vector<double>>& A); // 计算上Hessenberg分解

void doubleShiftQR(vector<vector<double>>& H); // 双重步位移的QR迭代

bool areEigenvaluesReal(vector<vector<double>>& A); // 判断特征值是否为实数

void zeroing(vector<vector<double>>& A, double u); // 将矩阵中满足条件的次对角元置零

bool isQuasi(vector<vector<double>>& A); // 判断二维矩阵块是否是准上三角阵的对角元

bool quasii(vector<vector<double>>& H, int m); // 判断矩阵块是否是拟上三角阵

int quasi(vector<vector<double>>& H); // 找到拟上三角阵H33的最大维数

bool isHessenberg(vector<vector<double>>& A); // 判断是否是不可约Hessenberg矩阵

int IrredHessenberg(vector<vector<double>>& H, int quasi); // 找到不可约Hessenberg矩阵H22的最大维数

vector<vector<double>> getHessenberg(vector<vector<double>>& H, int m); // 提取不可约Hessenberg矩阵H22

void setHessenberg(vector<vector<double>>& H, vector<vector<double>>& A, int m); // 把QR迭代后的H22赋值到矩阵的对应位置

void eigenvalues2D(vector<vector<double>>& A); // 求解二阶矩阵的复特征值

void prEigens(vector<vector<double>>& A); // 输出QR迭代后得到矩阵的特征值

int implicitQR(vector<vector<double>>& A); //隐式QR迭代, 返回迭代次数

//第七章的函数

vector<vector<double>> jacobiClassic(vector<vector<double>>& A, int p, int q); // 对矩阵A的第p/q的行/列做Jacobi方法, 返回Jacobi旋转矩阵J

double offDiagNorm(vector<vector<double>>& A); // 求矩阵A的非对角范数

bool passingThreshold(vector<vector<double>>& A, double delta); // 判断矩阵A是否过了delta的关

vector<vector<double>> thresholdJacobi(vector<vector<double>>& A); // 过关Jacobi方法

int reversals(vector<vector<double>>& A, double miu); // 计算变号数

double dichEigenvalue(vector<vector<double>>& A, int m); // 二分法求第m个特征值

vector<double> VectorAddition(vector<double>x, vector<double>y); //计算向量相加得到的向量

vector<double> inversePowerMethod(vector<vector<double>>& matrix, double lambda); // 反幂法求特征向量