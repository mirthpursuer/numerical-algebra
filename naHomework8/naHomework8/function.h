#pragma once
#include<iostream>
#include<vector>
using namespace std;

//���ڵ��Եĺ���

void PrintMatrix(vector<vector<double>>& matrix);//��ӡ����

void PrintVector(vector<double>& vector);//��ӡ����

// ��һ�µĺ���

void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>>& U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>>& U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void gauss_elim(vector<vector<double>>& A);//Gauss��ȥ��

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//ȫ��ԪGauss��ȥ��

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& b); //������ԪGauss��ȥ��ⷽ��

void vector_pb(vector<int>& u, vector<double>& b);//��������P*b����ѡ��

void vector_qb(vector<int>& v, vector<double>& b);//��������Q*b����ѡ��

void cholesky_decomp(vector<vector<double>>& A);//�Գ��������׼Cholesky�ֽ�

void transposeMatrix(vector<vector<double>>& inputMatrix, vector<vector<double>>& transposedMatrix);
//ʵ�־����ת�ò����ת�ú���¾���

void modified_cholesky_decomp(vector<vector<double>>& A);//�Ľ���ƽ������

void Hilbert_Matrix(vector<vector<double>>& A);

void matrix_DLT(vector<vector<double>>& A);//�������D*L^T����ѡ��

//�ڶ��µĺ���

vector<double> sign(vector<double>& w); //����w�ķ�������

int sign(double m); //m>0Ϊ1, m<0Ϊ-1, m=0Ϊ0

double InnerProduct(vector<double>& a, vector<double>& b);//���������ڻ�

double VectorInfinityNorm(vector<double>& vec); //���������������

double VectorOneNorm(vector<double>& vec); //����������һ����

vector<double> UnitVectorGenerating(vector<double>& vec, int n); //�����������Ӧ�±�ĵ�λ����

double MatrixOneNorm(int n, vector<vector<double>>& A);//��������һ����

double MatrixInfinityNorm(vector<vector<double>>& matrix); //�������������

vector<double> MatrixVectorMultiply(vector<vector<double>>& A, vector<double>& b);//�������A������b��˵õ�������

vector<double> VectorSubtraction(vector<double>x, vector<double>y); //������������õ�������

//�����µĺ���

void house(vector<double>& x, vector<double>& v, double& beta); //����Householder�任

void QRDecomposition(vector<vector<double>>& A, vector<double>& d); //����QR�ֽ�

vector<vector<double>> HouseholderMatrix(vector<vector<double>>& A, vector<double>& d, int k); //����Householder�任����

vector<double> QR_equation_solving(vector<vector<double>>& A, vector<double>& b); //����QR�ֽ���ⷽ����

void equation_generating1(vector<vector<double>>& A, vector<double>& b); //���ɵ�һ�µĵ�һ��������

void equation_generating2(vector<vector<double>>& A, vector<double>& x, vector<double>& b); //���ɵ�һ�µĵڶ���������

void equation_generating3(vector<vector<double>>& A, vector<double>& b); //���ɵ�һ�µĵ�����������

vector<vector<double>> MultiplyMatrices(vector<vector<double>>& matrix1, vector<vector<double>>& matrix2); //����������

void ReduceMatrix(vector<vector<double>>& Q, int n); //������Qֻ����ǰn��

vector<double> LS_proplem_solving(vector<vector<double>>& A, vector<double>& b); //���������С��������

double VectorTwoNorm(vector<double>& x); //���������Ķ�����

//�����µĺ���

vector<double> Jacobi_Iteration(vector<vector<double>>& A, vector<double>& b); // Jacobi������

vector<double> GS_Iteration(vector<vector<double>>& A, vector<double>& b); // G-S������

vector<double> SOR_Iteration(vector<vector<double>>& A, vector<double>& b, double omega); // SOR������

void Iterations(double epsilon); // �����ֵ���������ַ���

int SOR_Performance(vector<vector<double>>& A, vector<double>& b, double omega); //���ض���ĳ��omega��SOR�������ĵ�������

double BisearchOmega(vector<vector<double>>& A, vector<double>& b); // ��(1,2)������Ѱ��SOR�����������omgea

vector<vector<double>> MatrixSubtraction(vector<vector<double>> x, vector<vector<double>> y); // ������������֮��

vector<vector<double>> MatrixAddition(vector<vector<double>> x, vector<vector<double>> y); // ������������֮��

void Jacobi_Iteration2(vector<vector<double>>& u); // Jacobi������

void GS_Iteration2(vector<vector<double>>& u); // G-S������

void SOR_Iteration2(vector<vector<double>>& u, double omega); // SOR������

int SOR_Performance2(vector<vector<double>>& A, double omega); //���ض���ĳ��omega��SOR�������ĵ�������

double BisearchOmega2(vector<vector<double>>& A); // ��(1,2)������Ѱ��SOR�����������omgea

void Iterations2(int n); // �����ֵ���������ά��ַ���

//�����µĺ���

double dotProduct(vector<double>& v1, vector<double>& v2); //�������������ڻ�

void ConjugateGradient1(vector<vector<double>>& A, vector<double>& x, vector<double>& b); // �����ݶȷ���ⷽ����, ���ѽ��������Ͻ��������

void ConjugateGradient2(vector<vector<double>>& A, vector<double>& x, vector<double>& b); // �����ݶȷ���ⷽ����

vector<vector<double>> generateMatrixA(int n);

//�����µĺ���

vector<vector<double>> equation_to_matrix(vector<double>& equation); // ����һ����ʽת��Ϊ������������ʽ�ľ���

double maxModulus(vector<double>& vec); // �ҵ�������ģ������

double powerMethod(vector<vector<double>>& matrix); // ���ݷ��������ģ�������ֵ

vector<vector<double>> getSubMatrix(vector<vector<double>>& A, int startRow, int endRow, int startCol, int endCol); // ��ȡ�Ӿ���

void setSubMatrix(vector<vector<double>>& A, vector<vector<double>>& submatrix, int startRow, int startCol); // ��ֵ�Ӿ���

void hessenberg(vector<vector<double>>& A); // ������Hessenberg�ֽ�

void doubleShiftQR(vector<vector<double>>& H); // ˫�ز�λ�Ƶ�QR����

bool areEigenvaluesReal(vector<vector<double>>& A); // �ж�����ֵ�Ƿ�Ϊʵ��

void zeroing(vector<vector<double>>& A, double u); // �����������������ĴζԽ�Ԫ����

bool isQuasi(vector<vector<double>>& A); // �ж϶�ά������Ƿ���׼��������ĶԽ�Ԫ

bool quasii(vector<vector<double>>& H, int m); // �жϾ�����Ƿ�������������

int quasi(vector<vector<double>>& H); // �ҵ�����������H33�����ά��

bool isHessenberg(vector<vector<double>>& A); // �ж��Ƿ��ǲ���ԼHessenberg����

int IrredHessenberg(vector<vector<double>>& H, int quasi); // �ҵ�����ԼHessenberg����H22�����ά��

vector<vector<double>> getHessenberg(vector<vector<double>>& H, int m); // ��ȡ����ԼHessenberg����H22

void setHessenberg(vector<vector<double>>& H, vector<vector<double>>& A, int m); // ��QR�������H22��ֵ������Ķ�Ӧλ��

void eigenvalues2D(vector<vector<double>>& A); // �����׾���ĸ�����ֵ

void prEigens(vector<vector<double>>& A); // ���QR������õ����������ֵ

int implicitQR(vector<vector<double>>& A); //��ʽQR����, ���ص�������

//�����µĺ���

vector<vector<double>> jacobiClassic(vector<vector<double>>& A, int p, int q); // �Ծ���A�ĵ�p/q����/����Jacobi����, ����Jacobi��ת����J

double offDiagNorm(vector<vector<double>>& A); // �����A�ķǶԽǷ���

bool passingThreshold(vector<vector<double>>& A, double delta); // �жϾ���A�Ƿ����delta�Ĺ�

vector<vector<double>> thresholdJacobi(vector<vector<double>>& A); // ����Jacobi����

int reversals(vector<vector<double>>& A, double miu); // ��������

double dichEigenvalue(vector<vector<double>>& A, int m); // ���ַ����m������ֵ

vector<double> VectorAddition(vector<double>x, vector<double>y); //����������ӵõ�������

vector<double> inversePowerMethod(vector<vector<double>>& matrix, double lambda); // ���ݷ�����������

//�ڰ˴���ҵ�ĺ���

vector<vector<double>> intMatrix(int rows, int cols); // ������������

vector<vector<double>> randomMatrix(int rows, int cols); // �������������

vector<vector<double>> idMatrix(int n); // ���ɵ�λ��

void transpose(vector<vector<double>>& A); // ʵ�־����ת��

vector<vector<double>> Householder(vector<double>& v, double beta); // ����Household����

void bidiag(vector<vector<double>>& A, vector<vector<double>>& P, vector<vector<double>>& Q); // ʵ��m*n�׾���Ķ��Խǻ�

void zeroing2(vector<vector<double>>& A, double epsilon); //�����������������ĶԽ�Ԫ�ʹζԽ�Ԫ����

void pqq(vector<vector<double>>B, int& p, int& q); // �ҵ�����p����С��q

void givens(double a, double b, double& c, double& s); // ����õ�����a��b��cos��sinֵ

vector<vector<double>> Givens(int n, int p, int q, double c, double s); // ��cos��sinֵ���ɶ�Ӧ��givens�任����

void GivensSVD(vector<vector<double>>& B, vector<vector<double>>& B22, int p, int i, vector<vector<double>>& U); // ��B22�жԽ�ԪΪ0����Givens�任�Ѹ��д��0

void WilkSVD(vector<vector<double>>& B, vector<vector<double>>& U, vector<vector<double>>& V);

double absmax(vector<vector<double>> A); // ���ؾ������ֵ����Ԫ��

void singluar(vector<vector<double>> A); // ��С�����������ֵ