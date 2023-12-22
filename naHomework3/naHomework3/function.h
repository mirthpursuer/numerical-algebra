#pragma once
#include<iostream>
#include<vector>
using namespace std;

//���ڵ��Եĺ���

void PrintMatrix(vector<vector<double>>& matrix);//��ӡ����

void PrintVector(std::vector<double>& vector);//��ӡ����

// ��һ�µĺ���

void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>>& U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>>& U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void gauss_elim(vector<vector<double>>& A);//Gauss��ȥ��

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//ȫ��ԪGauss��ȥ��

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

std::vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& b); //������ԪGauss��ȥ��ⷽ��

void vector_pb(vector<int>& u, vector<double>& b);//��������P*b����ѡ��

void vector_qb(vector<int>& v, vector<double>& b);//��������Q*b����ѡ��

void cholesky_decomp(vector<vector<double>>& A);//�Գ��������׼Cholesky�ֽ�

void transposeMatrix(std::vector<std::vector<double>>& inputMatrix, std::vector<std::vector<double>>& transposedMatrix);
//ʵ�־����ת�ò����ת�ú���¾���

void modified_cholesky_decomp(vector<vector<double>>& A);//�Ľ���ƽ������

void Hilbert_Matrix(vector<vector<double>>& A);

void matrix_DLT(vector<vector<double>>& A);//�������D*L^T����ѡ��

//�ڶ��µĺ���

std::vector<double> sign(std::vector<double>& w); //����w�ķ�������

double InnerProduct(std::vector<double>& a, std::vector<double>& b);//���������ڻ�

double VectorInfinityNorm(std::vector<double>& vec); //���������������

double VectorOneNorm(std::vector<double>& vec); //����������һ����

std::vector<double> UnitVectorGenerating(std::vector<double>& vec, int n); //�����������Ӧ�±�ĵ�λ����

double MatrixOneNorm(int n, std::vector<std::vector<double>>& A);//��������һ����

double MatrixInfinityNorm(std::vector<std::vector<double>>& matrix); //�������������

std::vector<double> MatrixVectorMultiply(std::vector<std::vector<double>>& A, std::vector<double>& b);//�������A������b��˵õ�������

vector<double> VectorSubtraction(vector<double>x, vector<double>y); //������������õ�������

//�����µĺ���

void house(std::vector<double>& x, std::vector<double>& v, double& beta); //����Householder�任

void QRDecomposition(std::vector<std::vector<double>>& A, std::vector<double>& d); //����QR�ֽ�

std::vector<std::vector<double>> HouseholderMatrix(std::vector<std::vector<double>>& A, std::vector<double>& d, int k); //����Householder�任����

std::vector<double> QR_equation_solving(vector<vector<double>>& A, vector<double>& b); //����QR�ֽ���ⷽ����

void equation_generating1(vector<vector<double>>& A, vector<double>& b); //���ɵ�һ�µĵ�һ��������

void equation_generating2(vector<vector<double>>& A, vector<double>& x, vector<double>& b); //���ɵ�һ�µĵڶ���������

void equation_generating3(vector<vector<double>>& A, vector<double>& b); //���ɵ�һ�µĵ�����������

std::vector<std::vector<double>> MultiplyMatrices(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2); //����������

void ReduceMatrix(std::vector<std::vector<double>>& Q, int n); //������Qֻ����ǰn��

vector<double> LS_proplem_solving(vector<vector<double>>& A, vector<double>& b); //���������С��������

double VectorTwoNorm(std::vector<double>& x); //���������Ķ�����