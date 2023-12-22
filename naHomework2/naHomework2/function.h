#pragma once
#include<iostream>
#include<vector>
using namespace std;

//���ڵ��Եĺ���

void PrintMatrix(vector<vector<double>>& mat);//��ӡ����

void PrintVector(std::vector<double>& vec);//��ӡ����

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