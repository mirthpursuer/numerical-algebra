#include "function.h"
#include <cmath>

const double EPSILON = 1e-10;  // ����һ����С�ĳ����������жϸ������Ƿ�ӽ���0

void PrintMatrix(vector<vector<double>>& mat) {
    int n = mat.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
    std::cout << std::endl;
}

void PrintVector(std::vector<double>& vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";

        // ÿʮ��������һ����
        if ((i + 1) % 10 == 0) {
            std::cout << std::endl;
        }
    }
    // ��ӡ�����з�������
    std::cout << std::endl;
}

//�㷨1.1.1�����������η����飺ǰ������
//for j = 1:n - 1
//    b(j) = b(j) / L(j, j)
//    b(j+1:n) = b(j+1:n) - b(j)*L(j + 1:n, j)
//end
//b(n) = b(n) / L(n, n)

void forward_subs(vector<vector<double>>& L, vector<double>& b) {
    int n = L.size();

    for (int j = 0; j < n - 1; j++) {
        b[j] = b[j] / L[j][j];
        for (int i = j + 1; i < n; i++) {
            b[i] -= b[j] * L[i][j];
        }
    }
    b[n - 1] = b[n - 1] / L[n - 1][n - 1];
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{
    int n = L.size();

    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * b[j];
        }
        b[i] = (b[i] - sum);  // �Խ�ԪΪ1��ǰ��������
    }
}

//�㷨1.1.2�����������η����飺�ش�����
//for j = n :-1 : 2
//    y(j) = y(j) / U(j, j)
//    (1:j - 1) = y(1:j - 1) - y(j)*U(1:j - 1, j)
//end
//y(1) = y(1) / U(1, 1)

void back_subs(vector<vector<double>>& U, vector<double>& y)
{
    int n = U.size();

    for (int j = n - 1; j >= 1; j--) {
        y[j] = y[j] / U[j][j];
        for (int i = 0; i < j; i++) {
            y[i] -= y[j] * U[i][j];
        }
    }
    y[0] = y[0] / U[0][0];
}

void back_subs1(vector<vector<double>>& U, vector<double>& y)
{
    int n = U.size();

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * y[j];
        }
        y[i] = (y[i] - sum);  // �Խ�ԪΪ1�Ļش�������
    }
}

//�㷨1.1.3�������Ƿֽ⣺Gauss��ȥ��
//for k = 1:n - 1
//    A(k + 1:n, k) = A(k + 1:n, k) / A(k, k)
//    A(k + 1:n, k + 1 : n) = A(k + 1:n, k + 1 : n) - A(k + 1:n, k)*A(k, k + 1:n)
//end

void gauss_elim(vector<vector<double>>& A)
{
    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            A[k][i] = factor; // Store the factor in the lower part of A
        }
    }

}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{
    int n = A.size();
    u.resize(n);
    v.resize(n);

    for (int i = 0; i < n; i++) {
        u[i] = i;
        v[i] = i;
    }

    for (int k = 0; k < n; k++) {
        double maxElem = 0.0;
        int i_max = k, j_max = k;
        for (int i = k; i < n; i++) {
            for (int j = k; j < n; j++) {
                if (fabs(A[i][j]) > maxElem) {
                    maxElem = fabs(A[i][j]);
                    i_max = i;
                    j_max = j;
                }
            }
        }

        if (maxElem == 0) {
            cerr << "Singular matrix encountered. LU decomposition failed." << endl;
            return;
        }

        swap(u[k], u[i_max]);
        swap(v[k], v[j_max]);

        for (int i = 0; i < n; i++) {
            swap(A[i][k], A[i][j_max]);
        }

        for (int i = k + 1; i < n; i++) {
            A[i][k] = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }
}

//�㷨1.2.2����������Ԫ���Ƿֽ⣺����ԪGauss��ȥ����
//for k = 1:n - 1
//    ȷ��p(k��p��n), ʹ��
//    | A(p, k) | = max{ |A(i,k)|��i = k:n }
//    A(k, 1:n)<->A(p, 1:n)(������k�к͵�p��)
//    u(k) = p(��¼�û�����P_k)
//    if A(k, k)��0
//        A(k + 1:n, k) = A(k + 1:n, k) / A(k, k)
//        A(k+1��n, k + 1��n) = A(k+1��n, k + 1:n) - A(k + 1:n, k)*A(k, k + 1:n)
//    else
//        stop(��������)
//    end
//end

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
{
    int n = A.size();
    double max = 0;
    int p;
    vector<double> vec(n);
    for (int k = 0; k <= n - 2; k++)
    {
        //ȷ����Ԫ
        max = fabs(A[k][k]);
        p = k;
        for (int i = k; i <= n - 1; i++)
        {
            if (fabs(A[i][k]) > max) {
                max = fabs(A[i][k]);
                p = i;
            }
        }
        u[k] = p;
        //������k��p�У���k��q��
        vec = A[k];
        A[k] = A[p];
        A[p] = vec;
    }
    for (int k = 0; k <= n - 2; k++)
    {
        if (A[k][k] != 0) {
            for (int i = k + 1; i <= n - 1; i++)
            {
                A[i][k] = A[i][k] / A[k][k];
            }
            for (int i = k + 1; i <= n - 1; i++)
            {
                for (int j = k + 1; j <= n - 1; j++)
                {
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                }
            }
        }
        else {
            cout << "��������" << endl;
            break;
        }
    }
}

//std::vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& u)
//{
//    gauss_elim_col_pivoting2(A);
//    forward_subs1(A, u);
//    back_subs(A, u);
//    return u;
//}

void vector_pb(vector<int>& u, vector<double>& b)
{
    int n = u.size();
    double key;
    for (int i = 0; i <= n - 2; i++)
    {
        key = b[i];
        b[i] = b[u[i]];
        b[u[i]] = key;
    }
}

std::vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& u)
{
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = u;
    std::vector<int>v(u.size());
    gauss_elim_col_pivoting(A2, v);
    vector_pb(v, b2);
    forward_subs1(A2, b2);
    back_subs(A2, b2);
    return b2;
}

void vector_qb(vector<int>& v, vector<double>& b)
{
    // ��������Q*b
    int n = v.size();
    vector<double> temp_b(n, 0);

    for (int i = 0; i < n; i++) {
        temp_b[i] = b[v[i]];
    }

    for (int i = 0; i < n; i++) {
        b[i] = temp_b[i];
    }
}

//�㷨1.3.1������Cholesky�ֽ⣺ƽ��������
//for k 1:n
//    A(k.k) = A(k.k)
//    A(k + 1:n, k) = A(k + 1:n��k) / A(k, k)
//    for j = k + 1 : n
//        A(j : n, j�� = A(j : n, j�� - A(j : n, k) * A(j, k)
//    end
//end

void cholesky_decomp(std::vector<std::vector<double>>& A) {
    int n = A.size();

    for (int k = 0; k < n; ++k) {
        A[k][k] = std::sqrt(A[k][k]);

        for (int i = k + 1; i < n; ++i) {
            A[i][k] /= A[k][k];
        }

        for (int j = k + 1; j < n; ++j) {
            for (int i = j; i < n; ++i) {
                A[i][j] -= A[i][k] * A[j][k];
            }
        }
    }
}

void transposeMatrix(std::vector<std::vector<double>>& inputMatrix, std::vector<std::vector<double>>& transposedMatrix) {
    int rows = inputMatrix.size();
    int cols = inputMatrix[0].size();

    // Resize the transposed matrix
    transposedMatrix.resize(cols, std::vector<double>(rows, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposedMatrix[j][i] = inputMatrix[i][j];
        }
    }
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{
    int n = A.size();

    for (int i = 0; i < n; i++) {
        double sum = A[i][i];
        for (int k = 0; k < i; k++)
            sum -= A[i][k] * A[i][k];

        if (sum < 0) {
            cerr << "������������޷�����Cholesky�ֽ⡣" << endl;
            return;
        }

        A[i][i] = sqrt(sum);

        for (int j = i + 1; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += A[j][k] * A[i][k];
            A[j][i] = (A[j][i] - sum) / A[i][i];  // �Ľ���ƽ����������
        }
    }
}

void Hilbert_Matrix(std::vector<std::vector<double>>& A) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = 1.0 / (i + j + 1);
        }
    }
}

void matrix_DLT(vector<vector<double>>& A)
{
    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = 0;
            for (int k = 0; k <= i; k++)
                sum += A[i][k] * A[j][k];
            A[i][j] = sum;  // �������D*L^T
        }
    }
}

//�㷨 2.5.1 ���ƾ����1�������Ż���
//k = 1
//while k = 1
//    w = Bx; v = sign(w); z = BTv
//    if || z || �� �� zTx
//        �� = || w || 1
//        k = 0
//    else
//        x = ej, �����±�j���� |zj| = ||z||��
//        k = 1
//    end
//end

std::vector<double> sign(std::vector<double>& w) {
    std::vector<double> v(w.size());

    for (size_t i = 0; i < w.size(); i++) {
        if (w[i] > 0) {
            v[i] = 1.0;  // ����
        }
        else if (w[i] < 0) {
            v[i] = -1.0; // ����
        }
        else {
            v[i] = 0.0;  // ��
        }
    }

    return v;
}

double InnerProduct(std::vector<double>& a, std::vector<double>& b) {
    if (a.size() != b.size()) {
        // �������Ȳ�һ�£��޷������ڻ�
        throw std::invalid_argument("Vector sizes do not match");
    }

    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }

    return result;
}

double VectorInfinityNorm(std::vector<double>& vec) {
    double norm = 0.0;
    for (double& element : vec) {
        double absElement = std::abs(element);
        if (absElement > norm) {
            norm = absElement;
        }
    }
    return norm;
}

double VectorOneNorm(std::vector<double>& vec) {
    double norm = 0.0;
    for (double& element : vec) {
        norm += std::abs(element);
    }
    return norm;
}

std::vector<double> UnitVectorGenerating(std::vector<double>& vec, int n) {
    double norm = 0.0;
    int maxIndex = -1; // ��ʼ��Ϊ-1����ʾδ�ҵ����ֵ
    for (size_t i = 0; i < vec.size(); i++) {
        double absElement = std::abs(vec[i]);
        if (absElement > norm) {
            norm = absElement;
            maxIndex = static_cast<int>(i); // �������ֵ���±�
        }
    }
    std::vector<double> unitVector(n, 0.0);  // ��ʼ��Ϊ0
    unitVector[maxIndex] = 1.0;  // ����ָ���±����Ϊ1
    return unitVector;
}

double MatrixInfinityNorm(std::vector<std::vector<double>>& matrix) {
    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();
    double maxNorm = 0.0;
    for (size_t i = 0; i < numRows; i++) {
        //PrintVector(matrix[i]);
        if (VectorOneNorm(matrix[i]) > maxNorm) {
            maxNorm = VectorOneNorm(matrix[i]);
        }
    }
    return maxNorm;
}

double MatrixOneNorm(int n, std::vector<std::vector<double>>& A) {
    int k = 1;
    int t = 0;
    //double OneNorm;
    vector<vector<double>> AT(n, vector<double>(n));
    transposeMatrix(A, AT);
    std::vector<double> w(n);
    std::vector<double> v(n);
    std::vector<double> z(n);
    std::vector<double> x(n, 1.0 / n);
    while (k == 1) {
        t++;
        // ���� w = Bx
        w = gauss_equation_solving(AT, x);
        // ���� v = sign(w)
        v = sign(w);
        // ���� z = BTv
        z = gauss_equation_solving(A, v);
        /*PrintVector(w);
        PrintVector(v);
        PrintVector(z);*/
        // ���� ||z||�� �� zTx
        if (VectorInfinityNorm(z) <= InnerProduct(z, x)) {
            /*OneNorm = VectorOneNorm(w);
            std::cout << OneNorm << std::endl;*/
            k = 0;
        }
        else {
            // ���� x Ϊ ej
            x = UnitVectorGenerating(z, n);
            //PrintVector(x);
            k = 1;
        }
        if (t > 100) { break; }
    }
    return VectorOneNorm(w);
}

// ����һ������������� A ������ b �ĳ˻�
std::vector<double> MatrixVectorMultiply(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int numRows = A.size();
    int numCols = A[0].size();
    std::vector<double> result(numRows, 0.0); // ��ʼ���������Ϊ0
    if (numCols != b.size()) {
        // �����������������������Ĵ�С���޷����
        return result; // ����ȫ������
    }
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i] += A[i][j] * b[j];
        }
    }
    return result;
}

vector<double> VectorSubtraction(vector<double>x, vector<double>y)
{
    int n = x.size();
    for (int i = 0; i < n; i++)
    {
        x[i] -= y[i];
    }
    return x;
}

