#include "function.h"
#include <cmath>
#include <random>
#include <iomanip> 

const double EPSILON = 1e-10;  // ����һ����С�ĳ����������жϸ������Ƿ�ӽ���0

void PrintMatrix(vector<vector<double>>& mat) {
    // �����������Ϊ��λС��
    //std::cout << std::fixed << std::setprecision(4);
    int m = mat.size();
    int n = mat[0].size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
    std::cout << std::endl;
}

void PrintVector(std::vector<double>& vec) {
    // �����������Ϊ��λС��
    std::cout << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
        // ÿʮ��������һ����
        if ((i + 1) % 11 == 0) {
            std::cout << std::endl;
        }
    }
    // ��ӡ�����з�������
    std::cout << std::endl;
    //std::cout << std::endl;
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
    vector<double>a(n);
    for (int i = 0; i < n; i++)
    {
        a[i] = x[i] - y[i];
    }
    return a;
}

//�㷨3.2.1������Householder�任��
//function:[v,beta]=house(x)
//	n=length(x)(����x�ĳ���)
//	eta = ||x||�ޣ�x=x/eta
//	sigma = x(2:n)T ��  x(2:n)
//	v(2:n)=x(2:n)
//	if sigma =0
//		beta=0
//	else
//		alpha=(x(1)^2+sigma)^0.5
//		if x(1)��0
//			v(1)=x(1)-alpha
//		else
//			v(1)=-sigma/((x(1)+alpha)
//		end
//		beta=2v(1)^2/(sigma+v(1)^2);v=v/v(1)
//	end

void house(std::vector<double>& x, std::vector<double>& v, double& beta) {
    int n = x.size();
    double eta = 0.0;
    for (int i = 0; i < n; ++i) {
        eta = std::max(eta, std::abs(x[i]));
    }
    /*std::cout << "eta:" << std::endl;
    cout << eta << endl;*/
    for (int i = 0; i < n; ++i) {
        x[i] /= eta;
    }
    /*std::cout << "x:" << std::endl;
    PrintVector(x);*/
    double sigma = 0.0;
    for (int i = 1; i < n; ++i) {
        sigma += x[i] * x[i];
    }
    /*std::cout << "sigma:" << std::endl;
    cout << sigma << endl;*/
    v.resize(n, 0.0);
    for (int i = 1; i < n; ++i) {
        v[i] = x[i];
    }
    /*std::cout << "v:" << std::endl;
    PrintVector(v);*/
    if (abs(sigma) < EPSILON) {
        beta = 0.0;
        /*std::cout << "beta:" << std::endl;
        cout << beta << endl;*/
    }
    else {
        double alpha = std::sqrt(x[0] * x[0] + sigma);
        /*std::cout << "alpha:" << std::endl;
        cout << alpha << endl;*/
        if (x[0] <= 0) {
            v[0] = x[0] - alpha;
            /*std::cout << "v:" << std::endl;
            PrintVector(v);*/
        }
        else {
            v[0] = -sigma / (x[0] + alpha);
            /*std::cout << "v:" << std::endl;
            PrintVector(v);*/
        }
        beta = (2.0 * v[0] * v[0]) / (sigma + v[0] * v[0]);
        /*std::cout << "beta:" << std::endl;
        cout << beta << endl;*/
        for (int i = 1; i < n; ++i) {
            v[i] /= v[0];
        }
        v[0] = 1;
        /*std::cout << "v:" << std::endl;
        PrintVector(v);*/
    }
}

//�㷨3.3.1������QR�ֽ⣺Householder������
//for j=1:n
//	if j<m
//		[v,beta] = house(A(j:m,j))
//		A(j:m,j:n) = (I_{m-j+1} - beta*v*vT) * A(j:m,j:n) //���Ǿ���˷�
//		d(j) = beta
//		A(j+1:m,j) = v(2:m-j+1)
//	end
//end

void QRDecomposition(std::vector<std::vector<double>>& A, std::vector<double>& d) {
    int m = A.size();
    int n = A[0].size();
    for (int j = 0; j < n; ++j) {
        std::vector<double> v;
        double beta;
        if (j < m - 1) {
            //if (1) {
            std::vector<double> subvector;
            for (int i = j; i < m; ++i) {
                subvector.push_back(A[i][j]);
            }
            /*std::cout << "subvector:" << std::endl;
            PrintVector(subvector);*/
            house(subvector, v, beta);
            /*std::cout << "v:" << std::endl;
            PrintVector(v);
            std::cout << "beta:" << std::endl;
            cout << beta << endl;
            std::cout << "A:" << std::endl;
            PrintMatrix(A);*/
            // Define the submatrix A(j:m, j:n)
            std::vector<std::vector<double>> submatrix(m - j, std::vector<double>(n - j, 0.0));
            for (int i = j; i < m; ++i) {
                for (int k = j; k < n; ++k) {
                    submatrix[i - j][k - j] = A[i][k];
                }
            }
            /*std::cout << "submatrix:" << std::endl;
            PrintMatrix(submatrix);*/
            // Define the matrix (I - beta * v * v^T)
            std::vector<std::vector<double>> identity(m - j, std::vector<double>(m - j, 0.0));
            for (int i = 0; i < m - j; ++i) {
                for (int k = 0; k < m - j; ++k) {
                    identity[i][k] = (i == k) ? 1.0 - beta * v[i] * v[k] : -beta * v[i] * v[k];
                }
            }
            /*std::cout << "identity:" << std::endl;
            PrintMatrix(identity);*/
            // Multiply the submatrix by the identity matrix
            std::vector<std::vector<double>> result(m - j, std::vector<double>(n - j, 0.0));
            for (int i = 0; i < m - j; ++i) {
                for (int k = 0; k < n - j; ++k) {
                    for (int p = 0; p < m - j; ++p) {
                        result[i][k] += identity[i][p] * submatrix[p][k];
                    }
                }
            }
            /*std::cout << "result:" << std::endl;
            PrintMatrix(result);*/
            // Update the original matrix A with the result
            for (int i = j; i < m; ++i) {
                for (int k = j; k < n; ++k) {
                    A[i][k] = result[i - j][k - j];
                }
            }
            /*std::cout << "new A:" << std::endl;
            PrintMatrix(A);*/
            d[j] = beta;
            for (int i = j + 1; i < m; ++i) {
                A[i][j] = v[i - j];
            }
            /*std::cout << "newnew A:" << std::endl;
            PrintMatrix(A);*/
        }
    }
}

//�����任���Ļ�������Ϊ��
//(1)����A��QR�ֽ⣻
//(2)����c1 = QTb��
//(3)��������Ƿ�����Rx = c1��
//����A=QR=H1������HnR, ��ôc1= Hn������H1 b, 
//����Hk = diag(I_{k-1}, I_{m-k+1} - beta*v_k��v_k^T)Ϊm*mά����
//�ݹ����MatrixVectorMultiply�������ɼ���õ�c1,����ûش���back_subs����õ�x��

std::vector<std::vector<double>> HouseholderMatrix(std::vector<std::vector<double>>& A, std::vector<double>& d, int k) {
    int m = A.size();
    int n = A[0].size();
    std::vector<std::vector<double>> matrix(m, std::vector<double>(m, 0.0));
    for (int i = 0; i < k - 1; ++i) {
        matrix[i][i] = 1.0;
    }
    //PrintMatrix(matrix);
    std::vector<double> v = { 1.0 };
    // ���� k �еĺ� m-k ���������ӵ����� v ����
    for (int i = k; i < m; ++i) {
        v.push_back(A[i][k - 1]);
    }
    //PrintVector(v);
    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v.size(); ++j) {
            matrix[k + i - 1][k + j - 1] = -d[k - 1] * v[i] * v[j];
        }
    }
    //PrintMatrix(matrix);
    for (int i = k - 1; i < A.size(); ++i) {
        matrix[i][i] += 1;
    }
    return matrix;
}

std::vector<double> QR_equation_solving(vector<vector<double>>& A, vector<double>& b) {
    int m = A.size();
    int n = A[0].size();
    // QR �ֽ�
    std::vector<double> d(n, 0.0);
    QRDecomposition(A, d);
    /*PrintMatrix(A);
    PrintVector(d);*/
    // ���� c1 = Q^T * b
    std::vector<double> c1 = b;
    //PrintVector(c1);
    std::vector<std::vector<double>> H(m, std::vector<double>(m, 0.0));
    for (int k = 1; k <= n; k++) {
        H = HouseholderMatrix(A, d, k);
        c1 = MatrixVectorMultiply(H, c1);
        /*cout << "��" << k << "��Householder�任��" << endl;
        PrintMatrix(H);
        PrintVector(c1);*/
    }
    // �������Ƿ����� Rx = c1
    back_subs(A, c1);
    return c1;
}

void equation_generating1(vector<vector<double>>& A, vector<double>& b) {
    int N = A.size();
    for (int i = 0; i < N - 1; i++)
    {
        A[i][i] = 6;
        A[i + 1][i] = 8;
        A[i][i + 1] = 1;
        b[i] = 15;
    }
    A[N - 1][N - 1] = 6;
    b[0] = 7;
    b[N - 1] = 14;
    //cout << "����AΪ��" << endl;
    //PrintMatrix(A);
    //cout << "����bΪ��" << endl;
    //PrintVector(b);
}

void equation_generating2(vector<vector<double>>& A, vector<double>& x, vector<double>& b) {
    int N = A.size();
    for (int i = 0; i < N - 1; i++)
    {
        A[i][i] = 10;
        A[i + 1][i] = 1;
        A[i][i + 1] = 1;
    }
    A[N - 1][N - 1] = 1;
    // ��ʼ�������������
    std::random_device rd;  // �������
    std::mt19937 gen(rd());  // Mersenne Twister ������
    std::uniform_real_distribution<double> dis(0.0, 1.0);  // ���� [0.0, 1.0) ��Χ�ڵ���� double ��
    // ������� double �����������
    for (int i = 0; i < N; ++i) {
        x[i] = dis(gen);
    }
    b = MatrixVectorMultiply(A, x);
    /*cout << "����AΪ��" << endl;
    PrintMatrix(A);*/
    cout << "����xΪ��" << endl;
    PrintVector(x);
    /*cout << "����bΪ��" << endl;
    PrintVector(b);*/
}

void equation_generating3(vector<vector<double>>& A, vector<double>& b) {
    int N = A.size();
    Hilbert_Matrix(A);
    //PrintMatrix(A);
    for (int i = 0; i < N; ++i) {
        for (int j = 1; j <= N; ++j) {
            b[i] += 1.0 / (i + j);
        }
    }
    /*cout << "����AΪ��" << endl;
    PrintMatrix(A);
    cout << "����bΪ��" << endl;
    PrintVector(b);*/
}

std::vector<std::vector<double>> MultiplyMatrices(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();
    // ������������Ƿ�������
    if (cols1 != rows2) {
        std::cerr << "Error: Matrix dimensions are not compatible for multiplication." << std::endl;
        return std::vector<std::vector<double>>();
    }
    // ��ʼ���������
    std::vector<std::vector<double>> result(rows1, std::vector<double>(cols2, 0.0));
    // ���о������
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            for (int k = 0; k < cols1; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

void ReduceMatrix(std::vector<std::vector<double>>& Q, int n) {
    for (std::vector<double>& row : Q) {
        if (row.size() > n) {
            row.erase(row.begin() + n, row.end()); // ɾ���ӵ� n+1 �е����һ�е�Ԫ��
        }
    }
}

vector<double> LS_proplem_solving(vector<vector<double>>& A, vector<double>& b) {
    int m = A.size();
    int n = A[0].size();
    //PrintMatrix(A);
    std::vector<double> d(n, 0.0);
    QRDecomposition(A, d);
    //PrintMatrix(A);
    std::vector<std::vector<double>> Q(m, std::vector<double>(m, 0.0));
    std::vector<std::vector<double>> H(m, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; ++i) {
        Q[i][i] = 1.0;
    }
    //PrintMatrix(Q);
    for (int k = 1; k <= n; k++) {
        H = HouseholderMatrix(A, d, k);
        Q = MultiplyMatrices(Q, H);
        /*PrintMatrix(H);
        PrintMatrix(Q);
        cout << "\n" << endl;*/
    }
    ReduceMatrix(Q, n);
    //PrintMatrix(Q);
    std::vector<std::vector<double>> Q1T(n, std::vector<double>(m, 0.0));
    transposeMatrix(Q, Q1T);
    //PrintMatrix(Q1T);
    std::vector<double> c1(n, 0.0);
    c1 = MatrixVectorMultiply(Q1T, b);
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            R[i][j] = A[i][j];
        }
    }
    //PrintMatrix(R);
    back_subs(R, c1);
    return c1;
}

double VectorTwoNorm(std::vector<double>& x) {
    double norm = 0.0;
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        norm += x[i] * x[i];
    }
    return std::sqrt(norm);
}

vector<double> Jacobi_Iteration(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0.0); // ��ʼ��������x
    vector<double> x_new(n, 0.0); // ���ڴ洢�µĵ������
    vector<double> diff(n, 0.0);  // ���ڴ洢�����ͱ�׼������
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int k = 0; k < maxIterations; k++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }
        diff = VectorSubtraction(x_new, x);
        if (VectorInfinityNorm(diff) < tolerance) {
            cout << "�� " << k + 1 << " �ε���������." << endl;
            return x_new;
        }
        x = x_new;
    }
    cout << "��ָ���ĵ���������δ����." << endl;
    return x;
}

vector<double> GS_Iteration(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0.0); // ��ʼ��������x
    vector<double> x_new(n, 0.0); // ���ڴ洢�µĵ������
    vector<double> diff(n, 0.0);  // ���ڴ洢�����ͱ�׼������
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; // ���ڼ���ǰ�벿�ֵ����
            double sum2 = 0.0; // ���ڼ����벿�ֵ����
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x_new[j];
            }
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x[j];
            }
            x_new[i] = (1.0 / A[i][i]) * (b[i] - sum1 - sum2);
        }
        diff = VectorSubtraction(x_new, x);
        if (VectorInfinityNorm(diff) < tolerance) {
            cout << "�� " << iter + 1 << " �ε���������." << endl;
            return x_new;
        }
        x = x_new;
    }
    cout << "��ָ���ĵ���������δ����." << endl;
    return x;
}

vector<double> SOR_Iteration(vector<vector<double>>& A, vector<double>& b, double omega) {
    int n = A.size();
    vector<double> x(n, 0.0); // ��ʼ��������x
    vector<double> x_new(n, 0.0); // ���ڴ洢�µĵ������
    vector<double> diff(n, 0.0);  // ���ڴ洢�����ͱ�׼������
    //cout << "omega=" << omega << endl;
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; // ���ڼ���ǰ�벿�ֵ����
            double sum2 = 0.0; // ���ڼ����벿�ֵ����
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x_new[j];
            }
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x[j];
            }
            x_new[i] = (1.0 - omega) * x[i] + omega * ((1.0 / A[i][i]) * (b[i] - sum1 - sum2));
        }
        diff = VectorSubtraction(x_new, x);
        if (VectorInfinityNorm(diff) < tolerance) {
            cout << "�� " << iter + 1 << " �ε���������." << endl;
            //cout << iter + 1 ;
            return x_new;
        }
        x = x_new;
    }
    cout << "��ָ���ĵ���������δ����." << endl;
    return x;
}

int SOR_Performance(vector<vector<double>>& A, vector<double>& b, double omega) {
    int n = A.size();
    vector<double> x(n, 0.0); // ��ʼ��������x
    vector<double> x_new(n, 0.0); // ���ڴ洢�µĵ������
    vector<double> diff(n, 0.0);  // ���ڴ洢�����ͱ�׼������
    //cout << "omega=" << omega << endl;
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; // ���ڼ���ǰ�벿�ֵ����
            double sum2 = 0.0; // ���ڼ����벿�ֵ����
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x_new[j];
            }
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x[j];
            }
            x_new[i] = (1.0 - omega) * x[i] + omega * ((1.0 / A[i][i]) * (b[i] - sum1 - sum2));
        }
        diff = VectorSubtraction(x_new, x);
        if (VectorInfinityNorm(diff) < tolerance) {
            return iter + 1;
        }
        x = x_new;
    }
    return maxIterations;
}

double BisearchOmega(vector<vector<double>>& A, vector<double>& b) {
    double lowerBound = 1.0; // �����½�
    double upperBound = 2.0; // �����Ͻ�
    const int maxIterations = 20000; // ����������
    const double tolerance = 1e-2; // ������ֵ
    int bestIterations = maxIterations;
    double bestOmega = lowerBound;

    clock_t start_time;
    start_time = clock();

    double lower = SOR_Performance(A, b, lowerBound);
    double upper = SOR_Performance(A, b, upperBound);
    bestIterations = min(lower, upper);
    if (lower < upper) {
        bestOmega = lowerBound;
    }
    else {
        bestOmega = upperBound;
    }

    cout << "Ѱ�����omega��..." << endl;
    //cout << bestIterations << endl;
    
    while (upperBound - lowerBound > tolerance) {
        double mid = (lowerBound + upperBound) / 2;
        int performance = SOR_Performance(A, b, mid);
        //cout << performance << endl;

        if (performance <= bestIterations) {
            bestIterations = performance;
            bestOmega = mid;
            if (performance < SOR_Performance(A, b, mid + tolerance)) {
                upperBound = mid;
                upper = performance;
            }
            else {
                lowerBound = mid;
                lower = performance;
            }
        }
        else {
            if (lower < upper) {
                upperBound = mid;
            }
            else {
                lowerBound = mid;
            }
        }
    }
    cout << "���omegaΪ��" << bestOmega  << "\nѰ����ʱ:" << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
    return bestOmega;
}

void Iterations(double epsilon) {
    int n = 99; //��[0,1]����100�ȷ�
    double a = 0.5;
    double h = 1.0 / (n + 1);
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    //��ʼ��A��b
    for (int i = 0; i < n - 1; i++)
    {
        A[i][i] = -2 * epsilon - h;
        A[i + 1][i] = epsilon;
        A[i][i + 1] = epsilon + h;
        b[i] = a * h * h;
    }
    A[n - 1][n - 1] = -2 * epsilon - h;
    b[n - 1] = a * h * h - epsilon - h;

    clock_t start_time;  // ��¼��ʼʱ��

    cout << "epsilon = " << epsilon << endl;
    vector<double> accurate_y(n, 0.0);
    double temp = (1 - a) / (1 - exp(-1.0 / epsilon));
    for (int i = 0; i < n; i++)
    {
        accurate_y[i] = temp * (1 - exp(-(i + 1) * h / epsilon)) + a * (i + 1) * h;
    }
    cout << "��ȷ��Ϊ��" << endl;
    PrintVector(accurate_y);
    cout << endl;

    start_time = clock();
    cout << "Jacobi������" << endl;
    vector<double> y1 = Jacobi_Iteration(A, b);
    PrintVector(y1);
    cout << "Jacobi����������ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;

    start_time = clock();
    cout << "G-S������" << endl;
    vector<double> y2 = GS_Iteration(A, b);
    PrintVector(y2);
    cout << "G-S����������ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;

    cout << "SOR������" << endl;
    double omega = BisearchOmega(A, b);
    start_time = clock();
    vector<double> y3 = SOR_Iteration(A, b, omega);
    PrintVector(y3);
    cout << "SOR����������ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
}

vector<vector<double>> MatrixSubtraction(vector<vector<double>> x, vector<vector<double>> y){
    int n = x.size();
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            A[i][j] = x[i][j] - y[i][j];
        }
    }
    return A;
}

void Jacobi_Iteration2(vector<vector<double>>& u) {

    clock_t start_time;  // ��¼��ʼʱ��

    int n = u.size() - 1;
    double h = 1.0 / n;
    const int maxIterations = 20000;
    const double tolerance = 1e-7;
    vector<vector<double>> u1 = u;
    vector<vector<double>> diff = MatrixSubtraction(u, u1);

    start_time = clock();
    for (int k = 0; k < maxIterations; k++) {
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                u1[i][j] = u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1] + h * h * (i + j) * h;
                u1[i][j] /= 4 + h * h * exp(i * j * h * h);
            }
        }
        diff = MatrixSubtraction(u, u1);
        if (MatrixInfinityNorm(diff) < tolerance) {
            cout << "Jacobi������ " << k + 1 << " �ε���������." << endl;
            double mini = 1.0;
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (u1[i][j] < mini) {
                        mini = u1[i][j];
                    }
                }
            }
            //PrintMatrix(u1);
            cout << "�����С����Ϊ " << mini << "." << endl;
            cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
            break;
        }
        u = u1;
        if (k == maxIterations - 1) {
            cout << "��ָ���ĵ���������δ����." << endl;
        }
    }
}

void GS_Iteration2(vector<vector<double>>& u) {

    clock_t start_time;  // ��¼��ʼʱ��

    int n = u.size() - 1;
    double h = 1.0 / n;
    const int maxIterations = 20000;
    const double tolerance = 1e-7;
    vector<vector<double>> u1 = u;
    vector<vector<double>> diff = MatrixSubtraction(u, u1);

    start_time = clock();
    for (int k = 0; k < maxIterations; k++) {
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                u1[i][j] = u1[i - 1][j] + u1[i][j - 1] + u[i + 1][j] + u[i][j + 1] + h * h * (i + j) * h;
                u1[i][j] /= 4 + h * h * exp(i * j * h * h);
            }
        }
        diff = MatrixSubtraction(u, u1);
        if (MatrixInfinityNorm(diff) < tolerance) {
            cout << "G-S������ " << k + 1 << " �ε���������." << endl;
            double mini = 1.0;
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (u1[i][j] < mini) {
                        mini = u1[i][j];
                    }
                }
            }
            //PrintMatrix(u1);
            cout << "�����С����Ϊ " << mini << "." << endl;
            cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
            break;
        }
        u = u1;
        if (k == maxIterations - 1) {
            cout << "��ָ���ĵ���������δ����." << endl;
        }
    }
}

void SOR_Iteration2(vector<vector<double>>& u, double omega) {
    clock_t start_time;  // ��¼��ʼʱ��

    int n = u.size() - 1;
    double h = 1.0 / n;
    const int maxIterations = 20000;
    const double tolerance = 1e-7;
    vector<vector<double>> u1 = u;
    vector<vector<double>> diff = MatrixSubtraction(u, u1);

    start_time = clock();
    for (int k = 0; k < maxIterations; k++) {
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                double temp = 4 + h * h * exp(i * j * h * h);
                u1[i][j] = omega * (u1[i - 1][j] + u1[i][j - 1] + u[i + 1][j] + u[i][j + 1] + h * h * (i + j) * h) + (1 - omega) * temp * u[i][j];
                u1[i][j] /= temp;
            }
        }
        diff = MatrixSubtraction(u, u1);
        if (MatrixInfinityNorm(diff) < tolerance) {
            cout << "SOR������ " << k + 1 << " �ε���������." << endl;
            double mini = 1.0;
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (u1[i][j] < mini) {
                        mini = u1[i][j];
                    }
                }
            }
            //PrintMatrix(u1);
            cout << "�����С����Ϊ " << mini << "." << endl;
            cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
            break;
        }
        u = u1;
        if (k == maxIterations - 1) {
            cout << "��ָ���ĵ���������δ����." << endl;
        }
    }

}

int SOR_Performance2(vector<vector<double>>& u, double omega) {
    int n = u.size() - 1;
    double h = 1.0 / n;
    const int maxIterations = 20000;
    const double tolerance = 1e-7;
    vector<vector<double>> u1 = u;
    vector<vector<double>> u2 = u;
    //vector<vector<double>> diff = MatrixSubtraction(u, u1);
    vector<vector<double>> diff = MatrixSubtraction(u1, u2);

    for (int k = 0; k < maxIterations; k++) {
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                double temp = 4 + h * h * exp(i * j * h * h);
                /*u1[i][j] = omega * (u1[i - 1][j] + u1[i][j - 1] + u[i + 1][j] + u[i][j + 1] + h * h * (i + j) * h) + (1 - omega) * temp * u[i][j];*/
                u1[i][j] = omega * (u1[i - 1][j] + u1[i][j - 1] + u2[i + 1][j] + u2[i][j + 1] + h * h * (i + j) * h) + (1 - omega) * temp * u2[i][j];
                u1[i][j] /= temp;
            }
        }
        /*diff = MatrixSubtraction(u, u1);*/
        diff = MatrixSubtraction(u1, u2);
        if (MatrixInfinityNorm(diff) < tolerance) {
            return k + 1;
        }
        u2 = u1;
        if (k == maxIterations - 1) {
             return k;
        }
    }
}

double BisearchOmega2(vector<vector<double>>& A) {
    double lowerBound = 1.0; // �����½�
    double upperBound = 2.0; // �����Ͻ�
    const int maxIterations = 20000; // ����������
    const double tolerance = 1e-2; // ������ֵ
    int bestIterations = maxIterations;
    double bestOmega = lowerBound;

    clock_t start_time;
    start_time = clock();

    double lower = SOR_Performance2(A, lowerBound);
    double upper = SOR_Performance2(A, upperBound);
    bestIterations = min(lower, upper);
    if (lower < upper) {
        bestOmega = lowerBound;
    }
    else {
        bestOmega = upperBound;
    }

    cout << "Ѱ�����omega��..." << endl;
    //cout << bestIterations << endl;

    while (upperBound - lowerBound > tolerance) {
        double mid = (lowerBound + upperBound) / 2;
        int performance = SOR_Performance2(A, mid);
        //cout << mid << endl << performance << endl;

        if (performance <= bestIterations) {
            bestIterations = performance;
            bestOmega = mid;
            if (performance < SOR_Performance2(A, mid + tolerance)) {
                upperBound = mid;
                upper = performance;
            }
            else {
                lowerBound = mid;
                lower = performance;
            }
        }
        else {
            if (lower < upper) {
                upperBound = mid;
            }
            else {
                lowerBound = mid;
            }
        }
    }
    cout << "���omegaΪ��" << bestOmega << "\nѰ����ʱ:" << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
    return bestOmega;
}

void Iterations2(int n) {
    
    vector<vector<double>> u(n + 1, vector<double>(n + 1));
    //��ʼ��A��b
    for (int i = 0; i <= n; i++)
    {
        u[0][i] = 1;
        u[i][0] = 1;
        u[n][i] = 1;
        u[i][n] = 1;
    }
    vector<vector<double>> backup = u;
    cout << "N = " << n << endl;
    
    Jacobi_Iteration2(u);
    u = backup;
    GS_Iteration2(u);
    u = backup;
    double omega = BisearchOmega2(u);
    SOR_Iteration2(u, omega);
}