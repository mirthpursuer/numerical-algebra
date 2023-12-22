#include "function.h"
#include <cmath>
#include <random>
#include <iomanip> 

const double EPSILON = 1e-10;  // ����һ����С�ĳ����������жϸ������Ƿ�ӽ���0

void PrintMatrix(vector<vector<double>>& mat) {
    // �����������Ϊ��λС��
    cout << fixed << setprecision(4);
    int m = mat.size();
    int n = mat[0].size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void PrintVector(vector<double>& vec) {
    // �����������Ϊ��λС��
    cout << fixed << setprecision(4);
    for (size_t i = 0; i < vec.size(); i++) {
        cout << vec[i] << " ";
        // ÿʮ��������һ����
        if ((i + 1) % 10 == 0) {
            cout << endl;
        }
    }
    // ��ӡ�����з�������
    cout << endl;
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

//vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& u)
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

vector<double> gauss_equation_solving(vector<vector<double>>& A, vector<double>& u)
{
    vector<vector<double>> A2 = A;
    vector<double> b2 = u;
    vector<int>v(u.size());
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

void cholesky_decomp(vector<vector<double>>& A) {
    int n = A.size();

    for (int k = 0; k < n; ++k) {
        A[k][k] = sqrt(A[k][k]);

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

void transposeMatrix(vector<vector<double>>& inputMatrix, vector<vector<double>>& transposedMatrix) {
    int rows = inputMatrix.size();
    int cols = inputMatrix[0].size();

    // Resize the transposed matrix
    transposedMatrix.resize(cols, vector<double>(rows, 0.0));

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

void Hilbert_Matrix(vector<vector<double>>& A) {
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

vector<double> sign(vector<double>& w) {
    vector<double> v(w.size());

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

int sign(double m) {
    if (m > 0)return 1;
    else if (m < 0)return -1;
    else return 0;
}

double InnerProduct(vector<double>& a, vector<double>& b) {
    if (a.size() != b.size()) {
        // �������Ȳ�һ�£��޷������ڻ�
        throw invalid_argument("Vector sizes do not match");
    }

    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }

    return result;
}

double VectorInfinityNorm(vector<double>& vec) {
    double norm = 0.0;
    for (double& element : vec) {
        double absElement = abs(element);
        if (absElement > norm) {
            norm = absElement;
            //norm = element;
        }
    }
    return norm;
}

double VectorOneNorm(vector<double>& vec) {
    double norm = 0.0;
    for (double& element : vec) {
        norm += abs(element);
    }
    return norm;
}

vector<double> UnitVectorGenerating(vector<double>& vec, int n) {
    double norm = 0.0;
    int maxIndex = -1; // ��ʼ��Ϊ-1����ʾδ�ҵ����ֵ
    for (size_t i = 0; i < vec.size(); i++) {
        double absElement = abs(vec[i]);
        if (absElement > norm) {
            norm = absElement;
            maxIndex = static_cast<int>(i); // �������ֵ���±�
        }
    }
    vector<double> unitVector(n, 0.0);  // ��ʼ��Ϊ0
    unitVector[maxIndex] = 1.0;  // ����ָ���±����Ϊ1
    return unitVector;
}

double MatrixInfinityNorm(vector<vector<double>>& matrix) {
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

double MatrixOneNorm(int n, vector<vector<double>>& A) {
    int k = 1;
    int t = 0;
    //double OneNorm;
    vector<vector<double>> AT(n, vector<double>(n));
    transposeMatrix(A, AT);
    vector<double> w(n);
    vector<double> v(n);
    vector<double> z(n);
    vector<double> x(n, 1.0 / n);
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
            cout << OneNorm << endl;*/
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

vector<double> MatrixVectorMultiply(vector<vector<double>>& A, vector<double>& b) {
    int numRows = A.size();
    int numCols = A[0].size();
    vector<double> result(numRows, 0.0); // ��ʼ���������Ϊ0
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

void house(vector<double>& x, vector<double>& v, double& beta) {
    int n = x.size();
    double eta = 0.0;
    for (int i = 0; i < n; ++i) {
        eta = max(eta, abs(x[i]));
    }
    /*cout << "eta:" << endl;
    cout << eta << endl;*/
    for (int i = 0; i < n; ++i) {
        x[i] /= eta;
    }
    /*cout << "x:" << endl;
    PrintVector(x);*/
    double sigma = 0.0;
    for (int i = 1; i < n; ++i) {
        sigma += x[i] * x[i];
    }
    /*cout << "sigma:" << endl;
    cout << sigma << endl;*/
    v.resize(n, 0.0);
    for (int i = 1; i < n; ++i) {
        v[i] = x[i];
    }
    /*cout << "v:" << endl;
    PrintVector(v);*/
    if (abs(sigma) < EPSILON) {
        beta = 0.0;
        /*cout << "beta:" << endl;
        cout << beta << endl;*/
    }
    else {
        double alpha = sqrt(x[0] * x[0] + sigma);
        /*cout << "alpha:" << endl;
        cout << alpha << endl;*/
        if (x[0] <= 0) {
            v[0] = x[0] - alpha;
            /*cout << "v:" << endl;
            PrintVector(v);*/
        }
        else {
            v[0] = -sigma / (x[0] + alpha);
            /*cout << "v:" << endl;
            PrintVector(v);*/
        }
        beta = (2.0 * v[0] * v[0]) / (sigma + v[0] * v[0]);
        /*cout << "beta:" << endl;
        cout << beta << endl;*/
        for (int i = 1; i < n; ++i) {
            v[i] /= v[0];
        }
        v[0] = 1;
        /*cout << "v:" << endl;
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

void QRDecomposition(vector<vector<double>>& A, vector<double>& d) {
    int m = A.size();
    int n = A[0].size();
    for (int j = 0; j < n; ++j) {
        vector<double> v;
        double beta;
        if (j < m - 1) {
            vector<double> subvector;
            for (int i = j; i < m; ++i) {
                subvector.push_back(A[i][j]);
            }
            /*cout << "subvector:" << endl;
            PrintVector(subvector);*/
            house(subvector, v, beta);
            /*cout << "v:" << endl;
            PrintVector(v);
            cout << "beta:" << endl;
            cout << beta << endl;
            cout << "A:" << endl;
            PrintMatrix(A);*/
            // Define the submatrix A(j:m, j:n)
            vector<vector<double>> submatrix(m - j, vector<double>(n - j, 0.0));
            for (int i = j; i < m; ++i) {
                for (int k = j; k < n; ++k) {
                    submatrix[i - j][k - j] = A[i][k];
                }
            }
            // Define the matrix (I - beta * v * v^T)
            vector<vector<double>> identity(m - j, vector<double>(m - j, 0.0));
            for (int i = 0; i < m - j; ++i) {
                for (int k = 0; k < m - j; ++k) {
                    identity[i][k] = (i == k) ? 1.0 - beta * v[i] * v[k] : -beta * v[i] * v[k];
                }
            }
            /*cout << "identity:" << endl;
            PrintMatrix(identity);*/
            // Multiply the submatrix by the identity matrix
            vector<vector<double>> result(m - j, vector<double>(n - j, 0.0));
            for (int i = 0; i < m - j; ++i) {
                for (int k = 0; k < n - j; ++k) {
                    for (int p = 0; p < m - j; ++p) {
                        result[i][k] += identity[i][p] * submatrix[p][k];
                    }
                }
            }
            /*cout << "result:" << endl;
            PrintMatrix(result);*/
            // Update the original matrix A with the result
            for (int i = j; i < m; ++i) {
                for (int k = j; k < n; ++k) {
                    A[i][k] = result[i - j][k - j];
                }
            }
            /*cout << "new A:" << endl;
            PrintMatrix(A);*/
            d[j] = beta;
            for (int i = j + 1; i < m; ++i) {
                A[i][j] = v[i - j];
            }
            /*cout << "newnew A:" << endl;
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

vector<vector<double>> HouseholderMatrix(vector<vector<double>>& A, vector<double>& d, int k) {
    int m = A.size();
    int n = A[0].size();
    vector<vector<double>> matrix(m, vector<double>(m, 0.0));
    for (int i = 0; i < k - 1; ++i) {
        matrix[i][i] = 1.0;
    }
    //PrintMatrix(matrix);
    vector<double> v = { 1.0 };
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

vector<double> QR_equation_solving(vector<vector<double>>& A, vector<double>& b) {
    int m = A.size();
    int n = A[0].size();
    // QR �ֽ�
    vector<double> d(n, 0.0);
    QRDecomposition(A, d);
    /*PrintMatrix(A);
    PrintVector(d);*/
    // ���� c1 = Q^T * b
    vector<double> c1 = b;
    //PrintVector(c1);
    vector<vector<double>> H(m, vector<double>(m, 0.0));
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
    random_device rd;  // �������
    mt19937 gen(rd());  // Mersenne Twister ������
    uniform_real_distribution<double> dis(0.0, 1.0);  // ���� [0.0, 1.0) ��Χ�ڵ���� double ��
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

vector<vector<double>> MultiplyMatrices(vector<vector<double>>& matrix1, vector<vector<double>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();
    // ������������Ƿ�������
    if (cols1 != rows2) {
        cerr << "Error: Matrix dimensions are not compatible for multiplication." << endl;
        return vector<vector<double>>();
    }
    // ��ʼ���������
    vector<vector<double>> result(rows1, vector<double>(cols2, 0.0));
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

void ReduceMatrix(vector<vector<double>>& Q, int n) {
    for (vector<double>& row : Q) {
        if (row.size() > n) {
            row.erase(row.begin() + n, row.end()); // ɾ���ӵ� n+1 �е����һ�е�Ԫ��
        }
    }
}

vector<double> LS_proplem_solving(vector<vector<double>>& A, vector<double>& b) {
    int m = A.size();
    int n = A[0].size();
    //PrintMatrix(A);
    vector<double> d(n, 0.0);
    QRDecomposition(A, d);
    //PrintMatrix(A);
    vector<vector<double>> Q(m, vector<double>(m, 0.0));
    vector<vector<double>> H(m, vector<double>(m, 0.0));
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
    vector<vector<double>> Q1T(n, vector<double>(m, 0.0));
    transposeMatrix(Q, Q1T);
    //PrintMatrix(Q1T);
    vector<double> c1(n, 0.0);
    c1 = MatrixVectorMultiply(Q1T, b);
    vector<vector<double>> R(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            R[i][j] = A[i][j];
        }
    }
    //PrintMatrix(R);
    back_subs(R, c1);
    return c1;
}

double VectorTwoNorm(vector<double>& x) {
    double norm = 0.0;
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        norm += x[i] * x[i];
    }
    return sqrt(norm);
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
            //cout << "�� " << k + 1 << " �ε���������." << endl;
            cout << "��������Ϊ��" << k + 1 << endl;
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
            //cout << "�� " << iter + 1 << " �ε���������." << endl;
            cout << "��������Ϊ��" << iter + 1 << endl;
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
    cout << "���omegaΪ��" << bestOmega << "\nѰ����ʱ:" << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
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

vector<vector<double>> MatrixSubtraction(vector<vector<double>> x, vector<vector<double>> y) {
    int m = x.size();
    int n = x[0].size();
    vector<vector<double>> A(m, vector<double>(n));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++) {
            A[i][j] = x[i][j] - y[i][j];
        }
    }
    return A;
}

vector<vector<double>> MatrixAddition(vector<vector<double>> x, vector<vector<double>> y) {
    int m = x.size();
    int n = x[0].size();
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++) {
            A[i][j] = x[i][j] + y[i][j];
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

//�㷨5.2.1����Գ����������飺�����ݶȷ���
//x0 = ��ֵ
//r0 = b - Axo; k = 0
//while rk �� 0
//k = k + 1
//if k = 1
//p0 = r0
//else
//Beta_{ k - 2 } = r ^ T_{ k - 1 }r_{ k - 1 } / r ^ T_{ k - 2 }r_{ k - 1 }
//p_{ k - 1 } = r_{ k - 1 } + Beta_{ k - 2 }p{ k - 2 }
//end
//alpha_{ k - 1 } = r ^ T_{ k - 1 }r_{ k - 1 } / p ^ T_{ k - 1 }��A��p_{ k - 1 }
//xk = x_{ k - 1 } + alpha_{ k - 1 }p_{ k - 1 }
//rk = r_{ k - 1 } - alpha_{ k - 1 }��A��p_{ k - 1 }
//end
//x = xk
//�㷨5.3.1����Գ����������飺ʵ�ù����ݶȷ���
//x = ��ֵ
//k = 0;
//r = b - Ax;
//rou = rTr;
//while (sqrt(rou) > epsilon * ||b||2) and (k < k_max)
//    k = k + 1
//    if k = 1
//        p = r
//    else
//        Beta = rou / rou1;
//        p = r + Beta * p;
//    end
//    w = Ap ;
//    alpha = rou / pTw;
//    x = x + alpha * p;
//    r = r - alpha * w;
//    rou1 = rou;
//    rou = rTr;
//end

double dotProduct(vector<double>& v1, vector<double>& v2) {
    double result = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

void ConjugateGradient1(vector<vector<double>>& A, vector<double>& x, vector<double>& b) {
    clock_t start_time;
    start_time = clock();

    int sizeA = A.size(); //����AΪ(n-1)*(n-1)ά��
    int sizeu = int(sqrt(sizeA)); // �����sizeu��ָ������м�ľ���ά����Ϊn-1
    double h = 1.0 / (sizeu + 1);
    /*cout << sizeu << endl;*/
    vector<vector<double>> u(sizeu + 2, vector<double>(sizeu + 2, 0.0)); //��uΪ(n+1)*(n+1)ά��
    for (int i = 0; i < u.size(); i++)
    {
        u[0][i] = u[i][0] = i * i * h * h;
        u[u.size() - 1][i] = u[i][u.size() - 1] = i * i * h * h + (u.size() - 1) * (u.size() - 1) * h * h;
    }
    //PrintMatrix(u);
    double epsilon = 1e-7;
    int k_max = 1000;
    //vector<double>x(n, 0.0);
    vector<double> r(sizeA), p(sizeA), w(sizeA);
    double alpha, Beta, rou, rou1;
    alpha = 10000;
    r = b;
    rou = dotProduct(r, r);
    rou1 = rou;
    p = r;
    //cout << VectorInfinityNorm(p) * alpha << endl;
    int k = 0;
    while ((abs(alpha) * VectorInfinityNorm(p) > epsilon) && k < k_max) {
        //while (sqrt(rou) > epsilon * sqrt(dotProduct(b, b)) && k < k_max) {
        k++;
        if (k == 1) {
            p = r;
        }
        else {
            Beta = rou / rou1;
            for (int i = 0; i < sizeA; ++i) {
                p[i] = r[i] + Beta * p[i];
            }
        }
        for (int i = 0; i < sizeA; ++i) {
            w[i] = 0.0;
            for (int j = 0; j < sizeA; ++j) {
                w[i] += A[i][j] * p[j];
            }
        }
        alpha = rou / dotProduct(p, w);
        //cout << "alpha = " << alpha << endl;
        for (int i = 0; i < sizeA; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * w[i];
        }
        rou1 = rou;
        rou = dotProduct(r, r);
        //cout << "alpha * ||p||�� = " << VectorInfinityNorm(p) * alpha << endl;
    }
    cout << "�����ݶȷ�\n��������Ϊ��" << k << "\n������Ϊ��\n";
    //PrintVector(x);
    for (int i = 0; i < sizeu;i++) {
        for (int j = 0; j < sizeu; j++) {
            u[i + 1][j + 1] = x[i * sizeu + j];
            /*b[i * sizeB + j] = sin((i + 1) * (j + 1) * h * h);
            b[i * sizeB + j] *= h * h / 4;*/
        }
    }
    PrintMatrix(u);
    cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}

void ConjugateGradient2(vector<vector<double>>& A, vector<double>& x, vector<double>& b) {
    clock_t start_time;
    start_time = clock();

    int n = b.size();
    double epsilon = 1e-6;
    int k_max = 1000;
    //vector<double>x(n, 0.0);
    vector<double> r(n), p(n), w(n);
    double alpha, Beta, rou, rou1;
    r = b;
    rou = dotProduct(r, r);
    rou1 = rou;

    int k = 0;
    while (sqrt(rou) > epsilon * sqrt(dotProduct(b, b)) && k < k_max) {
        k++;
        if (k == 1) {
            p = r;
        }
        else {
            Beta = rou / rou1;
            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + Beta * p[i];
            }
        }
        for (int i = 0; i < n; ++i) {
            w[i] = 0.0;
            for (int j = 0; j < n; ++j) {
                w[i] += A[i][j] * p[j];
            }
        }
        alpha = rou / dotProduct(p, w);
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * w[i];
        }
        rou1 = rou;
        rou = dotProduct(r, r);
    }
    /*cout << "�����ݶȷ�\n��������Ϊ��" << k << "\n������Ϊ��\n";
    PrintVector(x);
    cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;*/
}

vector<vector<double>> generateMatrixA(int n) {
    double h = 1.0 / n;
    double diagValue = 1.0 + h * h / 4.0;
    double subDiagValue = -1.0 / 4.0;

    // Create the matrix A with the described properties
    vector<vector<double>> A(n - 1, vector<double>(n - 1, 0.0));

    for (int i = 0; i < n - 1; ++i) {
        A[i][i] = diagValue;  // Diagonal elements
        if (i > 0) {
            A[i][i - 1] = subDiagValue;  // Sub-diagonal elements
            A[i - 1][i] = subDiagValue;  // Super-diagonal elements
        }
    }

    return A;
}

vector<vector<double>> equation_to_matrix(vector<double>& equation) {
    int m = equation.size();
    vector<vector<double>> A(m, vector<double>(m));
    for (int i = 0; i < m; i++) {
        A[0][i] = -equation[i];
        if (i < m - 1) {
            A[i + 1][i] = 1;
        }
    }
    //PrintMatrix(A);
    return A;
}

double maxModulus(vector<double>& vec) {
    int n = vec.size();
    double elem = 0.0;
    for (int i = 0; i < n; i++) {
        if (abs(vec[i]) > abs(elem)) {
            elem = vec[i];
        }
    }
    return elem;
}

double powerMethod(vector<vector<double>>& matrix) {
    clock_t start_time;
    start_time = clock();

    int maxIterations = 1000; // ����������
    double epsilon = 1e-6; // �ݲ�
    int n = matrix.size();
    vector<double> x(n); 
    for (int i = 0; i < n; i++) {
        x[i] = 1.0 / (i + 1);
    } // ��ʼ����x, ��ֹx�պ�����������
    //PrintVector(x);
    vector<double> prevX(n);
    vector<double> diff(n);
    double lambda = 0.0;
    int iteration = 0;
    double error = 1.0;

    while (iteration < maxIterations && error > epsilon) {
        prevX = x;
        double normX = maxModulus(x);
        for (int i = 0; i < n; ++i) {
            x[i] = 0.0;
            for (int j = 0; j < n; ++j) {
                x[i] += matrix[i][j] * prevX[j];
            }
        }
        //PrintVector(x);
        lambda = maxModulus(x);
        for (int i = 0; i < n; ++i) {
            x[i] /= lambda; // ��һ������
        }
        diff = VectorSubtraction(x, prevX);
        error = VectorInfinityNorm(diff);
        iteration++;
    }

    cout << "��������Ϊ��" << iteration << endl;
    //PrintVector(x);
    //cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    return lambda;
}

//�㷨6.4.1��������Hessenberg�ֽ⣺Householder�任����
//for k = 1 : n - 2
//    [v, beta] = house(A(k + 1 : n, k))
//    A(k + 1: n , k : n) = (I - beta��vv^T)A(k + 1 : n, k : n)
//    A(1 : n, k + 1 : n) = A(1 : n, k + 1 : n)(I - beta��vv ^ T)
//end

vector<vector<double>> getSubMatrix(vector<vector<double>>& A, int startRow, int endRow, int startCol, int endCol) {
    if (startRow < 0) {
        startRow = 0;
    }
    if (endRow >= A.size()) {
        endRow = A.size() - 1;
    }
    if (startCol < 0) {
        startCol = 0;
    }
    if (endCol >= A.size()) {
        endCol = A.size() - 1;
    }
    vector<vector<double>> submatrix;
    for (int i = startRow; i <= endRow; ++i) {
        vector<double> row;
        for (int j = startCol; j <= endCol; ++j) {
            row.push_back(A[i][j]);
        }
        submatrix.push_back(row);
    }
    return submatrix;
}

void setSubMatrix(vector<vector<double>>& A, vector<vector<double>>& submatrix, int startRow, int startCol) {
    int numRows = submatrix.size();
    int numCols = submatrix[0].size();
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            A[startRow + i][startCol + j] = submatrix[i][j];
        }
    }
}

void hessenberg(vector<vector<double>>& A) {
    int n = A.size();
    /*cout << "ԭ����Ϊ��" << endl;
    PrintMatrix(A);*/
    double beta = 0.0;
    //int k = 0;
    for (int k = 0; k < n - 2; k++) {
        //cout << "k = " << k + 1 << endl;
        vector<double> x;
        vector<double> v;
        for (int i = k + 1; i < n;i++) {
            x.push_back(A[i][k]);
        }
        //PrintVector(x);
        // Ӧ�� Householder �任�õ����� v �� beta
        house(x, v, beta);
        //PrintVector(v);
        int m = x.size();
        //cout << "m = " << m << endl;
        vector<vector<double>> identity(m, vector<double>(m, 0.0));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                identity[i][j] = (i == j) ? 1.0 - beta * v[i] * v[j] : -beta * v[i] * v[j];
            }
        }
        //PrintMatrix(identity);
        vector<vector<double>> submatrix = getSubMatrix(A, k + 1, n - 1, k, n - 1);
        //PrintMatrix(submatrix);
        submatrix = MultiplyMatrices(identity, submatrix);
        //PrintMatrix(submatrix);
        setSubMatrix(A, submatrix, k + 1, k);
        /*cout << "A:" << endl;
        PrintMatrix(A);*/
        submatrix = getSubMatrix(A, 0, n - 1, k + 1, n - 1);
        /*PrintMatrix(submatrix);*/
        submatrix = MultiplyMatrices(submatrix, identity);
        /*PrintMatrix(submatrix);*/
        setSubMatrix(A, submatrix, 0, k + 1);
        /*cout << "A:" << endl;
        PrintMatrix(A);*/
    }
    /*cout << "��Hessenberg�ֽ�õ�����:" << endl;
    PrintMatrix(A);*/
}

//�㷨6.4.2��˫�ز�λ�Ƶ�QR������
//m = n - 1
//s = H(m, m) + H(n, n)
//t = H(m, m)H(n, n) - H(m, n)H(n, m)
//x = H(1, 1)H(1, 1) + H(1, 2)H(2, 1) - sH(1, 1) + t
//y = H(2, 1)(H(1, 1) + H(2, 2) - s)
//z = H(2, 1)H(3, 2)
//for k = 0 : n - 3
//    [v, beta] = house([x, y, z]T��
//    q = max{1,k}
//    H(k + 1:k + 3, q : n) = (I - beta��vvT)H(k + 1:k + 3, q : n)
//    r = min{ k + 4, n }
//    H(1:r, k + 1 : k + 3) = H(1:r, k + 1 : k + 3)(I - beta��vvT)
//    x = H(k + 2, k + 1)
//    y = H(k + 3, k + 1)
//    if k < n - 3
//        z = H(k + 4, k + 1)
//    end
//end
//[v, beta] = house([x, y]T)
//H(n - 1:n, n - 2 : n) = (I - beta��vvT)H(n - 1:n, n - 2 : n)
//H(1 : n, n - 1 : n) = H(1 : n, n - 1 : n)(I - beta��vvT)

void doubleShiftQR(vector<vector<double>>& H) {
    int n = H.size();
    if (n == 2) {
        double a = H[0][0];
        double b = H[0][1];
        double c = H[1][0];
        double d = H[1][1];

        // ������ת�Ƕ�
        double s = hypot(a, c);
        double cos_val = a / s;
        double sin_val = c / s;
        // Ӧ�� Givens �任��ԭʼ����
        for (int i = 0; i < n; ++i) {
            double temp1 = H[0][i];
            double temp2 = H[1][i];
            H[0][i] = cos_val * temp1 + sin_val * temp2;
            H[1][i] = -sin_val * temp1 + cos_val * temp2;
        }

    }
    else {
        int m = n - 1;
        double s, t, x, y, z;
        s = H[m - 1][m - 1] + H[n - 1][n - 1];
        t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1] * H[n - 1][m - 1];
        x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
        y = H[1][0] * (H[0][0] + H[1][1] - s);
        z = H[1][0] * H[2][1];
        for (int k = 0; k <= n - 3; ++k) {
            vector<double> temp, v;
            temp = { x, y, z };
            double beta;
            house(temp, v, beta);
            int q = max(1, k);
            int r = min(k + 4, n);

            int N = v.size();
            vector<vector<double>> identity(N, vector<double>(N, 0.0));
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    identity[i][j] = (i == j) ? 1.0 - beta * v[i] * v[j] : -beta * v[i] * v[j];
                }
            }
            vector<vector<double>> submatrix1 = getSubMatrix(H, k, k + 2, q - 1, n - 1);
            submatrix1 = MultiplyMatrices(identity, submatrix1);
            setSubMatrix(H, submatrix1, k, q - 1);
            vector<vector<double>> submatrix2 = getSubMatrix(H, 0, r - 1, k, k + 2);
            submatrix2 = MultiplyMatrices(submatrix2, identity);
            setSubMatrix(H, submatrix2, 0, k);
            x = H[k + 1][k];
            y = H[k + 2][k];
            if (k < n - 3) {
                z = H[k + 3][k];
            }
        }

        vector<double> temp, v;
        double beta;
        temp = { x, y };
        house(temp, v, beta);
        int N = v.size();
        vector<vector<double>> identity(N, vector<double>(N, 0.0));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                identity[i][j] = (i == j) ? 1.0 - beta * v[i] * v[j] : -beta * v[i] * v[j];
            }
        }

        vector<vector<double>> submatrix1 = getSubMatrix(H, n - 2, n - 1, n - 3, n - 1);
        submatrix1 = MultiplyMatrices(identity, submatrix1);
        setSubMatrix(H, submatrix1, n - 2, n - 3);
        vector<vector<double>> submatrix2 = getSubMatrix(H, 0, n - 1, n - 2, n - 1);
        submatrix2 = MultiplyMatrices(submatrix2, identity);
        setSubMatrix(H, submatrix2, 0, n - 2);
        /*cout << "˫�ز�λ�Ƶ�QR�����õ��ľ���Ϊ��" << endl;
        PrintMatrix(H);*/
    }
    

}

//����˫�ز�λ�Ƶ�QR����
//void doubleShiftQR(vector<vector<double>>& H) {
//    int n = H.size();
//    int m = n - 1;
//    double s, t, x, y, z;
//    s = H[m - 1][m - 1] + H[n - 1][n - 1];
//    t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1] * H[n - 1][m - 1];
//    x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
//    y = H[1][0] * (H[0][0] + H[1][1] - s);
//    z = H[1][0] * H[2][1];
//    cout << "m = " << m << endl;
//    cout << "n = " << n << endl;
//    cout << "s = " << s << endl;
//    cout << "t = " << t << endl;
//    cout << "x = " << x << endl;
//    cout << "y = " << y << endl;
//    cout << "z = " << z << endl;
//
//    //int k = n - 3;
//    for (int k = 0; k <= n - 3; ++k) {
//        cout << "\nk = " << k << endl;
//        vector<double> temp, v;
//        temp = { x, y, z };
//        cout << "temp = ";
//        PrintVector(temp);
//        double beta;
//        house(temp, v, beta);
//        cout << "v = ";
//        PrintVector(v);
//        cout << "beta: " << beta << endl;
//        int q = max(1, k);
//        int r = min(k + 4, n);
//        cout << "q = " << q << endl;
//        cout << "r = " << r << endl;
//
//        int N = v.size();
//        vector<vector<double>> identity(N, vector<double>(N, 0.0));
//        for (int i = 0; i < N; ++i) {
//            for (int j = 0; j < N; ++j) {
//                identity[i][j] = (i == j) ? 1.0 - beta * v[i] * v[j] : -beta * v[i] * v[j];
//            }
//        }
//        cout << "\nidentity = " << endl;
//        PrintMatrix(identity);
//        vector<vector<double>> submatrix1 = getSubMatrix(H, k, k + 2, q - 1, n - 1);
//        cout << "submatrix1 = " << endl;
//        PrintMatrix(submatrix1);
//        submatrix1 = MultiplyMatrices(identity, submatrix1);
//        cout << "submatrix1 = " << endl;
//        PrintMatrix(submatrix1);
//        setSubMatrix(H, submatrix1, k, q - 1);
//        cout << "H = " << endl;
//        PrintMatrix(H);
//
//        vector<vector<double>> submatrix2 = getSubMatrix(H, 0, r - 1, k, k + 2);
//        cout << "submatrix2 = " << endl;
//        PrintMatrix(submatrix2);
//        submatrix2 = MultiplyMatrices(submatrix2, identity);
//        cout << "submatrix2 = " << endl;
//        PrintMatrix(submatrix2);
//        setSubMatrix(H, submatrix2, 0, k);
//        cout << "H = " << endl;
//        PrintMatrix(H);
//        x = H[k + 1][k];
//        y = H[k + 2][k];
//        if (k < n - 3) {
//            z = H[k + 3][k];
//        }
//    }
//
//    vector<double> temp, v;
//    double beta;
//    temp = { x, y };
//    house(temp, v, beta);
//    int N = v.size();
//    vector<vector<double>> identity(N, vector<double>(N, 0.0));
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < N; ++j) {
//            identity[i][j] = (i == j) ? 1.0 - beta * v[i] * v[j] : -beta * v[i] * v[j];
//        }
//    }
//    cout << "\nidentity = " << endl;
//    PrintMatrix(identity);
//
//    vector<vector<double>> submatrix1 = getSubMatrix(H, n - 2, n - 1, n - 3, n - 1);
//    cout << "submatrix1 = " << endl;
//    PrintMatrix(submatrix1);
//    submatrix1 = MultiplyMatrices(identity, submatrix1);
//    cout << "submatrix1 = " << endl;
//    PrintMatrix(submatrix1);
//    setSubMatrix(H, submatrix1, n - 2, n - 3);
//    cout << "H = " << endl;
//    PrintMatrix(H);
//
//    vector<vector<double>> submatrix2 = getSubMatrix(H, 0, n - 1, n - 2, n - 1);
//    cout << "submatrix2 = " << endl;
//    PrintMatrix(submatrix2);
//    submatrix2 = MultiplyMatrices(submatrix2, identity);
//    cout << "submatrix2 = " << endl;
//    PrintMatrix(submatrix2);
//    setSubMatrix(H, submatrix2, 0, n - 2);
//    cout << "˫�ز�λ�Ƶ�QR�����õ��ľ���Ϊ��" << endl;
//    PrintMatrix(H);
//}

bool areEigenvaluesReal(vector<vector<double>>& A) {
    double a = A[0][0];
    double b = A[0][1];
    double c = A[1][0];
    double d = A[1][1];
    double trace = a + d;
    double det = a * d - b * c;
    double discriminant = trace * trace - 4 * det;
    //cout << discriminant << endl;
    /*if (discriminant >= 0) {
        cout << "��������ֵΪʵ��." << endl;
    }
    else {
        cout << "��������ֵΪ����." << endl;
    }*/
    return (discriminant >= 0); // ����б�ʽ���ڵ���0������ֵΪʵ��
}

void zeroing(vector<vector<double>>& A, double u) {
    //double u = 100;
    //double u = 1e-10;
    int n = A.size();
    for (int i = 1; i < n;i++) {
        if (abs(A[i][i - 1]) <= (abs(A[i][i]) + abs(A[i - 1][i - 1])) * u) {
            A[i][i - 1] = 0.0;
        }
    }
}

bool isQuasi(vector<vector<double>>& A) {
    double a = A[0][0];
    double b = A[0][1];
    double c = A[1][0];
    double d = A[1][1];
    if (c == 0 or !(areEigenvaluesReal(A))) {
        //cout << "�������׼��������Խ�Ԫ." << endl;
        return true;
    }
    else {
        //cout << "����鲻��׼��������Խ�Ԫ." << endl;
        return false;
    }
}

/*int quasi(vector<vector<double>>& H) {  //Ĭ������ľ�������Hessenberg����
    int n = H.size();
    //int hess = 0;
    vector<int> hess = { }; // �洢����������H33ά��
    vector<vector<double>> submatrix_up;
    vector<vector<double>> submatrix_down;
    for (int m = 0; m <= n; m++) {
        if (m == 0) {
            cout << "\nm = " << m << endl;
            submatrix_up = getSubMatrix(H, n - m - 2, n - m - 1, n - m - 2, n - m - 1);
            cout << "submatrix_up = " << endl;
            PrintMatrix(submatrix_up);
            //cout << "����λ�õ�ֵΪ��" << H[n - m][n - m - 1] << endl;
            if (!isQuasi(submatrix_up)) {
                hess.push_back(m);
                //return m;
            }
        }
        if (m == 1) {
            cout << "\nm = " << m << endl;
            submatrix_up = getSubMatrix(H, n - m - 2, n - m - 1, n - m - 2, n - m - 1);
            cout << "submatrix_up = " << endl;
            PrintMatrix(submatrix_up);
            cout << "����λ�õ�ֵΪ��" << H[n - m][n - m - 1] << endl;
            if (!isQuasi(submatrix_up) && H[n - m][n - m - 1] == 0) {
                hess.push_back(m);
                //return m;
            }
        }
        // ����1<m<n-1���������Ҫ�������������Ӿ���
        if (1 < m && m < n - 1) {
            cout << "\nm = " << m << endl;
            submatrix_up = getSubMatrix(H, n - m - 2, n - m - 1, n - m - 2, n - m - 1);
            cout << "submatrix_up = " << endl;
            PrintMatrix(submatrix_up);
            cout << isQuasi(submatrix_up) << endl;
            submatrix_down = getSubMatrix(H, n - m, n - m + 1, n - m, n - m + 1);
            cout << "submatrix_down = " << endl;
            PrintMatrix(submatrix_down);
            cout << isQuasi(submatrix_down) << endl;
            cout << "����λ�õ�ֵΪ��" << H[n - m][n - m - 1] << endl;


            if (!isQuasi(submatrix_up) && H[n - m][n - m - 1] == 0 && isQuasi(submatrix_down)) {
                hess.push_back(m);
                //return m;
            }
        }
        if (m == n - 1) {
            cout << "\nm = " << m << endl;
            //cout << "\nhess = " << hess << endl;
            submatrix_down = getSubMatrix(H, n - m, n - m + 1, n - m, n - m + 1);
            cout << "submatrix_down = " << endl;
            PrintMatrix(submatrix_down);
            cout << isQuasi(submatrix_down) << endl;
            cout << "����λ�õ�ֵΪ��" << H[n - m][n - m - 1] << endl;
            if (isQuasi(submatrix_down) && H[n - m][n - m - 1] == 0) {
                hess.push_back(m + 1);
                //return m + 1;
            }
        }
        if (m == n) {
            cout << "\nm = " << m << endl;
            //cout << "\nhess = " << hess << endl;
            submatrix_down = getSubMatrix(H, n - m, n - m + 1, n - m, n - m + 1);
            cout << "submatrix_down = " << endl;
            PrintMatrix(submatrix_down);
            cout << isQuasi(submatrix_down) << endl;
            //cout << "����λ�õ�ֵΪ��" << H[n - m][n - m - 1] << endl;
            if (isQuasi(submatrix_down) && H[n - m][n - m - 1] == 0) {
                hess.push_back(m);
                //return m + 1;
            }
        }

    }
    cout << "����������H33ά��: ";
    for (int i = 0; i < hess.size(); ++i) {
        cout << hess[i] << " "; // ������vectorԪ��
    }
    cout << endl;
    return n;
}*/

bool quasii(vector<vector<double>>& H, int m) {
    int n = H.size();
    vector<vector<double>> submatrix_down;
    //cout << "m = " << m << endl;
    //cout << H[n - m][n - m - 1] << endl;
    if (m < n && H[n - m][n - m - 1]) {
        return false;
    }

    if (m == 1) {
        if (H[n - m][n - m - 1]) {
            return false;
        }
        else {
            return true;
        }
    }
    if (m == 2) {
        if (H[n - m + 1][n - m]) {
            submatrix_down = getSubMatrix(H, n - m, n - m + 1, n - m, n - m + 1);
            PrintMatrix(submatrix_down);
            //cout << isQuasi(submatrix_down) << endl;
            if (isQuasi(submatrix_down)) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            return(quasii(H, m - 1));
        }
    }
    if (H[n - m + 1][n - m] == 0) {
        return(quasii(H, m - 1));
    }
    else {
        submatrix_down = getSubMatrix(H, n - m, n - m + 1, n - m, n - m + 1);
        //PrintMatrix(submatrix_down);
        //cout << isQuasi(submatrix_down) << endl;
        if (isQuasi(submatrix_down)) {
            return quasii(H, m - 2);
        }
        else {
            return false;
        }
    }    
}

int quasi(vector<vector<double>>& H) {
    int n = H.size();
    for (int k = n; k >= 1; k--) {
        if (quasii(H, k)) {
            //cout << "k = " << k << "ʱ������������" << endl;
            return k;
        }
    }
    return 0;
}

bool isHessenberg(vector<vector<double>>& A) {
    int n = A.size();
    double flag = 1.0;
    for (int i = 1; i < n; i++) {
        flag *= A[i][i - 1];
    }
    /*return(!(flag==0));*/
    return(flag);
}

int IrredHessenberg(vector<vector<double>>& H, int quasi) {
    int n = H.size();
    if (quasi == n) {
        return 0;
    }
    int end = n - 1 - quasi;
    vector<vector<double>> submatrix;
    //cout << end << endl;
    for (int k = 1; k <= end; k++) {
        submatrix = getSubMatrix(H, end - k, end, end - k, end);
        /*cout << "submatrix = " << endl;
        PrintMatrix(submatrix);*/
        if (!isHessenberg(submatrix)) {
            return submatrix.size() - 1;
        }
    }
    return submatrix.size();
}

vector<vector<double>> getHessenberg(vector<vector<double>>& H, int m) {
    int n = H.size();
    int l = IrredHessenberg(H, m);
    int end = n - m - 1;
    vector<vector<double>> submatrix;
    submatrix = getSubMatrix(H, end - l + 1, end, end - l + 1, end);
    return submatrix;
}

void setHessenberg(vector<vector<double>>& H, vector<vector<double>>& A, int m) {
    int n = H.size();
    int l = A.size();
    int end = n - m - 1;
    setSubMatrix(H, A, end - l + 1, end - l + 1);
}

void eigenvalues2D(vector<vector<double>>& A) {
    double a = A[0][0];
    double b = A[0][1];
    double c = A[1][0];
    double d = A[1][1];
    double trace = a + d;
    double det = a * d - b * c;
    double discriminant = trace * trace - 4 * det;
    if (discriminant > 0) {
        double lambda1 = (trace + sqrt(discriminant)) / 2;
        double lambda2 = (trace - sqrt(discriminant)) / 2;
        //cout << "����ֵΪ��" << lambda1 << " and " << lambda2 << endl;
        cout << lambda1 << "\n" << lambda2 << endl;
    }
    else if (discriminant == 0) {
        double lambda = trace / 2;
        //cout << "����ֵΪ��" << lambda << endl;
        cout << lambda << endl;
    }
    else {
        double realPart = trace / 2;
        double imaginaryPart = sqrt(-discriminant) / 2;
        //cout << "����ֵΪ��" << realPart << " �� " << imaginaryPart << "i" << endl;
        cout << realPart << " �� " << imaginaryPart << "i" << endl;
    }
}

void prEigens(vector<vector<double>>& A) {
    int n = A.size();
    int k = 0;
    vector<vector<double>> B;
    //cout << "���õ�������ֵΪ��" << endl;
    while (k < n) {
        //cout << k << endl;
        if (k == n-1) {
            cout << A[k][k] << endl;
            break;
        }
        if (A[k + 1][k]) {
            B = getSubMatrix(A, k, k+1, k, k + 1);
            //PrintMatrix(B);
            eigenvalues2D(B);
            k += 2;
        }
        else {
            cout << A[k][k] << endl;
            k++;
        }
    }
}

int implicitQR(vector<vector<double>>& A) {
    /*cout << "ԭ����Ϊ��" << endl;
    PrintMatrix(A);*/
    int n = A.size();
    //double u = 1e-6; // �ڶ�С����
    double u = 1e-4; // ����С����
    hessenberg(A); //������A��Hessenberg��
    int m = quasi(A); // H33��ά��
    int l = IrredHessenberg(A, m); // H22��ά��
    vector<vector<double>> H22;
    int k = 0;
    while (m != n && k <= 1000) {
        //cout << "k = " << k << endl;
        H22 = getHessenberg(A, m);
        //PrintMatrix(H22);
        doubleShiftQR(H22);
        //PrintMatrix(H22);
        setHessenberg(A, H22, m);
        //PrintMatrix(A);
        zeroing(A, u);
        //PrintMatrix(A);
        m = quasi(A);
        //cout << "m = " << m << endl;
        k++;
    }
    return k;
    //cout << "\n��������Ϊ��" << k << endl;
    /*cout << "��ʽQR�����õ�����Ϊ��" << endl;
    PrintMatrix(A);*/
    /*cout << "\n������Ľ�Ϊ��" << endl;
    prEigens(A);*/
}

vector<vector<double>> jacobiClassic(vector<vector<double>>& A, int p, int q) {
    double alpha_pp = A[p][p];
    double alpha_qq = A[q][q];
    double alpha_pq = A[p][q];
    vector<vector<double>> J(A.size(), vector<double>(A.size(), 0.0));
    for (int k = 0; k < A.size(); ++k) {
        J[k][k] = 1.0;
    }
    if (alpha_pq == 0) {
        A[p][q] = A[q][p] = 0.0;
        return J;
    }
    double tau = (alpha_qq - alpha_pp) / (2 * alpha_pq);
    double t = (tau >= 0) ? (1.0 / (tau + sqrt(1 + tau * tau))) : (-1.0 / (-tau + sqrt(1 + tau * tau)));
    double c = 1.0 / sqrt(1 + t * t);
    double s = t * c;
    for (int i = 0; i < A.size(); ++i) {
        if (i != p && i != q) {
            double temp = A[i][p];
            A[i][p] = c * temp - s * A[i][q];
            A[p][i] = A[i][p];
            A[i][q] = s * temp + c * A[i][q];
            A[q][i] = A[i][q];
        }
    }
    A[p][p] = c * c * alpha_pp - 2 * s * c * alpha_pq + s * s * alpha_qq;
    A[q][q] = s * s * alpha_pp + 2 * s * c * alpha_pq + c * c * alpha_qq;
    A[p][q] = A[q][p] = 0.0;
    J[p][p] = J[q][q] = c;
    J[p][q] = -s;
    J[q][p] = s;
    return J;
}

double offDiagNorm(vector<vector<double>>& A) {
    double norm = 0.0;
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                norm += A[i][j] * A[i][j];
            }
        }
    }
    return sqrt(norm);
}

bool passingThreshold(vector<vector<double>>& A, double delta) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                if (abs(A[i][j]) > delta) {
                    return true;
                }
            }
        }
    }
    return false;
}

vector<vector<double>> thresholdJacobi(vector<vector<double>>& A) {
    int n = A.size();
    double delta = offDiagNorm(A);
    double epsilon = 1e-6;
    vector<vector<double>> Q(n, vector<double>(n, 0.0));
    for (int k = 0; k < n; ++k) {
        Q[k][k] = 1.0;
    }
    int k = 0;
    while (offDiagNorm(A) > epsilon && k<100) {
        while (passingThreshold(A, delta)) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (i != j) {
                        if (abs(A[i][j]) > delta) {
                            //cout << k << endl;
                            vector<vector<double>> J = jacobiClassic(A, i, j);
                            Q = MultiplyMatrices(Q, J);
                        }
                    }
                }
            }
        }
        delta /= n;
        k++;
    }
    //cout << "ѭ������Ϊ" << k;
    //PrintMatrix(Q);
    return Q;
}

//�㷨7.4.1������������
//x = [a1, a2, ��, an]
//y = [0, b2, ��, bn]
//s = 0; q = x(1) - miu
//for k = 1 : n
//    if q < 0
//        s = s + 1
//    end
//    if k < n
//        if q = 0
//            q = |y(k + 1)|u
//        end
//        q = x(k + 1) - miu - y(k + 1)^2 / q
//    end
//end

int reversals(vector<vector<double>>& A, double miu) {
    int n = A.size();
    vector<double> x;
    for (int i = 0; i < n;i++) {
        x.push_back(A[i][i]);
    }
    //PrintVector(x);
    vector<double> y = { 0 };
    for (int i = 0; i < n - 1;i++) {
        y.push_back(A[i][i + 1]);
    }
    //PrintVector(y);
    int s = 0;
    double u = 1e-3;
    double q = x[0] - miu;
    for (int k = 0; k < n; ++k) {
        if (q < 0) {
            s += 1;
        }

        if (k < n - 1) {
            if (q == 0) {
                q = abs(y[k + 1]) * u;
            }
            q = x[k + 1] - miu - y[k + 1] * y[k + 1] / q;
        }
    }
    return s;
}

double dichEigenvalue(vector<vector<double>>& A, int m) {
    double l = -MatrixInfinityNorm(A);
    double u = -l;
    double r = (l + u) / 2;
    double tolerance = 1e-6;
    while (abs(l - r) >= tolerance) {
        if (reversals(A, r) >= m) {
            u = r;
            r = (l + u) / 2;
        }
        else {
            l = r;
            r = (l + u) / 2;
        }
    }
    return r;
}

vector<double> VectorAddition(vector<double>x, vector<double>y) {
    int n = x.size();
    vector<double>a(n);
    for (int i = 0; i < n; i++)
    {
        a[i] = x[i] + y[i];
    }
    return a;
}

vector<double> inversePowerMethod(vector<vector<double>>& A, double lambda) {
    clock_t start_time;
    start_time = clock();
    int maxIterations = 100; // ����������
    double epsilon = 1e-4; // �ݲ�
    int iteration = 0;
    double error1 = 1.0;
    double error2 = 1.0;

    vector<vector<double>> matrix = A;
    int n = matrix.size();
    //PrintMatrix(matrix);
    vector<double> v(n);
    vector<double> z(n);

    random_device rd;  // �������
    mt19937 gen(rd());  // Mersenne Twister ������
    uniform_real_distribution<double> dis(0.0, 1.0);  // ���� [0.0, 1.0) ��Χ�ڵ���� double ��
    for (int i = 0; i < n; i++) {
        matrix[i][i] -= lambda;
        //z[i] = 1.0 / (i + 1);
        z[i] = dis(gen);
    }
    //PrintMatrix(matrix);
    //PrintVector(z);
    double normZ = VectorTwoNorm(z);
    for (int i = 0; i < n; ++i) {
        z[i] = z[i] / normZ;
    }
    //PrintVector(z);
    vector<double> prevz(n);
    vector<double> diff1(n);
    vector<double> diff2(n);


    while (iteration < maxIterations && error1 > epsilon && error2 > epsilon) {
        prevz = z;
        ConjugateGradient2(matrix, v, z);
        /*cout << "v: ";
        PrintVector(v);*/
        double normV = VectorTwoNorm(v);
        for (int i = 0; i < n; ++i) {
            z[i] = v[i] / normV;
        }
        /*cout << "z: ";
        PrintVector(z);
        cout << "prevz: ";
        PrintVector(prevz);*/
        diff1 = VectorSubtraction(z, prevz);
        diff2 = VectorAddition(z, prevz);
        /*cout << "diff2: ";
        PrintVector(diff2);*/
        error1 = VectorInfinityNorm(diff1);
        error2 = VectorInfinityNorm(diff2);
        iteration++;
    }

    cout << "��������Ϊ��" << iteration << endl;
    //PrintVector(z);
    //cout << "����ʱ��: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    return z;
}

vector<vector<double>> intMatrix(int rows, int cols) {
    vector<vector<double>> matrix(rows, vector<double>(cols, 0));
    double value = 1.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = value++;
        }
    }
    return matrix;
}

vector<vector<double>> randomMatrix(int rows, int cols) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(1.0, 10.0); // Adjust the range as needed

    // Create a matrix with random values
    vector<vector<double>> matrix(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = dis(gen);
        }
    }
    return matrix;
}

vector<vector<double>> idMatrix(int n) {
    vector<vector<double>> matrix(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        matrix[i][i] = 1;
    }
    return matrix;
}

void transpose(vector<vector<double>>& A) {
    int rows = A.size();
    int cols = A[0].size();

    // Resize the transposed matrix
    vector<vector<double>> transposedMatrix(cols, vector<double>(rows, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposedMatrix[j][i] = A[i][j];
        }
    }
    A = transposedMatrix;
}

vector<vector<double>> Householder(vector<double>& v, double beta) {
    int m = v.size();
    vector<vector<double>> H(m, vector<double>(m, 0.0));
    for (int i = 0;i < m;i++) {
        for (int j = 0;j < m;j++) {
            H[i][j] = (i == j) ? 1.0 - beta * v[i] * v[j] : -beta * v[i] * v[j];
        }
    }
    return H;
}

void bidiag(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V) {
    int m = A.size();
    int n = A[0].size();
    U = idMatrix(m);
    V = idMatrix(n);
    for (int k = 0; k < n; ++k) {
        //[u, beta] = house(A��k : m, k))
        vector<double> u(n - k, 0.0);
        vector<double> v(m - k, 0.0);
        double beta;
        vector<double> column(m - k);
        for (int i = k; i < m; ++i) {
            column[i - k] = A[i][k];
        }
        /*cout << "A(k:m,k):" << endl;
        PrintVector(column);*/
        house(column, v, beta);
        /*cout << "v:" << endl;
        PrintVector(v);
        cout << "beta:" << beta << "\n" << endl;*/
        //uT = (beta��vT)A��k : m, k : n)
        vector<vector<double>> H1 = Householder(v, beta);
        vector<vector<double>> id1 = idMatrix(m);
        setSubMatrix(id1, H1, k, k);
        /*cout << "H1:" << endl;
        PrintMatrix(H1);
        cout << "id1:" << endl;
        PrintMatrix(id1);*/
        A = MultiplyMatrices(id1, A);
        U = MultiplyMatrices(id1, U);
        if (k < n - 2) {
            //[v, beta] = house(A(k, k + 1 : n))
            vector<double> column(n - k - 1);
            for (int i = k + 1; i < n; ++i) {
                column[i - k - 1] = A[k][i];
            }
            /*cout << "A(k,k+1:n):" << endl;
            PrintVector(column);*/
            house(column, v, beta);
            /*cout << "v:" << endl;
            PrintVector(v);
            cout << "beta:" << beta << "\n" << endl;*/
            vector<vector<double>> H2 = Householder(v, beta);
            vector<vector<double>> id2 = idMatrix(n);
            setSubMatrix(id2, H2, k + 1, k + 1);
            /*cout << "H2:" << endl;
            PrintMatrix(H2);
            cout << "id2:" << endl;
            PrintMatrix(id2);*/
            A = MultiplyMatrices(A, id2);
            V = MultiplyMatrices(V, id2);
            /*vector<vector<double>> C = MultiplyMatrices(A, id2);
            cout << "C:" << endl;
            PrintMatrix(C);*/
        }
    }
}

void zeroing2(vector<vector<double>>& A, double epsilon) {
    //epsilon = 0.3;
    int n = A.size();
    for (int i = 0; i < n - 1;i++) {
        if (fabs(A[i][i + 1]) <= (fabs(A[i][i]) + fabs(A[i + 1][i + 1])) * epsilon) {
            A[i][i + 1] = 0.0;
        }
    }
    double temp = MatrixInfinityNorm(A) * epsilon;
    for (int i = 0; i < n;i++) {
        if (fabs(A[i][i]) < temp) {
            A[i][i] = 0.0;
        }
    } 
}

void pqq(vector<vector<double>>B, int& p, int& q) {
    int n = B.size();
    p = q = 0;
    for (q = 0; q < n - 1; ++q) {
        if (B[n - 2 - q][n - 1 - q] != 0) {
            break;
        }
    }
    if (q != n - 1) {
        int l = 1;
        while (l < n - q - 1 && B[n - q - 2 - l][n - q - 1 - l] != 0) {
            ++l;
        }
        p = n - q - l - 1;
    }
    else {
        q = n;
    }
}

void givens(double a, double b, double& c, double& s) {
    double t;
    if (b == 0) {
        c = 1;
        s = 0;
    }
    else {
        if (abs(b) > abs(a)) {
            t = a / b;
            s = 1 / sqrt(1 + t * t);
            c = s * t;
        }
        else {
            t = b / a;
            c = 1 / sqrt(1 + t * t);
            s = c * t;
        }
    }
}

vector<vector<double>> Givens(int n, int p, int q, double c, double s) {
    vector<vector<double>> A = idMatrix(n);
    A[p][p] = A[q][q] = c;
    A[p][q] = s;
    A[q][p] = -s;
    return A;
}

void GivensSVD(vector<vector<double>>& B, vector<vector<double>>& B22, int p, int i, vector<vector<double>>& U) {
    int m = B.size();
    int r = B22.size();
    double cos, sin;
    for (int j = i + 1; j < r - 1; j++) {
        givens(B22[j][j], B22[i][j], cos, sin);
        double a = B[p + i][p + j];
        double b = B[p + j][p + j];
        double c = B[p + j][p + j + 1];
        B[p + i][p + j] = -sin * a + cos * b;
        B[p + i][p + j + 1] = cos * c;
        B[p + j][p + j] = cos * a + sin * b;
        B[p + j][p + j + 1] = sin * c;
        vector<vector<double>> G = Givens(m, p + i, p + j, cos, sin);
        G = MultiplyMatrices(G, U);
    }
    int j = r - 1;
    givens(B22[j][j], B22[i][j], cos, sin);
    double a = B22[i][j];
    double b = B22[j][j];
    B22[j][j] = cos * a + sin * b;
    B22[i][j] = -sin * a + cos * b;
    setSubMatrix(B, B22, p, p);
    vector<vector<double>> G = Givens(m, p + i, p + j, cos, sin);
    U = MultiplyMatrices(G, U);
}

void WilkSVD(vector<vector<double>>& B, vector<vector<double>>& U, vector<vector<double>>& V) {
    int m = U.size(), n = V.size(), p, q;
    pqq(B, p, q);
    int r = n - p - q;
    vector<vector<double>>B22 = getSubMatrix(B, p, p + r - 1, p, p + r - 1);
    if (r == 2) { // B22Ϊ2*2άʱ
        double a, b, c, d, cos, sin;
        a = B22[0][0];
        b = B22[0][1];
        d = B22[1][1];
        cos = (a + d) / sqrt((a + d) * (a + d) + b * b);
        sin = -b / sqrt((a + d) * (a + d) + b * b);
        vector<vector<double>> G = Givens(r, 0, 1, cos, sin);
        B22 = MultiplyMatrices(G, B22);
        G = Givens(U.size(), p, p + 1, cos, sin);
        U = MultiplyMatrices(G, U);
        vector<vector<double>> J = thresholdJacobi(B22);
        setSubMatrix(B, B22, p, p);
        a = J[0][0]; b = J[0][1]; c = J[1][0]; d = J[1][1];
        vector<vector<double>> submatrix = getSubMatrix(U, p, p + 1, 0, m - 1);
        transpose(J);
        submatrix = MultiplyMatrices(J, submatrix);
        setSubMatrix(U, submatrix, p, 0);
        submatrix = getSubMatrix(V, 0, n - 1, p, p + 1);
        transpose(J);
        submatrix = MultiplyMatrices(submatrix, J);
        setSubMatrix(V, submatrix, 0, p);
        a = sign(B22[0][0]); b = sign(B22[1][1]);
        J = { {a,0},{0,b} };
        submatrix = getSubMatrix(U, p, p + 1, 0, m - 1);
        submatrix = MultiplyMatrices(J, submatrix);
        setSubMatrix(U, submatrix, p, 0);
    }
    else { // Ӧ��Wilkinsonλ�Ƶ�SVD����
        vector<vector<double>> IU, IV;
        IU = IV = idMatrix(r);
        double alpha, delta, beta, miu, gama, y, z, c, s;
        alpha = B22[r - 1][r - 1] * B22[r - 1][r - 1] + B22[r - 2][r - 1] * B22[r - 2][r - 1];
        delta = (B22[r - 2][r - 2] * B22[r - 2][r - 2] + B22[r - 3][r - 2] * B22[r - 3][r - 2] - alpha) / 2;
        beta = B22[r - 2][r - 2] * B22[r - 2][r - 1];
        miu = alpha - beta * beta / (delta + sign(delta) * sqrt(delta * delta + beta * beta));
        y = B22[0][0] * B22[0][0] - miu;
        z = B22[0][0] * B22[0][1];
        for (int k = 0; k < r - 1; k++) {
            givens(y, z, c, s);
            vector<vector<double>> G = Givens(r, k, k + 1, c, -s);
            B22 = MultiplyMatrices(B22, G);
            IV = MultiplyMatrices(IV, G);
            y = B22[k][k];
            z = B22[k + 1][k];
            givens(y, z, c, s);
            G = Givens(r, k, k + 1, c, s);
            IU = MultiplyMatrices(G, IU);
            if (k < r - 2) {
                B22 = MultiplyMatrices(G, B22);
                y = B22[k][k + 1];
                z = B22[k][k + 2];
            }
            else {
                B22 = MultiplyMatrices(G, B22);
            }
        }
        vector<vector<double>> submatrix = getSubMatrix(U, p, p + r - 1, 0, m - 1);
        submatrix = MultiplyMatrices(IU, submatrix);
        setSubMatrix(U, submatrix, p, 0);
        submatrix = getSubMatrix(V, 0, n - 1, p, p + r - 1);
        submatrix = MultiplyMatrices(submatrix, IV);
        setSubMatrix(V, submatrix, 0, p);
        setSubMatrix(B, B22, p, p);
    }
}

double absmax(vector<vector<double>> A) {
    int m = A.size(), n = A[0].size();
    double temp = A[0][0];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (abs(A[i][j]) > abs(temp)) {
                temp = A[i][j];
            }
        }
    }
    return temp;
}

void singluar(vector<vector<double>> A) {
    int n = A.size();
    vector<double> singvalues;
    for (int i = 0; i < n;i++) {
        if (A[i][i] >= 0) {
            singvalues.push_back(A[i][i]);
        }
    }
    sort(singvalues.begin(), singvalues.end());
    cout << "����ֵ��С��������Ϊ��" << endl;
    for (int i = 0;i < singvalues.size();i++) {
        cout << singvalues[i] << " ";
    }
    cout << endl;
}