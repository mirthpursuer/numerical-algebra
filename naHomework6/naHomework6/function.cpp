#include "function.h"
#include <cmath>
#include <random>
#include <iomanip> 

const double EPSILON = 1e-10;  // 定义一个很小的常数，用于判断浮点数是否接近于0

void PrintMatrix(vector<vector<double>>& mat) {
    // 设置输出精度为四位小数
    std::cout << std::fixed << std::setprecision(4);
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
    // 设置输出精度为四位小数
    std::cout << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
        // 每十个分量换一次行
        if ((i + 1) % 10 == 0) {
            std::cout << std::endl;
        }
    }
    // 打印完所有分量后换行
    std::cout << std::endl;
}

//算法1.1.1（解下三角形方程组：前代法）
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
        b[i] = (b[i] - sum);  // 对角元为1的前代法计算
    }
}

//算法1.1.2（解上三角形方程组：回代法）
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
        y[i] = (y[i] - sum);  // 对角元为1的回代法计算
    }
}

//算法1.1.3计算三角分解：Gauss消去法
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

//算法1.2.2（计算列主元三角分解：列主元Gauss消去法）
//for k = 1:n - 1
//    确定p(k≤p≤n), 使得
//    | A(p, k) | = max{ |A(i,k)|：i = k:n }
//    A(k, 1:n)<->A(p, 1:n)(交换第k行和第p行)
//    u(k) = p(记录置换矩阵P_k)
//    if A(k, k)≠0
//        A(k + 1:n, k) = A(k + 1:n, k) / A(k, k)
//        A(k+1：n, k + 1：n) = A(k+1：n, k + 1:n) - A(k + 1:n, k)*A(k, k + 1:n)
//    else
//        stop(矩阵奇异)
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
        //确定主元
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
        //交换第k行p行，第k列q列
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
            cout << "矩阵奇异" << endl;
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
    // 计算向量Q*b
    int n = v.size();
    vector<double> temp_b(n, 0);

    for (int i = 0; i < n; i++) {
        temp_b[i] = b[v[i]];
    }

    for (int i = 0; i < n; i++) {
        b[i] = temp_b[i];
    }
}

//算法1.3.1（计算Cholesky分解：平方根法）
//for k 1:n
//    A(k.k) = A(k.k)
//    A(k + 1:n, k) = A(k + 1:n，k) / A(k, k)
//    for j = k + 1 : n
//        A(j : n, j） = A(j : n, j） - A(j : n, k) * A(j, k)
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
            cerr << "矩阵非正定，无法进行Cholesky分解。" << endl;
            return;
        }

        A[i][i] = sqrt(sum);

        for (int j = i + 1; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += A[j][k] * A[i][k];
            A[j][i] = (A[j][i] - sum) / A[i][i];  // 改进的平方根法计算
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
            A[i][j] = sum;  // 计算矩阵D*L^T
        }
    }
}

//算法 2.5.1 估计矩阵的1范数：优化法
//k = 1
//while k = 1
//    w = Bx; v = sign(w); z = BTv
//    if || z || ∞ ≤ zTx
//        ν = || w || 1
//        k = 0
//    else
//        x = ej, 其中下标j满足 |zj| = ||z||∞
//        k = 1
//    end
//end

std::vector<double> sign(std::vector<double>& w) {
    std::vector<double> v(w.size());

    for (size_t i = 0; i < w.size(); i++) {
        if (w[i] > 0) {
            v[i] = 1.0;  // 正数
        }
        else if (w[i] < 0) {
            v[i] = -1.0; // 负数
        }
        else {
            v[i] = 0.0;  // 零
        }
    }

    return v;
}

double InnerProduct(std::vector<double>& a, std::vector<double>& b) {
    if (a.size() != b.size()) {
        // 向量长度不一致，无法计算内积
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
            //norm = element;
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
    int maxIndex = -1; // 初始化为-1，表示未找到最大值
    for (size_t i = 0; i < vec.size(); i++) {
        double absElement = std::abs(vec[i]);
        if (absElement > norm) {
            norm = absElement;
            maxIndex = static_cast<int>(i); // 更新最大值的下标
        }
    }
    std::vector<double> unitVector(n, 0.0);  // 初始化为0
    unitVector[maxIndex] = 1.0;  // 设置指定下标分量为1
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
        // 计算 w = Bx
        w = gauss_equation_solving(AT, x);
        // 计算 v = sign(w)
        v = sign(w);
        // 计算 z = BTv
        z = gauss_equation_solving(A, v);
        /*PrintVector(w);
        PrintVector(v);
        PrintVector(z);*/
        // 计算 ||z||∞ 和 zTx
        if (VectorInfinityNorm(z) <= InnerProduct(z, x)) {
            /*OneNorm = VectorOneNorm(w);
            std::cout << OneNorm << std::endl;*/
            k = 0;
        }
        else {
            // 更新 x 为 ej
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
    std::vector<double> result(numRows, 0.0); // 初始化结果向量为0
    if (numCols != b.size()) {
        // 如果矩阵的列数不等于向量的大小，无法相乘
        return result; // 返回全零向量
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

//算法3.2.1（计算Householder变换）
//function:[v,beta]=house(x)
//	n=length(x)(向量x的长度)
//	eta = ||x||∞；x=x/eta
//	sigma = x(2:n)T ·  x(2:n)
//	v(2:n)=x(2:n)
//	if sigma =0
//		beta=0
//	else
//		alpha=(x(1)^2+sigma)^0.5
//		if x(1)≤0
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

//算法3.3.1（计算QR分解：Householder方法）
//for j=1:n
//	if j<m
//		[v,beta] = house(A(j:m,j))
//		A(j:m,j:n) = (I_{m-j+1} - beta*v*vT) * A(j:m,j:n) //这是矩阵乘法
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

//正交变换法的基本步骤为：
//(1)计算A的QR分解；
//(2)计算c1 = QTb；
//(3)求解上三角方程组Rx = c1。
//已有A=QR=H1···HnR, 那么c1= Hn···H1 b, 
//其中Hk = diag(I_{k-1}, I_{m-k+1} - beta*v_k·v_k^T)为m*m维矩阵。
//递归调用MatrixVectorMultiply函数即可计算得到c1,最后用回代法back_subs计算得到x。

std::vector<std::vector<double>> HouseholderMatrix(std::vector<std::vector<double>>& A, std::vector<double>& d, int k) {
    int m = A.size();
    int n = A[0].size();
    std::vector<std::vector<double>> matrix(m, std::vector<double>(m, 0.0));
    for (int i = 0; i < k - 1; ++i) {
        matrix[i][i] = 1.0;
    }
    //PrintMatrix(matrix);
    std::vector<double> v = { 1.0 };
    // 将第 k 列的后 m-k 个分量连接到向量 v 后面
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
    // QR 分解
    std::vector<double> d(n, 0.0);
    QRDecomposition(A, d);
    /*PrintMatrix(A);
    PrintVector(d);*/
    // 计算 c1 = Q^T * b
    std::vector<double> c1 = b;
    //PrintVector(c1);
    std::vector<std::vector<double>> H(m, std::vector<double>(m, 0.0));
    for (int k = 1; k <= n; k++) {
        H = HouseholderMatrix(A, d, k);
        c1 = MatrixVectorMultiply(H, c1);
        /*cout << "第" << k << "次Householder变换：" << endl;
        PrintMatrix(H);
        PrintVector(c1);*/
    }
    // 解上三角方程组 Rx = c1
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
    //cout << "矩阵A为：" << endl;
    //PrintMatrix(A);
    //cout << "向量b为：" << endl;
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
    // 初始化随机数生成器
    std::random_device rd;  // 随机种子
    std::mt19937 gen(rd());  // Mersenne Twister 生成器
    std::uniform_real_distribution<double> dis(0.0, 1.0);  // 生成 [0.0, 1.0) 范围内的随机 double 数
    // 生成随机 double 数并填充向量
    for (int i = 0; i < N; ++i) {
        x[i] = dis(gen);
    }
    b = MatrixVectorMultiply(A, x);
    /*cout << "矩阵A为：" << endl;
    PrintMatrix(A);*/
    cout << "向量x为：" << endl;
    PrintVector(x);
    /*cout << "向量b为：" << endl;
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
    /*cout << "矩阵A为：" << endl;
    PrintMatrix(A);
    cout << "向量b为：" << endl;
    PrintVector(b);*/
}

vector<vector<double>> MultiplyMatrices(vector<vector<double>>& matrix1, vector<vector<double>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();
    // 检查两个矩阵是否可以相乘
    if (cols1 != rows2) {
        cerr << "Error: Matrix dimensions are not compatible for multiplication." << endl;
        return vector<vector<double>>();
    }
    // 初始化结果矩阵
    vector<vector<double>> result(rows1, vector<double>(cols2, 0.0));
    // 进行矩阵相乘
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
            row.erase(row.begin() + n, row.end()); // 删除从第 n+1 列到最后一列的元素
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
    vector<double> x(n, 0.0); // 初始化解向量x
    vector<double> x_new(n, 0.0); // 用于存储新的迭代结果
    vector<double> diff(n, 0.0);  // 用于存储计算解和标准解的误差
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
            //cout << "在 " << k + 1 << " 次迭代后收敛." << endl;
            cout << "迭代次数为：" << k + 1 << endl;
            return x_new;
        }
        x = x_new;
    }
    cout << "在指定的迭代次数内未收敛." << endl;
    return x;
}

vector<double> GS_Iteration(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0.0); // 初始化解向量x
    vector<double> x_new(n, 0.0); // 用于存储新的迭代结果
    vector<double> diff(n, 0.0);  // 用于存储计算解和标准解的误差
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; // 用于计算前半部分的求和
            double sum2 = 0.0; // 用于计算后半部分的求和
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
            //cout << "在 " << iter + 1 << " 次迭代后收敛." << endl;
            cout << "迭代次数为：" << iter + 1 << endl;
            return x_new;
        }
        x = x_new;
    }
    cout << "在指定的迭代次数内未收敛." << endl;
    return x;
}

vector<double> SOR_Iteration(vector<vector<double>>& A, vector<double>& b, double omega) {
    int n = A.size();
    vector<double> x(n, 0.0); // 初始化解向量x
    vector<double> x_new(n, 0.0); // 用于存储新的迭代结果
    vector<double> diff(n, 0.0);  // 用于存储计算解和标准解的误差
    //cout << "omega=" << omega << endl;
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; // 用于计算前半部分的求和
            double sum2 = 0.0; // 用于计算后半部分的求和
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
            cout << "在 " << iter + 1 << " 次迭代后收敛." << endl;
            //cout << iter + 1 ;
            return x_new;
        }
        x = x_new;
    }
    cout << "在指定的迭代次数内未收敛." << endl;
    return x;
}

int SOR_Performance(vector<vector<double>>& A, vector<double>& b, double omega) {
    int n = A.size();
    vector<double> x(n, 0.0); // 初始化解向量x
    vector<double> x_new(n, 0.0); // 用于存储新的迭代结果
    vector<double> diff(n, 0.0);  // 用于存储计算解和标准解的误差
    //cout << "omega=" << omega << endl;
    const int maxIterations = 20000;
    const double tolerance = 1e-6;
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; // 用于计算前半部分的求和
            double sum2 = 0.0; // 用于计算后半部分的求和
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
    double lowerBound = 1.0; // 区间下界
    double upperBound = 2.0; // 区间上界
    const int maxIterations = 20000; // 最大迭代次数
    const double tolerance = 1e-2; // 收敛阈值
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

    cout << "寻找最佳omega中..." << endl;
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
    cout << "最佳omega为：" << bestOmega << "\n寻找用时:" << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
    return bestOmega;
}

void Iterations(double epsilon) {
    int n = 99; //将[0,1]区间100等分
    double a = 0.5;
    double h = 1.0 / (n + 1);
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    //初始化A和b
    for (int i = 0; i < n - 1; i++)
    {
        A[i][i] = -2 * epsilon - h;
        A[i + 1][i] = epsilon;
        A[i][i + 1] = epsilon + h;
        b[i] = a * h * h;
    }
    A[n - 1][n - 1] = -2 * epsilon - h;
    b[n - 1] = a * h * h - epsilon - h;

    clock_t start_time;  // 记录开始时间

    cout << "epsilon = " << epsilon << endl;
    vector<double> accurate_y(n, 0.0);
    double temp = (1 - a) / (1 - exp(-1.0 / epsilon));
    for (int i = 0; i < n; i++)
    {
        accurate_y[i] = temp * (1 - exp(-(i + 1) * h / epsilon)) + a * (i + 1) * h;
    }
    cout << "精确解为：" << endl;
    PrintVector(accurate_y);
    cout << endl;

    start_time = clock();
    cout << "Jacobi迭代：" << endl;
    vector<double> y1 = Jacobi_Iteration(A, b);
    PrintVector(y1);
    cout << "Jacobi迭代法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;

    start_time = clock();
    cout << "G-S迭代：" << endl;
    vector<double> y2 = GS_Iteration(A, b);
    PrintVector(y2);
    cout << "G-S迭代法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;

    cout << "SOR迭代：" << endl;
    double omega = BisearchOmega(A, b);
    start_time = clock();
    vector<double> y3 = SOR_Iteration(A, b, omega);
    PrintVector(y3);
    cout << "SOR迭代法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
}

vector<vector<double>> MatrixSubtraction(vector<vector<double>> x, vector<vector<double>> y) {
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

    clock_t start_time;  // 记录开始时间

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
            cout << "Jacobi迭代在 " << k + 1 << " 次迭代后收敛." << endl;
            double mini = 1.0;
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (u1[i][j] < mini) {
                        mini = u1[i][j];
                    }
                }
            }
            //PrintMatrix(u1);
            cout << "解的最小分量为 " << mini << "." << endl;
            cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
            break;
        }
        u = u1;
        if (k == maxIterations - 1) {
            cout << "在指定的迭代次数内未收敛." << endl;
        }
    }
}

void GS_Iteration2(vector<vector<double>>& u) {

    clock_t start_time;  // 记录开始时间

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
            cout << "G-S迭代在 " << k + 1 << " 次迭代后收敛." << endl;
            double mini = 1.0;
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (u1[i][j] < mini) {
                        mini = u1[i][j];
                    }
                }
            }
            //PrintMatrix(u1);
            cout << "解的最小分量为 " << mini << "." << endl;
            cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
            break;
        }
        u = u1;
        if (k == maxIterations - 1) {
            cout << "在指定的迭代次数内未收敛." << endl;
        }
    }
}

void SOR_Iteration2(vector<vector<double>>& u, double omega) {
    clock_t start_time;  // 记录开始时间

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
            cout << "SOR迭代在 " << k + 1 << " 次迭代后收敛." << endl;
            double mini = 1.0;
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (u1[i][j] < mini) {
                        mini = u1[i][j];
                    }
                }
            }
            //PrintMatrix(u1);
            cout << "解的最小分量为 " << mini << "." << endl;
            cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
            break;
        }
        u = u1;
        if (k == maxIterations - 1) {
            cout << "在指定的迭代次数内未收敛." << endl;
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
    double lowerBound = 1.0; // 区间下界
    double upperBound = 2.0; // 区间上界
    const int maxIterations = 20000; // 最大迭代次数
    const double tolerance = 1e-2; // 收敛阈值
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

    cout << "寻找最佳omega中..." << endl;
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
    cout << "最佳omega为：" << bestOmega << "\n寻找用时:" << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
    return bestOmega;
}

void Iterations2(int n) {

    vector<vector<double>> u(n + 1, vector<double>(n + 1));
    //初始化A和b
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

//算法5.2.1（解对称正定方程组：共轭梯度法）
//x0 = 初值
//r0 = b - Axo; k = 0
//while rk ≠ 0
//k = k + 1
//if k = 1
//p0 = r0
//else
//Beta_{ k - 2 } = r ^ T_{ k - 1 }r_{ k - 1 } / r ^ T_{ k - 2 }r_{ k - 1 }
//p_{ k - 1 } = r_{ k - 1 } + Beta_{ k - 2 }p{ k - 2 }
//end
//alpha_{ k - 1 } = r ^ T_{ k - 1 }r_{ k - 1 } / p ^ T_{ k - 1 }·A·p_{ k - 1 }
//xk = x_{ k - 1 } + alpha_{ k - 1 }p_{ k - 1 }
//rk = r_{ k - 1 } - alpha_{ k - 1 }·A·p_{ k - 1 }
//end
//x = xk
//算法5.3.1（解对称正定方程组：实用共轭梯度法）
//x = 初值
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

    int sizeA = A.size(); //矩阵A为(n-1)*(n-1)维的
    int sizeu = int(sqrt(sizeA)); // 这里的sizeu是指解矩阵中间的矩阵维数，为n-1
    double h = 1.0 / (sizeu + 1);
    /*cout << sizeu << endl;*/
    vector<vector<double>> u(sizeu + 2, vector<double>(sizeu + 2, 0.0)); //解u为(n+1)*(n+1)维的
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
        //while (std::sqrt(rou) > epsilon * sqrt(dotProduct(b, b)) && k < k_max) {
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
        //cout << "alpha * ||p||∞ = " << VectorInfinityNorm(p) * alpha << endl;
    }
    cout << "共轭梯度法\n迭代次数为：" << k << "\n解向量为：\n";
    //PrintVector(x);
    for (int i = 0; i < sizeu;i++) {
        for (int j = 0; j < sizeu; j++) {
            u[i + 1][j + 1] = x[i * sizeu + j];
            /*b[i * sizeB + j] = sin((i + 1) * (j + 1) * h * h);
            b[i * sizeB + j] *= h * h / 4;*/
        }
    }
    PrintMatrix(u);
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
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
    while (std::sqrt(rou) > epsilon * sqrt(dotProduct(b, b)) && k < k_max) {
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
    cout << "共轭梯度法\n迭代次数为：" << k << "\n解向量为：\n";
    PrintVector(x);
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}

std::vector<std::vector<double>> generateMatrixA(int n) {
    double h = 1.0 / n;
    double diagValue = 1.0 + h * h / 4.0;
    double subDiagValue = -1.0 / 4.0;

    // Create the matrix A with the described properties
    std::vector<std::vector<double>> A(n - 1, std::vector<double>(n - 1, 0.0));

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

    int maxIterations = 1000; // 最大迭代次数
    double epsilon = 1e-6; // 容差
    int n = matrix.size();
    vector<double> x(n); 
    for (int i = 0; i < n; i++) {
        x[i] = 1.0 / (i + 1);
    } // 初始向量x, 防止x刚好是特征向量
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
            x[i] /= lambda; // 归一化向量
        }
        diff = VectorSubtraction(x, prevX);
        error = VectorInfinityNorm(diff);
        iteration++;
    }

    cout << "迭代次数为：" << iteration << endl;
    //PrintVector(x);
    //cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    return lambda;
}

//算法6.4.1（计算上Hessenberg分解：Householder变换法）
//for k = 1 : n - 2
//    [v, beta] = house(A(k + 1 : n, k))
//    A(k + 1: n , k : n) = (I - beta·vv^T)A(k + 1 : n, k : n)
//    A(1 : n, k + 1 : n) = A(1 : n, k + 1 : n)(I - beta·vv ^ T)
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
    /*cout << "原矩阵为：" << endl;
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
        // 应用 Householder 变换得到向量 v 和 beta
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
    /*cout << "上Hessenberg分解得到矩阵:" << endl;
    PrintMatrix(A);*/
}

//算法6.4.2（双重步位移的QR迭代）
//m = n - 1
//s = H(m, m) + H(n, n)
//t = H(m, m)H(n, n) - H(m, n)H(n, m)
//x = H(1, 1)H(1, 1) + H(1, 2)H(2, 1) - sH(1, 1) + t
//y = H(2, 1)(H(1, 1) + H(2, 2) - s)
//z = H(2, 1)H(3, 2)
//for k = 0 : n - 3
//    [v, beta] = house([x, y, z]T）
//    q = max{1,k}
//    H(k + 1:k + 3, q : n) = (I - beta·vvT)H(k + 1:k + 3, q : n)
//    r = min{ k + 4, n }
//    H(1:r, k + 1 : k + 3) = H(1:r, k + 1 : k + 3)(I - beta·vvT)
//    x = H(k + 2, k + 1)
//    y = H(k + 3, k + 1)
//    if k < n - 3
//        z = H(k + 4, k + 1)
//    end
//end
//[v, beta] = house([x, y]T)
//H(n - 1:n, n - 2 : n) = (I - beta·vvT)H(n - 1:n, n - 2 : n)
//H(1 : n, n - 1 : n) = H(1 : n, n - 1 : n)(I - beta·vvT)

void doubleShiftQR(vector<vector<double>>& H) {
    int n = H.size();
    if (n == 2) {
        double a = H[0][0];
        double b = H[0][1];
        double c = H[1][0];
        double d = H[1][1];

        // 计算旋转角度
        double s = hypot(a, c);
        double cos_val = a / s;
        double sin_val = c / s;
        // 应用 Givens 变换到原始矩阵
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
        /*cout << "双重步位移的QR迭代得到的矩阵为：" << endl;
        PrintMatrix(H);*/
    }
    

}

//调试双重步位移的QR迭代
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
//    cout << "双重步位移的QR迭代得到的矩阵为：" << endl;
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
        cout << "矩阵特征值为实数." << endl;
    }
    else {
        cout << "矩阵特征值为复数." << endl;
    }*/
    return (discriminant >= 0); // 如果判别式大于等于0，特征值为实数
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
        //cout << "矩阵块是准上三角阵对角元." << endl;
        return true;
    }
    else {
        //cout << "矩阵块不是准上三角阵对角元." << endl;
        return false;
    }
}

/*int quasi(vector<vector<double>>& H) {  //默认输入的矩阵是上Hessenberg矩阵
    int n = H.size();
    //int hess = 0;
    vector<int> hess = { }; // 存储满足条件的H33维数
    vector<vector<double>> submatrix_up;
    vector<vector<double>> submatrix_down;
    for (int m = 0; m <= n; m++) {
        if (m == 0) {
            cout << "\nm = " << m << endl;
            submatrix_up = getSubMatrix(H, n - m - 2, n - m - 1, n - m - 2, n - m - 1);
            cout << "submatrix_up = " << endl;
            PrintMatrix(submatrix_up);
            //cout << "考察位置的值为：" << H[n - m][n - m - 1] << endl;
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
            cout << "考察位置的值为：" << H[n - m][n - m - 1] << endl;
            if (!isQuasi(submatrix_up) && H[n - m][n - m - 1] == 0) {
                hess.push_back(m);
                //return m;
            }
        }
        // 对于1<m<n-1的情况，需要考虑上下两个子矩阵
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
            cout << "考察位置的值为：" << H[n - m][n - m - 1] << endl;


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
            cout << "考察位置的值为：" << H[n - m][n - m - 1] << endl;
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
            //cout << "考察位置的值为：" << H[n - m][n - m - 1] << endl;
            if (isQuasi(submatrix_down) && H[n - m][n - m - 1] == 0) {
                hess.push_back(m);
                //return m + 1;
            }
        }

    }
    std::cout << "满足条件的H33维数: ";
    for (int i = 0; i < hess.size(); ++i) {
        std::cout << hess[i] << " "; // 逐个输出vector元素
    }
    std::cout << std::endl;
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
            //cout << "k = " << k << "时是拟上三角阵" << endl;
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
        //cout << "特征值为：" << lambda1 << " and " << lambda2 << endl;
        cout << lambda1 << "\n" << lambda2 << endl;
    }
    else if (discriminant == 0) {
        double lambda = trace / 2;
        //cout << "特征值为：" << lambda << endl;
        cout << lambda << endl;
    }
    else {
        double realPart = trace / 2;
        double imaginaryPart = sqrt(-discriminant) / 2;
        //cout << "特征值为：" << realPart << " ± " << imaginaryPart << "i" << endl;
        cout << realPart << " ± " << imaginaryPart << "i" << endl;
    }
}

void prEigens(vector<vector<double>>& A) {
    int n = A.size();
    int k = 0;
    vector<vector<double>> B;
    //cout << "求解得到的特征值为：" << endl;
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
    /*cout << "原矩阵为：" << endl;
    PrintMatrix(A);*/
    int n = A.size();
    //double u = 1e-6; // 第二小题用
    double u = 1e-4; // 第三小题用
    hessenberg(A); //将矩阵A上Hessenberg化
    int m = quasi(A); // H33的维数
    int l = IrredHessenberg(A, m); // H22的维数
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
    //cout << "\n迭代次数为：" << k << endl;
    /*cout << "隐式QR迭代得到矩阵为：" << endl;
    PrintMatrix(A);*/
    /*cout << "\n方程组的解为：" << endl;
    prEigens(A);*/
}

