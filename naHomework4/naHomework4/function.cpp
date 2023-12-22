#include "function.h"
#include <cmath>
#include <random>
#include <iomanip> 

const double EPSILON = 1e-10;  // 定义一个很小的常数，用于判断浮点数是否接近于0

void PrintMatrix(vector<vector<double>>& mat) {
    // 设置输出精度为四位小数
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
    // 设置输出精度为四位小数
    std::cout << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
        // 每十个分量换一次行
        if ((i + 1) % 11 == 0) {
            std::cout << std::endl;
        }
    }
    // 打印完所有分量后换行
    std::cout << std::endl;
    //std::cout << std::endl;
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

std::vector<std::vector<double>> MultiplyMatrices(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();
    // 检查两个矩阵是否可以相乘
    if (cols1 != rows2) {
        std::cerr << "Error: Matrix dimensions are not compatible for multiplication." << std::endl;
        return std::vector<std::vector<double>>();
    }
    // 初始化结果矩阵
    std::vector<std::vector<double>> result(rows1, std::vector<double>(cols2, 0.0));
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
            cout << "在 " << k + 1 << " 次迭代后收敛." << endl;
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
            cout << "在 " << iter + 1 << " 次迭代后收敛." << endl;
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
    cout << "最佳omega为：" << bestOmega  << "\n寻找用时:" << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
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