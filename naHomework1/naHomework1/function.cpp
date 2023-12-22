#include "function.h"
#include <cmath>

const double EPSILON = 1e-10;  // 定义一个很小的常数，用于判断浮点数是否接近于0

void print(vector<vector<double>>& mat) {
    int n = mat.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
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
/*
void forward_subs(vector<vector<double>>& L, vector<double>& b)
{
    int n = L.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            b[i] -= L[i][j] * b[j];
        }
        b[i] /= L[i][i];
    }

}*/

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

//void gauss_elim(std::vector<std::vector<double>>& A) {
//    int n = A.size();
//
//    for (int k = 0; k < n - 1; k++) {
//        for (int i = k + 1; i < n; i++) {
//            A[i][k] = A[i][k] / A[k][k];
//            for (int j = k + 1; j < n; j++) {
//                A[i][j] = A[i][j] - A[i][k] * A[k][j];
//            }
//        }
//    }
//}

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
//
void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
{
    int n = A.size();
    u.resize(n);

    for (int i = 0; i < n; i++) {
        u[i] = i;
    }

    for (int k = 0; k < n - 1; k++) {
        int pivot_row = k;
        double max_abs = std::abs(A[k][k]);

        // Find pivot element
        for (int i = k + 1; i < n; i++) {
            double abs_val = std::abs(A[i][k]);
            if (abs_val > max_abs) {
                max_abs = abs_val;
                pivot_row = i;
            }
        }

        // Check for matrix singularity
        if (max_abs == 0) {
            std::cout << "Matrix is singular. Aborting." << std::endl;
            return;
        }

        // Swap rows in A and u
        for (int i = 0; i < n; i++) {
            swap(A[i][k], A[i][pivot_row]);
        }
        std::swap(u[k], u[pivot_row]);


        for (int i = k + 1; i < n; i++) {
            A[i][k] = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }
}

void vector_pb(vector<int>& u, vector<double>& b)
{
    // 计算向量P*b
    int n = u.size();
    vector<double> temp_b(n, 0);

    for (int i = 0; i < n; i++) {
        temp_b[i] = b[u[i]];
    }

    for (int i = 0; i < n; i++) {
        b[i] = temp_b[i];
    }
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

//
//void cholesky_decomp(vector<vector<double>>& A)
//{
//    int n = A.size();
//
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < i; j++) {
//            double sum = 0;
//            for (int k = 0; k < j; k++)
//                sum += A[i][k] * A[j][k];
//            A[i][j] = (A[i][j] - sum) / A[j][j];  // 对称正定阵标准Cholesky分解计算
//        }
//
//        double sum = A[i][i];
//        for (int k = 0; k < i; k++)
//            sum -= A[i][k] * A[i][k];
//
//        if (sum < 0) {
//            cerr << "矩阵非正定，无法进行Cholesky分解。" << endl;
//            return;
//        }
//
//        A[i][i] = sqrt(sum);
//    }
//}

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

void transposeMatrix(const std::vector<std::vector<double>>& inputMatrix, std::vector<std::vector<double>>& transposedMatrix) {
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


