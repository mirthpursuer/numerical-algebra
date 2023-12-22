#include "exercise.h"
#include <random>
clock_t start, ende;

void exercise1_1()
{
    int N = 50; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);

    //初始化A和b
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
    //PrintMatrix(A);
    //PrintVector(b);

    clock_t start_time;  // 记录开始时间

    //不选主元
    //读取 A,b
    std::vector<std::vector<double>> A1 = A;
    vector<double> b1 = b;

    start_time = clock();  // 记录开始时间

    gauss_elim(A1);
    //PrintMatrix(A1);
    forward_subs1(A1, b1);
    back_subs(A1, b1);
    cout << "不选主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b1[i] << endl;
    //   cout << b1[i] << endl;
    cout << "不选主元的Gauss消去法解与精确解的差值：" << endl;
    for (int i = 0; i < N; i++)
        cout << "|x[" << i << "] -1| = " << abs(b1[i] - 1) << endl;
    //   cout << abs(b1[i] - 1) << endl;
    cout << "不选主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    // 全主元
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;
    std::vector<int> u(N);
    std::vector<int> v(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_full_pivoting(A2, u, v);
    //PrintMatrix(A2);
    forward_subs1(A2, b2);
    back_subs(A2, b2);
    cout << "全主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        //cout << b2[i] << endl;
        cout << "x[" << i << "] = " << b2[i] << endl;
    cout << "全主元的Gauss消去法解与精确解的差值：" << endl;
    for (int i = 0; i < N; i++)
        //cout << abs(b2[i] - 1) << endl;
        cout << "|x[" << i << "] -1| = " << abs(b2[i] - 1) << endl;

    cout << "全主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    //列主元
    std::vector<std::vector<double>> A3 = A;
    std::vector<double> b3 = b;
    std::vector<int> P(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_col_pivoting(A3, P);
    //PrintMatrix(A3);
    vector_pb(P, b3);
    forward_subs1(A3, b3);
    back_subs(A3, b3);
    //gauss_equation_solving(A3, b3);
    cout << "列主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        //cout << b3[i] << endl;
        cout << "x[" << i << "] = " << b3[i] << endl;
    cout << "列主元的Gauss消去法解与精确解的差值：" << endl;
    for (int i = 0; i < N; i++)
        //cout << abs(b3[i] - 1) << endl;
        cout << "|x[" << i << "] -1| = " << abs(b3[i] - 1) << endl;

    cout << "列主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}

void exercise1_2_1()
{
    //   初始化A和b
    int N = 100; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);

    clock_t start_time;  // 记录开始时间

    //初始化A和b
    for (int i = 0; i < N - 1; i++)
    {
        A[i][i] = 10;
        A[i + 1][i] = 1;
        A[i][i + 1] = 1;
    }
    A[N - 1][N - 1] = 1;
    //PrintMatrix(A);

    // 初始化随机数生成器
    std::random_device rd;  // 随机种子
    std::mt19937 gen(rd());  // Mersenne Twister 生成器
    std::uniform_real_distribution<double> dis(0.0, 1.0);  // 生成 [0.0, 1.0) 范围内的随机 double 数

    // 生成随机 double 数并填充向量
    for (int i = 0; i < N; ++i) {
        b[i] = dis(gen);
    }

    // 输出向量的内容
 /*   std::cout << "Randomly generated vector:" << std::endl;
    for (const auto& elem : b) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;*/

    //备份A和b，为后续改进平方根法所用
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;

    start_time = clock();  // 记录开始时间

    cholesky_decomp(A);
    //PrintMatrix(A);

    vector<vector<double>> AT(N, vector<double>(N));
    transposeMatrix(A, AT);
    //PrintMatrix(AT);

    forward_subs1(A, b);
    back_subs(AT, b);

    cout << "平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b[i] << endl;
    //cout <<  b[i] << endl;
    cout << "平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    start_time = clock();  // 记录开始时间

    modified_cholesky_decomp(A2);
    //PrintMatrix(A2);

    vector<vector<double>> A2T(N, vector<double>(N));
    transposeMatrix(A2, A2T);
    //PrintMatrix(A2T);

    forward_subs1(A2, b2);
    back_subs(A2T, b2);

    cout << "改进平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b2[i] << endl;
    //cout << b2[i] << endl;
    cout << "改进平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

}

void exercise1_2_2()
{
    int N = 10; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);
    Hilbert_Matrix(A);
    //PrintMatrix(A);
    for (int i = 0; i < N; ++i) {
        for (int j = 1; j <= N; ++j) {
            b[i] += 1.0 / (i + j);
        }
    }
    //备份A和b，为后续改进平方根法所用
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;

    clock_t start_time;  // 记录开始时间
    start_time = clock();  // 记录开始时间

    cholesky_decomp(A);
    //PrintMatrix(A);

    vector<vector<double>> AT(N, vector<double>(N));
    transposeMatrix(A, AT);
    //PrintMatrix(AT);

    forward_subs1(A, b);
    back_subs(AT, b);

    cout << "平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b[i] << endl;
    cout << "平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    start_time = clock();  // 记录开始时间

    modified_cholesky_decomp(A2);
    //PrintMatrix(A2);

    vector<vector<double>> A2T(N, vector<double>(N));
    transposeMatrix(A2, A2T);
    //PrintMatrix(A2T);

    forward_subs1(A2, b2);
    back_subs(A2T, b2);

    cout << "改进平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b2[i] << endl;
    cout << "改进平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

}

void exercise1_3_1()
{
    //   初始化A和b
    int N = 100; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);
    vector<double> x(N);
    clock_t start_time;  // 记录开始时间

    //初始化A和b
    for (int i = 0; i < N - 1; i++)
    {
        A[i][i] = 10;
        A[i + 1][i] = 1;
        A[i][i + 1] = 1;
    }
    A[N - 1][N - 1] = 1;
    //PrintMatrix(A);

    // 初始化随机数生成器
    std::random_device rd;  // 随机种子
    std::mt19937 gen(rd());  // Mersenne Twister 生成器
    std::uniform_real_distribution<double> dis(0.0, 1.0);  // 生成 [0.0, 1.0) 范围内的随机 double 数

    // 生成随机 double 数并填充向量
    for (int i = 0; i < N; ++i) {
        x[i] = dis(gen);
    }

    // 输出向量的内容
    /*std::cout << "Randomly generated vector:" << std::endl;
    for (const auto& elem : x) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;*/
    b = MatrixVectorMultiply(A, x);

    std::vector<std::vector<double>> A1 = A;
    vector<double> b1 = b;

    start_time = clock();  // 记录开始时间

    gauss_elim(A1);
    //PrintMatrix(A1);
    forward_subs1(A1, b1);
    back_subs(A1, b1);
    vector<double> d1 = VectorSubtraction(b1, x);
    cout << "不选主元的Gauss求解误差：" << VectorInfinityNorm(d1) << endl;
    //for (int i = 0; i < N; i++)
    //    //   cout << "x[" << i << "] = " << b1[i] << endl;
    //    cout << b1[i] << endl;

    cout << "不选主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    // 全主元
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;
    std::vector<int> u(N);
    std::vector<int> v(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_full_pivoting(A2, u, v);
    //PrintMatrix(A2);
    forward_subs1(A2, b2);
    back_subs(A2, b2);
    vector<double> d2 = VectorSubtraction(b2, x);
    cout << "全主元的Gauss求解误差：" << VectorInfinityNorm(d2) << endl;
    /*cout << "全主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << b2[i] << endl;*/
        //cout << "x[" << i << "] = " << b2[i] << endl;

    cout << "全主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    //列主元
    std::vector<std::vector<double>> A3 = A;
    std::vector<double> b3 = b;
    std::vector<int> P(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_col_pivoting(A3, P);
    //PrintMatrix(A3);
    vector_pb(P, b3);
    forward_subs1(A3, b3);
    back_subs(A3, b3);
    vector<double> d3 = VectorSubtraction(b3, x);
    cout << "列主元的Gauss求解误差：" << VectorInfinityNorm(d3) << endl;
    //cout << "列主元的Gauss消去法解：" << endl;
    //for (int i = 0; i < N; i++)
    //    cout << b3[i] << endl;
    //cout << "x[" << i << "] = " << b3[i] << endl;

    cout << "列主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

}

void exercise1_3_2()
{
    int N = 40; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);
    vector<double> x(N, 1.0);
    Hilbert_Matrix(A);
    //PrintMatrix(A);
    for (int i = 0; i < N; ++i) {
        for (int j = 1; j <= N; ++j) {
            b[i] += 1.0 / (i + j);
        }
    }
    /*for (const auto& elem : b) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;*/
    clock_t start_time;  // 记录开始时间
    start_time = clock();  // 记录开始时间

    //不选主元
    //读取 A,b
    std::vector<std::vector<double>> A1 = A;
    vector<double> b1 = b;
    gauss_elim(A1);
    //PrintMatrix(A1);
    forward_subs1(A1, b1);
    back_subs(A1, b1);
    vector<double> d1 = VectorSubtraction(b1, x);
    cout << "不选主元的Gauss求解误差：" << VectorInfinityNorm(d1) << endl;
    //cout << "不选主元的Gauss消去法解：" << endl;
    //for (int i = 0; i < N; i++)
    //    //   cout << "x[" << i << "] = " << b1[i] << endl;
    //    cout << b1[i] << endl;
    cout << "不选主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    // 全主元
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;
    std::vector<int> u(N);
    std::vector<int> v(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_full_pivoting(A2, u, v);
    //PrintMatrix(A2);
    forward_subs1(A2, b2);
    back_subs(A2, b2);
    //cout << "全主元的Gauss消去法解：" << endl;
    //for (int i = 0; i < N; i++)
    //    //   cout << "x[" << i << "] = " << b2[i] << endl;
    //    cout << b2[i] << endl;
    vector<double> d2 = VectorSubtraction(b2, x);
    cout << "不选主元的Gauss求解误差：" << VectorInfinityNorm(d2) << endl;
    cout << "全主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    //列主元
    std::vector<std::vector<double>> A3 = A;
    std::vector<double> b3 = b;
    std::vector<int> P(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_col_pivoting(A3, P);
    //PrintMatrix(A3);
    vector_pb(P, b3);
    forward_subs1(A3, b3);
    back_subs(A3, b3);
    //cout << "列主元的Gauss消去法解：" << endl;
    //for (int i = 0; i < N; i++)
    //    //   cout << "x[" << i << "] = " << b3[i] << endl;
    //    cout << b3[i] << endl;
    vector<double> d3 = VectorSubtraction(b3, x);
    cout << "不选主元的Gauss求解误差：" << VectorInfinityNorm(d3) << endl;
    cout << "列主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}

void exercise2_1_1(int N)
{
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);
    Hilbert_Matrix(A);
    //cout << "Hilbert矩阵为：" << endl;
    //PrintMatrix(A);
    //std::cout << std::endl;
    vector<vector<double>> AT(N, vector<double>(N));
    transposeMatrix(A, AT);
    //cout << "Hilbert矩阵转置为：" << endl;
    //PrintMatrix(AT);
    //std::cout << std::endl;
    //std::vector<double> x(N, 1.0/N);
    //;PrintVector(x);
    std::cout << N << "阶Hilbert矩阵的无穷范数条件数为" << MatrixInfinityNorm(A) * MatrixOneNorm(N, A) << std::endl;
    std::cout << std::endl;
    /*  调试函数
      std::vector<double> vector1 = { 1.0, -2.0, 3.0, -4.0, 5.0 };
      PrintVector(vector1);
      std::vector<double> vector2 = { 1.0, 2.0, 3.0, 4.0, 5.0 };
      PrintVector(vector2);
      std::vector<double> vector3 = sign(vector1);
      PrintVector(vector3);
      std::cout << InnerProduct(vector1, vector2) << std::endl;
      std::cout << VectorOneNorm(vector3) << std::endl;
      std::vector<double> vector4 = UnitVectorGenerating(vector1, N);
      PrintVector(vector4);
      std::cout << MatrixInfinityNorm(A) << std::endl;*/

}

void exercise2_1_2(int N)
{
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> x(N);
    //PrintVector(x);
    vector<double> x1(N);
    vector<double> x2(N);
    vector<double> b(N);
    for (int i = 0; i <= N - 1; i++)
    {
        A[i][i] = 1;
        for (int j = 0; j < i; j++)
        {
            A[i][j] = -1;
        }
        A[i][N - 1] = 1;
    }
    //PrintMatrix(A);
    // 初始化随机数生成器
    std::random_device rd;  // 随机种子
    std::mt19937 gen(rd());  // Mersenne Twister 生成器
    std::uniform_real_distribution<double> dis(0.0, 1.0);  // 生成 [0.0, 1.0) 范围内的随机 double 数
    // 生成随机 double 数并填充向量
    for (int i = 0; i < N; ++i) {
        x[i] = dis(gen);
    }
    std::cout << "n=" << N << "" << std::endl;
    std::cout << "随机选取的x为：" << std::endl;
    PrintVector(x);
    b = MatrixVectorMultiply(A, x);
    //PrintVector(b);
    x1 = gauss_equation_solving(A, b);
    x2 = VectorSubtraction(x, x1);
    //PrintMatrix(A);
    std::cout << "计算得到的x为：" << std::endl;
    PrintVector(x1);
    double estimated_accuracy = MatrixInfinityNorm(A) * MatrixOneNorm(N, A) * VectorInfinityNorm(x2) / VectorInfinityNorm(b);
    double actual_accuracy = VectorInfinityNorm(x2) / VectorInfinityNorm(x);
    std::cout << "估计精度为：" << estimated_accuracy << "\n真实精度为：" << actual_accuracy << "\n" << std::endl;
}

void exercise3_1_1()
{
    //用于测试Householder变换

    /*std::vector<double> x = { 0.0, 1.0, 0.0 };
    std::vector<double> v;
    double beta;
    house(x, v, beta);
    std::cout << "x:" << std::endl;
    PrintVector(x);
    std::cout << "v:" << std::endl;
    PrintVector(v);
    std::cout << "beta:" << std::endl;
    cout << beta << endl;*/

    //用于测试QRdecomposition

    /*std::vector<std::vector<double>> A = {
        {0.0, 4.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 3.0, 2.0}
    };
    PrintMatrix(A);
    std::vector<double> d = { 0.0, 0.0};
    QRDecomposition(A, d);
    //PrintMatrix(A);
    std::cout << "R Matrix:" << std::endl;
    // 打印 R 矩阵，它是 A 的上三角部分
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            if (j < i) {
                std::cout << "0.0\t";
            }
            else {
                std::cout << A[i][j] << "\t";
            }
        }
        std::cout << std::endl;
    }
    PrintVector(d);*/

    clock_t start_time;  // 记录开始时间

    cout << "求解方程组1：" << endl;
    std::vector<int> matrix_sizes1 = { 10, 30, 50, 55, 56, 84 };
    for (int i = 0; i < matrix_sizes1.size(); i++) {
        int N1 = matrix_sizes1[i];
        cout << "矩阵规模为：" << N1 << endl;
        vector<vector<double>> A1(N1, vector<double>(N1));
        vector<double> b1(N1);
        vector<double> x1(N1);
        //vector<double> accurate_x(N1);
        equation_generating1(A1, b1);
        start_time = clock();  // 记录开始时间
        x1 = QR_equation_solving(A1, b1);
        /*cout << "解x为：" << endl;
        PrintVector(x1);*/
        for (int j = 0; j < x1.size(); j++) {
            x1[j] -= 1;
        }
        cout << "与准确解的差值：" << VectorInfinityNorm(x1) << endl;
        cout << "用QR分解求解方程组1用时: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
        std::cout << std::endl;
    }

    cout << "求解方程组2：\n" << endl;
    std::vector<int> matrix_sizes2 = { 10, 50, 100 };
    //std::vector<int> matrix_sizes2 = { 10  };
    for (int i = 0; i < matrix_sizes2.size(); i++) {
        int N2 = matrix_sizes2[i];
        cout << "矩阵规模为：" << N2 << endl;
        vector<vector<double>> A2(N2, vector<double>(N2));
        vector<double> actual_x2(N2);
        vector<double> b2(N2);
        vector<double> x2(N2);
        vector<double> difference(N2);
        equation_generating2(A2, actual_x2, b2);
        start_time = clock();  // 记录开始时间
        x2 = QR_equation_solving(A2, b2);
        cout << "解x为：" << endl;
        PrintVector(x2);
        difference = VectorSubtraction(actual_x2, x2);
        cout << "与准确解的差值：" << VectorInfinityNorm(difference) << endl;
        cout << "用QR分解求解方程组2用时: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
        std::cout << std::endl;
    }

    cout << "求解方程组3：\n" << endl;
    std::vector<int> matrix_sizes3 = { 10, 13, 20, 40 };
    for (int i = 0; i < matrix_sizes3.size(); i++) {
        int N3 = matrix_sizes3[i];
        cout << "矩阵规模为：" << N3 << endl;
        vector<vector<double>> A3(N3, vector<double>(N3));
        vector<double> b3(N3);
        vector<double> x3(N3);
        //vector<double> accurate_x(N3);
        equation_generating3(A3, b3);
        start_time = clock();  // 记录开始时间
        x3 = QR_equation_solving(A3, b3);
        cout << "解x为：" << endl;
        PrintVector(x3);
        for (int j = 0; j < x3.size(); j++) {
            x3[j] -= 1;
        }
        cout << "与准确解的差值：" << VectorInfinityNorm(x3) << endl;
        cout << "用QR分解求解方程组3用时: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
        std::cout << std::endl;
    }
    //int N3 = 10;
    //vector<vector<double>> A3(N3, vector<double>(N3));
    //vector<double> b3(N3);
    //vector<double> x3(N3);
    //equation_generating3(A3, b3);
    //start_time = clock();  // 记录开始时间
    //x3 = QR_equation_solving(A3, b3);
    //cout << "解x为：" << endl;
    //PrintVector(x3);
    //cout << "用QR分解求解方程组3用时: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    //std::cout << std::endl;
}

void exercise3_1_2()
{
    /*//调试ReduceMatrix函数
    int m = 4; // 行数
    int n = 2; // 保留的列数
    std::vector<std::vector<double>> Q(m, std::vector<double>(m, 0.0));
    ReduceMatrix(Q, n);
    for (const std::vector<double>& row : Q) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }*/
    clock_t start_time;  // 记录开始时间
    int m = 7;
    int n = 3;
    vector<double> t = { -1,-0.75,-0.5,0,0.25,0.5,0.75 };
    vector<double> y = { 1,0.8125,0.75,1,1.3125,1.75,2.3125 };
    vector<double> x(n, 0.0);
    vector<double> r(m, 0.0);
    vector<vector<double>> A(m, vector<double>(n, 0.0));
    vector<vector<double>> A1 = A;
    for (int i = 0; i < m; i++) {
        A[i][0] = t[i] * t[i];
        A[i][1] = t[i];
        A[i][2] = 1;
    }
    start_time = clock();
    x = LS_proplem_solving(A, y);
    cout << "拟合多项式为：y = " << x[0] << "*t^2 +" << x[1] << "*t +" << x[0] << " *1." << endl;
    r = MatrixVectorMultiply(A1, x);
    r = VectorSubtraction(r, y);
    cout << "残向量的二范数为：" << VectorTwoNorm(r) << endl;
    // 输出运行时间
    cout << "运行时间： " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}

void exercise3_1_3()
{
    vector<vector<double>> A =
    { {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
    {1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
    {1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
    {1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
    {1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
    {1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
    {1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
    {1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
    {1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
    {1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
    {1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
    {1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
    {1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
    {1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
    {1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
    {1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
    {1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
    {1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
    {1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
    {1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
    {1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
    {1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
    {1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
    {1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
    {1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
    {1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
    {1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
    {1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
    vector<double> b =
    { 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
    28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
    30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
    37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 };
    int m = A.size();
    int n = A[0].size();
    clock_t start_time;  // 记录开始时间
    start_time = clock();  // 记录开始时间
    vector<double> x(n, 0.0);
    vector<double> r(m, 0.0);
    vector<vector<double>> A1 = A;
    x = LS_proplem_solving(A, b);
    PrintVector(x);
    cout << "\n房屋估价的拟合模型为：\ny = " << x[0] << "+";
    for (int i = 1; i < n - 1; i++) {
        cout << x[i] << "a" << i << " + ";
        if (i == 5) {
            cout << "\n   ";
        }
    }
    cout << x[n - 1] << "a" << n - 1 << endl;
    r = MatrixVectorMultiply(A1, x);
    r = VectorSubtraction(r, b);
    cout << "残向量的二范数为：" << VectorTwoNorm(r) << endl;
    // 输出运行时间
    cout << "运行时间： " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}

void exercise4_1()
{
    vector<double> epsilons = { 1.0, 0.1, 0.01, 0.0001 };
    for (int i = 0; i < epsilons.size(); i++) {
        Iterations(epsilons[i]);
    }
}

void exercise4_2()
{
    vector<int> partitions = { 20, 40, 60 };
    for (int i = 0; i < partitions.size(); i++) {
        Iterations2(partitions[i]);
    }
}

void exercise5_1_1() {
    cout << "exercise5.1:" << endl;
    int n = 20;
    int sizeA = (n - 1) * (n - 1);
    int sizeB = n - 1;
    double h = 1.0 / n; 
    double diagonalValue = 1.0 + h * h / 4;
    double offDiagonalValue = -1.0 / 4.0;

    vector<vector<double>> A(sizeA, vector<double>(sizeA, 0.0));
    vector<vector<double>> B(sizeB, vector<double>(sizeB, 0.0));
    vector<vector<double>> S(sizeB, vector<double>(sizeB, 0.0));
    for (int i = 0; i < B.size(); ++i) {
        B[i][i] = 1;
    }
    //PrintMatrix(B);
    for (int i = 0; i < S.size(); ++i) {
        S[i][i] = diagonalValue;
        if (i > 0) {
            S[i][i - 1] = S[i - 1][i] = offDiagonalValue;
        }
    }
    //PrintMatrix(S);
    for (int i = 0; i <= sizeB-1; i++) {
        int startRow = i * sizeB;  // 选取起始行
        int startCol = i * sizeB;  // 选取起始列
        for (int j = 0; j < sizeB; ++j) {
            for (int k = 0; k < sizeB; ++k) {
                A[startRow + j][startCol + k] = S[j][k];
                if (i <= sizeB - 2) {
                    A[startRow + j + sizeB][startCol + k] = B[j][k];
                    A[startRow + j][startCol + k + sizeB] = B[j][k];
                }
                
            }
        }
    }
    //PrintMatrix(A);
    
    vector<double> x(sizeA);
    vector<double> b(sizeA);
    for (int i = 0; i < sizeB;i++) {
        for (int j = 0; j < sizeB; j++) {
            b[i * sizeB + j] = sin((i + 1) * (j + 1) * h * h);
            b[i * sizeB + j] *= h * h / 4;
        }
    }
    //PrintVector(b);

    ConjugateGradient1(A, x, b);
    cout << endl;
}

void exercise5_1_2() {
    cout << "exercise5.2:" << endl;
    vector<int> matrixsize = { 20, 40, 60, 80 };
    for (int i = 0; i < matrixsize.size(); i++) {
        int N = matrixsize[i];
        cout << "N = " << N << endl;
        vector<vector<double>> A(N, vector<double>(N));
        vector<double> x(N);
        vector<double> b(N);
        Hilbert_Matrix(A);
        for (int i = 0; i < N; ++i) {
            for (int j = 1; j <= N; ++j) {
                b[i] += (1.0 / 3) / (i + j);
            }
        }
        ConjugateGradient2(A, x, b);
        cout << endl;
    }
}

void exercise5_1_3() {
    cout << "exercise5.3:" << endl;
    clock_t start_time;
    
    vector<vector<double>> A = {
        {10, 1,  2,  3,  4},
        {1,  9, -1,  2, -3},
        {2, -1,  7,  3, -5},
        {3,  2,  3, 12, -1},
        {4, -3, -5, -1, 15}
    };
    vector<double> b = { 12, -27, 14, -17, 12 };
    vector<double> x(b.size());

    start_time = clock();
    cout << "Jacobi迭代法" << endl;
    x = Jacobi_Iteration(A, b);
    cout << "解向量为：";
    PrintVector(x);
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    start_time = clock();
    cout << "G-S迭代法" << endl;
    x = GS_Iteration(A, b);
    cout << "解向量为：";
    PrintVector(x);
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    
    std::fill(x.begin(), x.end(), 0.0); // 把x重新归零
    ConjugateGradient2(A, x, b);

}