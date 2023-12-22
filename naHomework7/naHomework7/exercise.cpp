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
    for (int i = 0; i <= sizeB - 1; i++) {
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

void exercise6_1() {
    clock_t start_time;
    start_time = clock();
    vector<double> equation1 = { 1.0, -5.0, 3.0 };
    vector<double> equation2 = { 0.0, -3.0, -1.0 };
    vector<double> equation3 = { 101.0, 208.01, 10891.01, 9802.08, 79108.9, -99902.0, 790.0, -1000.0 };
    //PrintMatrix(equation_to_matrix(equation1));
    int m1 = equation1.size();
    vector<vector<double>> A1(m1, vector<double>(m1));
    A1 = equation_to_matrix(equation1);
    //PrintMatrix(A1);
    cout << "多项式方程为：x^3 + x^2 - 5x + 3 = 0" << endl;
    double lambda1 = powerMethod(A1);
    cout << "模最大根为：" << lambda1 << endl;
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
     
    int m2 = equation2.size();
    vector<vector<double>> A2(m2, vector<double>(m2));
    A2 = equation_to_matrix(equation2);
    //PrintMatrix(A2);
    cout << "多项式方程为：x^3 - 3x - 1 = 0" << endl;
    double lambda2 = powerMethod(A2);
    cout << "模最大根为：" << lambda2 << endl;
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;

    int m3 = equation3.size();
    vector<vector<double>> A3(m3, vector<double>(m3));
    A3 = equation_to_matrix(equation3);
    //PrintMatrix(A3);
    cout << "多项式方程为：x^8 + 101x^7 + 208.01x^6 + 10891.01x^5 + 9802.08x^4 + 79108.9x^3 - 99902x^2 + 790x - 1000 = 0" << endl;
    double lambda3 = powerMethod(A3);
    cout << "模最大根为：" << lambda3 << endl;
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
}

void exercise6_2_1() {
    // 测试上Hessenberg分解
    /*const int rows = 7;
    const int cols = 7;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(1.0, 10.0); // 修改范围以符合你的需求
    // 生成随机矩阵
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = dis(gen);
        }
    }
    hessenberg(matrix);
    doubleShiftQR(matrix);*/

    //测试双重步位移的QR迭代
    /*std::vector<std::vector<double>> A = {
        {1, 2, 3, 4, 5, 6, 7},
        {8, 9, 10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19, 20, 21},
        {22, 23, 24, 25, 26, 27, 28},
        {29, 30, 31, 32, 33, 34, 35},
        {36, 37, 38, 39, 40, 41, 42},
        {43, 44, 45, 46, 47, 48, 49}
    };
    PrintMatrix(A);
    doubleShiftQR(A);*/
    
    //测试判断特征值是否为实数
    /*vector<vector<double>> A = {
        {2, -2},
        {1, 9}
        //{6, 1},
        //{-8, 0}
    };
    bool realEigenvalues = areEigenvaluesReal(A);
    if (realEigenvalues) {
        cout << "矩阵特征值为实数." << endl;
    }
    else {
        cout << "矩阵特征值为复数." << endl;
    }*/

    //测试置零函数
    /*std::vector<std::vector<double>> A = {
        {1, 2, 3, 4, 5, 6, 7},
        {8, 9, 10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19, 20, 21},
        {22, 23, 24, 25, 26, 27, 28},
        {29, 30, 31, 32, 33, 34, 35},
        {36, 37, 38, 39, 40, 41, 42},
        {43, 44, 45, 46, 47, 48, 49}
    };
    PrintMatrix(A);
    zeroing(A, 100);
    PrintMatrix(A);*/

    //测试二阶矩阵块是否是拟上三角阵的对角元
    /*vector<vector<double>> A = {
        {2, -2},
        {1, 1}
    };
    cout << isQuasi(A) << endl;*/

    //测试寻找拟上三角阵
    /*vector<vector<double>> A1 = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {2, 3, 0, 0, 0, 0, 0, 0},
        {0, 0, 5, 0, 0, 0, 0, 0},
        {0, 0, 0, 6, 1, 0, 0, 0},
        {0, 0, 0,-8, 0,-2, 0, 0},
        {0, 0, 0, 0, 9, 2, 1, 0},
        {0, 0, 0, 0, 0, 0, 4, 1},
        {0, 0, 0, 0, 0, 0, 0, 5}
    };
    vector<vector<double>> A2 = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {0, 3, 4, 0, 0, 0, 0, 0},
        {0, 0, 5, 0, 0, 0, 0, 0},
        {0, 0, 0, 6, 7, 0, 0, 0},
        {0, 0, 0,-8, 0,-2, 0, 0},
        {0, 0, 0, 0, 0, 2, 1, 0},
        {0, 0, 0, 0, 0, 9, 4, 0},
        {0, 0, 0, 0, 0, 0, 1, 5}
    };
    vector<vector<double>> A3 = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {0, 3, 4, 0, 0, 0, 0, 0},
        {0, 0, 5, 0, 0, 0, 0, 0},
        {0, 0, 0, 6, 7, 0, 0, 0},
        {0, 0, 0,-8, 9,-2, 0, 0},
        {0, 0, 0, 0, 0, 2, 1, 0},
        {0, 0, 0, 0, 0,-9, 4, 0},
        {0, 0, 0, 0, 0, 0, 0, 5}
    };
    vector<vector<double>> A4 = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {0, 3, 4, 0, 0, 0, 0, 0},
        {0, 0, 5, 0, 0, 0, 0, 0},
        {0, 0, 0, 6, 7, 0, 0, 0},
        {0, 0, 0,-1, 0,-2, 0, 0},
        {0, 0, 0, 0, 0, 2, 1, 0},
        {0, 0, 0, 0, 0,-9, 4, 0},
        {0, 0, 0, 0, 0, 0, 0, 5}
    };

    cout << "原矩阵A1为：" << endl;
    PrintMatrix(A1);
    cout << "H33维度为" << quasi(A1) << endl;

    cout << "原矩阵A2为：" << endl;
    PrintMatrix(A2);
    cout << "H33维度为" << quasi(A2) << endl;

    cout << "原矩阵A3为：" << endl;
    PrintMatrix(A3);
    cout << "H33维度为" << quasi(A3) << endl;

    cout << "原矩阵A4为：" << endl;
    PrintMatrix(A4);
    cout << "H33维度为" << quasi(A4) << endl;*/
    
    // 测试判断是否是不可约上Hessenberg矩阵
    /*vector<vector<double>> A = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {0, 3, 4, 0, 0, 0, 0, 0},
        {0, 0, 5, 0, 0, 0, 0, 0},
        {0, 0, 2, 6, 7, 0, 0, 0},
        {0, 0, 0,-8, 9,-2, 0, 0}, 
        {0, 0, 0, 0, 0, 2, 1, 0},
        {0, 0, 0, 0, 0,-9, 4, 0},
        {0, 0, 0, 0, 0, 0, 0, 5}
    };
    cout << isHessenberg(A) << endl;
    vector<vector<double>> B = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {1, 3, 4, 0, 0, 0, 0, 0},
        {0, 3, 5, 0, 0, 0, 0, 0},
        {0, 0, 2, 6, 7, 0, 0, 0},
        {0, 0, 0,-8, 9,-2, 0, 0},
        {0, 0, 0, 0, 1, 2, 1, 0},
        {0, 0, 0, 0, 0,-9, 4, 0},
        {0, 0, 0, 0, 0, 0, -3, 5}
    };
    cout << isHessenberg(B) << endl;*/

    // 测试不可约Hessenberg矩阵
    /*vector<vector<double>> A = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {0, 3, 4, 0, 0, 0, 0, 0},
        {0, 0, 5, 0, 0, 0, 0, 0},
        {0, 0, 2, 6, 7, 0, 0, 0},
        {0, 0, 0,-8, 9,-2, 0, 0},
        {0, 0, 0, 0, 0, 2, 1, 0},
        {0, 0, 0, 0, 0,-9, 4, 0},
        {0, 0, 0, 0, 0, 0, 0, 5}
    };
    cout << IrredHessenberg(A, 3) << endl;
    vector<vector<double>> B = {
        {1, 2, 0, 0, 0, 0, 0, 0},
        {2, 3, 4, 0, 0, 0, 0, 0},
        {0, 5, 5, 0, 0, 0, 0, 0},
        {0, 0, 2, 6, 7, 0, 0, 0},
        {0, 0, 0,-8, 9,-2, 0, 0},
        {0, 0, 0, 0, 9, 2, 1, 0},
        {0, 0, 0, 0, 0,-9, 4, 0},
        {0, 0, 0, 0, 0, 0, 0, 5}
    };
    cout << IrredHessenberg(B, 3) << endl;
    vector<vector<double>> C;
    C = getHessenberg(A, 3);
    PrintMatrix(C);
    C = getHessenberg(B, 3);
    PrintMatrix(C);
    setHessenberg(A, C, 3);
    PrintMatrix(A);*/

    // 测试求解二阶矩阵的复特征值
    /*vector<vector<double>> A = {
        {-1.4625, -1.0557},
        { 1.0675, -1.3061}
    };*/
    //eigenvalues2D(A);

    // 测试隐式QR迭代
    /*vector<vector<double>> A = {
        {1.1908, -1.0565, -2.1707, 0.5913, 0.0000, 0.7310},
        {-1.2025, 1.4151, -0.0592, -0.6436, -0.3179, 0.5779},
        {-0.0198, -0.8051, -1.0106, 0.3803, 1.0950, 0.0403},
        {-0.1567, 0.5287, 0.6145, -1.0091, -1.8740, 0.6771},
        {-1.6041, 0.2193, 0.5077, -0.0195, 0.4282, 0.5689},
        {0.2573, -0.9219, 1.6924, -0.0482, 0.8956, -0.2556}
    };
    cout << "原矩阵为：" << endl;
    PrintMatrix(A);
    cout << "求解得到的特征值为：" << endl;
    implicitQR(A);*/
}

void exercise6_2_2() {
    // 测试求解方程程序
    //cout << "多项式方程为：(x-1)(x-2)(x-3)(x-4) = x^4-10x^3+35x^2-50x+24=0" << endl;
    //vector<double> equation = { -10, 35, -50, 24 };
    //int m = equation.size();
    //vector<vector<double>> A(m, vector<double>(m));
    //A = equation_to_matrix(equation);
    ////PrintMatrix(A);
    //implicitQR(A);
    clock_t start_time;
    start_time = clock();
    cout << "多项式方程为：x^41 + x^3 + 1 = 0" << endl;
    vector<double> equation = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 1.0
    };
    int m = equation.size();
    vector<vector<double>> A(m, vector<double>(m));
    A = equation_to_matrix(equation);
    //PrintMatrix(A);
    cout << "迭代次数:"<<implicitQR(A) << endl;
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    cout << "方程组的解为：" << endl;
    prEigens(A);
}

void exercise6_2_3() {
    /*//测试求解特征值
    vector<vector<double>> A = {
        { 9.1, 3.0, 2.6, 4.0 },
        { 4.2, 5.3, 4.7, 1.6 },
        { 3.2, 1.7, 9.4, x[0]},
        { 6.1, 4.9, 3.5, 6.2 }
    };
    PrintMatrix(A);
    cout << "迭代次数:" << implicitQR(A) << endl;
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
    PrintMatrix(A);
    cout << "矩阵的特征值为：" << endl;
    prEigens(A);*/
    clock_t start_time;
    start_time = clock();
    vector<double> x = { 0.9, 1.0, 1.1 };
    for (int i = 0; i < x.size(); i++) {
        cout << "x = " << x[i] << endl;
        vector<vector<double>> A = {
        { 9.1, 3.0, 2.6, 4.0 },
        { 4.2, 5.3, 4.7, 1.6 },
        { 3.2, 1.7, 9.4, x[i]},
        { 6.1, 4.9, 3.5, 6.2 }
        };
        //PrintMatrix(A);
        cout << "迭代次数:" << implicitQR(A) << endl;
        cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
        //PrintMatrix(A);
        cout << "矩阵的特征值为：" << endl;
        prEigens(A);
        cout << endl;
    }
}

void exercise7_1_1() {
    
    vector<vector<double>> A = {
        {4, 1,-2, 2},
        {1, 2, 0, 1},
        {-2,0, 3,-2},
        {2, 1,-2,-1}
    };
    PrintMatrix(A);
    //测试非对角范数和经典Jacobi方法
    /*cout << "A的非对角范数为：" << offDiagNorm(A) << endl;
    vector<vector<double>> J;
    J = jacobiClassic(A, 0, 2);
    PrintMatrix(A);
    PrintMatrix(J);*/
    //测试过关判断函数
    /*cout << passingThreshold(A, 1) << endl;
    cout << passingThreshold(A, 3) << endl;
    cout << passingThreshold(A, offDiagNorm(A)) << endl;*/
    //测试过关Jacobi方法
    vector<vector<double>> Q = thresholdJacobi(A);
    PrintMatrix(A);
    PrintMatrix(Q);
}

void exercise7_1_2(){
    clock_t start_time;
    vector<int> dimension = { 50, 60, 70, 80, 90, 100 };
    //vector<int> dimension = { 50 };
    for (int i = 0; i < dimension.size(); i++) {
        start_time = clock();
        int n = dimension[i];
        vector<vector<double>> A(n, vector<double>(n));
        for (int i = 0; i < n - 1; i++)
        {
            A[i][i] = 4;
            A[i + 1][i] = 1;
            A[i][i + 1] = 1;
        }
        A[n - 1][n - 1] = 4;
        //PrintMatrix(A);
        cout << "n = " << n << ", ";
        vector<vector<double>> Q = thresholdJacobi(A);
        cout << ", 运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
        //vector<double> eigenvalues;
        //for (int i = 0; i < n;i++) {
        //    eigenvalues.push_back(A[i][i]);
        //}
        ////PrintVector(eigenvalues);
        //sort(eigenvalues.begin(), eigenvalues.end());
        ////cout << "特征值为：" << endl;
        //PrintVector(eigenvalues);
        cout << "Ak= " << endl;
        PrintMatrix(A);
        cout << "Qk= " << endl;
        PrintMatrix(Q);
    }
}

void exercise7_2_1(){
    int n = 10;
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n - 1; i++)
    {
        A[i][i] = 2;
        A[i + 1][i] = -1;
        A[i][i + 1] = -1;
    }
    A[n - 1][n - 1] = 2;
    //PrintMatrix(A);

    // 测试特征值
    vector<vector<double>> Q = thresholdJacobi(A);
    cout << endl;
    vector<double> eigenvalues;
    for (int i = 0; i < n;i++) {
        eigenvalues.push_back(A[i][i]);
    }
    //PrintVector(eigenvalues);
    sort(eigenvalues.begin(), eigenvalues.end());
    //cout << "特征值为：" << endl;
    PrintVector(eigenvalues);
    

    for (int i = 0; i < n - 1; i++)
    {
        A[i][i] = 2;
        A[i + 1][i] = -1;
        A[i][i + 1] = -1;
    }
    A[n - 1][n - 1] = 2;
    //PrintMatrix(A);
    // 测试变号数
    //cout << reversals(A, 1.0) << endl;
    // 测试二分法求特征值
    for (int i = 1;i <= n;i++) {
        cout << dichEigenvalue(A, i)<<"\t";
    }
    cout << endl;
    // 测试反幂法求特征向量
    //cout << dichEigenvalue(A, 10) << endl; // 输出二分法求得的最大的特征值
    //vector<double> z = inversePowerMethod(A, dichEigenvalue(A, 10)); // 用反幂法求A最大的特征值对应的特征向量
    //vector<double> z = inversePowerMethod(A, 3.9190);
    for (int i = 1;i <= n;i++) {
        vector<double> z = inversePowerMethod(A, dichEigenvalue(A, i));
        PrintVector(z);
    }
}

void exercise7_2_2(){
    clock_t start_time;
    int n = 100;
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n - 1; i++)
    {
        A[i][i] = 2;
        A[i + 1][i] = -1;
        A[i][i + 1] = -1;
    }
    A[n - 1][n - 1] = 2;
    start_time = clock();
    cout << "最小特征值为：";
    cout << dichEigenvalue(A, 1) << endl;
    vector<double> z = inversePowerMethod(A, dichEigenvalue(A, 1));
    PrintVector(z);
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
    start_time = clock();
    cout << "最大特征值为：";
    cout << dichEigenvalue(A, n) << endl;
    z = inversePowerMethod(A, dichEigenvalue(A, n));
    PrintVector(z);
    cout << "运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds.\n" << endl;
}