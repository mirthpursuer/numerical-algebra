#include "exercise.h"
#include <random>
clock_t start, ende;

void exercise_1()
{
    int N = 84; //矩阵大小
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
    //print(A);
    for (int i = 0; i < N; i++)
        cout << "b[" << i << "] = " << b[i] << endl;

    clock_t start_time;  // 记录开始时间

    //不选主元
    //读取 A,b
    std::vector<std::vector<double>> A1 = A;
    vector<double> b1 = b;

    start_time = clock();  // 记录开始时间

    gauss_elim(A1);
    //print(A1);
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
    //print(A2);
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
    std::vector<int> P;

    start_time = clock();  // 重新记录开始时间

    gauss_elim_col_pivoting(A3, P);
    //print(A3);
    forward_subs1(A3, b3);
    back_subs(A3, b3);
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

void exercise_2_1()
{
    //   初始化A和b
    int N = 10; //矩阵大小
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
    //print(A);

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
    //print(A);

    vector<vector<double>> AT(N, vector<double>(N));
    transposeMatrix(A, AT);
    //print(AT);

    forward_subs1(A, b);
    back_subs(AT, b);

    cout << "平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b[i] << endl;
    //cout <<  b[i] << endl;
    cout << "平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    start_time = clock();  // 记录开始时间

    modified_cholesky_decomp(A2);
    //print(A2);

    vector<vector<double>> A2T(N, vector<double>(N));
    transposeMatrix(A2, A2T);
    //print(A2T);

    forward_subs1(A2, b2);
    back_subs(A2T, b2);

    cout << "改进平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b2[i] << endl;
    //cout << b2[i] << endl;
    cout << "改进平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

}

void exercise_2_2()
{
    int N = 10; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);
    Hilbert_Matrix(A);
    print(A);
    for (int i = 0; i < N; ++i) {
        for (int j = 1; j <= N; ++j) {
            b[i] += 1.0 / (i + j);
        }
    }
    /*for (const auto& elem : b) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;*/


    //备份A和b，为后续改进平方根法所用
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;

    clock_t start_time;  // 记录开始时间
    start_time = clock();  // 记录开始时间

    cholesky_decomp(A);
    print(A);


    vector<vector<double>> AT(N, vector<double>(N));
    transposeMatrix(A, AT);
    print(AT);

    forward_subs1(A, b);
    back_subs(AT, b);

    cout << "平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b[i] << endl;
    cout << "平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    start_time = clock();  // 记录开始时间

    modified_cholesky_decomp(A2);
    print(A2);

    vector<vector<double>> A2T(N, vector<double>(N));
    transposeMatrix(A2, A2T);
    print(A2T);

    forward_subs1(A2, b2);
    back_subs(A2T, b2);

    cout << "改进平方根法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i << "] = " << b2[i] << endl;
    cout << "改进平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

}


//void exercise_2_2()
//{
//    int N = 40; //矩阵大小
//    vector<vector<double>> A(N, vector<double>(N));
//    vector<double> b(N);
//    Hilbert_Matrix(A);
//    //print(A);
//    for (int i = 0; i < N; ++i) {
//        for (int j = 1; j <= N; ++j) {
//            b[i] += 1.0 / (i + j);
//        }
//    }
//    /*for (const auto& elem : b) {
//        std::cout << elem << " ";
//    }
//    std::cout << std::endl;*/
//
//
//    //备份A和b，为后续改进平方根法所用
//    std::vector<std::vector<double>> A2 = A;
//    std::vector<double> b2 = b;
//
//    clock_t start_time;  // 记录开始时间
//    start_time = clock();  // 记录开始时间
//
//    cholesky_decomp(A);
//    //print(A);
//
//
//    vector<vector<double>> AT(N, vector<double>(N));
//    transposeMatrix(A, AT);
//    //print(AT);
//
//    forward_subs1(A, b);
//    back_subs(AT, b);
//
//    cout << "平方根法解：" << endl;
//    for (int i = 0; i < N; i++)
//        //cout << "x[" << i << "] = " << b[i] << endl;
//        cout <<  b[i] << endl;
//    cout << "平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
//
//    start_time = clock();  // 记录开始时间
//
//    modified_cholesky_decomp(A2);
//    //print(A2);
//
//    vector<vector<double>> A2T(N, vector<double>(N));
//    transposeMatrix(A2, A2T);
//    //print(A2T);
//
//    forward_subs1(A2, b2);
//    back_subs(A2T, b2);
//
//    cout << "改进平方根法解：" << endl;
//    for (int i = 0; i < N; i++)
//        //cout << "x[" << i << "] = " << b2[i] << endl;
//        cout <<  b2[i] << endl; 
//    cout << "改进平方根法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
//
//}

void exercise_3_1()
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
    //print(A);

    // 初始化随机数生成器
    std::random_device rd;  // 随机种子
    std::mt19937 gen(rd());  // Mersenne Twister 生成器
    std::uniform_real_distribution<double> dis(0.0, 1.0);  // 生成 [0.0, 1.0) 范围内的随机 double 数

    // 生成随机 double 数并填充向量
    for (int i = 0; i < N; ++i) {
        b[i] = dis(gen);
    }

    // 输出向量的内容
    std::cout << "Randomly generated vector:" << std::endl;
    for (const auto& elem : b) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    //不选主元
    //读取 A,b
    std::vector<std::vector<double>> A1 = A;
    vector<double> b1 = b;

    start_time = clock();  // 记录开始时间

    gauss_elim(A1);
    //print(A1);
    forward_subs1(A1, b1);
    back_subs(A1, b1);
    cout << "不选主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        //   cout << "x[" << i << "] = " << b1[i] << endl;
        cout << b1[i] << endl;

    cout << "不选主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    // 全主元
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;
    std::vector<int> u(N);
    std::vector<int> v(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_full_pivoting(A2, u, v);
    //print(A2);
    forward_subs1(A2, b2);
    back_subs(A2, b2);
    cout << "全主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << b2[i] << endl;
    //cout << "x[" << i << "] = " << b2[i] << endl;

    cout << "全主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    //列主元
    std::vector<std::vector<double>> A3 = A;
    std::vector<double> b3 = b;
    std::vector<int> P;

    start_time = clock();  // 重新记录开始时间

    gauss_elim_col_pivoting(A3, P);
    //print(A3);
    forward_subs1(A3, b3);
    back_subs(A3, b3);
    cout << "列主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        cout << b3[i] << endl;
    //cout << "x[" << i << "] = " << b3[i] << endl;

    cout << "列主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

}

void exercise_3_2()
{
    int N = 40; //矩阵大小
    vector<vector<double>> A(N, vector<double>(N));
    vector<double> b(N);
    Hilbert_Matrix(A);
    //print(A);
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
    //print(A1);
    forward_subs1(A1, b1);
    back_subs(A1, b1);
    cout << "不选主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        //   cout << "x[" << i << "] = " << b1[i] << endl;
        cout << b1[i] << endl;
    cout << "不选主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    // 全主元
    std::vector<std::vector<double>> A2 = A;
    std::vector<double> b2 = b;
    std::vector<int> u(N);
    std::vector<int> v(N);

    start_time = clock();  // 重新记录开始时间

    gauss_elim_full_pivoting(A2, u, v);
    //print(A2);
    forward_subs1(A2, b2);
    back_subs(A2, b2);
    cout << "全主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        //   cout << "x[" << i << "] = " << b2[i] << endl;
        cout << b2[i] << endl;
    cout << "全主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;

    //列主元
    std::vector<std::vector<double>> A3 = A;
    std::vector<double> b3 = b;
    std::vector<int> P;

    start_time = clock();  // 重新记录开始时间

    gauss_elim_col_pivoting(A3, P);
    //print(A3);
    forward_subs1(A3, b3);
    back_subs(A3, b3);
    cout << "列主元的Gauss消去法解：" << endl;
    for (int i = 0; i < N; i++)
        //   cout << "x[" << i << "] = " << b3[i] << endl;
        cout << b3[i] << endl;
    cout << "列主元的Gauss消去法运行时间: " << double(clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl;
}