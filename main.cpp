#include <iostream>
#include <iomanip>
#include <fstream>
#include <windows.h>
#include <gmp.h>
#include <gmpxx.h>

#define M 100

int N = 0;
typedef mpq_class matrix[M][M + 1];

using namespace std;

void SwapRow(matrix mat, int i, int j)
{
    for (int k = 0; k <= N; k++)
    {
        mpq_class temp = mat[i][k];
        mat[i][k] = mat[j][k];
        mat[j][k] = temp;
    }
}

int ForwardElimination(matrix mat)
{
    for (int k = 0; k < N; k++)
    {
        int i_max = k;
        mpq_class v_max = mat[i_max][k];

        for (int i = k + 1; i < N; i++)
            if (abs(mat[i][k]) > abs(v_max))
                v_max = mat[i][k], i_max = i;

        if (!mat[i_max][k])
            return k;

        if (i_max != k)
            SwapRow(mat, k, i_max);

        for (int i = k + 1; i < N; i++)
        {
            mpq_class d = mat[i][k] / mat[k][k];
            for (int j = k + 1; j <= N; j++)
                mat[i][j] -= mat[k][j] * d;

            mat[i][k] = 0;
        }
    }
    return -1;
}

void BackSubstitution(matrix mat, matrix mat1)
{
    mpq_class x[N];
    unsigned int max_len = 0;
    for (int i = N - 1; i >= 0; i--)
    {
        x[i] = mat[i][N];

        for (int j = i + 1; j < N; j++)
            x[i] -= mat[i][j] * x[j];

        x[i] = x[i] / mat[i][i];
        if (max_len < x[i].get_str().length())
            max_len = x[i].get_str().length();
    }

    cout << setprecision(17);
    for (int i = 0; i < N; i++)
        cout << left << setw(max_len + 3) << x[i] << setw(30) << x[i].get_d() << endl;

    bool bError = false;
    for (int i = 0; i < N; i++)
    {
        mpq_class d = 0;
        for (int j = 0; j < N; j++)
            d += mat1[i][j] * x[j];
        if (d != mat1[i][N])
        {
            cout << "Error in equation ¹" << i + 1 << ' ' << d << ' ' << mat1[i][N] << endl;
            bError = true;
        }
    }

    if (!bError)
        cout << "No errors" << endl;
}

int main(int argc, char *argv[])
{
    SetConsoleOutputCP(1251);
    matrix mat, mat1;

    ifstream iFile;
    if (argc == 1)
        iFile.open("input.txt");
    else
        iFile.open(argv[1]);

    if (!iFile.is_open()) {
        cout << "Can't open input file" << endl;
        return 1;
    }

    iFile >> N;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N + 1; j++) {
            iFile >> mat[i][j];
            mat1[i][j] = mat[i][j];
        }
    iFile.close();

    int singular_flag = ForwardElimination(mat);
    if (singular_flag != -1)
    {
        cout << "Singular matrix" << endl;
        if (mat[singular_flag][N])
            cout << "Inconsistent system" << endl;
        else
            cout << "The system has an infinite number of solutions" << endl;
        return 1;
    }

    BackSubstitution(mat, mat1);
    return 0;
}
