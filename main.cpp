#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <regex>
using namespace std;

int PRECISION;

class Matrix {
protected:
    int n, m;
public:
    vector<vector<double>> matrixData;
    Matrix(int n, int m) : n(n), m(m), matrixData(vector<vector<double>>(n, vector<double>(m, 0))) {}

    friend istream& operator>>(istream& in, Matrix& matrixToRead) {
        for (int i = 0; i < matrixToRead.n; ++i) {
            for (int j = 0; j < matrixToRead.m; ++j) {
                in >> matrixToRead.matrixData[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const Matrix& matrix) {
        for (int i = 0; i < matrix.n; ++i) {
            for (int j = 0; j < matrix.m; ++j) {
                out << fixed << setprecision(PRECISION) << matrix.matrixData[i][j];
                if (j != matrix.m - 1) {
                    out << " ";
                }
            }
            out << endl;
        }
        return out;
    }
    Matrix& operator=(Matrix& other) {
        n = other.n;
        m = other.m;
        matrixData = other.matrixData;
        return *this;
    }

    virtual Matrix* operator+(Matrix& other) {
        Matrix* newMatrix = new Matrix(n, m);
        if (n == other.n && m == other.m) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    newMatrix->matrixData[i][j] = matrixData[i][j] + other.matrixData[i][j];
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual Matrix* operator-(Matrix& other) {
        Matrix* newMatrix = new Matrix(n, m);
        if (n == other.n && m == other.m) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    newMatrix->matrixData[i][j] = matrixData[i][j] - other.matrixData[i][j];
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual Matrix* operator*(Matrix& other) {
        Matrix* newMatrix = new Matrix(n, other.m);
        if (m == other.n) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < other.m; j++) {
                    for (int p = 0; p < m; p++) {
                        newMatrix->matrixData[i][j] += matrixData[i][p] * other.matrixData[p][j];
                    }
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual Matrix* transpose() {
        Matrix* newMatrix = new Matrix(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                newMatrix->matrixData[j][i] = matrixData[i][j];
            }
        }
        return newMatrix;
    }

    int getRows() const { return n; }

    int getColumns() const { return m; }
};


void makeBasicColumn(Matrix& table, int basic_row, int basic_col) {
    // Make basic element 1
    double basic_element = table.matrixData[basic_row][basic_col];
    for (int i = 0; i < table.getColumns(); i++) {
        table.matrixData[basic_row][i] /= basic_element;
    }


    // Nullify all elements except basic row
    for (int i = 0; i < table.getRows(); i++) {
        if (i == basic_row) {
            continue;
        }

        double coefficient = table.matrixData[i][basic_col];
        for (int j = 0; j < table.getColumns(); j++) {
            table.matrixData[i][j] -= coefficient * table.matrixData[basic_row][j];
        }

    }
}

// Basic test
//3 6 1
//2 3 0 -1 0 0 0
//2 -1 0 -2 1 0 16
//3 2 1 -3 0 0 18
//-1 3 0 4 0 1 24


class Answer {
    bool solver_sate;
    Matrix solution;
    float z;

    Answer(bool state, const Matrix& sol, float z): solver_sate(state), solution(sol), z(z) {}
};

Answer solveLLP(Matrix& C, Matrix& A, Matrix& b, int eps) {


}

Matrix* concatenate(Matrix& A, Matrix& B) {
    Matrix* newMatrix = new Matrix(A.getRows(), A.getColumns() + B.getColumns());
    for (int i = 0; i < A.getRows(); i++) {
        for (int j = 0; j < A.getColumns(); j++) {
            newMatrix->matrixData[i][j] = A.matrixData[i][j];
        }
        for (int j = 0; j < B.getColumns(); j++) {
            newMatrix->matrixData[i][A.getColumns()+j] = B.matrixData[i][j];
        }
    }
    return newMatrix;
}

int main() {
    int n_constrains, n_var;
    cin >> n_constrains >> n_var >> PRECISION;

    Matrix table(1 + n_constrains, n_var + 1);
    cin >> table;

    Make Z coeff negative

    for (int i = 0; i < n_var; i++) {
        table.matrixData[0][i] *= -1;
    }
    // cout << table;

    for (int _ = 0; _ < n_constrains; _++) {
        int basic_col = -1;
        int min_col = 0;
        for (int i = 0; i < n_var; i++) {
            if (table.matrixData[0][i] < min_col) {
                min_col = table.matrixData[0][i];
                basic_col = i;
            }
        }
        if (basic_col == -1) {
            break;
        };



        int basic_row = 0;
        int min_ratio = pow(10, 10);
        for (int i = 1; i < n_constrains+1; i++) {
            if (table.matrixData[i][n_var] / table.matrixData[i][basic_col] < min_ratio &&
                table.matrixData[i][n_var] / table.matrixData[i][basic_col] > 0) {
                basic_row = i;
                min_ratio = table.matrixData[i][n_var] / table.matrixData[i][basic_col];
                }
        }
        if (basic_row == 0) {
            cout << "The method is not applicable!\n";
            break;
        }
        // cout << basic_col << ' ' << basic_row << '\n';
        makeBasicColumn(table, basic_row, basic_col);
        // cout << table;
    }
    cout << table;
}