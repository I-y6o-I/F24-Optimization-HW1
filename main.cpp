#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
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
    Matrix operator=(Matrix& other) {
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


class Answer {
public:
    bool solver_sate;
    Matrix solution;
    double z;

public:
    Answer(bool state, const Matrix& sol, float z): solver_sate(state), solution(sol), z(z) {}
};


Matrix* concatenate(Matrix& A, Matrix& B, bool dim) {
    Matrix* newMatrix;
    if (dim){
        newMatrix = new Matrix(A.getRows(), A.getColumns() + B.getColumns());
        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                newMatrix->matrixData[i][j] = A.matrixData[i][j];
            }
            for (int j = 0; j < B.getColumns(); j++) {
                newMatrix->matrixData[i][A.getColumns()+j] = B.matrixData[i][j];
            }
        }
    }
    else{
        newMatrix = new Matrix(A.getRows()+B.getRows(), A.getColumns());
        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                newMatrix->matrixData[i][j] = A.matrixData[i][j];
            }
        }
        for (int i = 0; i < B.getRows(); i++) {
            for (int j = 0; j < B.getColumns(); j++){
                newMatrix->matrixData[i+A.getRows()][j] = B.matrixData[i][j];
            }
        }
    }

    return newMatrix;
}

Matrix getSolution(Matrix& A, Matrix& solution) {
    for (int j = 0; j < A.getColumns() - 1; j++) {
        int countOne = 0, oneInd = -1;
        bool flag = false;
        for (int i = 0; i < A.getRows(); i++) {
            if (A.matrixData[i][j] == 1 && countOne == 0) {
                countOne++;
                oneInd = i;
            }
            else if (A.matrixData[i][j] == 0) {
                continue;
            }
            else {
                flag = true;
            }
        }
        if (countOne == 1 && !flag) {
            solution.matrixData[0][j] = A.matrixData[oneInd][A.getColumns() - 1];
        }
    }
    return solution;
}



Answer solveLLP(Matrix& C, Matrix& A, Matrix& b, int eps) {
    PRECISION = eps;

    Matrix temp_col = *concatenate(*new Matrix(1, 1), b, false);
    Matrix table = *concatenate(*concatenate(C, A, false), temp_col, true);

    int n_var = C.getColumns();
    int n_constrains = A.getRows();



    for (int i = 0; i < n_var; i++) {
        table.matrixData[0][i] *= -1;
    }
//    cout << table << endl;

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
            return Answer(false, Matrix(1, n_var), -123123);
        }
        makeBasicColumn(table, basic_row, basic_col);

    }

    Matrix solution = getSolution(table, *new Matrix(1, n_var));
    for (int i = 0; i < solution.getColumns(); i++) {
        if (solution.matrixData[0][i] < 0) {
            return Answer(false, Matrix(1, n_var), -123123);
        }
    }

    return Answer(true, solution, table.matrixData[0][n_var]);

}

bool compareDoubleVectors(const std::vector<double>& vec1, const std::vector<double>& vec2, double epsilon) {
    if (vec1.size() != vec2.size()) {
        return false;
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        if (std::fabs(vec1[i] - vec2[i]) >= epsilon) {
            return false;
        }
    }

    return true;
}

bool compareDouble(double a, double b, double epsilon) {
    return std::fabs(a - b) < epsilon;
}



void runTests() {
//    Test 1

    Matrix C(1, 6);
    Matrix A(3, 6);
    Matrix b(3, 1);

    C.matrixData[0] = vector<double> {9, 10, 16, 0, 0, 0};
    A.matrixData[0] = vector<double> {18, 15, 12, 1, 0, 0};
    A.matrixData[1] = vector<double> {6, 4, 8, 0, 1, 0};
    A.matrixData[2] = vector<double> {5, 3, 3, 0, 0, 1};
    b.matrixData[0][0] = 360;
    b.matrixData[1][0] = 192;
    b.matrixData[2][0] = 180;

    Answer ans1 = solveLLP(C, A, b, 3);

    if (ans1.solver_sate
        && compareDoubleVectors(ans1.solution.matrixData[0], vector<double> {0, 8, 20, 0, 0, 96}, 1e-3)
        && compareDouble(ans1.z, 400, 1e-3)) {
        cout << "Test 1 passed\n";
    } else {
        cout << "Test 1 failed\n";
    }

// Test 2

    C = *(new Matrix(1, 6));
    A = *(new Matrix(3, 6));
    b = *(new Matrix(3, 1));

    C.matrixData[0] = vector<double> {2, 3, 0, -1, 0, 0};
    A.matrixData[0] = vector<double> {2, -1, 0, -2, 1, 0};
    A.matrixData[1] = vector<double> {3, 2, 1, -3, 0, 0};
    A.matrixData[2] = vector<double> {-1, 3, 0, 4, 0, 1};
    b.matrixData[0][0] = 16;
    b.matrixData[1][0] = 18;
    b.matrixData[2][0] = 24;

    Answer ans2 = solveLLP(C, A, b, 4);

    if (ans2.solver_sate
        && compareDoubleVectors(ans2.solution.matrixData[0], vector<double> {0.5455, 8.1818, 0, 0, 23.0909, 0}, 1e-4)
        && compareDouble(ans2.z, 25.6364, 1e-4)) {
        cout << "Test 2 passed\n";
    } else {
        cout << "Test 2 failed\n";
    }

//    Test 3
    C = *(new Matrix(1, 5));
    A = *(new Matrix(3, 5));
    b = *(new Matrix(3, 1));

    C.matrixData[0] = vector<double> {1, 1, 0, 0, 0};
    A.matrixData[0] = vector<double> {2, 4, 1, 0, 0};
    A.matrixData[1] = vector<double> {-4, 2, 0, 1, 0};
    A.matrixData[2] = vector<double> {1, 3, 0, 0, 1};
    b.matrixData[0][0] = 16;
    b.matrixData[1][0] = 8;
    b.matrixData[2][0] = 9;

    Answer ans3 = solveLLP(C, A, b, 4);
    if (ans3.solver_sate
        && compareDoubleVectors(ans3.solution.matrixData[0], vector<double> {8.0000, 0, 0, 40.0000, 1.0000}, 1e-4)
        && compareDouble(ans3.z, 8.0000, 1e-4)) {
        cout << "Test 3 passed\n";
    } else {
        cout << "Test 3 failed\n";
    }


//    Test 4

    C.matrixData[0] = vector<double> {1, 1, 0, 0, 0};
    A.matrixData[0] = vector<double> {-2, -2, 1, 0, 0};
    A.matrixData[1] = vector<double> {-5, 3, 0, 1, 0};
    A.matrixData[2] = vector<double> {4, 6, 0, 0, 1};
    b.matrixData[0][0] = -14;
    b.matrixData[1][0] = 15;
    b.matrixData[2][0] = 24;

    Answer ans4 = solveLLP(C, A, b, 4);

    if (!ans4.solver_sate) {
        cout << "Test 4 passed. The method is not applicable!\n";
    } else {
        cout << "Test 4 failed\n";
    }

//    Test 5

    C = *(new Matrix(1, 4));
    A = *(new Matrix(2, 4));
    b = *(new Matrix(2, 1));

    C.matrixData[0] = vector<double> {1, -3, 0, 0};
    A.matrixData[0] = vector<double> {1, 7, 1, 0};
    A.matrixData[1] = vector<double> {-6, 5, 0, 1};
    b.matrixData[0][0] = 14;
    b.matrixData[1][0] = 15;

    Answer ans5 = solveLLP(C, A, b, 2);

    if (ans5.solver_sate
        && compareDoubleVectors(ans5.solution.matrixData[0], vector<double> {14, 0, 0, 99}, 1e-2)
        && compareDouble(ans5.z, 14, 1e-2)) {
        cout << "Test 5 passed\n";
    } else {
        cout << "Test 5 failed\n";
    }

}


int main() {
    runTests();

    cout << "Input LLP" << endl;
    int n_constrains, n_var;
    cin >> n_constrains >> n_var >> PRECISION;
    Matrix C(1, n_var);
    Matrix A(n_constrains, n_var);
    Matrix b(n_constrains, 1);
    cin >> C >> A >> b;


    Answer ans = solveLLP(C, A, b, 2);
    if (ans.solver_sate) {
        cout << ans.solver_sate << endl;
        cout << ans.solution;
        cout << ans.z << endl;
    } else {
        cout << "The method is not applicable!" << endl;
    }
}