#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
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

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {}


    virtual SquareMatrix* operator+(Matrix& other) override {
        SquareMatrix* otherSquare = dynamic_cast<SquareMatrix*>(&other);

        SquareMatrix* newMatrix = new SquareMatrix(n);
        if (n == otherSquare->n) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    newMatrix->matrixData[i][j] = matrixData[i][j] + otherSquare->matrixData[i][j];
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual SquareMatrix* operator-(Matrix& other) override {
        SquareMatrix* otherSquare = dynamic_cast<SquareMatrix*>(&other);

        SquareMatrix* newMatrix = new SquareMatrix(n);
        if (n == otherSquare->n) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    newMatrix->matrixData[i][j] = matrixData[i][j] - otherSquare->matrixData[i][j];
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual SquareMatrix* operator*(Matrix& other) override {
        SquareMatrix* otherSquare = dynamic_cast<SquareMatrix*>(&other);

        SquareMatrix* newMatrix = new SquareMatrix(n);
        if (n == otherSquare->n) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        newMatrix->matrixData[i][j] += matrixData[i][k] * otherSquare->matrixData[k][j];
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

    virtual SquareMatrix* transpose() override {
        SquareMatrix* newMatrix = new SquareMatrix(n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                newMatrix->matrixData[j][i] = matrixData[i][j];
            }
        }

        return newMatrix;
    }
    SquareMatrix& operator=(SquareMatrix& other) {
        n = other.n;
        matrixData = other.matrixData;
        return *this;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            matrixData[i][i] = 1;
        }
    }
};
class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix(int n, int i, int j, SquareMatrix givenMatrix) : SquareMatrix(n) {
        for (int p = 0; p < n; p++) {
            matrixData[p][p] = 1;
        }
        matrixData[i - 1][j - 1] = -givenMatrix.matrixData[i - 1][j - 1] / givenMatrix.matrixData[j - 1][j - 1];
    }
};
class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int n, int i, int j) : SquareMatrix(n) {
        for (int p = 0; p < n; p++) {
            matrixData[p][p] = 1;
        }
        swap(matrixData[i - 1], matrixData[j - 1]);
    };
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


int main() {
    int n_constrains, n_var;
    cin >> n_constrains >> n_var >> PRECISION;

    Matrix table(1 + n_constrains, n_var);
    cin >> table;

    // Make Z coeff negative
    for (int i = 0; i < n_var; i++) {
        table.matrixData[0][i] *= -1;
    }

    makeBasicColumn(table, 3, 1);
    makeBasicColumn(table, 2, 0);
    cout << table;



}