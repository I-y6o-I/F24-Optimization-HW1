#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
using namespace std;

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
                out << fixed << setprecision(4) << matrix.matrixData[i][j];
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


int main()
{
    int n;
    cin >> n;
    SquareMatrix A(n);
    cin >> A;
    SquareMatrix Lower(n);
    SquareMatrix Upper(n);
    SquareMatrix Diagonal(n);
    int m;
    cin >> m;
    Matrix b(m, 1);
    cin >> b;
    double epsStart; cin >> epsStart;
    for (int i = 0; i < n; i++) {
        Diagonal.matrixData[i][i] = A.matrixData[i][i];
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            Upper.matrixData[i][j] = A.matrixData[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            Lower.matrixData[j][i] = A.matrixData[j][i];
        }
    }
    SquareMatrix A1 = Diagonal;
    int countPerm = 0;
    for (int i = 0; i < n; i++) {
        int IndMaxEl = i;
        double maxEl = abs(A1.matrixData[i][i]);
        for (int j = i + 1; j < n; j++) {
            if (abs(A1.matrixData[j][i]) > maxEl) {
                IndMaxEl = j;
                maxEl = abs(A1.matrixData[j][i]);
            }
        }
        if (IndMaxEl > i) {

            SquareMatrix P = PermutationMatrix(n, i + 1, IndMaxEl + 1);
            SquareMatrix* temp = P * A1;
            A1 = *temp;
            countPerm++;
        }
        for (int j = i + 1; j < n; j++) {
            if (A1.matrixData[j][i] != 0) {
                SquareMatrix E = EliminationMatrix(n, j + 1, i + 1, A1);
                SquareMatrix* temp = E * A1;
                A1 = *temp;
            }
        }
    }
    double ans = 1;
    for (int i = 0; i < n; i++) {
        ans *= A1.matrixData[i][i];
    }
    if (countPerm % 2 != 0) {
        ans *= -1;
    }
    if (ans == 0) {
        cout << "The method is not applicable";
        return 0;
    }
    if (n != m) {
        cout << "The method is not applicable";
        return 0;
    }


    Matrix AugmentedMatrix(n, 2 * n);
    IdentityMatrix I(n);
    for (int i = 0; i < n; i++) {
        AugmentedMatrix.matrixData[i][i] = Diagonal.matrixData[i][i];
        for (int j = n; j < 2 * n; j++) {
            AugmentedMatrix.matrixData[i][j] = I.matrixData[i][j - n];
        }
    }
    int k = 1;
    for (int i = 0; i < n; i++) {
        int IndMaxEl = i;
        double maxEl = abs(Diagonal.matrixData[i][i]);
        for (int j = i + 1; j < n; j++) {
            if (abs(Diagonal.matrixData[j][i]) > maxEl) {
                IndMaxEl = j;
                maxEl = abs(Diagonal.matrixData[j][i]);
            }
        }
        if (IndMaxEl > i) {
            SquareMatrix P = PermutationMatrix(n, i + 1, IndMaxEl + 1);
            Matrix P1 = P;
            SquareMatrix* temp = P * Diagonal;
            Matrix* temp1 = P1 * AugmentedMatrix;
            Diagonal = *temp;
            AugmentedMatrix = *temp1;
            k++;
        }
        for (int j = i + 1; j < n; j++) {
            if (Diagonal.matrixData[j][i] != 0) {
                SquareMatrix E = EliminationMatrix(n, j + 1, i + 1, Diagonal);
                Matrix E1 = E;
                Matrix* temp = E1 * AugmentedMatrix;
                SquareMatrix* temp1 = E * Diagonal;
                AugmentedMatrix = *temp;
                Diagonal = *temp1;
                k++;
            }
        }
    }
    for (int i = n - 1; i > -1; i--) {
        for (int j = i - 1; j > -1; j--) {
            if (Diagonal.matrixData[j][i] != 0) {
                SquareMatrix E = EliminationMatrix(n, j + 1, i + 1, Diagonal);
                Matrix E1 = E;
                Matrix* temp = E1 * AugmentedMatrix;
                SquareMatrix* temp1 = E * Diagonal;
                AugmentedMatrix = *temp;
                Diagonal = *temp1;
                k++;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double diag = AugmentedMatrix.matrixData[i][i];
        for (int j = 0; j < 2 * n; j++) {
            AugmentedMatrix.matrixData[i][j] = AugmentedMatrix.matrixData[i][j] / diag;
        }
    }
    SquareMatrix InverseDiagonal(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            InverseDiagonal.matrixData[i][j] = AugmentedMatrix.matrixData[i][j + n];
        }
    }
    Matrix* UplusL = Upper + Lower;
    Matrix* matrixAlpha = InverseDiagonal * *UplusL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (matrixAlpha->matrixData[i][j] != 0) {
                matrixAlpha->matrixData[i][j] = -matrixAlpha->matrixData[i][j];
            }
        }
    }
    cout << "alpha:\n";
    cout << *matrixAlpha;
    Matrix* matrixBeta = (Matrix)InverseDiagonal * b;
    cout << "beta:\n";
    cout << *matrixBeta;
    Matrix B(n,n);
    Matrix C(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            B.matrixData[j][i] = matrixAlpha->matrixData[j][i];
            C.matrixData[i][j] = matrixAlpha->matrixData[i][j];
        }
    }
    cout << "B:" << endl;
    cout << B;
    cout << "C:" << endl;
    cout << C;
    SquareMatrix I_B = *(SquareMatrix*)((Matrix)I - B);
    cout << "I-B:" << endl;
    cout << I_B;
    
    Matrix AugmentedMatrix1(n, 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = n; j < 2 * n; j++) {
            AugmentedMatrix1.matrixData[i][j - n] = I_B.matrixData[i][j - n];
            AugmentedMatrix1.matrixData[i][j] = I.matrixData[i][j - n];
        }
    }
    k = 1;
    for (int i = 0; i < n; i++) {
        int IndMaxEl = i;
        double maxEl = abs(I_B.matrixData[i][i]);
        for (int j = i + 1; j < n; j++) {
            if (abs(I_B.matrixData[j][i]) > maxEl) {
                IndMaxEl = j;
                maxEl = abs(I_B.matrixData[j][i]);
            }
        }
        if (IndMaxEl > i) {
            SquareMatrix P = PermutationMatrix(n, i + 1, IndMaxEl + 1);
            Matrix P1 = P;
            SquareMatrix* temp = P * I_B;
            Matrix* temp1 = P1 * AugmentedMatrix1;
            I_B = *temp;
            AugmentedMatrix1 = *temp1;
            k++;
        }
        for (int j = i + 1; j < n; j++) {
            if (I_B.matrixData[j][i] != 0) {
                SquareMatrix E = EliminationMatrix(n, j + 1, i + 1, I_B);
                Matrix E1 = E;
                Matrix* temp = E1 * AugmentedMatrix1;
                SquareMatrix* temp1 = E * I_B;
                AugmentedMatrix1 = *temp;
                I_B = *temp1;
                k++;
            }
        }
    }
    for (int i = n - 1; i > -1; i--) {
        for (int j = i - 1; j > -1; j--) {
            if (I_B.matrixData[j][i] != 0) {
                SquareMatrix E = EliminationMatrix(n, j + 1, i + 1, I_B);
                Matrix E1 = E;
                Matrix* temp = E1 * AugmentedMatrix1;
                SquareMatrix* temp1 = E * I_B;
                AugmentedMatrix1 = *temp;
                I_B = *temp1;
                k++;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double diag = AugmentedMatrix1.matrixData[i][i];
        for (int j = 0; j < 2 * n; j++) {
            AugmentedMatrix1.matrixData[i][j] = AugmentedMatrix1.matrixData[i][j] / diag;
        }
    }
    Matrix I_B_1(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            I_B_1.matrixData[i][j] = AugmentedMatrix1.matrixData[i][j + n];
        }
    }
    cout << "(I-B)_-1:" << endl;
    cout << I_B_1;




    int count = 1;
    Matrix xPrev = *matrixBeta;
    Matrix xNew = *(*(*(I_B_1*C)*xPrev)+*(I_B_1**matrixBeta));
    double epsilon = 0;
    cout << "x";
    for (int i = 0; i < m; i++) {
        epsilon += pow(xNew.matrixData[i][0] - xPrev.matrixData[i][0], 2);
    }
    if (sqrt(epsilon) < epsStart) {
        cout << "(" << count << "):" << endl;
        cout << xNew;
        cout << "e:" << sqrt(epsilon) << endl;
        cout << "x~:" << endl;
        cout << xNew;
        return 0;
    }
    else {
        cout << "(" << count << "):" << endl;
        cout << xNew;
        count++;
        cout << "e:" << sqrt(epsilon) << endl;
        epsilon = 0;
    }
    xPrev = xNew;
    while (true) {
        cout << "x";
        Matrix xNew = *(*(*(I_B_1 * C) * xPrev) + *(I_B_1 * *matrixBeta));
        epsilon = 0;
        for (int i = 0; i < m; i++) {
            epsilon += pow(xNew.matrixData[i][0] - xPrev.matrixData[i][0], 2);
        }
        if (sqrt(epsilon) < epsStart) {
            cout << "(" << count << "):" << endl;
            cout << xNew;
            cout << "e:" << sqrt(epsilon) << endl;
            cout << "x~:" << endl;
            cout << xNew;
            return 0;
        }
        else {
            cout << "(" << count << "):" << endl;
            cout << xNew;
            cout << "e:" << sqrt(epsilon) << endl;
            xPrev = xNew;
            count++;
        }
    }
}