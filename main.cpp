#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <random>

#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

using namespace std;

const double EPS = 1e-9;

class Matrix {
public:
    Matrix(vector<double> vector1) {
        rows = vector1.size();
        cols = 1;
        data = vector<vector<double>>(rows, vector<double>(cols));
        for (int i = 0; i < rows; i++) {
            data[i][0] = vector1[i];
        }

    }

    int rows{};
    int cols{};
    vector<vector<double>> data;

    Matrix() {}

    Matrix(int n, int m) : rows(n), cols(m), data(vector<vector<double>>(n, vector<double>(m))) {}

    friend istream &operator>>(istream &is, Matrix &mat) {
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                is >> mat.data[i][j];
            }
        }
        return is;
    }

    friend ostream &operator<<(ostream &os, const Matrix &mat) {
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                if(fabs(mat.data[i][j]) < EPS)
                    os << 0.0000;
                else
                    os << setprecision(4) << fixed << mat.data[i][j];
                if (j != mat.cols - 1) {
                    os << " ";
                }
            }
            os << "\n";
        }
        return os;
    }

    virtual Matrix operator=(const Matrix &mat) {
        rows = mat.rows;
        cols = mat.cols;
        data = mat.data;
        return *this;
    }

    friend Matrix operator+(const Matrix &mat1, const Matrix &mat2) {
        if (mat1.rows != mat2.rows || mat1.cols != mat2.cols) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        Matrix sum(mat1.rows, mat1.cols);
        for (int i = 0; i < mat1.rows; i++) {
            for (int j = 0; j < mat1.cols; j++) {
                sum.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
            }
        }
        return sum;
    }

    friend Matrix operator-(const Matrix &mat1, const Matrix &mat2) {
        if (mat1.rows != mat2.rows || mat1.cols != mat2.cols) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        Matrix diff(mat1.rows, mat1.cols);
        for (int i = 0; i < mat1.rows; i++) {
            for (int j = 0; j < mat1.cols; j++) {
                diff.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
            }
        }
        return diff;
    }

    friend Matrix operator*(const Matrix &mat1, const Matrix &mat2) {
        if (mat1.cols != mat2.rows) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        Matrix prod(mat1.rows, mat2.cols);
        for (int i = 0; i < mat1.rows; i++) {
            for (int j = 0; j < mat2.cols; j++) {
                for (int k = 0; k < mat1.cols; k++) {
                    prod.data[i][j] += mat1.data[i][k] * mat2.data[k][j];
                }
            }
        }
        return prod;
    }

    Matrix print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if(fabs(data[i][j]) < EPS)
                    cout << 0.0000;
                else
                    cout << setprecision(4) << fixed << data[i][j];
                if (j != cols - 1) {
                    cout << " ";
                }
            }
            cout << "\n";
        }
        return *this;
    }

    Matrix transpose() const {
        Matrix transposed(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed.data[j][i] = data[i][j];
            }
        }
        return transposed;
    }

    Matrix set(int i, int j, int val) {
        data[i][j] = val;
        return *this;
    }

    Matrix inverse() {
        if (rows != cols) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        Matrix inv(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (i == j) {
                    inv.data[i][j] = 1;
                } else {
                    inv.data[i][j] = 0;
                }
            }
        }
        for (int i = 0; i < rows; i++) {
            double pivot = data[i][i];
            for (int j = 0; j < cols; j++) {
                data[i][j] /= pivot;
                inv.data[i][j] /= pivot;
            }
            for (int j = 0; j < rows; j++) {
                if (i != j) {
                    double factor = data[j][i];
                    for (int k = 0; k < cols; k++) {
                        data[j][k] -= factor * data[i][k];
                        inv.data[j][k] -= factor * inv.data[i][k];
                    }
                }
            }
        }
        return inv;
    }
};
class ColumnVector {
private:
    vector<double> v;
    int size;
public:

    ColumnVector(int n) {
        size = n;
        v.resize(n);
    }

    friend istream &operator >>(istream &in, ColumnVector &v){
        for (int i = 0; i < v.size; i++) {
            in >> v[i];
        }
        return in;
    }
    friend ostream &operator <<(ostream &out, ColumnVector &v){
        for (int i = 0; i < v.size; i++) {
            if(fabs(v[i]) < EPS)
                out << 0.0000;
            else
                out << v[i];
            out << "\n";
        }
        return out;
    }

    double &operator[](int i) {
        return v[i];
    }

    ColumnVector operator+(const ColumnVector &other) const {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = v[i] + other.v[i];
        }
        return result;
    }

    ColumnVector operator*(double scalar) const {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = v[i] * scalar;
        }
        return result;
    }

    ColumnVector operator-(const ColumnVector &other) const {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = v[i] - other.v[i];
        }
        return result;
    }

    ColumnVector operator*(const Matrix &other) const {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                result[i] += v[j] * other.data[j][i];
            }
        }
        return result;
    }
};


int find_max_abs_element(vector<vector<double>> &matrix, int col) {
    int max_row = col;
    double max_val = matrix[col][col];
    for (int i = col + 1; i < matrix.size(); ++i) {
        if (fabs(matrix[i][col]) > fabs(max_val)) {
            max_row = i;
            max_val = matrix[i][col];
        }
    }
    return max_row;
}

void print_matrix(vector<vector<double>> &matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            cout << fixed << setprecision(4) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

class SquareMatrix : public Matrix {
public:
    SquareMatrix() {}

    SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix operator=(const SquareMatrix &mat) {
        rows = mat.rows;
        cols = mat.cols;
        data = mat.data;
        return *this;
    }

    friend SquareMatrix operator+(const SquareMatrix &mat1, const SquareMatrix &mat2) {
        if (mat1.rows != mat2.rows || mat1.cols != mat2.cols) {
            cout << "Error: the dimensional problem occurred\n";
            return SquareMatrix();
        }
        SquareMatrix sum(mat1.rows);
        for (int i = 0; i < mat1.rows; i++) {
            for (int j = 0; j < mat1.cols; j++) {
                sum.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
            }
        }
        return sum;
    }

    friend SquareMatrix operator-(const SquareMatrix &mat1, const SquareMatrix &mat2) {
        if (mat1.rows != mat2.rows || mat1.cols != mat2.cols) {
            cout << "Error: the dimensional problem occurred\n";
            return SquareMatrix();
        }
        SquareMatrix diff(mat1.rows);
        for (int i = 0; i < mat1.rows; i++) {
            for (int j = 0; j < mat1.cols; j++) {
                diff.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
            }
        }
        return diff;
    }

    friend SquareMatrix operator*(const SquareMatrix &mat1, const SquareMatrix &mat2) {
        if (mat1.cols != mat2.rows) {
            cout << "Error: the dimensional problem occurred\n";
            return SquareMatrix();
        }
        SquareMatrix prod(mat1.rows);
        for (int i = 0; i < mat1.rows; i++) {
            for (int j = 0; j < mat2.cols; j++) {
                int val = 0;
                for (int k = 0; k < mat1.cols; k++) {
                    val += mat1.data[i][k] * mat2.data[k][j];
                }
                prod.data[i][j] = val;
            }
        }
        return prod;
    }

    SquareMatrix transpose() const {
        SquareMatrix transposed(cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed.data[j][i] = data[i][j];
            }
        }
        return transposed;
    }

    void inverse_matrix(vector<vector<double>> &matrix, ColumnVector &b) {
        int n = matrix.size();
        int step = 1;

        for (int i = 0; i < n; ++i) {
            int max_row = find_max_abs_element(matrix, i);
            if (max_row != i) {
                swap(b[i], b[max_row]);
                swap(matrix[i], matrix[max_row]);
                cout << "step #" << step << ": permutation" << endl;
                print_matrix(matrix);
                cout << b;
                step++;
            }
            for (int j = i + 1; j < n; ++j) {
                double coef = matrix[j][i] / matrix[i][i];
                if(coef == 0.0)continue;
                for (int k = i; k < n; ++k) {
                    matrix[j][k] -= matrix[i][k] * coef;
                }
                b [j] -= b[i] * coef;
                cout << "step #" << step << ": elimination" << endl;
                print_matrix(matrix);
                cout << b;
                step++;
            }
        }

        for (int i = n - 1; i >= 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                double coef = matrix[j][i] / matrix[i][i];
                matrix[j][i] -= matrix[i][i] * coef;
                b[j] -= b[i] * coef;
                cout << "step #" << step << ": elimination" << endl;
                print_matrix(matrix);
                cout << b;
                step++;

            }
        }
        for (int i = 0; i < n; ++i) {
            b[i] /= matrix[i][i];
            matrix[i][i] = 1;
        }
        cout << "Diagonal normalization:\n";
        print_matrix(matrix);
        cout << b;

        cout << "result:" << endl;
        cout << b;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix() : SquareMatrix(3) {
        for (int i = 0; i < 3; i++) {
            data[i][i] = 1;
        }
    }
};

class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            data[i][i] = 1;
        }
    }

    void set(int i, int j, int val) {
        data[i][j] = -val;
    }
};

class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            data[i][i] = 1;
        }
    }

    void swap(int i, int j) {
        for (int k = 0; k < rows; k++) {
            int tmp = data[i][k];
            data[i][k] = data[j][k];
            data[j][k] = tmp;
        }
    }
};


int main() {
    #ifdef WIN32
        FILE* pipe = _popen(GNUPLOT_NAME, "w");
    #else
        FILE* pipe = popen(GNUPLOT_NAME, "w");
    #endif

    int m;
    cin >> m;
    vector<double> t(m);
    vector<double> b(m);
    for (int i = 0; i < m; i++) {
        cin >> t[i] >> b[i];
    }
    int n;
    cin >> n;
    Matrix A(m, n + 1);
    for (int i = 0; i < m; i++) {
        A.data[i][0] = 1;
        for (int j = 1; j <= n; j++) {
            A.data[i][j] = A.data[i][j - 1] * t[i];
        }
    }
    cout << "A:" << endl;
    cout << A;

    Matrix AT_A = A.transpose() * A;
    cout << "A_T*A:" << endl;
    cout << AT_A;

    Matrix AT_A_inv = AT_A.inverse();
    cout << "(A_T*A)^-1:" << endl;
    cout << AT_A_inv;

    Matrix AT_b = A.transpose() * b;
    cout << "A_T*b:" << endl;
    cout << AT_b;

    Matrix x = AT_A_inv * AT_b;
    cout << "x~:" << endl;
    cout << x;


    int X, Y;
    default_random_engine _random{random_device{}()};
    uniform_real_distribution<double> interval(X, Y);


    X = interval(_random);
    Y = interval(_random);


    if (n == 1) {
        double k1 = x.data[0][0];
        double k0 = x.data[1][0];
        fprintf(pipe, "plot [5 : 15] [0 : 5] %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", k1, k0);
        fprintf(pipe, "%f\t%f\n", X, Y);
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);

        #ifdef WIN32
            _pclose(pipe);
        #else
            pclose(pipe);
        #endif

        return 0;
    }

    if (n == 2) {
        double k2 = x.data[0][0];
        double k1 = x.data[1][0];
        double k0 = x.data[2][0];
        fprintf(pipe, "plot [5 : 15] [0 : 5] %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", k2, k1, k0);
        fprintf(pipe, "%f\t%f\n", X, Y);
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);

        #ifdef WIN32
            _pclose(pipe);
        #else
            pclose(pipe);
        #endif

        return 0;
    }

    if (n == 3) {
        double k3 = x.data[0][0];
        double k2 = x.data[1][0];
        double k1 = x.data[2][0];
        double k0 = x.data[3][0];
        fprintf(pipe, "plot [5 : 15] [0 : 5] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", k3, k2, k1, k0);
        fprintf(pipe, "%f\t%f\n", X, Y);
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);

        #ifdef WIN32
            _pclose(pipe);
        #else
            pclose(pipe);
        #endif

        return 0;
    }

    return 0;
}