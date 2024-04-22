#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
using namespace std;

int Gaus(vector<vector<int>>& matrix) {
    int n = matrix.size();
    if (n == 0) {
        throw invalid_argument("Matrix must be non-empty.");
    }
    int det = 1;
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(matrix[j][i]) > abs(matrix[pivot][i])) {
                pivot = j;
            }
        }
        swap(matrix[i], matrix[pivot]);
        if (matrix[i][i] == 0) {
            return 0;
        }
        det *= matrix[i][i];
        for (int j = i + 1; j < n; ++j) {
            int factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k < n; ++k) {
                matrix[j][k] -= matrix[i][k] * factor;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (matrix[i][i] == 0) {
            return 0;
        }
        det *= matrix[i][i];
    }
    return det;
}

void readDataFromFile(const string& filename, vector<int>& V, vector<int>& C, vector<int>& R) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Не удалось открыть файл " << filename << endl;
        return;
    }
    string line;
    int value;
    if (getline(file, line)) {
        istringstream iss(line);
        while (iss >> value) {
            V.push_back(value);
        }
    }
    if (getline(file, line)) {
        istringstream iss(line);
        while (iss >> value) {
            C.push_back(value);
        }
    }
    if (getline(file, line)) {
        istringstream iss(line);
        while (iss >> value) {
            R.push_back(value);
        }
    }
    file.close();
}
void printDenseMatrixFromCSR(const vector<int>& values, const vector<int>& column_indices, const vector<int>& row_indices, bool printDeterminant = false) {
    size_t rows = row_indices.size() - 1;
    size_t cols = 0;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = row_indices[i]; j < row_indices[i + 1]; ++j) {
            if (column_indices[j] > cols) {
                cols = column_indices[j];
            }
        }
    }
    cols++;
    vector<vector<int>> denseMatrix(rows, vector<int>(cols, 0));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = row_indices[i]; j < row_indices[i + 1]; ++j) {
            denseMatrix[i][column_indices[j]] = values[j];
        }
    }
    for (const auto& row : denseMatrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << std::endl;
    }
    if (printDeterminant) {
        auto start2 = chrono::high_resolution_clock::now();
        cout << "Determinant: " << Gaus(denseMatrix) << std::endl;
        auto end2 = chrono::high_resolution_clock::now();
        auto duration2 = chrono::duration_cast<chrono::microseconds>(end2 - start2).count();
        cout << "Время выполения алгоритма: " << duration2 << " microseconds" << endl;
    }
}

vector<vector<int>> generateRandomMatrix(int n, double zeroPercentage) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, 100);
    vector<vector<int>> matrix(n, vector<int>(n, 0));
    int numZeros = static_cast<int>(n * n * zeroPercentage / 100.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (numZeros > 0 && dis(gen) % 100 < zeroPercentage) {
                matrix[i][j] = 0;
                numZeros--; 
            }
            else {
                matrix[i][j] = dis(gen);
            }
        }
    }
    return matrix;
}

void convertToSparse(const vector<vector<int>>& matrix, vector<int>& values, vector<int>& rowIndices, vector<int>& colIndices) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size(); 
    values.clear();
    rowIndices.clear();
    colIndices.clear();
    rowIndices.push_back(0);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (matrix[i][j] != 0) { 
                values.push_back(matrix[i][j]); 
                colIndices.push_back(j); 
            }
        }
        rowIndices.push_back(values.size()); 
    }
}

void writeVectorsToFile(const vector<int>& Values, const vector<int>& RowIndices, const vector<int>& ColIndices) {
    ofstream outfile("matrix.txt");
    copy(Values.begin(), Values.end(), ostream_iterator<int>(outfile, " "));
    outfile << "\n"; 
    copy(RowIndices.begin(), RowIndices.end(), ostream_iterator<int>(outfile, " "));
    outfile << "\n"; 
    copy(ColIndices.begin(), ColIndices.end(), ostream_iterator<int>(outfile, " "));
    outfile << "\n"; 
    outfile.close();
}

int main() {
    setlocale(LC_ALL, "Russian");
    vector<int> V, C, R;
    int k;
    cout << "Выберите способ запуск программы:" << endl;
    cout << "Если вы хотите создать разреженную матрицу используя входные данные и найти ее определитель. Напишите: 1" << endl;
    cout << "Если вы хотите создать разреженную матрицу используя случайную генерацию, найти определитель и записать матрицу в файл. Напишите: 2" << endl;
    cin >> k;
    cout << endl;
    if (k == 1) {
        readDataFromFile("b.txt", V, C, R);
        cout << "Vector V: ";
        for (int i = 0; i < V.size(); i++) {
            cout << V.at(i) << ' ';
        }
        cout << endl;
        cout << "Vector C: ";
        for (int i = 0; i < C.size(); i++) {
            cout << C.at(i) << ' ';
        }
        cout << endl;
        cout << "Vector R: ";
        for (int i = 0; i < R.size(); i++) {
            cout << R.at(i) << ' ';
        }
        cout << endl;
        printDenseMatrixFromCSR(V, C, R, true);
        cout << endl;
    }
    if (k == 2) {
        cout << "N x N ( enter N)" << endl;
        int n;
        double zeroPercent;
        cin >> n;
        cout << "Enter the percentage of zeros (0-100): ";
        cin >> zeroPercent;
        vector<vector<int>> matrix = generateRandomMatrix(n, zeroPercent);
        vector<int> Values, RowIndices, ColIndices;
        convertToSparse(matrix, Values, RowIndices, ColIndices);
        cout << endl;
        if (n < 10) {
            for (const auto& row : matrix) {
                for (const auto& element : row) {
                    cout << element << " ";
                }
                cout << endl; 
            }
        }
        cout << endl;
        int dd;
        auto start3 = chrono::high_resolution_clock::now();
        dd = Gaus(matrix);
        auto end3 = chrono::high_resolution_clock::now();
        cout << "D = " << dd << endl;
        auto duration3 = chrono::duration_cast<chrono::microseconds>(end3 - start3).count();
        cout << "Время выполения алгоритма: " << duration3 << " microseconds" << endl;
        writeVectorsToFile(Values, ColIndices, RowIndices);
        cout << "Матрица была записана в файл matrix.txt" << endl;
    }
    return 0;
}