#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

int numRows; // m rows
int numColumns; // n columns

// function to print out the matrix
void printMatrix(const string& type, vector<int> matrix) {

    int k;

    if (type == "column") {
        k = numRows;
    } else {
        k = numColumns;
    }

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numColumns; j++) {
            if (k == numRows) {
                cout << matrix[(j * numRows) + i] << " ";
            } else {
                cout << matrix[(i * numColumns) + j] << " ";
            }
        }
        cout << endl;
    }
}

// function to create matrix x via a .txt file
vector<int> makeMatrix() {

    ifstream matrixFile;

    matrixFile.open("./matrix.txt");

    if (matrixFile.is_open()) {
        matrixFile >> numRows >> numColumns;
    }

    vector<int> matrix(numRows * numColumns);

    // user-input for matrix
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numColumns; j++) {
            matrixFile >> matrix[((j * numRows) + i)];
        }
    }

    matrixFile.close();

    cout << endl << "[MAIN MATRIX]" << endl;

    // prints the matrix
    printMatrix("column", matrix);

    return matrix;
}

// function to get the "transpose" of matrix x
vector<int> transpose(vector<int> matrix) {

    // creates an empty matrix for the transpose
    vector<int> transpose(numRows * numColumns);

    // switches around the row and column value for the transpose matrix
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numColumns; j++) {
            transpose[(i * numColumns) + j] = matrix[(j * numRows) + i];
        }
    }

    return transpose;
}

// function to perform bit-reversal on a given matrix
vector<int> bitReversal(vector<int> matrix) {

    vector<int> bitReverse(numRows * numColumns);

    int p = (int) log2(matrix.size());

    bitReverse[1] = (int) pow(2, p - 1);

    if (p > 1) {
        bitReverse[2] = bitReverse[1] / 2;
        bitReverse[3] = 3 * bitReverse[2];
        int nd0 = 3;

        for (int k = 3; k <= p; k++) {
            int nd = (int) pow(2, k) - 1;
            int nd1 = nd + 1;
            bitReverse[nd1 - 1] = bitReverse[nd0] + (int) pow(2, p - k);

            for (int l = 1; l <= nd0; l++) {
                bitReverse[nd1 - l - 1] = bitReverse[nd1 - 1] - bitReverse[l];
            }
            nd0 = nd;
        }
    }

    for (int i = 0; i < numRows * numColumns; i++) {
        if (bitReverse[i] != 0) {
            swap(matrix[i], matrix[bitReverse[i]]);
            bitReverse[bitReverse[i]] = 0;
            bitReverse[i] = 0;
        }
    }

    return matrix;
}

// function to perform either row or column FWHT method on a matrix
vector<int> FWHT(const string& type, vector<int> matrix, bool bitRev, bool bitRevTwice) {

    int k;

    if (type == "column") {
        k = numRows;
    } else {
        k = numColumns;
    }

    if ((k & k - 1) != 0 || k == 0) {
        cout << endl << "FWHT cannot be performed on this matrix. Please input a matrix where there are 2^m rows or columns.";
        exit(-1);
    }

    cout << endl << "[FWHT ANSWER]" << endl;

    if (bitRev && !bitRevTwice) {
        vector<int> newMatrix(numRows * numColumns);
        int p = 0;

        while (p < numRows * numColumns) {
            for (int j = 0; j <= k; j = j + 4) {
                for (int i = 0; i <= 1; i++) {
                    matrix[i + j] = matrix[i + j] + matrix[i + j + 2];
                    matrix[i + j + 2] = matrix[i + j] - 2 * matrix[i + j + 2];
                }
            }
            for (int j = 0; j < 2; j = j + 1) {
                for (int l = 0; l < k; l = l + 2) {
                    int u = matrix[j + l];
                    int v = matrix[j + l + 4];
                    newMatrix[p] = u + v;
                    newMatrix[p + 1] = u - v;
                    p += 2;
                }
            }
        }
        matrix = newMatrix;
    } else if (bitRevTwice) {
        int p = 2;
        while (p < k) {
            for (int j = 0; j <= k; j = j + 4) {
                for (int i = 0; i <= 1; i++) {
                    int u = matrix[i + j];
                    int v = matrix[i + j + p];
                    matrix[i + j] = u + v;
                    matrix[i + j + p] = u - v;
                }
            }
            p *= 2;
            for (int l = 0; l < k; l++) {
                int u = matrix[l];
                int v = matrix[l + p];
                matrix[l] = u + v;
                matrix[l + p] = u - v;
            }
        }
        matrix = bitReversal(matrix);
    } else {
        int h = 2;
        while (h <= k) {
            int hf = floor(h/2);
            for (int i = 0; i < numRows * numColumns; i = i + h) {
                for (int j = 0; j < hf; j++) {
                    int u = matrix[i + j];
                    int v = matrix[i + j + hf];
                    matrix[i + j] = u + v;
                    matrix[i + j + hf] = u - v;
                }
            }
            h *= 2;
        }
    }

    printMatrix(type, matrix);

    return matrix;
}

// function to measure the allotted time it takes for FWHT to perform
void timeFWHT(bool bitRev, bool bitRevTwice) {

    auto start = high_resolution_clock::now();
    vector<int> x = makeMatrix();

    if (numRows >= numColumns) {
        if (bitRev && !bitRevTwice) {
            FWHT("column", bitReversal(x), bitRev, bitRevTwice);
        } else if (bitRevTwice) {
            FWHT("column", bitReversal(x), bitRev, bitRevTwice);
        } else {
            FWHT("column", x, bitRev, bitRevTwice);
        }
    } else {
        if (bitRev && !bitRevTwice) {
            FWHT("row", bitReversal(transpose(x)), bitRev, bitRevTwice);
        } else if (bitRevTwice) {
            FWHT("row", bitReversal(transpose(x)), bitRev, bitRevTwice);
        } else {
            FWHT("row", transpose(x), bitRev, bitRevTwice);
        }
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    if (bitRev && !bitRevTwice) {
        cout << endl << "Allotted time with bit reversal: " << duration.count() << " microseconds" << endl;
    } else if (bitRevTwice) {
        cout << endl << "Allotted time with bit reversal twice: " << duration.count() << " microseconds" << endl;
    } else {
        cout << endl << "Allotted time without bit reversal: " << duration.count() << " microseconds" << endl;
    }
}

// main method
int main() {

    timeFWHT(false, false);
    timeFWHT(true, false);
    timeFWHT(true, true);

    return 0;
}