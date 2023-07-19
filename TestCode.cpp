#include <iostream>
#include <string>
#include <math.h>

#include "qbMatrix.h"

using namespace std;

template <class T>
void PrintMatrix(qbMatrix2<T> matrix) {
    int nRows = matrix.GetNumRows();
    int nCols = matrix.GetNumCols();
    for (int row = 0; row < nRows; ++row) {
        for (int col = 0; col < nCols; ++col) {
            cout << matrix.GetElement(row, col) << " ";
        }
        cout << endl;
    }
}

int main() {
    cout << "code to test qbMatrix2" << endl;
    double simpleData[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    qbMatrix2<double> testMatrix(3, 4, simpleData);

    cout << endl << "------------------------" << endl;
    cout << "3 x 4 matrix test" << endl;
    PrintMatrix(testMatrix);

    cout << endl << "------------------------"  << endl;
    cout << "test retrieval" << endl;
    cout << "element (0, 0): " << testMatrix.GetElement(0, 0) << endl;
    cout << "element (1, 0): " << testMatrix.GetElement(1, 0) << endl;
    cout << "element (2, 0): " << testMatrix.GetElement(2, 0) << endl;
    cout << "element (0, 1): " << testMatrix.GetElement(0, 1) << endl;
    cout << "element (1, 1): " << testMatrix.GetElement(1, 1) << endl;
    cout << "element (2, 1): " << testMatrix.GetElement(2, 1) << endl;
    cout << "element (5, 5): " << testMatrix.GetElement(5, 5) << endl;

    cout << endl << "------------------------"  << endl;
    cout << "test multiplication" << endl;
    double simpleData2[12] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    qbMatrix2<double> testMatrix2(4, 3, simpleData2);
    cout << "4 x 3 matrix (testMatrix2)" << endl;
    PrintMatrix(testMatrix2);
    cout << "multiplication (testMatrix * testMatrix2) result:" << endl;
    PrintMatrix(testMatrix * testMatrix2);

    cout << endl << "------------------------"  << endl;
    cout << "testMatrix2 * testMatrix: " << endl;
    PrintMatrix(testMatrix2 * testMatrix);

    cout << endl << "------------------------"  << endl;
    cout << "test multiplication column vector by matrix" << endl;
    double columnData[3] = {1.5, 2.5, 3.5};
    double squareData[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    qbMatrix2<double> testColumn(3, 1, columnData);
    qbMatrix2<double> squareMatrix(3, 3, squareData);
    cout << "column vector: " << endl;
    PrintMatrix(testColumn);
    cout << "square matrix: " << endl;
    PrintMatrix(squareMatrix);
    cout << "column vector * square matrix:" << endl;
    PrintMatrix(testColumn * squareMatrix);
    cout << "square matrix * column vector:" << endl;
    PrintMatrix(squareMatrix * testColumn);
    cout << "square matrix + 1.0: " << endl;
    PrintMatrix(squareMatrix + 1.0);
    cout << "(square matrix + 1.0) * columVector:" << endl;
    PrintMatrix((squareMatrix + 1.0) * testColumn);
    return 0;
}
