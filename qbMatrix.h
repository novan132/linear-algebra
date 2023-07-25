#ifndef QBMATRIX2_H
#define QBMATRIX2_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

template <class T>
class qbMatrix2 {
public:
    qbMatrix2();
    qbMatrix2(int nRows, int nCols);
    qbMatrix2(int nRows, int nCols, const T* inputData);
    qbMatrix2(const qbMatrix2<T>& inputMatrix);
    qbMatrix2(int nRows, int nCols, const std::vector<T>* inputData);
    ~qbMatrix2();

    bool Resize(int numRows, int numCols);
    void SetToIdentity();

    // access method
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    bool Inverse();

    // overload operator
    bool operator== (const qbMatrix2<T>& rhs);

    bool Compare(const qbMatrix2<T> matrix1, double tolerance);
    bool Separate(qbMatrix2<T>* matrix1, qbMatrix2<T>*, int colNum);

    template <class U> friend qbMatrix2<U> operator+ (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator+ (const U& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator+ (const qbMatrix2<U>& lhs, const U& rhs);

    template <class U> friend qbMatrix2<U> operator- (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator- (const U& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator- (const qbMatrix2<U>& lhs, const U& rhs);

    template <class U> friend qbMatrix2<U> operator* (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator* (const U& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator* (const qbMatrix2<U>& lhs, const U& rhs);

public:
    int Sub2Ind(int row, int col);
    bool IsSquare();
    bool CloseEnough(T f1, T f2);
    void SwapRow(int i, int j);
    void MultAdd(int i, int j, T multFactor);
    void MultRow(int i, T multFactor);
    bool Join(const qbMatrix2<T>& matrix2);
    int FindRowWithMaxElement(int colNumber, int startingRow);
    void PrintMatrix();

private:
    T* m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};

template <class T>
qbMatrix2<T>::qbMatrix2() {
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

template <class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols) {
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i)
        m_matrixData[i] = 0.0;
}

template <class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols, const T* inputData) {
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i)
        m_matrixData[i] = inputData[i];
}

template <class T>
qbMatrix2<T>::qbMatrix2(const qbMatrix2<T>& inputMatrix) {
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i)
        m_matrixData[i] = inputMatrix.m_matrixData[i];
}

template <class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols, const std::vector<T>* inputData) {
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i)
        m_matrixData[i] = inputData->at(i);
}

template <class T>
qbMatrix2<T>::~qbMatrix2() {
    if (m_matrixData != nullptr){
        delete[] m_matrixData;
    } 
}

template <class T>
bool qbMatrix2<T>::Resize(int numRows, int numCols) {
    m_nRows = numRows;
    m_nCols = numCols;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;

    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr)  {
        for (int i = 0; i < m_nElements; ++i)
            m_matrixData[i] = 0.0;
        return true;
    }
    return false;
}

template <class T>
void qbMatrix2<T>::SetToIdentity() {
    if (!IsSquare()) throw std::invalid_argument("cannot form identity from non square matrix");
    for (int row = 0; row < m_nRows; ++row) {
        for (int col = 0; col < m_nCols; ++col) {
            if (col == row) m_matrixData[Sub2Ind(row, col)] = 1.0;
            else m_matrixData[Sub2Ind(row, col)] = 0.0;
        }
    }
}

template <class T>
T qbMatrix2<T>::GetElement(int row, int col) {
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0) {
        return m_matrixData[linearIndex];
    }
    return 0.0;
}

template <class T>
bool qbMatrix2<T>::SetElement(int row, int col, T elementValue) {
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0) {
        m_matrixData[linearIndex] = elementValue;
        return true;
    }
    return false;
}

template <class T>
int qbMatrix2<T>::Sub2Ind(int row, int col) {
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0)) {
        return (row * m_nCols) + col;
    }
    return -1;
}

template <class T>
int qbMatrix2<T>::GetNumRows() {
    return m_nRows;
}

template <class T>
int qbMatrix2<T>::GetNumCols() {
    return m_nCols;
}

template <class T>
bool qbMatrix2<T>::Compare(const qbMatrix2<T>& matrix1, double tolerance) {
   int numRows1 = matrix1.m_nRows;
   int numCols1 = matrix1.m_nCols;
   if ((numRows1 != m_nRows) || (numCols1 != m_nCols)) {
        return false;
   }

   double cumulativeSum = 0.0;
   for (int i = 0; i < m_nElements; ++i) {
        T element1 = matrix1.m_matrixData[i];
        T element2 = m_matrixData[i];
        cumulativeSum += ((element1 - element22) * (element1 - element2));
   }
   double finalValue = sqrt(cumulativeSum / ((numRows1 * numCols1) - 1));

   if (finalValue < tolerance) return true;
   return false;
}

// + operator
template <class T>
qbMatrix2<T> operator+ (const qbMatrix2<T>& lhs, const qbMatrix2<T>& rhs) {
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

template <class T>
qbMatrix2<T> operator+ (const T& lhs, const qbMatrix2<T>& rhs) {
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs + rhs.m_matrixData[i];
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;

}

template <class T>
qbMatrix2<T> operator+ (const qbMatrix2<T>& lhs, const T& rhs) {
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs.m_matrixData[i] + rhs;
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;

}

// - operator
template <class T>
qbMatrix2<T> operator- (const qbMatrix2<T>& lhs, const qbMatrix2<T>& rhs) {
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

template <class T>
qbMatrix2<T> operator- (const T& lhs, const qbMatrix2<T>& rhs) {
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs - rhs.m_matrixData[i];
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

template <class T>
qbMatrix2<T> operator- (const qbMatrix2<T>& lhs, const T& rhs) {
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs.m_matrixData[i] - rhs;
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// * operator
template <class T>
qbMatrix2<T> operator* (const qbMatrix2<T>& lhs, const qbMatrix2<T>& rhs) {
    int r_numRows = rhs.m_nRows;
    int r_numCols = rhs.m_nCols;
    int l_numRows = lhs.m_nRows;
    int l_numCols = lhs.m_nCols;

    if (l_numCols == r_numRows) {
        T* tempResult = new T[lhs.m_nRows * rhs.m_nCols];
        for (int lhsRow = 0; lhsRow < l_numRows; ++lhsRow) {
            for (int rhsCol = 0; rhsCol < r_numCols; ++rhsCol) {
                T elementResult = 0.0;
                for (int lhsCol = 0; lhsCol < l_numCols; ++lhsCol) {
                    int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
                    int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;
                    elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]) ;
                }
                int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
                tempResult[resultLinearIndex] = elementResult;
            }
        }
        qbMatrix2<T> result(l_numRows, r_numCols, tempResult);
        delete[] tempResult;
        return result;
    }
    else {
        qbMatrix2<T> result(1, 1);
        return result;
    }
}

template <class T>
qbMatrix2<T> operator* (const T& lhs, const qbMatrix2<T>& rhs) {
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs * rhs.m_matrixData[i];
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

template <class T>
qbMatrix2<T> operator* (const qbMatrix2<T>& lhs, const T& rhs) {
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i) {
        tempResult[i] = lhs.m_matrixData[i] * rhs;
    }

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// == operator
template <class T>
bool qbMatrix2<T>::operator== (const qbMatrix2<T>& rhs) {
    if ((this->m_nRows != rhs.m_nRows) || (this->m_nCols != rhs.m_nCols)) return false;
    for (int i = 0; i < this->m_nElements; ++i) {
        // if (this->m_matrixData[i] != rhs.m_matrixData) return false;
        if (!CloseEnough(this->m_matrixData[i], rhs.m_matrixData[i])) return false;
    }
    return true;
}

// separate matrix
template <class T>
bool qbMatrix2<T> Separate(qbMatrix2<T>* matrix1, qbMatrix2<T>* matrix2, int colNum) {
    int numRows = m_nRows;
    int numCols1 = colNum;
    int numCol2 = m_nCols - colNum;

    matrix1->Resize(numRows, numCols1);
    matrix2->Resize(numRows, numCols1);

    for (int row = 0; row < m_nRows; ++row) {
        for (int col = 0; col < m_nCols; ++col) {
            if (col < colNum) {
                matrix1->SetElement(row, col, this->GetElement(row, col));
            }
            else {
                matrix2->SetElement(row, col-colNum, this->GetElement(row, col));
            }
        }
    }
    return true;
}

template <class T>
bool qbMatrix2<T>::CloseEnough(T f1, T f2) {
    return fabs(f1 - f2) < 1e-9;
}

#endif
