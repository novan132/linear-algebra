#ifndef QBMATRIX2_H
#define QBMATRIX2_H

template <class T>
class qbMatrix2 {
public:
    qbMatrix2();
    qbMatrix2(int nRows, int nCols);
    qbMatrix2(int nRows, int nCols, const T* inputData);
    qbMatrix2(const qbMatrix2<T>& inputMatrix);
    ~qbMatrix2();

    bool resize(int numRows, int numCols);

    // access method
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    // overload operator
    bool operator== (const qbMatrix2<T>& rhs);

    template <class U> friend qbMatrix2<U> operator+ (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator+ (const U& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator+ (const qbMatrix2<U>& lhs, const U& rhs);

    template <class U> friend qbMatrix2<U> operator- (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator- (const U& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator- (const qbMatrix2<U>& lhs, const U& rhs);

    template <class U> friend qbMatrix2<U> operator* (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator* (const U& lhs, const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator* (const qbMatrix2<U>& lhs, const U& rhs);

private:
    int Sub2Ind(int row, int col);

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
qbMatrix2<T>::~qbMatrix2() {
    if (m_matrixData != nullptr){
        delete[] m_matrixData;
    } 
}

template <class T>
bool qbMatrix2<T>::resize(int numRows, int numCols) {
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
        if (this->m_matrixData[i] != rhs.m_matrixData) return false;
    }
    return true;
}

#endif
