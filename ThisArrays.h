#pragma once
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdexcept>

class ThisArrays {
private:
    double* data;
    int rows, cols;

public:
    ThisArrays() : rows(0), cols(0), data(nullptr) {}
    ThisArrays(double* data,int rows, int cols):data(data),rows(rows),cols(cols){}
    ThisArrays(int rows, int cols) : rows(rows), cols(cols) {
        data = new double[rows * cols];
        memset(data, 0, rows * cols * sizeof(double));
    }
    ThisArrays(int rows) : ThisArrays(rows, 1) {
    }

    ThisArrays(const ThisArrays& other) : rows(other.rows), cols(other.cols) {
        data = new double[rows * cols];
        memcpy(data, other.data, rows * cols * sizeof(double));
    }

    ~ThisArrays() {
        delete[] data;
    }

    double& operator()(int x) {
        if (x < 0 || x >= rows) {
            throw std::out_of_range("Index out of range");
        }
        return data[x];
    }

    double& operator()(int x, int y) {
        if (x < 0 || x >= rows || y < 0 || y >= cols) {
            throw std::out_of_range("Index out of range");
        }
        return data[x * cols + y];
    }



    double& operator[](int x) {
        if (x < 0 || x >= rows * cols) {
            throw std::out_of_range("Index out of range");
        }
        return data[x];
    }

    double operator[](int x) const {
        if (x < 0 || x >= rows * cols) {
            throw std::out_of_range("Index out of range");
        }
        return data[x];
    }

    ThisArrays& operator=(const ThisArrays& other) {
        if (this != &other) {
            if (rows != other.rows || cols != other.cols) {
                delete[] data;
                rows = other.rows;
                cols = other.cols;
                data = new double[rows * cols];
            }
            memcpy(data, other.data, rows * cols * sizeof(double));
        }
        return *this;
    }

    ThisArrays operator+(const ThisArrays& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Array dimensions do not match");
        }
        ThisArrays result(rows, cols);
        for (int i = 0; i < rows * cols; ++i) {
            result[i] = data[i] + other[i];
        }
        return result;
    }



    ThisArrays operator-() const {
        ThisArrays result(rows, cols);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = -data[i];
        }
        return result;
    }
    ThisArrays& operator*(double factor)  {
       
        return *this;
    }

    
    ThisArrays operator-(const ThisArrays& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        ThisArrays result(rows, cols);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
    }

    ThisArrays& operator+=(const ThisArrays& other) {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Array dimensions do not match");
        }
        for (int i = 0; i < rows * cols; ++i) {
            data[i] += other[i];
        }
        return *this;
    }

    ThisArrays& operator*=(double factor) {
        for (int i = 0; i < rows * cols; ++i) {
            data[i] *= factor;
        }
        return *this;
    }

    ThisArrays col(int j) {
        return  ThisArrays(data+j*rows,rows,1);
    }


    double at(int x) const {
        if (x < 0 || x >= rows * cols) {
            throw std::out_of_range("Index out of range");
        }
        return data[x];
    }

    void save(const char* fileName) const {
        std::ofstream outFile(fileName, std::ios::out | std::ios::binary);
        outFile.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        outFile.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        outFile.write(reinterpret_cast<const char*>(data), rows * cols * sizeof(double));
        outFile.close();
    }
    friend ThisArrays abs(ThisArrays A);
    friend double accu(ThisArrays A);
    friend double max(ThisArrays A);
    friend ThisArrays linspace(double start, double end, int num);
    ThisArrays operator*(double scalar) const;
    friend ThisArrays operator*(double scalar, const ThisArrays& arr);
};
inline ThisArrays ThisArrays::operator*(double scalar) const {
    ThisArrays result(rows, cols);
    for (int i = 0; i < rows * cols; i++) {
        result.data[i] = data[i] * scalar;
    }
    return result;
}

inline ThisArrays operator*(double scalar, const ThisArrays& arr) {
    return arr * scalar;
}

inline ThisArrays abs(ThisArrays A) {
    ThisArrays result(A.rows, A.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            result(i, j) = std::abs(A(i, j));
        }
    }
    return result;
}

inline double accu(ThisArrays A) {
    double result=0;
    for (int i = 0; i < A.rows*A.cols; ++i) {
        
        result += A[i];
        
    }
    return result;
}
inline double max(ThisArrays A) {
    double result=A[0];
    for (int i = 0; i < A.rows*A.cols; ++i) {
        if (A[i]>result)
        {
            result = A[i];
        }
    }
    return result;
}



inline ThisArrays linspace(double start, double end, int num) {
    auto result = ThisArrays(num);
    double step = (end - start) / num;
    for (int i = 0; i < num; i++) {
        result(i) = start + i * step;
    }
    return result;
}