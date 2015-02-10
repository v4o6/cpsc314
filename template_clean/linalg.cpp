
#include <iostream>
#include <cmath>

#include "linalg.hpp"


Vector::Vector(double x, double y, double z, double w) {
    data_[0] = x;
    data_[1] = y;
    data_[2] = z;
    data_[3] = w;
}


Vector::Vector(Vector const &other) {
    for (int i = 0; i < 4; i++) {
        data_[i] = other.data_[i];
    }
}


void Vector::print() const {
    std::cout << '(';
    for (int i = 0; i < 4; i++) {
        std::cout << (i ? "," : "") << data_[i];
    }
    std::cout << ')';
}


double Vector::length() const {
    double sum = 0;
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        sum += data_[i] * data_[i];
    }
    return sqrt(sum);
}


double Vector::length2() const {
    double sum = 0;
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        sum += data_[i] * data_[i];
    }
    return sum;
}


void Vector::normalize() {
    double len = length();
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        data_[i] /= len;
    }
}


Vector Vector::normalized() const {
    Vector v(*this);
    v.normalize();
    return v;
}


Vector Vector::operator+(Vector const &other) const {
    Vector result;
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        result[i] = data_[i] + other.data_[i];
    }
    return result;
}


Vector Vector::operator-(Vector const &other) const {
    Vector result;
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        result[i] = data_[i] - other.data_[i];
    }
    return result;
}


Vector Vector::operator*(double value) const {
    Vector result;
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        result[i] = data_[i] * value;
    }
    return result;
}


Vector Vector::operator/(double value) const {
    Vector result;
    // Ignoring W.
	for (int i = 0; i < 3; i++) {
        result[i] = data_[i] / value;
    }
    return result;
}


Vector Vector::operator*(Vector const &other) const {
    Vector result;
    // Ignoring W.
    for (int i = 0; i < 3; i++) {
        result[i] = data_[i] * other.data_[i];
    }
    return result;
}


double Vector::dot(Vector const &other, bool homogeneous) const {
    double result = 0;
    for (int i = 0; i < (homogeneous ? 4 : 3); i++) {
        result += data_[i] * other.data_[i];
    }
    return result;
}


Vector Vector::cross(Vector const &other, double w) const {
    return Vector(
         data_[1] * other.data_[2] - other.data_[1] * data_[2],
        -data_[0] * other.data_[2] + other.data_[0] * data_[2],
         data_[0] * other.data_[1] - other.data_[0] * data_[1],
         w
    );
}

Vector Vector::operator*(Matrix const &other) const {
	return Vector(
		data_[0] * (other[0][0] + other[1][0] + other[2][0] + other[3][0]),
		data_[1] * (other[0][1] + other[1][1] + other[2][1] + other[3][1]),
		data_[2] * (other[0][2] + other[1][2] + other[2][2] + other[3][2]),
		data_[3] * (other[0][3] + other[1][3] + other[2][3] + other[3][3])
	);
}


Matrix::Matrix(Matrix const &other) {
    for (int i = 0; i < 4; i++) {
        data_[i] = other.data_[i];
    }
}


Matrix::Matrix(Vector const &a, Vector const &b, Vector const &c, Vector const &d) {
    data_[0] = a;
    data_[1] = b;
    data_[2] = c;
    data_[3] = d;
}


Matrix::Matrix(double a, double b, double c, double d,
         double e, double f, double g, double h,
         double i, double j, double k, double l,
         double m, double n, double o, double p) {
     data_[0] = Vector(a, b, c, d);
     data_[1] = Vector(e, f, g, h);
     data_[2] = Vector(i, j, k, l);
     data_[3] = Vector(m, n, o, p);
}


Matrix Matrix::identity() {
    return Matrix(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    );
}


Matrix Matrix::rotation(double angle, Vector const &axis) {

    double c_angle = cos(angle);
    double v_angle = 1.0 - c_angle;
    double s_angle = sin(angle);

    return Matrix(
        axis[0] * axis[0] * v_angle + c_angle,
        axis[1] * axis[0] * v_angle - axis[2] * s_angle,
        axis[2] * axis[0] * v_angle + axis[1] * s_angle,
        0,
        axis[0] * axis[1] * v_angle + axis[2] * s_angle,
        axis[1] * axis[1] * v_angle + c_angle,
        axis[2] * axis[1] * v_angle - axis[0] * s_angle,
        0,
        axis[0] * axis[2] * v_angle - axis[1] * s_angle,
        axis[1] * axis[2] * v_angle + axis[0] * s_angle,
        axis[2] * axis[2] * v_angle + c_angle,
        0,
        0, 0, 0, 1
    );
}

void Matrix::print(bool pretty) const {
    std::cout << '(';
    for (int r = 0; r < 4; r++) {
        std::cout << (r ? ";" : "") << (r && pretty ? "\n " : "");
        for (int c = 0; c < 4; c++) {
            std::cout << (c ? "," : "") << data_[r][c];
        }
    }
    std::cout << ')';
}


Matrix Matrix::transpose() const {
    Matrix result;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            result[r][c] = data_[c][r];
        }
    }
    return result;
}

Matrix Matrix::inverse() const {
	/* 4x4 version */
	//Vector data0 = Vector(	data_[1][1] * data_[2][2] * data_[3][3] -
	//						data_[1][1] * data_[2][3] * data_[3][2] -
	//						data_[2][1] * data_[1][2] * data_[3][3] +
	//						data_[2][1] * data_[1][3] * data_[3][2] +
	//						data_[3][1] * data_[1][2] * data_[2][3] -
	//						data_[3][1] * data_[1][3] * data_[2][2],

	//						-data_[0][1] * data_[2][2] * data_[3][3] +
	//						data_[0][1] * data_[2][3] * data_[3][2] +
	//						data_[2][1] * data_[0][2] * data_[3][3] -
	//						data_[2][1] * data_[0][3] * data_[3][2] -
	//						data_[3][1] * data_[0][2] * data_[2][3] + 
	//						data_[3][1] * data_[0][3] * data_[2][2],

	//						data_[0][1] * data_[1][2] * data_[3][3] - 
	//						data_[0][1] * data_[1][3] * data_[3][2] - 
	//						data_[1][1] * data_[0][2] * data_[3][3] + 
	//						data_[1][1] * data_[0][3] * data_[3][2] +
	//						data_[3][1] * data_[0][2] * data_[1][3] - 
	//						data_[3][1] * data_[0][3] * data_[1][2],

	//						-data_[0][1] * data_[2][2] * data_[2][3] +
	//						data_[0][1] * data_[2][3] * data_[2][2] +
	//						data_[1][1] * data_[0][2] * data_[2][3] -
	//						data_[1][1] * data_[0][3] * data_[2][2] -
	//						data_[2][1] * data_[0][2] * data_[1][3] +
	//						data_[2][3] * data_[0][3] * data_[1][2]);

	//Vector data1 = Vector(	-data_[1][0] * data_[2][2] * data_[3][3] +
	//						data_[1][0] * data_[2][3] * data_[3][2] +
	//						data_[2][0] * data_[1][2] * data_[3][3] - 
	//						data_[2][0] * data_[1][3] * data_[3][2] -
	//						data_[3][0] * data_[1][2] * data_[2][3] + 
	//						data_[3][0] * data_[1][3] * data_[2][2],

	//						data_[0][0] * data_[2][2] * data_[3][3] -
	//						data_[0][0] * data_[2][3] * data_[3][2] -
	//						data_[2][0] * data_[0][2] * data_[3][3] +
	//						data_[2][0] * data_[0][3] * data_[3][2] +
	//						data_[3][0] * data_[0][2] * data_[2][3] -
	//						data_[3][0] * data_[0][3] * data_[2][2],

	//						-data_[0][0] * data_[1][2] * data_[3][3] +
	//						data_[0][0] * data_[1][3] * data_[3][2] +
	//						data_[1][0] * data_[0][2] * data_[3][3] -
	//						data_[1][0] * data_[0][3] * data_[3][2] -
	//						data_[3][0] * data_[0][2] * data_[1][3] + 
	//						data_[3][0] * data_[0][3] * data_[1][2],

	//						data_[0][0] * data_[1][2] * data_[2][3] - 
	//						data_[0][0] * data_[1][3] * data_[2][2] - 
	//						data_[1][0] * data_[0][2] * data_[2][3] + 
	//						data_[1][0] * data_[0][3] * data_[2][2] +
	//						data_[2][0] * data_[0][2] * data_[1][3] - 
	//						data_[2][0] * data_[0][3] * data_[1][2]);

	//Vector data2 = Vector(	data_[1][0] * data_[2][1] * data_[3][3] -
	//						data_[1][0] * data_[2][3] * data_[3][1] -
	//						data_[2][0] * data_[1][1] * data_[3][3] +
	//						data_[2][0] * data_[1][3] * data_[3][1] +
	//						data_[3][0] * data_[1][1] * data_[2][3] -
	//						data_[3][0] * data_[1][3] * data_[2][1],

	//						-data_[0][0] * data_[2][1] * data_[3][3] +
	//						data_[0][0] * data_[2][3] * data_[3][1] +
	//						data_[2][0] * data_[0][1] * data_[3][3] -
	//						data_[2][0] * data_[0][3] * data_[3][1] -
	//						data_[3][0] * data_[0][1] * data_[2][3] + 
	//						data_[3][0] * data_[0][3] * data_[2][1],

	//						data_[0][0] * data_[1][1] * data_[3][3] - 
	//						data_[0][0] * data_[1][3] * data_[3][1] - 
	//						data_[1][0] * data_[0][1] * data_[3][3] + 
	//						data_[1][0] * data_[0][3] * data_[3][1] +
	//						data_[3][0] * data_[0][1] * data_[1][3] - 
	//						data_[3][0] * data_[0][3] * data_[1][1],

	//						-data_[0][0] * data_[2][1] * data_[2][3] +
	//						data_[0][0] * data_[2][3] * data_[2][1] +
	//						data_[1][0] * data_[0][1] * data_[2][3] -
	//						data_[1][0] * data_[0][3] * data_[2][1] -
	//						data_[2][0] * data_[0][1] * data_[1][3] +
	//						data_[2][0] * data_[0][3] * data_[1][1]);

	//Vector data3 = Vector(	-data_[1][0] * data_[2][1] * data_[3][2] +
	//						data_[1][0] * data_[2][2] * data_[3][1] +
	//						data_[2][0] * data_[1][1] * data_[3][2] - 
	//						data_[2][0] * data_[1][2] * data_[3][1] -
	//						data_[3][0] * data_[1][1] * data_[2][2] + 
	//						data_[3][0] * data_[1][2] * data_[2][1],

	//						data_[0][0] * data_[2][1] * data_[3][2] -
	//						data_[0][0] * data_[2][2] * data_[3][1] -
	//						data_[2][0] * data_[0][1] * data_[3][2] +
	//						data_[2][0] * data_[0][2] * data_[3][1] +
	//						data_[3][0] * data_[0][1] * data_[2][2] -
	//						data_[3][0] * data_[0][2] * data_[2][1],

	//						-data_[0][0] * data_[1][1] * data_[3][2] +
	//						data_[0][0] * data_[1][2] * data_[3][1] +
	//						data_[1][0] * data_[0][1] * data_[3][2] -
	//						data_[1][0] * data_[0][2] * data_[3][1] -
	//						data_[3][0] * data_[0][1] * data_[1][2] + 
	//						data_[3][0] * data_[0][2] * data_[1][1],

	//						data_[0][0] * data_[1][1] * data_[2][2] - 
	//						data_[0][0] * data_[1][2] * data_[2][1] - 
	//						data_[1][0] * data_[0][1] * data_[2][2] + 
	//						data_[1][0] * data_[0][2] * data_[2][1] +
	//						data_[2][0] * data_[0][1] * data_[1][2] - 
	//						data_[2][0] * data_[0][2] * data_[1][1]);

	//double det = data_[0][0] * data0[0] - data_[0][1] * data1[0] + data_[0][2] * data2[0] - data_[0][3] * data3[0];
	//det = 1.0/det;
	//return Matrix(data0 * det, data1 * det, data2 * det, data3 * det);

	/* 3x3 version; ignoring W. */
	Vector data0 = Vector(	data_[1][1] * data_[2][3] - data_[1][2] * data_[2][1],
							data_[0][2] * data_[2][1] - data_[0][1] * data_[2][3],
							data_[0][1] * data_[1][2] - data_[0][2] * data_[1][1],
							0.0		);
	Vector data1 = Vector(	data_[1][2] * data_[2][0] - data_[1][0] * data_[2][3],
							data_[0][0] * data_[2][3] - data_[0][2] * data_[2][0],
							data_[0][2] * data_[1][0] - data_[0][0] * data_[1][2],
							0.0		);
	Vector data2 = Vector(	data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0],
							data_[0][1] * data_[2][0] - data_[0][0] * data_[2][1],
							data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0],
							0.0		);
	Vector data3 = Vector(	0.0, 0.0, 0.0, 1.0	);

	double det = data_[0][0] * data0[0] + data_[0][1] * data1[0] + data_[0][2] * data2[0];
	det = 1.0/det;
	return Matrix(data0 * det, data1 * det, data2 * det, data3);
}

Vector Matrix::operator*(Vector const &other) const {
    return Vector(
        other.dot(data_[0], true),
        other.dot(data_[1], true),
        other.dot(data_[2], true),
        other.dot(data_[3], true)
    );
}


Matrix Matrix::operator*(Matrix const &other) const {
    Matrix transposed = other.transpose();
    Matrix result;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            result[r][c] = data_[r].dot(transposed.data_[c], true);
        }
    }
    return result;
}


void Matrix::operator*=(Matrix const &other) {
    Matrix copy = (*this) * other;
    for (int i = 0; i < 4; i++) {
        data_[i] = copy.data_[i];
    }
}






