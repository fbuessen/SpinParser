/**
 * @file Geometry.hpp
 * @author Finn Lasse Buessen
 * @brief Implementation of three- and four-dimensional vectors and matrices for geometric operations. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <cstring>
#include <cmath>

namespace geometry
{
	#pragma region vector
	/**
	 * @brief Three-dimensional vector. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 */
	template <typename T> struct Vec3
	{
		/**
		 * @brief Construct a new Vec3 object with uninitialized values. 
		 */
		Vec3() {};

		/**
		 * @brief Construct a new Vec3 object and initialize all numbers with the same value. 
		 * 
		 * @param xyz Value to initialize. 
		 */
		Vec3(const T &xyz) : x(xyz), y(xyz), z(xyz) {};

		/**
		 * @brief Construct a new Vec3 object and initialize numbers.  
		 * 
		 * @param x Value to initialize x component. 
		 * @param y Value to initialize y component. 
		 * @param z Value to initialize z component. 
		 */
		Vec3(const T &x, const T &y, const T &z) : x(x), y(y), z(z) {};

		T x; ///< x Component. 
		T y; ///< y Component. 
		T z; ///< z Component. 

		/**
		 * @brief Compute the Euclidean norm of a Vec3 object. 
		 * 
		 * @return T Value of the norm. 
		 */
		T norm() const
		{
			return std::sqrt(x * x + y * y + z * z);
		}

		/**
		 * @brief Normalize the Vec3 object with respect to its Euclidean norm. 
		 * 
		 * @return Vec3& Reference to self. 
		 */
		Vec3 &normalize()
		{
			T n = norm();
			x /= n;
			y /= n;
			z /= n;
			return *this;
		}
	};

	/**
	 * @brief Compute the dot product of two Vec3 objects. 
	 * 
	 * @tparam T Underlying fundamental type of vectors. 
	 * @param v Left operand. 
	 * @param w Right operand. 
	 * @return T Dot product of v and w. 
	 */
	template <typename T> T dot(const Vec3<T> &v, const Vec3<T> &w)
	{
		return v.x * w.x + v.y * w.y + v.z * w.z;
	}

	/**
	 * @brief Compute the cross producdt of two Vec3 objects. 
	 * 
	 * @tparam T Underlying fundamental type of vectors. 
	 * @param v Left operand. 
	 * @param w Right operand. 
	 * @return Vec3<T> Cross product of v and w. 
	 */
	template <typename T> Vec3<T> cross(const Vec3<T> &v, const Vec3<T> &w)
	{
		return Vec3<T>(v.y * w.z - v.z * w.y, v.z * w.x - v.x * w.z, v.x * w.y - v.y * w.x);
	}

	/**
	 * @brief Scalar left-multiplication of a Vec3 object. 
	 * 
	 * @tparam T Underlying fundamental type of the vector. 
	 * @param lhs Left operand. 
	 * @param rhs Right operand. 
	 * @return Vec3<T> Product lhs * rhs. 
	 */
	template <typename T> Vec3<T> operator*(const T lhs, const Vec3<T> &rhs)
	{
		return Vec3<T>(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
	}

	/**
	 * @brief Addition of two Vec3 objects. 
	 * 
	 * @tparam T Underlying fundamental type of the vector. 
	 * @param lhs Left operand.
	 * @param rhs Right operand. 
	 * @return Vec3<T> Sum lhs + rhs. 
	 */
	template <typename T> Vec3<T> operator+(const Vec3<T> &lhs, const Vec3<T> &rhs)
	{
		return Vec3<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

	/**
	 * @brief Subtraction of two Vec3 objects. 
	 * 
	 * @tparam T Underlying fundamental type of the vector. 
	 * @param lhs Left operand.
	 * @param rhs Right operand. 
	 * @return Vec3<T> Difference lhs - rhs. 
	 */
	template <typename T> Vec3<T> operator-(const Vec3<T> &lhs, const Vec3<T> &rhs)
	{
		return Vec3<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

	/**
	 * @brief Inversion of a Vec3 object with respect to addition. 
	 * 
	 * @tparam T Underlying fundamental type of the vector. 
	 * @param rhs Vector object. 
	 * @return Vec3<T> Inverted vector -rhs. 
	 */
	template <typename T> Vec3<T> operator-(const Vec3<T> &rhs)
	{
		return Vec3<T>(-rhs.x, -rhs.y, -rhs.z);
	}

	/**
	 * @brief Comparison of two Vec3 objects. 
	 * 
	 * @tparam T Underlying fundamental type of the vector. 
	 * @param lhs Left operand. 
	 * @param rhs Right operand. 
	 * @return bool Returns true if lhs equls rhs. Returns false otherwise. 
	 */
	template <typename T> bool operator==(const Vec3<T> &lhs, const Vec3<T> &rhs)
	{
		return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
	}

	/**
	 * @brief Four-dimensional vector. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 */
	template <typename T> struct Vec4
	{
		/**
		 * @brief Construct a new Vec4 object with uninitialized components.
		 */
		Vec4() {};

		/**
		 * @brief Construct a new Vec4 object and initialize the x, y, and z components. The w component is set to unity. 
		 * 
		 * @param xyz Value to initialize the x, y, and z component. 
		 */
		Vec4(const T &xyz) : x(xyz), y(xyz), z(xyz), w(1) {};
		
		/**
		 * @brief Construct a new Vec4 object and initialize the x, y, and z components. The w component is set to unity.  
		 * 
		 * @param x Value to initialize the x component. 
		 * @param y Value to initialize the y component. 
		 * @param z Value to initialize the z component. 
		 */
		Vec4(const T &x, const T &y, const T &z) : x(x), y(y), z(z), w(1) {};

		/**
		 * @brief Construct a new Vec4 object and initialize the components. 
		 * 
		 * @param x Value to initialize the x component. 
		 * @param y Value to initialize the y component. 
		 * @param z Value to initialize the z component. 
		 * @param w Value to initialize the w component. 
		 */
		Vec4(const T &x, const T &y, const T &z, const T &w) : x(x), y(y), z(z), w(w) {};

		T x; ///< x Component. 
		T y; ///< y Component. 
		T z; ///< z Component. 
		T w; ///< w Component. 

		/**
		 * @brief Conversion to Vec3 object by dropping the w component. 
		 * 
		 * @return Vec3<T> Reduced Vec3 object with w component dropped. 
		 */
		operator Vec3<T>() const
		{
			return Vec3<T>(x, y, z);
		}
	};
	#pragma endregion

	#pragma region matrix
	/**
	 * @brief 3x3 dimensional matrix. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 */
	template <typename T> struct Mat3
	{
		/**
		 * @brief Construct a new Mat3 with uninitialized values
		 */
		Mat3() {};

		/**
		 * @brief Construct a new Mat3 and initialize all entries with the same value. 
		 *
		 * @param m Value to initialize entries. 
		 */
		Mat3(const T &m)
		{
			data[0][0] = m;
			data[0][1] = m;
			data[0][2] = m;
			data[1][0] = m;
			data[1][1] = m;
			data[1][2] = m;
			data[2][0] = m;
			data[2][1] = m;
			data[2][2] = m;
		};

		/**
		 * @brief Construct a new Mat3 object and initialize the column vectors. 
		 * 
		 * @param a0 First column vector. 
		 * @param a1 Second column vector. 
		 * @param a2 Third column vector. 
		 */
		Mat3(const Vec3<T> &a0, const Vec3<T> &a1, const Vec3<T> &a2)
		{
			data[0][0] = a0.x;
			data[0][1] = a1.x;
			data[0][2] = a2.x;
			data[1][0] = a0.y;
			data[1][1] = a1.y;
			data[1][2] = a2.y;
			data[2][0] = a0.z;
			data[2][1] = a1.z;
			data[2][2] = a2.z;
		}

		/**
		 * @brief Construct a new Mat3 object and initialize entries. 
		 * 
		 * @param m00 Top left entry. 
		 * @param m01 Top middle entry. 
		 * @param m02 Top right entry. 
		 * @param m10 Middle left entry. 
		 * @param m11 Middle middle entry.
		 * @param m12 Middle right entry. 
		 * @param m20 Bottom left entry. 
		 * @param m21 Bottom middle entry. 
		 * @param m22 Bottom right entry. 
		 */
		Mat3(const T &m00, const T &m01, const T &m02, const T &m10, const T &m11, const T &m12, const T &m20, const T &m21, const T &m22)
		{
			data[0][0] = m00;
			data[0][1] = m01;
			data[0][2] = m02;
			data[1][0] = m10;
			data[1][1] = m11;
			data[1][2] = m12;
			data[2][0] = m20;
			data[2][1] = m21;
			data[2][2] = m22;
		}

		T data[3][3]; ///< Matrix entries in row-major order. 

		/**
		 * @brief Generate a 3x3 identity matrix. 
		 * 
		 * @return Mat3<T> Identity matrix. 
		 */
		static Mat3<T> identity()
		{
			return Mat3<T>(1, 0, 0, 0, 1, 0, 0, 0, 1);
		}

		/**
		 * @brief Calculate the determinant of a Mat3 object. 
		 * 
		 * @return T Value of the determinant. 
		 */
		T determinant() const
		{
			return -data[0][2] * data[1][1] * data[2][0] + data[0][1] * data[1][2] * data[2][0] + data[0][2] * data[1][0] * data[2][1] - data[0][0] * data[1][2] * data[2][1] - data[0][1] * data[1][0] * data[2][2] + data[0][0] * data[1][1] * data[2][2];
		}

		/**
		 * @brief Calculate the inverse of a matrix. 
		 * 
		 * @return Mat3<T> Inverse matrix. 
		 */
		Mat3<T> inverse() const
		{
			return (1 / this->determinant()) * Mat3<T>(-data[1][2] * data[2][1] + data[1][1] * data[2][2], data[0][2] * data[2][1] - data[0][1] * data[2][2], -data[0][2] * data[1][1] + data[0][1] * data[1][2], data[1][2] * data[2][0] - data[1][0] * data[2][2], -data[0][2] * data[2][0] + data[0][0] * data[2][2], data[0][2] * data[1][0] - data[0][0] * data[1][2], -data[1][1] * data[2][0] + data[1][0] * data[2][1], data[0][1] * data[2][0] - data[0][0] * data[2][1], -data[0][1] * data[1][0] + data[0][0] * data[1][1]);
		}
	};

	/**
	 * @brief Scalar left-multiplication of a Mat3 object. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 * @param lhs Scalar operand. 
	 * @param rhs Matrix operand. 
	 * @return Mat3<T> Value of lhs * rhs. 
	 */
	template <typename T> Mat3<T> operator*(const T &lhs, const Mat3<T> &rhs)
	{
		Mat3<T> m(T(0));
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				m.data[i][j] = lhs * rhs.data[i][j];
			}
		}
		return m;
	}

	/**
	 * @brief 4x4 dimensional matrix. 
	 * 
	 * @tparam T 
	 */
	template <typename T> struct Mat4
	{
		/**
		 * @brief Construct a new Mat4 object with uninitialized entries. 
		 */
		Mat4() {};

		/**
		 * @brief Construct a new Mat4 object and initialize all entries with the same value.
		 */
		Mat4(const T &m) 
		{
			data[0][0] = m;
			data[0][1] = m;
			data[0][2] = m;
			data[0][3] = m;
			data[1][0] = m;
			data[1][1] = m;
			data[1][2] = m;
			data[1][3] = m;
			data[2][0] = m;
			data[2][1] = m;
			data[2][2] = m;
			data[2][3] = m;
			data[3][0] = m;
			data[3][1] = m;
			data[3][2] = m;
			data[3][3] = m;
		};

		/**
		 * @brief Construct a new Mat 4 object and initialize entires. 
		 * 
		 * @param m00 ///< Matrix entry m[0][0]. 
		 * @param m01 ///< Matrix entry m[0][1]. 
		 * @param m02 ///< Matrix entry m[0][2]. 
		 * @param m03 ///< Matrix entry m[0][3]. 
		 * @param m10 ///< Matrix entry m[1][0]. 
		 * @param m11 ///< Matrix entry m[1][1]. 
		 * @param m12 ///< Matrix entry m[1][2]. 
		 * @param m13 ///< Matrix entry m[1][3]. 
		 * @param m20 ///< Matrix entry m[2][0]. 
		 * @param m21 ///< Matrix entry m[2][1]. 
		 * @param m22 ///< Matrix entry m[2][2]. 
		 * @param m23 ///< Matrix entry m[2][3]. 
		 * @param m30 ///< Matrix entry m[3][0]. 
		 * @param m31 ///< Matrix entry m[3][1]. 
		 * @param m32 ///< Matrix entry m[3][2]. 
		 * @param m33 ///< Matrix entry m[3][3]. 
		 */
		Mat4(const T &m00, const T &m01, const T &m02, const T &m03, const T &m10, const T &m11, const T &m12, const T &m13, const T &m20, const T &m21, const T &m22, const T &m23, const T &m30, const T &m31, const T &m32, const T &m33)
		{
			data[0][0] = m00;
			data[0][1] = m01;
			data[0][2] = m02;
			data[0][3] = m03;
			data[1][0] = m10;
			data[1][1] = m11;
			data[1][2] = m12;
			data[1][3] = m13;
			data[2][0] = m20;
			data[2][1] = m21;
			data[2][2] = m22;
			data[2][3] = m23;
			data[3][0] = m30;
			data[3][1] = m31;
			data[3][2] = m32;
			data[3][3] = m33;
		}

		T data[4][4]; ///< Matrix entries in row-major format. 

		/**
		 * @brief Generate a 4x4 identity matrix.
		 * 
		 * @return Mat4<T> Identity matrix. 
		 */
		static Mat4<T> identity()
		{
			return Mat4<T>(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
		}

		/**
		 * @brief Generate a spatial inversion operation, acting on a three-dimensional vector with fourth entry of unity. 
		 * 
		 * @return Mat4<T> Matrix describing spatial inversion. 
		 */
		static Mat4<T> inversion()
		{
			return Mat4<T>(-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1);
		}

		/**
		 * @brief Generate a spatial translation operation, acting on a three-dimensional vector with fourth entry of unity. 
		 * 
		 * @param v Translation vector. 
		 * @return Mat4<T> Matrix describing the spatial translation. 
		 */
		static Mat4<T> translation(const Vec3<T> &v)
		{
			return Mat4<T>(1, 0, 0, v.x, 0, 1, 0, v.y, 0, 0, 1, v.z, 0, 0, 0, 1);
		}

		/**
		 * @brief Generate a spatial rotation around a specific axis by a specific angle, acting on a three-dimensional vector with fourth entry of unity. 
		 * 
		 * @param axis Rotation axis. 
		 * @param angle Rotation angle. 
		 * @return Mat4<T> Matrix describing the spatial rotation operation. 
		 */
		static Mat4<T> rotation(const Vec3<T> &axis, T angle)
		{
			Vec3<T> a = axis;
			a.normalize();
			return Mat4<T>(
				a.x * a.x * (1 - cos(angle)) + cos(angle), 
				a.x * a.y * (1 - cos(angle)) - a.z * sin(angle),
				a.x * a.z * (1 - cos(angle)) + a.y * sin(angle),
				0,

				a.y * a.x * (1 - cos(angle)) + a.z * sin(angle),
				a.y * a.y * (1 - cos(angle)) + cos(angle),
				a.y * a.z * (1 - cos(angle)) - a.x * sin(angle),
				0,

				a.z * a.x * (1 - cos(angle)) - a.y * sin(angle),
				a.z * a.y * (1 - cos(angle)) + a.x * sin(angle),
				a.z * a.z * (1 - cos(angle)) + cos(angle), 
				0,

				0, 
				0, 
				0, 
				1);
		}

		/**
		 * @brief Generate a spatial rotation around a specific rotation center and axis by a specific angle, acting on a three-dimensional vector with fourth entry of unity. 
		 * 
		 * @param axis Rotation axis. 
		 * @param point Rotation center. 
		 * @param angle Rotation angle. 
		 * @return Mat4<T> Matrix describing the spatial rotation operation. 
		 */
		static Mat4<T> rotation(const Vec3<T> &axis, const Vec3<T> &point, const T &angle)
		{
			return translation(point) * rotation(axis, angle) * translation(-point);
		}
	};

	/**
	 * @brief Multiplication of two Mat4 objects. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 * @param lhs Left operand. 
	 * @param rhs Right operand. 
	 * @return Mat4<T> Value of the product lhs * rhs. 
	 */
	template <typename T> Mat4<T> operator*(const Mat4<T> &lhs, const Mat4<T> &rhs)
	{
		Mat4<T> m(T(0));
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				for (int k = 0; k < 4; ++k) m.data[i][j] += lhs.data[i][k] * rhs.data[k][j];
			}
		}
		return m;
	}
	#pragma endregion

	#pragma region matrixvectoroperations
	/**
	 * @brief Matrix vector multiplication of 3x3 matrix and a three-dimensional vector. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 * @param lhs Left operand. 
	 * @param rhs Right operand. 
	 * @return Vec3<T> Value of lhs * rhs. 
	 */
	template <typename T> Vec3<T> operator*(const Mat3<T> &lhs, const Vec3<T> &rhs)
	{
		return Vec3<T>(
			lhs.data[0][0] * rhs.x + lhs.data[0][1] * rhs.y + lhs.data[0][2] * rhs.z,
			lhs.data[1][0] * rhs.x + lhs.data[1][1] * rhs.y + lhs.data[1][2] * rhs.z,
			lhs.data[2][0] * rhs.x + lhs.data[2][1] * rhs.y + lhs.data[2][2] * rhs.z
			);
	}

	/**
	 * @brief Matrix vector multiplication of 4x4 matrix and a three-dimensional vector by padding the vector with 1 in the fourth dimension and returning only the first three entries. 
	 * 
	 * @tparam T Underlying fundamental data type. 
	 * @param lhs Left operand. 
	 * @param rhs Right operand. 
	 * @return Vec3<T> Result of lhs * rhs. 
	 */
	template <typename T> Vec3<T> operator*(const Mat4<T> &lhs, const Vec3<T> &rhs)
	{
		return Vec3<T>(
			lhs.data[0][0] * rhs.x + lhs.data[0][1] * rhs.y + lhs.data[0][2] * rhs.z + lhs.data[0][3],
			lhs.data[1][0] * rhs.x + lhs.data[1][1] * rhs.y + lhs.data[1][2] * rhs.z + lhs.data[1][3],
			lhs.data[2][0] * rhs.x + lhs.data[2][1] * rhs.y + lhs.data[2][2] * rhs.z + lhs.data[2][3]
			);
	}
	#pragma endregion
}