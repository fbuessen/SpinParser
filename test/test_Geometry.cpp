#define BOOST_TEST_MODULE "GeometryTest"
#include <boost/test/included/unit_test.hpp>
#define _USE_MATH_DEFINES
#include <math.h>
#include "lib/Geometry.hpp"

BOOST_TEST_DONT_PRINT_LOG_VALUE(geometry::Vec3<double>);

typedef geometry::Vec3<double> v3;
typedef geometry::Vec4<double> v4;
typedef geometry::Mat3<double> m3;
typedef geometry::Mat4<double> m4;

BOOST_AUTO_TEST_SUITE(GeometryTest);

BOOST_AUTO_TEST_CASE(Vec3Constructor)
{
	v3 x;
	
	x = v3(1.0);
	BOOST_CHECK_EQUAL(x.x, 1.0);
	BOOST_CHECK_EQUAL(x.y, 1.0);
	BOOST_CHECK_EQUAL(x.z, 1.0);

	x = v3(1.0, 2.0, 3.0);
	BOOST_CHECK_EQUAL(x.x, 1.0);
	BOOST_CHECK_EQUAL(x.y, 2.0);
	BOOST_CHECK_EQUAL(x.z, 3.0);
};

BOOST_AUTO_TEST_CASE(Vec3Normalization)
{
	v3 x = v3(1.0, 2.0, 3.0);
	BOOST_CHECK_EQUAL(x.norm(), sqrt(14.0));
	x.normalize();
	BOOST_CHECK_EQUAL(x.x, 1.0/sqrt(14.0));
	BOOST_CHECK_EQUAL(x.y, 2.0/sqrt(14.0));
	BOOST_CHECK_EQUAL(x.z, 3.0/sqrt(14.0));
}

BOOST_AUTO_TEST_CASE(Vec4Constructor)
{
	v4 x;

	x = v4(1.0, 2.0, 3.0);
	BOOST_CHECK_EQUAL(x.x, 1.0);
	BOOST_CHECK_EQUAL(x.y, 2.0);
	BOOST_CHECK_EQUAL(x.z, 3.0);
	BOOST_CHECK_EQUAL(x.w, 1.0);

	x = v4(1.0, 2.0, 3.0, 4.0);
	BOOST_CHECK_EQUAL(x.x, 1.0);
	BOOST_CHECK_EQUAL(x.y, 2.0);
	BOOST_CHECK_EQUAL(x.z, 3.0);
	BOOST_CHECK_EQUAL(x.w, 4.0);

	v3 y(x);
	BOOST_CHECK_EQUAL(y.x, 1.0);
	BOOST_CHECK_EQUAL(y.y, 2.0);
	BOOST_CHECK_EQUAL(y.z, 3.0);
};

BOOST_AUTO_TEST_CASE(Mat3Constructor)
{
	m3 x;

	x = m3(1.0);
	for (int i = 0; i < 9; ++i) BOOST_CHECK_EQUAL((&x.data[0][0])[i], 1.0);

	v3 y(1.0, 2.0, 3.0);
	x = m3(y, y, y);
	for (int i = 0; i < 3; ++i)
	{
		BOOST_CHECK_EQUAL(x.data[0][i], 1.0);
		BOOST_CHECK_EQUAL(x.data[1][i], 2.0);
		BOOST_CHECK_EQUAL(x.data[2][i], 3.0);
	}

	x = m3(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);
	for (int i = 0; i < 9; ++i) BOOST_CHECK_EQUAL((&x.data[0][0])[i], double(i));

	x = m3::identity();
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (i == j) BOOST_CHECK_EQUAL(x.data[i][j], 1.0);
			else BOOST_CHECK_EQUAL(x.data[i][j], 0.0);
		}
	}
}

BOOST_AUTO_TEST_CASE(Mat3Operations)
{
	m3 x(1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);

	BOOST_CHECK_EQUAL(x.determinant(), -3.0);

	x = x.inverse();
	BOOST_CHECK_EQUAL(x.data[0][0], 1.0);
	BOOST_CHECK_EQUAL(x.data[0][1], -2.0);
	BOOST_CHECK_EQUAL(x.data[0][2], 1.0);
	BOOST_CHECK_EQUAL(x.data[1][0], -2.0);
	BOOST_CHECK_EQUAL(x.data[1][1], 4.0/3.0);
	BOOST_CHECK_EQUAL(x.data[1][2], -1.0/3.0);
	BOOST_CHECK_EQUAL(x.data[2][0], 1.0);
	BOOST_CHECK_EQUAL(x.data[2][1], 1.0/3.0);
	BOOST_CHECK_EQUAL(x.data[2][2], -1.0/3.0);
}

BOOST_AUTO_TEST_CASE(Mat4Constructor)
{
	m4 x;

	x = m4(1.0);
	for (int i = 0; i < 16; ++i) BOOST_CHECK_EQUAL((&x.data[0][0])[i], 1.0);

	x = m4(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0);
	for (int i = 0; i < 16; ++i) BOOST_CHECK_EQUAL((&x.data[0][0])[i], double(i));

	x = m4::identity();
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (i == j) BOOST_CHECK_EQUAL(x.data[i][j], 1.0);
			else BOOST_CHECK_EQUAL(x.data[i][j], 0.0);
		}
	}

	v3 v;
	v = m4::inversion() * v3(1.0, 2.0, 3.0);
	BOOST_CHECK_CLOSE(v.x, -1.0, 1e-10);
	BOOST_CHECK_CLOSE(v.y, -2.0, 1e-10);
	BOOST_CHECK_CLOSE(v.z, -3.0, 1e-10);

	v = m4::translation(v3(1.0, 2.0, 3.0)) * v3(1.0, 2.0, 3.0);
	BOOST_CHECK_CLOSE(v.x, 2.0, 1e-10);
	BOOST_CHECK_CLOSE(v.y, 4.0, 1e-10);
	BOOST_CHECK_CLOSE(v.z, 6.0, 1e-10);

	v = m4::rotation(v3(1.0, 1.0, 1.0), 2 * M_PI / 3.0) * v3(1.0, 0.0, 0.0);
	BOOST_CHECK_SMALL(v.x, 1e-10);
	BOOST_CHECK_CLOSE(v.y, 1.0, 1e-10);
	BOOST_CHECK_SMALL(v.z, 1e-10);

	v = m4::rotation(v3(0.0, 0.0, 1.0), v3(1.0, 0.0, 0.0), M_PI) * v3(2.0, 0.0, 0.0);
	BOOST_CHECK_SMALL(v.x, 1e-10);
	BOOST_CHECK_SMALL(v.y, 1e-10);
	BOOST_CHECK_SMALL(v.z, 1e-10);
}

BOOST_AUTO_TEST_CASE(VectorArithmetics)
{
	double x = geometry::dot(v3(1.0, 2.0, 3.0), v3(2.0, 3.0, 4.0));
	BOOST_CHECK_CLOSE(x, 20.0, 1e-10);

	v3 y = geometry::cross(v3(1.0, 2.0, 3.0), v3(2.0, 3.0, 4.0));
	BOOST_CHECK_CLOSE(y.x, -1.0, 1e-10);
	BOOST_CHECK_CLOSE(y.y, 2.0, 1e-10);
	BOOST_CHECK_CLOSE(y.z, -1.0, 1e-10);

	y = v3(1.0, 2.0, 3.0) + v3(1.0, 2.0, 3.0);
	BOOST_CHECK_CLOSE(y.x, 2.0, 1e-10);
	BOOST_CHECK_CLOSE(y.y, 4.0, 1e-10);
	BOOST_CHECK_CLOSE(y.z, 6.0, 1e-10);

	y = v3(1.0, 2.0, 3.0) - v3(2.0, 3.0, 4.0);
	BOOST_CHECK_CLOSE(y.x, -1.0, 1e-10);
	BOOST_CHECK_CLOSE(y.y, -1.0, 1e-10);
	BOOST_CHECK_CLOSE(y.z, -1.0, 1e-10);

	y = -v3(1.0, 2.0, 3.0);
	BOOST_CHECK_CLOSE(y.x, -1.0, 1e-10);
	BOOST_CHECK_CLOSE(y.y, -2.0, 1e-10);
	BOOST_CHECK_CLOSE(y.z, -3.0, 1e-10);

	bool b = v3(1.0, 2.0, 3.0) == v3(1.0, 2.0, 3.0);
	BOOST_CHECK(b == true);
	b = v3(1.0, 2.0, 3.0) == v3(2.0, 2.0, 3.0);
	BOOST_CHECK(b == false);
}

BOOST_AUTO_TEST_CASE(MatrixArithmetics)
{
	m3 x = 2.0 * m3::identity();
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (i == j) BOOST_CHECK_EQUAL(x.data[i][j], 2.0);
			else BOOST_CHECK_EQUAL(x.data[i][j], 0.0);
		}
	}

	v3 y = m4::rotation(v3(1.0, 1.0, 1.0), M_PI) * m4::rotation(v3(1.0, 1.0, 1.0), M_PI) * v3(2.0, 3.0, 4.0);
	BOOST_CHECK_CLOSE(y.x, 2.0, 1e-10);
	BOOST_CHECK_CLOSE(y.y, 3.0, 1e-10);
	BOOST_CHECK_CLOSE(y.z, 4.0, 1e-10);

	y = m3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0) * v3(2.0, 3.0, 4.0);
	BOOST_CHECK_CLOSE(y.x, 20.0, 1e-10);
	BOOST_CHECK_CLOSE(y.y, 47.0, 1e-10);
	BOOST_CHECK_CLOSE(y.z, 74.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END();