#pragma once

#include <ostream>
#include <iostream>

struct Point
{
	double x;
	double y;
	double z;

	friend std::ostream& operator<< ( std::ostream& out, const Point& point);
	friend Point operator+(const Point& a, const Point& b);
	friend Point operator-(const Point& a, const Point& b);
	friend Point operator/(const Point& a, const Point& b);
	friend Point operator*(const Point& a, const Point& b);
	friend Point operator*(const Point& a, double d);
	friend Point operator/(const Point& a, double d);

	double norm();
	double scalar_product(const Point& b);
	static Point cross(const Point& a, const Point& b);
};

struct Points
{
	Points() = default;

	Points(Point a, Point b, Point c)
	{
		A = a;
		B = b;
		C = c;
	}
	
	Points(Point a, Point b, Point c, Point d)
	{
		A = a;
		B = b;
		C = c;
		D = d;
	}
	
	Point A{};
	Point B{};
	Point C{};
	Point D{};
};