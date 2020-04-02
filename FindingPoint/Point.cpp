#include "Point.h"


using namespace std;

double Point::norm() const
{
	return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

double Point::scalar_product(const Point& b) const
{
	return x * b.x + y * b.y + z * b.z;
}

Point Point::cross(const Point& a, const Point& b)
{
	Point result;

	result.x = a.y * b.z - a.z * b.y;
	result.y = a.z * b.x - a.x * b.z;
	result.z = a.x * b.y - a.y * b.x;

	return result;
}

std::ostream& operator<< (std::ostream& out, const Point& n)
{
	out << "\n\tX - " << n.x << "\n\tY - " << n.y << "\n\tZ - " << n.z << endl;
	return out;
}

Point operator+(const Point& a, const Point& b)
{
	Point result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;

	return result;
}

Point operator-(const Point& a, const Point& b)
{
	Point result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;

	return result;
}

Point operator/(const Point& a, const Point& b)
{
	Point result;
	result.x = a.x / b.x;
	result.y = a.y / b.y;
	result.z = a.z / b.z;

	return result;
}

Point operator*(const Point& a, const Point& b)
{
	Point result;
	result.x = a.x * b.x;
	result.y = a.y * b.y;
	result.z = a.z * b.z;

	return result;
}

Point operator*(const Point& a, double d)
{
	Point result;
	result.x = a.x * d;
	result.y = a.y * d;
	result.z = a.z * d;

	return result;
}

Point operator/(const Point& a, double d)
{
	Point result;
	result.x = a.x / d;
	result.y = a.y / d;
	result.z = a.z / d;

	return result;
}