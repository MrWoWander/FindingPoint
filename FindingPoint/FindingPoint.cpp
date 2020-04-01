#define _USE_MATH_DEFINES
#include <iostream>
#include "omp.h"
#include <cmath>
#include <vector>
#include "Point.h"

using namespace std;

vector<Point> point;
vector<Points> PointsArray;
vector<Points> PointDArray;

double h = 3;
float AB = 1.52f;
float BC = 2.6f;
float AC = sqrt(AB * AB + BC * BC);

int quality = 10;
float stepR = 0.1f;
int stepV = 360;

void findXYZ()
{
	for (double v = 0; v < 2 * M_PI; v += (2 * M_PI) / stepV)
		for (double r = 3.1; r < 7.1; r += stepR)
		{
			Point p;
			p.x = r * cos(v);
			p.y = r * sin(v);
			p.z = (h / (2 * M_PI)) * v;
			
			point.push_back(p);
		}
}

void findNormalPoint()
{
	cout << "Старт поиска точек" << endl;

#pragma omp parallel for
	for (int i = 0; i < point.size() - 1; i++)
	{
		for (int j = i + 1; j < point.size(); j++)
		{
			float ab = pow(point[i].x - point[j].x, 2) +
				pow(point[i].y - point[j].y, 2) +
				pow(point[i].z - point[j].z, 2);

			if (round(ab * quality) / quality == round(pow(AB, 2) * quality) / quality)
			{
				for (int k = 0; k < point.size(); k++)
				{
					float bc = pow(point[j].x - point[k].x, 2) +
						pow(point[j].y - point[k].y, 2) +
						pow(point[j].z - point[k].z, 2);

					if (round(bc * quality) / quality == round(pow(BC, 2) * quality) / quality)
					{
						float abc = pow(point[i].x - point[k].x, 2) +
							pow(point[i].y - point[k].y, 2) +
							pow(point[i].z - point[k].z, 2);

						if (round(abc * quality) / quality == round((pow(AC, 2) * quality) / quality))
						{
							
							PointsArray.emplace_back(point[i], point[j], point[k]);
						}
					}
				}
			}
		}
	}
}

void findDPoint(Point a, Point b, Point c)
{
	Point temp1 = b - a;
	double t = temp1.norm();

	Point e_x = temp1 / t;

	Point temp2 = c - a;
	double i = e_x.scalar_product(temp2);

	Point temp3 = temp2 - e_x * i;
	Point e_y = temp3 / temp3.norm();

	Point e_z = Point::cross(e_x, e_y);

	double j = e_y.scalar_product(temp2);

	const double r1 = BC;
	const double r2 = AC;
	const double r3 = AB;
	
	double x = (r1 * r1 - r2 * r2 + t * t) / (2 * t);
	double y = (r1 * r1 - r3 * r3 - 2 * i * x + i * i + j * j) / (2 * j);

	double temp4 = r1 * r1 - x * x - y * y;

	if (temp4 < 0)
	{
		std::cout << "temp4 < 0" << endl;
		return;
	}

	double z = sqrt(temp4);

	//Point p_12_a = a + e_x * x + e_y * y + e_z * z;
	Point d = a + e_x * x + e_y * y - e_z * z;

	
	
	PointDArray.emplace_back(a, b, c, d);
}

void findExcessZ()
{
	cout << "Поиск превышения Z" << endl;

	vector<double> excessZ;

	
	for (int i = 0; i < point.size(); i++)
	{
		for (int j = 0; j < PointDArray.size(); j++)
		{
			if (round(PointDArray[j].D.x * 10) / 10 == round(point[i].x * 10) / 10)
				if (round(PointDArray[j].D.y * 10) / 10 == round(point[i].y * 10) / 10)
				{
					double z = round(PointDArray[j].D.z * 10) / 10 - round(point[i].z * 10) / 10;
					excessZ.push_back(z);

					if (z > 1.5)
					{
						cout << "A: " << PointDArray[j].A << endl;
						cout << "B: " << PointDArray[j].B << endl;
						cout << "C: " << PointDArray[j].C << endl;
						cout << "D: " << PointDArray[j].D << endl;
						cout << "Point: " << point[i] << endl;
						
						cout << "Z у D: " << round(PointDArray[j].D.z * 10) / 10 << endl;
						cout << "Z у стандартной точки: " << round(point[i].z * 10) / 10 << endl;
						cout << "Z итоговая: " << z << endl;
						cout << "\n";
					}
				}
		}

		//cout << "Осталось " << point.size() - i + 1 << " точек" << endl;
	}

	cout << "Количество Z: " << excessZ.size() << endl;

	double maxExcessZ = 0;
	double sumExcessZ = 0;

	for (auto b : excessZ)
	{
		if (b > maxExcessZ)
			maxExcessZ = b;

		sumExcessZ += b;
	}

	cout << "Максимальный Z: " << maxExcessZ << endl;
	cout << "Средний Z: " << sumExcessZ / excessZ.size() << endl;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Hello World!\n";

	findXYZ();

	cout << point.size() << endl;

	findNormalPoint();
	cout << "Пар: " << PointsArray.size() << endl;


	if (!PointsArray.empty())
		for (int i = 0; i < PointsArray.size(); i++)
			findDPoint(PointsArray[i].A, PointsArray[i].B, PointsArray[i].C);

	if (!PointDArray.empty())
		findExcessZ();


	system("pause");
}
