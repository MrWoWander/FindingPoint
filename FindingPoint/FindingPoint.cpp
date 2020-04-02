#define _USE_MATH_DEFINES
#include <iostream>
#include "omp.h"
#include <cmath>
#include <vector>
#include <fstream>
#include "Point.h"

using namespace std;

void designation_variable();

vector<Point> point;
vector<Points> PointsArray;
vector<Points> PointDArray;

double h = 3;
float stepR = 0.1f;
int stepV = 360;
float start_r = 3.1f;
float end_r = 7.1f;

float AB = 1.52f;
float BC = 2.6f;
float AC = sqrt(AB * AB + BC * BC);

void find_direct_helicoid()
{
	cout << "Поиск прямого геликоида" << endl;

	designation_variable();

	for (double v = 0; v < 2 * M_PI; v += 2 * M_PI / stepV)
		for (double r = start_r; r < end_r; r += stepR)
		{
			Point p;
			p.x = r * cos(v);
			p.y = r * sin(v);
			p.z = h / (2 * M_PI) * v;

			point.push_back(p);
		}
}

void pseudo_developable_helicoid()
{
	cout << "Поиск псевдо-развёртывающего геликоида" << endl;

	designation_variable();

	float a = 0.f;

	cout << "Задай a: ";
	cin >> a;

	for (double v = 0; v < 2 * M_PI; v += 2 * M_PI / stepV)
		for (double r = start_r; r < end_r; r += stepR)
		{
			Point p;
			p.x = a * cos(v) - r * sin(v);
			p.y = a * sin(v) + r * cos(v);
			p.z = h / (2 * M_PI) * v;

			point.push_back(p);
		}
}

void developable_helicoid()
{
	cout << "Поиск развёртывающего геликоида" << endl;

	designation_variable();

	float a = 0.f;
	float y = 0.f;

	cout << "Задай a: ";
	cin >> a;

	for (double v = 0; v < 2 * M_PI; v += 2 * M_PI / stepV)
		for (double r = start_r; r < end_r; r += stepR)
		{
			Point p;

			const double m = sqrt(pow(a, 2) + pow(h / (2 * M_PI), 2));
			
			p.x = a * cos(v) - a / m * r * sin(v);
			p.y = a * sin(v) + a / m * r * cos(v);
			p.z = h / (2 * M_PI) * v + r * (h / (2 * M_PI)) / m;

			point.push_back(p);
		}
}

void convolute_helicoid()
{
	cout << "Поиск конволютного геликоида" << endl;

	designation_variable();

	float a = 0.f;
	float y = 0.f;

	cout << "Задай a: ";
	cin >> a;

	cout << "Задай y: ";
	cin >> y;

	for (double v = 0; v < 2 * M_PI; v += 2 * M_PI / stepV)
		for (double r = start_r; r < end_r; r += stepR)
		{
			Point p;
			const double RadY = y * (180 / M_PI);

			p.x = a * cos(v) - r * sin(v) * sin(RadY);
			p.y = a * sin(v) + r * cos(v) * sin(RadY);
			p.z = h / (2 * M_PI) * v + r * cos(RadY);

			point.push_back(p);
		}
}

void find_oblique_helicoid()
{
	cout << "Поиск косого геликоида" << endl;

	designation_variable();

	float k = 0.f;

	cout << "Задай k: ";
	cin >> k;

	for (double v = 0; v < 2 * M_PI; v += 2 * M_PI / stepV)
		for (double r = start_r; r < end_r; r += stepR)
		{
			Point p;
			p.x = r * cos(v);
			p.y = r * sin(v);
			p.z = h / (2 * M_PI) * v + k * r;

			point.push_back(p);
		}
}

const int quality = 10;

void findNormalPoint()
{
	cout << "Старт поиска точек" << endl;

	omp_lock_t myLock;
	omp_init_lock(&myLock);

	int count = 0;

#pragma omp parallel for shared(point, PointsArray, count)
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
							omp_set_lock(&myLock);
							Points points = Points(point[i], point[j], point[k]);

							PointsArray.push_back(points);

							omp_unset_lock(&myLock);
						}
					}
				}
			}
		}
		omp_set_lock(&myLock);
		count++;
		cout << "Осталось " << point.size() - count << " точек" << endl;
		omp_unset_lock(&myLock);
	}
}

void findDPoint(Point a, Point b, Point c)
{
	const Point temp1 = b - a;
	const double t = temp1.norm();

	const Point e_x = temp1 / t;

	const Point temp2 = c - a;
	const double i = e_x.scalar_product(temp2);

	const Point temp3 = temp2 - e_x * i;
	const Point e_y = temp3 / temp3.norm();

	const Point e_z = Point::cross(e_x, e_y);

	const double j = e_y.scalar_product(temp2);

	const double r1 = BC;
	const double r2 = AC;
	const double r3 = AB;

	const double x = (r1 * r1 - r2 * r2 + t * t) / (2 * t);
	const double y = (r1 * r1 - r3 * r3 - 2 * i * x + i * i + j * j) / (2 * j);

	const double temp4 = r1 * r1 - x * x - y * y;

	if (temp4 < 0)
	{
		cout << "temp4 < 0 " << endl;

		cout << "temp4: " << temp4 << endl;

		cout << "A: " << a << endl;
		cout << "B: " << b << endl;
		cout << "C: " << c << endl;
		cout << endl;

		return;
	}

	double z = sqrt(temp4);

	//Point p_12_a = a + e_x * x + e_y * y + e_z * z;
	Point d = a + e_x * x + e_y * y - e_z * z;


	PointDArray.emplace_back(a, b, c, d);
	
	//Point d;
	//d.x = 2 * a.x - b.x + c.x;
	//d.y = 2 * a.y - b.y + c.y;
	//d.z = 2 * a.z - b.z + c.z;

	//PointDArray.emplace_back(a, b, c, d);
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

void designation_variable()
{
	cout << "Задайте h: ";
	cin >> h;

	cout << "Задайте минимальный r: ";
	cin >> start_r;

	cout << "Задайте максимальный r: ";
	cin >> end_r;
}

int main()
{
	setlocale(LC_ALL, "Russian");

	bool isMode = true;
	int tryMode = 0;
	int mode = 0;

	while (isMode)
	{
		cout << "Выберите режим:\n"
			<< "\t1 - Прямой геликоид\n"
			<< "\t2 - Косой геликоид\n"
			<< "\t3 - Псевдо-развёртывающийся геликоид\n"
			<< "\t4 - Конволютный  геликоид\n"
			<< "\t5 - Развёртывающийся геликоид\n"
			<< endl;
		cout << "Режим: ";
		cin >> mode;
		cout << endl;
		
		switch (mode) {
		case 1:

			isMode = false;
			find_direct_helicoid();

			break;


		case 2:

			isMode = false;
			find_oblique_helicoid();

			break;

		case 3:

			isMode = false;
			pseudo_developable_helicoid();

			break;

		case 4:

			isMode = false;
			convolute_helicoid();

			break;

		case 5:

			isMode = false;
			developable_helicoid();

			break;

		default:

			if (tryMode < 2)
			{
				tryMode++;
				cout << "Такого режима нет!\n" << endl;
				continue;
			}

			cout << "Иди нахуй!" << endl;
			return 0;
		}
	}

	cout << "\nВсего точек: " << point.size() << endl;

	findNormalPoint();
	cout << "Пар: " << PointsArray.size() << endl;


	if (!PointsArray.empty())
		for (int i = 0; i < PointsArray.size(); i++)
			findDPoint(PointsArray[i].A, PointsArray[i].B, PointsArray[i].C);

	if (!PointDArray.empty())
		findExcessZ();


	system("pause");
}
