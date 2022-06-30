#include <vector>
#include <random>
#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <cstring>

//#define dimension 1
//#define number_of_dots 200
//#define N1 7
//#define N2 7

//#define dimension 2
//#define number_of_dots 100000
//#define N1 40
//#define N2 100

#define dimension 3
#define number_of_dots 2000000
#define N1 100
#define N2 100

int read_input (double values[])
{
	FILE *finp;
	char param[6];
	double value;
	int i = 0;
	if ((finp = fopen ("cluster.inp", "rt")) != NULL)
	{
		//std::cout << "File " << fname << " has been opened successfully!" << std::endl;
	}
	else
	{
		//std::cout << "File " << fname << " couldn't be opened. Program finished with error." << std::endl;
		return 1;
	}
	while (1)
	{
		if (fscanf (finp, "%5c %lf\n", &param, &value) != 2)
		{
			break;
		}
		values[i] = value;
//		std::cout << values[i] << std::endl;
		i++;
	}

	fclose (finp);
	return 0;
}

int generate_dots (double main_array[], double average_radius, double sigma, double cdse, int seed)
{
	double r_min, r_max;
	if (cdse > 0.5)
	{
		r_min = 1.080611;
		r_max = 2.415654;
	}
	else
	{
		r_min = 1.487189;
		r_max = 1.996367;
	}
	sigma *= average_radius*0.01;
	srand(time(NULL)+seed);
	std::default_random_engine dre(time(NULL)+100*seed);
	for (int i = 2; i < number_of_dots*6+2; i+=6)
	{
		std::normal_distribution<double> norm_distr(average_radius, sigma);
		main_array[i] = norm_distr(dre);
		if (main_array[i] < r_min || main_array[i] > r_max
				|| main_array[i] > average_radius+3*sigma
				|| main_array[i] < average_radius-3*sigma)
		{
			i-=6;
			continue;
		}
	}

	return 0;
}

int init1 (double main_array[], double shell, int neighbors[][N1]);
int init2 (double main_array[], double shell, int neighbors[][N1], int neighbors2[][N2]);
int init3 (double main_array[], double shell, int neighbors[][N1], int neighbors2[][N2]);
int add1 (double main_array[], double shell, int neighbors[][N1]);
int add2 (double main_array[], double average_radius, double shell, int neighbors[][N1], int neighbors2[][N2]);
int add3 (double main_array[], double average_radius, double shell, int neighbors[][N1], int neighbors2[][N2]);
int generate_cluster (double main_array[], double average_radius, double shell, int neighbors[][N1], int neighbors2[][N2], double values[])
{
//	std::cout << "generate_cluster" << std::endl;
	double t0 = time(NULL); //время начала работы функции
	for (int i = 0; i < number_of_dots; i++)
	{
		neighbors[i][0] = 0;
		neighbors2[i][0] = 0;
		for (int j = 1; j < N1; j++)
		{
			neighbors[i][j] = -1;
		}
		for (int j = 1; j < N2; j++)
		{
			neighbors2[i][j] = -1;
		}
	}

	if (dimension == 1)
	{
//		std::cout << "dimension = 1" << std::endl;
		init1 (main_array, shell, neighbors);
		while (add1 (main_array, shell, neighbors));
		return 0;
	}
	if (dimension == 2)
	{
//		std::cout << "dimension = 2" << std::endl;
		init2 (main_array, shell, neighbors, neighbors2);
		while (add2 (main_array, average_radius, shell, neighbors, neighbors2));
	}
	if (dimension == 3)
	{
//		std::cout << "dimension = 3" << std::endl;
		init3 (main_array, shell, neighbors, neighbors2);
		while (add3 (main_array, average_radius, shell, neighbors, neighbors2));
	}
	FILE *fpdb;
	std::ostringstream name_file_pdb_o;
	name_file_pdb_o << "pdb-dim=" << dimension
			<< "-r=" << values[0]
			<< "-sig=" << values[1]
			<< "-l=" << values[2]
			<< "-k=" << values[3]
			<< "-phi=" << values[4]
			<< "-vphi=" << values[5] << ".pdb";
	std::string name_file_pdb = name_file_pdb_o.str();
	if ((fpdb = fopen (name_file_pdb.c_str(), "wt")) != NULL)
	{
		for (int i = 0; i < main_array[0]; i++)
		{
			fprintf (fpdb, "ATOM  *****  CCC XXX     1    %8.3lf%8.3lf%8.3lf  1.00  0.00\n",0.1*main_array[i*6+3], 0.1*main_array[i*6+4], 0.1*main_array[i*6+5]);
		}
	}
	fclose (fpdb);
	std::cout << "Кластер создан за " << time(NULL)-t0 << " сек." << std::endl;
	return 0;
}

int init1 (double main_array[], double shell, int neighbors[][N1])
{
	main_array[3] = 0.0; //Х координата первой точки -0- расположение точек
	main_array[4] = 0.0; //Y координата первой точки
	main_array[5] = 0.0; //Z координата первой точки
	main_array[9] = main_array[2] + main_array[8] + shell*2.; //X координата второй точки
	main_array[10] = 0.; //Y координата второй точки -0-1-
	main_array[11] = 0.; //Z координата второй точки
	main_array[15] = - main_array[2] - main_array[14] - shell*2.; //X третьей точки -2-0-1-
	main_array[16] = 0.0; //Y координата третьей точки
	main_array[17] = 0.0; //Z координата третьей точки
	main_array[21] = main_array[9] + main_array[8] + main_array[20] + shell*2.; //-2-0-1-3-
	main_array[22] = 0.0; //Y координата четвёртой точки
	main_array[23] = 0.0; //Z координата четвёртой точки
	main_array[27] = main_array[15] - main_array[14] - main_array[26] - shell*2.; //-4-2-0-1-3-
	main_array[28] = 0.0; //Y координата пятой точки
	main_array[29] = 0.0; //Z координата пятой точки
	main_array[33] = main_array[21] + main_array[20] + main_array[32] + shell*2.; //-4-2-0-1-3-5-
	main_array[34] = 0.0; //Y координата шестой точки
	main_array[35] = 0.0; //Z координата шестой точки
	main_array[39] = main_array[27] - main_array[26] - main_array[38] - shell*2.; //-6-4-2-0-1-3-5-
	main_array[40] = 0.0; //Y координата седьмой точки
	main_array[41] = 0.0; //Z координата седьмой точки
	main_array[0] = 7; //Общее количество точек на текущий момент, для которых вычислены координаты
	neighbors[0][0] = 6; //У первой частицы 6 соседей
	neighbors[0][1] = 1; //Номер соседа
	neighbors[0][2] = 2; //Номер соседа
	neighbors[0][3] = 3; //Номер соседа
	neighbors[0][4] = 4; //Номер соседа
	neighbors[0][5] = 5; //Номер соседа
	neighbors[0][6] = 6; //Номер соседа
	neighbors[1][0] = 5; //У второй частицы 5 соседей
	neighbors[1][1] = 0; //Номер соседа
	neighbors[1][2] = 3; //Номер соседа
	neighbors[1][3] = 2; //Номер соседа
	neighbors[1][4] = 5; //Номер соседа
	neighbors[1][5] = 4; //Номер соседа
	neighbors[2][0] = 5; //У третьей частицы 5 соседей
	neighbors[2][1] = 0; //Номер соседа
	neighbors[2][2] = 1; //Номер соседа
	neighbors[2][3] = 3; //Номер соседа
	neighbors[2][4] = 4; //Номер соседа
	neighbors[2][5] = 6; //Номер соседа
	neighbors[3][0] = 4; //У четвёртой частицы 4 соседа
	neighbors[3][1] = 1; //Номер соседа
	neighbors[3][2] = 5; //Номер соседа
	neighbors[3][3] = 0; //Номер соседа
	neighbors[3][4] = 2; //Номер соседа
	neighbors[4][0] = 4; //У пятой частицы 4 соседа
	neighbors[4][1] = 2; //Номер соседа
	neighbors[4][2] = 6; //Номер соседа
	neighbors[4][3] = 0; //Номер соседа
	neighbors[4][4] = 1; //Номер соседа
	neighbors[5][0] = 3; //У шестой частицы 3 соседа
	neighbors[5][1] = 3; //Номер соседа
	neighbors[5][2] = 1; //Номер соседа
	neighbors[5][3] = 0; //Номер соседа
	neighbors[6][0] = 3; //У седьмой частицы 3 соседа
	neighbors[6][1] = 4; //Номер соседа
	neighbors[6][2] = 2; //Номер соседа
	neighbors[6][3] = 0; //Номер соседа
	return 0;
}

int init2 (double main_array[], double shell, int neighbors[][N1], int neighbors2[][N2])
{
//	std::cout << "init2" << std::endl;
	main_array[3] = 0.0;
	main_array[4] = 0.0;
	main_array[5] = 0.0;
	main_array[9] = 0.0;
	main_array[10] = main_array[2]+main_array[8]+shell*2.;
	main_array[11] = 0.0;
	main_array[0] = 2;
	neighbors[0][0] = 1;
	neighbors[0][1] = 1;
	neighbors[1][0] = 1;
	neighbors[1][1] = 0;
	neighbors2[0][0] = 1;
	neighbors2[0][1] = 1;
	neighbors2[1][0] = 1;
	neighbors2[1][1] = 0;
	return 0;
}

int init3 (double main_array[], double shell, int neighbors[][N1], int neighbors2[][N2])
{
	main_array[3] = 0.0;
	main_array[4] = 0.0;
	main_array[5] = 0.0;
	main_array[9] = 0.0;
	main_array[10] = main_array[2]+main_array[8]+shell*2.;
	main_array[11] = 0.0;
	main_array[16] = main_array[2]+shell+
				(main_array[14]+shell)*(main_array[2]-main_array[8])/main_array[10];
	main_array[15] = sqrt((main_array[2]+main_array[14]+shell*2.0)
				* (main_array[2]+main_array[14]+shell*2.0)
				- main_array[16]*main_array[16]);
	main_array[17] = 0.0;
	main_array[0] = 3;
	neighbors[0][0] = 2;
	neighbors[0][1] = 1;
	neighbors[0][2] = 2;
	neighbors[1][0] = 2;
	neighbors[1][1] = 0;
	neighbors[1][2] = 2;
	neighbors[2][0] = 2;
	neighbors[2][1] = 0;
	neighbors[2][2] = 1;
	neighbors2[0][0] = 2;
	neighbors2[0][1] = 1;
	neighbors2[0][2] = 2;
	neighbors2[1][0] = 2;
	neighbors2[1][1] = 0;
	neighbors2[1][2] = 2;
	neighbors2[2][0] = 2;
	neighbors2[2][1] = 0;
	neighbors2[2][2] = 1;
	return 0;
}

int add1 (double main_array[], double shell, int neighbors[][N1])
{
	if (main_array[0] == number_of_dots)
	{
		return 0;
	}
	int i = main_array[0];
	//Добавляется сразу 2 частицы
	//Первая из добавляемых частиц пристраивается справа от самой правой частицы
	main_array[i*6+3] = main_array[(i-2)*6+3]+main_array[(i-2)*6+2]+main_array[i*6+2]+shell*2.;
	main_array[i*6+4] = 0.0;
	main_array[i*6+5] = 0.0;
	main_array[0]++;
	neighbors[i-2][0]++;
	neighbors[i-4][0]++;
	neighbors[i-6][0]++;
	neighbors[i-2][4] = i;
	neighbors[i-4][5] = i;
	neighbors[i-6][6] = i;
	neighbors[i][0] = 3;
	neighbors[i][1] = i-2;
	neighbors[i][2] = i-4;
	neighbors[i][3] = i-6;
	if (main_array[0] == number_of_dots)
	{
		return 0;
	}
	i = main_array[0];
	//Вторая из добавляемых частиц пристраивается слева от самой левой частицы
	main_array[i*6+3] = main_array[(i-2)*6+3]-main_array[(i-2)*6+2]-main_array[i*6+2]-shell*2.;
	main_array[i*6+4] = 0.0;
	main_array[i*6+5] = 0.0;
	main_array[0]++;
	neighbors[i-2][0]++;
	neighbors[i-4][0]++;
	neighbors[i-6][0]++;
	neighbors[i-2][4] = i;
	neighbors[i-4][5] = i;
	neighbors[i-6][6] = i;
	neighbors[i][0] = 3;
	neighbors[i][1] = i-2;
	neighbors[i][2] = i-4;
	neighbors[i][3] = i-6;
	return 1;
}

int sphere2 (int i, int j, double xyz[6], double main_array[], double shell);
int add2 (double main_array[], double average_radius, double shell, int neighbors[][N1], int neighbors2[][N2])
{
/*	std::cout << "Добавляем частицу с индексом: " << main_array[0] << std::endl;
	std::cout << "Координаты всех имеющихся частиц с соседями2:" << std::endl;
	for (int i = 0; i < main_array[0]; i++)
	{
		std::cout << i << ": (" << main_array[i*6+3] << ", " << main_array[i*6+4] << ") ";
		for (int j = 1; j <= neighbors2[i][0]; j++)
		{
			std::cout << neighbors2[i][j] << " ";
		}
		std::cout << std::endl;
	}
*/	if (main_array[0] == number_of_dots)
	{
		return 0;
	}
	//srand(time(0));
	int i, j, k ,l, l2, m, n, p, q = 0;
	int par1, par2;
	int ii, jj;
	double xyz[6], vars[50000][3];
	p = main_array[0]; //глобальный номер добавляемой частицы, он же число имеющихся частиц
	double rp = main_array[p*6+2];
	int s = (p>10000?round(5.*sqrt(p)):500);
	i = p>s?p-s:0;
	int step = p-i;
	while (1)
	{
		l = neighbors[i][0]; //число соседей у i-ой частицы из поверхностного слоя
		l2 = neighbors2[i][0]; //число соседей у i-ой частицы из поверхностного слоя
		for (j = 1; j <= l; j++)
		{
			m = neighbors[i][j]; //глобальные номера всех соседей i частицы из поверхностного слоя
			if (sphere2(m, i, xyz, main_array, shell))
			{
				par1 = 1;
				par2 = 1;
				for (jj = 1; jj <= l2; jj++) //пробегаем по всем частицам из ближайшего окружения и проверяем, не пересекает ли их новая сфера
				{
					ii = neighbors2[i][jj];
					double rii = main_array[ii*6+2];
					double xii = main_array[ii*6+3];
					double yii = main_array[ii*6+4];
/*					if (p == 3)
					{
						std::cout << "Специальный вывод: " << std::endl;
						std::cout << "xyz[0], xii, xyz[1], yii, rp, rii, shell: " << xyz[0] << ", " << xii << ", " << xyz[1] << ", " << yii << ", " << rp << ", " << rii << ", " << shell << std::endl;
						std::cout << "Контрольная разность 1: " << (xyz[0]-xii)*(xyz[0]-xii)+(xyz[1]-yii)*(xyz[1]-yii)-((rp+rii+2.*shell)*(rp+rii+2.*shell)) << std::endl;
						std::cout << "xyz[3], xii, xyz[4], yii, rp, rii, shell: " << xyz[3] << ", " << xii << ", " << xyz[4] << ", " << yii << ", " << rp << ", " << rii << ", " << shell << std::endl;
						std::cout << "Контрольная разность 2: " << (xyz[3]-xii)*(xyz[3]-xii)+(xyz[4]-yii)*(xyz[4]-yii)-((rp+rii+2.*shell)*(rp+rii+2.*shell)) << std::endl;
					}
*/					if (((xyz[0]-xii)*(xyz[0]-xii)+(xyz[1]-yii)*(xyz[1]-yii))<((rp+rii+2.*shell)*(rp+rii+2.*shell)-0.1))
					{
						par1 = 0;
					}
					if (((xyz[3]-xii)*(xyz[3]-xii)+(xyz[4]-yii)*(xyz[4]-yii))<((rp+rii+2.*shell)*(rp+rii+2.*shell)-0.1))
					{
						par2 = 0;
					}
				}
				if (par1)
				{
					vars[q][0] = xyz[0];
					vars[q][1] = xyz[1];
					vars[q][2] = xyz[2];
					q++;
				}
				if (par2)
				{
					vars[q][0] = xyz[3];
					vars[q][1] = xyz[4];
					vars[q][2] = xyz[5];
					q++;
				}
			}
		}
		step = step>1?step/2:1;
		if ((step == 1 && q) || q > 20)
		{
			break;
		}
		if (q == 0)
		{
			i += step;
		}
		if (q && i > step)
		{
			i -= step;
		}
	}
	ii = 0;
	for (i = 1; i < q; i++)
	{
		if ((vars[i][0]*vars[i][0]+vars[i][1]*vars[i][1])<(vars[ii][0]*vars[ii][0]+vars[ii][1]*vars[ii][1]))
		{
			ii = i;
		}
	}
	main_array[p*6+3] = vars[ii][0];
	main_array[p*6+4] = vars[ii][1];
	main_array[p*6+5] = vars[ii][2];
	main_array[0]++;
	double xp = main_array[p*6+3];
	double yp = main_array[p*6+4];
	if (p%1000 == 0)
	{
		//std::cout << "(x,y)[" << p << "] = ("<< xp << ", " << yp << "), q = " << q << std::endl;
	}
	for (i = 0; i < p; i++)
	{
		double xi = main_array[i*6+3];
		double yi = main_array[i*6+4];
		double ri = main_array[i*6+2];
		//В список ближайших соседей попадают частицы, удалённые от данной на расстояние,
		//не превышающее сумму радиусов наших двух частиц и средний диаметр с учётом оболочки
		double test_dist = rp+ri+shell*3.+average_radius;
		if ((xi-xp)*(xi-xp)+(yi-yp)*(yi-yp) <= test_dist * test_dist)
		{
			neighbors[i][0]++;
			neighbors[p][0]++;
			if (neighbors[i][0] > N1-1)
			{
				std::cout << "Error N1!" << std::endl;
				return 0;
			}
			neighbors[i][neighbors[i][0]] = p;
			neighbors[p][neighbors[p][0]] = i;
		}
		//test_dist = rp+ri+shell*4.+2*average_radius;
		if ((xi-xp)*(xi-xp)+(yi-yp)*(yi-yp) <= test_dist * test_dist)
		{
			neighbors2[i][0]++;
			neighbors2[p][0]++;
			if (neighbors2[i][0] > N2-1)
			{
				std::cout << "Error N2!" << std::endl;
				return 0;
			}
			neighbors2[i][neighbors2[i][0]] = p;
			neighbors2[p][neighbors2[p][0]] = i;
		}
	}
	return 1;
}

int sphere3 (int i, int j, int k, double xyz[6], double main_array[], double shell);
int add3 (double main_array[], double average_radius, double shell, int neighbors[][N1], int neighbors2[][N2])
{
/*	std::cout << "Добавляем частицу с индексом: " << main_array[0] << std::endl;
	std::cout << "Координаты всех имеющихся частиц с соседями2:" << std::endl;
	for (int i = 0; i < main_array[0]; i++)
	{
		std::cout << i << ": (" << main_array[i*6+3] << ", " << main_array[i*6+4] << ", "
			<< main_array[i*6+5] << ") ";
		std::cout << neighbors[i][0] << " -- " << neighbors2[i][0];
		for (int j = 1; j <= neighbors2[i][0]; j++)
		{
			std::cout << neighbors2[i][j] << " ";
		}
		std::cout << std::endl;
	}
*/	if (main_array[0] == number_of_dots)
	{
		return 0;
	}
	//srand(time(0));
	int i, j, k ,l, l2, m, n, p, q = 0;
	int par1, par2;
	int ii, jj;
	double xyz[6], vars[50000][3];
	p = main_array[0]; //глобальный номер добавляемой частицы, он же число имеющихся частиц
	double rp = main_array[p*6+2];
	int s = p>1000?round(5.*pow(p,0.67)):p;
	i = p>s?p-s:0;
	int step = p-i;
	while (1)
	{
		l = neighbors[i][0]; //малое число соседей у i-ой частицы из поверхностного слоя
		l2 = neighbors2[i][0]; //большое число соседей у i-ой частицы из поверхностного слоя
		for (j = 2; j <= l; j++)
		{
			m = neighbors[i][j]; //глобальные номера всех соседей i частицы из поверхностного слоя
			for (k = 1; k < j; k++)
			{
				n = neighbors[i][k];
				if (sphere3(i, m, n, xyz, main_array, shell))
				{
					par1 = 1;
					par2 = 1;
					for (jj = 1; jj <= l2; jj++) //пробегаем по всем частицам из ближайшего окружения и проверяем, не пересекает ли их новая сфера
					{
						ii = neighbors2[i][jj];
						double rii = main_array[ii*6+2];
						double xii = main_array[ii*6+3];
						double yii = main_array[ii*6+4];
						double zii = main_array[ii*6+5];
/*						if (p == 3)
						{
							std::cout << "Специальный вывод: " << std::endl;
							std::cout << "xyz[0], xii, xyz[1], yii, rp, rii, shell: " << xyz[0] << ", " << xii << ", " << xyz[1] << ", " << yii << ", " << rp << ", " << rii << ", " << shell << std::endl;
							std::cout << "Контрольная разность 1: " << (xyz[0]-xii)*(xyz[0]-xii)+(xyz[1]-yii)*(xyz[1]-yii)-((rp+rii+2.*shell)*(rp+rii+2.*shell)) << std::endl;
							std::cout << "xyz[3], xii, xyz[4], yii, rp, rii, shell: " << xyz[3] << ", " << xii << ", " << xyz[4] << ", " << yii << ", " << rp << ", " << rii << ", " << shell << std::endl;
							std::cout << "Контрольная разность 2: " << (xyz[3]-xii)*(xyz[3]-xii)+(xyz[4]-yii)*(xyz[4]-yii)-((rp+rii+2.*shell)*(rp+rii+2.*shell)) << std::endl;
						}
	*/					if ((xyz[0]-xii)*(xyz[0]-xii)+(xyz[1]-yii)*(xyz[1]-yii)+(xyz[2]-zii)*(xyz[2]-zii)<(rp+rii+2.*shell)*(rp+rii+2.*shell)-0.1)
						{
							par1 = 0;
						}
						if ((xyz[3]-xii)*(xyz[3]-xii)+(xyz[4]-yii)*(xyz[4]-yii)+(xyz[5]-zii)*(xyz[5]-zii)<(rp+rii+2.*shell)*(rp+rii+2.*shell)-0.1)
						{
							par2 = 0;
						}
					}
					if (par1)
					{
						vars[q][0] = xyz[0];
						vars[q][1] = xyz[1];
						vars[q][2] = xyz[2];
						q++;
					}
					if (par2)
					{
						vars[q][0] = xyz[3];
						vars[q][1] = xyz[4];
						vars[q][2] = xyz[5];
						q++;
					}
				}
			}
		}
		step = step>1?step/2:1;
		if ((step == 1 && q) || q > 20)
		{
			break;
		}
		if (q == 0)
		{
			i += step;
		}
		if (q && i > step)
		{
			i -= step;
		}
	}
	ii = 0;
	for (i = 1; i < q; i++)
	{
		if ((vars[i][0]*vars[i][0]+vars[i][1]*vars[i][1])+vars[i][2]*vars[i][2]
			<(vars[ii][0]*vars[ii][0]+vars[ii][1]*vars[ii][1]+vars[ii][2]*vars[ii][2]))
		{
			ii = i;
		}
	}
	main_array[p*6+3] = vars[ii][0];
	main_array[p*6+4] = vars[ii][1];
	main_array[p*6+5] = vars[ii][2];
	main_array[0]++;
	double xp = main_array[p*6+3];
	double yp = main_array[p*6+4];
	double zp = main_array[p*6+5];
	if (p%100 == 0)
	{
		//std::cout << "(x,y,z)[" << p << "] = ("<< xp << ", " << yp << ", " << zp << "), q = " << q << std::endl;
	}
	for (i = 0; i < p; i++)
	{
		double xi = main_array[i*6+3];
		double yi = main_array[i*6+4];
		double zi = main_array[i*6+5];
		double ri = main_array[i*6+2];
		double test_dist = rp+ri+shell*3.+average_radius;
		if ((xi-xp)*(xi-xp)+(yi-yp)*(yi-yp)+(zi-zp)*(zi-zp) <= test_dist*test_dist)
		{
			neighbors[i][0]++;
			neighbors[p][0]++;
			if (neighbors[i][0] > N1-1)
			{
				std::cout << "Error N1!" << std::endl;
				std::cout << "p = " << p << std::endl;
				std::cout << "N1 = " << N1 << std::endl;
				std::cout << neighbors[i][0] << std::endl;
				std::cout << "average_radius = " << average_radius << std::endl;
				std::cout << "shell = " << shell << std::endl;
				std::cout << "S = " << (xi-xp)*(xi-xp)+(yi-yp)*(yi-yp)+(zi-zp)*(zi-zp) << std::endl;
				std::cout << "xi = " << xi << std::endl;
				std::cout << "yi = " << yi << std::endl;
				std::cout << "zi = " << zi << std::endl;
				std::cout << "xp = " << xp << std::endl;
				std::cout << "yp = " << yp << std::endl;
				std::cout << "zp = " << zp << std::endl;
				exit(1);
			}
			neighbors[i][neighbors[i][0]] = p;
			neighbors[p][neighbors[p][0]] = i;
		}
		//test_dist = rp+ri+shell*4.+average_radius*2;
		if ((xi-xp)*(xi-xp)+(yi-yp)*(yi-yp)+(zi-zp)*(zi-zp) <= test_dist*test_dist)
		{
			neighbors2[i][0]++;
			neighbors2[p][0]++;
			if (neighbors2[i][0] > N2-1)
			{
				std::cout << "Error N2!" << std::endl;
				exit(1);
			}
			neighbors2[i][neighbors2[i][0]] = p;
			neighbors2[p][neighbors2[p][0]] = i;
		}
	}
	return 1;
}

int sphere2 (int i, int j, double xyz[], double main_array[], double shell)
{
	if (main_array[0] == number_of_dots)
	{
		return 0;
	}
	double a, b, znam, sum, p, q, aa, bb, cc, dd;
	int k = main_array[0];
	double ri = main_array[i*6+2];
	double rj = main_array[j*6+2];
	double rk = main_array[k*6+2];
	double xi = main_array[i*6+3];
	double xj = main_array[j*6+3];
	double yi = main_array[i*6+4];
	double yj = main_array[j*6+4];

	a = ri + rk + shell*2.0;
	b = rj + rk + shell*2.0;
	znam = yj - yi;
	sum = a*a - b*b + xj*xj - xi*xi + yj*yj - yi*yi;
	if (fabs(znam) < 1.e-5)
	{
		xyz[0] = 0.5*sum/(xj-xi);
		xyz[3] = xyz[0];
		dd = a*a - (xyz[0]-xi)*(xyz[0]-xi);
		if (dd < -1e-5)
		{
			return 0;
		}
		if (dd < 0)
		{
			dd = 0.;
		}
		dd = sqrt(dd);
		xyz[1] = yi+dd;
		xyz[4] = yi-dd;
	}
	else
	{
		p = 0.5*sum/znam - yi;
		q = (xi-xj)/znam;
		aa = 1 + q*q;
		bb = 2*p*q - 2*xi;
		cc = xi*xi + p*p - a*a;
		dd = bb*bb - 4*aa*cc;
		if (dd < -1e-5)
		{
			return 0;
		}
		if (dd < 0)
		{
			dd = 0.;
		}
		dd = sqrt(dd);
		xyz[0] = 0.5*(-bb-dd)/aa;
		xyz[3] = 0.5*(-bb+dd)/aa;
		xyz[1] = 0.5*sum/znam + (xi-xj)/znam*xyz[0];
		xyz[4] = 0.5*sum/znam + (xi-xj)/znam*xyz[3];
	}
	xyz[2] = 0.;
	xyz[5] = 0.;
	return 1;
}

int sphere3 (int i, int j, int k, double xyz[6], double main_array[], double shell)
{
	if (main_array[0] == number_of_dots)
	{
		return 0;
	}
	double a, b, c, m, n, p, q, s, t, znam, aa, bb, cc, dd;
	int l = main_array[0];
	double ri = main_array[i*6+2];
	double rj = main_array[j*6+2];
	double rk = main_array[k*6+2];
	double rl = main_array[l*6+2];
	double xi = main_array[i*6+3];
	double xj = main_array[j*6+3];
	double xk = main_array[k*6+3];
	double yi = main_array[i*6+4];
	double yj = main_array[j*6+4];
	double yk = main_array[k*6+4];
	double zi = main_array[i*6+5];
	double zj = main_array[j*6+5];
	double zk = main_array[k*6+5];

	a = ri + rl + shell*2.;
	b = rj + rl + shell*2.;
	c = rk + rl + shell*2.;
	m = a*a - (xi*xi+yi*yi+zi*zi) - b*b + (xj*xj+yj*yj+zj*zj);
	n = a*a - (xi*xi+yi*yi+zi*zi) - c*c + (xk*xk+yk*yk+zk*zk);
	znam = 2.*(xi-xj)*(yi-yk)-2.*(xi-xk)*(yi-yj);

	if (fabs(znam) < 1.e-5)
	{
		return 0;
	}

	p = 2.*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))/znam;
	q = ((yi-yj)*n-(yi-yk)*m)/znam;
	s = 2.*((xi-xk)*(zi-zj)-(xi-xj)*(zi-zk))/znam;
	t = ((xi-xk)*m-(xi-xj)*n)/znam;

	aa = p*p+s*s+1.;
	bb = 2.*(p*q+s*t-p*xi-s*yi-zi);
	cc = q*q+t*t-2.*(q*xi+t*yi)+xi*xi+yi*yi+zi*zi-a*a;
	dd = bb*bb-4.*aa*cc;

	if (dd < -1e-5)
	{
		return 0;
	}
	if (dd < 0)
	{
		dd = 0.;
	}

	dd = sqrt(dd);
	xyz[2] = 0.5*(-bb+dd)/aa;
	xyz[5] = 0.5*(-bb-dd)/aa;
	xyz[0] = p*xyz[2]+q;
	xyz[3] = p*xyz[5]+q;
	xyz[1] = s*xyz[2]+t;
	xyz[4] = s*xyz[5]+t;

	return 1;
}

int read_rf6 (double rf6[], double cdse)
{
	FILE *finp;
	double value1, value2, value3;
	int i = 0;
	if (cdse > 0.5)
	{
		finp = fopen ("rf_cdse.txt", "rt");
	}
	else
	{
		finp = fopen ("rf_inp.txt", "rt");
	}

	while (1)
	{
		if (fscanf (finp, "%lf %lf %lf\n", &value1, &value2, &value3) != 3)
		{
			break;
		}
		rf6[i] = value1;
		rf6[i+1] = value2;
		rf6[i+2] = value3;
		i+=3;
	}

	fclose (finp);
	return 0;
}

int init_probabilities (double main_array[], int p_initial)
{
	main_array[1] = p_initial;
	for (int i = 0; i < main_array[0]; i++)
	{
		main_array[6*i+6] = 0;
	}
	main_array[6*p_initial+6] = 1;
	return 0;
}

int init_black (double main_array[], double vphi)
{
	for (int i = 0; i < main_array[0]; i++)
	{
		if ((double)rand()/RAND_MAX > vphi)
		{
			//эти частицы красим в чёрный цвет
			main_array[6*i+7] = 0;
		}
		else
		{
			//эти частицы красим в белый цвет
			main_array[6*i+7] = 1;
		}
	}
	return 0;
}

double step_probabilities (double main_array[], int neighbors[][N1], double rf6[], double k0, double phi, double cdse);
void go_probabilities (double main_array[], int neighbors[][N1], double rf6[], double all_ksi2[], double k0, double phi, double cdse, double all_delta_peak_fg[])
{
	//Вычисляем положение пика спектра поглощения всего нанокластера
	//Суммируем пики поглощения всех НЧ и усредняем
	//Пик поглощения НЧ находим исходя из её радиуса по формуле
	//Пик вычисляется 1 раз для данного набора радиусов
	//Размерность нанометры
	double peak_f = 0.;
	for (int i = 2; i < number_of_dots*6+2; i+=6)
	{
		peak_f += 1035/(1 + 2.414/2/main_array[i]);
	}
	peak_f /= number_of_dots;

	//Вычисляем положение пика спектра люминесценции нанокластера
	//В каждый момент времени (из 5000)! возбуждение находится на одной НЧ, поэтому
	//спектр нанокластера определяется радиусом только этой частицы, поэтому нужен массив
	double peak_g[5000];//!
	double std_peak_g[10000], std_delta_peak_fg[10000];

	double ksi2[5000], time[5000], time_step;//!
	int i, j = 1, np;
	ksi2[0] = 0.;
	time[0] = 0.;
	for (i = 1; i < 5000; i++)//!
	{
		ksi2[i] = 0.;
		time[i] = 1000000.+i;
	}
	double std_time[10000], std_ksi2[10000], std_step = 5.e-12;
	for (i = 0; i < 10000; i++)
	{
		std_time[i] = i*std_step;
		std_ksi2[i] = 0.;
	}
	for (i = 1; i < 5000; i++)//!
	{
		time_step = step_probabilities (main_array, neighbors, rf6, k0, phi, cdse);
		//if (i < 5)
			//std::cout << time_step << std::endl;
		if (time_step < 1e-20)
		{
			while (i < 5000)//!
			{
				ksi2[i] = ksi2[i-1];
				i++;
			}
			break;
		}
		time[i] = time[i-1] + time_step;
		np = main_array[1];//Узнаём номер частицы, на которой сейчас возбуждение
		ksi2[i] = main_array[6*np+3]*main_array[6*np+3]//Вычисляем расстояние частицы 
			+ main_array[6*np+4]*main_array[6*np+4]//от первой частицы
			+ main_array[6*np+5]*main_array[6*np+5];//Размерность кв.нанометры
		//Здесь вычисляем пик люминесценции НЧ, на которой сейчас возбуждение
		//Размерность нанометры
		peak_g[i] = 1/((1 + 2.414/2/main_array[6*np+2]) / 1035 - 0.000039);
	}
	//time[999] = 10000000.;
	for (i = 1; i < 10000; i++)
	{
		while (std_time[i] > time[j])
		{
			j++;
			if (j == 5000)//!
			{
				//std::cout << std_time[i] << "; " << time[j-1] << std::endl;
				break;
			}
		}
		if (j == 5000)//!
		{
			break;
		}
		std_ksi2[i] = ksi2[j-1]+(std_time[i]-time[j-1])*(ksi2[j]-ksi2[j-1])/(time[j]-time[j-1]);
		all_ksi2[i] += std_ksi2[i];

		std_peak_g[i] = peak_g[j-1]
			+(std_time[i]-time[j-1])*(peak_g[j]-peak_g[j-1])/(time[j]-time[j-1]);
		std_delta_peak_fg[i] = std_peak_g[i] - peak_f;
		all_delta_peak_fg[i] += std_delta_peak_fg[i];
		//std::cout << "time_j = " << time[j] << std::endl;
	}
}

double step_probabilities (double main_array[], int neighbors[][N1], double rf6[], double k0, double phi, double cdse)
{
	int n = main_array[1]; //Узнаём номер частицы, на которой сейчас возбуждение
	if (main_array[6*n+7] == 0) //Если эта частица чёрная, сразу выходим из функции
	{
		return 0;
	}
	int m = neighbors[n][0]; //Узнаём число ближайших соседей этой частицы
	double k_ij[N1], k_ji[N1];
	int i, j;
	double r0 = main_array[6*n+2]; //Радиус частицы, на которой возбуждение
	double ri; //Радиусы всех соседей
	double rij6; //Расстояние до каждого соседа сразу в 6 степени
	double x0, y0, z0, xi, yi, zi;
	x0 = main_array[6*n+3]; //Координаты частицы, на которой возбуждение
	y0 = main_array[6*n+4];
	z0 = main_array[6*n+5];
	double rf6_ij, rf6_ji; //Радиусы фёрстера в 6 степени
	double rmin, rmax;
	int nsteps;
	if (cdse > 0.5)
	{
		rmin = 1.080611;
		rmax = 2.415654;
		nsteps = 400;
	}
	else
	{
		rmin = 1.487189;
		rmax = 1.996367;
		nsteps = 200;
	}
	double delta = (rmax - rmin) / nsteps;
	int n1, n2, s;
	double sum_ij = 0., sum_ji = 0.;
	for (i = 0; i < m; i++) //Перебираем всех соседей в списке ближайших соседей
	{
		j = neighbors[n][i+1]; //Номер текущего соседа в глобальном списке
		ri = main_array[6*j+2]; //Радиус текущего соседа
		xi = main_array[6*j+3]; //Координаты текущего соседа
		yi = main_array[6*j+4];
		zi = main_array[6*j+5];
		rij6 = (x0-xi)*(x0-xi)+(y0-yi)*(y0-yi)+(z0-zi)*(z0-zi);
		rij6 = rij6 * rij6 * rij6;
		n1 = (r0 - rmin) / delta;
		n2 = (ri - rmin) / delta;
		s = (nsteps+1)*n1+n2;
		rf6_ij = rf6[3*s+2];
		s = (nsteps+1)*n2+n1;
		rf6_ji = rf6[3*s+2];
		k_ij[i] = phi*rf6_ij/rij6;
		k_ji[i] = phi*rf6_ji/rij6;
		sum_ij += k_ij[i];
		sum_ji += k_ji[i];
	}
	sum_ij += k0 * 5.e7;
	sum_ji += k0 * 5.e7;
	double tmp_sum_ij = 0.;
	double tmp_sum_ji = 0.;
	double p = (double)rand()/RAND_MAX;
	for (i = 0; i <= m; i++)
	{
		if (i == m)
		{
			return 0.;
		}
		tmp_sum_ij += k_ij[i];
		if (p <= tmp_sum_ij/sum_ij)
		{
			j = neighbors[n][i+1];
			main_array[6*j+6] = 1;
			main_array[1] = j;
			main_array[6*n+6] = 0;
			break;
		}
	}
	if (sum_ij)
	{
		//std::cout << "time = " << 1./sum_ij << std::endl;
		return 1./sum_ij;
	}
	else
	{
		return 0;
	}
}

///////////////////////////////////////////////////////////////////
//Набор функций для получения файла с радиусами Фёрстера для CdSe//
///////////////////////////////////////////////////////////////////
int read_spectrum (char fname[], double v[401][1502])
{
	FILE *finp;
	double value1, value2;
	int i = 0, j = 0;
	while (v[i][0] > 1)
	{
		i++;
	}
	if ((finp = fopen (fname, "rt")) != NULL)
	{
		//std::cout << "File " << fname << " has been opened successfully!" << std::endl;
	}
	else
	{
		//std::cout << "File " << fname << " couldn't be opened. Program finished with error." << std::endl;
		return 1;
	}
	while (1)
	{
		if (fscanf (finp, "%lf %lf\n", &value1, &value2) != 2)
		{
			break;
		}
		v[i][j] = value2;
		//std::cout << values[i] << std::endl;
		j++;
	}

	fclose (finp);
	return 0;
}

void aToSigma (double vA[401][1502])
{
	int i, j, n;
	double diameter, peak_ext_value, coeff;
	for (i = 0; i < 8; i++)
	{
		//Размерность диаметра и длины волны, приходящейся на пик поглощения, в нанометрах
		//vA[i][0] -- положение пика поглощения для восьми исходных спектров CdSe
		diameter = (1.6122e-9)*pow(vA[i][0],4) - (2.6575e-6)*pow(vA[i][0],3) + (1.6242e-3)*pow(vA[i][0],2) - 0.4277*vA[i][0] + 41.57;
        	peak_ext_value = 5857*pow(diameter, 2.65);
		n = (vA[i][0]-400)*5+1;
		coeff = peak_ext_value / vA[i][n] * log(10) * 0.1661e-23;
		for (j = 1; j < 1502; j++)
		{
			//Размерность сечения поглощения в квадратных метрах
			//Спектр сечения поглощения от длины волны
			vA[i][j] *= coeff;
		}
	}
}

void bToNormilize (double vA[401][1502], double vB[401][1502])
{
	int i, j;
	double Sum = 0.;
	for (i = 0; i < 8; i++)
	{
		for (j = 1; j < 1502; j++)
		{
			Sum += vB[i][j];
		}
		//Шаг интегрирования в метрах
		Sum *= 2.e-10;
		for (j = 1; j < 1502; j++)
		{
			//Размерность нормированного спектра люминесценции в обратных метрах
			vB[i][j] /= Sum;
		}
	}
}

void spreadAB (double vA[401][1502], double vB[401][1502], int p[8])
{
	int i;
	double d0, di, d7;
	d0 = (1.6122e-9)*pow(vA[0][0],4) - (2.6575e-6)*pow(vA[0][0],3) + (1.6242e-3)*pow(vA[0][0],2) - 0.4277*vA[0][0] + 41.57;
	d7 = (1.6122e-9)*pow(vA[7][0],4) - (2.6575e-6)*pow(vA[7][0],3) + (1.6242e-3)*pow(vA[7][0],2) - 0.4277*vA[7][0] + 41.57;
	for (i = 0; i < 8; i++)
	{
		di = (1.6122e-9)*pow(vA[i][0],4) - (2.6575e-6)*pow(vA[i][0],3) + (1.6242e-3)*pow(vA[i][0],2) - 0.4277*vA[i][0] + 41.57;
		p[i] = 400*(di-d0)/(d7-d0);
	}
	for (i = 0; i < 1502; i++)
	{
		vA[p[7]][i] = vA[7][i];
		vA[7][i] = 0;
		vA[p[6]][i] = vA[6][i];
		vA[6][i] = 0;
		vA[p[5]][i] = vA[5][i];
		vA[5][i] = 0;
		vA[p[4]][i] = vA[4][i];
		vA[4][i] = 0;
		vA[p[3]][i] = vA[3][i];
		vA[3][i] = 0;
		vA[p[2]][i] = vA[2][i];
		vA[2][i] = 0;
		vA[p[1]][i] = vA[1][i];
		vA[1][i] = 0;
		vB[p[7]][i] = vB[7][i];
		vB[7][i] = 0;
		vB[p[6]][i] = vB[6][i];
		vB[6][i] = 0;
		vB[p[5]][i] = vB[5][i];
		vB[5][i] = 0;
		vB[p[4]][i] = vB[4][i];
		vB[4][i] = 0;
		vB[p[3]][i] = vB[3][i];
		vB[3][i] = 0;
		vB[p[2]][i] = vB[2][i];
		vB[2][i] = 0;
		vB[p[1]][i] = vB[1][i];
		vB[1][i] = 0;
	}
}

void fillAB (double vA[401][1502], double vB[401][1502], int p[8])
{
	int i, j, k, j1, j2;
	int t1A, t1B, t, t2A, t2B;
	double w1, w2;
	for (k = 0; k < 7; k++)
	{
		t1A = (vA[p[k]][0]-400)*5+1;
		t1B = (vB[p[k]][0]-400)*5+1;
		t2A = (vA[p[k+1]][0]-400)*5+1;
		t2B = (vB[p[k+1]][0]-400)*5+1;
		for (i = p[k]+1; i < p[k+1]; i++)
		{
			w1 = (p[k+1]-(double)i)/(p[k+1]-p[k]);
			w2 = 1-w1;
			vA[i][0] = round((vA[p[k]][0]*w1+vA[p[k+1]][0]*w2)*5)/5;
			vB[i][0] = round((vB[p[k]][0]*w1+vB[p[k+1]][0]*w2)*5)/5;
			t = (vA[i][0]-400)*5+1;
			w1 = (t2A-(double)t)/(t2A-t1A);
			w2 = 1-w1;
			for (j = 1; j < 1502; j++)
			{
				if (j+(t1A-t) < 1)
				{
					j1 = 1;
				}
				else if (j+(t1A-t) > 1501)
				{
					j1 = 1501;
				}
				else
				{
					j1 = j+(t1A-t);
				}
				if (j+(t2A-t) < 1)
				{
					j2 = 1;
				}
				else if (j+(t2A-t) > 1501)
				{
					j2 = 1501;
				}
				else
				{
					j2 = j+(t2A-t);
				}
				vA[i][j] = vA[p[k]][j1]*w1+vA[p[k+1]][j2]*w2;
			}
			t = (vB[i][0]-400)*5+1;
			w1 = (t2B-(double)t)/(t2B-t1B);
			w2 = 1-w1;
			for (j = 1; j < 1502; j++)
			{
				if (j+(t1B-t) < 1)
				{
					j1 = 1;
				}
				else if (j+(t1B-t) > 1501)
				{
					j1 = 1501;
				}
				else
				{
					j1 = j+(t1B-t);
				}
				if (j+(t2B-t) < 1)
				{
					j2 = 1;
				}
				else if (j+(t2B-t) > 1501)
				{
					j2 = 1501;
				}
				else
				{
					j2 = j+(t2B-t);
				}
				vB[i][j] = vB[p[k]][j1]*w1+vB[p[k+1]][j2]*w2;
			}
		}
	}
}

void printAB (double vA[401][1502], double vB[401][1502])
{
	FILE *vAfile, *vBfile;
	int i, j;
	vAfile = fopen ("vA.txt", "wt");
	vBfile = fopen ("vB.txt", "wt");
	for (j = 0; j < 1502; j++)
	{
		for (i = 0; i < 401; i++)
		{
			if (j == 0)
			{
				fprintf (vAfile, "%8.1lf", vA[i][j]);
				fprintf (vBfile, "  %8.1lf", vB[i][j]);
			}
			else
			{
				fprintf (vAfile, " %1.5lf", vA[i][j]*1.e20);
				fprintf (vBfile, " %9.5lf", vB[i][j]);
			}
		}
		fprintf (vAfile, "\n");
		fprintf (vBfile, "\n");
	}
	fclose (vAfile);
	fclose (vBfile);
}

void calcRF (double vA[401][1502], double vB[401][1502])
{
	FILE *RFfile;
	int i, j, k;
	double Integral;
	double r[401];
	int m;
	RFfile = fopen ("rf_cdse.txt", "wt");

	for (i = 0; i < 401; i++)
	{
		r[i] = ((1.6122e-9)*pow(vA[i][0],4) - (2.6575e-6)*pow(vA[i][0],3) + (1.6242e-3)*pow(vA[i][0],2) - 0.4277*vA[i][0] + 41.57)*0.5;
	}
	//Пробегаем по всем донорам i
	for (i = 0; i < 401; i++)
	{
		//Пробегаем по всем акцепторам j при фиксированном доноре i
		for (j = 0; j < 401; j++)
		{
			Integral = 0.;
			for (k = 1; k < 1502; k++)
			{
				Integral += vB[i][k]*vA[j][k]*pow((0.2*(k-1)+400)*1.e-9,4);
			}
			Integral *= 0.2*1.e-9*9*2./3./128./pow(4*atan(1),5)/pow(1.5,4);
			//Умножаем на 1/tau
			Integral *= 5.e7;
			if (Integral)
			{
				m = floor(log10(Integral));
			}
			else
			{
				m = 0;
			}
			//Умножаем на 1е54, чтобы потом можно было использовать R^6 в нм
			fprintf (RFfile, "%lf %lf %.15lfe+%d\n", r[i], r[j], Integral/pow(10,m), m+54);
		}
	}
	//Конечные радиусы Фёрстера домножаем на 1е46 для удобного отображения на печати
	fclose (RFfile);
	return;
}

////////////////////////////
//Функции для работы с InP//
////////////////////////////
int read_fg (double f[2251], double g[2251])
{
	FILE *finp;
	//f -- это значения поглощения, value2
	//g -- это значения люминесценции, value3
	double value1, value2, value3;

	int i;
	for (i = 0; i < 2251; i++)
	{
		f[i] = 0.;
		g[i] = 0.;
	}
	i = 1000;

	finp = fopen ("fg.txt", "rt");

	while (1)
	{
		if (fscanf (finp, "%lf %lf %lf\n", &value1, &value2, &value3) != 3)
		{
			break;
		}
		f[i] = value2;
		g[i] = value3;
		i++;
	}

	fclose (finp);
	return 0;
}

int read_rho (double rho[665])
{
	FILE *finp;
	double value1, value2;
	int i = 0;
	if ((finp = fopen ("rho.txt", "rt")) != NULL)
	{
		//std::cout << "File " << fname << " has been opened successfully!" << std::endl;
	}
	else
	{
		//std::cout << "File " << fname << " couldn't be opened. Program finished with error." << std::endl;
		return 1;
	}
	while (1)
	{
		if (fscanf (finp, "%lf %lf\n", &value1, &value2) != 2)
		{
			break;
		}
		rho[i] = value2;
		i++;
	}
	//Нормируем спектр rho, единица измерения метры
	double sum = 0.;
	for (i = 0; i < 665; i++)
	{
		sum += rho[i];
	}
	//Переводим шаг интегрирования из 1 миллионной доли обратного нанометра в 1000 обратных метров
	//Поэтому домножаем не на 1е-6, а на 1е3, т.к. в обратных метрах шаг интегрирования равен 1000
	sum *= 1.e3;
	for (i = 0; i < 665; i++)
	{
		rho[i] /= sum;
	}

	fclose (finp);
	return 0;
}

int read_eas (double eas[2160])
{
	FILE *finp;
	double value1, value2;
	int i = 0;
	if ((finp = fopen ("eas.txt", "rt")) != NULL)
	{
		//std::cout << "File eas.txt has been opened successfully!" << std::endl;
	}
	else
	{
		//std::cout << "File eas.txt couldn't be opened. Program finished with error." << std::endl;
		return 1;
	}
	while (1)
	{
		if (fscanf (finp, "%lf %lf\n", &value1, &value2) != 2)
		{
			break;
		}
		eas[i] = value2;
		i++;
	}

	fclose (finp);
	return 0;
}

//С новыми файлами Товстуна эта функция не нужна. Присланный f уже является сечением поглощения 
void fToSigma (double f[2251], double rho[665], double eas[2160])
{
	int i, j;
	double svertka_f_rho[2251];

	for (i = 0; i < 2251; i++)
	{
		svertka_f_rho[i] = 0.;
		for (j = 0; j < 665; j++)
		{
			if ((i+300-j >= 0) && (i+300-j < 2251))//!!!!!!!!!!!!!!!!!!!!!!!
			{
				svertka_f_rho[i] += f[i+300-j]*rho[j]*1.e3;//!!!!!!!!!!!!!!!!!
			}
		}
		std::cout << i << " " << svertka_f_rho[i] << std::endl;
	}
	for (i = 0; i < 2160; i++)
	{
		eas[i] *= (37890*pow((2.414/(1035./560.-1.)),3)*0.6156/3.9296)*log(10)/6.022e24;//!!!!!
		//std::cout << 300+i << " " << eas[i] << std::endl;
	}
	for (i = 1000; i < 1251; i++)
	{
		f[i] *= eas[260]/svertka_f_rho[1363];//!!!!!!!!!!!!!!!!!!!!
		//std::cout << i << " " << f[i] << std::endl;
	}
}

void gToNormilize (double g[2251])
{
	int i;
	double sum = 0.;
	for (i = 1000; i < 1251; i++)
	{
		sum += 0.5e3*(g[i]+g[i+1]);
	}
	for (i = 1000; i < 1251; i++)
	{
		//Размерность нормированного спектра люминесценции в метрах
		g[i] /= sum;
	}
}

void calcRF_inp (double f[2251], double g[2251])
{
	FILE *RFfile;
	int i, j, k;
	double Integral;
	double r[201];
	int m;
	RFfile = fopen ("rf_inp.txt", "wt");

	for (i = 0; i < 201; i++)
	{
		r[200-i] = 0.5*(2.414/(1035*(0.0015173364+(33+i)*1.e-6)-1.));
	}
	//Пробегаем по всем донорам i
	for (i = 0; i < 201; i++)
	{
		//Пробегаем по всем акцепторам j при фиксированном доноре i
		for (j = 0; j < 201; j++)
		{
			Integral = 0.;
			for (k = 900; k < 1351; k++)
			{
				Integral += f[k+j-100]*g[k+i-100]*pow((0.0015173364+(k-1000)*1.e-6)*1.e9,-4)*1.e3;//!!!!!!!!!!!!!!!!!!!!
			}
			//std::cout << Integral << std::endl;
			Integral *= 9*2./3./128./pow(4*atan(1),5)/pow(1.5,4);
			//Умножаем на 1/tau
			Integral *= 5.e7;
			if (Integral)
			{
				m = floor(log10(Integral));
			}
			else
			{
				m = 0;
			}
			//Умножаем на 1е54, чтобы потом можно было использовать R^6 в нм
			fprintf (RFfile, "%lf %lf %.15lfe+%d\n", r[i], r[j], Integral/pow(10,m), m+54);
		}
	}
	//Конечные радиусы Фёрстера домножаем на 1е46 для удобного отображения на печати
	fclose (RFfile);
	return;
}
