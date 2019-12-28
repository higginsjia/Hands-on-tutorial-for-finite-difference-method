#ifndef _TEMP_SOLVER_IMPLICT_
#define _TEMP_SOLVER_IMPLICT_

#include <algorithm>
#include "math.h"
#include "stdio.h"
#include "time.h"
#include "omp.h"
#pragma warning( disable : 4996 4305 4244)
using namespace std;


class implictSolver
{
public:
	implictSolver(){};
	~implictSolver(){};


	int nx;
	int ny;
	int cellNum;

	float Tbound = 25.0f;

	float* phi_cur;
	float* phi_nex;
	bool* bstatus;

	float* vB;//
	float* vX;//t=t
	float* vXn;//t=t+dt

	float*coeff_a, *coeff_b, *coeff_c, *coeff_d, *coeff_e;

	void set_init_value()
	{
		nx = 200;
		ny = 200;
		cellNum = nx*ny;

		printf("nx=%d ny=%d N=%d\n", nx, ny, cellNum);

	}


	void run()
	{
		clock_t sta = clock();

		const int n_out = 1000;
		for (int ITER = 0; ITER < 10000; ++ITER)
		{
			if (ITER%n_out==0)
				printf("T=%5d \n",ITER);
			
			if (jacobiSolver(1000, 1.0e-6f, phi_cur, phi_nex) < 0)
			{
				printf("Not converge!\n");
				break;
			}
			
			for (int cur = 0; cur < cellNum; ++cur)
			{
				phi_cur[cur] = phi_nex[cur];
			}


			if (ITER % n_out== 0)
			{
				int nframe = ITER / n_out;
				writeTec(nx, ny, nframe, phi_nex);
			}
		}

		float t_used = (float)(clock() - sta) / CLOCKS_PER_SEC;

		printf("time used=%.3f (s)\n", t_used);
	}


	int jacobiSolver(int max_iter,float eps, float *in_var,float* out_var)
	{
		//AX=B
		const float k0 = 0.1f;//dt/dx^2

#pragma omp parallel for
		for (int x = 0; x < nx; ++x)
		for (int y = 0; y < ny; ++y)
		{
			const int cur = cellIndex(x, y);
			bool bound = (x == 0) || (x == nx - 1) || (y == 0) || (y == ny - 1);

			if (bound)
			{
				vB[cur] = Tbound;
				coeff_a[cur] = 0.0f;	coeff_b[cur] = 0.0f;coeff_c[cur] = 0.0f; coeff_d[cur] = 0.0f;
				coeff_e[cur] = 1.0f;
				vX[cur] = Tbound;
			}
			else
			{
				vB[cur] = in_var[cur];
				coeff_a[cur] = k0;	coeff_b[cur] = k0;	coeff_c[cur] = k0;	coeff_d[cur] = k0;
				coeff_e[cur] = 1.0f + 4.0f*k0;
				vX[cur] = in_var[cur];
			}
		}

		//int max_iter = 100;
		for (int iter = 0; iter <= max_iter; ++iter)
		{
#pragma omp parallel for
			for (int cur = 0; cur < cellNum; ++cur)
			{
				if (!bstatus[cur])
				{
					vXn[cur] =
						(vB[cur] + coeff_a[cur] * vX[cur - ny] + coeff_b[cur] * vX[cur - 1] +
						coeff_c[cur] * vX[cur + 1] + coeff_d[cur] * vX[cur + ny]) / coeff_e[cur];//
				}
				else
				{
					vXn[cur] = vB[cur];
				}
			}

			float err_sum = 0.0f;
			for (int cur = 0; cur < cellNum; ++cur)
			{
				err_sum += fabs(vXn[cur] - vX[cur]);	
			}
			err_sum /= (float)cellNum;

			//printf("IT=%6d Err=%e\n", iter, err_sum);

			if (err_sum < eps)
			{
				//printf("IT=%6d Err=%e Converged!\n", iter, err_sum);		
				for (int cur = 0; cur < cellNum; ++cur)
					out_var[cur] = vXn[cur];

				return 0;
			}

#pragma omp parallel for
			for (int cur = 0; cur < cellNum; ++cur)
				vX[cur] = vXn[cur];
		}

		return -1;//
	}



	void init_field()
	{
		
		for (int x = 0; x < nx;++x)
		for (int y = 0; y < ny; ++y)
		{
			const int cur = cellIndex(x, y);
			bool bound = (x == 0) || (x == nx - 1) || (y == 0) || (y == ny - 1);

			if (bound)
			{
				phi_cur[cur] = phi_nex[cur] = Tbound;
				bstatus[cur] = true;//boundary layer
			}
			else
			{
				phi_cur[cur] = phi_nex[cur] = 100.0f;
				bstatus[cur] = false;
			}
		}
	}
	
	void allo_mem()
	{
		phi_cur = (float*)malloc(sizeof(float)*cellNum);
		phi_nex = (float*)malloc(sizeof(float)*cellNum);
		vB = (float*)malloc(sizeof(float)*cellNum);
		vX = (float*)malloc(sizeof(float)*cellNum);
		vXn = (float*)malloc(sizeof(float)*cellNum);
		bstatus = (bool*)malloc(sizeof(bool)*cellNum);


		coeff_a = (float*)malloc(sizeof(float)*cellNum);
		coeff_b = (float*)malloc(sizeof(float)*cellNum);
		coeff_c = (float*)malloc(sizeof(float)*cellNum);
		coeff_d = (float*)malloc(sizeof(float)*cellNum);
		coeff_e = (float*)malloc(sizeof(float)*cellNum);

	}



	int cellIndex(int x, int y)
	{
		return x*ny + y;
	}


	void writeTec(int mx, int my,int nframe, float* var)
	{
		int i, j, z;
		FILE *fp;
		char filename[100];
		sprintf(filename, "field_%ld.plt", 10000+nframe);
		fp = fopen(filename, "w");
		fprintf(fp, "VARIABLES=\"X\",\"Y\",\"Var\"\n");
		fprintf(fp, "ZONE T=\"1\"\nI=%d,J=%d,F=POINT\n", mx, my);
		
		for (i = 0; i < mx; i++)
		for (j = 0; j < my; j++)
		{
			z = i*my + j;
			fprintf(fp, "%d %d %f\n", i, j, var[z]);
			fprintf(fp, "\n");
		}
		fclose(fp);
			
	}

private:

};



#endif