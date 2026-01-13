#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <math.h> 
#ifdef _OPENMP
#include <omp.h>
#endif 
#define PARA_PI_Forward_Numerical_Simulation 3.14159265 
#define PARA_ORDER_Forward_Numerical_Simulation 4
#define Nx_Forward_Numerical_Simulation 300
#define Ny_Forward_Numerical_Simulation 300
#define Nz_Forward_Numerical_Simulation 300
#define dx_Forward_Numerical_Simulation 5
#define dy_Forward_Numerical_Simulation 5
#define dz_Forward_Numerical_Simulation 5
#define NTIME_Forward_Numerical_Simulation 1000
#define DETATIME_Forward_Numerical_Simulation 0.0005
#define FRE_Forward_Numerical_Simulation 30          
#define NAB_Forward_Numerical_Simulation 100              
#define NPML_Forward_Numerical_Simulation 50      
#define NPMLRPARA_Forward_Numerical_Simulation 1000000.0   
#define VVMAXPARA_Forward_Numerical_Simulation 3000
#define NQPARA_Forward_Numerical_Simulation 3
#define PARALY_Forward_Numerical_Simulation  300
#define PARAPAOX_Forward_Numerical_Simulation 150
#define PARAPAOY_Forward_Numerical_Simulation 150
#define PARAPAOZ_Forward_Numerical_Simulation 150
#define PARAPAORES1X_Forward_Numerical_Simulation 0
#define PARAPAORES1Y_Forward_Numerical_Simulation 0
float*** space3d(int nx, int ny, int nz);
float** space2d(int nr, int nc);
void wfile(char filename[], float*** data, int nx, int ny, int nz);
void wfile3d(char filename[], float*** data, int nx, int ny, int nz);
void ColPivot(float* c, int n, float x[]);
void main()
{
	printf("¡¾Forward¡¿");
	int FORWARD_SIZE_PARA_nx1;
	int FORWARD_SIZE_PARA_nx2;
	int FORWARD_SIZE_PARA_ny1;
	int FORWARD_SIZE_PARA_ny2;
	int FORWARD_SIZE_PARA_nz1;
	int FORWARD_SIZE_PARA_nz2;
	int FORWARD_SIZE_PARA_jx1;
	int FORWARD_SIZE_PARA_jx2;
	int FORWARD_SIZE_PARA_jy1;
	int FORWARD_SIZE_PARA_jy2;
	int FORWARD_SIZE_PARA_jz1;
	int FORWARD_SIZE_PARA_jz2;
	int FORWARD_SIZE_PARA_NNx;
	int FORWARD_SIZE_PARA_NNy;
	int FORWARD_SIZE_PARA_NNz;
	int FORWARD_SIZE_PARA_bx1;
	int FORWARD_SIZE_PARA_bx2;
	int FORWARD_SIZE_PARA_by1;
	int FORWARD_SIZE_PARA_by2;
	int FORWARD_SIZE_PARAix;
	int FORWARD_SIZE_PARAiy;
	int FORWARD_SIZE_PARAiz; 
	int FORWARD_TIMESIZEit; 
	int FORWARD_TIMESIZEk; 
	int FORWARD_PARA_ORDER;
	float PARAtime;
	float PARAsum1;
	float PARAsum2;
	float PARAsum3;
	int PARAipx;
	int PARAipy;
	int PARASX;
	int PARASY;
	int PARASZ;
	int PARATIME;
	FORWARD_SIZE_PARA_NNx = Nx_Forward_Numerical_Simulation + 2 * NPML_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_NNy = Ny_Forward_Numerical_Simulation + 2 * NPML_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_NNz = Nz_Forward_Numerical_Simulation + 2 * NPML_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_nx1 = PARA_ORDER_Forward_Numerical_Simulation; 
	FORWARD_SIZE_PARA_nx2 = FORWARD_SIZE_PARA_NNx - PARA_ORDER_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_ny1 = PARA_ORDER_Forward_Numerical_Simulation; 
	FORWARD_SIZE_PARA_ny2 = FORWARD_SIZE_PARA_NNy - PARA_ORDER_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_nz1 = PARA_ORDER_Forward_Numerical_Simulation; 
	FORWARD_SIZE_PARA_nz2 = FORWARD_SIZE_PARA_NNz - PARA_ORDER_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_jx1 = NPML_Forward_Numerical_Simulation; 
	FORWARD_SIZE_PARA_jx2 = FORWARD_SIZE_PARA_NNx - NPML_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_jy1 = NPML_Forward_Numerical_Simulation; 
	FORWARD_SIZE_PARA_jy2 = FORWARD_SIZE_PARA_NNy - NPML_Forward_Numerical_Simulation;
	FORWARD_SIZE_PARA_jz1 = NPML_Forward_Numerical_Simulation; 
	FORWARD_SIZE_PARA_jz2 = FORWARD_SIZE_PARA_NNz - NPML_Forward_Numerical_Simulation;
	float*** FORWARD_PARA_V;
	float*** FORWARD_PARA_FAI;
	float***FORWARD_PARA_PHO;
	float*** FORWARD_PARA_VElX;
	float*** FORWARD_PARA_VElY;
	float *** FORWARD_PARA_VEZ;
	float*** FORWARD_PARA_VELOCITY;
	float*** FORWARD_PARA_ESP;
	float*** FORWARD_PARA_QUA;
	float*** FORWARD_PARA_DEN;
	float*** FORWARD_PARA_DDEN;
	float*** FORWARD_PARA_K;
	float*** FORWARD_PARA_QQUA;
	float*** FORWARD_PARA_VV1;
	float*** FORWARD_PARA_VV2;
	float*** FORWARD_PARA_FAI1;
	float*** FORWARD_PARA_FAI2;
	float*** PARAMX;
	float*** PARASIA;
	float*** PARAesp;
	float*** PARAESP;
	float*** PARAdeta;
	float*** PARADETA;
	float* PARAC, * PARAWAVELET;
	float* PARATAO, * PARAD;
	float*** DATAPaoS;
	float*** PARA_theta;
	float*** PARA_Theta;
	float*** PARAfi;
	float*** PARAFI;
	float*** PARG1;
	float*** PARG2;
	float*** PARG3;
	float*** FORWARD_DATA;
	FORWARD_PARA_V = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation); 
	FORWARD_PARA_DEN = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation); 
	FORWARD_PARA_ESP = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation);
	PARAesp = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation); 
	PARAdeta = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation); 
	PARA_theta = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation);
	PARAfi = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation);
	FORWARD_PARA_PHO = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	FORWARD_PARA_VElX = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_VElY = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_VEZ = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	FORWARD_PARA_VELOCITY = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	FORWARD_PARA_QQUA = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_K = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_DDEN = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_QUA = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	PARA_Theta = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	PARAFI = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_FAI = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	PARAESP = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	PARADETA = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	PARAMX = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	PARASIA = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_VV1 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_VV2 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_FAI1 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_PARA_FAI2 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	PARG1 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	PARG2 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz); 
	PARG3 = space3d(FORWARD_SIZE_PARA_NNx, FORWARD_SIZE_PARA_NNy, FORWARD_SIZE_PARA_NNz);
	FORWARD_DATA = space3d(Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation);;
	DATAPaoS = space3d(PARALY_Forward_Numerical_Simulation, Nx_Forward_Numerical_Simulation, NTIME_Forward_Numerical_Simulation);
	PARAC = (float*)calloc(PARA_ORDER_Forward_Numerical_Simulation, sizeof(float));
	PARAWAVELET = (float*)calloc(NTIME_Forward_Numerical_Simulation, sizeof(float));
	PARATAO = (float*)calloc(NQPARA_Forward_Numerical_Simulation, sizeof(float));
	PARAD = (float*)calloc(NQPARA_Forward_Numerical_Simulation, sizeof(float));
   #ifdef _OPENMP
   #pragma omp parallel for
   #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_V[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_DEN[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_ESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				PARAesp[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				PARAdeta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				PARA_theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
			{
				PARAfi[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_FAI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_VElY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_VV1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_VV2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_FAI1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_FAI2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_VELOCITY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_DDEN[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				FORWARD_PARA_QUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARADETA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				PARAMX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARASIA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				PARG1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARG2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
				PARG3[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
			}
		}
	}
	PARAD[0] = 1.74278563;
	PARAD[1] = 1.41221145;
	PARAD[2] = 1.72183357;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < PARALY_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_TIMESIZEk = 0; FORWARD_TIMESIZEk < NTIME_Forward_Numerical_Simulation; FORWARD_TIMESIZEk++)
			{
				DATAPaoS[FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAix][FORWARD_TIMESIZEk] = 0.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_TIMESIZEit = 0; FORWARD_TIMESIZEit < NTIME_Forward_Numerical_Simulation; FORWARD_TIMESIZEit++)
	{
		PARAWAVELET[FORWARD_TIMESIZEit] = 0.0;
	}
	PARATAO[0] = 0.00159180; 
	PARATAO[1] = 0.01450505;
	PARATAO[2] = 0.1402500;
	float S = 0;
	for (int i = 0; i < NQPARA_Forward_Numerical_Simulation; i++)
	{
		S = S + PARAD[i];
	}
	float A[PARA_ORDER_Forward_Numerical_Simulation][PARA_ORDER_Forward_Numerical_Simulation + 1];
	for (int i = 1; i <= PARA_ORDER_Forward_Numerical_Simulation; i++)
	{
		for (int j = 1; j <= PARA_ORDER_Forward_Numerical_Simulation; j++)
		{
			A[i - 1][j - 1] = pow(2.0 * j - 1, 2.0 * i - 1);
		}
	}
	for (int i = 1; i < PARA_ORDER_Forward_Numerical_Simulation; i++)
	{
		A[i][PARA_ORDER_Forward_Numerical_Simulation] = 0;
	}
	A[0][PARA_ORDER_Forward_Numerical_Simulation] = 1.0;
	ColPivot(A[0], PARA_ORDER_Forward_Numerical_Simulation, PARAC);
	for (FORWARD_TIMESIZEit = 0; FORWARD_TIMESIZEit < NTIME_Forward_Numerical_Simulation; FORWARD_TIMESIZEit++)
	{
		PARAtime = (FORWARD_TIMESIZEit - NAB_Forward_Numerical_Simulation) * DETATIME_Forward_Numerical_Simulation;
		PARAWAVELET[FORWARD_TIMESIZEit] = (1 - 2 * pow((PARA_PI_Forward_Numerical_Simulation * FRE_Forward_Numerical_Simulation * PARAtime), 2)) * exp(-pow((PARA_PI_Forward_Numerical_Simulation * FRE_Forward_Numerical_Simulation * PARAtime), 2));
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_VELOCITY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 4000.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
 	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_DDEN[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 1.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_QUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 50.0;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.25;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				PARADETA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.05;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 30 * PARA_PI_Forward_Numerical_Simulation / 180;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 30 * PARA_PI_Forward_Numerical_Simulation / 180;
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_K[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = FORWARD_PARA_DDEN[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_VELOCITY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_VELOCITY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz];
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_VV1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = (2 - DETATIME_Forward_Numerical_Simulation * FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) / (2 + DETATIME_Forward_Numerical_Simulation * FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]);
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_VV2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 2 * DETATIME_Forward_Numerical_Simulation / (FORWARD_PARA_DDEN[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * (FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * DETATIME_Forward_Numerical_Simulation + 2));
			}
		}
	}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_FAI1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = (4 - DETATIME_Forward_Numerical_Simulation * (FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] + FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz])) / (4 + DETATIME_Forward_Numerical_Simulation * (FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] + FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]));
			}
		}
	}
	FILE* fpdata;
	fpdata = fopen("data.dat", "rb");
	for (FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
	{
		for (FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
		{
			fread(&FORWARD_DATA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][0], Nz_Forward_Numerical_Simulation, sizeof(float), fpdata);
		}
	}
	fclose(fpdata);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
	{
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
		{
			for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
			{
				FORWARD_PARA_FAI2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 4 * DETATIME_Forward_Numerical_Simulation / (4 + DETATIME_Forward_Numerical_Simulation * (FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] + FORWARD_PARA_QQUA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]));
			}
		}
	}
	FILE* file;
	char filename[30]; 
			PARASX = PARAPAOX_Forward_Numerical_Simulation + NPML_Forward_Numerical_Simulation;
			PARASY = PARAPAOY_Forward_Numerical_Simulation + NPML_Forward_Numerical_Simulation;
			PARASZ = PARAPAOZ_Forward_Numerical_Simulation + NPML_Forward_Numerical_Simulation;
			FORWARD_SIZE_PARA_bx1 = PARAPAORES1X_Forward_Numerical_Simulation + NPML_Forward_Numerical_Simulation;
			FORWARD_SIZE_PARA_bx2 = FORWARD_SIZE_PARA_bx1 + Nx_Forward_Numerical_Simulation;
			FORWARD_SIZE_PARA_by1 = PARAPAORES1Y_Forward_Numerical_Simulation + NPML_Forward_Numerical_Simulation;
			FORWARD_SIZE_PARA_by2 = FORWARD_SIZE_PARA_by1 + PARALY_Forward_Numerical_Simulation;
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
			for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
			{
				for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
				{
					for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
					{
						FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0; 
					}
				}
			}
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
			for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
			{
				for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
				{
					for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
					{
						FORWARD_PARA_FAI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
					}
				}
			}
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
			for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
			{
				for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
				{
					for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
					{
						FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
					}
				}
			}
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
			for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
			{
				for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
				{
					for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
					{
						FORWARD_PARA_VElY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
					}
				}
			}
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
			for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_NNx; FORWARD_SIZE_PARAix++)
			{
				for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_NNy; FORWARD_SIZE_PARAiy++)
				{
					for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_NNz; FORWARD_SIZE_PARAiz++)
					{
						FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0.0;
					}
				}
			}
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
			for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < Nx_Forward_Numerical_Simulation; FORWARD_SIZE_PARAix++)
			{
				for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < Ny_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiy++)
				{
					for (int FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < Nz_Forward_Numerical_Simulation; FORWARD_SIZE_PARAiz++)
					{
						FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = FORWARD_DATA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz];
					}
				}
			}
			for (FORWARD_TIMESIZEk = 1; FORWARD_TIMESIZEk < NTIME_Forward_Numerical_Simulation - 1; FORWARD_TIMESIZEk++)
			{
				FORWARD_PARA_PHO[PARASX][PARASY][PARASZ] = FORWARD_PARA_PHO[PARASX][PARASY][PARASZ] + PARAWAVELET[FORWARD_TIMESIZEk];
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
				for (int FORWARD_SIZE_PARAix = FORWARD_SIZE_PARA_nx1; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_nx2; FORWARD_SIZE_PARAix++)
				{
					for (int FORWARD_SIZE_PARAiy = FORWARD_SIZE_PARA_ny1; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_ny2; FORWARD_SIZE_PARAiy++)
					{
						for (int FORWARD_SIZE_PARAiz = FORWARD_SIZE_PARA_nz1; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_nz2; FORWARD_SIZE_PARAiz++)
						{
							float Variable1 = 0.0;
							float Variable2 = 0.0;
							float Variable3= 0.0;
							float Variable4 = 0.0;
							float Variable5 = 0.0;
							float Variable6 = 0.0;
							float Variable7 = 0.0;
							float Variable8 = 0.0;
							float Variable9 = 0.0;
							float Variable10 = 0.0;
							float Variable11 = 0.0;
							for (int FORWARD_PARA_ORDER = 0; FORWARD_PARA_ORDER < PARA_ORDER_Forward_Numerical_Simulation; FORWARD_PARA_ORDER++)
							{
								Variable1 = Variable1 + PARAC[FORWARD_PARA_ORDER] * (FORWARD_PARA_PHO[FORWARD_SIZE_PARAix + FORWARD_PARA_ORDER + 1][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] - FORWARD_PARA_PHO[FORWARD_SIZE_PARAix - FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) / dx_Forward_Numerical_Simulation;
							}
							for (int FORWARD_PARA_ORDER = 0; FORWARD_PARA_ORDER < PARA_ORDER_Forward_Numerical_Simulation; FORWARD_PARA_ORDER++)
							{
								Variable2 = Variable2 + PARAC[FORWARD_PARA_ORDER] * (FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy + FORWARD_PARA_ORDER + 1][FORWARD_SIZE_PARAiz] - FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy - FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiz]) / dy_Forward_Numerical_Simulation;
							}
							for (int FORWARD_PARA_ORDER = 0; FORWARD_PARA_ORDER < PARA_ORDER_Forward_Numerical_Simulation; FORWARD_PARA_ORDER++)
							{
								Variable3= Variable3+ PARAC[FORWARD_PARA_ORDER] * (FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz + FORWARD_PARA_ORDER + 1] - FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz - FORWARD_PARA_ORDER]) / dz_Forward_Numerical_Simulation;
							}
							Variable4 = pow(PARG1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz], 2);
							Variable5 = pow(PARG2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz], 2);
							Variable6 = pow(PARG3[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz], 2);
							Variable10 = -2 * (PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] - PARADETA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (Variable4 + Variable5) * Variable6 * (Variable4 + Variable5 + Variable6);
							Variable11 = ((1 + 2 * PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (Variable4+ Variable5) + Variable6) * ((1 + 2 * PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz])+ 2 * (1 + PARADETA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (Variable4 + Variable5) * Variable6);
							if (Variable11 != 0)
							{
								PARASIA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = Variable10 / Variable11;
							}
							else
							{
								PARASIA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = 0;
							}
						}
					}
				}
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
				for (int FORWARD_SIZE_PARAix = FORWARD_SIZE_PARA_nx1; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_nx2; FORWARD_SIZE_PARAix++)
				{
					for (int FORWARD_SIZE_PARAiy = FORWARD_SIZE_PARA_ny1; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_ny2; FORWARD_SIZE_PARAiy++)
					{
						for (int FORWARD_SIZE_PARAiz = FORWARD_SIZE_PARA_nz1; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_nz2; FORWARD_SIZE_PARAiz++)
						{
							FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = FORWARD_PARA_VV1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] + FORWARD_PARA_VV2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * PARG3[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * (1 + PARASIA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]);
						}
					}
				}
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
				for (int FORWARD_SIZE_PARAix = FORWARD_SIZE_PARA_nx1; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_nx2; FORWARD_SIZE_PARAix++)
				{
					for (int FORWARD_SIZE_PARAiy = FORWARD_SIZE_PARA_ny1; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_ny2; FORWARD_SIZE_PARAiy++)
					{
						for (int FORWARD_SIZE_PARAiz = FORWARD_SIZE_PARA_nz1; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_nz2; FORWARD_SIZE_PARAiz++)
						{
							FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = FORWARD_PARA_VV1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] + FORWARD_PARA_VV2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * PARG1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * (1 + 2 * PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (1 + PARASIA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]);
						}
					}
				}
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
				for (int FORWARD_SIZE_PARAix = FORWARD_SIZE_PARA_nx1; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_nx2; FORWARD_SIZE_PARAix++)
				{
					for (int FORWARD_SIZE_PARAiy = FORWARD_SIZE_PARA_ny1; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_ny2; FORWARD_SIZE_PARAiy++)
					{
						for (int FORWARD_SIZE_PARAiz = FORWARD_SIZE_PARA_nz1; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_nz2; FORWARD_SIZE_PARAiz++)
						{
							FORWARD_PARA_VElY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = FORWARD_PARA_VV1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_VElY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] + FORWARD_PARA_VV2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * PARG2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * (1 + 2 * PARAESP[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (1 + PARASIA[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]);
						}
					}
				}
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
				for (int FORWARD_SIZE_PARAix = FORWARD_SIZE_PARA_nx1; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_nx2; FORWARD_SIZE_PARAix++)
				{
					for (int FORWARD_SIZE_PARAiy = FORWARD_SIZE_PARA_ny1; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_ny2; FORWARD_SIZE_PARAiy++)
					{
						for (int FORWARD_SIZE_PARAiz = FORWARD_SIZE_PARA_nz1; FORWARD_SIZE_PARAiz < FORWARD_SIZE_PARA_nz2; FORWARD_SIZE_PARAiz++)
						{
							float Variable12 = 0.0;
							float Variable13 = 0.0;
							float Variable14 = 0.0;
							float Variable15 = 0.0;
							for (int FORWARD_PARA_ORDER = 0; FORWARD_PARA_ORDER < PARA_ORDER_Forward_Numerical_Simulation; FORWARD_PARA_ORDER++)
							{
								Variable12 = Variable12 + PARAC[FORWARD_PARA_ORDER] * (cos(PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * cos(PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VElX[FORWARD_SIZE_PARAix + FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] - FORWARD_PARA_VElX[FORWARD_SIZE_PARAix - FORWARD_PARA_ORDER - 1][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) + cos(PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * sin(PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy + FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiz] - FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy - FORWARD_PARA_ORDER - 1][FORWARD_SIZE_PARAiz]) - sin(PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz + FORWARD_PARA_ORDER] - FORWARD_PARA_VElX[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz - FORWARD_PARA_ORDER - 1])) / dx_Forward_Numerical_Simulation;
							}
							for (int FORWARD_PARA_ORDER = 0; FORWARD_PARA_ORDER < PARA_ORDER_Forward_Numerical_Simulation; FORWARD_PARA_ORDER++)
							{
								Variable13 = Variable13 + PARAC[FORWARD_PARA_ORDER] * (-sin(PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VElY[FORWARD_SIZE_PARAix + FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] - FORWARD_PARA_VElY[FORWARD_SIZE_PARAix - FORWARD_PARA_ORDER - 1][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) + cos(PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VElY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy + FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiz] - FORWARD_PARA_VElY[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy - FORWARD_PARA_ORDER - 1][FORWARD_SIZE_PARAiz])) / dy_Forward_Numerical_Simulation;
							}
							for (int FORWARD_PARA_ORDER = 0; FORWARD_PARA_ORDER < PARA_ORDER_Forward_Numerical_Simulation; FORWARD_PARA_ORDER++)
							{
								Variable14 = Variable14 + PARAC[FORWARD_PARA_ORDER] * (sin(PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * cos(PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix + FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] - FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix - FORWARD_PARA_ORDER - 1][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) + sin(PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * sin(PARAFI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy + FORWARD_PARA_ORDER][FORWARD_SIZE_PARAiz] - FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy - FORWARD_PARA_ORDER - 1][FORWARD_SIZE_PARAiz]) + cos(PARA_Theta[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz]) * (FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz + FORWARD_PARA_ORDER] - FORWARD_PARA_VEZ[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz - FORWARD_PARA_ORDER - 1])) / dz_Forward_Numerical_Simulation;
							}
							for (int in = 0; in < NQPARA_Forward_Numerical_Simulation; in++)
							{
								FORWARD_PARA_FAI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] = FORWARD_PARA_FAI1[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_FAI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] - FORWARD_PARA_FAI2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * FORWARD_PARA_FAI[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] / PARATAO[in] - FORWARD_PARA_FAI2[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz] * PARAD[in] * (Variable12 + Variable13 + Variable14) / PARATAO[in];
							}
						}
					}
				}
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
				for (int FORWARD_SIZE_PARAix = FORWARD_SIZE_PARA_bx1; FORWARD_SIZE_PARAix < FORWARD_SIZE_PARA_bx2; FORWARD_SIZE_PARAix++)
				{
					for (int FORWARD_SIZE_PARAiy = FORWARD_SIZE_PARA_by1; FORWARD_SIZE_PARAiy < FORWARD_SIZE_PARA_by2; FORWARD_SIZE_PARAiy++)
					{
						DATAPaoS[FORWARD_SIZE_PARAiy - FORWARD_SIZE_PARA_by1][FORWARD_SIZE_PARAix - FORWARD_SIZE_PARA_bx1][FORWARD_TIMESIZEk - 1] = FORWARD_PARA_PHO[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][PARASZ];
					}
				}
			}
	printf("\n\n¡¾Output_File¡¿===>¡¾Success¡¿");
	char pname01[] = "FORWARD_PARA_PHO.dat";
	wfile3d(pname01, FORWARD_PARA_PHO, Nx_Forward_Numerical_Simulation, Ny_Forward_Numerical_Simulation, Nz_Forward_Numerical_Simulation);
	printf("SUCCESS");
}
float** space2d(int nr, int nc)
{
	float** a;
	int i;
	a = (float**)calloc(nr, sizeof(float*));
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int i = 0; i < nr; i++)
		a[i] = (float*)calloc(nc, sizeof(float));
	return a;
}
float*** space3d(int nx, int ny, int nz)
{
	float*** a;
	int FORWARD_SIZE_PARAix, FORWARD_SIZE_PARAiy;
	a = (float***)calloc(nx, sizeof(float**));
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	for (int FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < nx; FORWARD_SIZE_PARAix++)
	{
		a[FORWARD_SIZE_PARAix] = (float**)calloc(ny, sizeof(float*));
		for (int FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < ny; FORWARD_SIZE_PARAiy++)
		{
			a[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy] = (float*)calloc(nz, sizeof(float));
		}
	}
	return a;
}
void wfile(char filename[], float*** data, int nx, int ny, int nz)
{
	int FORWARD_SIZE_PARAix, FORWARD_SIZE_PARAiy, FORWARD_SIZE_PARAiz;
	FILE* fp = fopen(filename, "wb");
	for (FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < ny; FORWARD_SIZE_PARAiy++)
		for (FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < nx; FORWARD_SIZE_PARAix++)
			fwrite(&data[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][0], nz, sizeof(float), fp);
	fclose(fp);
}
void wfile3d(char filename[], float*** data, int nx, int ny, int nz)
{
	int FORWARD_SIZE_PARAix, FORWARD_SIZE_PARAiy, FORWARD_SIZE_PARAiz;
	int ix1, iy1, iz1;
	FILE* fp = fopen(filename, "wb");
	for (FORWARD_SIZE_PARAiy = 0; FORWARD_SIZE_PARAiy < ny ; FORWARD_SIZE_PARAiy++)
	{
		for (FORWARD_SIZE_PARAix = 0; FORWARD_SIZE_PARAix < nx ; FORWARD_SIZE_PARAix++)
		{
			for (FORWARD_SIZE_PARAiz = 0; FORWARD_SIZE_PARAiz < nz ; FORWARD_SIZE_PARAiz++)
			{
				fwrite(&data[FORWARD_SIZE_PARAix][FORWARD_SIZE_PARAiy][FORWARD_SIZE_PARAiz], 1, sizeof(float), fp);
			}
		}
	}
	fclose(fp);
}
void ColPivot(float* c, int n, float x[])
{
	int i, j, PARAtime, FORWARD_TIMESIZEk;
	float p;
	for (i = 0; i <= n - 2; i++)
	{
		FORWARD_TIMESIZEk = i;
		for (j = i + 1; j <= n - 1; j++)
			if (fabs(*(c + j * (n + 1) + i)) > (fabs(*(c + FORWARD_TIMESIZEk * (n + 1) + i)))) FORWARD_TIMESIZEk = j;
		if (FORWARD_TIMESIZEk != i)
			for (j = i; j <= n; j++)
			{
				p = *(c + i * (n + 1) + j);
				*(c + i * (n + 1) + j) = *(c + FORWARD_TIMESIZEk * (n + 1) + j);
				*(c + FORWARD_TIMESIZEk * (n + 1) + j) = p;
			}
		for (j = i + 1; j <= n - 1; j++)
		{
			p = (*(c + j * (n + 1) + i)) / (*(c + i * (n + 1) + i));
			for (PARAtime = i; PARAtime <= n; PARAtime++) *(c + j * (n + 1) + PARAtime) -= p * (*(c + i * (n + 1) + PARAtime));
		}
	}
	for (i = n - 1; i >= 0; i--)
	{
		for (j = n - 1; j >= i + 1; j--)
			(*(c + i * (n + 1) + n)) -= x[j] * (*(c + i * (n + 1) + j));
		x[i] = *(c + i * (n + 1) + n) / (*(c + i * (n + 1) + i));
	}
}