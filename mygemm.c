#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double R = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                R += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = R;
        }
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 2)
    {
        for (j = 0; j < n; j += 2)
        {
            register double C_0_0 = C[i * n + j];
            register double C_0_1 = C[i * n + (j + 1)];
            register double C_1_0 = C[(i + 1) * n + j];
            register double C_1_1 = C[(i + 1) * n + (j + 1)];

            for (k = 0; k < n; k += 2)
            {
                register double A_0_0 = A[i * n + k];
                register double A_0_1 = A[i * n + (k + 1)];               
                register double A_1_0 = A[(i + 1) * n + k];
                register double A_1_1 = A[(i + 1) * n + (k + 1)];
                register double B_0_0 = B[k * n + j];
                register double B_0_1 = B[k * n + (j + 1)];                
                register double B_1_0 = B[(k + 1) * n + j];
                register double B_1_1 = B[(k + 1) * n + (j + 1)];

                C_0_0 += A_0_0 * B_0_0 + A_0_1 * B_1_0;
                C_0_1 += A_0_0 * B_0_1 + A_0_1 * B_1_1;                
                C_1_0 += A_1_0 * B_0_0 + A_1_1 * B_1_0;
                C_1_1 += A_1_0 * B_0_1 + A_1_1 * B_1_1;
            }
    
            C[i * n + j] = C_0_0;
            C[i * n + (j + 1)] = C_0_1;
            C[(i + 1) * n + j] = C_1_0;            
            C[(i + 1) * n + (j + 1)] = C_1_1;
        }
    }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 3)
    {
        for (j = 0; j < n; j += 3)
        {
            register double C_0_0 = C[i * n + j];
            register double C_0_1 = C[i * n + (j + 1)];
            register double C_0_2 = C[i * n + (j + 2)];
            register double C_1_0 = C[(i + 1) * n + j];
            register double C_1_1 = C[(i + 1) * n + (j + 1)];
            register double C_1_2 = C[(i + 1) * n + (j + 2)];
            register double C_2_0 = C[(i + 2) * n + j];
            register double C_2_1 = C[(i + 2) * n + (j + 1)];
            register double C_2_2 = C[(i + 2) * n + (j + 2)];

            for (k = 0; k < n; k++)
            {
                register double A_0 = A[i * n + k];
                register double A_1 = A[(i + 1) * n + k];
                register double A_2 = A[(i + 2) * n + k];
                register double B_0 = B[k * n + j];
                register double B_1 = B[k * n + (j + 1)];
                register double B_2 = B[k * n + (j + 2)];

                C_0_0 += A_0 * B_0;
                C_0_1 += A_0 * B_1;
                C_0_2 += A_0 * B_2;
                C_1_0 += A_1 * B_0;
                C_1_1 += A_1 * B_1;
                C_1_2 += A_1 * B_2;
                C_2_0 += A_2 * B_0;
                C_2_1 += A_2 * B_1;
                C_2_2 += A_2 * B_2;
            }

            C[i * n + j] = C_0_0;
            C[i * n + (j + 1)] = C_0_1;
            C[i * n + (j + 2)] = C_0_2;
            C[(i + 1) * n + j] = C_1_0;
            C[(i + 1) * n + (j + 1)] = C_1_1;
            C[(i + 1) * n + (j + 2)] = C_1_2;    
            C[(i + 2) * n + j] = C_2_0;
            C[(i + 2) * n + (j + 1)] = C_2_1;
            C[(i + 2) * n + (j + 2)] = C_2_2;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{

}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jik(const double *A, const double *B, double *C, const int n) 
{

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kij(const double *A, const double *B, double *C, const int n) 
{

}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{

}


void ikj(const double *A, const double *B, double *C, const int n) 
{

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jki(const double *A, const double *B, double *C, const int n) 
{

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kji(const double *A, const double *B, double *C, const int n) 
{

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{

}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}
