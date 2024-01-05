/*  Projeto Final para HPC2
*   Autor: Luis Armando Quintanilla Villon
*   Data : 22/09/2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h> 
#include <unistd.h> 
#include "mpi.h"

/*
* IMPORTANTE: 
* - É necessário 2 ou mais processos para executar este programa.
* - É necessário que o número de processos seja múltiplo de Nx.
*/

// Altura das sub-regiões.
const int  Nx = 384;
const int  Ny = 384;

// Constantes físicas do problema.
const double k = 150.0;
const double h = 450.0; 
const double Tinf = 0.0;
const double H = 0.1;
const double L = 0.1;
const double Sp = -12.0;
const double Sc = 12500.0;
    
// Erro inicial e precisão. 
double erro = 1.0;
const double precisao = 0.00001; 

// Função que imprime a matriz final
void print_matrix(int x, int y, double *m)
{     
    int i,j;
    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            printf("%5.2f ", m[(y * i) + j]);
        }
        printf("\n");
    }
}

// Função que imprime a matriz auxiliar com processos ghosts
void mpi_print_matrix_aux(int x, int y, double *m)
{    
    int i,j;
    for (i = 0; i < x; i++)
    {
        if ((i == 0) || (i == x - 1)) printf("Ghost: ");
        else printf("       ");
        for (j = 0; j < y; j++)
        {
            printf("%5.2f ", m[(y * i) + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Função que imprime a matriz com ghost
void mpi_print_matrix(int x, int y, double *m, int myid, int nprocs, char name)
{   
    if (0 == myid)
    {
        printf("----------------- Matriz %c: ----------------\n", name);
        double buffer[x * y];
        mpi_print_matrix_aux(x, y, m);
        //Recebendo valores dos outros processos e imprimindo
        int id;
        for (id = 1; id < nprocs; id++)
        {
            MPI_Recv(&buffer[0], x * y, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpi_print_matrix_aux(x, y, buffer);
        }
    }
    else
    {
        MPI_Send(&m[0], x * y, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

// Função que armazena a matriz.
void mpi_save_matrix_aux(int x, int y, double *m)
{    
    int i,j;
    FILE *f1;
    f1 = fopen("matriz-384-mpi.csv","a");
    for (i = 1; i <= x - 2; i++)
    {
        for (j = 0; j < y; j++)
        {
            fprintf(f1,"%5.2f ", m[(y * i) + j]);
        }
        fprintf(f1,"\n");
    }
    fclose(f1);
}

// Função que armazena a matriz real sem ghost.
void mpi_save_matrix(int x, int y, double *m, int myid, int nprocs)
{   
    if (0 == myid)
    {
        double buffer[x * y];
        mpi_save_matrix_aux(x, y, m);
        int id;
        for (id = 1; id < nprocs; id++)
        {
            MPI_Recv(&buffer[0], x * y, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpi_save_matrix_aux(x, y, buffer);
        }
    }
    else
    {
        MPI_Send(&m[0], x * y, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

// Função que armazena o tempo de execução na seção paralela
void mpi_save_time(int nprocs, double t)
{
    FILE *f1;
    f1 = fopen("matriz-384.txt","a+");
    fprintf(f1,"%d -- %.6f ", nprocs, t);
    fprintf(f1,"\n");
    fclose(f1);
}

// Condições iniciais
void mpi_initial_conditions(int x, int y, double *m, int myid, int nprocs)
{
    int i,j;
    for(i = 0; i < x; i++)
    {
        for(j = 0; j < y; j++)
        {
            if ((i == 0) || (i == x - 1))
            {
                m[(y * i) + j] = -1;
            }
            else
            {
                m[(y * i) + j] = 75.0;
            }
        }
    }
}

// Comunicações na região A.
void mpi_comunications_A(int x, int y, double *m, int myid, int nprocs)
{
    /*
      'myid' envia a primeira fila real para o último 'ghost' de 'myid - 1'.
      'myid' recebe a ultima fila real de 'myid - 1' e o armazena no 
      primeiro 'ghost'.
      
      'myid' envia sua ultima fila real para o primeiro 'ghost' de 'myid + 1'.
      'myid' recebe a primeira fila real de 'myid + 1' e o armazena no 
      último 'ghost'.
    */
    MPI_Request request;
    if (0 != myid)
    {
        MPI_Isend(&m[y + Ny], Ny, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, &request);
        MPI_Recv(&m[0 + Ny], Ny, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
    if (myid != nprocs - 1)
    {
        MPI_Isend(&m[((x - 2) * y) + Ny], Ny, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Recv(&m[((x - 1) * y) + Ny], Ny, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

// Comunicações entre as regiões A e B.
void mpi_comunications_A_B(int x, int y, double *mA, double *mB, int myid, int nprocs)
{
    /*
      'myid' envia a primeira fila real para o último 'ghost' de 'myid - 1'.
      'myid' recebe a ultima fila real de 'myid - 1' e o armazena no 
      primeiro 'ghost'.
      
      'myid' envia sua ultima fila real para o primeiro 'ghost' de 'myid + 1'.
      'myid' recebe a primeira fila real de 'myid + 1' e o armazena no 
      último 'ghost'.
    */
    MPI_Request request;
    if (nprocs > 1)
    {
        if (0 == myid)
        {
            MPI_Isend(&mB[y + Ny], Ny, MPI_DOUBLE, nprocs - 1, 0, MPI_COMM_WORLD, &request);
            MPI_Recv(&mB[0 + Ny], Ny, MPI_DOUBLE, nprocs - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
        if (myid == nprocs - 1)
        {
            MPI_Isend(&mA[((x - 2) * y) + Ny], Ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
            MPI_Recv(&mA[((x - 1) * y) + Ny], Ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
    }
}

// Comunicações na região B.
void mpi_comunications_B(int x, int y, double *m, int myid, int nprocs)
{
    /*
      'myid' envia a primeira fila real para o último 'ghost' de 'myid - 1'.
      'myid' recebe a ultima fila real de 'myid - 1' e o armazena no 
      primeiro 'ghost'.
      
      'myid' envia sua ultima fila real para o primeiro 'ghost' de 'myid + 1'.
      'myid' recebe a primeira fila real de 'myid + 1' e o armazena no 
      último 'ghost'.
    */
    MPI_Request request;
    if (0 != myid)
    {
        MPI_Isend(&m[y], y, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, &request);
        MPI_Recv(&m[0], y, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
    if (myid != nprocs - 1)
    {
        MPI_Isend(&m[(x - 2) * y], y, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, &request);
        MPI_Recv(&m[(x - 1) * y], y, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
}

// Condições de contorno constantes.
void mpi_boundary_constant(int x, int y, double *mA, double *mB, int myid, int nprocs)
{ 
    int i,j;
    //---------------Condições de contorno com valores constantes fixos------------
    //--------------------------------------Contorno 1.
    if (nprocs > 1)
    {
        j = 0;
        if(0 == myid)
        {
            for(i = 2; i <= x - 2; i++)
            {
                mB[(y * i) + j] = 100.0;
            }
        }
        else
        {
            for(i = 1; i <= x - 2; i++)
            {
                mB[(y * i) + j] = 100.0;
            }
        }
    }
    //--------------------------------------Contorno 4.
    if (myid == 0)
    {
        i = 1;
        for(j = Ny + 1; j <= 2*Ny-2; j++)
        {
            mA[(y * i) + j] = 50.0;
        } 
    }
    //--------------------------------------Contorno 6.
    if (myid == nprocs - 1)
    {
        i = x - 2;
        for(j = 1; j <= 2*Ny-2; j++)
        {
            mB[(y * i) + j] = 50.0;
        } 
    }
}

//--------------------------------------Contorno variável 2.
void mpi_boundary_2(int x, int y, double *mB, double *mBnova, int myid, int nprocs)
{
    int i, j;
    if (0 == myid)
    {
        i = 1;
        for(j = 0; j <= Ny-1; j++)
        {
            mB[(y * i) + j] = mB[(y * (i + 1)) + j];
            mBnova[(y * i) + j] = mBnova[(y * (i + 1)) + j];
        } 
    }
}

//--------------------------------------Contorno variável 3.
void mpi_boundary_3(int x, int y, double *mA, double *mAnova, int myid, int nprocs)
{
    int i, j = Ny;
    for(i = 1; i <= x - 2; i++)
    {
        mA[(y * i) + j] = mA[(y * i) + j + 1];
        mAnova[(y * i) + j] = mAnova[(y * i) + j + 1];
    } 
}

//--------------------------------------Contorno variável 5.
void mpi_boundary_5(int x, int y, double *mA, double *mB, double *mAnova, double *mBnova, int myid, int nprocs, double c7, double c8)
{
    int i, j = 2*Ny-1;
    for(i = 1; i <= x - 2; i++)
    {
        mA[(y * i) + j] = c8 * mA[(y * i) + j - 1] + c7;
        mB[(y * i) + j] = c8 * mB[(y * i) + j - 1] + c7;
        mAnova[(y * i) + j] = c8 * mAnova[(y * i) + j - 1] + c7;
        mBnova[(y * i) + j] = c8 * mBnova[(y * i) + j - 1] + c7;
    } 
}

//-----------------------Cálculo dos elementos da matriz A e Anova.
void mpi_heat_A(int x, int y, double *mA, double *mAnova, int myid, int nprocs, double c3, double c4, double c5, double c6)
{
    erro = 0.0;
    int i, j;
    for(i = 1; i <= x - 2; i++)
	{
	    for(j = Ny + 1; j <= 2*Ny-2; j++)
	    {
            if (0 == myid && 1 == i) continue;
	        mAnova[(y * i) + j] = c4 + c3*mA[(y * i) + j - 1] + c3*mA[(y * i) + j + 1] + c5*mA[(y * (i - 1)) + j] + c6*mA[(y * (i + 1)) + j];
		}
	}
    for(i = 1; i <= x - 2; i++)
	{
	    for(j = Ny + 1; j <= 2*Ny-2; j++)
	    {
	        // Calculando a diferença entre a matriz nova e velha.
            if (0 == myid && 1 == i) continue;
	        erro = fmax(fabs(mAnova[(y * i) + j] - mA[(y * i) + j]), erro);
	        mA[(y * i) + j] = mAnova[(y * i) + j];
	    }
	}
}

//-----------------------Cálculo dos elementos da matriz B e Bnova.
void mpi_heat_B(int x, int y, double *mB, double *mBnova, int myid, int nprocs, double c3, double c4, double c5, double c6)
{
    int i, j;
    for(i = 1; i <= x - 2; i++)
	{
	    for(j = 1; j <= 2*Ny-2; j++)
	    {
            if ((0 == myid) && (1 == i) && (j <= Ny - 1)) continue;
            if ((nprocs - 1 == myid) && (x - 2 == i)) continue;
	        mBnova[(y * i) + j] = c4 + c3*mB[(y * i) + j - 1] + c3*mB[(y * i) + j + 1] + c5*mB[(y * (i - 1)) + j] + c6*mB[(y * (i + 1)) + j];
		}
	}
    for(i = 1; i <= x - 2; i++)
	{
	    for(j = 1; j <= 2*Ny-2; j++)
	    {
	        // Calculando a diferença entre a matriz nova e velha.
            if ((0 == myid) && (1 == i) && (j <= Ny - 1)) continue;
            if ((nprocs - 1 == myid) && (x - 2 == i)) continue;  
	        erro = fmax(fabs(mBnova[(y * i) + j] - mB[(y * i) + j]), erro);
	        mB[(y * i) + j] = mBnova[(y * i) + j];
	    }
	}
    double erro_B = erro;
    MPI_Allreduce(&erro_B, &erro, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
    // Variaveis auxiliares para o cálculo dos elementos da matriz.
    const double dx = L/((double)Nx-1.0);
    const double dy = H/((double)Ny-1.0);
    const int Nx_2 = 2*Nx;
    const int Ny_2 = 2*Ny;
    const double alpha = dx/dy;
    const double An = alpha*alpha;
    const double As = alpha*alpha;
    const double Ap = -(2+2*alpha*alpha)+Sp*dx*dx;
    const double Su = -Sc*dx*dx;
    const double c1 = h*dx*Tinf/k; 
    const double c2 = 1+h*dx/k;
    const double c3 = -(1/Ap);
    const double c4 = -(1/Ap)*(-Su);
    const double c5 = -(1/Ap)*An;
    const double c6 = -(1/Ap)*As;
    const double c7 = c1/c2;
    const double c8 = 1/c2;

    // Inicializando variáveis e parametros.
    // Variáveis do ambiente MPI
    int myid, nproc;
    double start_time, max_time, min_time, avg_time, local_time ;
    MPI_Status rstatus;
    
    MPI_Init(&argc, &argv); // Inicialização dos processos
    MPI_Comm_rank( MPI_COMM_WORLD, &myid ); // Itentificação de cada processo os processos mediante 'myid'
    MPI_Comm_size( MPI_COMM_WORLD, &nproc); // Retorno do número de processos mediante 'nproc' 
    
    // Criação das matrizes locales correspondentes as rigóes A e B.
    char name[] = {'A', 'B'}; // Rótulo das matrizes A e B. 
    char name_novo[] = {'C', 'D'}; // Rótulo das matrizes novas A(C) e B(D). 
    int i,j;
    int Nl = Nx/nproc + 2; // Altura da matriz local, ique inclui dois processos ghosts. 
    double Ml_A[Nl * Ny_2]; // Matriz local na região A
    double Ml_B[Nl * Ny_2]; // Matriz local na região B
    double Ml_Anova[Nl * Ny_2]; // Matriz local nova na região A 
    double Ml_Bnova[Nl * Ny_2]; // Matriz local nova na região B
    
    // Inizialização das condições inicias 
    mpi_initial_conditions(Nl, Ny_2, Ml_A, myid, nproc);
    mpi_initial_conditions(Nl, Ny_2, Ml_B, myid, nproc);
    mpi_initial_conditions(Nl, Ny_2, Ml_Anova, myid, nproc);
    mpi_initial_conditions(Nl, Ny_2, Ml_Bnova, myid, nproc);
    
    // Inizialização das condições de contorno constantes. 
    mpi_boundary_constant(Nl, Ny_2, Ml_A, Ml_B, myid, nproc);
    mpi_boundary_constant(Nl, Ny_2, Ml_Anova, Ml_Bnova, myid, nproc); 
    
    MPI_Barrier(MPI_COMM_WORLD);  // Sinchronização inicial dos processes
    
    // Inizialização da comunicação mediante mpi. 
    mpi_comunications_A(Nl, Ny_2, Ml_A, myid, nproc);
    mpi_comunications_A_B(Nl, Ny_2, Ml_A, Ml_B, myid, nproc);
    mpi_comunications_B(Nl, Ny_2, Ml_B, myid, nproc);
    mpi_comunications_A(Nl, Ny_2, Ml_Anova, myid, nproc);
    mpi_comunications_A_B(Nl, Ny_2, Ml_Anova, Ml_Bnova, myid, nproc);
    mpi_comunications_B(Nl, Ny_2, Ml_Bnova, myid, nproc);
    
    start_time = MPI_Wtime();  // Medição inicial do tempo.
    //--------------------------------------Loop principal para resolver o sistema--------------------.    
    while ( precisao < erro )
    {    
        // ---------- Atualização dos contornos variáveis.-----------
        //-------------------------------------- Contorno 2.	
		mpi_boundary_2(Nl, Ny_2, Ml_B, Ml_Bnova, myid, nproc);
		//-------------------------------------- Contorno 3.
        mpi_boundary_3(Nl, Ny_2, Ml_A, Ml_Anova, myid, nproc);
	    //-------------------------------------- Contorno 5.
        mpi_boundary_5(Nl, Ny_2, Ml_A, Ml_B, Ml_Anova, Ml_Bnova, myid, nproc, c7, c8);
        // Atualização da comunicação das duas matrizes mediante mpi. 
        mpi_comunications_A(Nl, Ny_2, Ml_Anova, myid, nproc);
        mpi_comunications_A_B(Nl, Ny_2, Ml_Anova, Ml_Bnova, myid, nproc);
        mpi_comunications_B(Nl, Ny_2, Ml_Bnova, myid, nproc);
        mpi_comunications_A(Nl, Ny_2, Ml_A, myid, nproc);
        mpi_comunications_A_B(Nl, Ny_2, Ml_A, Ml_B, myid, nproc);
        mpi_comunications_B(Nl, Ny_2, Ml_B, myid, nproc);
        MPI_Barrier(MPI_COMM_WORLD);   
        //----------------------------------Cálculo da região A.
        mpi_heat_A(Nl, Ny_2, Ml_A, Ml_Anova, myid, nproc, c3, c4, c5, c6);

		//----------------------------------Cálculo da região B.
        //mpi_print_matrix(Nl, Ny_2, Ml_B, myid, nproc, name[0]);
        mpi_heat_B(Nl, Ny_2, Ml_B, Ml_Bnova, myid, nproc, c3, c4, c5, c6);
    }
    //printf("error: %f, process: %d\n", erro, myid);

    local_time = MPI_Wtime() - start_time;  // Contagem final do tempo por processo
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0,MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

    /*
    //Impressão da martriz final.
    //mpi_print_matrix(Nl, Ny_2, Ml_Anova, myid, nproc, name[0]);
    //mpi_print_matrix(Nl, Ny_2, Ml_Bnova, myid, nproc, name[1]);
    */

    /* 
    //Armazenamento da martriz final.
    //mpi_save_matrix(Nl, Ny_2, Ml_Anova, myid, nproc);
    //mpi_save_matrix(Nl, Ny_2, Ml_Bnova, myid, nproc);
    */

    if (myid == 0)
    {
        double tempoh = max_time/3600;
        printf("Tempo em segundos: %f \n", max_time);
        printf("Tempo em horas: %f \n", tempoh);
        mpi_save_time(nproc, max_time);  // Armazanando o tempo
    }
    
    MPI_Finalize(); //Finalização do ambiente do MPI
    return 0;
}
