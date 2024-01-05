/*  Programa otimizado para HPC1
*   Autor: Luis Armando Quintanilla Villon
*   Data: 05/05/2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h> 
#include <unistd.h> 

#define  Nx 720
#define  Ny 720

//Inicializando as matrizes necessárias. 
static double Mvelha[1439][1439];
static double Mnova[1439][1439];

// Função que armazena a matriz final
void savematrix(size_t a, size_t b, double m[a][b])
{     
     int i,j;
     FILE *f1;
     f1=fopen("matriz-1440.csv","w");
     for (j=a-1;j>=0;j--)
     {
         for (i=0;i<b;i++)
         {
             fprintf(f1,"%5.2f ", m[i][j]);
             }
             fprintf(f1,"\n");
     }
     fclose(f1);
}

void contorno_constante(int l1, int l2, int l4, int l5, int maxx, int maxy)
{
    int i,j;
    for(i = 0; i<=maxx-1; i++)
    {
       for(j = 0; j<=maxy-1; j++)
       {
           Mvelha[i][j] = 75.0;
           Mnova[i][j] = 75.0;
       }
    }
    //--------------------------Condiciones de contorno com valores constantes fixos-------------------
    //--------------------------------------Contorno 1.
    for(j = 0; j<=l1; j++)
    {
       Mvelha[0][j] = 100.0;
       Mnova[0][j] = 100.0;
    }
    //--------------------------------------Contorno 4.
    for(i = l2; i<=l4; i++)
    {
       Mvelha[i][l5] = 50.0;
       Mnova[i][l5] = 50.0;
    }
    //--------------------------------------Contorno 6.
    for(i = 1; i<=l4; i++)
    {
       Mvelha[i][0] = 50.0;
       Mnova[i][0] = 50.0;
    }    
}

int main()
{ 
    //Inicializando variáveis e parametros.
    double tempo = 0.0;
    clock_t ini = clock();
    double  k = 150.0;
    double h =450.0; 
    double  Tinf = 0.0;
    double  H = 0.1;
    double  L = 0.1;
    double  Sp = -12.0;
    double  Sc = 12500.0;
    double dx = L/((double)Nx-1.0);
    double dy = H/((double)Ny-1.0);
       
    const int maxx = 2*Nx-1;
    const int maxy = 2*Ny-1;
    
    const double alpha = dx/dy;
    const double An = alpha*alpha;
    const double As = alpha*alpha;
    const double Ap = -(2+2*alpha*alpha)+Sp*dx*dx;
    const double Su = -Sc*dx*dx;
    
    //Variaveis auxiliares e liomites dos contornos:
    const int l1 = Ny-1;
    const int l2 = Nx-1;
    const int l3 = 2*Ny-3;
    const int l4 = 2*Nx-2;
    const int l5 = 2*Ny-2;
    const int l6 = Ny-2;
    const int l7 = 2*Nx-3;
    const int l8 = maxx-2;
    const int l9 = maxy-2;
    const int l10 = l8*l6+(Nx-2)*l1;
    
    const double c1 = h*dx*Tinf/k; 
    const double c2 = 1+h*dx/k;
    const double c3 = -(1/Ap);
    const double c4 = -(1/Ap)*(-Su);
    const double c5 = -(1/Ap)*An;
    const double c6 = -(1/Ap)*As;
    const double c7 = c1/c2;
    const double c8 = 1/c2;
    
    //Definindo o erro e a precisão. 
    double erro = 1.0;
    double precisao = 0.00001; 
    //int p = omp_get_max_threads();
    int i,j;
    contorno_constante(l1,l2,l4,l5,maxx,maxy);

    //---------------------------------------Loop principal para resolver o sistema--------------------.    
    while ( precisao < erro )
    {
                //---------------------------------------Contorno 2.	
		for(i = 1; i<=l2; i++)
		{
		    Mvelha[i][l1] = Mvelha[i][l6];
		    Mnova[i][l1] = Mnova[i][l6];
		}
		
		//---------------------------------------Contorno 3.
		for(j = Ny; j <=l3; j++)
		{
		    Mvelha[l2][j] = Mvelha[Nx][j];
		    Mnova[l2][j] = Mnova[Nx][j];
		}
		
		//----------------------------------------Contorno 5.
		for(j = 1; j<=l3; j++)
		{
		    Mvelha[l4][j] = c8*Mvelha[l7][j]+c7;
		    Mnova[l4][j] = c8*Mvelha[l7][j]+c7;
		}

		//----------------------------------Cálculo da primeira parte da matriz nova.
		erro = 0.0;
		for(i = 1; i<=l8; i++)
		{
		    for(j = 1; j<=l6; j++)
		    {
		         Mnova[i][j] = c4+c3*Mvelha[i-1][j]+c3*Mvelha[i+1][j]+c5*Mvelha[i][j+1]+c6*Mvelha[i][j-1];
		    }
		}
		
		for(i = 1; i<=l8; i++)
		{
		    for(j = 1; j<=l6; j++)
		    {
		         // Calculando a diferença entre a matriz nova e velha.
		         erro = fmax(fabs(Mnova[i][j]-Mvelha[i][j]),erro);
		         Mvelha[i][j] = Mnova[i][j];
		    }
		}
                 
		//----------------------------------Cálculo da segunda ṕarte da matriz nova.
		for(i = Nx; i<=l8; i++)
		{
		    for(j = l1; j<=l9; j++)
		    {
		         Mnova[i][j] = c4+c3*Mvelha[i-1][j]+c3*Mvelha[i+1][j]+c5*Mvelha[i][j+1]+c6*Mvelha[i][j-1];
		    }
		}
		
		for(i = Nx; i<=l8; i++)
		{
		    for(j = l1; j<=l9; j++)
		    {
		         // Calculando a diferença entre a matriz nova e velha.
		         erro = fmax(fabs(Mnova[i][j]-Mvelha[i][j]),erro);
		         Mvelha[i][j] = Mnova[i][j];
		    }
		}
     }
   //Contagem final do tempo.
   clock_t fim = clock();
   tempo += (double)(fim - ini) / CLOCKS_PER_SEC;
   double tempoh = tempo/3600;
   //Armazanando a matriz
   savematrix(maxx, maxy, Mnova);  
   printf("Tempo em segundos: %f \n", tempo);
   printf("Tempo em horas: %f \n", tempoh);
   return 0;
}
