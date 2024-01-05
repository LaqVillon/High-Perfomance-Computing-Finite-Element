/*  Programa sequencial base para HPC1
*   Autor: Luis Armando Quintanilla Villon
*   Data: 05/05/2021
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h> 
#include <unistd.h> 


// Função que calcula o valor máximo de uma matriz. 
double valormax(size_t a, size_t b, double m[a][b])
{
    double max = m[0][0]; 
    for(int i = 0; i<a; i++)
        {
            for( int j = 0; j<b; j++)
            {
                 if (m[i][j] > max)
			 max = m[i][j];
            }
        }
    return max;
}

// Função que armazena o valor de uma matriz em um novo arquivo.
/*
void savematrix(size_t a, size_t b, double m[a][b])
{
     
     int i,j;
     FILE *f1;
     f1=fopen("LuisArmando_matrix.csv","w");
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
*/


int main(int argc, char** argv)
{
    // Inicio da contagem de tempo.
    double tempo = 0.0;
    clock_t ini = clock();
    
    //Inicializando variáveis e parametros.
    double  k = 150.0;
    double h =450.0; 
    double  Tinf = 0.0;
    double  H = 0.1;
    double  L = 0.1;
    double  Sp = -12.0;
    double  Sc = 12500.0;
    
    //Dominío de solução
    int  Nx = 720;
    int  Ny = 720; 
    double dx = L/((double)Nx-1.0);
    double dy = H/((double)Ny-1.0);
       
    const int maxx = 2*Nx-1;
    const int maxy = 2*Ny-1;
    
    double alpha = dx/dy;
    double An = alpha*alpha;
    double As = alpha*alpha;
    double Ap = -(2+2*alpha*alpha)+Sp*dx*dx;
    double Su = -Sc*dx*dx;

    //Inicializando as matrizes necessárias. 
    static double Mvelha[1439][1439];
    static double Mnova[1439][1439];
    static double desv[1439][1439];
   
    //Introduzindo elementos em cada matriz.
    for(int i = 0; i<=maxx-1; i++)
        {
            for( int j = 0; j<=maxy-1; j++)
            {
                 Mvelha[i][j] = 75.0;
                 Mnova[i][j] = 75.0;
		 desv[i][j] = 0.0;
            }
        }
    
    //Definindo o erro e a precisão. 
    double erro = 1.0;
    double precisao = 0.00001;
    
    //Loop principal para resolver o sistema.
    while ( precisao < erro )
    {
	//Contorno 1.
	for(int j = 0; j<=Ny-1; j++)
        {
            Mvelha[0][j] = 100.0;
            Mnova[0][j] = 100.0;
        }
        
        //Contorno 2.
	for(int i = 1; i<=Nx-1; i++)
        {
            Mvelha[i][Ny-1] = Mvelha[i][Ny-2];
            Mnova[i][Ny-1] = Mnova[i][Ny-2];
        }
        
        //Contorno 3.
	for(int j = Ny; j <=2*Ny-3; j++)
        {
            Mvelha[Nx-1][j] = Mvelha[Nx][j];
            Mnova[Nx-1][j] = Mnova[Nx][j];
        }
        
        //Contorno 4.
	for(int i = Nx-1; i<=2*Nx-2; i++)
        {
            Mvelha[i][2*Ny-2] = 50.0;
            Mnova[i][2*Ny-2] = 50.0;
        }
        
        //Contorno 5.
	for(int j = 1; j<=2*Ny-3; j++)
        {
            Mvelha[2*Nx-2][j] = (Mvelha[2*Nx-3][j]+h*dx*Tinf/k)/(1+h*dx/k);
            Mnova[2*Nx-2][j] = (Mvelha[2*Nx-3][j]+h*dx*Tinf/k)/(1+h*dx/k);
        }
	
	//Contorno 6.
	for(int i = 1; i<=2*Nx-2; i++)
        {
            Mvelha[i][0] = 50.0;
            Mnova[i][0] = 50.0;
        }
        
        //Cálculo da primeira parte da matriz nova.
	for(int i = 1; i<=maxx-2; i++)
        {
            for( int j = 1; j<=Ny-2; j++)
            {
                 Mnova[i][j] = -(1/Ap)*(-Su+Mvelha[i-1][j]+Mvelha[i+1][j]+An*Mvelha[i][j+1]+As*Mvelha[i][j-1]);
            }
        }
        
        for(int i = 1; i<=maxx-2; i++)
        {
            for( int j = 1; j<=Ny-2; j++)
            {
                 // Calculando a diferença entre a matriz nova e velha.
                 double aux = Mnova[i][j]-Mvelha[i][j];
                 desv[i][j] = fabs(aux);
                 Mvelha[i][j] = Mnova[i][j];
            }
        }

        //Cálculo da segunda ṕarte da matriz nova.
	for(int i = Nx; i<=maxx-2; i++)
        {
            for( int j = Ny-1; j<=maxy-2; j++)
            {
                 Mnova[i][j] = -(1/Ap)*(-Su+Mvelha[i-1][j]+Mvelha[i+1][j]+An*Mvelha[i][j+1]+As*Mvelha[i][j-1]);
            }
        }
        
        for(int i = Nx; i<=maxx-2; i++)
        {
            for( int j = Ny-1; j<=maxy-2; j++)
            {
                 // Calculando a diferença entre a matriz nova e velha.
                 double aux = Mnova[i][j]-Mvelha[i][j];
                 desv[i][j] = fabs(aux);
                 Mvelha[i][j] = Mnova[i][j];
            }
        } 
        //atualização do erro
        erro = valormax(maxx, maxy, desv);
   }
   //Armazanando a matriz
   /*savematrix(maxx, maxy, Mnova);*/
   //Contagem final do tempo.
   clock_t fim = clock();
   tempo += (double)(fim - ini) / CLOCKS_PER_SEC;
   double tempoh = tempo/3600;
   printf("Tempo em segundos: %f \n", tempo);
   printf("Tempo em horas: %f \n", tempoh);
   return 0;
}
