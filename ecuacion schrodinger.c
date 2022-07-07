/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// ECUACUACIÓN DE SCHRODINGER                                                  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

//BIBLIOTECAS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"

//CONSTANTES DEL PROGRAMA
#define N 100
#define T 1000      //tiempo máximo
#define nciclos 10  //Esta variable va de 1 a N/4 como máximo
#define lambda 3.0 
#define PI 3.14159265359

//VARIABLES EXTERNAS
//Variables reales
double k0,s,norma;
double V[N]; 
//Variables complejas (N+1 para tener en cuenta las condiciones de contorno)
fcomplex phi[N];
fcomplex chi[N];
fcomplex alpha[N-1];
fcomplex beta[N-1];
fcomplex A0[N];
//Ficheros
FILE *f1,*f2,*f3;

//DECLARACIÓN DE FUNCIONES
int CalculoAlpha();
int CalculoBeta();
int CalculoChi();
int MostrarCondicionesIniciales();

//FUNCIÓN PRINCIPAL........................................................................................................

int main()
{
    /*
    ...........................................................................................................
    .  La idea es calcular primero las condiciones iniciales para t=0                                         .
    .  COSAS A CALCULAR:                                                                                      .       
    .    -Constantes k0 y s                                                                                   .
    .    -Potencial V                                                                                         .
    .    -Vector phi (para representarlo se usa el módulo cuadrado, que es un vector de números reales)       .
    .    -Vector de alphas (independiente del tiempo)                                                         .
    ...........................................................................................................
    */

   //Abrir los ficheros
   f1=fopen("Potencial.txt","w");
   f2=fopen("Schrodinger.txt","w");
   f3=fopen("Norma.txt","w");

   //Dar valores a los parámetros iniciales
   k0=2*PI*nciclos/N;
   s=0.25/pow(k0,2);

   //Fijar condiciones de contorno
   phi[0]=phi[N-1]=Complex(0.0,0.0);
   chi[0]=chi[N-1]=Complex(0.0,0.0);

   for(int i=0;i<N;i++)
   {
       //Cálculo del potencial (al cuadrado)
       if(i>=(2*N/5-1) && i<=(3*N/5-1)) V[i]=lambda*pow(k0,2);
       else V[i]=0;

       //Cálculo del vector A0
       A0[i]=Complex(-2-V[i],2/s);
   }

   for(int i=0;i<T;i++)
   {
       //Guardar el valor del potencial para cada tiempo en el fichero
       for(int j=0;j<N;j++) fprintf(f1,"%i\t&lf\n",j,V[j]);

       fprintf(f1,"\n\n");
   }

   CalculoAlpha();
   MostrarCondicionesIniciales();


   //BUCLE TEMPORAL: se calcula phi para cada unidad de tiempo, que es lo que al final se representa

   for(int i=0;i<T;i++)
   {
        CalculoBeta();
        CalculoChi();

        //Cálculo y guardado de la phi en el fichero
        for(int j=0;j<N;j++) phi[j]=Csub(chi[j],phi[j]);
        for(int j=0;j<N;j++) fprintf(f2,"%i\t%lf\n",j,pow(Cabs(phi[j],2)));

        fprintf(f3,"%i\t%lf\n",i,norma);
   }

   //Cerrar los ficheros
   fclose(f1);
   fclose(f2);
   fclose(f3);

   return 0;
}

//DEFINICIÓN DE FUNCIONES..................................................................................................

int CalculoAlpha() //Me calcula los valores del vector alpha
{
    //Inicializo la última componente del vector a 0
    alpha[N-2]=Complex(0.0,0.0);
    //Bucle para completar todos los valores de alpha
    for (int i=N-2;i>0;i--) alpha[i-1]=Cdiv(Complex(-1.0,0.0),Cadd(A0[i],alpha[i]));

    return 0;
}

int CalculoBeta()  //Me calcula los valores del vector beta
{
    //Inicializo la última componente del vector a 0
    beta[N-2]=Complex(0.0,0.0);
    //Bucle para completar todos los valores de alpha
    for (int i=N-2;i>0;i--) beta[i-1]=Cdiv(Csub(Cmul(Complex(0.0,4.0/s),phi[i]),beta[i]),Cadd(A0[i],alpha[i]));

    return 0;
}

int CalculoChi()   //Me calcula los valores del vector phi
{
    for(int i=1;i<N-1;i++) chi[i]=Cadd(Cmul(alpha[i-1],chi[i-1]),beta[i-1]);

    return 0;
}

int MostrarCondicionesIniciales()  //Me muestra por pantalla todas las condiciones iniciales
{
    //Vector de potencial
    printf("VECTOR DE POTENCIAL\n\n");
    for(int i=0,i<N;i++) printf("V[%i]=%lf\n",i,V[i]);

    printf("\n\n\n");

    //Vector alpha
    printf("VECTOR alpha\n\n");
    for(int i=0,i<N-1;i++) 
    {
        printf("Alpha[%i]=%lf+%lf i\n",i,alpha[i].r,alpha[i].i);
    }

    printf("\n\n\n");

    //Vector phi
    printf("VECTOR PHI\n\n");
    for(int i=0,i<N;i++) printf("Phi[%i]=%lf+%lf i\n",i,phi[i].r,phi[i].i);

    printf("\n\n\n");

    //Norma total, k0 y s
    printf("NORMA TOTAL: %lf \n\n",norma);
    printf("k0: %lf \n\n",k0);
    printf("s: %lf \n\n",s);

    printf("\n\n\n");

    return 0;
}