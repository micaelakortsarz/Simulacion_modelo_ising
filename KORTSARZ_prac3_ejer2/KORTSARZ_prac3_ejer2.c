#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define it_MC 10000
float ran2(long *idum);
//-------------------------------------------------------------------------------------------------------------
//Estructura
typedef struct
{
    int L;          //Tamaño del sistema
    double T;       //Temperatura
    double delta;    //delta/J
    double **s;        //Partículas de la red
    double Magnetizacion;
    double Energia;
    int it;             //Nro de iteraciones
    long semilla;
} Red_t;

double Calcula_delta_E(Red_t *r,double p,int i, int j);
//-------------------------------------------------------------------------------------------------------------
//Manejo de memoria de matrices
void liberaMatriz(double ** pmat, int nfil)
{
    int i;
    // primero debo liberar los elementos de pmat
    for( i=0; i<nfil; i++)
        free(pmat[i]);
    // y finalmente libero pmat
    free(pmat);
    return;
}

double ** allocaMatriz(int nfil, int ncol)
{
    int i;
    double  **pmat;
    // pmat apunta a un arreglo de nfil punteros
    pmat = (double **) malloc(nfil * sizeof(double *));
    if( pmat == NULL )
        return NULL;
    for( i=0; i<nfil; i++)
    {
        pmat[i]=(double *) malloc(ncol * sizeof(double ));
        if( pmat[i] == NULL )
        {
            liberaMatriz(pmat,i);
            return NULL;
        }
    }
    return pmat;
}
//----------------------------------------------------------------------------------------------------------------
//Funciones para calcular magnitudes iniciales e inicializar el sistema
int Cond_Con(int i,int L)           //Condiciones de contorno periodicas
{
    if(i==-1)
        return L-1;
    if(i==L)
        return 0;
    return i;
}

void Energia(Red_t *r)
{
    int i,j;
    double E=0;
    for(i=0; i<r->L; i++)
    {
        for(j=0; j<r->L; j++){
            E+=-0.5*r->s[i][j]*(r->s[Cond_Con(i+1,r->L)][j]+r->s[Cond_Con(i-1,r->L)][j]+r->s[i][Cond_Con(j+1,r->L)]+r->s[i][Cond_Con(j-1,r->L)]);
            E+=r->delta*r->s[i][j]*r->s[i][j];
        }
    }
    r->Energia=E;
}

void Magnetizacion(Red_t *r)
{
    int i,j;
    double M=0;
    for(i=0; i<r->L; i++)
    {
        for(j=0; j<r->L; j++)
            M+=r->s[i][j];
    }
    r->Magnetizacion=M;
}

void Inicializa_Sistema(Red_t *r, int L, double T,double delta)
{
    int i,j;
    r->it=0;
    r->L=L;
    r->T=T;
    r->delta=delta;
    r->s=allocaMatriz(L,L);
    r->semilla=56789;

    for(i=0; i<L; i++)
    {
        for(j=0; j<L; j++)
            r->s[i][j]=1;
    }
    Energia(r);
    Magnetizacion(r);
}

double Calcula_delta_E(Red_t *r,double p,int i, int j)
{
    double dE=0;
    dE=(-p + r->s[i][j])*(r->s[Cond_Con(i+1,r->L)][j]+r->s[Cond_Con(i-1,r->L)][j]+r->s[i][Cond_Con(j+1,r->L)]+r->s[i][Cond_Con(j-1,r->L)]);
    dE+=r->delta*(p*p - r->s[i][j]*r->s[i][j]);
    return dE;
}
//---------------------------------------------------------------------------------------------------------------
//Histogramas

void Arma_Histograma(Red_t *r, double *P,int nroit,int Nbin)
{
    int N=r->L*r->L,i;
    double dx=2.0/Nbin, E=(r->Energia)/N;
    for(i=0; i<Nbin+1; i++)
    {
        if(-1.1+(i*dx)<=E && (-1.1+(i+1)*dx)>=E)
            P[i]++;
    }
    return;
}

//-----------------------------------------------------------------------------------------------------------------
//Resolucion por monte carlo
double Elige_Propuesta_s(Red_t *r,int i, int j)
{
   double eta= ran2(&(r->semilla)),n;
   if(r->s[i][j]==0)
   {
      if(eta<=0.5)
        n=-1;
      else n=1;
      return n;
   }

   else
   {
      if(eta<=0.5)
        n=-1*r->s[i][j];
      else n=0;
      return n;
   }
}

void Asigna_Probabilidades(Red_t *r,double dE,double prop,int i, int j)
{

    if(dE<=0)
    {
        r->Magnetizacion+=prop-r->s[i][j];
        r->s[i][j]=prop;
        r->Energia+=dE;

        return;
    }

    double eta= ran2(&r->semilla);
    double p=exp(-dE/r->T);
    if(eta<p)
    {
        r->Magnetizacion+=prop-r->s[i][j];
        r->s[i][j]=prop;
        r->Energia+=dE;
    }
    return;
}

void Monte_Carlo(Red_t *r)
{
    int i,j;
    double dE,p;
    for(i=0; i<r->L; i++)
    {
        for(j=0; j<r->L; j++)
        {
            p=Elige_Propuesta_s(r,i,j);
            dE=Calcula_delta_E(r,p,i,j);
            Asigna_Probabilidades(r,dE,p,i,j);

        }
    }
    r->it++;
}

void Resuelve_Monte_Carlo(FILE *f,double T,int NroIt, Red_t *r,int imprime)
{
    srand(time(NULL));
    int i,j,N=r->L*r->L;
     r->T=T;
    double Etot=0,Etot2=0,Mtot=0,Mabs=0,Mtot2=0;
    double x,C,xx;
  int Nbin=5.0*(1.0+3.322*log(NroIt));

    double *P=(double *) malloc((NroIt+1) * sizeof(double ));



    double norm=(1.0/(N*(0.9*NroIt))),dx=2.0/Nbin;

    for(i=0; i<NroIt+1; i++)
        P[i]=0;

    do
    {
        if(imprime==1)
        fprintf(f,"%d %lf %lf\n",r->it,r->Energia/N,(r->Magnetizacion)/N);

         if(imprime==3)
            Arma_Histograma(r,P,NroIt,Nbin);

        Monte_Carlo(r);


         if(r->it>0.1*NroIt)
        {
            Etot+=r->Energia;
            Etot2+=r->Energia*r->Energia;
            Mtot+=r->Magnetizacion;
            Mtot2+=r->Magnetizacion*r->Magnetizacion;
            Mabs+=abs(r->Magnetizacion);
        }
    }
    while(r->it<NroIt);

    double Etotp=Etot*norm;
    double Etot2p=Etot2*norm;
    double  Mtotp=Mtot*norm;
    double  Mtot2p=Mtot2*norm;
    double Mabsp=Mabs*norm;
    x=(Mtot2p-(Mtotp*Mtotp)*N)/(r->T);      //Susceptibilidad por sitio
    xx=(Mtot2p-(Mabsp*Mabsp)*N)/(r->T);     //Suceptibilidad por sitio calculada con el valor medio del modulo de la magnetizacion
    C=(Etot2p-(Etotp*Etotp)*N)/(r->T*r->T);  //Calor especifico por sitio

    if(imprime==2)
        fprintf(f,"%lf %lf %lf %lf %lf\n",(r->T),Mabsp,xx,Etotp,C);

    if(imprime==3)
    { int n=0;
        for(i=0; i<Nbin+1; i++)
        {
            if(P[i]!=0)
                n++;

        }
        double dx2=2.0/n;
        for(i=0; i<Nbin+1; i++)
        {
            if(P[i]!=0)
            {
                fprintf(f,"%lf %lf %lf\n",-1.1+1.0*i*dx,P[i]/(NroIt*dx2),-r->T*log(P[i]/(NroIt*dx2)));
                P[i]=P[i]/(NroIt*dx2);
            }
        }
    }
    liberaMatriz(&P,1);

     r->it=0;
}

//--------------------------------------------------------------------------------------------------------------------------
//Ejercicios
void Ejercicio_2_a()
{
    FILE *f;
    Red_t r;

    f=fopen("2a T=0.4.txt","w");
    Inicializa_Sistema(&r,20,0.4,1.975);
    Resuelve_Monte_Carlo(f,0.4,5*it_MC,&r,1);
    fclose(f);
    liberaMatriz(r.s,r.L);

     f=fopen("2a T=0.585.txt","w");
    Inicializa_Sistema(&r,20,0.585,1.975);
    Resuelve_Monte_Carlo(f,0.585,5*it_MC,&r,1);
    fclose(f);
    liberaMatriz(r.s,r.L);


    f=fopen("2a T=1.txt","w");
    Inicializa_Sistema(&r,20,1.0,1.975);
    Resuelve_Monte_Carlo(f,1.0,5*it_MC,&r,1);
    fclose(f);

    liberaMatriz(r.s,r.L);

}

void Ejercicio_2_b()
{
    double T;
    Red_t r;
    FILE *f;

   /* f=fopen("2b L=20 delta=1.975.txt","w");
     Inicializa_Sistema(&r,20,0.1,1.975);
    for(T=10; T<100; T++)
        Resuelve_Monte_Carlo(f,0.01*T,100000,&r,2);
    for(T=100; T<300; T++)
        Resuelve_Monte_Carlo(f,0.01*T,it_MC,&r,2);
    fclose(f);
    liberaMatriz(r.s,r.L);*/

    f=fopen("2b L=100 delta=1.975.txt","w");
     Inicializa_Sistema(&r,100,0.1,1.975);
    for(T=10; T<53; T++)
        Resuelve_Monte_Carlo(f,0.01*T,10000,&r,2);
    for(T=53; T<60; T++)
        Resuelve_Monte_Carlo(f,0.01*T,1000000,&r,2);
    for(T=60; T<300; T++)
        Resuelve_Monte_Carlo(f,0.01*T,10000,&r,2);
    fclose(f);
    liberaMatriz(r.s,r.L);

  /*  f=fopen("2b L=20 delta=1.txt","w");
     Inicializa_Sistema(&r,20,0.1,1);
    for(T=10; T<100; T++)
        Resuelve_Monte_Carlo(f,0.01*T,10000,&r,2);
    for(T=100; T<300; T++)
        Resuelve_Monte_Carlo(f,0.01*T,100000,&r,2);
    fclose(f);
    liberaMatriz(r.s,r.L);*/

 /*  f=fopen("2b L=100 delta=1.txt","w");
     Inicializa_Sistema(&r,100,0.1,1);
    for(T=10; T<100; T++)
        Resuelve_Monte_Carlo(f,0.01*T,10000,&r,2);
    for(T=100; T<300; T++)
        Resuelve_Monte_Carlo(f,0.01*T,100000,&r,2);
    fclose(f);
    liberaMatriz(r.s,r.L);*/
}

void Ejercicio_2_c()
{
    double T;
    Red_t r;
    FILE *f;

    f=fopen("2c1 delta=1.975.txt","w");
     Inicializa_Sistema(&r,30,0.2,1.975);
    for(T=10; T<100; T++){
        if(T<=50)
        Resuelve_Monte_Carlo(f,0.02*T,50000,&r,2);
        else Resuelve_Monte_Carlo(f,(2- 0.02*T),50000,&r,2);
    }
    fclose(f);
    liberaMatriz(r.s,r.L);

   /*f=fopen("2c delta=1.95.txt","w");
     Inicializa_Sistema(&r,20,0.2,1.95);
    for(T=20; T<200; T++){
        if(T<=100)
        Resuelve_Monte_Carlo(f,0.01*T,100000,&r,2);
        else Resuelve_Monte_Carlo(f,(2- 0.01*T),100000,&r,2);
    }
    fclose(f);
    liberaMatriz(r.s,r.L);
*/
 /*  f=fopen("2c delta=1.txt","w");
     Inicializa_Sistema(&r,20,0.2,1.0);
    for(T=100; T<320; T++){
        if(T<=200)
        Resuelve_Monte_Carlo(f,0.01*T,100000,&r,2);
        else Resuelve_Monte_Carlo(f,(4.0- 0.01*T),100000,&r,2);
    }
    fclose(f);
    liberaMatriz(r.s,r.L);*/
}

void Ejercicio_2_d()
{
    double T;
    Red_t r;
    FILE *f;

   f=fopen("2d T=0.581 delta=1.975.txt","w");
     Inicializa_Sistema(&r,20,0.581,1.975);
        Resuelve_Monte_Carlo(f,0.581,10000000,&r,3);
    fclose(f);
    liberaMatriz(r.s,r.L);


}

void Evolucion_Sistema_T()
{

   int T,j,i;
   int NroIt=10000;
    FILE *f;
    Red_t r;
    char str[80];
    Inicializa_Sistema(&r,100,0.1,1.975);

    for(T=10; T<250; T++)
    {
       sprintf(str, "EvT delta =1.975 %d.txt",T);
        f=fopen(str,"w");
        r.T=0.01*T;
       do
    {
        Monte_Carlo(&r);
    }
    while(r.it<NroIt);
    for(i=0;i<r.L;i++)
    {
        for(j=0;j<r.L;j++)
            fprintf(f,"%d %d %lf\n",i,j,r.s[i][j]);
    }
        fclose(f);
    r.it=0;
   }
 liberaMatriz(r.s,r.L);


 Inicializa_Sistema(&r,100,0.1,1.0);

    for(T=10; T<250; T++)
    {
       sprintf(str, "EvT delta =1 %d.txt",T);
        f=fopen(str,"w");
        r.T=0.01*T;
       do
    {
        Monte_Carlo(&r);
    }
    while(r.it<NroIt);
    for(i=0;i<r.L;i++)
    {
        for(j=0;j<r.L;j++)
            fprintf(f,"%d %d %lf\n",i,j,r.s[i][j]);
    }
        fclose(f);
    r.it=0;
   }
 liberaMatriz(r.s,r.L);
}

void Histogramas_en_T()
{
    FILE *f;
    Red_t r;
    int T;
    char str[80];
    int NroIt=5000000;
    int i,n=0;

    double *P =(double *) malloc((NroIt+1) * sizeof(double ));

    int Nbin=100.0*(1.0+3.322*log(NroIt));
    double dx=2.2/Nbin,M;

Inicializa_Sistema(&r,20,0.3,1.975);

    for(T=300; T<750; T++)
    {

        sprintf(str, "H 1.975 %d.txt",T);
        f=fopen(str,"w");
        r.T=0.001*T;
    for(i=0; i<NroIt+1; i++)
        P[i]=0;

       do
    {
        Monte_Carlo(&r);
        M=r.Energia/(r.L*r.L);

    for(i=0; i<Nbin+2; i++)
    {
        if((-1.0+1.0*i*dx)<=M && (-1.0+(1.0*i+1)*dx)>=M)
            P[i]++;
    }

    } while(r.it<NroIt);

        for(i=0; i<Nbin+1; i++)
        {
            if(P[i]!=0)
                n++;

        }
        double dx2=2.2/n;
        for(i=0; i<Nbin+1; i++)
        {
            if(P[i]!=0)
            {
                fprintf(f,"%lf %lf %lf\n",-1.0+1.0*i*dx,P[i]/(NroIt*dx2),-r.T*log(P[i]/(NroIt*dx2)));
                P[i]=P[i]/(NroIt*dx2);
            }

    }
            printf("holis\n");
    r.it=0;
        fclose(f);
    }

     liberaMatriz(r.s,r.L);
       free(P);

       Inicializa_Sistema(&r,20,1.4,1.0);

    for(T=1400; T<1850; T++)
    {

        sprintf(str, "H 1 %d.txt",T);
        f=fopen(str,"w");
        r.T=0.001*T;
    for(i=0; i<NroIt+1; i++)
        P[i]=0;

       do
    {
        Monte_Carlo(&r);
        M=r.Energia/(r.L*r.L);

    for(i=0; i<Nbin+2; i++)
    {
        if((-1.0+1.0*i*dx)<=M && (-1.0+(1.0*i+1)*dx)>=M)
            P[i]++;
    }

    } while(r.it<NroIt);

        for(i=0; i<Nbin+1; i++)
        {
            if(P[i]!=0)
                n++;

        }
        double dx2=2.2/n;
        for(i=0; i<Nbin+1; i++)
        {
            if(P[i]!=0)
            {
                fprintf(f,"%lf %lf %lf\n",-1.0+1.0*i*dx,P[i]/(NroIt*dx2),-r.T*log(P[i]/(NroIt*dx2)));
                P[i]=P[i]/(NroIt*dx2);
            }

    }
            printf("holis\n");
    r.it=0;
        fclose(f);
    }

     liberaMatriz(r.s,r.L);
       free(P);
}





//--------------------------------------------------------------------------------------------------------------------------
int main()
{
  //  Ejercicio_2_a();
//Ejercicio_2_b();
//Ejercicio_2_c();
//Ejercicio_2_d();
//Evolucion_Sistema_T();
Histogramas_en_T();
    return 0;
}
