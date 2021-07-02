#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define it_MC 10000

float ran2(long *idum);
//-----------------------------------------------------------------------------------------------

typedef struct
{
    int L;          //Tamaño del sistema
    double T;       //Temperatura
    double **s;        //Partículas de la red
    double Magnetizacion;
    double Energia;
    int it;             //Nro de iteraciones
} Red_t;
//-------------------------------------------------------------------------------------------------
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

int Cond_Con(int i,int L)
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
        for(j=0; j<r->L; j++)
            E+=-1.0*r->s[i][j]*(r->s[Cond_Con(i+1,r->L)][j]+r->s[Cond_Con(i-1,r->L)][j]+r->s[i][Cond_Con(j+1,r->L)]+r->s[i][Cond_Con(j-1,r->L)]);
    }
    r->Energia=E*0.5;
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

void Inicializa_Sistema(Red_t *r, int L, double T)
{
    int i,j;
    r->it=0;
    r->L=L;
    r->T=T;
    r->s=allocaMatriz(L,L);

    for(i=0; i<L; i++)
    {
        for(j=0; j<L; j++)
            r->s[i][j]=1;
    }
    Energia(r);
    Magnetizacion(r);
}

double Calcula_delta_E(Red_t *r, int i, int j)
{
    double dE=0;
    dE=2.0*r->s[i][j]*(r->s[Cond_Con(i+1,r->L)][j]+r->s[Cond_Con(i-1,r->L)][j]+r->s[i][Cond_Con(j+1,r->L)]+r->s[i][Cond_Con(j-1,r->L)]);
    return dE;
}

void Arma_Histograma(Red_t *r, double P[],int nroit,int Nbin)
{
    int N=r->L*r->L,i;
    double dx=2.0/Nbin, M=(r->Magnetizacion)/N;

    for(i=0; i<Nbin+2; i++)
    {
        if((-1.0+1.0*i*dx)<=M && (-1.0+(1.0*i+1)*dx)>=M)
            P[i]++;
    }
    return;
}

double Halla_Ms(int nrobin,double P[])
{
    int i;
    double dx=2.0/nrobin,max=-100;
    double mi;
    for(i=nrobin*0.5; i<nrobin+1; i++)
    {
        if(max<P[i])
        {
            max=P[i];
            mi=1.0*i;
        }
    }
    return -1.0+mi*dx;
}

double Calcula_F(double M,int nrobin,double P[],double T)
{
    int i;
    double dx=2.0/nrobin;
    for(i=0; i<nrobin+1; i++)
    {
        if((-1.0+1.0*i*dx)<=M && (-1.0+(1.0*i+1)*dx)>M)
            return  -T*log(P[i]);
    }
    return 0;
}


void Asigna_Probabilidades(Red_t *r,double dE,int i, int j)
{
    long int ran=rand();
    double eta= ran2(&ran);
    double p=exp(-dE/r->T);
    if(dE<=0)
    {
        r->s[i][j]=-1.0*(r->s[i][j]);
        r->Energia+=dE;
        r->Magnetizacion+=2*r->s[i][j];
        return;
    }
    if(eta<p)
    {
        r->s[i][j]=-1.0*(r->s[i][j]);
        r->Energia+=dE;
        r->Magnetizacion+=2*r->s[i][j];
    }
    return;
}

void Monte_Carlo(Red_t *r)
{
    int i,j;
    double dE;
    for(i=0; i<r->L; i++)
    {
        for(j=0; j<r->L; j++)
        {
            dE=Calcula_delta_E(r,i,j);
            Asigna_Probabilidades(r,dE,i,j);

        }
    }
    r->it++;
}

void Resuelve_Monte_Carlo(FILE *f,int L, double T,int NroIt,int imprime,Red_t *r,double *P)
{
    srand(time(NULL));

    double Etot=0,Etot2=0,Mtot=0,Mabs=0,Mtot2=0;
    double x,C,xx,f0,fms,ms;
    int i,n,N;
    N=L*L;
    int Nbin=100*(1.0+3.322*log(NroIt));


    double norm=(1.0/(N*(NroIt-100))),dx=2.0/Nbin;

    for(i=0; i<NroIt+1; i++)
        P[i]=0;

    r->T=T;

    do
    {
        if(imprime==1)
            fprintf(f,"%d %lf %lf\n",r->it,r->Energia/N,(r->Magnetizacion)/N);


        if(imprime==4)
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
        fprintf(f,"%lf %lf %lf\n",r->T,Mabsp,xx);



    if(imprime==4)
    {
        int n=0;
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
                fprintf(f,"%lf %lf %lf\n",-1+1.0*i*dx,P[i]/(NroIt*dx2),-r->T*log(P[i]/(NroIt*dx2)));
                P[i]=P[i]/(NroIt*dx2);
            }
        }
       /* ms=Halla_Ms(Nbin,P);
        f0=Calcula_F(0,Nbin,P,r->T);
        fms=Calcula_F(ms,Nbin,P,r->T);
        printf("deltaF= %lf\n",f0-fms);*/
    }


    r->it=0;
}

void Ejercicio_1_a()
{
    //FILE *f;
    /*f=fopen("1a T=2.5 L=10.txt","w");
    Resuelve_Monte_Carlo(f,10,2.5,it_MC,1,0);
    fclose(f);*/

    /* f=fopen("1a T=2.1 L=10.txt","w");
     Resuelve_Monte_Carlo(f,10,2.1,it_MC,1);
     fclose(f);*/

    /*  f=fopen("1a T=2.5 L=30.txt","w");
      Resuelve_Monte_Carlo(f,30,2.5,it_MC,1,0);
      fclose(f);

      f=fopen("1a T=2.1 L=30.txt","w");
      Resuelve_Monte_Carlo(f,30,2.1,it_MC,1,0);
      fclose(f);*/
}

void Ejercicio_1_b()
{
    /*   double T;
       FILE *f;
       f=fopen("1b3 L=10.txt","w");
       for(T=10; T<200; T++)
           Resuelve_Monte_Carlo(f,10,0.01*T,it_MC,2);
       for(T=200; T<250; T++)
           Resuelve_Monte_Carlo(f,10,0.01*T,10*it_MC,2);
       for(T=250; T<500; T++)
           Resuelve_Monte_Carlo(f,10,0.01*T,it_MC,2);
       fclose(f);

       f=fopen("1b3 L=20.txt","w");
       for(T=10; T<200; T++)
           Resuelve_Monte_Carlo(f,20,0.01*T,it_MC,2);
       for(T=200; T<250; T++)
           Resuelve_Monte_Carlo(f,20,0.01*T,10*it_MC,2);
       for(T=250; T<500; T++)
           Resuelve_Monte_Carlo(f,20,0.01*T,it_MC,2);
       fclose(f);

       f=fopen("1b3 L=30.txt","w");
       for(T=10; T<200; T++)
           Resuelve_Monte_Carlo(f,30,0.01*T,it_MC,2);
       for(T=200; T<250; T++)
           Resuelve_Monte_Carlo(f,30,0.01*T,10*it_MC,2);
       for(T=250; T<500; T++)
           Resuelve_Monte_Carlo(f,30,0.01*T,it_MC,2);
       fclose(f);

       f=fopen("1b3 L=50.txt","w");
       for(T=10; T<200; T++)
           Resuelve_Monte_Carlo(f,50,0.01*T,it_MC,2);
       for(T=200; T<250; T++)
           Resuelve_Monte_Carlo(f,50,0.01*T,10*it_MC,2);
       for(T=250; T<500; T++)
           Resuelve_Monte_Carlo(f,50,0.01*T,it_MC,2);
       fclose(f);

       f=fopen("1b3 L=100.txt","w");
       for(T=10; T<200; T++)
           Resuelve_Monte_Carlo(f,100,0.01*T,it_MC,2);
       for(T=200; T<250; T++)
           Resuelve_Monte_Carlo(f,100,0.01*T,10*it_MC,2);
       for(T=250; T<500; T++)
           Resuelve_Monte_Carlo(f,100,0.01*T,it_MC,2);
       fclose(f);*/
}

void Ejercicio_1_c()
{
    double T;
    FILE *f;
    Red_t r;
    f=fopen("1c L=100.txt","w");

    Inicializa_Sistema(&r,100,2.2);

    for(T=2200; T<2255; T+=2)
//        Resuelve_Monte_Carlo(f,100,0.001*T,100000,2,&r);
    liberaMatriz(r.s,r.L);
    fclose(f);
}

void Ejercicio_1_d()
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

Inicializa_Sistema(&r,10,1.7);

    for(T=170; T<350; T++)
    {

        sprintf(str, "H %d.txt",T);
        f=fopen(str,"w");
        r.T=0.01*T;
    for(i=0; i<NroIt+1; i++)
        P[i]=0;

       do
    {
        Monte_Carlo(&r);
        M=r.Magnetizacion/(r.L*r.L);

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

void Evolucion_Sistema_T()
{

   int T,j,i;
   int NroIt=10000;
    FILE *f;
    Red_t r;
    char str[80];
    Inicializa_Sistema(&r,100,1.0);

    for(T=100; T<350; T++)
    {
       sprintf(str, "EvT %d.txt",T);
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

void Evolucion_Sistema_iT()
{
int NroIt=10000,i,j;
    FILE *f;
    Red_t r;
    char str[80];
     Inicializa_Sistema(&r,100,1.0);
       do
    {
           if(r.it%50==0){
                sprintf(str, "EviT T=1 %d.txt",r.it);
        f=fopen(str,"w");
        for(i=0;i<r.L;i++)
    {
        for(j=0;j<r.L;j++)
            fprintf(f,"%d %d %lf\n",i,j,r.s[i][j]);
    }}
        Monte_Carlo(&r);
        if(r.it%50==0)
        fclose(f);
    }
    while(r.it<NroIt);
    r.it=0;
 liberaMatriz(r.s,r.L);

  Inicializa_Sistema(&r,100,2.2692);
       do
    {
           if(r.it%50==0){
                sprintf(str, "EviT T=3.5 %d.txt",r.it);
        f=fopen(str,"w");
        for(i=0;i<r.L;i++)
    {
        for(j=0;j<r.L;j++)
            fprintf(f,"%d %d %lf\n",i,j,r.s[i][j]);
    }}
        Monte_Carlo(&r);
        if(r.it%50==0)
        fclose(f);
    }
    while(r.it<NroIt);
    r.it=0;
 liberaMatriz(r.s,r.L);

  Inicializa_Sistema(&r,100,3.5);
      do
    {
           if(r.it%50==0){
                sprintf(str, "EviT T=Tc %d.txt",r.it);
        f=fopen(str,"w");
        for(i=0;i<r.L;i++)
    {
        for(j=0;j<r.L;j++)
            fprintf(f,"%d %d %lf\n",i,j,r.s[i][j]);
    }}
        Monte_Carlo(&r);
        if(r.it%50==0)
        fclose(f);
    }
    while(r.it<NroIt);
    r.it=0;
 liberaMatriz(r.s,r.L);

}
//-------------------------------------------------------------------------------------------------
int main()
{
    // Ejercicio_1_a();
    // Ejercicio_1_b();
//   Ejercicio_1_c();
    Ejercicio_1_d();
 //Evolucion_Sistema_T();
    return 0;
}
