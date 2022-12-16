//Program to calculate the performance parameters of the blade: 

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define SIZE 6
#define PI 3.14159
#define R2D (180/PI)

double calculate_sigma_prime(double B,double c,double r);
double calculate_beta1(double lemda);
double calculate_a(double beta,double sigma_prime,double CL,double Q,double CD);
double calculate_a_prime1(double a);
double calculate_i(double gamma,double beta);
double calculate_beta(double lemda,double a_prime,double a);
double calculate_a_prime(double sigma_prime,double CL,double lemda,double beta,double a,double Q);
double calculate_c(double r,double beta, double B, double lemda_r);
double calculate_error(double a1, double a2);
double Calculate_CP(double lemda_r[],double a[],double a_prime[], int size1,double CL[],double CD[],double beta[]);
double calculate_tip_loss(double beta,double r,double R,double B);
double NACA02CL(double ALPHA);
double NACA23018CL(double alpha);
double NACA25112CL(double alpha);
double NACA25112CD(double alpha);
double NACA0012CD(double alpha);
double NACA23018CD(double alpha);
double NACA0015CL(double alpha);
double NACA0015CD(double alpha);
double NACA64210CL(double alpha);
double NACA64210CD(double alpha);
double NACA23024CL(double alpha);
double NACA23024CD(double alpha);
double NACA63210CL(double alpha);
double NACA63210CD(double alpha);
double NACA2415CL(double alpha);
double NACA2415CD(double alpha);


int main()
{
    /* List of all input */
    double r[SIZE]={0.2,1,2,3,4,5};
    double gamma[SIZE]={61,74.3,84.9,89.1,91.3,92.6};
    double c[SIZE]={.7,.71,.44,.30,.23,.19};
    double lemda=8;
    /* End of input */
    /* Start of variable declaration */
    double P[SIZE]={0};
    double Cp[SIZE]={0};
    double eta[SIZE]={0};
    double R[SIZE]={0};
    double V[SIZE]={0};
    double B[SIZE]={0};
    double CD[SIZE]={0};
    double CL[SIZE]={0};
    double a[SIZE]={0};
    double a_new[SIZE]={0};
    double a_prime[SIZE]={0};
    double a_prime_new[SIZE]={0};
    double beta[SIZE]={0};
    double sigma_prime[SIZE]={0};
    double temp=0;
    double lemda_r[SIZE]={0};
    double I[SIZE]={0};
    double beta_new[SIZE]={0};
    double e_a[SIZE]={0};
    double e_a_prime[SIZE]={0};
    double Q[SIZE]={0};
    /* End of variable declaration */
    /* Start of additional variable declaration */
    int i;
    int j=0;
    double CP=0;
    FILE *fp1;
    FILE *fp;
    fp=fopen("NACA2415.csv","w");
    fp1=fopen("NACA2415Result.csv","w");
    fprintf(fp1,"j,r,c,gamma,a,a_prime,i,beta\n");
    fprintf(fp,"j,a[i],a_prime[i],beta[i],I[i],CL[i],e_a[i],e_a_prime[i]\n");
    /* End of additional variable declaration */


    /** Starting of main part **/

    printf("i+1\tr[i]\t\tc[i]\t\tgamma[i]\ta[i]\t\ta_prime[i]\tI[i]\t\tbeta[i]\n");
    for(i=0;i<SIZE;i++)
    {
    B[i]=3;
    lemda_r[i]=(r[i]*lemda)/(r[SIZE-1]);
    sigma_prime[i]=calculate_sigma_prime(B[i],c[i],r[i]);

    for(j=0;j<50;j++)
    {
        if(j!=0)
    {
       a_new[i]=a[i];
       a_prime_new[i]=a_prime[i];
    }
    if(j==0)
        beta[i]=calculate_beta1(lemda_r[i]);
    else
        beta[i]=calculate_beta(lemda_r[i],a_prime[i],a[i]);
    I[i]=calculate_i(gamma[i],beta[i]);

    CL[i]=NACA2415CL(I[i]);
    CD[i]=NACA2415CD(I[i]);

    Q[i]=calculate_tip_loss(beta[i],r[i],r[SIZE-1],B[i]);

    a[i]=calculate_a(beta[i],sigma_prime[i],CL[i],Q[i],CD[i]);


    if(j==0)
        a_prime[i]=calculate_a_prime1(a[i]);
    else
        a_prime[i]=calculate_a_prime(sigma_prime[i],CL[i],lemda_r[i],beta[i],a[i],Q[i]);
    if(j>=15 && a_prime>1)
    {
            a_prime[i]==a_prime_new[i];
            break;
    }
    if(j !=0 )
    {
        e_a[i]=calculate_error(a_new[i],a[i]);
        e_a_prime[i]=calculate_error(a_prime_new[i],a_prime[i]);
    }
    if(i==SIZE-1)
    {
        if(e_a[i]<=5 && e_a_prime[i]<=5 && j!=0)
        {
            break;
        }
    }
    if(e_a[i]<=0.0001 && e_a_prime[i]<=0.0001&& j!=0)
    {
        break;
    }
    fprintf(fp,"%d, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",j,a[i],a_prime[i],beta[i],I[i],CL[i],e_a[i],e_a_prime[i]);
    }
    c[i]=calculate_c(r[i],beta[i],B[i],lemda_r[i]);
    printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i+1,r[i],c[i],gamma[i],a[i],a_prime[i],I[i],beta[i],Q[i]);
    fprintf(fp1,"%d, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",i+1,r[i],c[i],gamma[i],a[i],a_prime[i],I[i],beta[i]);

    }
    fclose(fp);
    fclose(fp1);
    CP=Calculate_CP(lemda_r,a,a_prime,SIZE,CL,CD,beta);



}

double calculate_sigma_prime(double B,double c,double r)
{
    /* Program to calculate local solidity (sigma prime) */
    double sigma_prime=(B*c)/(2*PI*r);
    return sigma_prime;
}

double calculate_beta1(double lemda)
{
    double temp1 = atan(1/lemda)*R2D;
    double temp=(temp1*2)/3;
    double beta = 90 - temp;
    if(beta<0)
            {
                beta=-beta;
            }
    return beta;
}

double calculate_a(double beta,double sigma_prime,double CL,double Q,double CD)
{
    double temp= 1+(4*cos(beta/R2D)*cos(beta/R2D)*Q)/((sigma_prime*(CL*sin(beta/R2D)+CD*cos(beta/R2D))));
    double a=1/temp;
    return a;
}

double calculate_a_prime1(double a)
{
    double a_prime = (1-3*a)/(4*a-1);
    return a_prime;
}

double calculate_tip_loss(double beta,double r,double R,double B)
{
    double temp;
    double Q;
    temp=((B/2)*(1-(r/R)))/((r/R)*cos(beta/R2D));
    Q=(2/PI)*acos(exp(-temp));
    if(Q<=0)
    {
        Q=1;
    }
    return Q;
}

double calculate_i(double gamma,double beta)
{
    double i=gamma-beta;
    if(i<0)
    {
        i=-i;
    }
    return i;
}

double calculate_beta(double lemda,double a_prime,double a)
{
    double temp= (lemda*(1+a_prime))/(1-a);
    double beta= atan(temp)*R2D;
    if(beta<0)
            {
                beta=-beta;
            }

    return beta;
}

double calculate_a_prime(double sigma_prime,double CL,double lemda,double beta,double a,double Q)
{
    double temp=(sigma_prime*CL)/(4*lemda*cos(beta/R2D)*Q);
    double a_prime= temp*(1-a);
    return a_prime;
}

double calculate_c(double r,double beta, double B, double lemda_r)
{
    double c = (8*PI*r*cos(beta/R2D))/(3*B*lemda_r);
    return c;
}

double calculate_error(double a1, double a2)
{
    double temp= a2-a1;
    double e= (temp/a1)*100;
    if(e<0)
    {
        e=-e;
    }
    return e;
}

double NACA02CL(double ALPHA)
{
    double CCL3[5]={.7419,-.078,.0205,-.0008,0};
    double CCL1[5]={-2.3054,-1.2862,-.285,-.0247,-.0008};
    double CCL2[5]={.1208,.1025,.0016,.0003,0};
    int VAR1=0;
    double CL=0;

    if(ALPHA <= -5)
    {
        for(VAR1=0;VAR1<5;VAR1++)
        {
            CL=CL+CCL1[VAR1]*pow(ALPHA,VAR1);
        }
        if(CL<-1.5)
        {
            CL=-1.5;
        }
    }
    else if(ALPHA <= 5.5)
    {
        for(VAR1=0;VAR1<5;VAR1++)
        {
            CL=CL+CCL2[VAR1]*pow(ALPHA,VAR1);
        }
    }
    else
    {
        for(VAR1=0;VAR1<5;VAR1++)
        {
            CL=CL+CCL3[VAR1]*pow(ALPHA,VAR1);
        }
    }
    return CL;

}

double NACA25112CL(double alpha)
{
    double c1[5]={0.08473,0.1376,0.0001129,-0.0002109,0};
    int i;
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    return CL;
}

double NACA25112CD(double alpha)
{
    double c1[5]={0.002,0.0043,-0.0008,0.00004,0};
    int i;
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}

double NACA0012CD(double alpha)
{
    double c1[5]={0.01081,-0.005385,0.001663,-0.0001592,5.042*pow(10,-6)};
    int i;
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}

double NACA0015CL(double alpha)
{
    double c1[5]={0.02451,0.07537,0.009332,-0.0007413,1.226*pow(10,-5)};
    int i;
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    return CL;
}
double NACA2415CL(double alpha)
{
    double c1[5]={0.2455,0.1199,- 0.0004,-0.0001,0};
    int i;
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    return CL;
}

double NACA2415CD(double alpha)
{
    double c1[5]={0.0037,0.0028,- 0.0005,0.00003,0};
    int i;
    if (alpha<0)
    {
        alpha=-alpha;
    }
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}

double NACA0015CD(double alpha)
{
    double c1[5]={0.008467,-0.002984,0.001008,-9.957*pow(10,-5),3.455*pow(10,-6)};
    int i;
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}
double NACA63210CL(double alpha)
{
    double c1[5]={0.1688,0.1309,-0.008934,0.001017,-4.245*pow(10,-5)};
    int i;
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    return CL;
}
double NACA63210CD(double alpha)
{
    double c1[5]={0.009519,-0.00731,0.002802,-0.0002977,1.032*pow(10,-5)};
    int i;
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}

double NACA64210CD(double alpha)
{
    double c1[5]={0.006415,-0.003523,0.001852,-0.0002616,1.304*pow(10,-5)};
    int i;
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}

double NACA23018CD(double alpha)
{
    double c1[5]={0.004,0.0033,-0.0005,0.00003,0};
    int i;
    if (alpha<0)
    {
        alpha=-alpha;
    }
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}
double NACA23024CL(double alpha)
{
    double c1[5]={0.08799,-0.008518,0.01961,-0.001152,1.656*pow(10,-5)};
    int i;
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    return CL;
}
double NACA23024CD(double alpha)
{
    double c1[5]={0.009363,-0.0006452,0.0002301,-1.686*pow(10,-5),8.307*pow(10,-7)};
    int i;
    double CD=0;
    for(i=0;i<5;i++)
        {
            CD=CD+c1[i]*pow(alpha,i);
        }
    return CD;
}


double NACA64210CL(double alpha)
{
    double c1[5]={0.1577,0.1581,-0.02288,0.003336,-0.0001634};
    int i;
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    if(CL<0)
    {
        CL=-CL;
    }
    if(CL>2)
    {
        CL=1.933;
    }
    return CL;
}

double NACA23018CL(double alpha)
{
    double c1[5]={0.1243,0.1033,0.0027,-0.0002,0};
    int i;
    if (alpha<0)
    {
        alpha=-alpha;
    }
    double CL=0;
    for(i=0;i<5;i++)
        {
            CL=CL+c1[i]*pow(alpha,i);
        }
    return CL;
}





double Calculate_CP(double lemda_r[],double a[],double a_prime[], int size1,double CL[],double CD[],double beta[])
{

    double f_x[size1];
    double lemda_d[size1];
    double f_x_i[size1];
    double sum=0;
    double cp=0;
    int i=0;
    for(i=0;i<size1;i++)
    {
        f_x[i]=pow(lemda_r[i],3)*a_prime[i]*(1-a[i])*(1-(CD[i]/CL[i])*tan(beta[i]));

    }
    for(i=0;i<size1;i++)
    {
        if(i==size1-1)
        {
            lemda_d[i]=0;
            f_x_i[i]=0;
        }
        else
        {
            lemda_d[i]=(.5)*(lemda_r[i+1]-lemda_r[i]);
            f_x_i[i]=lemda_d[i]*(f_x[i]+f_x[i+1]);
        }
        sum=sum+f_x_i[i];

    }
    cp=8*sum*(1/(lemda_r[size1-1]*lemda_r[size1-1]));
    printf("\nCP = %lf",cp);
    return cp;
 
}
