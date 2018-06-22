#include <stdio.h>
#include <stdlib.h>
#include<math.h>

/* Random number generator in [0, 1] */
#define Rand()		((double)rand()/RAND_MAX)
double optimum_value;
#define NMax 1000
#define RMax 2000	
FILE *f_run1,*f_run;
int total_feval;
double acc_err;
struct landscape funct;
double GlobalMins[RMax]; //r_max=total number of runs
/* Structure of a particle */
typedef struct {
	double *x;
	double *v;
	double f;
	double pbest;
	double *x_star;
} ParticleRec,*Particle;
struct landscape {int N;double x[NMax];double fx[NMax];}; // 1D-landscape
/*
	Parameters for PSO
*/
int Pr;
/* Number of particles.this is the swarm size. Good if this lies between 10 and 30 */
#define Nparticles	50
/* Maximum number of iterations */
#define T_MAX		100
/* The value of inertia weight at t=0 (W_0) and t=T_MAX (W_T) . THis is increasing inertia.*/
#define W_0		0.9
#define W_T		0.4
#define MAX_V		2.0
/* The cognitive parameter (c1) and the social parameter (c2) */
#define c1		2.0
#define c2		2.0

/*
	Definitions for a problem
*/

/* Number of variables: problem dependent . */
#define Nvariables	10

/* Objective function for minimization: problem dependent */
#define better(y1, y2)	(y1<y2)

void fun(Particle P)
   {
     
    
    long  double f=0;
    
    
    switch (Pr)
    {
      case 0: // Schwefel prob 1.2 (in list 1)
	{ 

	int i,j;
       double s1=0,f;
	//optimum_value=0;
	P->f = 0.0;                                    
       double s2=0;
       for(i=0;i < Nvariables;i++)
       {
          for(j=0;j < i;j++)
          {
                 s1+=P->x[j];
          }
          P->f+=s1*s1;
	if((P->f)!=0)
			total_feval++;
       }
       break;
       }
	case 1:	//step
	{
	int i;
	double f;
	//optimum_value=0;
	P->f = 0.0;
	for(i=0;i < Nvariables;i++){
	P->f+=pow(floor(P->x[i]+0.5),(double)2);
	}
	if((P->f)!=0)
			total_feval++;
	break;
	}  
	case 2:
	 {
       double s1=0;                                    
         int i;
	optimum_value=-418.9829*Nvariables;
       for(i=0;i < Nvariables;i++)
       {
          s1=s1+P->x[i]*sin(sqrt(fabs(P->x[i])));
          
       }
       P->f=-s1;
	if((P->f)!=0)
			total_feval++;
       break;
     }
	case 3:
	 {
	int k = 10;
	int i,top = 0;
	//optimum_value=0;
	double xd,pi;
        //f = 0;
        for(i=0;i < Nvariables;i++)
        {
          xd = P->x[i];
          P->f += P->x[i] *P->x[i] - k * cos( 2 * pi * P->x[i] );
        }
        P->f += P->x[i] * k;
	if((P->f)!=0)
			total_feval++;
      break;
	 } 
	case 4:
	{	                                  
   
double s =0;
int i;
//optimum_value=0;
P->f = 0.0;
for(i=0;i < Nvariables-1;i++){
s+=((1+0.25*(P->x[i]+1))-1)*((1+0.25*(P->x[i]+1))-1)*(1+10*sin(3.14*(1+0.25*(P->x[i+1]+1)))*sin(3.14*(1+0.25*(P->x[i+1]+1))));
P->f = (3.14/(Nvariables))*(10*sin(3.14*(1+0.25*(P->x[0]+1)))*sin(3.14*(1+0.25*(P->x[0]+1)))+s+((1+0.25*(P->x[i-1]+1))-1)*((1+0.25*(P->x[i-1]+1))-1));
}
if((P->f)!=0)
			total_feval++;
break;		
}
	case 5:
	{
	int i;	                                                        
double s =0;

P->f = 0.0;
for(i=0;i < Nvariables-1;i++){
          s+=(P->x[i]-1)*(P->x[i]-1)*(1+sin(3*3.14*P->x[i+1])*sin(3*3.14*P->x[i+1]));
         P->f = 0.1*(sin(3*3.14*P->x[0])*sin(3*3.14*P->x[0])+s+(P->x[Nvariables-1]-1)*(P->x[Nvariables-1]-1)*(1+sin(2*3.14*P->x[Nvariables-1])*sin(2*3.14*P->x[Nvariables-1])));
 	//return top;	 
		}
	if((P->f)!=0)
			total_feval++;
	break;
	}	
	case 6:
	{ 	
	//optimum_value=-10.1532;
		static double c[5] =
            { 		
		0.1,0.2,0.2,0.4
//,0.3,0.7,0.5,0.5,0.4
                   
            };	
		
		static double a[5][Nvariables] = {
				3.0,7.0,2.0,9.0,5.0,45,44,0.2,2.5
//,3.0,8.0,1.0,6.0,2.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = P->x[0];
			float x2 = P->x[1];
			float x3 = P->x[2];
			float x4 = P->x[3];
			//float x5 = x[4];
			//float x6 = x[5];
		int s= 0, i,j;
		
		for(j=0;j<5;j++){
			int p = 0;
		for(i=0;i<Nvariables;i++){
			p = p + ((x1*x2*x3*x4)-a[j][i])*((x1*x2*x3*x4)-a[j][i]);
		}
		s = s + 1/(p+c[j]);
		}
			P->f = -s;
			if((P->f)!=0)
			total_feval++;	
		break;
	}
case 7:
{
            int x1=P->x[0],i;
            int x2=P->x[1];
		//optimum_value=0.000307486;
		
		P->f = 0.0;	
	      int x3=P->x[2];
            int x4=P->x[3];
           	double a[]={0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235,0.0246};
           	double b[]={4.0, 2.0,1.0, 0.5, 0.25, 0.1667, 0.125, 0.1, 0.0833,  0.07143, 0.0625};
            for(i=0;i<Nvariables;i++)
            {
                P->f += pow((P->x[i]-(x1*(b[i]*b[i]+b[i]*x2))/(b[i]*b[i]+b[i]*x3+x4)),(double)2);
            }
           if((P->f)!=0)
			total_feval++;
        break;
          
		  }
case 8:
{
          
	//optimum_value=-1.0316;
	P->f = 0.0;
           float x1=P->x[0];
           float x2=P->x[1];
           P->f=4*x1*x1-2.1*pow(x1,4)+(pow(x1,6)/3.0)+x1*x2-4*x2*x2+4*pow(x2,4);
		if((P->f)!=0)
			total_feval++;
              break;
       }
case 9:
{	double pi=3.14;
           	//double a,b,c,d,e,fav;
           	double a=1; 
		//optimum_value=0.3979;
		//long double r = ;
           	double b=5.1/(4*pi*pi);
           	double c=5/(pi);
           	double d=6;
           	double e=10;
           	double fav=1/(8*(pi));
		//P->f = 0.0;
          int x1=P->x[0];
      int  x2=P->x[1];
            
            P->f = (1*pow((x2-(b*x1*x1)+(c*x1)-6),2)+(10*(1-fav)*cos(x1)+10));
	if((P->f)!=0)
			total_feval++;
           break;
    }
case 10:
{
        double f;
	optimum_value=3;
	P->f = 0.0;
          float x1=P->x[0];
          float x2=P->x[1];
           P->f=(1+pow((x1+x2+1),(double)2)*(19-14*x1+3*x1*x1-14*x2+6*x1*x2+3*x2*x2))* (30+pow((2*x1-3*x2),(double)2)*(18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2));
	if((P->f)!=0)
			total_feval++;
           	break;
       }
case 11:
{ 
//optimum_value=-3.86278;
		static double p[4][3] =
            { 		0.36890,0.11700,0.26730,
			0.46990,0.43870,0.74700,0.10910,0.87320,0.55470,
			0.03815,0.57430,0.88280
                   
            };	
		static double c[4]=
		{
			1.0,1.2,3.0,3.2

		};
		static double a[4][3] = {
				3.0,0.1,30.0,35.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = P->x[0];
			float x2 = P->x[1];
			float x3 = P->x[2];
		int s= 0, i,j;
		double f;
		for(i=0;i<4;i++){
			int sm = 0;
		for(j=0;j<3;j++){
			sm = sm + a[i][j]*(((x1*x2*x3)-p[i][j])*((x1*x2*x3)-p[i][j]));
		}
		s = s + c[i] * exp(-sm);
		}
			P->f = -s;
		if((P->f)!=0)
			total_feval++;	
		break;	
	}
case 12:
{ 	optimum_value=-10.4029;
		static double c[7] =
            { 		
		0.1,0.2,0.2,0.4,0.3,0.7,0.5
//,0.5,0.4
                
            };	
		
		static double a[7][4] = {
				3.0,7.0,2.0,9.0,5.0,3.0,8.0
//,1.0,6.0,2.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = P->x[0];
			float x2 = P->x[1];
			float x3 = P->x[2];
			float x4 = P->x[3];
			//float x5 = x[4];
			//float x6 = x[5];
		int s= 0, i,j;
		double f;
		for(j=0;j<7;j++){
			int p = 0;
		for(i=0;i<4;i++){
			p = p + ((x1*x2*x3*x4)-a[j][i])*((x1*x2*x3*x4)-a[j][i]);
		}
		s = s + 1/(p+c[j]);
		}
			P->f = -s;
		if((P->f)!=0)
			total_feval++;	
		break;
	}
case 13:
{ 	
//optimum_value=-10.5364;
		static double c[10] =
            {		
		0.1,0.2,0.2,0.4,0.3,0.7,0.5,0.5,0.4,0.8
                
            };	
		
		static double a[10][4] = {
				3.0,7.0,2.0,9.0,5.0,3.0,8.0,1.0,6.0,2.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = P->x[0];
			float x2 = P->x[1];
			float x3 = P->x[2];
			float x4 = P->x[3];
			//float x5 = x[4];
			//float x6 = x[5];
		int s= 0, i,j;
		double f;
		for(j=0;j<7;j++){
			int p = 0;
		for(i=0;i<4;i++){
			p = p + ((x1*x2*x3*x4)-a[j][i])*((x1*x2*x3*x4)-a[j][i]);
		}
		s = s + 1/(p+c[j]);
		}
			P->f = -s;
		if((P->f)!=0)
			total_feval++;	
		break;	
		
	}
case 14:
{ 	
		//optimum_value=-3.32237;
		static double p[4][7] =
            { 		0.36890,0.11700,0.26730,
			0.46990,0.43870,0.74700,0.10910,0.87320,0.55470,
			0.03815,0.57430,0.88280
                   
            };	
		static double c[4]=
		{
			1.0,1.2,3.0,3.2

		};
		static double b[4][3] = {
				3.0,0.1,30.0,35.0,25,33
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = P->x[0];
			float x2 = P->x[1];
			float x3 = P->x[2];
			float x4 = P->x[3];
			float x5 = P->x[4];
			float x6 = P->x[5];
		int o= 0, i,j;
		
		for(i=0;i<4;i++){
			int sm1 = 0;
		for(j=0;j<3;j++){
			sm1 = sm1 + b[i][j]*(((x1*x2*x3*x4*x5*x6)-p[i][j])*((x1*x2*x3*x4*x5*x6)-p[i][j]));
		}
		o = o + c[i] * exp(-sm1);
		//printf("%d \n",o);
		}
			P->f = -o;
		if((P->f)!=0)
			total_feval++;	
		break;	
	}
case 15:
	{

P->f = 0.0;
double s=0;
int i;          //Cigar (madam paper)
for(i=1;i < Nvariables;i++){
s+=P->x[i]*P->x[i];
P->f=P->x[0]*P->x[0]+100000*s;
}
if((P->f)!=0)
			total_feval++;
break;
}
case 16:
{
        double f;
	int i;
	optimum_value=0;
	P->f = 0.0;
        for ( i = 0; i < Nvariables; i++ )
        {
           P->f +=(i+1)*P->x[i] * P->x[i];
        }
	if((P->f)!=0)
			total_feval++;
      break;
 }
case 17:
{		
      	double f;
		optimum_value=0;
           	float x1=P->x[0];
           float  x2=P->x[1];	
	
        
            	P->f = pow((1.5-x1*(1-x2)),(double)2)+pow((2.25-x1*(1-x2*x2)),(double)2)+pow((2.625-x1*(1-x2*x2*x2)),(double)2);
		if((P->f)!=0)
			total_feval++;
          break;
    }
case 18:
{	
//optimum_value=-450;
          	static double offset_0[30] =
            { 
                   -3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
                   -8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000, 
                   -1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
                   6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001, 
                   3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001, 
                   -6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
            };		
		
		P->f = -450.0;
		double xd;
		int i;
	    for (i=0;i<Nvariables;i++)
        {
			xd = P->x[i]-offset_0[i];
			P->f += xd * xd;  
		}
	if((P->f)!=0)
			total_feval++;
		break;
    }
case 19:
{		
//optimum_value=-450;
       // Shifted Schwefel (F2 CEC 2005)
       static double offset_4[30] =
       { 
              3.5626700e+001, -8.2912300e+001, -1.0642300e+001, -8.3581500e+001,  8.3155200e+001,
              4.7048000e+001, -8.9435900e+001, -2.7421900e+001,  7.6144800e+001, -3.9059500e+001,
              4.8885700e+001, -3.9828000e+000, -7.1924300e+001,  6.4194700e+001, -4.7733800e+001,
              -5.9896000e+000 ,-2.6282800e+001, -5.9181100e+001,  1.4602800e+001, -8.5478000e+001,
              -5.0490100e+001,  9.2400000e-001,  3.2397800e+001,  3.0238800e+001, -8.5094900e+001,
              6.0119700e+001, -3.6218300e+001, -8.5883000e+000, -5.1971000e+000,  8.1553100e+001 
       };

		for (int i=0;i<Nvariables;i++)
		{
			P->x[i]=P->x[i]-offset_4[i];
		}

    
	P->f = -450.0;
	//int i,j;
    for (int i=0;i<Nvariables;i++)
    {
       float  sum2 = 0.0;
        for (int j=0; j<=i; j++)
        {
            sum2 += P->x[j];
        }
        P->f += sum2*sum2;
    }
	if((P->f)!=0)
			total_feval++;
		break;
  }

case 20:
{		//optimum_value=-180;
       // Shifted Griewank (CEC 2005)
       static double offset_5[30] =
       { 
         -2.7626840e+002, -1.1911000e+001, -5.7878840e+002, -2.8764860e+002, -8.4385800e+001,
         -2.2867530e+002, -4.5815160e+002, -2.0221450e+002, -1.0586420e+002, -9.6489800e+001,
         -3.9574680e+002, -5.7294980e+002, -2.7036410e+002, -5.6685430e+002, -1.5242040e+002,
         -5.8838190e+002, -2.8288920e+002, -4.8888650e+002, -3.4698170e+002, -4.5304470e+002,
         -5.0658570e+002, -4.7599870e+002, -3.6204920e+002, -2.3323670e+002, -4.9198640e+002,
         -5.4408980e+002, -7.3445600e+001, -5.2690110e+002, -5.0225610e+002, -5.3723530e+002 
       };
	double xd;
	 float  sum1 = 0.0;
	 float  sum2 = 1.0;
	 double  f;
		P->f = 180.0;
	   for (int i=0;i<Nvariables;i++)
		 {
		    xd=P->x[i]-offset_5[i];
				sum1 += xd*xd;
        sum2 *= cos(xd/sqrt(1.0+i));
 		 }
    P->f +=1.0 + sum1/4000.0 - sum2;
	if((P->f)!=0)
			total_feval++;
		break;
  }

case 21:
{	//optimum_value=-140;
  // Shifted Ackley (CEC 2005)
  static double offset_6[30] =
  { 
    -1.6823000e+001,  1.4976900e+001,  6.1690000e+000,  9.5566000e+000,  1.9541700e+001,
    -1.7190000e+001, -1.8824800e+001,  8.5110000e-001, -1.5116200e+001,  1.0793400e+001,
    7.4091000e+000,  8.6171000e+000, -1.6564100e+001, -6.6800000e+000,  1.4543300e+001,
    7.0454000e+000, -1.8621500e+001,  1.4556100e+001, -1.1594200e+001, -1.9153100e+001,
    -4.7372000e+000,  9.2590000e-001,  1.3241200e+001, -5.2947000e+000,  1.8416000e+000,
    4.5618000e+000, -1.8890500e+001,  9.8008000e+000, -1.5426500e+001,  1.2722000e+000
  };
    double xd,E;
double pi = 3.14;
	P->f = 140.0;
   float sum1 = 0.0;
   float sum2 = 0.0;
    for (int i=0;i<Nvariables;i++)
    {
        xd = P->x[i]-offset_6[i];
				sum1 += xd*xd;
        sum2 += cos(2.0*pi*xd);
    }
    sum1 = -0.2*sqrt(sum1/Nvariables);
    sum2 /= Nvariables;
    P->f +=20.0 + E - 20.0*exp(sum1) - exp(sum2);
if((P->f)!=0)
			total_feval++;
		break;
}
case 22:
		
       {	double pi = 3.14;
              
		 //optimum_value=-1;
		P->f = 0.0;
             float  x1=P->x[0];
             float  x2=P->x[1];
               P->f=-cos(x1)*cos(x2)*exp((-pow((x1-pi),(double)2)-pow((x2-pi),(double)2)));
if((P->f)!=0)
			total_feval++;
              break;
       }
case 23:
		
       {
              
		//optimum_value=-24777;
              float x1=P->x[0];
              float x2=P->x[1];
               P->f=1.0e5*x1*x1+x2*x2-pow((x1*x1+x2*x2),2)+1.0e-5*pow((x1*x1+x2*x2),4);
if((P->f)!=0)
			total_feval++;
               break;
       }

}}
    

void initialize_params(int Pr)
{
 
 if (Pr==0)
        {
           
           optimum_value=0;
           acc_err=1.0e-3;
           
        }
       
        if (Pr==1)
        {
          
           optimum_value=0;
           acc_err=1.0e-3;
           
        }
       if (Pr==2)
        {    
           optimum_value=-418.9829*Nvariables;
           acc_err=1.0e-3;
        }
	if (Pr==3)
	  {         
           optimum_value=0;
           acc_err=1.0e-3;          
        }
       if (Pr==4)
	  {         
           optimum_value=0;
           acc_err=1.0e-3;          
        }
	if (Pr==5)
	  {         
           optimum_value=0;
           acc_err=1.0e-3;          
        }
	if (Pr==6)
	  {         
           optimum_value=0.998;
           acc_err=1.0e-3;          
        }
	if (Pr==7)
	  {         
           optimum_value=0.0003075;
           acc_err=1.0e-3;          
        }
       if (Pr==8)
	  {         
           optimum_value=-1.0316;
           acc_err=1.0e-3;          
        }
	if (Pr==9)
	  {         
           optimum_value=0.397887;
           acc_err=1.0e-3;          
        }
    if (Pr==10)
	  {         
           optimum_value=3;
           acc_err=1.0e-3;          
        }
	if (Pr==11)
	  {         
           optimum_value=-3.86278;
           acc_err=1.0e-3;          
        }
	if (Pr==12)
	  {         
           optimum_value=-3.32237;
           acc_err=1.0e-3;          
        }
	if (Pr==13)
	  {         
           optimum_value=-3.32237;
           acc_err=1.0e-3;          
        }
	if (Pr==14)
	  {         
           optimum_value=-10.1532;
           acc_err=1.0e-3;          
        }
	if (Pr==15)
	  {         
           optimum_value=-10.4029;
           acc_err=1.0e-3;          
        }
	if (Pr==16)
	  {         
           optimum_value=-10.5364;
           acc_err=1.0e-3;          
        }
	if (Pr==17)
	  {         
           optimum_value=0;
           acc_err=1.0e-5;          
        }
	if (Pr==18)
	  {         
           optimum_value=0;
           acc_err=1.0e-5;          
        }
	if (Pr==19)
	  {         
           optimum_value=0;
           acc_err=1.0e-5;          
        }
	if (Pr==20)
	  {         
           optimum_value=-450;
           acc_err=1.0e-5;          
        }
	if (Pr==21)
	  {         
           optimum_value=-450;
           acc_err=1.0e-5;          
        }
	if (Pr==22)
	  {         
           optimum_value=-180;
           acc_err=1.0e-5;          
        }
	if (Pr==23)
	  {         
           optimum_value=-2477;
           acc_err=1.0e-5;          
        }
	
}

/* update pbest */
void UpdateBest(Particle P)
{
	int j;

	for(j=0; j<Nvariables; j++) 
	P->x_star[j]=P->x[j];
	P->pbest=P->f;
	
}

/* Initialization of particles: problem dependent */
/* The function returns the index of the best particle */
int Initialize(Particle P, int n)
{
	int i, j;
	int G;		/* the index of the best particle */

	G=0;
	for(i=0; i<n; i++) {
		for(j=0; j<Nvariables; j++) {
			P[i].x[j]=Rand();	/* problem dependent */
			P[i].v[j]=0.0;		/* problem dependent */
		}
int p;
		
		fun(&P[i]);
		UpdateBest(&P[i]);
		if(better(P[i].f, P[G].f)) G=i;
		
	}
	return G;
	
}

/*
	Definitions for PSO
*/

/* allocate new data structures.purane velocities ko rakhna zaroori hai cauz if the next value is not good then it retains the previous value*/
#define New(type, n, msg)	(type *)NewCell(sizeof(type), n, msg)

void *NewCell(int size, int n, char *msg)
{
	void *new;

	if((new=malloc(size*n))==NULL) {
		fprintf(stderr, "Cannot allocate memory for %d %s\n", n, msg);
		exit(1);

	}
	return new;
}

/* allocate "n" new particles */
Particle NewParticles(int n)
{
	int i;
	Particle P;

	P=New(ParticleRec, n, "particles");
	for(i=0; i<n; i++) {
		P[i].x=New(double, Nvariables, "x");
		P[i].v=New(double, Nvariables, "v");
		P[i].x_star=New(double, Nvariables, "x*");
	
	}
	return P;
}

/* Print a particle */
void Print(Particle P)
{
	int j;

	//for(j=0; j<Nvariables; j++)
		//printf("%f ", P->x_star[j]);
	//printf(" = %e", P->pbest);
	printf("%e \n", P->f);
	
}

/* Particle Swarm Optimization */
void main()
{
f_run1 = fopen( "PSO1_op.ods", "w" ); 
f_run = fopen( "PSO2_op.ods", "w" );       
fprintf(f_run,"  Mean Fun Val   \t Std_Dev \t   Error  \t RUN \t Succ Rate \t Acceptable error\n  ");
fprintf(f_run1,"  RUN  \t Fun Val \t   Error \n  ");

	int t, i, j;
	Particle P;
	int G,p;
	double w;
	
	for(Pr=0;Pr<24;Pr++)
	{
		initialize_params(Pr);
	double mean = 0.0,error=0,mean_error=0,var=0.0,sd=0.0;
	int total_feval;
	double mean_feval=0;
	int succ_rate;
	P=NewParticles(Nparticles);
	G=Initialize(P, Nparticles);
	w=W_0;
	double feval=0;
	for(t=0; t<T_MAX; t++) {
		for(i=0; i<Nparticles; i++) {
			for(j=0; j<Nvariables; j++) {
				P[i].v[j]=w*P[i].v[j]+c1*Rand()*(P[i].x_star[j]-P[i].x[j])+c2*Rand()*(P[G].x_star[j]-P[i].x[j]);
				if(P[i].v[j]<-MAX_V)
					P[i].v[j]=-MAX_V;
				else if(P[i].v[j]>MAX_V)
					P[i].v[j]=MAX_V;
				P[i].x[j]+=P[i].v[j];
			}
			//printf("enter the function you want");
			//scanf("%c",&function);


			fun(&P[i]);
			if(better(P[i].f, P[i].pbest)) {       /*each particle ka pbest*/
				if(better(P[i].f, P[G].pbest)) G=i; /* best particle ke index ke saath compare karra hai*/
				UpdateBest(&P[i]);
				
			}
			
	if(fabs((P->f)-optimum_value)<=acc_err)
    {
        succ_rate+=1;
        mean_feval=mean_feval+feval;
	
        break;
        
    }//if GlobalMin*/
}
	error=fabs((P->f)-optimum_value);
	mean_error=mean_error+error;
	mean=mean+(P->f);
	
fprintf(f_run1,"%d \t  %e \t %e \n",t+1,(P->f),mean_error);
printf("function %d Dimension %d  Run %d  F_val %e Error %e \n", Pr,Nvariables, t,(P->f),mean_error);
		

		w-=(W_0-W_T)/T_MAX;
		
	}
 GlobalMins[t]=P->f;
mean=mean/T_MAX;
mean_error=mean_error/T_MAX;

 for(int k=0;k<T_MAX;k++)
    var=var+pow(GlobalMins[k]-mean,2);
    var=var/T_MAX;
    sd=sqrt(var);
printf("Means of %d runs: %e \n",T_MAX,mean);
fprintf(f_run," %e \t %e \t %e   \t %d \t %d  \t %e \n", mean,sd, mean_error,T_MAX, succ_rate,optimum_value );
}
fclose(f_run);
fclose(f_run1);
printf("%d total f_eval",total_feval);
}



