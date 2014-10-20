
#include <stdio.h>
#include <complex.h>
#include <math.h>

#define N	100
#define pi	3.141592654
#define h	657.9276

#define z0	complex(0.0,0.0)
#define zi	complex(0.0,1.0)

/* using namespace std; */


double T2;
double chi0, delta_t, detuning, de ; /* parameters for laser */
double gg[N+1][N+1];
complex pt[N+1][2];

main()
{

 /*  declare variables  */

     int j, k, ks, nt,nt0;
     double t, dt, ti, tmax;
     double aB, Er, emax, g_density, g_hf, g_pola;
     double sumf;
     complex sump;
     complex p[N+1][2];
     double can[N+1];


	  void ff(double t, complex p[][2], complex pt[][2]);
     void RK4(double t, double dt, complex p[][2],
    		     void (ff)(double t, complex p[][2], complex pt[][2]));


     FILE *outputp;
     FILE *outputf;

    /* parameters */
   chi0=0.1;       /* for intensity of pulse */
   delta_t= 20.0;	/*width of laser pulse  */
   detuning=0.0 ;       /*Delta_T */
   aB=125.0e-8; 	/* Bohr radius (cm) */
   Er=4.2;        /* energy  */
   emax=300.0;     /* max energy */
   T2=200.0;  /* dephasing time */
   nt0=5;     /* time output step  */

/*  define initial values  */
     ti=-3*delta_t; /* initial time or t_0 */
     dt=2.0;    /* time step */
     tmax=12*delta_t;

/* computing parameters */

   de=emax/N;
   g_hf=sqrt(Er*emax/N)/pi/h; /* for hartree-fock term */
   g_density=0.5*(de*sqrt(de))/(Er*sqrt(Er))/((aB*aB*aB)*(pi*pi))/10e+16;
   g_pola=de*sqrt(de)/(Er*sqrt(Er))/pi/pi;
   delta_t=delta_t*0.84932 ;
   chi0=sqrt(pi)*chi0/delta_t ;

   de=de/h;
   detuning=detuning/h;



 /* open and print out to file  */

     outputp=fopen("fpRK", "w");
     outputf=fopen("deRK", "w");

   for (k=1;k<=N;k++)   can[k]= sqrt(k);


/* pre-computing g(k,k') in Hartree-Fock terms */

   for (k=1;k<=N;k++){
     for (ks=1;ks<=N;ks++){
        if(k==ks)
        {
        	 gg[k][ks]=0.0;
        }else if(ks<k)
        {
          gg[k][ks]= (g_hf/can[k])*
          			log(double ( (can[k]+can[ks])/ (can[k]-can[ks]) ) );
        }else {
          gg[k][ks]=  (g_hf/can[k])*
  				log(double ( (can[k]+can[ks])/ (can[ks]-can[k]) ) );
        }
     }
   }


/* at t=t_0=ti */

     t=ti;           
     for (k=1;k<=N;k++) p[k][0]=p[k][1]=z0;

		fprintf(outputf, "%10.5f\n", tmax);
     for (k=1;k<=N;k=k+4){
     	fprintf(outputp, "%10.5f %10.5f %10.5f %10.5f %10.5f\n",
      	t,k*de*h, real(p[k][0]),imag(p[k][0]), abs(p[k][1]) );
     }
      fprintf(outputp, "\n");


/*  loop step time by dt  */
/*  call RK4 to step x[i] */
/*  print out to file     */
        nt=nt0;
	for (j=1; t<tmax; j++) {

         RK4(t, dt, p, ff);
         t=j*dt+ti;

			sump=z0;
         sumf=0.0;
       for (k=1;k<=N;k++){
         sump += can[k]*p[k][1]*g_pola;
         sumf += can[k]*real(p[k][0]);
       }
         fprintf(outputf, "%10.5f %10.5f %10.5f\n",
          			t, abs(sump),g_density*sumf);

   	if(j==nt){
       for (k=1;k<=N;k=k+4){
          fprintf(outputp, "%10.5f %10.5f %10.5f %10.5f %10.5f\n",
          			t,k*de*h, real(p[k][0]),imag(p[k][0]),abs(p[k][1]) );
       }
	       fprintf(outputp, "\n");
				nt+=nt0;
      }
	}

/* close output file and end */
     fclose(outputp);
     fclose(outputf);
}

/*------------------------------------------------------------*/
/*                                                            */
/*               FUNCTIONs for RUNGE-KUTTA                     */
/*                                                            */
/*               solves p' = dp/dt = pt(t,x)                   */
/*               for p=p[i] with i=1,1,...,N                */
/*                                                            */
/*------------------------------------------------------------*/

void RK4(double t, double dt, complex p[][2],
			void (*ff)(double t, complex p[][2], complex pt[][2]))
{
     int i;
     double dt2=dt/2.0;
     complex ptemp[N+1][2], kp1[N+1][2], kp2[N+1][2];
/*     void ff(double t, complex p[], complex pt[]); */

     (*ff)(t,p,pt);
     for (i=1;i<=N;i++)                  /*  k1 calculation  */
	  {
	     		kp1[i][0]= dt*pt[i][0];
	     		kp1[i][1]= dt*pt[i][1];
   	      ptemp[i][0] = p[i][0]+ 0.5*kp1[i][0];
   	      ptemp[i][1] = p[i][1]+0.5*kp1[i][1];
	  }

	  (*ff)(t+dt2,ptemp,pt);
     for (i=1;i<=N;i++)                  /*  k2 calculation  */
	  {
				kp2[i][0]= dt*pt[i][0];
				kp2[i][1]= dt*pt[i][1];
   			ptemp[i][0] = p[i][0]+ 0.5*kp2[i][0];
   			ptemp[i][1] = p[i][1]+ 0.5*kp2[i][1];
            kp1[i][0] += 2*kp2[i][0];
            kp1[i][1] += 2*kp2[i][1];
	  }

     (*ff)(t+dt2,ptemp,pt);
     for (i=1;i<=N;i++)                  /*  k3 calculation  */
	  {
				kp2[i][0] = dt*pt[i][0];
				kp2[i][1] = dt*pt[i][1];
   			ptemp[i][0]  = p[i][0]+ 0.5*kp2[i][0];
   			ptemp[i][1]  = p[i][1]+ 0.5*kp2[i][1];
            kp1[i][0] += 2*kp2[i][0];
            kp1[i][1] += 2*kp2[i][1];
	  }

	  (*ff)(t+dt,ptemp,pt);
     for (i=1;i<=N;i++)                  /*  k4 calculation  */
	  {
            kp1[i][0] += dt*pt[i][0];
            kp1[i][1] += dt*pt[i][1];
     }

     for (i=1;i<=N;i++)                      /*  step x  */
	  {
             p[i][0] += kp1[i][0]/6.0;
             p[i][1] += kp1[i][1]/6.0;
     }

}

/*------------------------------------------------------------*/
/*                                                            */
/*               FUNCTION for DERIVATIVES                     */
/*                                                            */
/*               defines pt(t,p)= dp/dt                       */
/*               for p=p[i] with i=1,1,...,N                  */
/*                                                            */
/*------------------------------------------------------------*/


void ff(double t, complex p[][2], complex pt[][2])
{
	int k1,k2;
	double chi, ek;
	complex rabi, e_renorm;
	complex rabi0, f_hf, p_hf;

		chi=chi0*exp(-(t*t)/(delta_t*delta_t));	/* chi=chi0/cosh(t/delta_t) */
		rabi0=zi*0.5*chi;

	for (k1=1;k1<=N;k1++) {
	     	f_hf=z0;
			p_hf=z0;
			/*	---- H-F --- */
         for (k2=1; k2<=N; k2++){
	          f_hf+= gg[k1][k2]*p[k2][0];
   		    p_hf+= gg[k1][k2]*p[k2][1];
         }

       ek = de*k1-(real(f_hf)+imag(f_hf))-detuning;
       e_renorm = -zi*ek*p[k1][1];
       rabi= rabi0 + p_hf;

     pt[k1][0]=-2*complex(imag(rabi*conj(p[k1][1])),imag(rabi*conj(p[k1][1])));
     pt[k1][1]= e_renorm + zi*(1.0-real(p[k1][0])-imag(p[k1][0]))*rabi-
     				 p[k1][1]/T2;
   }

}


