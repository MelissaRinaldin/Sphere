/******************************************************************
**                                                               **
**                           SPHERE                              **
**                                                               **
** Program for sphere (April 2016)                               **
**                                                               **
*******************************************************************/

/* x=theta (0:PI) and y= phi (0:2PI)  */


#include "membrane.h"

#define  sign(x) ((x < 0) ? -1 : 1)
#define  min(x,y) ((x < y) ? x : y)
#define  mod(x,y) ((( x & (y - 1)) + y) & ( y - 1 ))

time_t t1, t2;

long num_of_meshpoint_theta, num_of_meshpoint_phi;

long num_of_iteration;
long method;

double lagrange;
double epsilon=0.1;
double sigma;

long seed;

double **mu;
double **phi;
double **energy;

double ***k0;
double ***k1;
double ***k2;
double ***k3;
double ***k4;

double dtheta;
double dphi;
double DT;

//double get_area();
double ran2(long *);

double DD(double **, long, long);
double DX(double **, long, long);
double DY(double **, long, long);
double DXX(double **, long, long);
double DYY(double **, long, long);
double DXY(double **, long, long);

void rk2();
void rk4();
void run();
void end();
void init();
void debug();
void euler();
void one_step();
//void get_energy();
void init_random();
void init_import();
void get_chemical();
void allocate_memory();
void progress_bar(long);
void get_rhs(double ***);
void write_hi(FILE *, long);
void export_conf(long, long);
void get_time(time_t, CPU_Time *);

/*******************************************************************/

main() 
{	
	init();
	time(&t1);
	run();
	end();
	
  	return EXIT_SUCCESS;	
}

/*******************************************************************/

void init()
{
	double run_time ;
	char f_name[32];
	long cflag;

	

	printf("Input run time "); 
	scanf("%lg",&run_time); //total time for the simulation

  	printf("Input grid size theta ");
    scanf("%ld",&num_of_meshpoint_theta);
    
    printf("Input grid size phi (only half of them) ");
    scanf("%ld",&num_of_meshpoint_phi);
    
    
    printf("Input iterations ");
  	scanf("%ld",&num_of_iteration);
    
    dtheta= PI/num_of_meshpoint_theta;
    dphi=PI/(num_of_meshpoint_phi);
    
    DT = run_time/num_of_iteration;

	printf("Choose integration method "); 
  	printf("\n");
  	printf("\n\t1. Euler");
  	printf("\n\t2. Runge-Kutta (2 steps)");
  	printf("\n\t3. Runge-Kutta (4 steps)");
  	printf("\n\n");
  	printf("Input flag ");
  	scanf("%ld",&method);
  
	printf("Input sigma "); 
  	scanf("%lg",&sigma);

	printf("Choose intial configuration "); 
  	printf("\n");
  	printf("\n\t1. Random");
	printf("\n\t2. Import");
	printf("\n\t3. Debug");
  	printf("\n\n");
  	printf("Input flag ");
  	scanf("%ld",&cflag);

  	switch(cflag){
  	case 1:
	    init_random();
	    break;
  	case 2:
		init_import();
		break;
	case 3:
		debug();
		break;
  	}
}

/*******************************************************************/

void allocate_memory()
{
	long i, j;
	
	mu = (double **) malloc(num_of_meshpoint_theta*sizeof(double));
	phi = (double **) malloc(num_of_meshpoint_theta*sizeof(double));
	energy = (double **) malloc(num_of_meshpoint_theta*sizeof(double));
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		mu[i] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
		phi[i] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
	}
	
	k0 = (double ***) malloc(4*num_of_meshpoint_theta*sizeof(double));
	k1 = (double ***) malloc(4*num_of_meshpoint_theta*sizeof(double));
	k2 = (double ***) malloc(4*num_of_meshpoint_theta*sizeof(double));
	k3 = (double ***) malloc(4*num_of_meshpoint_theta*sizeof(double));
	k4 = (double ***) malloc(4*num_of_meshpoint_theta*sizeof(double));
	
	for (i=0; i<4; i++){
		k0[i] = (double **) malloc(num_of_meshpoint_theta*num_of_meshpoint_theta*sizeof(double));
		k1[i] = (double **) malloc(num_of_meshpoint_theta*num_of_meshpoint_theta*sizeof(double));
		k2[i] = (double **) malloc(num_of_meshpoint_theta*num_of_meshpoint_theta*sizeof(double));
		k3[i] = (double **) malloc(num_of_meshpoint_theta*num_of_meshpoint_theta*sizeof(double));
		k4[i] = (double **) malloc(num_of_meshpoint_theta*num_of_meshpoint_theta*sizeof(double));
	}
	
	for (i=0; i<4; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			k0[i][j] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
			k1[i][j] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
			k2[i][j] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
			k3[i][j] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
			k4[i][j] = (double *) malloc(2*num_of_meshpoint_phi*sizeof(double));
		}
	}
}

/*******************************************************************/

void init_random()
{
	double ns=0.1, area;
	long i, j;
	
	printf("Input area percentage ");
	scanf("%lg",&area);
	
	printf("Input random seed ");
    scanf("%ld",&seed);

	allocate_memory();
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = (2*area-1)+ns*(1-2*ran2(&seed));
		}
	}
}

/******************************************************************/

void init_import()
{
	printf("This function is not yet implemented\n");
	exit(0);
}

/*******************************************************************/

void get_chemical()
{
	long i, j, num_of_gridpoint=2*num_of_meshpoint_theta*num_of_meshpoint_phi;
	double u;
		
	lagrange = 0;	
		
    for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			u = phi[i][j]*(phi[i][j]*phi[i][j]-1)/(epsilon*epsilon);
			mu[i][j] = sigma*(u-DD(phi,i,j));
			lagrange += u/num_of_gridpoint;
            printf("Derivata %f %f %f %f \n",DD(phi,i,j),DX(phi,i,j),DXX(phi,i,j),DYY(phi,i,j));
		}
	}
}

/*******************************************************************/

void get_rhs(double ***rhs)
{
	long i, j;
	
	get_chemical();
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			rhs[0][i][j] = lagrange-mu[i][j];
		}
	}
}

/*******************************************************************/

void run()
{
	long t, export=1;
	
	FILE *f_hi, *f_st;
	
	f_hi = fopen("hi.dat","w");	
	
	for (t=0; t<num_of_iteration; t++){
		export_conf(t,export);
		write_hi(f_hi,t);
		progress_bar(t);
		one_step();
	}

	fclose(f_hi);
}

/******************************************************************
 
 
 double get_area()
{
	double area=0;
	long i, j;
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			area += 0.5*DS*DS*(phi[i][j]+1);
		}
	}
	return area;
}

/******************************************************************
 
 void get_energy()
{
	long i, j;
	
	for (i=0; i<num_of_meshpoint; i++){
		for (j=0; j<num_of_meshpoint; j++){
			
			// To be implemented
			
			energy[i][j] = 0;
		}
	}
}

/*******************************************************************/

//Laplace operator on a spehere

double DD(double **f, long i, long j)
{	
	return DXX(f,i,j)+ 1/tan(i*dtheta)*DX(f,i,j)+1/(sin(i*dtheta)*sin(i*dtheta))*DYY(f,i,j); // since tan(0)=tan(pi)=0 the for loops have to start from i=1 and stop in i=num_of_meshpoint_theta-1
}

/*******************************************************************/

double DX(double **f, long i, long j)
{
	long prev_i,prev_j, next_i, next_j;
	
    if (i==1) {
        prev_i = 1;
        prev_j = mod(j+ num_of_meshpoint_phi, 2*num_of_meshpoint_phi);}
    
    
    else {prev_i = mod(i-1, num_of_meshpoint_theta);
          prev_j = j;}
        
        
    
    if (i==2*num_of_meshpoint_theta) {
        next_i = 2*num_of_meshpoint_theta;
        next_j = mod(j+ num_of_meshpoint_phi, 2*num_of_meshpoint_phi);}
    
    
    else {next_i = mod(i+1, num_of_meshpoint_theta);
        next_j = j;}
    
	
	return (f[next_i][next_j]-f[prev_i][prev_j])/(2*dtheta);
}

/*******************************************************************
 
 double DY(double **f, long i, long j)
{
	long prev, next;
	
	prev = mod(j-1,num_of_meshpoint);
	next = mod(j+1,num_of_meshpoint);
	
	return (f[i][next]-f[i][prev])/(2*DS);
}

/*******************************************************************/

double DXX(double **f, long i, long j)
{

    
    long prev_i,prev_j, next_i, next_j;
    
    if (i==1) {
        prev_i = 1;
        prev_j = mod(j+ num_of_meshpoint_phi, 2*num_of_meshpoint_phi);}
    
    
    else {prev_i = mod(i-1, num_of_meshpoint_theta);
        prev_j = j;}
    
    
    
    if (i==2*num_of_meshpoint_theta) {
        next_i = 2*num_of_meshpoint_theta;
        next_j = mod(j+num_of_meshpoint_phi, 2*num_of_meshpoint_phi);}
    
    
    else {next_i = mod(i+1, num_of_meshpoint_theta);
        next_j = j;}
	
	return (f[next_i][next_j]-2*f[i][j]+f[prev_i][prev_j])/(dtheta*dtheta);
}

/*******************************************************************/

double DYY(double **f, long i, long j)
{
	long prev, next;
	
	prev = mod(j-1,2*num_of_meshpoint_phi);
	next = mod(j+1,2*num_of_meshpoint_phi);
	
	return (f[i][next]-2*f[i][j]+f[i][prev])/(dphi*dphi);
}

/*******************************************************************
 
 double DXY(double **f, long i, long j)
{
	long prev_i, prev_j, next_i, next_j;
	
	prev_i = mod(i-1,num_of_meshpoint);
	prev_j = mod(j-1,num_of_meshpoint);
	next_i = mod(i+1,num_of_meshpoint);
	next_j = mod(j+1,num_of_meshpoint);
	
	return (f[next_i][next_j]-f[next_i][prev_j]
		-f[prev_i][next_j]+f[prev_i][prev_j])/(4*DS*DS);
}

/*******************************************************************/

void one_step()
{
	switch(method){
  	case 1:
		euler();
    	break;
  	case 2:
		rk2();
    	break;
  	case 3:
		rk4();
    	break;
	}
}

/*******************************************************************/

void euler()
{
	long i, j;
	double h=DT;
	
	get_rhs(k0);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] += h*k0[0][i][j];
		}
	}
}

/*******************************************************************/

void rk2()
{
	double h=DT/2;
	long i, j;
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			k0[0][i][j] = phi[i][j]; 
		}
	}
	
	get_rhs(k1);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = k0[0][i][j]+DT*k1[0][i][j];
		}
	}
	
	get_rhs(k2);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = k0[0][i][j]+h*(k1[0][i][j]+k2[0][i][j]);
		}
	}
}

/*******************************************************************/

void rk4()
{
	double h=DT/6;
	long i, j;
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			k0[0][i][j] = phi[i][j]; 
		}
	}
	
	get_rhs(k1);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = k0[0][i][j]+0.5*DT*k1[0][i][j];
		}
	}
	
	get_rhs(k2);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = k0[0][i][j]+0.5*DT*k2[0][i][j];
		}
	}
	
	get_rhs(k3);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = k0[0][i][j]+DT*k3[0][i][j];
		}
	}
	
	get_rhs(k4);
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi[i][j] = k0[0][i][j]+h*(k1[0][i][j]+2*k2[0][i][j]+2*k3[0][i][j]+k4[0][i][j]);
		}
	}
}

/*******************************************************************/

void progress_bar(long t)
{
	long i, ticks, percent, num_of_ticks=20;
	char progress[32];

	percent = 100*t/num_of_iteration;
	ticks = num_of_ticks*percent/100;
	
	for (i=0; i<20; i++){
		progress[i] = ((i<=ticks)?'#':' ');
	}
	progress[20] = '\0';
	
	printf("Progress: [%s] %ld%%\r",progress,percent);
	fflush(stdout);
}

/*******************************************************************/

void end()
{	
	long defect[2];
	CPU_Time cpu_time;
	
	FILE *f_ou;
	
	time(&t2);
  	get_time(t2-t1,&cpu_time);
	
	export_conf(-1,-1);	
		
	printf("\n\n");
	printf("\tGrid spacing in theta %g and in phi %g\n",dtheta,dphi);
	printf("\tTime step %g\n",DT);
	printf("\tSigma %g\n",sigma);
	//printf("\tArea %g\n",get_area());
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
  	printf("\n");
	
	f_ou = fopen("ou.dat","w");	
	fprintf(f_ou,"Time step %g\n",DT);
	fprintf(f_ou,"Sigma %g\n",sigma);
	//fprintf(f_ou,"Area %.10f\n",get_area());
	fprintf(f_ou,"CPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
	fclose(f_ou);

}

/*******************************************************************/

void write_hi(FILE *f_ou, long t)
{
	double phi_tot=0;
	long m, i, j;
 	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			phi_tot += phi[i][j];
		}
	}
	
	m = num_of_meshpoint_theta/2;
	
	fprintf(f_ou,"%.10f\t%.10f\t%.10f\n",t*DT,phi[m][m],phi_tot);
	
	// Legend
	
	// 1) Time
	// 2) Phi
	// 3) Integral of Phi
}

/*******************************************************************/

void export_conf(long t, long period)
{
	char f_na[32];
	double x, y;
	long i, j;

	FILE *f_ou;

	if (t%period!=0) return;
	sprintf(f_na,"s-t%ld.dat",t);
	if (t<0) sprintf(f_na,"last.dat");
	
	f_ou = fopen(f_na,"w");
	
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			x = i*dtheta;
			y = j*dphi;
			fprintf(f_ou,"%.10f\t%.10f\t%.10f\n",x,y,phi[i][j]);
		}
	}
	fclose(f_ou);
}


/*******************************************************************/

void get_time(time_t t, CPU_Time *cpu_time)
{
  long sec_in_day, sec_in_ora, sec_in_min;

  FILE *f_ti;

  /*
  f_ti = fopen("ti.dat","w");
  fprintf(f_ti,"%ld",t);
  fclose(f_ti);
  */

  sec_in_day = 86400;
  sec_in_ora = 3600;
  sec_in_min = 60;

  cpu_time->d = t/sec_in_day; t = t%sec_in_day;
  cpu_time->h = t/sec_in_ora; t = t%sec_in_ora;
  cpu_time->m = t/sec_in_min; t = t%sec_in_min;
  cpu_time->s = t;
}

/*******************************************************************/

double ran2(long *idum)
{
  	int j;
  	long k;
  	static long idum2=123456789;
  	static long iy=0;
  	static long iv[NTAB];
  	double temp;
 
  	if (*idum <= 0) { 
    	if (-(*idum) < 1) *idum=1; 
    	else *idum = -(*idum);
    	idum2 = (*idum);
    	for (j=NTAB+7; j>=0; j--) {
      		k = (*idum)/IQ1;
      		*idum = IA1*(*idum-k*IQ1)-k*IR1;
      		if (*idum < 0) *idum += IM1;
      		if (j < NTAB) iv[j] = *idum;
    	}
    	iy = iv[0];
  	}
  	k = (*idum)/IQ1; 
  	*idum = IA1*(*idum-k*IQ1)-k*IR1; 
  	if (*idum < 0) *idum += IM1; 
  	k = idum2/IQ2;
  	idum2 = IA2*(idum2-k*IQ2)-k*IR2; 
  	if (idum2 < 0) idum2 += IM2;
  	j = iy/NDIV; 
  	iy = iv[j]-idum2;
  	iv[j] = *idum; 
  	if (iy < 1) iy += IMM1;
  	if ((temp = AM*iy) > RNMX) return RNMX; 
  	else return temp;
}

/*******************************************************************/

void debug()
{
	long i, j, k, l;
	double x, y;
	
	FILE *f_ou;
	
	seed=123456789;
	
	init_random();

	i = 0;
	j = 32;
	k = 0;
	l = 10;
	
	exit(0);
	
	f_ou = fopen("debug-phi.dat","w");
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			x = i*dtheta;
			y = j*dphi;
			fprintf(f_ou,"%g %g %g\n",x,y,phi[i][j]);
		}
	}
	fclose(f_ou);
	
	get_rhs(k0);

	f_ou = fopen("debug-rhs1.dat","w");
	for (i=1; i<=num_of_meshpoint_theta-1; i++){
		for (j=1; j<2*num_of_meshpoint_phi-1; j++){
			x = i*dtheta;
			y = j*dphi;
			fprintf(f_ou,"%g %g %g\n",x,y,k0[0][i][j]);
		}
	}
	fclose(f_ou);
	
	exit(0);
}
