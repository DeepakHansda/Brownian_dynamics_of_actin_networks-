```c
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


//calculate random force
 #define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum);
float randforce(float sum,float w);
long int id=-9084099491-PARA;

//long int id= -9081093757;
// it is ks = 100 simulation
int main()
{
FILE *energy,*energy1,*fp1,*fp2,*fp4,*fp3,*fp7;
float mu=0.0002,kb=300.0,ks=200.0,b= 1.0, pi=3.14159265359,pol_x,pol_y,pol_z,E_lj_pw,rc_pw = b*b,ini_dis_x = 14.0,ini_dis_y = 0.0;
int N=2,timep=0,M = 20,bead_count,pol_len,time_step = 10000,t_max = pow(10,8),nothing,v,ip,po_i,length_pol[N],pol_num;
///long id=-1829307645;
//M = (t_max/time_step) - 5;
//long int id=-9084099497;
//b is linear distance between two wall bead

/*
N=2;
M=30;
ini_dis_x = 13.0;
ini_dis_y = 0.0;


*/



float L,hL,r2,f,r6i,dx1,dy1,dz1,dx2,dy2,dz2,r_nn2,F_LJ_pw_x[N],F_LJ_pw_y[N],F_LJ_pw_z[N];
float x[N],y[N],z[N],fx[N],fy[N],fz[N],fsX[N],fsY[N],fsZ[N],F_LJx[N],F_LJy[N],F_LJz[N],px[N],py[N],pz[N],pol_wall_lj_totx[N],pol_wall_lj_toty[N],pol_wall_lj_totz[N],y_dis_pol[N];

float  rantemp, sum, rantemp1, sum1, sum2, w=1.0,cn, dis[N], end_dis, epsilon=1.0,sigma=1.0,rc= b*b,fnnX1[N],fnnY1[N],fnnZ1[N],dx[N-1],dy[N-1],dz[N-1],rsq[N-1],rsq_root[N-1],inv_rsq_root[N-1],mmrsq[N-2],mmrsq_root[N-2],mmrsq_root2[N-2],mdx[N-2],mdy[N-2],mdz[N-2],min_disx,min_disy,min_disz,dis_pw[N],max_dis,rq_d_sq = 4.0;
int t,i,j,S2,t1;
cn=sqrt(2*mu);



char Num[100],Num1[100],Num2[100],num3[100],num4[100],Num7[100];


for(i=0;i<=(N-1);i++)
{
x[i] = 0.0;
y[i] = 0.0;
z[i] = 0.0;
fx[i] = 0.0;
fy[i] = 0.0;
fz[i] = 0.0;
fsX[i] = 0.0;
fsY[i] = 0.0;
fsZ[i] = 0.0;
F_LJx[i] = 0.0;
F_LJy[i] = 0.0;
F_LJz[i] = 0.0;
px[i] = 0.0;
py[i] = 0.0;
pz[i] = 0.0;
y_dis_pol[i] = 0.0;
length_pol[i] = 0;
}


for (i=0;i<=(N-1);i++)
{
F_LJ_pw_x[i]=0.0;
F_LJ_pw_y[i]=0.0;
F_LJ_pw_z[i]=0.0;
//F_LJ_pw_x2[i]=0.0; //define
//F_LJ_pw_y2[i]=0.0;
//F_LJ_pw_z2[i]=0.0;
pol_wall_lj_totx[i] = 0.0;
pol_wall_lj_toty[i] = 0.0;
pol_wall_lj_totz[i] = 0.0;
}



for(i=0;i<(N-1);i++)
{
dx[i] = 0.0;
dy[i] = 0.0;
dz[i] = 0.0;
rsq[i] = 0.0;
rsq_root[i] = 0.0;
inv_rsq_root[i] = 0.0;
}


for(i=0;i<(N-2);i++)
{
mmrsq[i] = 0.0;
mmrsq_root[i] = 0.0;
mmrsq_root2[i] = 0.0;
mdx[i] = 0.0;
mdy[i] = 0.0;
mdz[i] = 0.0;
}

t1=0;

for(i=0;i<=(N-1);i++)
{
x[i] = ini_dis_x - t1;
y[i] = ini_dis_y;
z[i] = 0.0;
t1 = t1 + 1.0;
}


/*
if(i==0)
{
x[i] = ini_dis_x;
y[i] = ini_dis_y;
z[i] = 0.0;
}
else
{
if(i%2 == 0)
{
x[i] = ini_dis_x;
y[i] = y[i-1] + 0.5;
z[i] = 0.0;
}
else
{
x[i] = x[i-1] + cos((3.14/180)*30);
y[i] = y[i-1] + sin((3.14/180)*30);
z[i] = 0.0;
}
}
*/



/*
x[0] = 13.0;
y[0] = -25.0;
z[0] = 0.0;
for(i=0;i<=(N-2);i++)
{
x[i+1] = x[i] - cos(pi/100*i);
y[i+1] = y[i] - sin(pi/100*i);
z[i+1] = 0.0;
}
*/
sprintf( Num7 ,"initial_config.xyz");
fp7=fopen(Num7, "w");
fprintf(fp7," %d \n",N);
fprintf(fp7,"timesptep %d \n",0);
for(i=0;i<=(N-1);i++)
{
fprintf(fp7," %d %f %f %f \n",1,x[i],y[i],z[i]);
}
fclose(fp7);


length_pol[1] = 5;
//length_pol[2] = 5; // define
y_dis_pol[1] = 0.0;
//y_dis_pol[2] =  1.0;

sprintf( Num2 ,"Bead_length");
fp1=fopen(Num2, "w");


for(t=0;t < t_max;t++) // time loop
{

	for(i=0;i<(N-1);i++)
	{
	    dx[i]=x[i]-x[i+1];
        dy[i]=y[i]-y[i+1];
        dz[i]=z[i]-z[i+1];




        rsq[i]= dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i];
        rsq_root[i]=sqrt(rsq[i]);
        inv_rsq_root[i]=1.0/rsq_root[i];
	}


/*

	for(i=0;i<(N-2);i++)

{
mmrsq[i]=(rsq[i])*(rsq[i+1]);
mmrsq_root[i]=1.0/sqrt(mmrsq[i]);
mmrsq_root2[i]=1.0/(mmrsq[i]*sqrt(mmrsq[i]));
mdx[i]=(dx[i])*(dx[i+1]);
mdy[i]=(dy[i])*(dy[i+1]);
mdz[i]=(dz[i])*(dz[i+1]);
}

*/
	// bending energy





	i=0;

/*
fx[i] = kb*(-mmrsq_root[i]*(dx[i+1])+(dx[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fy[i] = kb*(-mmrsq_root[i]*(dy[i+1])+(dy[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fz[i] = kb*(-mmrsq_root[i]*(dz[i+1])+(dz[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

*/
fsX[i]=2*ks*(-b+rsq_root[i])*(dx[i])*inv_rsq_root[i];

fsY[i]=2*ks*(-b+rsq_root[i])*(dy[i])*inv_rsq_root[i];

fsZ[i]=2*ks*(-b+rsq_root[i])*(dz[i])*inv_rsq_root[i];



/*


i=1;
fx[i]= kb*(-mmrsq_root[i-1]*(dx[i-1]-dx[i])-mmrsq_root[i]*(dx[i+1])+((dx[i]*2.0)*(rsq[i-1])-(dx[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)+(dx[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fy[i] =kb*( -mmrsq_root[i-1]*(dy[i-1]-dy[i])-mmrsq_root[i]*(dy[i+1])+((dy[i]*2.0)*(rsq[i-1])-(dy[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)+(dy[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fz[i]= kb*(-mmrsq_root[i-1]*(dz[i-1]-dz[i])-mmrsq_root[i]*(dz[i+1])+((dz[i]*2.0)*(rsq[i-1])-(dz[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)+(dz[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fsX[i]=-2*ks*(-b+rsq_root[i-1])*(dx[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dx[i])*inv_rsq_root[i];

fsY[i]=-2*ks*(-b+rsq_root[i-1])*(dy[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dy[i])*inv_rsq_root[i];

fsZ[i]=-2*ks*(-b+rsq_root[i-1])*(dz[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dz[i])*inv_rsq_root[i];


for(i=2;i<=(N-3);i++)
{
fx[i] = kb*(-mmrsq_root[i-1]*(dx[i-1]-dx[i])+mmrsq_root[i-2]*(dx[i-2])-mmrsq_root[i]*(dx[i+1])+((dx[i]*2.0)*(rsq[i-1])-(dx[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)-(dx[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5)+(dx[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));


fy[i] = kb*(-mmrsq_root[i-1]*(dy[i-1]-dy[i])+mmrsq_root[i-2]*(dy[i-2])-mmrsq_root[i]*(dy[i+1])+((dy[i]*2.0)*(rsq[i-1])-(dy[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)-(dy[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5)+(dy[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));


fz[i] = kb*(-mmrsq_root[i-1]*(dz[i-1]-dz[i])+mmrsq_root[i-2]*(dz[i-2])-mmrsq_root[i]*(dz[i+1])+((dz[i]*2.0)*(rsq[i-1])-(dz[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)-(dz[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5)+(dz[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fsX[i]=-2*ks*(-b+rsq_root[i-1])*(dx[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dx[i])*inv_rsq_root[i];

fsY[i]=-2*ks*(-b+rsq_root[i-1])*(dy[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dy[i])*inv_rsq_root[i];

fsZ[i]=-2*ks*(-b+rsq_root[i-1])*(dz[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dz[i])*inv_rsq_root[i];

}
i=(N-2);
fx[i]= kb*(-mmrsq_root[i-1]*(dx[i-1]-dx[i])+mmrsq_root[i-2]*(dx[i-2])+((dx[i]*2.0)*(rsq[i-1])-(dx[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)-(dx[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));


fy[i] = kb*(-mmrsq_root[i-1]*(dy[i-1]-dy[i])+mmrsq_root[i-2]*(dy[i-2])+((dy[i]*2.0)*(rsq[i-1])-(dy[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)-(dy[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));


fz[i]= kb*(-mmrsq_root[i-1]*(dz[i-1]-dz[i])+mmrsq_root[i-2]*(dz[i-2])+((dz[i]*2.0)*(rsq[i-1])-(dz[i-1]*2.0)*(rsq[i]))*mmrsq_root2[i-1]*(mdx[i-1]+mdy[i-1]+mdz[i-1])*(0.5)-(dz[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));

fsX[i]=-2*ks*(-b+rsq_root[i-1])*(dx[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dx[i])*inv_rsq_root[i];

fsY[i]=-2*ks*(-b+rsq_root[i-1])*(dy[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dy[i])*inv_rsq_root[i];

fsZ[i]=-2*ks*(-b+rsq_root[i-1])*(dz[i-1])*inv_rsq_root[i-1] +2*ks*(-b+rsq_root[i])*(dz[i])*inv_rsq_root[i];



*/


i=(N-1);

/*
fx[i] = kb*(mmrsq_root[i-2]*(dx[i-2])-(dx[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));


fy[i]= kb*(mmrsq_root[i-2]*(dy[i-2])-(dy[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));


fz[i]= kb*(mmrsq_root[i-2]*(dz[i-2])-(dz[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));
*/

// spring force

fsX[i]=-2*ks*(-b+rsq_root[i-1])*(dx[i-1])*inv_rsq_root[i-1];


fsY[i]=-2*ks*(-b+rsq_root[i-1])*(dy[i-1])*inv_rsq_root[i-1];

fsZ[i]=-2*ks*(-b+rsq_root[i-1])*(dz[i-1])*inv_rsq_root[i-1];


//===================Leonard Jones potential start===================

for (i=0;i<=(N-1);i++)
	{
	F_LJx[i]=0.0;
	F_LJy[i]=0.0;
	F_LJz[i]=0.0;
	}


   for (i=0;i<=(N-2);i++) {
     for (j=i+1;j<=(N-1);j++) {
		dx1=x[i]-x[j];
                dy1=y[i]-y[j];
                dz1=z[i]-z[j];



                r2=dx1*dx1+dy1*dy1+dz1*dz1;

	if (r2 < rc) {// here rc is also the "distance squared" between the wall beads
	  r6i   = pow(rc,3)/(r2*r2*r2);

	  f     = -12.0*epsilon*(r6i*r6i-r6i);
	  F_LJx[i] += dx1*f/r2;
	  F_LJx[j] -= dx1*f/r2;
	  F_LJy[i] += dy1*f/r2;
	  F_LJy[j] -= dy1*f/r2;
	  F_LJz[i] += dz1*f/r2;
	  F_LJz[j] -= dz1*f/r2;



	}
     }


}


//***************end**********Leonard******************************** Jones*************************
for (i=0;i<=(N-1);i++)
{
pol_wall_lj_totx[i]=0.0; // define variable
pol_wall_lj_toty[i]=0.0;
pol_wall_lj_totz[i]=0.0;
}



//============================================================
for(pol_num = 1; pol_num <=1;pol_num++) // define pol_num
{


for (i=0;i<=(N-1);i++)
{
px[i] = 0.0;
py[i] = 0.0;
pz[i] = 0.0;
}

for(po_i = 1; po_i <= length_pol[pol_num]; po_i++) // define length_pol[pol_num]
{
px[po_i] = po_i;
py[po_i] = y_dis_pol[pol_num]; // define
pz[po_i] = 0.0;

for (i=0;i<=(N-1);i++)
{
F_LJ_pw_x[i]=0.0;
F_LJ_pw_y[i]=0.0;
F_LJ_pw_z[i]=0.0;

}


//ip = 1;
//E_lj_pw =0.0;

for (i=0;i<=(N-1);i++)
	{
                dx1=-px[po_i]+x[i];
                dy1=-py[po_i]+y[i];
                dz1=-pz[po_i]+z[i];
                r2=dx1*dx1+dy1*dy1+dz1*dz1;

	if (r2 < rc_pw  && r2 > 0){ // "rc_pw" is the "distance squared" between the wall and the polymer beads
	  r6i   = pow(rc_pw , 3)/(r2*r2*r2);
	   //E_lj_pw    = E_lj_pw + (r6i*r6i - 2.0*r6i);
	  f     = -12.0*epsilon*(r6i*r6i-r6i);

	  F_LJ_pw_x[i]= dx1*f/r2;

	  F_LJ_pw_y[i]= dy1*f/r2;

	  F_LJ_pw_z[i]= dz1*f/r2;
	}
	pol_wall_lj_totx[i] = pol_wall_lj_totx[i] + F_LJ_pw_x[i] ;
	pol_wall_lj_toty[i] = pol_wall_lj_toty[i] + F_LJ_pw_y[i] ;
	pol_wall_lj_totz[i] = pol_wall_lj_totz[i] + F_LJ_pw_z[i] ;
	}


}

}
//==================================================






for(i=1;i<=(N-1);i++)
{

sum=0.0;
sum1=0.0;
//sum2=0.0;
	for(j=0;j<12;j++)
	{
         rantemp=ran2(&id);
         rantemp1=ran2(&id);
         //rantemp2=ran2(&id);
         sum=sum+rantemp;
         sum1=sum1+rantemp1;
         //sum2=sum2+rantemp2;

	}
x[i]=x[i]-mu*(fx[i]+fsX[i]+F_LJx[i]+ pol_wall_lj_totx[i])+cn*(sum-6.0);
y[i]=y[i];//-mu*(fy[i]+fsY[i]+F_LJy[i]+pol_wall_lj_toty[i])+cn*(sum1-6.0);
z[i]=z[i] ;//-mu*(fz[i]+fsZ[i]+F_LJz[i]+F_LJ_pw_z[i])+cn*randforce(sum2,w);
}


for(pol_num =1; pol_num <= 1;pol_num++)
{

max_dis = 0.0;
for(i=0;i<=(N-1);i++)
{
dis_pw[i] = 0.0;
}

for(i=0;i<=(N-1);i++)
{
 min_disx = x[i] - length_pol[pol_num];
 min_disy = y[i] - y_dis_pol[pol_num];
 min_disz = z[i] - 0.0;
dis_pw[i] = min_disx*min_disx + min_disy*min_disy + min_disz*min_disz ;


if(dis_pw[i] > max_dis)
{
max_dis = dis_pw[i];
}
}

for(i=0;i<=(N-1);i++)
{
if(dis_pw[i] < max_dis)
{
max_dis = dis_pw[i];
}
}
if(sqrt(max_dis) >= sqrt(rq_d_sq))
{
length_pol[pol_num] = length_pol[pol_num] + 1;
}


}

//end_dis=sqrt((x[0] - x[N-1])*(x[0] - x[N-1]) +(y[0] - y[N-1])*(y[0] - y[N-1])+(z[0] - z[N-1])*(z[0] - z[N-1]));
//if(length_pol[1]>=length_pol[2])
pol_len = length_pol[1];
//else
//pol_len = length_pol[2];

if( t >= timep )
{
fprintf(fp1, "%d %d\n",timep,pol_len );

sprintf( Num ,"Wall_coordinate%d.xyz",timep);
energy=fopen(Num, "w");
fprintf(energy, "%d\n",(N+M));
fprintf(energy, "timestep %d\n",timep );
bead_count = 0;
for(i=0;i<=(N-1);i++)
{

	fprintf(energy,"%d  %f %f %f \n",1,x[i],y[i],z[i]);
bead_count = bead_count + 1;
}
/*
if(pol_len >= 11)
{
pol_len = 11;
}
else
{
nothing = 0;
}
*/
for(pol_num = 1;pol_num<=1;pol_num++)
{
for(i=1;i<= length_pol[pol_num];i++)
{
pol_x = i;
pol_y = y_dis_pol[pol_num];
pol_z = 0.0;
fprintf(energy,"%d  %f %f %f \n",2,pol_x,pol_y,pol_z);
bead_count = bead_count + 1;
}

}
do
{
pol_x = 0.0;
pol_y = 0.0;
pol_z = 0.0;

fprintf(energy,"%d  %f %f %f \n",2,pol_x,pol_y,pol_z);
bead_count = bead_count + 1;
}while(bead_count < (N+M));


timep = timep + time_step;
/*
if(timep >= timex*v)
{
pol_len = pol_len + 1;
v = v + 1;
}
*/
fclose(energy);
}

} // time loop
fclose(fp1); //  file bead_count
//fclose(fp3); // file closing polymer and wall
//fclose(fp4); // file closing bead and bead

return 0;

}





float ran2(long *idum)

{

	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;

	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;

		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;

			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;

	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;


}

#undef IM1
#undef IM2
#undef AM

#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1

#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* (C) Copr. 1986-92 Numerical Recipes Software -0)+'. */

```
