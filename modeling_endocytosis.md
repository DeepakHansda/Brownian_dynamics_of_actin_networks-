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
long int id=-9084097497;

//long int id= -9081093757;
// it is ks = 100 simulation -PARA
main()
{
FILE *energy,*energy1,*fp1,*fp2,*fp4,*fp3,*fp7,*fp10,*fp11;
float mu=0.0002,mu2=0.000025, pi=3.14159265359,pol_x,pol_y,pol_z,E_lj_pw,ini_dis_x = 13.0,ini_dis_y = -25.0;
int mlt=3,N=161,timep=0,M =200,bead_count,pol_len,time_step = 10000,t_max = pow(10,8),nothing,v,ip,po_i,pol_num,ji1,ji2,num_wall_bead = 20,t_flag,lower_len,upper_len,tm_store,additional_bead,middle_dis;

mu = 0.000025;
float L,hL,r2,f,r6i,dx1,dy1,dz1,dx2,dy2,dz2;
float x[mlt*N],y[mlt*N],z[mlt*N],fx[mlt*N],fy[mlt*N],fz[mlt*N],fsX[mlt*N],fsY[mlt*N],fsZ[mlt*N],F_LJx[mlt*N],F_LJy[mlt*N],F_LJz[mlt*N];

float  rantemp, sum, rantemp1, sum1, sum2, w=1.0,cn,cn2, dis[mlt*N], end_dis, epsilon=1.0,sigma=1.0,dx[mlt*N-1],dy[mlt*N-1],dz[mlt*N-1],rsq[mlt*N-1],rsq_root[mlt*N-1],inv_rsq_root[mlt*N-1],mmrsq[mlt*N-2],mmrsq_root[mlt*N-2],mmrsq_root2[mlt*N-2],mdx[mlt*N-2],mdy[mlt*N-2],mdz[mlt*N-2],min_disx,min_disy,min_disz,dis_pw[mlt*N],max_dis, bend_start_low_y = 20.0;
int t,i,j,S2,t1;
cn=sqrt(2*mu);
cn2 = sqrt(2*mu2);
// pol_wall_lj_toty[N],pol_wall_lj_totx[N],pol_wall_lj_totz[N], F_LJ_pw_x[N],F_LJ_pw_y[N],F_LJ_pw_z[N], y_front_bead[N]
//=== start filament section definition

int num_fil =7,sum_num_fil,br_bead_one,br_bead_two,beyond_lower_len; // "num_fil" number will be 1 greater than the number of filament in the system

float fl_x[num_fil][mlt*M], fl_y[num_fil][mlt*M], fl_z[num_fil][mlt*M], fl_fx[num_fil][mlt*M], fl_fy[num_fil][mlt*M], fl_fz[num_fil][mlt*M], fl_fsX[num_fil][mlt*M],
fl_fsY[num_fil][mlt*M], fl_fsZ[num_fil][mlt*M], fl_F_LJx[num_fil][mlt*M], fl_F_LJy[num_fil][mlt*M], fl_F_LJz[num_fil][mlt*M], fl_dx[num_fil][mlt*M], fl_dy[num_fil][mlt*M],
fl_dz[num_fil][mlt*M], fl_rsq[num_fil][mlt*M], fl_rsq_root[num_fil][mlt*M],fl_inv_rsq_root[num_fil][mlt*M],fl_mmrsq[num_fil][mlt*M],fl_mmrsq_root[num_fil][mlt*M],
fl_mmrsq_root2[num_fil][mlt*M], fl_mdx[num_fil][mlt*M], fl_mdy[num_fil][mlt*M], fl_mdz[num_fil][mlt*M];

float pol_and_fl_dx,pol_and_fl_dy,pol_and_fl_dz,pol_and_fl_rsq,pol_and_fl_rsq_root,
pol_and_fl_inv_rsq_root,pol_and_fl_fsX,pol_and_fl_fsY,pol_and_fl_fsZ,
Sf_on_fl_bead_due_to_polX[num_fil][mlt*M], Sf_on_fl_bead_due_to_polY[num_fil][mlt*M],Sf_on_fl_bead_due_to_polZ[num_fil][mlt*M],


b=1.0, rc = b*b, fl_b = 2.0, fl_rc = fl_b*fl_b, rc_pw = (0.5*b + 0.5*fl_b)*(0.5*b + 0.5*fl_b), ks=200.0,
fl_ks = 200.0, kb=60.0, fl_kb = 1000.0, pol_and_fl_b = fl_b, pol_and_fl_kb, pol_and_fl_ks = 200.0,
rq_d_sq = (0.5*b + fl_b + 0.5*fl_b)*(0.5*b + fl_b + 0.5*fl_b), fl_fl_rsq = fl_b*fl_b, poly_and_fl_b = 0.5*b + 0.5*fl_b,
poly_and_fl_ks = 70.0,

pol_wall_lj_totx[num_fil][mlt*M], pol_wall_lj_toty[num_fil][mlt*M], pol_wall_lj_totz[num_fil][mlt*M],
fl_pol_wall_lj_totx[num_fil][mlt*M],fl_pol_wall_lj_toty[num_fil][mlt*M],fl_pol_wall_lj_totz[num_fil][mlt*M],
F_LJ_pw_x[num_fil][mlt*M],F_LJ_pw_y[num_fil][mlt*M],F_LJ_pw_z[num_fil][mlt*M],
fil_F_LJ_pw_x[num_fil][mlt*M],fil_F_LJ_pw_y[num_fil][mlt*M],fil_F_LJ_pw_z[num_fil][mlt*M],

y_front_bead[mlt*N],x_front_bead[mlt*N],bead_gap,ji3,sum_x[mlt*N],sum_y[mlt*N],px[mlt*N],py[mlt*N],pz[mlt*N],y_dis_pol[mlt*N],x_dis_pol[mlt*N],

min_dis,length_pol[mlt*N], bead1_y_limit ,bead2_y_limit,


	fl_fl_lj_totx[num_fil][mlt*M],
	fl_fl_lj_toty[num_fil][mlt*M],
	fl_fl_lj_totz[num_fil][mlt*M],

fl_fl_LJ_pw_x[mlt*N],
fl_fl_LJ_pw_y[mlt*N],
fl_fl_LJ_pw_z[mlt*N],

Sf_on_pol_bead_due_to_flX[mlt*N],
Sf_on_pol_bead_due_to_flY[mlt*N],
Sf_on_pol_bead_due_to_flZ[mlt*N],
fl_fl_prevent_collaps_ks = 1000.0,
fl_fl_prevent_collaps_fl_b,
collaps_dis,fl_fl_prevent_collaps_fl_b_other_side;

int pol_N,j_pol,start_bead_num,j_temp,fl_bead_num,increased_bead,pol_N_of_fl[num_fil],shortest_dis_bead,ini_fil,pol_bd_num;

//pol_and_fl_b = ?decide here what should be the resting distance between that "filament bead" and "polymer bead"
//fl_ks =? spring constant between two filament beads
//fl_b = ? resting distance between the two filament beads
// fl_rc = ? square of the resting distance between filament beads



//=== end filament section definition
// what is "epsilon" "sigma" "w" and "b"
/*
variables
rq_d_sq = 3.0625
b= 0.5
rc_pw = 0.5625
mu=0.0002
rc = 0.25
*/

char Num[100],Num1[100],Num2[100],num3[100],num4[100],Num7[100],Num11[100];

for(i=0;i< num_fil;i++)
{
pol_N_of_fl[i] = 0;

}

for(i=0;i<=(mlt*N-1);i++)
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
Sf_on_pol_bead_due_to_flX[i] = 0.0;
Sf_on_pol_bead_due_to_flY[i] = 0.0;
Sf_on_pol_bead_due_to_flZ[i] = 0.0;
//=== below start filament variable

y_front_bead[i] = 0.0;
x_front_bead[i] = 0.0;
length_pol[i] = 0;
y_dis_pol[i] = 0.0;
x_dis_pol[i] = 0.0;
sum_x[i] = 0.0;
sum_y[i] = 0.0;
fl_fl_LJ_pw_x[i] = 0.0;
fl_fl_LJ_pw_y[i] = 0.0;
fl_fl_LJ_pw_z[i] = 0.0;

//=== end filament variable
}


for(i=0;i<=(mlt*N-1);i++)
{
px[i] = 0.0;
py[i] = 0.0;
pz[i] = 0.0;
}



//=== start filaments section: initializing the forces arrays to 0.0
for(j=1;j< num_fil;j++)
{
 	for(i=0;i<=(mlt*M-1);i++)
	{
fl_x[j][i] = 0.0;
fl_y[j][i] = 0.0;
fl_z[j][i] = 0.0;
fl_fx[j][i] = 0.0;
fl_fy[j][i] = 0.0;
fl_fz[j][i] = 0.0;
fl_fsX[j][i] = 0.0;
fl_fsY[j][i] = 0.0;
fl_fsZ[j][i] = 0.0;
fl_F_LJx[j][i] = 0.0;
fl_F_LJy[j][i] = 0.0;
fl_F_LJz[j][i] = 0.0;
Sf_on_fl_bead_due_to_polX[j][i] = 0.0;
Sf_on_fl_bead_due_to_polY[j][i] = 0.0;
Sf_on_fl_bead_due_to_polZ[j][i] = 0.0;
pol_wall_lj_totx[j][i] = 0.0;
pol_wall_lj_toty[j][i] = 0.0;
pol_wall_lj_totz[j][i] = 0.0;
fl_pol_wall_lj_totx[j][i] = 0.0;
fl_pol_wall_lj_toty[j][i] = 0.0;
fl_pol_wall_lj_totz[j][i] = 0.0;
F_LJ_pw_x[j][i] = 0.0;
F_LJ_pw_y[j][i] = 0.0;
F_LJ_pw_z[j][i] = 0.0;
fil_F_LJ_pw_x[j][i] = 0.0;
fil_F_LJ_pw_y[j][i] = 0.0;
fil_F_LJ_pw_z[j][i] = 0.0;

fl_fl_lj_totx[j][i] = 0.0;
fl_fl_lj_toty[j][i] = 0.0;
fl_fl_lj_totz[j][i] = 0.0;
	}
}
//===end filament section

for(i=0;i<(mlt*N-1);i++)
{
dx[i] = 0.0;
dy[i] = 0.0;
dz[i] = 0.0;
rsq[i] = 0.0;
rsq_root[i] = 0.0;
inv_rsq_root[i] = 0.0;
}

// start of filament section
for(j=1;j< num_fil;j++)
{

	for(i=0;i<(mlt*M-1);i++)
	{
	fl_dx[j][i] = 0.0;
	fl_dy[j][i] = 0.0;
	fl_dz[j][i] = 0.0;
	fl_rsq[j][i] = 0.0;
	fl_rsq_root[j][i] = 0.0;
	fl_inv_rsq_root[j][i] = 0.0;
	}
}
// end of filament section

for(i=0;i<(mlt*N-2);i++)
{
mmrsq[i] = 0.0;
mmrsq_root[i] = 0.0;
mmrsq_root2[i] = 0.0;
mdx[i] = 0.0;
mdy[i] = 0.0;
mdz[i] = 0.0;
}

//=== start of filament section

for(j=1;j< num_fil;j++)
{

	for(i=0;i<(mlt*M-2);i++)
	{
	fl_mmrsq[j][i] = 0.0;
	fl_mmrsq_root[j][i] = 0.0;
	fl_mmrsq_root2[j][i] = 0.0;
	fl_mdx[j][i] = 0.0;
	fl_mdy[j][i] = 0.0;
	fl_mdz[j][i] = 0.0;

	}
}
//=== end of filament section


//=== start of section to put initial configuration for polymer (membrane)
additional_bead = 55;
lower_len = 20;
middle_dis = 10;
upper_len = lower_len + middle_dis;
bead1_y_limit = ini_dis_y + (additional_bead - 15)*1.0 ,bead2_y_limit = ini_dis_y + (additional_bead + middle_dis + 15)*1.0,
// for i =0
x[0] = ini_dis_x;
y[0] = ini_dis_y;
z[0] = 0.0;

for(i=1;i <= additional_bead; i++)
// "bend_start_low_y" is  the "y" co-ordinate value till which the wall goes, starting from the least value of "y" co-ordinate
{
x[i] = ini_dis_x;
y[i] = ini_dis_y + i*1.0;
z[i] = 0.0;
}



j_temp=1;
for(i= (additional_bead + 1); i<=(additional_bead + lower_len); i++)
{
x[i] =  x[additional_bead] + j_temp*1.0 ;
y[i] =  y[additional_bead];
z[i] = 0.0;
j_temp = j_temp + 1;
}

ji1 = 1.0;
for(i=( additional_bead + lower_len + 1); i<=( additional_bead + lower_len + middle_dis); i++)
{
x[i] = x[additional_bead + lower_len];
y[i] = y[additional_bead + lower_len] + ji1;
z[i] = 0.0;
ji1 = ji1 +1.0;
}


ji2 = 1.0;
for(i= (additional_bead + lower_len + middle_dis + 1); i<=(additional_bead + lower_len+ middle_dis + lower_len); i++)
{
x[i] = x[additional_bead+lower_len+middle_dis] - ji2;
y[i] = y[additional_bead+lower_len+middle_dis];
z[i] = 0.0;
ji2 = ji2 + 1.0;
}



ji3 = 1.0;
for(i=(additional_bead+lower_len+middle_dis+lower_len + 1); i<N; i++ )
{
x[i] = x[additional_bead+lower_len+middle_dis+lower_len];
y[i] = y[additional_bead+lower_len+middle_dis+lower_len] + ji3;
z[i] = 0.0;
ji3 = ji3 + 1.0;
}
//=== end of section to put initial configuration for polymer (membrane)

//=== start of section to put initial configuration of filaments

beyond_lower_len = 0;

x_dis_pol[1] = ini_dis_x + lower_len + beyond_lower_len;//x-coordinate of the "first" bead, but here y-coordinate of all the beads in the filament will ramain as stated in the right hand side of the equation
y_dis_pol[1] = bead1_y_limit;
x_dis_pol[4] = ini_dis_x + lower_len + beyond_lower_len;//x-coordinate of the "first" bead, but here y-coordinate of all the beads in the filament will ramain as stated in the right hand side of the equation
y_dis_pol[4] = bead2_y_limit;


pol_N_of_fl[1] = (lower_len + beyond_lower_len)/2; // check every time this value and this value should not exceed the array size i.e 100

pol_N_of_fl[2] = 3;
pol_N_of_fl[3] = 3;
pol_N_of_fl[4] = (lower_len + beyond_lower_len)/2; // check every time this value and this value should not exceed the array size i.e 100
pol_N_of_fl[5] = 3;
pol_N_of_fl[6] = 3;

 

br_bead_one = 3 + beyond_lower_len/2; // always make sure that the these "br_bead_one and br_bead_two" number are within br_bead_one
br_bead_two = 6 + beyond_lower_len/2;

br_bead_one_given = pol_N_of_fl[1];
br_bead_one_new_value = pol_N_of_fl[1];

br_bead_four_given = pol_N_of_fl[4] ;
br_bead_four_new_value = pol_N_of_fl[4];

j_pol = 1;
start_bead_num = pol_N_of_fl[j_pol];
bead_gap = 0.0;
		for(i=0; i <=(start_bead_num - 1); i++)
		{
		fl_x[j_pol][i] = x_dis_pol[j_pol] - bead_gap;
		fl_y[j_pol][i] = y_dis_pol[j_pol];
		fl_z[j_pol][i] = 0.0;
		bead_gap = bead_gap + 2.0;
		}


j_pol = 2;
start_bead_num = pol_N_of_fl[j_pol];
for(i=0;i<=(start_bead_num - 1);i++)
{
                if(i==0)
		{
		
		fl_x[j_pol][i] = fl_x[1][br_bead_one]; // do some thing for the array fl_x[i][j]; so that it becomes automatic
		fl_y[j_pol][i] = fl_y[1][br_bead_one] + fl_b;
		fl_z[j_pol][i] = 0.0;
		}
		else
		{
		fl_x[j_pol][i] = fl_x[j_pol][i-1] - (2.0*cos((70.0*3.14)/180.0));
		fl_y[j_pol][i] = fl_y[j_pol][i-1] + (2.0*sin((70.0*3.14)/180.0));
		fl_z[j_pol][i] = 0.0;

		}
}


j_pol = 3;
start_bead_num = pol_N_of_fl[j_pol];

for(i=0; i<=(start_bead_num - 1) ;i++)
{

		if(i==0)
		{
		fl_x[j_pol][i] = fl_x[1][br_bead_two]; // do some thing for the array fl_x[i][j]; so that it becomes automatic
		fl_y[j_pol][i] = fl_y[1][br_bead_two] + fl_b;
                fl_z[j_pol][i] = 0.0;
		}
		else
		{
		fl_x[j_pol][i] = fl_x[j_pol][i-1] - (2.0*cos((70.0*3.14)/180.0));
		fl_y[j_pol][i] = fl_y[j_pol][i-1] + (2.0*sin((70.0*3.14)/180.0));
		fl_z[j_pol][i] = 0.0;

		}

}


j_pol = 4;
start_bead_num = pol_N_of_fl[j_pol] ;
bead_gap = 0.0;
		for(i=0; i<=(start_bead_num - 1);i++)
		{
		fl_x[j_pol][i] = x_dis_pol[j_pol] - bead_gap;
		fl_y[j_pol][i] = y_dis_pol[j_pol];
		fl_z[j_pol][i] = 0.0;
		bead_gap = bead_gap + 2.0;
		}


j_pol = 5;
start_bead_num = pol_N_of_fl[j_pol];

for(i=0; i<=(start_bead_num - 1);i++)
{
		if(i==0)
		{
		fl_x[j_pol][i] = fl_x[4][br_bead_one];// do some thing for the array fl_x[i][j]; so that it becomes automatic
		fl_y[j_pol][i] = fl_y[4][br_bead_one] - fl_b;
		fl_z[j_pol][i] = 0.0;
		}
		else
		{
		fl_x[j_pol][i] = fl_x[j_pol][i-1] - (2.0*cos((70.0*3.14)/180.0));
		fl_y[j_pol][i] = fl_y[j_pol][i-1] - (2.0*sin((70.0*3.14)/180.0));
		fl_z[j_pol][i] = 0.0;

		}
}


j_pol = 6;
start_bead_num = pol_N_of_fl[j_pol];

for(i=0; i<=(start_bead_num - 1);i++)
{
		if(i==0)
		{
		fl_x[j_pol][i] = fl_x[4][br_bead_two];// do some thing for the array fl_x[i][j]; so that it becomes automatic
		fl_y[j_pol][i] = fl_y[4][br_bead_two] - fl_b;
		fl_z[j_pol][i] = 0.0;
		}
		else
		{
		fl_x[j_pol][i] = fl_x[j_pol][i-1] - (2.0*cos((70.0*3.14)/180.0));
		fl_y[j_pol][i] = fl_y[j_pol][i-1] - (2.0*sin((70.0*3.14)/180.0));
		fl_z[j_pol][i] = 0.0;

		}
}



collaps_dis = (fl_x[1][br_bead_one+2] - fl_x[2][1])*(fl_x[1][br_bead_one+2] - fl_x[2][1]) + (fl_y[1][br_bead_one+2] - fl_y[2][1])*(fl_y[1][br_bead_one+2] - fl_y[2][1]);
fl_fl_prevent_collaps_fl_b = sqrt(collaps_dis);



collaps_dis = (fl_x[1][br_bead_one-2] - fl_x[2][1])*(fl_x[1][br_bead_one-2] - fl_x[2][1]) + (fl_y[1][br_bead_one-2] - fl_y[2][1])*(fl_y[1][br_bead_one-2] - fl_y[2][1]);
fl_fl_prevent_collaps_fl_b_other_side = sqrt(collaps_dis);

//=== end of section to put initial configuration of filaments

/*
for(i=0;i<=(N-1);i++)
{
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

}
*/

/*
x[0] = 13.0;
y[0] = -25.0;
z[0] = 0.0;
for(i=0;i<=(N-2);i++)
{
x[i+1] = x[i] - cos((pi/100)*i);
y[i+1] = y[i] - sin((pi/100)*i);
z[i+1] = 0.0;
}
*/

sprintf( Num7 ,"initial_config.xyz");
fp7=fopen(Num7, "w");
sum_num_fil = 0;
for(i=1; i< num_fil; i++)
{
sum_num_fil = sum_num_fil + pol_N_of_fl[i];
}

fprintf(fp7," %d \n",N + sum_num_fil);
fprintf(fp7,"timesptep %d \n",0);
for(i=0;i<=(N-1);i++)
{
fprintf(fp7," %d %f %f %f \n",1,x[i],y[i],z[i]);
}

for(j_pol = 1;j_pol < num_fil; j_pol++)
{
start_bead_num = pol_N_of_fl[j_pol];

	for(i=0;i<=(start_bead_num-1);i++)
	{
	fprintf(fp7," %d %f %f %f \n",2,fl_x[j_pol][i],fl_y[j_pol][i],fl_z[j_pol][i]);
	}
}
fclose(fp7);





//v = 1;
sprintf( Num2 ,"Bead_length");
fp1=fopen(Num2, "w");
sprintf(Num11, "tau");
fp11 = fopen(Num11, "w");
t_flag = 0;
//sprintf(num3,"pw_force"); // for printing polymer and wall force
//fp3 = fopen(num3,"w");

//sprintf(num4,"bb_force"); // for printing bead and bead force
//fp4 = fopen(num4,"w");
for(t=0;t < t_max;t++)
{


/*
if(br_bead_len_given[1] > br_bead_len_new_value[1])
{

br_bead[1][br_bead_one] = br_bead[1][br_bead_one] + 1;

br_bead_len_given[1] = br_bead_len_new_value[1];
}

if(br_bead_four_given > br_bead_four_new_value)
{
br_bead_four = br_bead_four + 1;

br_bead_four_given = br_bead_four_new_value;
}

/==========================================
br_bead_one have to change this value
if(br_bead_one_given > br_bead_one_new_value)
{
br_bead_one = br_bead_one + 1;

br_bead_one_given = br_bead_one_new_value;
}

br_bead_two have to change this value
if(br_bead_two_given > br_bead_two_new_value)
{
br_bead_two = br_bead_two + 1;
br_bead_two_given = br_bead_two_new_value;

}
*/
	for(i=0;i<(N-1);i++)
	{
	dx[i]=x[i]-x[i+1];
        dy[i]=y[i]-y[i+1];
        dz[i]=z[i]-z[i+1];




	rsq[i]= dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i];
        rsq_root[i]=sqrt(rsq[i]);
	inv_rsq_root[i]=1.0/rsq_root[i];


	}


//=== start of filament section
for(j=1;j< num_fil;j++)
{
pol_N = pol_N_of_fl[j];
for(i=0;i<(pol_N - 1);i++)
	{
	fl_dx[j][i]=fl_x[j][i]-fl_x[j][i+1];
        fl_dy[j][i]=fl_y[j][i]-fl_y[j][i+1];
        fl_dz[j][i]=fl_z[j][i]-fl_z[j][i+1];




	fl_rsq[j][i]= fl_dx[j][i]*fl_dx[j][i] + fl_dy[j][i]*fl_dy[j][i] + fl_dz[j][i]*fl_dz[j][i];
        fl_rsq_root[j][i]=sqrt(fl_rsq[j][i]);
	fl_inv_rsq_root[j][i]=1.0/fl_rsq_root[j][i];
	}
}
//=== end of filament section

	for(i=0;i<(N-2);i++)

{
mmrsq[i]=(rsq[i])*(rsq[i+1]);
mmrsq_root[i]=1.0/sqrt(mmrsq[i]);
mmrsq_root2[i]=1.0/(mmrsq[i]*sqrt(mmrsq[i]));
mdx[i]=(dx[i])*(dx[i+1]);
mdy[i]=(dy[i])*(dy[i+1]);
mdz[i]=(dz[i])*(dz[i+1]);
}

//=== start of filament section
for(j=1;j< num_fil;j++)
{
pol_N = pol_N_of_fl[j];
	for(i=0;i<(pol_N - 2);i++)
	{
	fl_mmrsq[j][i]=(fl_rsq[j][i])*(fl_rsq[j][i+1]);
	fl_mmrsq_root[j][i]=1.0/sqrt(fl_mmrsq[j][i]);
	fl_mmrsq_root2[j][i]=1.0/(fl_mmrsq[j][i]*sqrt(fl_mmrsq[j][i]));
	fl_mdx[j][i]=(fl_dx[j][i])*(fl_dx[j][i+1]);
	fl_mdy[j][i]=(fl_dy[j][i])*(fl_dy[j][i+1]);
	fl_mdz[j][i]=(fl_dz[j][i])*(fl_dz[j][i+1]);
	}
}
//=== end of filament section


	// bending energy


//=====start of (spring+bending) potential of membrane
i=0;

fx[i] = kb*(-mmrsq_root[i]*(dx[i+1])+(dx[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fy[i] = kb*(-mmrsq_root[i]*(dy[i+1])+(dy[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fz[i] = kb*(-mmrsq_root[i]*(dz[i+1])+(dz[i]*2.0)*mmrsq_root2[i]*(rsq[i+1])*(mdx[i]+mdy[i]+mdz[i])*(0.5));

fsX[i]=2*ks*(-b+rsq_root[i])*(dx[i])*inv_rsq_root[i];

fsY[i]=2*ks*(-b+rsq_root[i])*(dy[i])*inv_rsq_root[i];

fsZ[i]=2*ks*(-b+rsq_root[i])*(dz[i])*inv_rsq_root[i];



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




i=(N-1);
fx[i] = kb*(mmrsq_root[i-2]*(dx[i-2])-(dx[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));


fy[i]= kb*(mmrsq_root[i-2]*(dy[i-2])-(dy[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));


fz[i]= kb*(mmrsq_root[i-2]*(dz[i-2])-(dz[i-1]*2.0)*mmrsq_root2[i-2]*(rsq[i-2])*(mdx[i-2]+mdy[i-2]+mdz[i-2])*(0.5));

// spring force
fsX[i]=-2*ks*(-b+rsq_root[i-1])*(dx[i-1])*inv_rsq_root[i-1];
fsY[i]=-2*ks*(-b+rsq_root[i-1])*(dy[i-1])*inv_rsq_root[i-1];
fsZ[i]=-2*ks*(-b+rsq_root[i-1])*(dz[i-1])*inv_rsq_root[i-1];
//=====end  of (spring+bending) potential of membrane



//=== start of calculation for the (spring+bending) potential "filament"
for(j=1;j< num_fil;j++)
{
pol_N = pol_N_of_fl[j];

i=0;

if(j==1 && i==0)
{
 pol_bd_num = additional_bead + lower_len; // it shows the bead number on polymer to which the first bead of this particular "j" filament would interact


        pol_and_fl_dx = fl_x[j][i]-x[pol_bd_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-y[pol_bd_num];
        pol_and_fl_dz = fl_z[j][i]-z[pol_bd_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;

//===============================================



pol_and_fl_fsX = 2*poly_and_fl_ks*(-poly_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*poly_and_fl_ks*(-poly_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*poly_and_fl_ks*(-poly_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_pol_bead_due_to_flX[pol_bd_num] = -pol_and_fl_fsX;
Sf_on_pol_bead_due_to_flY[pol_bd_num] = -pol_and_fl_fsY;
Sf_on_pol_bead_due_to_flZ[pol_bd_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dx[j][i+1])+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fy[j][i] = fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dy[j][i+1])+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fz[j][i]=  fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dz[j][i+1])+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=2*fl_ks*(-fl_b+fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i]+pol_and_fl_fsX;

fl_fsY[j][i]=2*fl_ks*(-fl_b+fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i]+pol_and_fl_fsY;

fl_fsZ[j][i]=2*fl_ks*(-fl_b+fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i]+pol_and_fl_fsZ;
}
else if(j==4 && i==0)
{
pol_bd_num = (additional_bead + lower_len + middle_dis) ; // it shows the bead number on polymer to which the first bead of this particular "j" filament would interact
		// do some thing to automate it

	pol_and_fl_dx = fl_x[j][i]-x[pol_bd_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-y[pol_bd_num];
        pol_and_fl_dz = fl_z[j][i]-z[pol_bd_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;

//===============================================



pol_and_fl_fsX = 2*poly_and_fl_ks*(-poly_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*poly_and_fl_ks*(-poly_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*poly_and_fl_ks*(-poly_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_pol_bead_due_to_flX[pol_bd_num] = -pol_and_fl_fsX;
Sf_on_pol_bead_due_to_flY[pol_bd_num] = -pol_and_fl_fsY;
Sf_on_pol_bead_due_to_flZ[pol_bd_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dx[j][i+1])+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fy[j][i] = fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dy[j][i+1]) + (fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fz[j][i]=fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dz[j][i+1])+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX;

fl_fsY[j][i]=2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}
else
{
fl_fx[j][i] = fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dx[j][i+1])+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fy[j][i] = fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dy[j][i+1]) + (fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fz[j][i]=fl_kb*(-(fl_mmrsq_root[j][i])*(fl_dz[j][i+1])+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=2*fl_ks*(-fl_b+fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i];

fl_fsY[j][i]=2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i];

fl_fsZ[j][i]=2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i];
}


i=1;

fl_fx[j][i]= fl_kb*(-(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1] - fl_dx[j][i])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fy[j][i] = fl_kb*(-(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));

fl_fz[j][i] = fl_kb*(-(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] +2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i];

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i];

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i];


for(i=2;i<=(pol_N - 3);i++)
{
/*
//======part for the two beads from two different filaments; see the diagram (spring between # and *)
-------#------
        *
	 *
	  *
*/

if(j==1 && i== br_bead_one) // comparing the filament number, so that we can calculate the spring potential with branched filament
{
// "j" shows the mother filament
j_pol = 2; // it shows the branched filament
fl_bead_num = 0; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;

//===============================================



pol_and_fl_fsX = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;

fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;

}

else if(j==1 && i==(br_bead_one-2))
{
// "j" shows the mother filament
j_pol = 2; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}

else if(j==1 && i==(br_bead_two-2))
{
// "j" shows the mother filament
j_pol = 3; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}


else if(j==4 && i==(br_bead_one-2))
{
// "j" shows the mother filament
j_pol = 5; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}
else if(j==4 && i==(br_bead_two-2))
{
// "j" shows the mother filament
j_pol = 6; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b_other_side + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}

else if(j==1 && i== br_bead_two)
{

// "j" shows the mother filament
j_pol = 3; // it shows the branched filament
fl_bead_num = 0; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;

//===============================================



pol_and_fl_fsX = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;

fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;



}
else if(j==4 && i== br_bead_one)
{

// "j" shows the mother filament
j_pol = 5; // it shows the branched filament
fl_bead_num = 0; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;

//===============================================



pol_and_fl_fsX = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;

fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;





}
else if(j==4 && i== br_bead_two)
{

// "j" shows the mother filament
j_pol = 6; // it shows the branched filament
fl_bead_num = 0; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;

//===============================================



pol_and_fl_fsX = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*pol_and_fl_ks*(-pol_and_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;

fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}

else if(j==1 && i== br_bead_one+2 )
{
// "j" shows the mother filament
j_pol = 2; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;

}

else if(j==1 && i== br_bead_two+2 )
{
// "j" shows the mother filament
j_pol = 3; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}


else if(j==4 && i== br_bead_one+2)
{
// "j" shows the mother filament
j_pol = 5; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}

else if(j==4 && i== br_bead_two+2)
{
// "j" shows the mother filament
j_pol = 6; // it shows the branched filament
fl_bead_num = 1; // it shows the first bead of the branched filament


	pol_and_fl_dx = fl_x[j][i]-fl_x[j_pol][fl_bead_num]; // define the variables
        pol_and_fl_dy = fl_y[j][i]-fl_y[j_pol][fl_bead_num];
        pol_and_fl_dz = fl_z[j][i]-fl_z[j_pol][fl_bead_num];




pol_and_fl_rsq = pol_and_fl_dx*pol_and_fl_dx + pol_and_fl_dy*pol_and_fl_dy + pol_and_fl_dz*pol_and_fl_dz;
pol_and_fl_rsq_root = sqrt(pol_and_fl_rsq);
pol_and_fl_inv_rsq_root = 1.0/pol_and_fl_rsq_root;





pol_and_fl_fsX = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dx)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsY = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dy)*pol_and_fl_inv_rsq_root;
pol_and_fl_fsZ = 2*fl_fl_prevent_collaps_ks*(-fl_fl_prevent_collaps_fl_b + pol_and_fl_rsq_root)*(pol_and_fl_dz)*pol_and_fl_inv_rsq_root;

Sf_on_fl_bead_due_to_polX[j_pol][fl_bead_num] = -pol_and_fl_fsX;
Sf_on_fl_bead_due_to_polY[j_pol][fl_bead_num] = -pol_and_fl_fsY;
Sf_on_fl_bead_due_to_polZ[j_pol][fl_bead_num] = -pol_and_fl_fsZ;


fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsX ;

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsY;

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i] + pol_and_fl_fsZ;
}

else
{
fl_fx[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])-fl_mmrsq_root[j][i]*(fl_dx[j][i+1])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dx[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fy[j][i]=fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])-fl_mmrsq_root[j][i]*(fl_dy[j][i+1])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dy[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+
fl_mdz[j][i])*(0.5));


fl_fz[j][i] = fl_kb*(-1.0*(fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])-fl_mmrsq_root[j][i]*(fl_dz[j][i+1])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5)+(fl_dz[j][i]*2.0)*fl_mmrsq_root2[j][i]*(fl_rsq[j][i+1])*(fl_mdx[j][i]+fl_mdy[j][i]+fl_mdz[j][i])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i];

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i];

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i];
}


}

i=(pol_N - 2);

fl_fx[j][i]= fl_kb*((-fl_mmrsq_root[j][i-1])*(fl_dx[j][i-1]-fl_dx[j][i])+fl_mmrsq_root[j][i-2]*(fl_dx[j][i-2])+((fl_dx[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dx[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-((fl_dx[j][i-1])*2.0)*(fl_mmrsq_root2[j][i-2])*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5));


fl_fy[j][i] = fl_kb*((-fl_mmrsq_root[j][i-1])*(fl_dy[j][i-1]-fl_dy[j][i])+fl_mmrsq_root[j][i-2]*(fl_dy[j][i-2])+((fl_dy[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dy[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5));


fl_fz[j][i]= fl_kb*((-fl_mmrsq_root[j][i-1])*(fl_dz[j][i-1]-fl_dz[j][i])+fl_mmrsq_root[j][i-2]*(fl_dz[j][i-2])+((fl_dz[j][i]*2.0)*(fl_rsq[j][i-1])-(fl_dz[j][i-1]*2.0)*(fl_rsq[j][i]))*fl_mmrsq_root2[j][i-1]*(fl_mdx[j][i-1]+fl_mdy[j][i-1]+fl_mdz[j][i-1])*(0.5)-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5));


fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dx[j][i])*fl_inv_rsq_root[j][i];

fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dy[j][i])*fl_inv_rsq_root[j][i];

fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1] + 2*fl_ks*(-fl_b + fl_rsq_root[j][i])*(fl_dz[j][i])*fl_inv_rsq_root[j][i];




i=(pol_N-1);


fl_fx[j][i] = fl_kb*((fl_mmrsq_root[j][i-2])*(fl_dx[j][i-2])-(fl_dx[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5));
fl_fy[j][i]= fl_kb*((fl_mmrsq_root[j][i-2])*(fl_dy[j][i-2])-(fl_dy[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5));
fl_fz[j][i]= fl_kb*((fl_mmrsq_root[j][i-2])*(fl_dz[j][i-2])-(fl_dz[j][i-1]*2.0)*fl_mmrsq_root2[j][i-2]*(fl_rsq[j][i-2])*(fl_mdx[j][i-2]+fl_mdy[j][i-2]+fl_mdz[j][i-2])*(0.5));


// spring force
fl_fsX[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dx[j][i-1])*fl_inv_rsq_root[j][i-1];
fl_fsY[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dy[j][i-1])*fl_inv_rsq_root[j][i-1];
fl_fsZ[j][i]=-2*fl_ks*(-fl_b + fl_rsq_root[j][i-1])*(fl_dz[j][i-1])*fl_inv_rsq_root[j][i-1];


}

//=== end of calculation for the (spring+bending) potential filament



//=================== Start of membrane(polymer) Lennard Jones potential===================

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
//=================== End of membrane(polymer) Leonard Jones potential===================


//========== Start of filament Lennard Jones potential --- within the beads of each filaments===========
for(j_pol=1;j_pol< num_fil;j_pol++)
{
pol_N = pol_N_of_fl[j_pol];
for (i=0;i<=(pol_N - 1);i++)
	{
	fl_F_LJx[j_pol][i]=0.0;
	fl_F_LJy[j_pol][i]=0.0;
	fl_F_LJz[j_pol][i]=0.0;
	}


   for (i=0;i<=(pol_N - 2);i++) {
     for (j=i+1;j<=(pol_N - 1);j++) {
		dx1=fl_x[j_pol][i]-fl_x[j_pol][j];
                dy1=fl_y[j_pol][i]-fl_y[j_pol][j];
                dz1=fl_z[j_pol][i]-fl_z[j_pol][j];



                r2=dx1*dx1+dy1*dy1+dz1*dz1;

	if (r2 < fl_rc) {// here fl_rc is also the "distance squared" between "filaments" beads
	  r6i   = pow(fl_rc,3)/(r2*r2*r2);

	  f     = -12.0*epsilon*(r6i*r6i-r6i);
	  fl_F_LJx[j_pol][i] += dx1*f/r2;
	  fl_F_LJx[j_pol][j] -= dx1*f/r2;
	  fl_F_LJy[j_pol][i] += dy1*f/r2;
	  fl_F_LJy[j_pol][j] -= dy1*f/r2;
	  fl_F_LJz[j_pol][i] += dz1*f/r2;
	  fl_F_LJz[j_pol][j] -= dz1*f/r2;



	}
     }


}
}
//========== End of filament Lennard Jones potential --- within the beads of each filaments ============


//======= Start of Filament-Membrane Lennard Jones potential calculation ========
for(j_pol=1;j_pol< num_fil;j_pol++)
{
	for(i=0;i<=(N-1);i++)
	{

	pol_wall_lj_totx[j_pol][i]=0.0; // total L-J potential on the polymer(membrane) beads due to the filaments beads
	pol_wall_lj_toty[j_pol][i]=0.0;// total L-J potential on the polymer(membrane) beads due to the filaments beads
	pol_wall_lj_totz[j_pol][i]=0.0;// define these variable and check the array limits
	}
}


for(j_pol=1;j_pol< num_fil;j_pol++)
{
	for(i=0;i<=(M-1);i++)
{
	fl_pol_wall_lj_totx[j_pol][i] = 0.0; // total L-J potential on the "filament beads" due to the "membrane beads"
	fl_pol_wall_lj_toty[j_pol][i] = 0.0;//i think initialization part needs to be kept out side "i.e. before the start of this loop"
	fl_pol_wall_lj_totz[j_pol][i] = 0.0;
}
}




for(j_pol=1;j_pol< num_fil;j_pol++)
{

pol_N = pol_N_of_fl[j_pol];




for(i=0;i<=(M-1);i++) // initialization is done for the "number of beads" in the "filaments" but the loop is done till the number of the number of beads in the filaments; there is no harm in doing so;
{
px[i] = 0.0;
py[i] = 0.0;
pz[i] = 0.0;
}



for(po_i = 0; po_i <= (pol_N -1); po_i++)
{
px[po_i] = fl_x[j_pol][po_i];
py[po_i] = fl_y[j_pol][po_i];
pz[po_i] = 0.0;

for (i=0;i<=(N-1);i++)
{
F_LJ_pw_x[j_pol][i]=0.0;
F_LJ_pw_y[j_pol][i]=0.0;
F_LJ_pw_z[j_pol][i]=0.0;

fil_F_LJ_pw_x[j_pol][i] = 0.0;
fil_F_LJ_pw_y[j_pol][i] = 0.0;
fil_F_LJ_pw_z[j_pol][i] = 0.0;
}


//E_lj_pw =0.0;

for (i=0;i<=(N-1);i++) // this loop will run till "N-1" because calculation is done for all the beads of membrane
	{
		dx1=-px[po_i]+x[i];
                dy1=-py[po_i]+y[i];
                dz1=-pz[po_i]+z[i];
              r2=dx1*dx1+dy1*dy1+dz1*dz1;

	if (r2 < rc_pw){ // "rc_pw" is the "distance squared" between the "filament" and the "polymer" beads
	  r6i   = pow(rc_pw , 3)/(r2*r2*r2);
	   //E_lj_pw    = E_lj_pw + (r6i*r6i - 2.0*r6i);
	  f     = -12.0*epsilon*(r6i*r6i-r6i);

	  F_LJ_pw_x[j_pol][i]  = dx1*f/r2;
      fil_F_LJ_pw_x[j_pol][i]  = dx1*f/r2;
	  F_LJ_pw_y[j_pol][i]  = dy1*f/r2;
      fil_F_LJ_pw_y[j_pol][i]  = dy1*f/r2;
	  F_LJ_pw_z[j_pol][i]  = dz1*f/r2;
      fil_F_LJ_pw_z[j_pol][i]  = dz1*f/r2;

	}

pol_wall_lj_totx[j_pol][i] = pol_wall_lj_totx[j_pol][i] + F_LJ_pw_x[j_pol][i] ; // L-J force on "polymer bead" due to the each "bead of filaments"
pol_wall_lj_toty[j_pol][i] = pol_wall_lj_toty[j_pol][i] + F_LJ_pw_y[j_pol][i] ;
pol_wall_lj_totz[j_pol][i] = pol_wall_lj_totz[j_pol][i] + F_LJ_pw_z[j_pol][i] ;


fl_pol_wall_lj_totx[j_pol][po_i] = fl_pol_wall_lj_totx[j_pol][po_i] - fil_F_LJ_pw_x[j_pol][i];//L-J force on "filament bead" due to the each "bead of polymer"
fl_pol_wall_lj_toty[j_pol][po_i] = fl_pol_wall_lj_toty[j_pol][po_i] - fil_F_LJ_pw_y[j_pol][i];
fl_pol_wall_lj_totz[j_pol][po_i] = fl_pol_wall_lj_totz[j_pol][po_i] - fil_F_LJ_pw_z[j_pol][i];


	}


}


}
//======= End of Filament-Membrane Lennard Jones potential calculation ========


//======== start of filament-filament Lennard Jones potential calculation
for(j_pol = 1; j_pol < num_fil; j_pol++)
{


	for(i=0;i<=(M-1); i++)
	{
	fl_fl_lj_totx[j_pol][i] = 0.0;
	fl_fl_lj_toty[j_pol][i] = 0.0;
	fl_fl_lj_totz[j_pol][i] = 0.0;
	}
}

for(ini_fil = 1; ini_fil < (num_fil-1); ini_fil++ ) // outer loop will run from filament number 1 to filament number 5
{

for(j_pol= (ini_fil + 1);j_pol< num_fil;j_pol++)
{

pol_N = pol_N_of_fl[j_pol];




for(i=0;i<=(M-1);i++) // initialization is done for the "number of beads" in the "filaments" but the loop is done till the maximum number of bead a filament can accomodate (array size of filaments beads); there is no harm in doing so;
{
px[i] = 0.0;
py[i] = 0.0;
pz[i] = 0.0;
}



for(po_i = 0; po_i <= (pol_N -1); po_i++)
{
px[po_i] = fl_x[j_pol][po_i];
py[po_i] = fl_y[j_pol][po_i];
pz[po_i] = 0.0;

for (i=0;i<=(pol_N_of_fl[ini_fil]-1);i++)
{
fl_fl_LJ_pw_x[i]=0.0; // check at the start for the declaration of this variable
fl_fl_LJ_pw_y[i]=0.0;
fl_fl_LJ_pw_z[i]=0.0;


}


//E_lj_pw =0.0;

for (i=0;i<=(pol_N_of_fl[ini_fil]-1);i++) // this loop will run till "N-1" because calculation is done for all the beads of membrane
	{
		dx1= -px[po_i]+fl_x[ini_fil][i];
                dy1= -py[po_i]+fl_y[ini_fil][i];
                dz1= -pz[po_i]+fl_z[ini_fil][i];
              r2=dx1*dx1+dy1*dy1+dz1*dz1;

	if (r2 < fl_fl_rsq ){ // "some_distance_square" is the "distance squared" between the "filament" and the "filament" beads
	  r6i   = pow(fl_fl_rsq,3)/(r2*r2*r2);
	   //E_lj_pw    = E_lj_pw + (r6i*r6i - 2.0*r6i);
	  f     = -12.0*epsilon*(r6i*r6i-r6i);

	  fl_fl_LJ_pw_x[i]  = dx1*f/r2;

	  fl_fl_LJ_pw_y[i]  = dy1*f/r2;

	  fl_fl_LJ_pw_z[i]  = dz1*f/r2;


	}

fl_fl_lj_totx[ini_fil][i] = fl_fl_lj_totx[ini_fil][i] + fl_fl_LJ_pw_x[i] ; // L-J force on "first filament bead" due to the each "bead of filaments"
fl_fl_lj_toty[ini_fil][i] = fl_fl_lj_toty[ini_fil][i] + fl_fl_LJ_pw_y[i] ;
fl_fl_lj_totz[ini_fil][i] = fl_fl_lj_totz[ini_fil][i] + fl_fl_LJ_pw_z[i] ;


fl_fl_lj_totx[j_pol][po_i] = fl_fl_lj_totx[j_pol][po_i] - fl_fl_LJ_pw_x[i] ;//L-J force on "filament bead" due to the each "bead of polymer"
fl_fl_lj_toty[j_pol][po_i] = fl_fl_lj_toty[j_pol][po_i] - fl_fl_LJ_pw_y[i];
fl_fl_lj_totz[j_pol][po_i] = fl_fl_lj_totz[j_pol][po_i] - fl_fl_LJ_pw_z[i];


	}


}


}

}
//======== end of filament-filament Lennard Jones potential calculation


//=== Start of Updating of membrane bead co-ordinates
for(i=1;i<=(N-2);i++)
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
sum_x[i] = pol_wall_lj_totx[1][i] + pol_wall_lj_totx[2][i] + pol_wall_lj_totx[3][i] + pol_wall_lj_totx[4][i] + pol_wall_lj_totx[5][i] + pol_wall_lj_totx[6][i]; // define sum_x[i] and sum_y[i] have to redine here(because increment of number of filament)
sum_y[i] = pol_wall_lj_toty[1][i] + pol_wall_lj_toty[2][i] + pol_wall_lj_toty[3][i] + pol_wall_lj_toty[4][i] + pol_wall_lj_toty[5][i] + pol_wall_lj_toty[6][i];


x[i]=x[i]-mu*(fx[i]+fsX[i]+F_LJx[i]+sum_x[i] + Sf_on_pol_bead_due_to_flX[i])+cn*(sum-6.0);
y[i]=y[i]-mu*(fy[i]+fsY[i]+F_LJy[i]+sum_y[i] + Sf_on_pol_bead_due_to_flY[i])+cn*(sum1-6.0);
z[i]=z[i];//-mu*(fz[i]+fsZ[i]+F_LJz[i]+F_LJ_pw_z[i])+cn*randforce(sum2,w);



}
//=== End of Updating of membrane bead co-ordinates





//=== Start Of updating filaments bead co-ordinates
for(j_pol = 1; j_pol < num_fil; j_pol++)
{
pol_N = pol_N_of_fl[j_pol];

for(i=0;i<=(pol_N - 1); i++)
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


	fl_x[j_pol][i]=fl_x[j_pol][i]-mu2*(fl_fx[j_pol][i] + fl_fsX[j_pol][i] + fl_F_LJx[j_pol][i] +
        fl_pol_wall_lj_totx[j_pol][i] + fl_fl_lj_totx[j_pol][i] + Sf_on_fl_bead_due_to_polX[j_pol][i]) + cn2*(sum - 6.0);

	fl_y[j_pol][i]=fl_y[j_pol][i]-mu2*(fl_fy[j_pol][i] + fl_fsY[j_pol][i] + fl_F_LJy[j_pol][i] +
         fl_pol_wall_lj_toty[j_pol][i]+fl_fl_lj_toty[j_pol][i] + Sf_on_fl_bead_due_to_polY[j_pol][i]) + cn2*(sum1-6.0);

	fl_z[j_pol][i]=fl_z[j_pol][i];//-mu*(fz[i]+fsZ[i]+F_LJz[i]+F_LJ_pw_z[i])+cn*randforce(sum2,w);

	}
}

//=== End  Of updating filaments bead co-ordinates



pol_N = pol_N_of_fl[1];

for(i=0; i<=(pol_N-1); i++)
{
	if(fl_y[1][i] < bead1_y_limit)
	{
	fl_y[1][i] = bead1_y_limit ;
	}
}

pol_N = pol_N_of_fl[4];
for(i=0; i<=(pol_N-1); i++)
{
	if(fl_y[4][i] > bead2_y_limit)
	{
	fl_y[4][i] = bead2_y_limit ;
	}
}






//===== start of section to figure out if the bead can be added to the "filaments"
for(j_pol = 1; j_pol < num_fil; j_pol++)
{

pol_N = pol_N_of_fl[j_pol];
	min_dis = 5000.0;
	for(i=0; i<=(N-1); i++)
	{
	dis_pw[i] = 0.0;
	}


x_front_bead[j_pol] = fl_x[j_pol][pol_N-1];
y_front_bead[j_pol] = fl_y[j_pol][pol_N-1];

for(i=0;i<=(N-1);i++)
{
min_disx = x[i] - fl_x[j_pol][pol_N-1];
min_disy = y[i] - fl_y[j_pol][pol_N-1];
min_disz = z[i] - 0.0;

dis_pw[i] = min_disx*min_disx + min_disy*min_disy + min_disz*min_disz ;


if(sqrt(dis_pw[i]) < min_dis)
{
min_dis = sqrt(dis_pw[i]);
shortest_dis_bead = i; // this is to remember that where the shortest distance is?
}
}

if( min_dis >= sqrt(rq_d_sq))
{
pol_N_of_fl[j_pol] = pol_N_of_fl[j_pol] + 1;
increased_bead = pol_N_of_fl[j_pol];


fl_x[j_pol][increased_bead - 1] = (( min_dis - 2.0)*fl_x[j_pol][pol_N-1] + 2.0*(x[shortest_dis_bead]))/min_dis ;

fl_y[j_pol][increased_bead - 1] = ((min_dis  - 2.0)*fl_y[j_pol][pol_N-1] + 2.0*(y[shortest_dis_bead]))/min_dis ;
fl_z[j_pol][increased_bead - 1] = 0.0;
}
else
nothing = 0;



if(j_pol==1)
{
br_bead_one_new_value = pol_N_of_fl[j_pol] ;
}

if(j_pol==4)
{
 br_bead_four_new_value = pol_N_of_fl[j_pol];
}

/*

br_bead_one have to change this value
if(br_bead_one_given > br_bead_one_new_value)
{
br_bead_one = br_bead_one + 1;

br_bead_one_given = br_bead_one_new_value;
}

if(br_bead_four_given > br_bead_four_new_value)
{
br_bead_four = br_bead_four + 1;

br_bead_four_given = br_bead_four_new_value;
}

*/

}
//===== end of section to figure out if the bead can be added to the "filaments"





/*
end_dis=sqrt((x[0] - x[N-1])*(x[0] - x[N-1]) +(y[0] - y[N-1])*(y[0] - y[N-1])+(z[0] - z[N-1])*(z[0] - z[N-1]));
if(t_flag ==0)
{
if(end_dis <= 1.0)
{fprintf(fp11, "%d \n", t);
t_flag = -1;}
}
*/


//printf("%d %d %d %d %d %d %d\n",timep,pol_N_of_fl[1],pol_N_of_fl[2],pol_N_of_fl[3],pol_N_of_fl[4],pol_N_of_fl[5],pol_N_of_fl[6]);
if( t >= timep )
{
//length_pol[1] = sqrt((x[36] - x[20])*(x[36] - x[20]) +(y[36] - y[20])*(y[36] - y[20])+(z[36] - z[20])*(z[36] - z[20]));
//length_pol[2] = sqrt((x[140] - x[124])*(x[140] - x[124]) +(y[140] - y[124])*(y[140] - y[124])+(z[140] - z[124])*(z[140] - z[124]));
//fprintf(fp1, "%d %f %f %d %d\n",timep,length_pol[1],length_pol[2],pol_N_of_fl[1],pol_N_of_fl[2]);
fprintf( fp1,"%d %d %d %d %d %d %d\n",timep,pol_N_of_fl[1],pol_N_of_fl[2],pol_N_of_fl[3],pol_N_of_fl[4],pol_N_of_fl[5],pol_N_of_fl[6]);
sprintf( Num ,"Wall_coordinate%d.xyz",timep);
fp10=fopen(Num, "w");
fprintf(fp10, "%d\n",(N+6*M));
fprintf(fp10, "timestep %d\n",timep );
bead_count = 0;
for(i=0;i<=(N-1);i++)
{

fprintf(fp10,"%d  %f %f %f \n",1,x[i],y[i],z[i]);
bead_count = bead_count + 1;
}


for(pol_num = 1;pol_num < num_fil;pol_num++)
{
pol_N = pol_N_of_fl[pol_num];
//j=0;
for(i=0;i<=(pol_N - 1);i++)
{

pol_x = fl_x[pol_num][i];
pol_y = fl_y[pol_num][i];
pol_z = 0.0;
//j = j+1;
fprintf(fp10,"%d  %f %f %f \n",2,pol_x,pol_y,pol_z);
bead_count = bead_count + 1;
}

}


do
{
pol_x = 0.0;
pol_y = 0.0;
pol_z = 0.0;

fprintf(fp10,"%d  %f %f %f \n",2,pol_x,pol_y,pol_z);
bead_count = bead_count + 1;
}while(bead_count < (N+6*M));


timep = timep + time_step;

fclose(fp10);
}

} // time loop
fclose(fp1); //  file bead_count
fclose(fp11);
//fclose(fp3); // file closing polymer and wall
//fclose(fp4); // file closing bead and bead




return 0;

}


float randforce(float sum,float w)
 {
 float randforce;
 randforce=(sum-6.0)/(sqrt(w));

return randforce;
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
