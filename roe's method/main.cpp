#include<bits/stdc++.h>
#include<iostream>
#include<fstream>
#include<cmath>

#define gamma 1.4
#define n 400
#define CFL 0.8
#define length_of_tube 1
#define enfix 0.000001

using namespace std;


double rho_l=1.0,rho_r=0.125,u_l=0,u_r=0,s_l,s_r,e_l,e_r,p_l=1.0,p_r=0.1,H_l,H_r;
double u_tilda[n],H_tilda[n],lambda[4],F1_l,F2_l,F3_l,F1_r,F2_r,F3_r,dummy1=0,dummy2=0,dummy3=0;
double delU1,delU2,delU3,alpha[4],k1[4],k2[4],k3[4];
double F1[n],F2[n],F3[n],U1[n],U2[n],U3[n];
double U_star_l1,U_star_l2,U_star_l3,U_star_r1,U_star_r2,U_star_r3;

double flux_1(double rho, double u){ return rho*u;}

double flux_2(double rho, double u,double p){return (rho*u*u+p);}

double flux_3(double rho, double u,double p)
{
    double e;
    e=(u*u)/2+(p/rho)*(1/(gamma-1));
    return u*(rho*e+p);
}

double maximum(double A[], int size)
{
	int i;
	double t;
	t=A[1];
	for(i=2;i<size-1;i++){
    	if(A[i]>t)
		t=A[i];
	}
	return(t);
}


int main()

{

    int i,j,k,count=0;

    double l=length_of_tube,delt,delx,a_l,a_r,a_tilda[n],rho[n],u[n],p[n],e[n],ie[n],x_node[n],x_center[n],lambda_max[n],z=0,t=0,part2=0;

    delx=(double)l/n;

  ofstream file,file2,file3,file4;
  file.open("densityPlotRoe.dat");
  file2.open("velocityPlotRoe.dat");
  file3.open("pressurePlotRoe.dat");
  file4.open("energyPlotRoe.dat");

    for (i=0;i<n;i++)
    {
        x_node[0]=(0);
        x_node[i+1]=x_node[i]+delx;

        x_center[i]=x_node[i]+delx/2;
        x_center[0]=delx/2;

        if(x_center[i]<(l/2))
        {
            rho[i]=rho_l;
            p[i]=p_l;
            u[i]=u_l;
            e[i]=(u[i]*u[i])/2+(p[i]/rho[i])*(1/(gamma-1));
        }

        else
        {
            rho[i]=rho_r;
            p[i]=p_r;
            u[i]=u_r;
            e[i]=(u[i]*u[i])/2+(p[i]/rho[i])*(1/(gamma-1));
        }
    }

    //cout<<"\t"<<"i"<<"\t"<<"x[i]"<<"\t"<<"rho[i]"<<"\t"<<"u[i]"<<"\t"<<"p[i]"<<endl;
    for (i=0;i<n;i++){

        //cout<<"\t"<<i<<"\t"<<x_center[i]<<"\t"<<rho[i]<<"\t"<<u[i]<<"\t"<<e[i]<<endl;
    }

    while(t<0.25)
    {
        for(i=1;i<n;i++)
        {

            rho_l=rho[i-1];     rho_r=rho[i];
            u_l=u[i-1];         u_r=u[i];
            p_l=p[i-1];         p_r=p[i];

            e_l=(u_l*u_l)/2+(p_l/rho_l)*(1/(gamma-1));
            e_r=(u_r*u_r)/2+(p_r/rho_r)*(1/(gamma-1));
            H_l=e_l+(p_l/rho_l);    H_r=e_r+(p_r/rho_r);



            U1[i-1]=rho_l;          U1[i]=rho_r;
            U2[i-1]=rho_l*u_l;      U2[i]=rho_r*u_r;
            U3[i-1]=rho_l*e_l;      U3[i]=rho_r*e_r;


            delU1=U1[i]-U1[i-1];
            delU2=U2[i]-U2[i-1];
            delU3=U3[i]-U3[i-1];


            u_tilda[i]=(sqrt(rho_l)*u_l+sqrt(rho_r)*u_r)/(sqrt(rho_l)+sqrt(rho_r));
            H_tilda[i]=(sqrt(rho_l)*H_l+sqrt(rho_r)*H_r)/(sqrt(rho_l)+sqrt(rho_r));
            a_tilda[i]=sqrt((gamma-1)*(H_tilda[i]-(pow(u_tilda[i],2))/2));


            lambda[1]=u_tilda[i] - a_tilda[i];
            lambda[2]=u_tilda[i];
            lambda[3]=u_tilda[i] + a_tilda[i];


            for(j=1;j<=3;j++){
                if(fabs(lambda[j])<enfix)
                lambda[j]=(((lambda[j]*lambda[j])/enfix)+enfix)/2;
            }
            //cout<<lambda[1]<<"  "<<lambda[2]<<"  "<<lambda[3]<<endl;
            //cout<<i<<"  "<<delU1<<"\t"<<delU2<<"  "<<delU3<<endl;


            k1[1]=1;    k2[1]=lambda[1];    k3[1]=H_tilda[i]-(u_tilda[i]*a_tilda[i]);
            k1[2]=1;    k2[2]=lambda[2];    k3[2]=(u_tilda[i]*u_tilda[i])/2;
            k1[3]=1;    k2[3]=lambda[3];    k3[3]=H_tilda[i]+(u_tilda[i]*a_tilda[i]);


            alpha[2]=(gamma-1)*(delU1*(H_tilda[i]-(u_tilda[i]*u_tilda[i]))+(u_tilda[i]*delU2)-delU3)/(a_tilda[i]*a_tilda[i]);
            alpha[1]=(delU1*(u_tilda[i]+a_tilda[i])-delU2-(a_tilda[i]*alpha[2]))/(2*a_tilda[i]);
            alpha[3]=delU1-(alpha[1]+alpha[2]);


            F1_l=flux_1(rho_l,u_l);         F1_r=flux_1(rho_r,u_r);
            F2_l=flux_2(rho_l,u_l,p_l);     F2_r=flux_2(rho_r,u_r,p_r);
            F3_l=flux_3(rho_l,u_l,p_l);     F3_r=flux_3(rho_r,u_r,p_r);


            dummy1=0;dummy2=0;dummy3=0;

            for(j=1;j<=3;j++){
                if(lambda[j]<=0){
                    dummy1+=(alpha[j]*lambda[j]*k1[j]);
                    dummy2+=(alpha[j]*lambda[j]*k2[j]);
                    dummy3+=(alpha[j]*lambda[j]*k3[j]);
                }
            }
            //cout<<i<<"  "<<dummy1<<"\t"<<dummy2<<"  "<<dummy3<<endl;

            F1[i]=F1_l+dummy1;
            F2[i]=F2_l+dummy2;
            F3[i]=F3_l+dummy3;

           // cout<<i<<"  "<<F1[i]<<"\t"<<F2[i]<<"\t"<<F3[i]<<endl;
        }

       	for(i=1;i<=n-1;i++){

			lambda_max[i]=fabs(u_tilda[i])+a_tilda[i];
	    }

		delt=(CFL*(delx)/maximum(lambda_max, n));
        t+=delt;

        for(i=1;i<n-1;i++){

            U1[i]=U1[i]-(delt*(F1[i+1]-F1[i])/delx);
            U2[i]=U2[i]-(delt*(F2[i+1]-F2[i])/delx);
            U3[i]=U3[i]-(delt*(F3[i+1]-F3[i])/delx);

        }

        for(i=1;i<n-1;i++){

            p[i]=U3[i]*(gamma-1)-(U2[i]*U2[i]*(gamma-1)/(2*U1[i]));
            rho[i]=U1[i];
            u[i]=U2[i]/U1[i];

        }

        count++;
        t+=delt;

    }

    for(i=0;i<n;i++){
     ie[i]= p[i]/((gamma-1)*rho[i]);
    }

     for(i=0;i<n;i++){

           file<<x_center[i]<<"\t"<<rho[i]<<endl;
           file2<<x_center[i]<<"\t"<<u[i]<<endl;
           file3<<x_center[i]<<"\t"<<p[i]<<endl;
           file4<<x_center[i]<<"\t"<<ie[i]<<endl;

        }
    return 0;
    }
