clc;clear all;close all;
%% Data
mesh=25;
N_r=10;
L=7;%m
T_inner_fluid=15; %celcius
P_inner=10*10^2;      %kp

T_outer_fluid=175;   %

T_inf=20;        %Air Temprature

D1_inner=8*10^-3;%m
D2_inner=9.5*10^-3;

D1_outer=14.5*10^-3;%m
D2_outer=16*10^-3;
D3_outer=36*10^-3;    %Insulating thickness

G_inner=300; %kg/s/m2
G_outer=100;

epsilon=1.5*10^-6; %m
R=8.314;

%% hydraulic diameter-permimeter-area calculation
D_h_inner=D1_inner;
D_h_outer=(D1_outer^2-D2_inner^2)/D2_inner;

P1_inner=pi*D1_inner;
P2_inner=pi*D2_inner;
P1_outer=pi*D1_outer;
P2_outer=pi*D2_outer;
P3_outer=pi*D3_outer;

A1_inner=pi*D_h_inner^2/4;
A1_outer=pi*(D_h_outer)^2/4;
A_inner_wall=pi*(D2_inner^2-D1_inner^2)/4;
A_outer_wall=pi*(D2_outer^2-D1_outer^2)/4;

%% laoding spline data
load spline_hl_T_coeff splined_hl_T_coeff
load spline_T_P_coeff splined_T_P_coeff
load spline_P_T_coeff splined_P_T_coeff
load spline_T_vl_coeff splined_T_vl_coeff
load spline_T_vg_coeff splined_T_vg_coeff
load spline_T_hl_coeff splined_T_hl_coeff
load spline_T_hg_coeff splined_T_hg_coeff

%% assumed data
k_inner=0.585; %W/(m*c)
k_outer_g=2.77; %W/(m*c)
k_outer_l=.585;

ro_inner=985; %kg/m3 average
mu_inner=0.8 *10^-3; %centipoise to SI

ro_outer_l=985; %kg/m3 average
mu_outer_l=0.2 *10^-3; %centipoise to SI
ro_outer_g=5.147; %kg/m3
mu_outer_g=0.015*10^-3; %centipoise to SI

thermal_cond_wall=55;  %w/(m c) for chrome steel ;
thermal_cond_ins=0.1;
    
inner_initial_entalpy=ppval(splined_T_hl_coeff,T_inner_fluid); 
outer_initial_entalpy=ppval(splined_T_hg_coeff,T_outer_fluid);

%% initial velocities
v_inner=G_inner/ro_inner;
v_outer=G_outer/ro_outer_g;

%% making variables
dz=L/mesh;
inner_f.P=ones(1,mesh+1)*P_inner;
inner_f.T=ones(1,mesh+1)*T_inner_fluid;
inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;
inner_f.hConv=ones(1,mesh);   

inner_w.T=zeros(2,mesh);

outer_f.P=ones(1,mesh+1)*ppval(splined_T_P_coeff,T_outer_fluid);
outer_f.T=ones(1,mesh+1)*T_outer_fluid;
outer_f.h=ones(1,mesh+1)*outer_initial_entalpy;
outer_f.x=ones(1,mesh+1);
outer_f.void=zeros(1,mesh);
outer_f.hConv=ones(1,mesh); 

outer_w.T=zeros(N_r,mesh);          %temperature
outer_w.hNat=ones(1,mesh)*10;      % This value is assumed for air natural convection
outer_T=linspace(T_outer_fluid-10,T_inf+10,N_r);
for i=1:N_r
    outer_w.T(i,:)=outer_T(i);
end

%% starting the program
inner_f.cp=Cp_inner(inner_f.T(1));
outer_f.cp=Cp_outer(outer_f.T(1),1);

Re_inner=ro_inner*v_inner*D_h_inner/mu_inner;
Re_outer=ro_outer_g*v_outer*D_h_outer/mu_outer_g;

Pr_inner=inner_f.cp*mu_inner/k_inner;
Pr_outer=outer_f.cp*mu_outer_g/k_outer_g;

Nu_inner = Gnieliski(Re_inner,Pr_inner(1));
Nu_outer = Monrad(Re_outer(1),Pr_outer(1),D1_outer,D2_inner);

inner_f.hConv=Nu_inner*k_inner/D_h_inner;
outer_f.hConv=Nu_outer*k_outer_g/D_h_outer;

temp=(inner_f.hConv*P1_inner*inner_f.T+outer_f.hConv*P2_inner*outer_f.T)/(inner_f.hConv*P1_inner+outer_f.hConv*P2_inner);
inner_w.T(1,:)=0.8*temp(1:mesh);
inner_w.T(2,:)=temp(1:mesh);

%% Solving the PDEs of Heat Exchanger simultaneoutsly (successive substitution method)
CRITERIA=1;
TEMP.inner_f_h=zeros(1,mesh+1);
TEMP.outer_f_h=zeros(1,mesh+1);
TEMP.inner_w_T=zeros(1,mesh);

while CRITERIA > 10^-5
    
    %% inner pipe simulation
    TEMP.inner_f_h=inner_f.h;
    f_inner=Churchill(Re_inner,epsilon/D1_inner);   %this line is for assumption.
    for i=1:mesh
       inner_f.P(i+1)=inner_f.P(i)-1e-3*dz*P1_inner/A1_inner*((f_inner/4) * G_inner^2/(2*ro_inner));
       ERROR=1;
       while ERROR>1e-5
            temp=inner_f.h(i+1);
            T_mean=0.5*(inner_f.T(i)+inner_f.T(i+1));
            Pr_inner=Cp_inner(T_mean)*mu_inner/k_inner;
            Nu=Gnieliski(Re_inner,Pr_inner);
            inner_f.hConv(i)=k_inner*Nu/D_h_inner;
            inner_f.h(i+1)=inner_f.h(i)+ 1e-3*(inner_f.hConv(i)*P1_inner*dz*(inner_w.T(i)-T_mean)/(G_inner*A1_inner));
            inner_f.T(i+1)=ppval(splined_hl_T_coeff,inner_f.h(i+1));
            ERROR=abs((inner_f.h(i+1)-temp)/inner_f.h(i+1));
        end
    end

    %% outer pipe simulation
    TEMP.outer_f_h=outer_f.h;
    TEMP.outer_f_T=outer_f.T;
    for i=1:mesh
        ERROR=1;
        outer_f.T(i+1)=outer_f.T(i);
        outer_f.x(i+1)=outer_f.x(i)  - 0.001;  
        while ERROR>10^-5

            T_mean=0.5*(outer_f.T(i)+outer_f.T(i+1));
            x_mean=0.5*(outer_f.x(i)+outer_f.x(i+1));

            outer_f.P(i+1)=ppval(splined_T_P_coeff,outer_f.T(i+1));

            all_x=[outer_f.x(i) x_mean outer_f.x(i+1)];
            all_ro_l=ppval(splined_T_vl_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
            all_ro_g=ppval(splined_T_vg_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;


            hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
            hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

            outer_f.h(i+1)=outer_f.x(i+1)*hg_outer+(1-outer_f.x(i+1))*hl_outer;
            previous_h=outer_f.h(i+1);
            previous_x=outer_f.x(i+1);

            alpha = Smith(all_x,all_ro_l,all_ro_g);
            mu_outer_ave=mu_outer_g; 
            Re_outer=G_outer*D_h_outer/mu_outer_ave;
            outer_f.cp=Cp_outer(T_mean,x_mean);
            k_outer_ave=x_mean*k_outer_g+(1-x_mean)*k_outer_l;
            Pr_outer=outer_f.cp*mu_outer_ave/k_outer_ave;
            Nu_outer = Monrad(Re_outer,Pr_outer,D1_outer,D2_inner);
            outer_f.hConv(i)=Nu_outer*k_outer_ave/D_h_outer;

            phi_friedel=1;
            G_l_outer=(1-all_x).*G_outer;
            all_vl=G_l_outer./all_ro_l;
            G_g_outer=all_x.*G_outer;
            all_vg=G_g_outer./all_ro_g;

            ro_tp=alpha(2)*all_ro_g(2)+(1-alpha(2))*all_ro_l(2);
            outer_f.void(i)=alpha(2);
            all_kapa=all_x.*all_vg+(1-all_x).*all_vl;
            f_outer=Churchill(Re_outer,epsilon/D_h_outer);

            outer_f.P(i+1)=outer_f.P(i)-1e-3*(   dz*P2_inner/A1_outer*((phi_friedel*f_outer/4) * (G_outer^2/(2*ro_tp)) )...
                                                + G_outer*(all_kapa(3)-all_kapa(1)));
            outer_f.h(i+1)=outer_f.h(i)+ 1e-3*(  outer_f.hConv(i)*P2_inner*dz*(inner_w.T(2,i)-T_mean)/(G_outer*A1_outer)...
                                               + outer_f.hConv(i)*P1_outer*dz*(outer_w.T(1,i)-T_mean)/(G_outer*A1_outer)...
                                               - 0.5*(all_kapa(3)-all_kapa(1)) );

            outer_f.T(i+1)=ppval(splined_P_T_coeff,outer_f.P(i+1));

            hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
            hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

            outer_f.x(i+1)=(outer_f.h(i+1)-hl_outer)/(hg_outer-hl_outer);
            ERROR=abs((outer_f.h(i+1)-previous_h)/outer_f.h(i+1))+abs((outer_f.x(i+1)-previous_x)/outer_f.x(i+1));
        end   
    end

    %% inner wall simulation
    TEMP.inner_w_T=inner_w.T;
    thermal_cond_wall=55;  %w/(m c) for chrome steel ;
    dr=(D2_inner-D1_inner)/4;
    wall_matrix=zeros(2*mesh,2*mesh);wall_vector=zeros(2*mesh,1);%wall_Temp=zeros(2*mesh,1);
    for i=1:2*mesh
        if i>mesh
             Coe=0;
             if i>mesh+1
                wall_matrix(i,i-1)=thermal_cond_wall*dr*(D2_inner^2-(D1_inner+dr)^2);
                Coe=Coe+1;
             end
             if i<2*mesh
                wall_matrix(i,i+1)=thermal_cond_wall*dr*(D2_inner^2-(D1_inner+dr)^2);
                Coe=Coe+1;
             end
             wall_matrix(i,i)=-( outer_f.hConv(i-mesh)*D2_inner*dz^2*dr + Coe*dr*(-(D1_inner+dr)^2+D2_inner^2) + 0.5*thermal_cond_wall*(D1_inner+dr)*dz^2);
             wall_matrix(i,i-mesh)= 0.5*thermal_cond_wall*(D1_inner+dr)*dz^2;                
             wall_vector(i)=-0.5*outer_f.hConv(i-mesh)*D2_inner*dz^2*dr*(outer_f.T(i-mesh)+outer_f.T(i+1-mesh));
        else
             Coe=0;
             if i>1
                wall_matrix(i,i-1)=thermal_cond_wall*dr*((D1_inner+dr)^2-D1_inner^2);
                Coe=Coe+1;
             end
             if i<mesh
                wall_matrix(i,i+1)=thermal_cond_wall*dr*((D1_inner+dr)^2-D1_inner^2);
                Coe=Coe+1;
             end
             wall_matrix(i,i)=-( inner_f.hConv(i)*D1_inner*dz^2*dr + Coe*dr*((D1_inner+dr)^2-D1_inner^2) +0.5* thermal_cond_wall*(D1_inner+dr)*dz^2);
             wall_matrix(i,i+mesh)=0.5* thermal_cond_wall*(D1_inner+dr)*dz^2;  
             wall_vector(i)=-0.5*inner_f.hConv(i)*D1_inner*dz^2*dr*(inner_f.T(i)+inner_f.T(i+1));
        end
        
        
    end
    wall_Temp=(wall_matrix\wall_vector)';
    inner_w.T(1,:)=wall_Temp(1:mesh);
    inner_w.T(2,:)=wall_Temp(mesh+1:2*mesh);
    
    
    %% outer wall simulation
    RaD=(9.81*D3_outer^3)*(outer_w.T(N_r,:)-T_inf*ones(1,mesh))/(15.68*1e-6*2.216*1e-5*T_inf);
    outer_w.hNat=(0.02624/D3_outer)* ( 0.36+ (0.518*RaD.^0.25)/( 1+(0.559/0.708)^(9/16) )^(4/9) );
    
    dr1=0.5*(D2_outer-D1_outer);
    dr2=(D3_outer-D2_outer)/(2*(N_r-1));
    
    wall_matrix=zeros(N_r*mesh,N_r*mesh);wall_vector=zeros(N_r*mesh,1);
    for i=1:mesh
        if i>1
            wall_matrix(i,i-1)=-thermal_cond_wall*A_outer_wall/dz;
        end
        wall_matrix(i,i)=2*thermal_cond_wall*A_outer_wall/dz + outer_f.hConv(i)*P1_outer*dz + thermal_cond_ins*P2_outer*dz/dr1;
        wall_matrix(i,i+mesh)= -thermal_cond_ins*P2_outer*dz/dr1;
        
        if i<mesh
            wall_matrix(i,i+1)=-thermal_cond_wall*A_outer_wall/dz;
        end
        wall_vector(i)=0.5*(outer_f.T(i)+outer_f.T(i+1))*outer_f.hConv(i)*P1_outer*dz;
    end
    j=2;  counter=1;
    for i=mesh+1:(N_r-1)*mesh
        Ds=D2_outer+(j-2)*dr2*2;
        Dn=D2_outer+(j-1)*dr2*2;
        A_ew=pi*(Dn^2-Ds^2)/4;
        A_s=pi*Ds*dz;
        A_n=pi*Dn*dz;
        if counter>1
            wall_matrix(i,i-1)=-thermal_cond_ins*A_ew/dz;
        end
        wall_matrix(i,i)= thermal_cond_ins* ( 2*A_ew/dz + (A_s+A_n)/dr2 );   
        wall_matrix(i,i-mesh)=-thermal_cond_ins*A_s/dr2;
        wall_matrix(i,i+mesh)=-thermal_cond_ins*A_n/dr2;
        if counter<mesh
            wall_matrix(i,i+1)=-thermal_cond_ins*A_ew/dz;
        end
        wall_vector(i)=0;  
        counter=counter+1;
        if counter>mesh
            counter=1;
            j=j+1;
        end
    end
    for i=(N_r-1)*mesh+1:N_r*mesh
        Ds=D2_outer+(j-2)*dr2*2;
        Dn=D2_outer+(j-1)*dr2*2;
        A_ew=pi*(Dn^2-Ds^2)/4;
        A_s=pi*Ds*dz;
        A_n=pi*Dn*dz;
        if counter>1
            wall_matrix(i,i-1)=-0.5*thermal_cond_ins*A_ew/dz;
        end
        wall_matrix(i,i)=thermal_cond_ins*A_ew/dz + outer_w.hNat(i-(j-1)*mesh)*P3_outer*dz + thermal_cond_ins*A_s/dr1;
        wall_matrix(i,i-mesh)= -thermal_cond_ins*A_s/dr1;
        if counter<mesh
            wall_matrix(i,i+1)=-0.5*thermal_cond_ins*A_ew/dz;
        end
        wall_vector(i)=outer_w.hNat(i-(j-1)*mesh)*P3_outer*dz*T_inf; 
        counter=counter+1;
    end
    T_vector=(wall_matrix\wall_vector)';
    for j=1:N_r
        outer_w.T(j,:)=T_vector((j-1)*mesh+1:j*mesh);
    end

    
    %% Defining the criteria for convergence of all PDEs
    CRITERIA1=norm(abs(inner_f.h-TEMP.inner_f_h)./inner_f.h);
    CRITERIA2=norm(abs(outer_f.h-TEMP.outer_f_h)./outer_f.h);
    CRITERIA=CRITERIA1+CRITERIA2;
end

figure('name','Inner Wall Temprature')
plot((0:mesh-1)*dz+(dz/2),inner_w.T(1,:),'.')
hold on
title('changing temperature in the inner pipe wall')
plot((0:mesh-1)*dz+(dz/2),inner_w.T(2,:),'r.')
xlabel('length(m)');ylabel('Inner wall Temperature(C)');

figure('name','Outer Wall Temprature')
for i=1:N_r
   plot((0:mesh-1)*dz+(dz/2),outer_w.T(i,:),'.')
   hold on
end
xlabel('length(m)');ylabel('Outer wall Temperature(C)');
title('trend of changing temperature in outer pipe and isolation')    