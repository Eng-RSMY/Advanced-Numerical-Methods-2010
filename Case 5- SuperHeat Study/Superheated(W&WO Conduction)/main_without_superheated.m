clc;clear all;close all;
%% Data
mesh=25;
L=3;%m
T_inner_fluid=22; %celcius
P_inner=5*10^2;      %kp

T_outer_fluid=400; %
P_outer_fluid=5*10^2;

D1_inner=9.5*10^-3;%m
D2_inner=11*10^-3;

D1_outer=14.5*10^-3;%m
D2_outer=16*10^-3;

G_inner=180; %kg/s/m2
G_outer=90;

epsilon=1.5*10^-6; %m
R=8.314;
tic

%% hydraulic diameter-permimeter-area calculation
D_h_inner=D1_inner;
D_h_outer=(D1_outer^2-D2_inner^2)/D2_inner;

P1_inner=pi*D1_inner;
P2_inner=pi*D2_inner;
P1_outer=pi*D1_outer;

A1_inner=pi*D_h_inner^2/4;
A1_outer=pi*(D_h_outer)^2/4;
A_inner_wall=pi*(D2_inner^2-D1_inner^2)/4;

%% laoding  data
load spline_hl_T_coeff splined_hl_T_coeff
load spline_T_P_coeff splined_T_P_coeff
load spline_P_T_coeff splined_P_T_coeff
load spline_T_vl_coeff splined_T_vl_coeff
load spline_T_vg_coeff splined_T_vg_coeff
load spline_T_hl_coeff splined_T_hl_coeff
load spline_T_hg_coeff splined_T_hg_coeff

load hSuperheated
load VSuperheated
load viscosity

T_vector=hSuperheated(2:end,1)';
P_vector=hSuperheated(1,2:end);
h_matrix=hSuperheated(2:end,2:end)';
V_matrix=VSuperheated(2:end,2:end)';
T_vis=viscosity(2:end,1)';
P_vis=viscosity(1,2:end);
mo_matrix=viscosity(2:end,2:end)';

%% assumed data
k_inner=0.585; %W/(m*c)
k_outer_g=2.77; %W/(m*c)
k_outer_l=.585;

ro_inner=985; %kg/m3 average
mu_inner=0.8 *10^-3; %centipoise to SI

ro_outer_g=1/interp2(T_vector,P_vector,V_matrix,T_outer_fluid,P_outer_fluid); %kg/m3
mu_outer_g= 1e-3*interp2(T_vis,P_vis,mo_matrix,T_outer_fluid,P_outer_fluid);        %0.015*10^-3; %centipoise to SI

%% initials
inner_initial_entalpy=ppval(splined_T_hl_coeff,T_inner_fluid); %92.207 for saturated liquid at 3kpa 
outer_initial_entalpy=interp2(T_vector,P_vector,h_matrix,T_outer_fluid,P_outer_fluid);
v_inner=G_inner/ro_inner;
v_outer=G_outer/ro_outer_g;

%% making variables
dz=L/mesh;
inner_f.P=ones(1,mesh+1)*P_inner;
inner_f.T=ones(1,mesh+1)*T_inner_fluid;
inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;
inner_f.hConv=ones(1,mesh);   

inner_w.T=zeros(1,mesh);

outer_f.P=ones(1,mesh+1)*P_outer_fluid;          
outer_f.T=ones(1,mesh+1)*T_outer_fluid;
outer_f.h=ones(1,mesh+1)*outer_initial_entalpy;
outer_f.hConv=ones(1,mesh);              
outer_f.vg=ones(1,mesh+1);
outer_f.mo=ones(1,mesh+1);

outer_w.T=zeros(1,mesh);%temperature

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
inner_w.T=linspace(temp(1),0.6*(T_inner_fluid+T_outer_fluid),mesh);
outer_w.T=outer_f.T(1:mesh);

%% Solving the PDEs of Heat Exchanger simultaneoutsly (successive substitution method)
CRITERIA=1;
TEMP.inner_f_h=zeros(1,mesh+1);
TEMP.outer_f_h=zeros(1,mesh+1);
TEMP.inner_w_T=zeros(1,mesh);

while CRITERIA > 10^-3
    
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
        outer_f.T(i+1)=0.998*outer_f.T(i);
        outer_f.P(i+1)=0.999*outer_f.P(i);   
        while ERROR>10^-5

            T_mean=0.5*(outer_f.T(i)+outer_f.T(i+1));
            P_mean=0.5*(outer_f.P(i)+outer_f.P(i+1));

            all_T=[outer_f.T(i) T_mean outer_f.T(i+1)];
            all_P=[outer_f.P(i) P_mean outer_f.P(i+1)];
            outer_f.h(i+1)=interp2(T_vector,P_vector,h_matrix,outer_f.T(i+1),outer_f.P(i+1));
            all_ro_g=interp2(T_vector,P_vector,V_matrix,all_T,all_P).^(-1);
            all_mo_g=1e-3*interp2(T_vis,P_vis,mo_matrix,all_T,all_P);
            outer_f.mo(i)=all_mo_g(1);outer_f.mo(i+1)=all_mo_g(3);
            
            previous_h=outer_f.h(i+1);
            previous_T=outer_f.T(i+1);
            mu_outer_ave=all_mo_g(2);                     
            Re_outer=G_outer*D_h_outer/mu_outer_ave;
            outer_f.cp=Cp_outer(T_mean,1);
            k_outer_ave=k_outer_g;                      
            Pr_outer=outer_f.cp*mu_outer_ave/k_outer_ave;
            Nu_outer = Monrad(Re_outer,Pr_outer,D1_outer,D2_inner);
            outer_f.hConv(i)=Nu_outer*k_outer_ave/D_h_outer;

            phi_friedel=1;
            all_vg=G_outer./all_ro_g;
            outer_f.vg(i)=all_vg(1);outer_f.vg(i+1)=all_vg(3);
                
            f_outer=Churchill(Re_outer,epsilon/D_h_outer);

            outer_f.P(i+1)=outer_f.P(i)- (1e-3)*(dz/A1_outer)*(   (phi_friedel*f_outer/4)*P2_inner * (G_outer^2/(2*all_ro_g(2))) ...
                                                + G_outer*A1_outer*(all_vg(3)-all_vg(1))/dz     )  ;
                                            
            outer_f.h(i+1)=outer_f.h(i)+ 1e-3*(  outer_f.hConv(i)*P2_inner*dz*(inner_w.T(i)-T_mean)/(G_outer*A1_outer)...
                                                 - 0.5*(all_vg(3)^2-all_vg(1)^2) );

            outer_f.T(i+1)=outer_f.T(i+1)*(outer_f.h(i+1)/previous_h)^2;   % A method for increasing convergence speed
            ERROR=abs((outer_f.h(i+1)-previous_h)/outer_f.h(i+1));   
        end   
    end

    %% inner wall simulation
    TEMP.inner_w_T=inner_w.T;
    thermal_cond_wall=55;  %w/(m c) for chrome steel ;
    wall_matrix=zeros(mesh,mesh);wall_vector=zeros(mesh,1);
    for i=1:mesh
        if i>1
            wall_matrix(i,i-1)=-thermal_cond_wall*A_inner_wall/dz;
        end
        wall_matrix(i,i)=2*thermal_cond_wall*A_inner_wall/dz+(inner_f.hConv(i)*P1_inner+outer_f.hConv(i)*P2_inner)*dz;
        if i<mesh
            wall_matrix(i,i+1)=-thermal_cond_wall*A_inner_wall/dz;
        end
        wall_vector(i)=(inner_f.hConv(i)*P1_inner*0.5*(inner_f.T(i)+inner_f.T(i+1))+outer_f.hConv(i)*P2_inner*0.5*(outer_f.T(i)+outer_f.T(i+1)))*dz;
    end
    inner_w.T=(wall_matrix\wall_vector)';
    
    %% Defining the criteria for convergence of all PDEs
    CRITERIA1=norm(abs(inner_f.h-TEMP.inner_f_h)./inner_f.h);
    CRITERIA2=norm(abs(outer_f.h-TEMP.outer_f_h)./outer_f.h);
    CRITERIA=CRITERIA1+CRITERIA2;
end
time=toc
plot((0:mesh)*dz,inner_f.T,'-s');
hold on
plot((0:mesh)*dz,outer_f.T,'r-o');
hold on
plot((0:mesh-1)*dz+dz/2,inner_w.T,'m-^');
hold on
legend('water Temperatuer','Steam Temperature','pipe wall Temperature')
xlabel('length(m)');ylabel('Temperature (c)');
title('Changing Temperature of superheated outer steam,subcooled inner water,pipe wall');
