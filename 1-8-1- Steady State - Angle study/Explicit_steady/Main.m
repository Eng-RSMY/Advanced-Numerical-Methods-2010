clc;clear all;close all;tic
%% laoding spline data
load spline_hl_T_coeff splined_hl_T_coeff
load spline_T_P_coeff splined_T_P_coeff
load spline_P_T_coeff splined_P_T_coeff
load spline_T_vl_coeff splined_T_vl_coeff
load spline_T_vg_coeff splined_T_vg_coeff
load spline_T_hl_coeff splined_T_hl_coeff
load spline_T_hg_coeff splined_T_hg_coeff

%% Data
mesh=20;
L=7;%m

T_inner_fluid=22; %celcius
P_inner=10*10^2;      %kp

T_outer_fluid=ppval(splined_P_T_coeff,P_inner);   

D1_inner=8*10^-3;%m
D2_inner=10*10^-3;
D1_outer=14*10^-3;%m
D2_outer=16*10^-3;

G_inner=250; %kg/s/m2
G_outer=80;

epsilon=1.5*10^-6; %m
R=8.314;

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

%Mu_L=@(T)( 0.000000032451574*T.^4  -0.000009061289916 *T.^3+  0.000984509345712 *T.^2 -0.055211010387099*T +  1.778370453327188)*10^-3;
Mu_L=@(T)( 0.8 )*10^-3;
%% hydraulic diameter-permimeter-area-velocity calculation
D_h_inner=D1_inner;
D_h_outer=(D1_outer^2-D2_inner^2)/D2_inner;

P1_inner=pi*D1_inner;
P2_inner=pi*D2_inner;
P1_outer=pi*D1_outer;

A1_inner=pi*D_h_inner^2/4;
A1_outer=pi*(D_h_outer)^2/4;
A_inner_wall=pi*(D2_inner^2-D1_inner^2)/4;

v_inner=G_inner/ro_inner;
v_outer=G_outer/ro_outer_g;

%% making variables
inner_initial_entalpy=ppval(splined_T_hl_coeff,T_inner_fluid);  
outer_initial_entalpy=ppval(splined_T_hg_coeff,T_outer_fluid);

dz=L/mesh;
inner_f.P=ones(1,mesh+1)*P_inner;inner_f.T=ones(1,mesh+1);inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;inner_f.hConv=ones(1,mesh);   
inner_f.vl=ones(1,mesh+1)*v_inner;

inner_w.T=zeros(1,mesh);

outer_f.P=ones(1,mesh+1)*ppval(splined_T_P_coeff,T_outer_fluid);outer_f.T=ones(1,mesh+1);outer_f.h=ones(1,mesh+1)*outer_initial_entalpy;
outer_f.x=ones(1,mesh+1);outer_f.void=zeros(1,mesh);outer_f.hConv=ones(1,mesh); 
outer_f.vg=zeros(1,mesh+1);outer_f.vl=zeros(1,mesh+1);outer_f.ro_g=zeros(1,mesh+1);outer_f.ro_tp=ones(1,mesh);
outer_w.T=zeros(1,mesh);%temperature

%% Wall temperature prediction
inner_f.T=T_inner_fluid*ones(1,mesh+1);
outer_f.T=T_outer_fluid*ones(1,mesh+1);

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
inner_w.T=temp(1:mesh);

outer_w.T=outer_f.T(1:mesh);

%% Solving the PDEs of Heat Exchanger (successive substitution method)
CRITERIA=1;
TEMP.inner_f_h=zeros(1,mesh+1);
TEMP.outer_f_h=zeros(1,mesh+1);
TEMP.inner_w_T=zeros(1,mesh);

while CRITERIA > 10^-5
    
    %% inner pipe simulation
    TEMP.inner_f_h=inner_f.h;
    f_inner=Churchill(Re_inner,epsilon/D1_inner);   %this line is for assumption.
    for i=1:mesh
                            % this line has been assumed for constant mu.
       inner_f.P(i+1)=inner_f.P(i)-1e-3*dz*P1_inner/A1_inner*((f_inner/4) * G_inner^2/(2*ro_inner));
       ERROR=1;
       while ERROR>1e-5
            temp=inner_f.h(i+1);
            T_mean=0.5*(inner_f.T(i)+inner_f.T(i+1));
            Pr_inner=Cp_inner(T_mean)*Mu_L(inner_f.T(i))/k_inner;
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
            outer_f.ro_g(i)=all_ro_g(1);outer_f.ro_g(i+1)=all_ro_g(3);
            hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
            hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

            outer_f.h(i+1)=outer_f.x(i+1)*hg_outer+(1-outer_f.x(i+1))*hl_outer;
            previous_h=outer_f.h(i+1);
            previous_x=outer_f.x(i+1);

            alpha = Smith(all_x,all_ro_l,all_ro_g);
            mu_outer_ave=mu_outer_g;
            Re_outer=G_outer*x_mean*D_h_outer/mu_outer_ave;
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
            outer_f.vg(i)=all_vg(1);outer_f.vg(i+1)=all_vg(3);
            outer_f.vl(i)=all_vl(1);outer_f.vl(i+1)=all_vl(3);
            ro_tp=alpha(2)*all_ro_g(2)+(1-alpha(2))*all_ro_l(2);
              outer_f.ro_tp(i)=ro_tp;
            outer_f.void(i)=alpha(2);

            all_kapa=all_x.*all_vg+(1-all_x).*all_vl;
            f_outer=Churchill(Re_outer,epsilon/D_h_outer);

            outer_f.P(i+1)=outer_f.P(i)-1e-3*(   dz*P2_inner/A1_outer*((phi_friedel*f_outer/4) * (G_outer^2/(2*ro_tp)) )...
                                                + G_outer*(all_kapa(3)-all_kapa(1)));
            outer_f.h(i+1)=outer_f.h(i)+ 1e-3*(  outer_f.hConv(i)*P2_inner*dz*(inner_w.T(i)-T_mean)/(G_outer*A1_outer)...
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
toc
%%  Showing the results

figure('name','showing the liquid film')
hold on;
fill([[0:mesh-1]*dz,[mesh-1:-1:0]*dz],[ones(1,mesh)*D1_inner/2*1e3,zeros(1,mesh)],'o-b','markersize',3,'markerfacecolor','b')
fill([[0:mesh-1]*dz,[mesh-1:-1:0]*dz],[ones(1,mesh)*D2_inner/2*1e3,ones(1,mesh)*D1_inner/2*1e3],'o-k','markersize',3,'markerfacecolor','k')
fill([[0:mesh-1]*dz,[mesh-1:-1:0]*dz],[(1-outer_f.void)*(D1_outer-D2_inner)/2*1e3+D2_inner/2*1e3,ones(1,mesh)*D2_inner/2*1e3],'o-c','markersize',3,'markerfacecolor','c')
fill([[0:mesh-1]*dz,[mesh-1:-1:0]*dz],[(1-outer_f.void)*(D1_outer-D2_inner)/2*1e3+D2_inner/2*1e3,ones(1,mesh)*D1_outer/2*1e3],'o-r','markersize',3,'markerfacecolor','r')
fill([[0:mesh-1]*dz,[mesh-1:-1:0]*dz],[ones(1,mesh)*D2_outer/2*1e3,ones(1,mesh)*D1_outer/2*1e3],'o-k','markersize',3,'markerfacecolor','k')
xlabel('lenght of pipe (m)');ylabel('Radius of heat exchagner * 10^3');
xlim([0,(mesh-1)*dz]);
title(' (axissymmetric view) Plot of water film  vs. length of pipe');
legend('Inner fluid(Water)','Inner pipe(chrome steel)','condensed water','Outer fluid(steam)','outer pipe(chrome steel)','location','southwest');

figure('name','inner fluid properties')
[ax,h1,h2]=plotyy([0:mesh]*dz,inner_f.T,[0:mesh]*dz,inner_f.P,'plot');
xlabel('length of Pipe')
ylabel('Temperature (c)','color','r');
grid on;
title('properties of inner fluid')
v=axis;
set(h2,'color','b','marker','s','markerfacecolor','b','markersize',5,'markeredgecolor','b')
set(h1,'color','r','marker','^','markerfacecolor','r','markersize',6,'markeredgecolor','r')
%set(get(ax(1),'Ylabel'),'String','')
set(get(ax(2),'Ylabel'),'String','Pressure (kpa)')
%text(v(2)*1.07,(v(4)-v(3))/2.3,'pressure (kpa)','rotation',90,'color','b')
%text(5*dz,inner_f.T(5)-5,'\leftarrow  temperature (c)','color','r');
text(5*dz,inner_f.T(5)+45,'pressure (kpa) \rightarrow','color','b')
legend('Temperature','Pressure')

figure('name','outer fluid properties')
[ax,h1,h2]=plotyy([0:mesh]*dz,outer_f.T,[0:mesh]*dz,outer_f.P);
xlabel('length of Pipe')
ylabel('Temperature (c)','color','r');
grid on;
v=axis;
title('properties of outer fluid')
set(h2,'color','b','marker','s','markerfacecolor','b','markersize',5,'markeredgecolor','b')
set(h1,'color','r','marker','^','markerfacecolor','r','markersize',6,'markeredgecolor','r')
set(get(ax(2),'Ylabel'),'String','Pressure (kpa)')
%text(v(2)*1.08,(v(3)+v(4))/2,'pressure (kpa)','rotation',90,'color','b')
%text(5*dz,outer_f.T(10)+.5,'\leftarrow  temperature (c)','color','r');
%text(5*dz,outer_f.T(5)-3,'pressure (kpa) \rightarrow','color','b')
legend('Temperature','Pressure')
%get(hObject,'color')

figure('name','Pipe temperature')
plot([0:mesh-1]*dz,inner_w.T(1:end),'-ob','markerfacecolor','b','markersize',4.5)
title('plot of Pipe Temperature vs. pipe length');
xlabel('pipe length (m)');ylabel('Temperature (c)');
legend('Pipe Temperature');grid on;

figure('name','Pressure contour')
y_mesh=200;
x=[0:mesh-1]*dz;
y=[1:y_mesh];
[X,Y]=meshgrid(x,y);
inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
for i=inner_fluid_mesh
    Z(i,1:mesh)=inner_f.T(1:end-1);
end
for i=inner_wall_mesh
    Z(i,1:mesh)=mean(inner_f.T);
end
for i=outer_fluid_mesh
    Z(i,1:mesh)=mean(inner_f.T(1:end-1));
end
for i=outer_wall_mesh
    Z(i,1:mesh)=inner_f.T(1:end-1);
end
surf(X,Y,Z);title('Temperature Contour map');
view(0,90);xlim([0,(mesh-1)*dz])
colorbar;ylim([1,floor(D1_inner*y_mesh/(D2_outer))]);xlabel('length of pipe (m)');
ylabel(' radius (axissymmetric view)');

figure('name','Pressure contour map')
y_mesh=200;
x=[0:mesh-1]*dz;
y=[1:y_mesh];
[X,Y]=meshgrid(x,y);
inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
for i=inner_fluid_mesh
    Z(i,1:mesh)=inner_f.P(1:end-1);
end
for i=inner_wall_mesh
    Z(i,1:mesh)=P_inner+2;
end
for i=outer_fluid_mesh
    Z(i,1:mesh)=outer_f.P(1:end-1);
end
for i=outer_wall_mesh
    Z(i,1:mesh)=P_inner+2;
end
surf(X,Y,Z);title('Pressure Contour map');
xlim([0,(mesh-1)*dz]);ylim([1,y_mesh]);
view(0,90);xlabel('length of pipe (m)');
colorbar;ylabel(' radius (axissymmetric view)');

figure('name','water velcity map')
y_mesh=200;
x=[0:mesh-1]*dz;
y=[1:y_mesh];
[X,Y]=meshgrid(x,y);
inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
for i=inner_fluid_mesh
    Z(i,1:mesh)=mean(outer_f.vl(1:end-1));%inner_f.vl(1:end-1);
end
for i=inner_wall_mesh
    Z(i,1:mesh)=mean(outer_f.vl(1:end-1));
end
for i=outer_fluid_mesh
    Z(i,1:mesh)=outer_f.vl(1:end-1);
end
for i=outer_wall_mesh
    Z(i,1:mesh)=mean(outer_f.vl(1:end-1));
end
surf(X,Y,Z);title('water velocity map');
xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
view(0,90);xlabel('length of pipe (m)');ylabel('radius*10^3')
colorbar

figure('name','steam velocity map')
y_mesh=200;
x=[0:mesh-1]*dz;
y=[1:y_mesh];
[X,Y]=meshgrid(x,y);
inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
for i=inner_fluid_mesh
    Z(i,1:mesh)=mean(outer_f.vg(1:end-1));%inner_f.vl(1:end-1);
end
for i=inner_wall_mesh
    Z(i,1:mesh)=mean(outer_f.vg(1:end-1));
end
for i=outer_fluid_mesh
    Z(i,1:mesh)=outer_f.vg(1:end-1);
end
for i=outer_wall_mesh
    Z(i,1:mesh)=mean(outer_f.vg(1:end-1));
end
surf(X,Y,Z);title('steam velocity map');
xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
view(0,90);xlabel('length of pipe (m)');ylabel('radius*10^3')
colorbar

figure('name','Hconv')
subplot(2,1,1)
plot([0:mesh-1]*dz,inner_f.hConv)
xlabel('length (m)');ylabel('hConv of subcooled water (J/m2/c)');
title('Plot of hConv vs. lenght of pipe')
subplot(2,1,2)
%figure('name','Hconv of steam')
plot([0:mesh-1]*dz,outer_f.hConv)
title('Plot of hConv vs. lenght of pipe')
xlabel('length (m)');ylabel('hConv of steam(J/m2/c)');

figure('name','velocity')
plotyy([0:mesh]*dz,inner_f.vl,[0:mesh]*dz,outer_f.vg,'plot')
