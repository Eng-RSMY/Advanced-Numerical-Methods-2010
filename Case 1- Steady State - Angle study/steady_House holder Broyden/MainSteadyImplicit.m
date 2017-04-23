clc;clear all;close all;tic
%% Data
mesh=20;
L=7;             %m
T_inner_fluid=22; %celcius
P_inner=10*10^2;  %kp

T_outer_fluid=180;   %

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

sigma=@(T)(-0.000000264568765*T.^2 - 0.000142361305361*T + 0.075698601398601);                     %T=degree c;interfacial tension N/m

%% hydraulic diameter-permimeter-area calculation
D_h_inner=D1_inner;
D_h_outer=(D1_outer^2-D2_inner^2)/D2_inner;

P1_inner=pi*D1_inner;
P2_inner=pi*D2_inner;
P1_outer=pi*D1_outer;

A1_inner=pi*D_h_inner^2/4;
A1_outer=pi*(D_h_outer)^2/4;
A_inner_wall=pi*(D2_inner^2-D1_inner^2)/4;

%% laoding spline data
load spline_hl_T_coeff splined_hl_T_coeff
load spline_T_P_coeff splined_T_P_coeff
load spline_P_T_coeff splined_P_T_coeff
load spline_T_vl_coeff splined_T_vl_coeff
load spline_T_vg_coeff splined_T_vg_coeff
load spline_T_hl_coeff splined_T_hl_coeff
load spline_T_hg_coeff splined_T_hg_coeff

global outer_f inner_f
%% initial velocities and entalpies
v_inner=G_inner/ro_inner;
v_outer=G_outer/ro_outer_g;

inner_initial_entalpy=ppval(splined_T_hl_coeff,T_inner_fluid); %92.207 for saturated liquid at 3kpa 
outer_initial_entalpy=ppval(splined_T_hg_coeff,T_outer_fluid);

%% making variables
dz=L/mesh;
inner_f.P=ones(1,mesh+1)*P_inner;inner_f.T=ones(1,mesh+1);inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;inner_f.hConv=ones(1,mesh);   

inner_w.T=zeros(1,mesh);

outer_f.P=ones(1,mesh+1)*ppval(splined_T_P_coeff,T_outer_fluid);outer_f.T=ones(1,mesh+1);outer_f.h=ones(1,mesh+1)*outer_initial_entalpy;
outer_f.x=ones(1,mesh+1);outer_f.void=zeros(1,mesh);outer_f.hConv=ones(1,mesh); 

outer_w.T=zeros(1,mesh);%temperature

%% starting the program
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

%% Defining the Varargin
varargin={splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
          splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
          L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
          G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
          ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall};



inner_temp_vector0=[inner_f.P(2:end);inner_f.T(2:end)];
outer_temp_vector0=[outer_f.T(2:end);outer_f.x(2:end)];
Vector0=[inner_temp_vector0(:)'  outer_temp_vector0(:)' inner_w.T ];

%% This is the core part of this m.file.

    AnswerVector=newtons_HHB('EquCalc',Vector0,1e-6,100,varargin{:});

%% getting the answer
inner_f.P=AnswerVector(1:2:2*mesh);
inner_f.T=AnswerVector(2:2:2*mesh);

outer_f.T=AnswerVector(2*mesh+1:2:4*mesh);
outer_f.x=AnswerVector(2*mesh+2:2:4*mesh);

inner_w.T=AnswerVector(4*mesh+1:5*mesh);

%% adding the data
inner_f.P=[P_inner inner_f.P];
inner_f.T=[T_inner_fluid inner_f.T];

outer_f.T=[T_outer_fluid outer_f.T];
outer_f.P=ppval(splined_T_P_coeff,outer_f.T);
outer_f.x=[1 outer_f.x];
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
[ax,h1,h2]=plotyy([0:mesh]*dz,inner_f.T,[0:mesh]*dz,inner_f.P);
xlabel('length of Pipe')
ylabel('Temperature (c)','color','r');
grid minor;
title('properties of inner fluid')
v=axis;
set(h2,'color','b','marker','s','markerfacecolor','b','markersize',5,'markeredgecolor','b')
set(h1,'color','r','marker','^','markerfacecolor','r','markersize',6,'markeredgecolor','r')
%text(v(2)*1.07,(v(4)-v(3))/2.3,'pressure (kpa)','rotation',90,'color','b')
%text(5*dz,inner_f.T(5)-5,'\leftarrow  temperature (c)','color','r');
text(5*dz,inner_f.T(5)+45,'pressure (kpa) \rightarrow','color','b')
legend('Temperature','Pressure')

figure('name','outer fluid properties')
[ax,h1,h2]=plotyy([0:mesh]*dz,outer_f.T,[0:mesh]*dz,outer_f.P);
xlabel('length of Pipe')
ylabel('Temperature (c)','color','r');
grid minor;
v=axis;
title('properties of outer fluid')
set(h2,'color','b','marker','s','markerfacecolor','b','markersize',5,'markeredgecolor','b')
set(h1,'color','r','marker','^','markerfacecolor','r','markersize',6,'markeredgecolor','r')
text(v(2)*1.08,(v(3)+v(4))/2,'pressure (kpa)','rotation',90,'color','b')
%text(5*dz,outer_f.T(10)+.5,'\leftarrow  temperature (c)','color','r');
%text(5*dz,outer_f.T(5)-3,'pressure (kpa) \rightarrow','color','b')
legend('Temperature','Pressure')
%get(hObject,'color')

figure('name','Pipe temperature')
plot([0:mesh-1]*dz,inner_w.T,'-ob','markerfacecolor','b','markersize',4.5)
title('plot of Pipe Temperature vs. pipe length');
xlabel('pipe length (m)');ylabel('Temperature (c)');
legend('Pipe Temperature');grid minor;

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
    Z(i,1:mesh)=inner_w.T;
end
for i=outer_fluid_mesh
    Z(i,1:mesh)=outer_f.T(1:end-1);
end
for i=outer_wall_mesh
    Z(i,1:mesh)=outer_f.T(1:end-1);
end
surf(X,Y,Z);title('Temperature Contour map');
view(0,90);xlim([0,(mesh-1)*dz])
colorbar;ylim([1,y_mesh]);xlabel('length of pipe (m)');
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