clc;clear all;close all;
%% Data
L=10;%m
mesh=50;

dt=0.3;
TimeStep=300;

%                         initial condition
ini_cond='pulse';
PULSE=3;

%                        boundary condition
T_inner_fluid=22; %celcius
P_inner=10*10^2;      %kp


T_outer_fluid=180;   

D1_inner=8*10^-3;%m
D2_inner=10*10^-3;

D1_outer=14*10^-3;%m
D2_outer=16*10^-3;

G_inner=300; %kg/s/m2
G_outer=300;

epsilon=1.5*10^-6; %m
R=8.314;

disp('--------------------------------------')
disp(['    please wait for ' num2str(TimeStep) ' Steps'])
disp('--------------------------------------')

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

thermal_cond_wall=55;
ro_wall=7865;
cp_wall=0.46;

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
A_annulus=(D1_outer^2-D2_inner^2)*pi/4;
%% laoding spline data
load spline_hl_T_coeff splined_hl_T_coeff
load spline_T_P_coeff splined_T_P_coeff
load spline_P_T_coeff splined_P_T_coeff
load spline_T_vl_coeff splined_T_vl_coeff
load spline_T_vg_coeff splined_T_vg_coeff
load spline_T_hl_coeff splined_T_hl_coeff
load spline_T_hg_coeff splined_T_hg_coeff

global outer_f TimeStruc inner_f inner_w

%% initial velocities and entalpies
v_inner=G_inner/ro_inner;
v_outer=G_outer/ro_outer_g;

inner_initial_entalpy=ppval(splined_T_hl_coeff,T_inner_fluid); %92.207 for saturated liquid at 3kpa 
outer_initial_entalpy=ppval(splined_T_hg_coeff,T_outer_fluid);

%% making variables
dz=L/mesh;
inner_f.P=ones(1,mesh+1)*P_inner;inner_f.T=ones(1,mesh+1)*T_inner_fluid;inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;
inner_f.hConv=ones(1,mesh);   

inner_w.T=ones(1,mesh);

outer_f.P=ones(1,mesh+1)*ppval(splined_T_P_coeff,T_outer_fluid);outer_f.T=ones(1,mesh+1);outer_f.h=ones(1,mesh+1)*outer_initial_entalpy;
outer_f.x=ones(1,mesh+1);outer_f.void=zeros(1,mesh);outer_f.hConv=ones(1,mesh); 
outer_f.ro_l=ones(1,mesh+1)*ro_outer_l;outer_f.ro_g=ones(1,mesh+1)*ro_outer_g;outer_f.m=ones(1,mesh+1)*G_outer*A_annulus;
outer_f.vl=zeros(1,mesh+1);outer_f.vg=ones(1,mesh+1)*v_outer;outer_f.ro_tp=ones(1,mesh)*ro_outer_g;

outer_w.T=ones(1,mesh);%temperature

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

%% Steady state condition for initial guess

if strcmp(ini_cond,'pulse')
    G_inner=G_inner*PULSE;
    varargin={splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
              splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
              L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
              G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
              ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall,thermal_cond_wall...
              thermal_cond_wall};


    inner_f.T=fliplr(inner_f.T);inner_f.P=fliplr(inner_f.P);
    inner_temp_vector0=[inner_f.P(1:end-1);inner_f.T(1:end-1)];
    outer_temp_vector0=[outer_f.T(2:end);outer_f.x(2:end)];
    Vector0=[inner_temp_vector0(:)'  outer_temp_vector0(:)' fliplr(inner_w.T) ];

    AnswerVector=newtons('steady_EquCalc',Vector0,1e-6,100,varargin{:});

    inner_f.T=PULSE*inner_f.T;
    TimeStruc.inner=inner_f;
    TimeStruc.outer=outer_f;
    TimeStruc.wall=inner_w;
    G_inner=G_inner/PULSE;
end


%% Defining the Varargin
varargin={dt,splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
          splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
          L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
          G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
          ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall,A_annulus ...
          ,thermal_cond_wall,ro_wall,cp_wall};

TimeStruc.inner=inner_f;
TimeStruc.outer=outer_f;
TimeStruc.wall=inner_w;
Temp_Time{1}=TimeStruc;
tic
%% Core part 
inner_temp_vector0=[inner_f.P(1:end-1);inner_f.T(1:end-1)];
outer_temp_vector0=[outer_f.m(2:end);outer_f.T(2:end);outer_f.x(2:end)];
Vector0=[inner_temp_vector0(:)'  outer_temp_vector0(:)' inner_w.T ];
for i=1:TimeStep  
    AnswerVector=newtons('unsteady_EquCalc',Vector0,1e-6,100,varargin{:});
    TimeStruc.inner=inner_f;
    TimeStruc.outer=outer_f;
    TimeStruc.wall=inner_w;
    Temp_Time{i}=TimeStruc;
    disp(i)
    toc
    Vector0=AnswerVector;
end

%% getting the answer
inner_f.P=AnswerVector(1:2:2*mesh);
inner_f.T=AnswerVector(2:2:2*mesh);

outer_f.m=AnswerVector(2*mesh+1:3:5*mesh);
outer_f.T=AnswerVector(2*mesh+2:3:5*mesh);
outer_f.x=AnswerVector(2*mesh+3:3:5*mesh);

inner_w.T=AnswerVector(5*mesh+1:6*mesh);

%% adding the data
inner_f.P=[inner_f.P P_inner];
inner_f.T=[inner_f.T T_inner_fluid];

outer_f.m=[G_outer*A_annulus outer_f.m];
outer_f.T=[T_outer_fluid outer_f.T];
outer_f.x=[1 outer_f.x];
outer_f.P=ppval(splined_T_P_coeff,outer_f.T);


save TS300_M20_P5 Temp_Time 

%%  Showing the results
figure('name','Water countou map')
z=zeros(1,mesh);
for ii=1:TimeStep
    y_mesh=200;
    x=[0:mesh-1]*dz;
    y=[1:y_mesh];
    [X,Y]=meshgrid(x,y);
    inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
    inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
    outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
    outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
    for i=inner_fluid_mesh
        Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
    end
    for i=inner_wall_mesh
        Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
    end
    for i=outer_fluid_mesh
        Z(i,1:mesh)=Temp_Time{ii}.inner.T(1:end-1);
    end
    for i=outer_wall_mesh
        Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
    end
    surf(X,Y,Z);title('subcooled water temperature vs. length');
    xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
    view(0,90);xlabel('length of pipe (m)');ylabel('temperature(c)')
    colorbar 
    disp(ii)
    legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
    pause(dt);
end

figure('name','water T Animation')
for i=1:TimeStep
    A=rand(1,3);
    plot(1:mesh+1,Temp_Time{i}.inner.T,'color',A);
    xlabel('length (m)');ylabel('temperatur of water(c)');
    legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
    pause(dt);hold on;
    disp(i)
end