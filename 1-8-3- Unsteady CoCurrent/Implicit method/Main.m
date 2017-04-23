clc;clear all;close all;tic
%% Data
L=7;%m
mesh=20;

dt=0.3;
TimeStep=100;

%                         initial condition
ini_cond='uniform';
%ini_cond='pulse';
%PULSE=10;

%                        boundary condition
T_inner_fluid=10; %celcius
P_inner=10*10^2;      %kp

T_outer_fluid=180;   

D1_inner=8*10^-3;%m
D2_inner=10*10^-3;

D1_outer=14*10^-3;%m
D2_outer=16*10^-3;

G_inner=350; %kg/s/m2
G_outer=150;

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

Mu_L=@(T)( 0.000000032451574*T.^4  -0.000009061289916 *T.^3+  0.000984509345712 *T.^2 -0.055211010387099*T +  1.778370453327188)*10^-3;
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
inner_f.P=ones(1,mesh+1)*P_inner;inner_f.T=ones(1,mesh+1);inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;inner_f.hConv=ones(1,mesh);   

inner_w.T=zeros(1,mesh);

outer_f.P=ones(1,mesh+1)*ppval(splined_T_P_coeff,T_outer_fluid);outer_f.T=ones(1,mesh+1);outer_f.h=ones(1,mesh+1)*outer_initial_entalpy;
outer_f.x=ones(1,mesh+1);outer_f.void=zeros(1,mesh);outer_f.hConv=ones(1,mesh); 
outer_f.ro_l=zeros(1,mesh+1);outer_f.ro_g=ones(1,mesh+1)*ro_outer_g;outer_f.m=ones(1,mesh+1)*G_outer*A_annulus;
outer_f.vl=zeros(1,mesh+1);outer_f.vg=ones(1,mesh+1)*v_outer;outer_f.ro_tp=ones(1,mesh)*ro_outer_g;

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

%% Steady state condition for initial guess.

if strcmp(ini_cond,'pulse')
        G_inner=G_inner/PULSE;
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
                    Pr_inner=Cp_inner(T_mean)*Mu_L(inner_f.T(i+1))/k_inner;
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
        %% PULSE
        G_inner=G_inner*PULSE;
        outer_f.x=ones(1,mesh+1);
end
%% Defining the Varargin
varargin={dt,splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
          splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
          L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
          G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
          ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall,A_annulus...
          ,thermal_cond_wall,ro_wall,cp_wall};

TimeStruc.inner=inner_f;
TimeStruc.outer=outer_f;
TimeStruc.wall=inner_w;
Temp_Time{1}=TimeStruc;

%% Core part 
inner_temp_vector0=[inner_f.P(2:end);inner_f.T(2:end)];
outer_temp_vector0=[outer_f.m(2:end);outer_f.T(2:end);outer_f.x(2:end)];
Vector0=[inner_temp_vector0(:)'  outer_temp_vector0(:)' inner_w.T ];
for i=2:TimeStep  
    AnswerVector=newtons('EquCalc',Vector0,1e-3,100,varargin{:});
    TimeStruc.inner=inner_f;
    TimeStruc.outer=outer_f;
    TimeStruc.wall=inner_w;
    Temp_Time{i}=TimeStruc;
    disp(i)
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
inner_f.P=[P_inner inner_f.P];
inner_f.T=[T_inner_fluid inner_f.T];

outer_f.m=[G_outer*A_annulus outer_f.m];
outer_f.T=[T_outer_fluid outer_f.T];
outer_f.P=ppval(splined_T_P_coeff,outer_f.T);
outer_f.x=[1 outer_f.x];
toc

save Pulse_mesh20 Temp_Time

%%  Showing the results just for water Temperature
figure('name','Water countou map')
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
    text
end