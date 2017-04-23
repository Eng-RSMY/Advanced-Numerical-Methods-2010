clc;clear all;close all;
%% Data
mesh=40;
L=5;%m
T_inner_fluid=22; %celcius
P_inner=10*10^2;      %kp

T_outer_fluid=180;   %

D1_inner=7.9*10^-3;%m
D2_inner=9.5*10^-3;

D1_outer=14*10^-3;%m
D2_outer=16*10^-3;

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

inner_initial_entalpy=ppval(splined_T_hl_coeff,T_inner_fluid); %92.207 for saturated liquid at 3kpa 
outer_initial_entalpy=ppval(splined_T_hg_coeff,T_outer_fluid);

sigma=@(T)(-0.000000264568765*T.^2 - 0.000142361305361*T + 0.075698601398601);                     %T=degree c;interfacial tension N/m

%% initial velocities
v_inner=G_inner/ro_inner;
v_outer=G_outer/ro_outer_g;

%% making variables
dz=L/mesh;
inner_f.P=ones(1,mesh+1)*P_inner;inner_f.T=ones(1,mesh+1);inner_f.h=ones(1,mesh+1)*inner_initial_entalpy;inner_f.hConv=ones(1,mesh);   

inner_w.T=zeros(2,mesh);

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
inner_w.T(1,:)=0.8*temp(1:mesh);
inner_w.T(2,:)=temp(1:mesh);

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
                            % this line has been assumed for constant mu.
       inner_f.P(i+1)=inner_f.P(i)-1e-3*dz*P1_inner/A1_inner*((f_inner/4) * G_inner^2/(2*ro_inner));
       %inner_f.T(i+1)=(1.)*inner_f.T(i); %the guess is assumed 1.015 of previous T.
       ERROR=1;
       while ERROR>1e-5
            temp=inner_f.h(i+1);
            T_mean=0.5*(inner_f.T(i)+inner_f.T(i+1));
            Pr_inner=Cp_inner(T_mean)*mu_inner/k_inner;
            Nu=Gnieliski(Re_inner,Pr_inner);
            inner_f.hConv(i)=k_inner*Nu/D_h_inner;
            inner_f.h(i+1)=inner_f.h(i)+ 1e-3*(inner_f.hConv(i)*P1_inner*dz*(inner_w.T(1,i)-T_mean)/(G_inner*A1_inner));
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
        outer_f.x(i+1)=outer_f.x(i)  - 0.001;   %this line must be modified.
        while ERROR>10^-4

            T_mean=0.5*(outer_f.T(i)+outer_f.T(i+1));
            x_mean=0.5*(outer_f.x(i)+outer_f.x(i+1));

            outer_f.P(i+1)=ppval(splined_T_P_coeff,outer_f.T(i+1));

            all_x=[outer_f.x(i) x_mean outer_f.x(i+1)];
            all_ro_l=ppval(splined_T_vl_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
            all_ro_g=ppval(splined_T_vg_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
            all_sigma=sigma([outer_f.T(i) T_mean outer_f.T(i+1)]);

            hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
            hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

            outer_f.h(i+1)=outer_f.x(i+1)*hg_outer+(1-outer_f.x(i+1))*hl_outer;
            previous_h=outer_f.h(i+1);
            previous_x=outer_f.x(i+1);

            if x_mean>0.5       %This condition must be set.
                alpha = Smith(all_x,all_ro_l,all_ro_g);
                mu_outer_ave=alpha(2)*mu_outer_g+(1-alpha(2))*mu_outer_l;
                mu_outer_ave=mu_outer_g; %attention
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
            else
                alpha = Smith(all_x,all_ro_l,all_ro_g);           
                Re_outer=G_outer*D_h_outer/mu_outer_g;
                outer_f.cp=Cp_inner(T_mean);
                Re_l=G_outer*D_h_outer*(1-all_x(2))/mu_outer_l;
                Re_g=G_outer*D_h_outer*all_x(2)/mu_outer_g;
                Pr_l=outer_f.cp*mu_outer_l/k_outer_l;

                Nu_outer = Dobson(Re_l,Pr_l,ro_outer_g/ro_outer_l,mu_outer_l/mu_outer_g,all_x(2));
                outer_f.hConv(i)=Nu_outer*k_outer_l/D_h_outer;

                f_l=Churchill(Re_l,epsilon/D_h_outer);
                f_g=Churchill(Re_g,epsilon/D_h_outer);
                phi_friedel=Friedel(G_outer,D_h_outer,all_x(2),alpha(2),all_ro_l(2),all_ro_g(2),mu_outer_l,mu_outer_g,all_sigma(2),f_g,f_l);

                all_vl=G_outer.*(1-all_x)./(all_ro_l.*(1-alpha));
                all_vg=G_outer.*all_x./(all_ro_g.*alpha);

                ro_tp=alpha(2)*all_ro_g(2)+(1-alpha(2))*all_ro_l(2);
                outer_f.void(i)=alpha(2);
            end

            all_kapa=all_x.*all_vg+(1-all_x).*all_vl;
            f_outer=Churchill(Re_outer,epsilon/D_h_outer);

            outer_f.P(i+1)=outer_f.P(i)-1e-3*(   dz*P2_inner/A1_outer*((phi_friedel*f_outer/4) * (G_outer^2/(2*ro_tp)) )...
                                                + G_outer*(all_kapa(3)-all_kapa(1)));
                                            
            outer_f.h(i+1)=outer_f.h(i)+ 1e-3*(  outer_f.hConv(i)*P2_inner*dz*(inner_w.T(2,i)-T_mean)/(G_outer*A1_outer)...
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
        
        
%         if i==1
%             wall_matrix(i,i-1)=-thermal_cond_wall*A_inner_wall/dz;
%         else if i<mesh
%                       wall_matrix(i,i)=2*thermal_cond_wall*A_inner_wall/dz+(inner_f.hConv(i)*P1_inner+outer_f.hConv(i)*P2_inner)*dz;
%             else if i==mesh
%                       wall_matrix(i,i+1)=-thermal_cond_wall*A_inner_wall/dz;
%                 end
%             end
%         end
%         wall_vector(i)=(inner_f.hConv(i)*P1_inner*0.5*(inner_f.T(i)+inner_f.T(i+1))+outer_f.hConv(i)*P2_inner*0.5*(outer_f.T(i)+outer_f.T(i+1)))*dz;
    end
    %inner_w.T=(wall_matrix\wall_vector)';
    wall_Temp=(wall_matrix\wall_vector)';
    inner_w.T(1,:)=wall_Temp(1:mesh);
    inner_w.T(2,:)=wall_Temp(mesh+1:2*mesh);
    
    %% Defining the criteria for convergence of all PDEs
    CRITERIA1=norm(abs(inner_f.h-TEMP.inner_f_h)./inner_f.h);
    CRITERIA2=norm(abs(outer_f.h-TEMP.outer_f_h)./outer_f.h);
    %CRITERIA3=norm(abs(inner_w.T-TEMP.inner_w_T)./inner_w.T);
    CRITERIA=CRITERIA1+CRITERIA2;%+CRITERIA3;
end

inner_f.hConv;
outer_f.hConv;
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

figure('name','Temperature contour')
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
for i=inner_wall_mesh(1:floor(length(inner_wall_mesh)))
    Z(i,1:mesh)=inner_w.T(1,1:mesh);
end
for i=inner_wall_mesh(ceil(length(inner_wall_mesh)):end)
    Z(i,1:mesh)=inner_w.T(2,1:mesh);
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
ylabel=(' radius (axissymmetric view)');

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
colorbar;ylabel=(' radius (axissymmetric view)');