function y = EquCalc(X,varargin)
global outer_f inner_f
%% interpretting the coming varargin
[splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall]=varargin{:};

%% interpretting the comming X
mesh=floor(length(X)/5);

inner_f.P=X(1:2:2*mesh);
inner_f.T=X(2:2:2*mesh);

outer_f.T=X(2*mesh+1:2:4*mesh);
outer_f.x=X(2*mesh+2:2:4*mesh);

inner_w.T=X(4*mesh+1:5*mesh);

%% adding boundary condition
inner_f.P=[P_inner inner_f.P];
inner_f.T=[T_inner_fluid inner_f.T];

outer_f.T=[T_outer_fluid outer_f.T];
outer_f.x=[1 outer_f.x];

%% initial velocities
v_inner=G_inner/ro_inner;

%% making variables
inner_f.h=ppval(splined_T_hl_coeff,inner_f.T);
outer_f.h = outer_f.x .* ppval(splined_T_hg_coeff,outer_f.T)+ ...
          (1-outer_f.x).* ppval(splined_T_hl_coeff,outer_f.T);

%% function starter
inner_f.cp=Cp_inner(inner_f.T);
outer_f.cp=Cp_outer(outer_f.T,outer_f.x);
Re_inner=ro_inner*v_inner*D_h_inner/mu_inner;
Pr_inner=inner_f.cp*mu_inner/k_inner;
Nu_inner = Gnieliski(Re_inner,Pr_inner);
inner_f.hConv=Nu_inner*k_inner/D_h_inner;

%% simulating inner fulid equations
f_inner=Churchill(Re_inner,epsilon/D1_inner);
for i=1:mesh
    j= 2*i - 1;
    k= 2*i;
    
    y(j)=-inner_f.P(i+1)+inner_f.P(i)-1e-3*dz*P1_inner/A1_inner*((f_inner/4) * G_inner^2/(2*ro_inner));
    T_mean=0.5*(inner_f.T(i)+inner_f.T(i+1));
    Pr_inner=Cp_inner(T_mean)*mu_inner/k_inner;
    Nu=Gnieliski(Re_inner,Pr_inner);
    inner_f.hConv(i)=k_inner*Nu/D_h_inner;
    inner_f.h(i+1)=inner_f.h(i)+ 1e-3*(inner_f.hConv(i)*P1_inner*dz*(inner_w.T(i)-T_mean)/(G_inner*A1_inner));
    y(k)=-inner_f.T(i+1)+ppval(splined_hl_T_coeff,inner_f.h(i+1));
end

%% simulating outer fluid equations
for i=1:mesh
    j=2*mesh + 2*i -1;
    k=2*mesh + 2*i;
    
    T_mean=0.5*(outer_f.T(i)+outer_f.T(i+1));
    x_mean=0.5*(outer_f.x(i)+outer_f.x(i+1));
    
    outer_f.P(i)=ppval(splined_T_P_coeff,outer_f.T(i));

    all_x=[outer_f.x(i) x_mean outer_f.x(i+1)];
    all_ro_l=ppval(splined_T_vl_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
    all_ro_g=ppval(splined_T_vg_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
 

    hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
    hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

    outer_f.h(i+1)=outer_f.x(i+1)*hg_outer+(1-outer_f.x(i+1))*hl_outer;
    
    alpha = Smith(all_x,all_ro_l,all_ro_g);
    %mu_outer_ave=alpha(2)*mu_outer_g+(1-alpha(2))*mu_outer_l;
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

    ro_tp=alpha(2)*all_ro_g(2)+(1-alpha(2))*all_ro_l(2);
    outer_f.void(i)=alpha(2);

    all_kapa=all_x.*all_vg+(1-all_x).*all_vl;
    f_outer=Churchill(Re_outer,epsilon/D_h_outer);

    outer_f.P(i+1)=outer_f.P(i)-1e-3*(   dz*P2_inner/A1_outer*((phi_friedel*f_outer/4) * (G_outer^2/(2*ro_tp)) )...
                                        + G_outer*(all_kapa(3)-all_kapa(1)));
    outer_f.h(i+1)=outer_f.h(i)+ 1e-3*(  outer_f.hConv(i)*P2_inner*dz*(inner_w.T(i)-T_mean)/(G_outer*A1_outer)...
                                         - 0.5*(all_kapa(3)-all_kapa(1)) );

    y(j)=-outer_f.T(i+1)+ppval(splined_P_T_coeff,outer_f.P(i+1));

    hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
    hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

    y(k)=-outer_f.x(i+1)+(outer_f.h(i+1)-hl_outer)/(hg_outer-hl_outer);

end

%% simulating wall equations
thermal_cond_wall=55;
a=zeros(1,mesh);d=zeros(1,mesh);
for i=1:mesh
    j=4*mesh+i;
      
      a(i)=2*thermal_cond_wall*A_inner_wall/dz+(inner_f.hConv(i)*P1_inner+outer_f.hConv(i)*P2_inner)*dz;
      b=thermal_cond_wall*A_inner_wall/dz;
      c=thermal_cond_wall*A_inner_wall/dz;
      d(i)=(inner_f.hConv(i)*P1_inner*0.5*(inner_f.T(i)+inner_f.T(i+1))+outer_f.hConv(i)*P2_inner*0.5*(outer_f.T(i)+outer_f.T(i+1)))*dz;
      if i==1
        y(j)=-a(i)*inner_w.T(i)+b*inner_w.T(i+1)+d(i);
      elseif i==mesh
        y(j)=-a(i)*inner_w.T(i)+c*inner_w.T(i-1)+d(i);
      else
        y(j)=-a(i)*inner_w.T(i)+b*inner_w.T(i+1)+c*inner_w.T(i-1)+d(i);       
      end
      
end
%keyboard