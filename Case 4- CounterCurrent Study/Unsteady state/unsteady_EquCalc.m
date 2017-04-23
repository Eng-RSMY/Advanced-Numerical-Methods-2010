function y = unsteady_EquCalc(X,varargin)
global TimeStruc inner_f outer_f inner_w
%% interpretting the coming varargin
[dt,splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall,A_annulus,...
 thermal_cond_wall,ro_wall,cp_wall]=varargin{:};

%% interpretting the comming X
mesh=floor(length(X)/6);

inner_f.P=X(1:2:2*mesh);
inner_f.T=X(2:2:2*mesh);

outer_f.m=X(2*mesh+1:3:5*mesh);
outer_f.T=X(2*mesh+2:3:5*mesh);
outer_f.x=X(2*mesh+3:3:5*mesh);

inner_w.T=X(5*mesh+1:6*mesh);
%keyboard
%% adding initial condition
inner_f.P=[inner_f.P P_inner];
inner_f.T=[inner_f.T T_inner_fluid];

outer_f.m=[G_outer*A_annulus outer_f.m];
outer_f.T=[T_outer_fluid outer_f.T ];
outer_f.x=[1 outer_f.x];

%% assumed data
sigma=@(T)(-0.000000264568765*T.^2 - 0.000142361305361*T + 0.075698601398601);                     %T=degree c;interfacial tension N/m

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
m_dot=G_inner*A1_inner;
for i=1:mesh
    j= 2*i - 1;
    k= 2*i;
    landa=mesh+2-i;
    
    y(j)=-inner_f.P(landa-1)+inner_f.P(landa)-1e-3*dz*P1_inner/A1_inner*((f_inner/4) * G_inner^2/(2*ro_inner));
    T_mean=0.5*(inner_f.T(landa)+inner_f.T(landa-1));
    Pr_inner=Cp_inner(T_mean)*mu_inner/k_inner;
    Nu=Gnieliski(Re_inner,Pr_inner);
    inner_f.hConv(landa)=k_inner*Nu/D_h_inner;
    c_coeff1=(0.5*(inner_f.P(landa)+inner_f.P(landa-1)))- 0.5*(TimeStruc.inner.P(landa)+TimeStruc.inner.P(landa-1));
    c_coeff2=inner_f.h(landa)- (TimeStruc.inner.h(landa)+TimeStruc.inner.h(landa-1));
    c_coeff=2*c_coeff1-ro_inner*c_coeff2;

    numerator=2*1e-3*inner_f.hConv(landa)*dz*P1_inner*(inner_w.T(landa-1)-T_mean) + 2 * m_dot*inner_f.h(landa) + c_coeff*A1_inner*dz/(dt);
    inner_f.h(landa-1)=numerator/(2*m_dot+ro_inner*A1_inner*dz/(dt));
    y(k)=-inner_f.T(landa-1)+ppval(splined_hl_T_coeff,inner_f.h(landa-1));
end

%% simulating outer fluid equations

for i=1:mesh
    j=2*mesh + 3*i -2;
    k=2*mesh + 3*i-1;
    kk=2*mesh + 3*i;
    
    T_mean=0.5*(outer_f.T(i)+outer_f.T(i+1));
    x_mean=0.5*(outer_f.x(i)+outer_f.x(i+1));
    
    outer_f.P(i)=ppval(splined_T_P_coeff,outer_f.T(i));

    all_x=[outer_f.x(i) x_mean outer_f.x(i+1)];
    all_ro_l=ppval(splined_T_vl_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
    all_ro_g=ppval(splined_T_vg_coeff,[outer_f.T(i) T_mean outer_f.T(i+1)]).^-1;
    outer_f.ro_g(i)=all_ro_g(1);outer_f.ro_g(i+1)=all_ro_g(3);
    all_sigma=sigma([outer_f.T(i) T_mean outer_f.T(i+1)]);

    hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
    hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

    outer_f.h(i+1)=outer_f.x(i+1)*hg_outer+(1-outer_f.x(i+1))*hl_outer;
    
  
    alpha = Smith(all_x,all_ro_l,all_ro_g);
    mu_outer_ave=alpha(2)*mu_outer_g+(1-alpha(2))*mu_outer_l;
    %mu_outer_ave=mu_outer_g;
    Re_outer=G_outer*D_h_outer/mu_outer_ave;
    outer_f.cp=Cp_outer(T_mean,x_mean);
    k_outer_ave=x_mean*k_outer_g+(1-x_mean)*k_outer_l;
    Pr_outer=outer_f.cp*mu_outer_ave/k_outer_ave;
    Nu_outer = Monrad(Re_outer,Pr_outer,D1_outer,D2_inner);
    outer_f.hConv(i)=Nu_outer*k_outer_ave/D_h_outer;
    for kij=1:3
        if alpha(kij)==1,alpha(kij)=.999;end
    end
    phi_friedel=1;
    G_l_outer=(1-all_x).*G_outer;
    all_vl=G_l_outer./(all_ro_l.*(1-alpha)); %this was changed.
    outer_f.vl(i)=all_vl(1);outer_f.vl(i+1)=all_vl(3);

    G_g_outer=all_x.*G_outer;
    all_vg=G_g_outer./(all_ro_g.*(alpha));   %this was changed.
    outer_f.vg(i)=all_vg(1);outer_f.vg(i+1)=all_vg(3);

    ro_tp=alpha(2)*all_ro_g(2)+(1-alpha(2))*all_ro_l(2);
    outer_f.ro_tp(i)=ro_tp;
    outer_f.void(i)=alpha(2);
    
    y(j)=-outer_f.m(i+1)+outer_f.m(i)-(A_annulus*dz/dt)*(ro_tp-TimeStruc.outer.ro_tp(i));
    
    m_dot_mean=0.5*(outer_f.m(i)+outer_f.m(i+1));
    all_m_dot=[outer_f.m(i) m_dot_mean outer_f.m(i+1)];
    all_kapa=all_m_dot.*(all_x.*all_vg+(1-all_x).*all_vl);
    f_outer=Churchill(Re_outer,epsilon/D_h_outer);

    outer_f.P(i+1)=outer_f.P(i)- 1e-3*dz/A_annulus*(phi_friedel*f_outer*m_dot_mean^2*P2_inner/(8*ro_tp*A_annulus^2)+...
                                                     (m_dot_mean-0.5*(TimeStruc.outer.m(i)+TimeStruc.outer.m(i+1)))/dt+ ...
                                                     (all_kapa(3)-all_kapa(1))/dz);
                                    
                                    
    ab_coeff=all_kapa.^2./all_m_dot;    %ab_coeff represent the second and the third praremeters in energy EQ.                                 
    c_coeff1=(0.5*(outer_f.P(i)+outer_f.P(i+1)))- 0.5*(TimeStruc.outer.P(i)+TimeStruc.outer.P(i+1));
    c_coeff2=outer_f.h(i)- (TimeStruc.outer.h(i)+TimeStruc.outer.h(i+1));
    c_coeff3=0.125*(outer_f.ro_g(i)+outer_f.ro_g(i+1))*(outer_f.vg(i)+outer_f.vg(i+1))^2 ...
                                      -0.125*(TimeStruc.outer.ro_g(i)+TimeStruc.outer.ro_g(i+1))*...
                                             (TimeStruc.outer.vg(i)+TimeStruc.outer.vg(i+1))^2;%we are not sure about this term.

    c_coeff=2*c_coeff1-TimeStruc.outer.ro_tp(i)*c_coeff2-c_coeff3*1e-3;

    numerator=2*1e-3*outer_f.hConv(i)*dz*P2_inner*(inner_w.T(i)-T_mean) - ab_coeff(3)+ ab_coeff(1)...
        +outer_f.m(i)*outer_f.h(i)+outer_f.m(i+1)*outer_f.h(i)+ c_coeff*A_annulus*dz/dt; 
    outer_f.h(i+1)=(numerator)/(all_m_dot(3)+all_m_dot(1)+ TimeStruc.outer.ro_tp(i)*A_annulus*dz/dt);

    y(k)=-outer_f.T(i+1)+ppval(splined_P_T_coeff,outer_f.P(i+1));

    hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
    hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));

    y(kk)=-outer_f.x(i+1)+(outer_f.h(i+1)-hl_outer)/(hg_outer-hl_outer);

end

%% simulating wall equations
a=zeros(1,mesh);d=zeros(1,mesh);
for i=1:mesh
    j=5*mesh+i;
      
      a(i)=2*thermal_cond_wall*A_inner_wall/dz+(inner_f.hConv(i)*P1_inner+outer_f.hConv(i)*P2_inner)*dz+A_inner_wall*dz/dt*ro_wall*cp_wall;
      b=thermal_cond_wall*A_inner_wall/dz;
      c=thermal_cond_wall*A_inner_wall/dz;
      d(i)=(inner_f.hConv(i)*P1_inner*0.5*(inner_f.T(i)+inner_f.T(i+1))+outer_f.hConv(i)*P2_inner*0.5*(outer_f.T(i)+outer_f.T(i+1)))*dz ...
                        +A_inner_wall*dz/dt*ro_wall*cp_wall*TimeStruc.wall.T(i);
      if i==1
        y(j)=-a(i)*inner_w.T(i)+b*inner_w.T(i+1)+d(i);
      elseif i==mesh
        y(j)=-a(i)*inner_w.T(i)+c*inner_w.T(i-1)+d(i);
      else
        y(j)=-a(i)*inner_w.T(i)+b*inner_w.T(i+1)+c*inner_w.T(i-1)+d(i);       
      end
      
end
