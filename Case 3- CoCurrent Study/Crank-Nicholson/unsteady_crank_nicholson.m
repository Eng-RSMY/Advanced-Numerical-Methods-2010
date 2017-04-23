clc;clear all;close all;tic
%% Data
mesh=20;
L=5;%m
dt=0.15;
TimeStep=3;
T_inner_fluid=22; %celcius
P_inner=10*10^2;  %kp

T_outer_fluid=180;   %

D1_inner=7.9*10^-3;%m
D2_inner=9.5*10^-3;

D1_outer=14*10^-3;%m
D2_outer=16*10^-3;

G_inner=300; %kg/s/m2
G_outer=250;

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
A_annulus=(D1_outer^2-D2_inner^2)*pi/4;
%% laoding spline data
load spline_hl_T_coeff splined_hl_T_coeff
load spline_T_P_coeff splined_T_P_coeff
load spline_P_T_coeff splined_P_T_coeff
load spline_T_vl_coeff splined_T_vl_coeff
load spline_T_vg_coeff splined_T_vg_coeff
load spline_T_hl_coeff splined_T_hl_coeff
load spline_T_hg_coeff splined_T_hg_coeff

global outer_f TimeStruc inner_f inner_w PreTime2
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
outer_f.x=ones(1,mesh+1);outer_f.void=zeros(1,mesh);outer_f.hConv=ones(1,mesh);outer_f.m=ones(1,mesh+1)*G_outer*A_annulus;

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

%% assuming the initial condition from steady state 
%%inner pipe
inner_f.T=[22,26.0818747181113,30.0601044819556,33.9382572251068,37.7143088571204,41.3906158376963,44.9718033699756,48.4598147136582,51.8573150736165,55.1664147124909,58.3885883546634,61.5252554682912,64.5787825770147,67.5522685399878,70.4482307198701,73.2684510127428,76.0136072285466,78.6852964987835,81.2863361973816,83.8190196717551,86.2850089401161];
inner_f.h=[92.2707022798782,109.361655391093,126.001296381425,142.204511282342,157.984402431458,173.353162116275,188.322140896773,202.902071226984,217.103273828569,230.935705403399,244.409057804796,257.532817709575,270.316232326001,282.768228738290,294.897412113275,306.712142299319,318.220637730835,329.430988266149,340.351040404522,350.988360571979,361.349865197254];
inner_f.hConv=[66.9394800898939,66.9636565732146,66.9871871955091,67.0100792641601,67.0323363774762,67.0539818775620,67.0750376965749,67.0955197733785,67.1154444123752,67.1348242594134,67.1536698076814,67.1719941216160,67.1898150801444,67.2071509081781,67.2240158303910,67.2404185565854,67.2563671826063,67.2718760291557,67.2869613739170,67.3016360120300];
inner_f.P=[1000,999.740155660654,999.480311321308,999.220466981961,998.960622642615,998.700778303269,998.440933963922,998.181089624576,997.921245285230,997.661400945884,997.401556606537,997.141712267191,996.881867927845,996.622023588499,996.362179249152,996.102334909806,995.842490570460,995.582646231113,995.322801891767,995.062957552421,994.803113213075];

%%inner wall
inner_w.T=[175.317726045318,175.299907005022,175.316189112561,175.351301934116,175.397130712350,175.449468796310,175.505901807139,175.564885891479,175.625393986785,175.686697210918,175.748261016797,175.809702372385,175.870760120970,175.931245513365,175.991002243507,176.049885923577,176.107779536248,176.164621500381,176.220379418645,176.271198503969];
%%outer pipe
outer_f.T=[180,179.832765574476,179.721344137839,179.637728876890,179.570758284917,179.514792130801,179.466652473912,179.424352550946,179.386568443665,179.352374245229,179.321097423056,179.292234385791,179.265398478742,179.240286519166,179.216656491730,179.194312237631,179.173092667872,179.152863962165,179.133513786022,179.114946960070,179.097085970221];
outer_f.P=[1002.60000000000,998.749096975703,996.189990362899,994.272986563544,992.739720489643,991.459851078594,990.360021246211,989.394416943555,988.532535715189,987.753061355055,987.040521307678,986.383336790713,985.772622585165,985.201415159207,984.664158475877,984.156351192085,983.674298210882,983.214931078889,982.775674911624,982.354348773131,981.949176211319];
outer_f.h=[2777.20000000000,2764.27644315862,2751.70526106499,2739.46393899760,2727.55488595004,2715.95369499735,2704.65287849430,2693.64489200396,2682.92225956473,2672.47765543790,2662.30388451198,2652.39386551940,2642.74067269197,2633.33760520493,2624.17819403456,2615.25614761436,2606.56527560160,2598.09948065767,2589.85284612426,2581.81966422803,2573.98382836640];
outer_f.hConv=[1832.69984341286,1831.19357448881,1829.67152430923,1828.12656515592,1826.57861445929,1825.01984962137,1823.45333992131,1821.88191616412,1820.30819188252,1818.73457564696,1817.16327738017,1815.59632000612,1814.03555889193,1812.48269955097,1810.93930520533,1809.40679617127,1807.88645169439,1806.37942554486,1804.88676534829,1803.40842666318];
outer_f.x=[1,0.993657368383068,0.987468062047627,0.981431949243281,0.975554557043161,0.969825935242419,0.964243563712325,0.958804434213093,0.953505319366834,0.948342916070421,0.943313891493983,0.938414907011735,0.933642658558250,0.928993923628891,0.924465572972208,0.920054549278705,0.915757833314577,0.911572442669383,0.907495477304326,0.903524135071056,0.899650451918572];
outer_f.vg=[38.7720000000000,38.6677455232434,38.5210302262852,38.3559558568072,38.1824121234358,38.0049255536566,37.8261830488606,37.6478171098017,37.4708732652423,37.2960454618440,37.1238037401567,36.9544688681390,36.7882594473061,36.6253231369232,36.4657567486241,36.3096189215181,36.1569383327208,36.0077205673317,35.8619550571403,35.7196195557039,35.5814767128765];
outer_f.vl=[0,0.00142984728983053,0.00282474921078557,0.00418489034504459,0.00550909846544289,0.00679965814069606,0.00805716461726457,0.00928231488147990,0.0104758502376748,0.0116385257289396,0.0127711007973529,0.0138743345349122,0.0149489770148766,0.0159957589940080,0.0170153894854884,0.0180085607108856,0.0189759558415206,0.0199182493653469,0.0208360968921977,0.0217301317103518,0.0226051628283005]*1.05;
outer_f.ro_g=[5.15836170432271,5.13946367928679,5.12690370037829,5.11749441420374,5.10996819105821,5.10368548873059,5.09828635084221,5.09354596266069,5.08931463975940,5.08548777398199,5.08198943242232,5.07876278974642,5.07576423883581,5.07295960314595,5.07032161347969,5.06782817670088,5.06546115651528,5.06320549208224,5.06104854494618,5.05897961013852,5.05677064928321];
outer_f.ro_tp=[6.82597757143925,10.1529191926507,13.4085620087455,16.6056276614265,19.7121638300850,22.7485433335881,25.7151746297387,28.6129005566485,31.4427872556539,34.2060238060229,36.9038765987267,39.5376562088909,42.1086836759209,44.6182677030781,47.0677024205128,49.4582789577347,51.7912924659401,54.0680278670106,56.2897412135661,58.4591105398507];

%outer_f.x=ones(1,mesh+1);
%% Defining the Varargin
varargin={dt,splined_hl_T_coeff,splined_T_P_coeff,splined_P_T_coeff,splined_T_vl_coeff,...
          splined_T_vg_coeff,splined_T_hl_coeff,splined_T_hg_coeff,...
          L,dz,T_inner_fluid,P_inner,T_outer_fluid,D1_inner,D2_inner,D1_outer,D2_outer,...
          G_inner,G_outer,epsilon,R,k_inner,k_outer_l,k_outer_g,ro_inner,mu_inner,ro_outer_l,mu_outer_l,...
          ro_outer_g,mu_outer_g,D_h_inner,D_h_outer,P1_inner,P2_inner,P1_outer,A1_inner,A1_outer,A_inner_wall,A_annulus};

TimeStruc.inner=inner_f;
TimeStruc.outer=outer_f;
TimeStruc.wall=inner_w;
Temp_Time{1}=TimeStruc;

%PreTime2.hl_outer=ppval(splined_T_hl_coeff,outer_f.T(i+1));
%PreTime2.hg_outer=ppval(splined_T_hg_coeff,outer_f.T(i+1));
%% Core part of this m file 
inner_temp_vector0=[inner_f.P(2:end);inner_f.T(2:end)];
outer_temp_vector0=[outer_f.m(2:end);outer_f.T(2:end);outer_f.x(2:end)];
Vector0=[inner_temp_vector0(:)'  outer_temp_vector0(:)' inner_w.T ];
for i=2:TimeStep  
    AnswerVector=newtons('EquCalc',Vector0,1e-3,5,varargin{:});
    
    TimeStruc.inner=inner_f;
    TimeStruc.outer=outer_f;
    TimeStruc.wall=inner_w;
    Temp_Time{i}=TimeStruc;
    
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