function cp=Cp_outer(T,x,R)
if nargin<3,R=8.314;end
T=T+273.15;
cp_water_R=8.712+1.25e-3.*T-.18e-6.*T.^2;
cp_steam_R=3.47+1.45e-3.*T+1.21e4.*T.^(-2);
cp=R*(x.*cp_steam_R+(1-x).*cp_water_R)/18.015;