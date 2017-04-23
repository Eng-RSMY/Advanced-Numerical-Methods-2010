function cp=Cp_inner(T,R)
if nargin<2,R=8.314;end
T=T+273.15; %celcius to kelvin
cp_R=8.712+1.25e-3*T-.18e-6*T.^2;
cp=R*cp_R/18.015;
