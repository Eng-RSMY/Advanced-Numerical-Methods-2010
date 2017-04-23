function f = Churchill(Re,e_D)
if nargin<2,e_D=0;end
temp=-2*log10(e_D/3.7+(7/Re)^.9);
f=temp^-2;