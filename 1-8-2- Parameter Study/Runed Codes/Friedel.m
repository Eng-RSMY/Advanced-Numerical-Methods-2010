function phi= Friedel(G_total,D_h,x,alpha,ro_l,ro_g,mu_l,mu_g,IFT,f_g,f_l,g)
if nargin<12, g=9.8;end
ro_tp=ro_g*alpha +(1-alpha)*ro_l;
A1=(1-x)^2 + x^2 * (ro_l* f_g)/(ro_g*f_l);
A2=x^0.78 * (1-x)^0.224;
A3=(ro_l/ro_g)^0.91 *(mu_g/mu_l)^0.19 *(1-mu_g/mu_l)^.7;
Fr=G_total^2/(g*D_h*ro_tp);
We=G_total^2*D_h/(ro_tp*IFT);
phi2=A1+3.24*A2*A3/(Fr^.045*We^.035);
phi=phi2^.5;

