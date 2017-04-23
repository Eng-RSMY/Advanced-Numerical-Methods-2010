function Nu=Gnieliski(Re,Pr)
f=(1.82*log10(Re)-1.64)^-2;
% Nu=.023*Re^.8*Pr.^.4   %this is for bolter.
Nu=( (f/8)*(Re-1000).*Pr )./ ( 1+12.7.*(f/8)^.5*(Pr.^(2/3)-1) );
