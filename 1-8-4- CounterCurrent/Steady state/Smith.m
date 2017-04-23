function alpha = Smith(x,ro_l,ro_g)
K=0.4;
PI=ro_g./ro_l;
S=K+(1-K)*((1./PI)+K.*(1-x)./x) ./ (1+K.*(1-x)./x).^.5;
alpha=1./(1+(1./x-1).* S .* (ro_g./ro_l));