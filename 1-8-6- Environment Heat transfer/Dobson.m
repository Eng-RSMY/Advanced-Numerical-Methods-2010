function Nu = Dobson(Re_l,Pr_l,ro_g_l,mu_l_g,x)
X= ( (1-x)/x ) ^ .9 * ( ro_g_l )^.5 * ( mu_l_g )^.1;
Nu=.023 * Re_l^.8 * Pr_l ^.4 * (1+2.22/X^.89);
