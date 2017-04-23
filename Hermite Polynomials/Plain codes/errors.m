function err=error_(NN,root)
aa=her_POL(NN);
% root=hermite_roots(NN);
for f=1:NN

erro=0;
for er=1:NN+1
         erro=erro+aa(er)*(root(f))^(NN+1-er);
end
  error(f)=0-erro;
end
err=error;