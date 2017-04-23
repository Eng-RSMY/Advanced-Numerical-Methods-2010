function err=error_main(NN,hroot)
aa=her_POL(NN);
if NN<17, delta=100;
else delta=2000;
end
power = min(15,NN-1);
for count=1:NN

mini=500000;
index=0;
for f=-delta:delta
num=hroot(count)+f*10^(-power);              
erro=0;
for er=1:NN+1
         erro=erro+aa(er)*(num)^(NN+1-er);
end
    if abs(0-erro)<mini 
        mini=abs(0-erro);
        index=num;
    end;
end
    hroot(count)=index;
end
err=hroot;