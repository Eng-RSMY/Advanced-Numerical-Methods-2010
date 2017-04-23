%% angle animation
load AngleStudy
mesh=20;clc;close all;
B='g-or-db-*k->c:dy-+b-sr-xg-or->b-^k:<y:dg-or-db-*k->c:dm-+bsr-xg-or->b-^k:<y:d';

flag=1;
while flag==1
hold off; 
counter=1;
colorIter=1;   
for i=-90:10:90

subplot(2,2,[2 4])
u=linspace(0,2*pi,41);
v=linspace(-2,2,41);
[U V]=meshgrid(u,v);
% vertical cylinder of radius 1
surf(7.9*cos(U),V,7.9*sin(U));hold on
surf(9*cos(U),V,9*sin(U))
surf(14.5*cos(U),V,14.5*sin(U))
surf(16*cos(U),V,16*sin(U))
colormap copper
view(0,i)
title('Effect of changing Angle on SteadyState condition');
subplot(2,2,[1,3])
plot(1:mesh+1,AngleTempTime{counter}{200}.outer.T,B(colorIter:colorIter+2));
hold on;
legend('-90','-80','-70','-60','-50','-40','-30','-20','-10','0','10','20','30','40','50','60','70','80','90','location','northeastoutside');
xlabel('length(m)');ylabel('temperature of water (c)');
colorIter=colorIter+3;
clc
pause(.5)
counter=counter+1;

end
end


