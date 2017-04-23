clc
mesh=20;
figure('name','T wall Animation')
for i=2:300
A=rand(1,3);
plot(1:mesh+1,Temp_Time{i}.outer.m,'color',A);
pause(.1);hold on;
end