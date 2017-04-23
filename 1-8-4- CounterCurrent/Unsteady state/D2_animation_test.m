clc
figure('name','T wall Animation')
for i=1:5
A=rand(1,3);
plot(1:mesh+1,Temp_Time{i}.outer.x,'color',A);
pause(1);hold on;
end