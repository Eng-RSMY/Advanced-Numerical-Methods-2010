%% another animation
load AngleStudy
mesh=20;clc;close all;
subplot(3,2,1)
outlet_collector=ones(1,19);
for i=1:19
    outlet_collector(i)=AngleTempTime{i}{200}.inner.P(end);
end
plot(-90:10:90,outlet_collector,'r','linewidth',2);hold on;
plot(-90:10:90,outlet_collector,'bo','markerfacecolor','b');
title('Plot of outlet Pressure of water vs. Angle');xlim([-90 90]);
xlabel('Angle (degree)');ylabel('Pressrue(kpa)');
grid on;
xlim([-100 100])
subplot(3,2,2)
outlet_collector=ones(1,19);
for i=1:19
    outlet_collector(i)=AngleTempTime{i}{200}.outer.P(end);
end
plot(-90:10:90,outlet_collector,'r','linewidth',2);hold on;
plot(-90:10:90,outlet_collector,'bo','markerfacecolor','b');
title('Plot of outlet Pressure of steam vs. Angle');xlim([-90 90]);
xlabel('Angle (degree)');ylabel('Pressrue(kpa)');
grid on;
xlim([-100 100])
subplot(3,2,3)
outlet_collector=ones(1,19);
for i=1:19
    outlet_collector(i)=AngleTempTime{i}{200}.inner.T(end);
end
plot(-90:10:90,outlet_collector,'r','linewidth',2);hold on;
plot(-90:10:90,outlet_collector,'bo','markerfacecolor','b');
title('Plot of outlet Temperature of water vs. Angle');xlim([-90 90]);
xlabel('Angle (degree)');ylabel('Temperature (c)');
grid on;
xlim([-100 100])
subplot(3,2,4)
for i=1:19
    outlet_collector(i)=AngleTempTime{i}{200}.outer.T(end);
end
plot(-90:10:90,outlet_collector,'r','linewidth',2);hold on;
plot(-90:10:90,outlet_collector,'bo','markerfacecolor','b');
title('Plot of outlet Temperature of Steam vs. Angle');xlim([-90 90]);
xlabel('Angle (degree)');ylabel('Temperature (c)');
grid on;
xlim([-100 100])
subplot(3,2,5)
for i=1:19
    outlet_collector(i)=1-AngleTempTime{i}{200}.outer.void(end);
end
plot(-90:10:90,outlet_collector,'r','linewidth',2);hold on;
plot(-90:10:90,outlet_collector,'bo','markerfacecolor','b');
title('Plot of outlet steam velocity vs. Angle');xlim([-90 90]);
xlabel('Angle (degree)');ylabel('Void');
grid;
xlim([-100 100])
subplot(3,2,6)
for i=1:19
    outlet_collector(i)=AngleTempTime{i}{200}.outer.vg(end);
end
plot(-90:10:90,outlet_collector,'r','linewidth',2);hold on;
plot(-90:10:90,outlet_collector,'bo','markerfacecolor','b');
title('Plot of steam velocity vs. Angle');xlim([-90 90]);
xlabel('Angle (degree)');ylabel('velocity(m/s)');
grid;
xlim([-100 100])