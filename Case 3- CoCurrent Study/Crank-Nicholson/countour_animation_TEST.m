close all;clc;
figure('name','Pressure contour map')
for ii=1:3
y_mesh=200;
x=[0:mesh-1]*dz;
y=[1:y_mesh];
[X,Y]=meshgrid(x,y);
inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
for i=inner_fluid_mesh
    Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));%inner_f.vl(1:end-1);
end
for i=inner_wall_mesh
    Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
end
for i=outer_fluid_mesh
    Z(i,1:mesh)=Temp_Time{ii}.inner.T(1:end-1);
end
for i=outer_wall_mesh
    Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
end
surf(X,Y,Z);title('steam velocity map');
xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
view(0,90);xlabel('length of pipe (m)');
colorbar
pause(1);%hold off;
end

