function hermite_plot(n,min_x,max_x)

if nargin<2
    min_x=-2;
end
if nargin<3
    max_x=2;
end
x=min_x:.01:max_x;
M=plot(x,polyval(her_POL(n),x),'.-','color',rand(1,3));A=num2str(n);legend(M,[' N= ' A]);
hold on;
plot([min_x,max_x],[0,0],'-*r','linewidth',2);
grid on;