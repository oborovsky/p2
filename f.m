x = 0:0.01:1;
m = 0.6;
color = ['r','g','b','y','c'];
for i=1:5
    m = m - 0.1
    y = U(x,m);
    plot(x,y,color(i));
    hold on;
end
hold off;

function y=U(x,m)
%     y = (exp(-x/m) - 1)/(exp(-1/m) - 1);
%      y = (exp(-x/m) - exp(-1/m))/(1- exp(-1/m));
    y = (2*m^2 - m - 2/3)*(exp(-x/m) - 1)/(exp(-1/m)-1) + 1 - 2*m^2*x + m*x.^2 -x.^3/3;
%     C2 = (1 + 6*exp((sqrt(1+m)-1)/m)/(m+1) - 3*exp(1)/(m+1))/(exp(-(sqrt(m+1)+1)/m) - exp((sqrt(1+m)-1)/m));
%     C1 =  -C2 - 6/(m+1);
%     y = C1*exp(x*(sqrt(1+m)-1)/m) + C2*exp(-x*(sqrt(m+1)+1)/m) + 3*(2-x).*exp(x)/(m+1);
    
end