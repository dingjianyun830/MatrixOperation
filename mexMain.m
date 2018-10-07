cd C:\Users\robin\Desktop\MartixMEX
mex FitLS.cpp
N=8;
x = [1:N];
for i =1:N
    y(i) = -2*x(i)*x(i)+2*x(i)+1;
end
m = 2;
[a,dt]=FitLS(x,y,m)


