function y = heaviside(x)

y = zeros(size(x));

y(x==0)=0.5;
y(x>0)= 1;
y(x<0)=0;

