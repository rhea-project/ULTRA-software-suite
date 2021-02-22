function y=CSTairfoil(A,x)


N1 = 0.5;
N2 = 1;

C = ((x.^N1)).*(1-x).^N2;

%% create Bernstein polynomial

n = length(A);

for v = 0:n-1
    Sx(v+1,:) = nchoosek(n-1,v)*x.^v.*(1-x).^(n-1-v);
end

% yb = zeros(1,length(x));

% for i = 1:n
%     yb(1,:) = yb(1,:) + A(i).*Sx(i,:);
% end

B1 = A(1).*Sx(1,:);
B2 = A(2).*Sx(2,:);
B3 = A(3).*Sx(3,:);
B4 = A(4).*Sx(4,:);
B5 = A(5).*Sx(5,:);

yb = B1+B2+B3+B4+B5;

y = C.*yb;