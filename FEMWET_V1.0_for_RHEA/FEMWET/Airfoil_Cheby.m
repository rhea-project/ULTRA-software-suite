function [X, Y] = Airfoil_Cheby(A,X0,Y0,ud)

% for upper surface X = 1..0  and for lower surface X = 0..1

Po = 0;
% 
n = Chebyshev2(A,X0)';
% 
X2 = meshgrid(X0,1:size(Y0,2))';
% 
% dX = (X2(2:end,:)-X2(1:end-1,:));
% if dX > 0
%     t = atan((Y0(2:end,:)-Y0(1:end-1,:))./(X2(2:end,:)-X2(1:end-1,:)));
% else
%     t = myatan2((Y0(2:end,:)-Y0(1:end-1,:)),(X2(2:end,:)-X2(1:end-1,:)));
%     t = reshape(t,size(X2,1)-1,size(X2,2));
% end
% 
% 
% t(end+1,:) = t(end,:);
% 
% t2 = t-pi/2;
% 
% dY = n.*abs(sin(t2));
% dX = n.*abs(cos(t2));
% 
% X = X2 + dX;
% Y = Y0 + dY;
% 
% % X = X0;
% % Y = Y0 + n;


t = atan((Y0(2:end,:)-Y0(1:end-1,:))./(X2(2:end,:)-X2(1:end-1,:)));
t(end+1,:) = t(end,:);

% n = Chebyshev(A,Snorm);

dY = n.*cos(t);
dX = n.*sin(t);

if ud ==1
%     X = Xs - dX;
%     Y = interp1(X0,Y0,Xs) + dY;
    X = X2 -dX;
    Y = Y0 + dY;
elseif ud ==2
%     X = Xs + dX;
%     Y = interp1(X0,Y0,Xs) - dY;
    X = X2 + dX;
    Y = Y0 -dY;
end


if Po ==1
    figure
    hold on
    plot(X0,Y0,'-b');
    plot(X0,n,'--k');
    plot(X,Y,'-r');
    hold off
end

end


function y = Chebyshev2(A,S)


SWT = 2;
 
SW = SWT*S ./ (1.0 + (SWT-1.0)*S);

X = 1.0 - 2.0*SW;

% X = min( 1.0 , X );
% X = max(-1.0 , X );

X = X.*heaviside(1-X+1e-9) + ones(size(X)).*heaviside(X-1-1e-9);
X = X.*heaviside(X+1+1e-9) - ones(size(X)).*heaviside(-1-X-1e-9);


THETA = acos(X);

y = zeros(size(A,1),length(S));
for i=1:size(A,2)
    RF = i+1;
    if mod(i,2) ==0
        y = y + meshgrid(A(:,i),1:length(S))'.*meshgrid((X-cos(RF*THETA))/RF,1:size(A,1));
    else
        y = y + meshgrid(A(:,i),1:length(S))'.*meshgrid((1-cos(RF*THETA))/RF,1:size(A,1));
    end
end

end