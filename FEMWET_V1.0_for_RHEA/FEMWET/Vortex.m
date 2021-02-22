function [u,v,w] = Vortex(x,y,z,x1,y1,z1,x2,y2,z2,Gamma)
    % define vortex function
    % used for calculating the induced velocity due to a vortex element
    % as taken from Katz and Plotkin page: 584
    %
    % Induced velocities are calculated in point [x,y,z]
    % due to a vortex element with coordinates [x1,y1,z1]-[x2,y2,z2] and vortex strength Gam.
    % function [u,v,w] = vortex(x,y,z,x1,y1,z1,x2,y2,z2,Gam);    
    

    % Simulate infinity
    t11 = find(x1<0);
    t12 = find(x2<0);
    
    x1(isinf(x1)) = 100000;
    x2(isinf(x2)) = 100000;
    
    t21 = find(x1 == 100000);
    t22 = find(x2 == 100000);
    
    x1(t11-t21 == 0) = -100000;
    x2(t12-t22 == 0) = -100000;
    
    %calculation of R1 x R2
    
    R1R2X = (y-y1).*(z-z2)-(z-z1).*(y-y2);
    R1R2Y = -((x-x1).*(z-z2)-(z-z1).*(x-x2));
    R1R2Z = (x-x1).*(y-y2)-(y-y1).*(x-x2);
    
    square = R1R2X.^2 + R1R2Y.^2 + R1R2Z.^2;
    
    % Calculation of R0(R1/R(R1)-R2/R(R2))
    R1 = sqrt((x-x1).^2+(y-y1).^2+(z-z1).^2);
    R2 = sqrt((x-x2).^2+(y-y2).^2+(z-z2).^2);

    R0R1 = (x2-x1).*(x-x1)+(y2-y1).*(y-y1)+(z2-z1).*(z-z1);
    R0R2 = (x2-x1).*(x-x2)+(y2-y1).*(y-y2)+(z2-z1).*(z-z2);

    COEF = Gamma ./(4*pi*square).*(R0R1./R1-R0R2./R2);

    u = R1R2X.*COEF;
    v = R1R2Y.*COEF;
    w = R1R2Z.*COEF;

    % remove NaN element due to collocation point on vortex
    u(isnan(u)) = 0;
    v(isnan(v)) = 0;
    w(isnan(w)) = 0;
return