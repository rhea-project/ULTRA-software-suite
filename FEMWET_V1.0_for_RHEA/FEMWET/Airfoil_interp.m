function [yu yl x] = Airfoil_interp(Coord,Position,Y)

n = length(Position);
nx = length(Coord(1,:,1));

Xi = Coord(1,:,1);
Yi = Position;
Zui = zeros(nx,n);
Zli = zeros(nx,n);

for i = 1:n
    for j = 1:nx
        Zui(j,i) = Coord(i,j,2);
        Zli(j,i) = Coord(i,j,3);
    end
end

 x = Coord(1,:,1);
 yu = interp2(Yi,Xi,Zui,Y,x);
 yl = interp2(Yi,Xi,Zli,Y,x);