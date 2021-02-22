function Wing_plot(AC,nNode,cNode,u,v,w,px,py,pz,ModePlot,ModeNum,ModeFr)

cNode = [cNode(1:nNode) cNode(nNode+1:2*nNode) cNode(2*nNode+1:3*nNode)];

%% Intial geometry

% airfoils

% Xa  = [Geo.Wing.Airfoil_coord(:,:,1)'; flipud(Geo.Wing.Airfoil_coord(:,:,1)')].*(ones(2*size(Geo.Wing.Airfoil_coord(:,:,1),2),1)*Geo.Wing.Chord);
% Xa = Xa + (ones(2*size(Geo.Wing.Airfoil_coord(:,:,1),2),1)*(Geo.Wing.Xle-Geo.Wing.Xle(1)));
% Za  = [Geo.Wing.Airfoil_coord(:,:,2)'; flipud(Geo.Wing.Airfoil_coord(:,:,3)')].*(ones(2*size(Geo.Wing.Airfoil_coord(:,:,1),2),1)*Geo.Wing.Chord);
% Za = Za + (ones(2*size(Geo.Wing.Airfoil_coord(:,:,1),2),1)*(Geo.Wing.Zle-Geo.Wing.Zle(1)));
% Ya  = (ones(2*size(Geo.Wing.Airfoil_coord(:,:,1),2),1)*(Geo.Wing.Yle-Geo.Wing.Yle(1)));

for i=1:size(AC.Wing.Geom,1)
    [yu yl x] = Airfoil_interp(AC.Wing.Airfoil_coord,AC.Wing.eta,AC.Wing.Geom(i,2)/AC.Wing.Geom(end,2));
    Airfoils(i,:,1) = x;
    Airfoils(i,:,2) = yu;
    Airfoils(i,:,3) = yl;
end

Xa  = [ Airfoils(:,:,1)'; flipud( Airfoils(:,:,1)')].*(ones(2*size( Airfoils(:,:,1),2),1)*AC.Wing.Geom(:,4)');
Xa = Xa + (ones(2*size( Airfoils(:,:,1),2),1)*(AC.Wing.Geom(:,1)-AC.Wing.Geom(1,1))');
Za  = [ Airfoils(:,:,2)'; flipud( Airfoils(:,:,3)')].*(ones(2*size( Airfoils(:,:,1),2),1)*AC.Wing.Geom(:,4)');
Za = Za + (ones(2*size( Airfoils(:,:,1),2),1)*(AC.Wing.Geom(:,3)-AC.Wing.Geom(1,3))');
Ya  = (ones(2*size( Airfoils(:,:,1),2),1)*(AC.Wing.Geom(:,2)-AC.Wing.Geom(1,2))');

Y = cNode(:,2)';
X = interp1(Ya(1,:)',Xa',Y')';
Z = interp1(Ya(1,:)',Za',Y')';


% twist
if sum(AC.Wing.Geom(:,5)) ~=0
    T = -interp1(AC.Wing.Geom(:,2),AC.Wing.Geom(:,5),Y);
else 
    T = zeros(length(Y),1)';
end

Y = ones(2*size( Airfoils(:,:,1),2),1)*Y;

for i=1:size(X,2);
    xt = min(X(:,i))+0.25*(max(X(:,i))-min(X(:,i)));
    yt = Z(1,i);
    [X(:,i), Z(:,i)] = twist(X(:,i),Z(:,i),T(i),xt,yt);
end


% planform

[Xlp, nl] = min(X);
Ylp = diag(Y(nl,:));
Zlp = diag(Z(nl,:));

[Xtp, nt] = max(X);
Ytp = diag(Y(nt,:));
Ztp = diag(Z(nt,:));


%% Deformed wing


Xn = X + ones(size(X,1),1)*u';
Yn = Y + ones(size(X,1),1)*v';
Zn = Z + ones(size(X,1),1)*w';

cNode_new = cNode + [u v w];

for i=1:size(Xn,2);
    xt = cNode_new(i,2);
    yt = cNode_new(i,3);
    [Yn(:,i), Zn(:,i)] = twist(Yn(:,i),Zn(:,i),px(i),xt,yt);
end

for i=1:size(Xn,2);
    xt = cNode_new(i,1);
    yt = cNode_new(i,3);
    [Xn(:,i), Zn(:,i)] = twist(Xn(:,i),Zn(:,i),-py(i),xt,yt);
end


[Xlpn, nln] = min(Xn);
Ylpn = diag(Yn(nln,:));
Zlpn = diag(Zn(nln,:));

[Xtpn, ntn] = max(Xn);
Ytpn = diag(Yn(ntn,:));
Ztpn = diag(Zn(ntn,:));

%%

if ModePlot ==0
    figure
    hold on
    plot3(Xlp,Ylp,Zlp,'-b');
    plot3(Xtp,Ytp,Ztp,'-b');
    plot3(Xlpn,Ylpn,Zlpn,'-r');
    plot3(Xtpn,Ytpn,Ztpn,'-r');

    for i=1: size(X,2)
        plot3(X(:,i),Y(:,i),Z(:,i),'-b');
        plot3(Xn(:,i),Yn(:,i),Zn(:,i),'-r');
    end
    % plot3(cNode(:,1),cNode(:,2),cNode(:,3),'o-k');
    axis equal
    view(3)
    title('Wing deformation under 2.5g static load')
    hold off
elseif ModePlot ==1
    figure
    hold on
    plot3(Xlpn,Ylpn,Zlpn,'-r');
    plot3(Xtpn,Ytpn,Ztpn,'-r');
    for i=1: size(X,2)
        plot3(Xn(:,i),Yn(:,i),Zn(:,i),'-r');
    end
    axis equal
    title(['Mode number ' num2str(ModeNum) ', Frequency ' num2str(ModeFr) ' Hz'] )
    view(3)
    hold off
%     print('-r300','-djpeg',['Mode' num2str(ModeNum) '3D'])
    
    figure
    hold on
    plot3(cNode(:,1),cNode(:,2),cNode(:,3),'-k','linewidth',1.5);
    plot3(cNode_new(:,1),cNode_new(:,2),cNode_new(:,3),'-r','linewidth',1.5);
    title(['Mode number ' num2str(ModeNum) ', Frequency ' num2str(ModeFr) ' Hz'] )
    view([-90 0])
    axis equal
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off
%     print('-r300','-djpeg',['Mode' num2str(ModeNum) '2D'])
end



end

function [X, Y] = twist(x,y,twist,x_twist,y_twist)

xt = x - x_twist;
yt = y - y_twist;

X = xt*cosd(twist) - yt*sind(twist) + x_twist;
Y = xt*sind(twist) + yt*cosd(twist) + y_twist;
end


function [yu, yl, x] = Airfoil_interp(Coord,Position,Y)

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
 
end

