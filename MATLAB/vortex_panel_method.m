
alpha = 0;
chord_length=2;
cruise_speed=20;
nodes = 100;
clear eps nu epsnu

% function [Cl,Cmle] = vortex_panel_method(xu,xl,yu,yl,alpha,chord_length,cruise_speed,nodes)
dtheta = 2*pi/(nodes-1);
i = linspace(1,nodes/2,nodes/2);
xoc=ceil(2000*.5*(1-cos((i-.5)*dtheta))); %cosine clustering 
xn1=xu(xoc);
xn2=xl(xoc);

diff(xn1);
diff(xn2);    % checks if any 2 nodes are the same, should print error
% diff(yn1);
% diff(yn2);

yn1=yu(xoc);
yn2=yl(xoc);
xN=[xn2 xn1];   %node coordinates
yN=[yn2 yn1];
scatter(xN,yN)
hold on
xlim([0 1])
ylim([-.5,.5])
xC=(xN(1:(nodes-1))+xN(2:end))/2;    % Control points
yC=(yN(1:(nodes-1))+yN(2:end))/2 ;    %
pl= sqrt((xN(2:end)-xN(1:(nodes-1))).^2+(xN(2:end)-xN(1:(nodes-1))).^2);
x = sym('x');
y = sym('y');
j=1;
k=2;
while j<nodes
  epsnu(:,j) =  (1./pl(j)).*[(xN(k)-xN(j)) (yN(k)-yN(j)); -(yN(k)-yN(j))... 
      (xN(k)-xN(j))]*[(x-xN(j)); (y-yN(j))] ;
j=j+1;
k = k+1;

end
    eps=epsnu(1,:);
    nu=epsnu(2,:);
phi=atan2(nu.*pl,(eps.^2+nu.^2-eps.*pl));
psi=.5*log((eps.^2+nu.^2)./((eps-pl).^2+nu.^2));













%  end