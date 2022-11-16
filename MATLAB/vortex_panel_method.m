
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
xNode=[xn2 xn1];   %node coordinates
yNode=[yn2 yn1];
scatter(xNode,yNode)
hold on
xlim([0 1])
ylim([-.5,.5])
xControl=(xNode(1:(nodes-1))+xNode(2:end))/2;    % Control points
yControl=(yNode(1:(nodes-1))+yNode(2:end))/2 ;    %
panelLength= sqrt((xNode(2:end)-xNode(1:(nodes-1))).^2+(xNode(2:end)-xNode(1:(nodes-1))).^2);

j=1;
k=2;
while j<nodes
  epsnu(:,j) =  (1./panelLength(j)).*[(xNode(k)-xNode(j)) (yNode(k)-yNode(j)); -(yNode(k)-yNode(j))... 
      (xNode(k)-xNode(j))]*[(xControl(j)-xNode(j)); (yControl(j)-yNode(j))] ;
j=j+1;
k = k+1;

end
    eps=epsnu(1,:);
    nu=epsnu(2,:);
phi=atan2(nu.*panelLength,(eps.^2+nu.^2-eps.*panelLength));
psi=.5*log((eps.^2+nu.^2)./((eps-panelLength).^2+nu.^2)); 

j=1; 
k=2;
AirfoilCoefficientMatrix=zeros(nodes);
AirfoilCoefficientMatrix(nodes,1)=1;
AirfoilCoefficientMatrix(nodes,nodes)=1;
while k<nodes
P=[(xNode(k)-xNode(j)) -(yNode(k)-yNode(j)); (yNode(k)-yNode(j))...
    (xNode(k)-xNode(j))]*[(panelLength(j)-eps(j)).*phi(j)+nu(j).*psi(j) ...
    eps(j).*phi(j)-nu(j).*psi(j); nu(j).*phi(j)-(panelLength(j)-eps(j)).*psi(j)-panelLength(j) ...
    (-nu(j).*phi(j)-eps(j).*psi(j)+panelLength(j))]
AirfoilCoefficientMatrix(j,j)=AirfoilCoefficientMatrix(j,j)+(xControl(k)-xControl(j))/panelLength(j)*P(2,1)...
    -(yControl(k)-yControl(j))/panelLength(j)*P(1,1);
AirfoilCoefficientMatrix(j,k)=AirfoilCoefficientMatrix(j,k)+(xControl(k)-xControl(j))/panelLength(j)*P(2,2)...
    -(yControl(k)-yControl(j))/panelLength(j)*P(1,2);



j=j+1;
k = k+1;

end
AirfoilCoefficientMatrix





%  end