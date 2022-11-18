
alpha = 0;
chord_length=2;
cruise_speed=20;
nodes = 10;
clear eps eta epseta

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

RearNodeCounter=1;
ForwardNodeCounter=2;
ControlPointCounter=1;
epsilon=zeros(nodes);
eta=zeros(nodes);
while RearNodeCounter<nodes
while ControlPointCounter<nodes
  epseta=(1./panelLength(RearNodeCounter)).*[(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)) (yNode(ForwardNodeCounter)-yNode(RearNodeCounter)); -(yNode(ForwardNodeCounter)-yNode(RearNodeCounter))... 
      (xNode(ForwardNodeCounter)-xNode(RearNodeCounter))]*[(xControl(ControlPointCounter)-xNode(RearNodeCounter)); (yControl(ControlPointCounter)-yNode(RearNodeCounter))]; 
  epsilon(RearNodeCounter,ControlPointCounter)=epsilon(RearNodeCounter,ControlPointCounter)+epseta(1,1);
 
  eta(RearNodeCounter,ControlPointCounter)=eta(RearNodeCounter,ControlPointCounter)+epseta(2,1);
ControlPointCounter=ControlPointCounter+1;


end
RearNodeCounter=RearNodeCounter+1;
ForwardNodeCounter=ForwardNodeCounter+1;
ControlPointCounter=1;
end
epsilon=epsilon(1:(nodes-1),1:(nodes-1));
eta=eta(1:(nodes-1),1:(nodes-1));

    RearNodeCounter=1;
while RearNodeCounter<nodes
phi=atan2(eta.*panelLength,(epsilon(RearNodeCounter,:).^2+eta.^2-epsilon(RearNodeCounter,:).*panelLength(RearNodeCounter)));
psi=.5*log((epsilon(RearNodeCounter,:).^2+eta.^2)./((epsilon(RearNodeCounter)-panelLength(RearNodeCounter)).^2+eta(RearNodeCounter).^2)); 
RearNodeCounter=RearNodeCounter+1;
end
RearNodeCounter=1; 
ForwardNodeCounter=2;
AirfoilCoefficientMatrix=zeros(nodes);
AirfoilCoefficientMatrix(nodes,1)=1;
AirfoilCoefficientMatrix(nodes,nodes)=1;
ControlPointCounter=1;
while ForwardNodeCounter<nodes
    while ControlPointCounter<nodes
P=(1/2/pi/panelLength(RearNodeCounter)^2)*[(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)) -(yNode(ForwardNodeCounter)-yNode(RearNodeCounter)); (yNode(ForwardNodeCounter)-yNode(RearNodeCounter))...
    (xNode(ForwardNodeCounter)-xNode(RearNodeCounter))]*[(panelLength(RearNodeCounter)-epsilon(RearNodeCounter,ControlPointCounter)).*phi(RearNodeCounter,ControlPointCounter)+eta(RearNodeCounter,ControlPointCounter).*psi(RearNodeCounter,ControlPointCounter) ...
    epsilon(RearNodeCounter,ControlPointCounter).*phi(RearNodeCounter,ControlPointCounter)-eta(RearNodeCounter,ControlPointCounter).*psi(RearNodeCounter,ControlPointCounter); eta(RearNodeCounter,ControlPointCounter).*phi(RearNodeCounter,ControlPointCounter)-(panelLength(RearNodeCounter)-epsilon(RearNodeCounter,ControlPointCounter)).*psi(RearNodeCounter,ControlPointCounter)-panelLength(RearNodeCounter) ...
    (-eta(RearNodeCounter,ControlPointCounter).*phi(RearNodeCounter,ControlPointCounter)-epsilon(RearNodeCounter,ControlPointCounter).*psi(RearNodeCounter,ControlPointCounter)+panelLength(RearNodeCounter))]
AirfoilCoefficientMatrix(RearNodeCounter,ControlPointCounter)=AirfoilCoefficientMatrix(RearNodeCounter,RearNodeCounter)+(xControl(ForwardNodeCounter)-xControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(2,1)...
    -(yControl(ForwardNodeCounter)-yControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(1,1);
AirfoilCoefficientMatrix(RearNodeCounter,ForwardNodeCounter)=AirfoilCoefficientMatrix(RearNodeCounter,ForwardNodeCounter)+(xControl(ForwardNodeCounter)-xControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(2,2)...
    -(yControl(ForwardNodeCounter)-yControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(1,2);
ControlPointCounter=ControlPointCounter+1;

    end
RearNodeCounter=RearNodeCounter+1;
ForwardNodeCounter = ForwardNodeCounter+1;

end
AirfoilCoefficientMatrix





%  end