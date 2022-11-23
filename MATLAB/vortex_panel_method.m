
alpha = 0;
chord_length=.23;
cruise_speed=20;
nodes = 10;
clear eps eta epseta AirfoilCoefficientMatrix 

% function [Cl,Cmle] = vortex_panel_method(xu,xl,yu,yl,alpha,chord_length,cruise_speed,nodes)

%% Cosine Clustering
dtheta = 2*pi/(nodes-1);
i = linspace(1,nodes/2,nodes/2);
xoc=ceil(2000*.5*(1-cos((i-.5)*dtheta))); %cosine clustering 
xn1=xu(xoc);
xn2=xl(xoc);

diff(xn1);
diff(xn2);    % checks if any 2 nodes are the same, should print error
% diff(yn1);
% diff(yn2);

%% Node Coordinates
yn1=yu(xoc);
yn2=yl(xoc);
xNode=[xn2 xn1];   %node coordinates
yNode=[yn2 yn1];
scatter(xNode,yNode)
hold on
xlim([0 1])
ylim([-.5,.5])
%% Control Point Coordinates
xControl=(xNode(1:(nodes-1))+xNode(2:end))/2;    % Control points
yControl=(yNode(1:(nodes-1))+yNode(2:end))/2 ;    %


%% Panel Length
panelLength= sqrt((xNode(2:end)-xNode(1:(nodes-1))).^2+(xNode(2:end)-xNode(1:(nodes-1))).^2);



%% Zeta and Eta matrices
RearNodeCounter=1;
ForwardNodeCounter=2;
ControlPointCounter=1;
zeta=zeros(nodes-1);
eta=zeros(nodes-1);
while ForwardNodeCounter<(nodes+1)
while ControlPointCounter<nodes
  zetaeta=(1./panelLength(RearNodeCounter)).*[(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)) (yNode(ForwardNodeCounter)-yNode(RearNodeCounter)); -(yNode(ForwardNodeCounter)-yNode(RearNodeCounter))... 
      (xNode(ForwardNodeCounter)-xNode(RearNodeCounter))]*[(xControl(ControlPointCounter)-xNode(RearNodeCounter)); (yControl(ControlPointCounter)-yNode(RearNodeCounter))]; 
  zeta(RearNodeCounter,ControlPointCounter)=zeta(RearNodeCounter,ControlPointCounter)+zetaeta(1,1);
 
  eta(RearNodeCounter,ControlPointCounter)=eta(RearNodeCounter,ControlPointCounter)+zetaeta(2,1);
ControlPointCounter=ControlPointCounter+1;


end
RearNodeCounter=RearNodeCounter+1;
ForwardNodeCounter=ForwardNodeCounter+1;
ControlPointCounter=1;
end

eta=eta(1:(nodes-1),1:(nodes-1));



%% Phi and Psi Matrics
    RearNodeCounter=1;
    phi=zeros(nodes-1);
    psi=zeros(nodes-1);
    
while RearNodeCounter<nodes
phi(RearNodeCounter,:)=phi(RearNodeCounter,:)+ atan2(eta(RearNodeCounter,:).*panelLength(RearNodeCounter),(zeta(RearNodeCounter,:).^2+eta(RearNodeCounter,:).^2-zeta(RearNodeCounter,:).*panelLength(RearNodeCounter)))
psi(RearNodeCounter,:)=psi(RearNodeCounter,:)+.5*log((zeta(RearNodeCounter,:).^2+eta(RearNodeCounter,:).^2)./((zeta(RearNodeCounter)-panelLength(RearNodeCounter)).^2+eta(RearNodeCounter,:).^2))
RearNodeCounter=RearNodeCounter+1;
end



%% Airfoil Coefficient Matrix
RearNodeCounter=1; 
ForwardNodeCounter=2;
AirfoilCoefficientMatrix=zeros(nodes);
AirfoilCoefficientMatrix(nodes,1)=1;
AirfoilCoefficientMatrix(nodes,nodes)=1;
ControlPointCounter=1;
ForwardControlPointCounter=2;
while ForwardNodeCounter<nodes
    while ForwardControlPointCounter<nodes
P=(1/2/pi/panelLength(RearNodeCounter)^2)*[(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)) -(yNode(ForwardNodeCounter)-yNode(RearNodeCounter)); (yNode(ForwardNodeCounter)-yNode(RearNodeCounter)) (xNode(ForwardNodeCounter)-xNode(RearNodeCounter))]*[(panelLength(RearNodeCounter)-zeta(RearNodeCounter,ControlPointCounter)).*phi(RearNodeCounter,ControlPointCounter)+eta(RearNodeCounter,ControlPointCounter).*psi(RearNodeCounter,ControlPointCounter) ...
    zeta(RearNodeCounter,ControlPointCounter).*phi(RearNodeCounter,ControlPointCounter)-eta(RearNodeCounter,ControlPointCounter).*psi(RearNodeCounter,ControlPointCounter); eta(RearNodeCounter,ControlPointCounter).*phi(RearNodeCounter,ControlPointCounter)-(panelLength(RearNodeCounter)-zeta(RearNodeCounter,ControlPointCounter)).*psi(RearNodeCounter,ControlPointCounter)-panelLength(RearNodeCounter) ...
    (-eta(RearNodeCounter,ControlPointCounter).*phi(RearNodeCounter,ControlPointCounter)-zeta(RearNodeCounter,ControlPointCounter).*psi(RearNodeCounter,ControlPointCounter)+panelLength(RearNodeCounter))];
AirfoilCoefficientMatrix(RearNodeCounter,ControlPointCounter)=AirfoilCoefficientMatrix(RearNodeCounter,RearNodeCounter)+(xControl(ForwardNodeCounter)-xControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(2,1)...
    -(yControl(ForwardNodeCounter)-yControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(1,1);
AirfoilCoefficientMatrix(ForwardNodeCounter,ControlPointCounter)=AirfoilCoefficientMatrix(ForwardNodeCounter,ControlPointCounter)+(xControl(ControlPointCounter)-xControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(2,2)...
    -(yControl(ForwardNodeCounter)-yControl(RearNodeCounter))/panelLength(RearNodeCounter)*P(1,2);
ControlPointCounter=ControlPointCounter+1;

    end
RearNodeCounter=RearNodeCounter+1;
ForwardNodeCounter = ForwardNodeCounter+1;
ControlPointCounter=1;
end
AirfoilCoefficientMatrix

PanelSolution=zeros(1,nodes);
RearNodeCounter=1;
ForwardNodeCounter=2;
while RearNodeCounter<nodes
    PanelSolution(RearNodeCounter)=PanelSolution(RearNodeCounter)+cruise_speed*((yNode(ForwardNodeCounter)-yNode(RearNodeCounter))*cos(alpha)...
        -(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)*sin(alpha)))/panelLength(RearNodeCounter);
RearNodeCounter=RearNodeCounter+1;
ForwardNodeCounter = ForwardNodeCounter+1;

end
PanelSolution=PanelSolution';
VortexStrengths=PanelSolution\AirfoilCoefficientMatrix;


VortexStrengthCounter=1;
ForwardVortexStrengthCounter=1;
InducedLifts=zeros(1,(nodes-1));
while VortexStrengthCounter<nodes
    InducedLifts(VortexStrengthCounter)=InducedLifts(VortexStrengthCounter)...
        +panelLength(VortexStrengthCounter)/chord_length*(VortexStrengths(VortexStrengthCounter)...
        +VortexStrengths(ForwardVortexStrengthCounter))/cruise_speed;
    

    
VortexStrengthCounter=VortexStrengthCounter+1;
ForwardVortexStrengthCounter=ForwardVortexStrengthCounter+1;
    
    


end
LiftCoeff=sum(InducedLifts)

%% Plot for testing alpha vs Cl

alphacounter=-.4
degree2Rad=180/pi;
LiftCoefftest=zeros(1,20)
testcounter=1
while alphacounter<.57
    clear PanelSolution InducedLifts
PanelSolution=zeros(1,nodes);
RearNodeCounter=1;
ForwardNodeCounter=2;
alpha=alphacounter;
while RearNodeCounter<nodes
    PanelSolution(RearNodeCounter)=PanelSolution(RearNodeCounter)+cruise_speed*((yNode(ForwardNodeCounter)-yNode(RearNodeCounter))*cos(alpha)...
        -(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)*sin(alpha)))/panelLength(RearNodeCounter);
RearNodeCounter=RearNodeCounter+1;
ForwardNodeCounter = ForwardNodeCounter+1;

end
PanelSolution=PanelSolution';
VortexStrengths=PanelSolution\AirfoilCoefficientMatrix;


VortexStrengthCounter=1;
ForwardVortexStrengthCounter=1;
InducedLifts=zeros(1,(nodes-1));
while VortexStrengthCounter<nodes
    InducedLifts(VortexStrengthCounter)=InducedLifts(VortexStrengthCounter)...
        +panelLength(VortexStrengthCounter)/chord_length*(VortexStrengths(VortexStrengthCounter)...
        +VortexStrengths(ForwardVortexStrengthCounter))/cruise_speed;
    

    
VortexStrengthCounter=VortexStrengthCounter+1;
ForwardVortexStrengthCounter=ForwardVortexStrengthCounter+1;
    
    


end
LiftCoefftest(testcounter)=sum(InducedLifts);
alphacounter=alphacounter+.05;
testcounter=testcounter+1;

end
alphavec=linspace(.24,32,20);
figure
plot(alphavec,LiftCoefftest)
%  end