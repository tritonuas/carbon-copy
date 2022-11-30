
alpha = 0;
chord_length=.23;
cruise_speed=20;
nodes = 50;
alpha=0;
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
xn2=flip(xn2);
yn2=flip(yn2);
xNode=[xn2 xn1];   %node coordinates
yNode=[yn2 yn1];
figure
scatter(xNode,yNode)
xlim([0 1])
ylim([-.5,.5])
% hold on
% for i=1:nodes
% scatter(xNode(i),yNode(i),'r')
% % hold on
% end
% xlim([0 1])
% ylim([-.5,.5])
%% Control Point Coordinates
% xNode=zeros(1,nodes);
% yControl=zeros(1,nodes);
for NodeCounter=1:(nodes-1)
xControl(NodeCounter)=(xNode(NodeCounter+1)+xNode(NodeCounter))/2;    % Control points
yControl(NodeCounter)=(yNode(NodeCounter+1)+yNode(NodeCounter))/2;   %

end
xControl=[xControl 0];
yControl=[yControl 0];
%% Panel Length
panelLength= sqrt((xNode(2:end)-xNode(1:(nodes-1))).^2+(yNode(2:end)-yNode(1:(nodes-1))).^2);



%% Zeta and Eta matrices
RearNodeCounter=1;
% ForwardNodeCounter=2;
% ControlPointCounter=1;
zeta=zeros(nodes-1);
eta=zeros(nodes-1);
for ForwardNodeCounter=2:nodes
    ForwardNodeCounter
    for ControlPointCounter=1:(nodes-1)
        ControlPointCounter
  zetaeta=(1./panelLength(RearNodeCounter)).*[(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)) (yNode(ForwardNodeCounter)-yNode(RearNodeCounter)); -(yNode(ForwardNodeCounter)-yNode(RearNodeCounter))... 
      (xNode(ForwardNodeCounter)-xNode(RearNodeCounter))]*[(xControl(ControlPointCounter)-xNode(RearNodeCounter)); (yControl(ControlPointCounter)-yNode(RearNodeCounter))];
  zeta(RearNodeCounter,ControlPointCounter)=zeta(RearNodeCounter,ControlPointCounter)+zetaeta(1,1);
 
  eta(RearNodeCounter,ControlPointCounter)=eta(RearNodeCounter,ControlPointCounter)+zetaeta(2,1);
% ControlPointCounter=ControlPointCounter+1;


    end
RearNodeCounter=RearNodeCounter+1
% ForwardNodeCounter=ForwardNodeCounter+1;
% ControlPointCounter=1;
end

eta=eta(1:(nodes-1),1:(nodes-1));



%% Phi and Psi Matrics
%     RearNodeCounter=1;
    phi=zeros(nodes-1);
    psi=zeros(nodes-1);
    
for PanelCounter=1:(nodes-1)
    PanelCounter
phi(PanelCounter,:)=phi(PanelCounter,:)+ atan2(eta(PanelCounter,:).*panelLength(PanelCounter),((zeta(PanelCounter,:).^2+eta(PanelCounter,:).^2-zeta(PanelCounter,:).*panelLength(PanelCounter))));
psi(PanelCounter,:)=psi(PanelCounter,:)+.5*log((zeta(PanelCounter,:).^2+eta(PanelCounter,:).^2)./((zeta(PanelCounter)-panelLength(PanelCounter)).^2+eta(PanelCounter,:).^2));
% RearNodeCounter=RearNodeCounter+1;
end
% zeta=zeta';
% eta=eta';
% phi=phi';
% psi=psi';


%% Airfoil Coefficient Matrix

% ForwardNodeCounter=2;
AirfoilCoefficientMatrix=zeros(nodes);
ControlPointCounter=1;
ForwardControlPointCounter=2;
for ForwardNodeCounter= 2:nodes
    ForwardNodeCounter
    for ControlPointCounter = 1:(nodes-1)
        ControlPointCounter
P=((1/(2*pi*((panelLength(ForwardNodeCounter-1))^2))))*[(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1)) -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1)); (yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1)) (xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))]*[(panelLength(ForwardNodeCounter-1)-zeta(ForwardNodeCounter-1,ControlPointCounter))*phi(ForwardNodeCounter-1,ControlPointCounter)+eta(ForwardNodeCounter-1,ControlPointCounter)*psi(ForwardNodeCounter-1,ControlPointCounter) ...
    zeta(ForwardNodeCounter-1,ControlPointCounter)*phi(ForwardNodeCounter-1,ControlPointCounter)-eta(ForwardNodeCounter-1,ControlPointCounter)*psi(ForwardNodeCounter-1,ControlPointCounter); eta(ForwardNodeCounter-1,ControlPointCounter).*phi(ForwardNodeCounter-1,ControlPointCounter)-(panelLength(ForwardNodeCounter-1)-zeta(ForwardNodeCounter-1,ControlPointCounter))*psi(ForwardNodeCounter-1,ControlPointCounter)-panelLength(ForwardNodeCounter-1) ...
    (-eta(ForwardNodeCounter-1,ControlPointCounter).*phi(ForwardNodeCounter-1,ControlPointCounter)-zeta(ForwardNodeCounter-1,ControlPointCounter).*psi(ForwardNodeCounter-1,ControlPointCounter)+panelLength(ForwardNodeCounter-1))];
AirfoilCoefficientMatrix(ControlPointCounter,ForwardNodeCounter-1)=AirfoilCoefficientMatrix(ForwardNodeCounter-1,ForwardNodeCounter-1)+(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))/panelLength(ControlPointCounter)*P(2,1)...
    -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1))/panelLength(ControlPointCounter)*P(1,1);
AirfoilCoefficientMatrix(ControlPointCounter,ForwardNodeCounter)=AirfoilCoefficientMatrix(ForwardNodeCounter,ControlPointCounter)+(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))/panelLength(ControlPointCounter)*P(2,2)...
    -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1))/panelLength(ControlPointCounter)*P(1,2);
% ControlPointCounter=ControlPointCounter+1;

    end

% ForwardNodeCounter = ForwardNodeCounter+1;
% ControlPointCounter=1;
end
% AirfoilCoefficientMatrix=AirfoilCoefficientMatrix'
AirfoilCoefficientMatrix(nodes,nodes)=1
AirfoilCoefficientMatrix(nodes,1)=1



PanelSolution=zeros(1,nodes);
ForwardNodeCounter=2;
for RearNodeCounter = 1:(nodes-1)
    RearNodeCounter
    PanelSolution(RearNodeCounter)=PanelSolution(RearNodeCounter)+cruise_speed*((yNode(ForwardNodeCounter)-yNode(RearNodeCounter))*cos(alpha)...
        -((xNode(ForwardNodeCounter)-xNode(RearNodeCounter))*sin(alpha)))/panelLength(RearNodeCounter);

ForwardNodeCounter = ForwardNodeCounter+1;

end
PanelSolution=PanelSolution';
VortexStrengths=PanelSolution\AirfoilCoefficientMatrix;

% 
% VortexStrengthCounter=1;
ForwardVortexStrengthCounter=2;
InducedLifts=zeros(1,(nodes-1));
for VortexStrengthCounter=1:(nodes-1)
    InducedLifts(VortexStrengthCounter)=InducedLifts(VortexStrengthCounter)...
        +(panelLength(VortexStrengthCounter)/chord_length)*(VortexStrengths(VortexStrengthCounter)...
        +VortexStrengths(ForwardVortexStrengthCounter))/cruise_speed;
    

    
ForwardVortexStrengthCounter=ForwardVortexStrengthCounter+1;
    
    


end
LiftCoeff=sum(InducedLifts)

%% Plot for testing alpha vs Cl

% alphacounter=-.4;
% degree2Rad=180/pi;
% LiftCoefftest=zeros(1,20);
% testcounter=1;
% while alphacounter<.57
%     clear PanelSolution InducedLifts
% PanelSolution=zeros(1,nodes);
% RearNodeCounter=1;
% ForwardNodeCounter=2;
% alpha=alphacounter;
% for RearNodeCounter=1:(nodes-1)
%     PanelSolution(RearNodeCounter)=PanelSolution(RearNodeCounter)+cruise_speed*((yNode(ForwardNodeCounter)-yNode(RearNodeCounter))*cos(alpha)...
%         -(xNode(ForwardNodeCounter)-xNode(RearNodeCounter)*sin(alpha)))/panelLength(RearNodeCounter);
% % RearNodeCounter=RearNodeCounter+1;
% ForwardNodeCounter = ForwardNodeCounter+1;
% 
% end
% PanelSolution=PanelSolution';
% VortexStrengths=PanelSolution\AirfoilCoefficientMatrix;
% 
% 
% VortexStrengthCounter=1;
% ForwardVortexStrengthCounter=1;
% InducedLifts=zeros(1,(nodes-1));
% for VortexStrengthCounter=1:(nodes-1)
%     InducedLifts(VortexStrengthCounter)=InducedLifts(VortexStrengthCounter)...
%         +panelLength(VortexStrengthCounter)/chord_length*(VortexStrengths(VortexStrengthCounter)...
%         +VortexStrengths(ForwardVortexStrengthCounter))/cruise_speed;
%     
% 
%     
% % VortexStrengthCounter=VortexStrengthCounter+1;
% ForwardVortexStrengthCounter=ForwardVortexStrengthCounter+1;
%     
%     
% 
% 
% end
% LiftCoefftest(testcounter)=sum(InducedLifts);
% alphacounter=alphacounter+.05;
% testcounter=testcounter+1;
% 
% end
% alphavec=linspace(24,32,20);
% figure
% plot(alphavec,LiftCoefftest)
% %  end