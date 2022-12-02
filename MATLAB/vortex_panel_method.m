clear all
close all
clc
alpha = 0;
chord_length=.23;
cruise_speed=20;
nodes = 6;



% function [Cl,Cmle] = vortex_panel_method(xu,xl,yu,yl,alpha,chord_length,cruise_speed,nodes)
% 
run get_foil_shape.m
%% Cosine Clustering
dtheta = 2*pi/(nodes-1);
i = linspace(1,nodes/2,nodes/2);
spacing_coeff=length(xu);
xoc=ceil(spacing_coeff*.5*(1-cos((i-.5)*dtheta))); %cosine clustering 
xn1=xu(xoc);
xn2=xl(xoc);

diff(xn1);
diff(xn2);    % checks if any 2 nodes are the same, creating panel length 
              % of 0 should print error
% diff(yn1);
% diff(yn2);

%% Node Coordinates
yn1=yu(xoc);
yn2=yl(xoc);
xn2=flip(xn2);
yn2=flip(yn2);
xNode=[xn2 xn1];   %node coordinates beginning from the trailing edge,
yNode=[yn2 yn1];   %clockwise about the airfoil
figure
scatter(xNode,yNode)
xlim([0 1])
ylim([-.5,.5])

% Step through here to see node coordinate order around airfoil

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

% Control points are established at the centers of each panel
% Order is the same convention as control points

for NodeCounter=1:(nodes-1)
xControl(NodeCounter)=(xNode(NodeCounter+1)+xNode(NodeCounter))/2;    % Control points
yControl(NodeCounter)=(yNode(NodeCounter+1)+yNode(NodeCounter))/2;   %

end

%Possibly necessary for airfoil matrix using "n" control points ? See line
% xControl=[xControl 0];
% yControl=[yControl 0];


% Step through here to see Control point order around airfoil

% figure
% hold on
% xlim([0 1])
% ylim([-.5,.5])
% for i=1:nodes
% scatter(xControl(i),yControl(i),'r')
% % hold on
% end
% xlim([0 1])
% ylim([-.5,.5])
%% Panel Length
% panelLength= sqrt((xNode(2:end)-xNode(1:(nodes-1))).^2+(yNode(2:end)-yNode(1:(nodes-1))).^2);
% Panel length is the distance measured between two consecutive nodes
% Follows the same clockwise order as nodes and control points

for NodeCounter=2:nodes
panelLength(NodeCounter-1)= sqrt((xNode(NodeCounter)-xNode(NodeCounter-1))^2+(yNode(NodeCounter)-yNode(NodeCounter-1))^2);
    
    
end

%% Zeta and Eta matrices

% Zeta and Eta are the parallel and normal directions with respect to
% the panel line
% the for loop creates a matrix evaluating the zeta and eta components of
% induced velocity at each control point


% The rows represent the relevant panel number and the columns represent the
% control point number (panel#,control#)



zeta=zeros(nodes-1);
eta=zeros(nodes-1);
for ForwardNodeCounter=2:nodes
    ForwardNodeCounter
    for ControlPointCounter=1:(nodes-1)
        ControlPointCounter
  zetaeta=(1/panelLength(ForwardNodeCounter-1))*[(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1)) (yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1)); -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1))... 
      (xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))]*[(xControl(ControlPointCounter)-xNode(ForwardNodeCounter-1)); (yControl(ControlPointCounter)-yNode(ForwardNodeCounter-1))];
  zeta(ForwardNodeCounter-1,ControlPointCounter)=zetaeta(1,1);
 
  eta(ForwardNodeCounter-1,ControlPointCounter)=zetaeta(2,1);



    end


end





%% Phi and Psi Matrics
%     RearNodeCounter=1;
    phi=zeros(nodes-1);
    psi=zeros(nodes-1);
    
    % Phi represents the angle from the control point to each nodes in a
    % specific panel and should be an angle between -pi and pi. 
    % Should be equal to -+pi when the control point
    % iteration is the same as the panel iteration
    % Psi is the stream function where dpsi/dy=u and -dpsi/dx=v
    % 
    % Rows Should represent panel numbers and columns represent control points
    % Same as zeta and eta (panel#,control#)
 
% zeta=zeta'; Probably not necessary as long as I iterate properly but
            % would make them control point number x panel number
% eta=eta';

% panelLength=panelLength';
for PanelCounter=1:(nodes-1)
    PanelCounter
phi(PanelCounter,:) = atan2(eta(PanelCounter,:)*panelLength(PanelCounter),((zeta(PanelCounter,:).^2+eta(PanelCounter,:).^2-zeta(PanelCounter,:).*panelLength(PanelCounter))));
psi(PanelCounter,:) = .5*log((zeta(PanelCounter,:).^2+eta(PanelCounter,:).^2)./((zeta(PanelCounter,:)-panelLength(PanelCounter)).^2+eta(PanelCounter,:).^2));

end



%% Airfoil Coefficient Matrix

% ACM is populated by P matrix determined by velocity induced at a control
% point 1:n-1 by a panel 1:n-1
% 
% The ACM should have rows representing the node number and the columns 
% representing the control point number (node#,control#) with n rows and 
% the nth column being the kutta condition




AirfoilCoefficientMatrix=zeros(nodes);
for ForwardNodeCounter= 2:nodes
    ForwardNodeCounter
    for ControlPointCounter = 1:(nodes-1)
        ControlPointCounter
P=((1/(2*pi*((panelLength(ForwardNodeCounter-1))^2))))*[(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1)) -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1)); (yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1)) (xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))]*[(panelLength(ForwardNodeCounter-1)-zeta(ForwardNodeCounter-1,ControlPointCounter))*phi(ForwardNodeCounter-1,ControlPointCounter)+eta(ForwardNodeCounter-1,ControlPointCounter)*psi(ForwardNodeCounter-1,ControlPointCounter) ...
    zeta(ForwardNodeCounter-1,ControlPointCounter)*phi(ForwardNodeCounter-1,ControlPointCounter)-eta(ForwardNodeCounter-1,ControlPointCounter)*psi(ForwardNodeCounter-1,ControlPointCounter); eta(ForwardNodeCounter-1,ControlPointCounter)*phi(ForwardNodeCounter-1,ControlPointCounter)-(panelLength(ForwardNodeCounter-1)-zeta(ForwardNodeCounter-1,ControlPointCounter))*psi(ForwardNodeCounter-1,ControlPointCounter)-panelLength(ForwardNodeCounter-1) ...
    (-eta(ForwardNodeCounter-1,ControlPointCounter)*phi(ForwardNodeCounter-1,ControlPointCounter)-zeta(ForwardNodeCounter-1,ControlPointCounter)*psi(ForwardNodeCounter-1,ControlPointCounter)+panelLength(ForwardNodeCounter-1))];
AirfoilCoefficientMatrix(ForwardNodeCounter-1,ControlPointCounter)=AirfoilCoefficientMatrix(ForwardNodeCounter-1,ControlPointCounter)+(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))/panelLength(ForwardNodeCounter-1)*P(2,1)...
    -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1))/panelLength(ForwardNodeCounter-1)*P(1,1);
AirfoilCoefficientMatrix(ForwardNodeCounter-1,ControlPointCounter+1)=AirfoilCoefficientMatrix(ForwardNodeCounter-1,ControlPointCounter+1)+(xNode(ForwardNodeCounter)-xNode(ForwardNodeCounter-1))/panelLength(ForwardNodeCounter-1)*P(2,2)...
    -(yNode(ForwardNodeCounter)-yNode(ForwardNodeCounter-1))/panelLength(ForwardNodeCounter-1)*P(1,2);


    end


end

AirfoilCoefficientMatrix(nodes,nodes)=1;
AirfoilCoefficientMatrix(nodes,1)=1;
% AirfoilCoefficientMatrix=AirfoilCoefficientMatrix'


PanelSolution=zeros(1,nodes);

for RearNodeCounter = 1:(nodes-1)
    RearNodeCounter
    PanelSolution(RearNodeCounter)=((yNode(RearNodeCounter+1)-yNode(RearNodeCounter))*cos(alpha)...
        -((xNode(RearNodeCounter+1)-xNode(RearNodeCounter))*sin(alpha)))/panelLength(RearNodeCounter);
PanelSolution;


end
PanelSolution=cruise_speed.*PanelSolution';
VortexStrengths=AirfoilCoefficientMatrix\PanelSolution

% 
% VortexStrengthCounter=1;

InducedLifts=zeros(1,(nodes-1));
for VortexStrengthCounter=1:(nodes-1)
    InducedLifts(VortexStrengthCounter)=(panelLength(VortexStrengthCounter)/chord_length)*(VortexStrengths(VortexStrengthCounter)...
        +VortexStrengths(VortexStrengthCounter+1))/cruise_speed
    

    

    
    


end
LiftCoeff=sum(InducedLifts)

%% Plot for testing alpha vs Cl

alphacounter=-.4;
degree2Rad=180/pi;
LiftCoefftest=zeros(1,20);
testcounter=1;
while alphacounter<.57
    clear PanelSolution InducedLifts
PanelSolution=zeros(nodes,1);
RearNodeCounter=1;
ForwardNodeCounter=2;
alpha=alphacounter;
for RearNodeCounter=1:(nodes-1)
     RearNodeCounter
    PanelSolution(RearNodeCounter)=PanelSolution(RearNodeCounter)+cruise_speed*((yNode(RearNodeCounter+1)-yNode(RearNodeCounter))*cos(alpha)...
        -((xNode(RearNodeCounter+1)-xNode(RearNodeCounter))*sin(alpha)))/panelLength(RearNodeCounter);
PanelSolution;

end

VortexStrengths=AirfoilCoefficientMatrix\PanelSolution


VortexStrengthCounter=1;
ForwardVortexStrengthCounter=1;
InducedLifts=zeros(1,(nodes-1));
for VortexStrengthCounter=1:(nodes-1)
   InducedLifts(VortexStrengthCounter)=InducedLifts(VortexStrengthCounter)...
        +(panelLength(VortexStrengthCounter)/chord_length)*(VortexStrengths(VortexStrengthCounter)...
        +VortexStrengths(VortexStrengthCounter+1))/cruise_speed
    
    


end
LiftCoefftest(testcounter)=sum(InducedLifts);
alphacounter=alphacounter+.05;
testcounter=testcounter+1;

end
alphavec=linspace(-.4,.54,20);
figure
plot(alphavec,LiftCoefftest)
LiftCoeff
%  end