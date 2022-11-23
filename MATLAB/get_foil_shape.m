%@param max_camber maximum value of the camber divided by 100
%@param max_camber_position Position of maximum camber in percent along
%chord
%@param Air foil thickness divided by 100
%@return xu the vector of x coordinates for the upper surface
%@return xl the vector of x coordinates for the lower surface
%@return yu the vector of y coordinates for the upper surface
%@return yl the vector of y coordinates for the lower surface
%@return yc the vector of y coordinates for the camber line
%@return dycdx the gradient of the camber line with respect to chord
%position

%* could return yt(thickness along chord) and a figure of the chord
% could give .txt file for solidworks

%For function testing:
max_camber=.02; 
max_camber_position=.4;
max_thickness=.12;
close all
%function [xu,xl,yu,yl,yc,dycdx]= get_foil_shape(max_camber,max_camber_position,max_thickness)
x1 = linspace(0,max_camber_position,1000);
x2 = linspace(max_camber_position,1,1000);
yc1=(max_camber/max_camber_position^2)*(2*max_camber_position.*x1-x1.^2);   %Camber from 0 to P
yc2=(max_camber/(1-max_camber_position)^2)*(1-2*max_camber_position+2*max_camber_position.*x2-x2.^2); %Camber from P to 1
yc=[yc1,yc2];         %Camber line vector
dycdx1=2*max_camber/(max_camber_position^2)*(max_camber_position-x1); %Gradient from 0 to P
dycdx2=(2*max_camber/(1-max_camber_position)^2)*(max_camber_position-x2); %Gradient from P to 1
x=[x1,x2];
dycdx=[dycdx1 dycdx2]; %gradient vector dyc/dx
yt = max_thickness/.2*(.2969.*x.^.5-.126.*x-.3516.*x.^2+.2843.*x.^3-.1015.*x.^4);
theta=atan(dycdx);
xc = linspace(0,1,2000); 
xu = xc-yt.*sin(theta);
yu=yc+yt.*cos(theta);
xl = xc-yt.*sin(theta);
yl = yc-yt.*cos(theta);
% For function testing or returning figure
plot(xc,yc)

hold on 
plot(xu,yu)
plot(xl,yl)
xlim([0 1])
ylim([-.5,.5])
%end
