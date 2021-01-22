function [Ytotal] = HoSparNonSymmAirfoil(x,a4,c,NACA)
%For the Triton UAS Team
% by: Ryan Beneduce
%   For the determination of max height of spar at desired location
% Source equations found at: airfoiltools.com/airfoil/naca4digit
% Non-Symmetrical airfoils only
% for symmetrical airfoils see other function
%inputs of function include: 
%'x'- location along chord length (range from 0-c)
%'a4'- value for closed tail edge= -.1036, normal tail edge= -.1015
%'c'-length of chord
% NACA= [ enter in 4-digit airfoil with spaces between]
% example: NACA=[2 0 14] %make sure not to have a space btwn last 2 digits

M=NACA(1)/100; %max camber
P=NACA(2)/10; %position of max camber
T=NACA(3)/100; %thickness
%Yt calculation
ao=.2969;
a1=-.126;
a2=-.3516;
a3=.2843;
yt= T*5*c*(ao*(x/c)^.5+a1*x/c+a2*(x/c)^2+a3*(x/c)^3+a4*(x/c)^4);
%Yc calculation
if x<(P*c)
    yc=((M*x)/P^2)*(2*P-(x/c)); %camber
    dycdx=(2*M/P^2)*(P-(x/c)); %gradient
elseif x>=(P*c)
    yc=(M*(c-x))/(1-P)^2*(1+(x/c)-2*P); %camber
    dycdx=(2*M/(1-P)^2)*(P-(x/c)); %gradient
end 
%theta and upper/lower surface calcs
theta=atan(dycdx);
xu=x-yt*sin(theta); %upper bound x
xl=x+yt*sin(theta);%lower bound x
yu=yc+yt*cos(theta); %upper bound y
yl=yc-yt*cos(theta); %lower bound y
Ytotal=yu-yl;
fprintf('At position x= %d the max spar height is = %d',x,Ytotal)
