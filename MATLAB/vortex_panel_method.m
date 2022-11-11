
alpha = 0;
chord_length=2;
cruise_speed=20;
nodes = 20
% function [Cl,Cmle] = vortex_panel_method(xu,xl,yu,yl,alpha,chord_length,cruise_speed,nodes)
dtheta = 2*pi/(nodes-1)
i = linspace(1,nodes/2,nodes/2)
xoc=.5*(1-cos(i-.5)*dtheta)













%  end