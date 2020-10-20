%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program is a MATLAB version of Cd = Cdi + Cd0

%%It is just saying that the dra coefficient is the sum of the 
%%drag induced by lift, and the drag not induced by lift.
%%You can look at each individual method to see how those are calculated.

%%Source: MAE 2 lectures slides

function [cd] = getDragCoeff(cd0,cdi)

cd = cd0 + cdi;


