%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program is for finding the induced drag

%%This is equation is the induced drag lift coefficient equation.

function [cdi] = getInducedDragCoeff(k,cl)

cdi = k*cl^2;