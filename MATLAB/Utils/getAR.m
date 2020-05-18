%@param wingSpan the total wing span, in meters
%@param S the reference area, in m^2
%@return AR the aspect ratio

%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program for finding the apsect ratio, AR.

%%Source: MAE2 Lecture Slides, and also MAE 155A (both Anderson)
function [AR] = getAR(wingSpan, S)

AR = wingSpan^2./S;