%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program for finding the supposed optimal taper ratio from a given
%%sweep angle

%This program follows uses an equation from a research paper that can be
%found here:
%https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
%I have reasonable trust for this source as they have seem to have done
%their research, and what logic I picked through seemed to make sense.
%Additionally, I was decently careful to only use equations that contained
%the same assumptions (or fewer assumptions) than I was making.

%@param sweepAngle the angle of the sweep at quarter chord
%@return taperRatio the taper ratio (chordTip/chordRoot)
function taperRatio = getTaper(sweepAngle)
%According to the paper, the optimal taper ratio is given by the equation:
taperRatio = 0.45*exp(-0.0375*sweepAngle); %sweep of quarter chord
%This gives an optimal taper ratio = 0.45 for no sweep (seems reasonable
%and the logic and result compare very closely to what I have read before)
%clarification: lambda = taper ratio in paper