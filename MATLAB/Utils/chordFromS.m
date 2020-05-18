%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program for finding the chord given wingspan and wing area

%%This may seem ridiculous, but I am someone who makes many stupid errors,
%%and this should be a good way of reducing them. Especially since I seem
%%to do this computation fairly often.
%%The reason why this is a thing is because there is often a wing area that
%%is required in order to maintain the desired amount of lift (given a lift
%%coefficient, not implying that we are trying to minimize wing area at the
%%expense of a higher lift coefficient. The actual computations are much
%%more complicated and can be found in the MAE 155A and maybe MAE 155B
%%notes). However, we are usually trying to max out our wingspan because,
%%as far as I understand because it's just better (higher AR for a given
%%wing area, and relatively low wing material weight anyway), so the 
%%design should be set with the wingspan being as high as possible given
%%constraints (space in car, on runway, in box, etc.)

%@param s wing area
%@param wingspan
%@return error the vector for the difference between f(x) and true f(x)
function [chord] = chordFromS(s, wingspan)

chord = s/wingspan;
