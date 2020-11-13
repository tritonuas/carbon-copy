%%This function gets the S compliance matrix from material properties
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param All the engineering moduli
%@return S the full compliance matrix
function S_planar = get_S_planar_from_full(S)

S_planar = [S(1,1), S(1,2), S(1,6);
            S(2,1), S(2,2), S(2,6);
            S(6,1), S(6,2), S(6,6)];