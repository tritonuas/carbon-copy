%%This function gets the planar version of the S compliance matrix 
%%from material properties.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param E1 the youngs modulus in the fiber direction
%@param E2 the youngs modulus in the orthogonal direction
%@param G12 the shear modulus
%@param Nu12 the poissons ratio
%@return S the full compliance matrix
function S_planar = get_S_planar_ortho(mat_props)

E1 = mat_props(1);
E2 = mat_props(2);
G12 = mat_props(3);
Nu12 = mat_props(4);

S_planar = [1/E1, -Nu12/E1, 0;
            -Nu12/E1,    1/E2,   0;
            0,    0,        1/G12];