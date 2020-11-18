%%Finds the ABD matrices of a laminate (planar, ortho)
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param thetas the angles of each ply
%@param thicknesses the vector of thicknesses for each lamina
%@param E1 the Young's modulus along the material 1 axis
%@param E2 the Young's modulus along the material 2 axis
%@param G12 the shear modulus
%@param Nu12 the Poisson's ratio
%@return ABD the complete ABD matrix for the laminate
%@return A the A (extensional) matrix for the laminate
%@return B the B (in-plnae/out-of-plane coupling) matrix for the laminate
%@return D the D (bending stiffness) matrix for the laminate
%@return AlphaBetaDelta the inverse of the ABD matrix for the laminate
%@return Alpha the Alpha matrix in the AlphaBetaDelta matrix
%@return Beta the Beta matrix in the AlphaBetaDelta matrix
%@return Delta the Delta matrix in the AlphaBetaDelta matrix
function [ABD, A, B, D, Ex, Ey, Gxy, Nuxy, AlphaBetaDelta, Alpha, Beta,...
    Delta] = laminate_to_ABD_planar_ortho(thetas, rad_or_deg, ...
    thicknesses, E1, E2, G12, Nu12)
                                                
%Preallocation
planar_size = 3;
S = zeros(planar_size, planar_size, length(thetas));
Q = zeros(planar_size, planar_size, length(thetas));
Qbar = zeros(planar_size, planar_size, length(thetas));

for k = 1:length(thetas)
    S(:,:,k) = get_S_planar_ortho(E1, E2, G12, Nu12);
    Q(:,:,k) = inv(S(:,:,k));
    Qbar(:,:,k) = get_Q_bar_planar(Q(:,:,k), thetas(k), rad_or_deg);
end

[ABD, A, B, D, Ex, Ey, Gxy, Nuxy] = QbartoABD(Qbar, thicknesses);
AlphaBetaDelta = inv(ABD);
Alpha = AlphaBetaDelta(1:3, 1:3);
Beta = AlphaBetaDelta(1:3, 4:6);
Delta = AlphaBetaDelta(4:6, 4:6);