%%This function gets the S compliance matrix from material properties
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param All the engineering moduli
%@return S the full compliance matrix
function S = get_S(Ex, Ey, Ez, Gyz, Gxz, Gxy, Nuxy, Nuxz, Nuyz, Nuyx, ...
    Nuzx, Nuzy, Etayz_x, Etaxz_x, Etaxy_x, Etayz_y, Etaxz_y, Etaxy_y, ...
    Etayz_z, Etaxz_z ,Etaxy_z, Etax_yz, Etay_yz, Etaz_yz, Etax_xz,  ...
    Etay_xz, Etaz_xz, Etax_xy, Etay_xy, Etaz_xy, Muxz_yz, Muxy_yz, ...
    Muyz_xz, Muxy_xz, Muyz_xy, Muxz_xy)

S = [1/Ex, -Nuyx/Ey, -Nuzx/Ez, Etayz_x/Gyz, Etaxz_x/Gxz, Etaxy_x/Gxy;
     -Nuxy/Ex, 1/Ey, -Nuzy/Ez, Etayz_y/Gyz, Etaxz_y/Gxz, Etaxy_y/Gxy;
     -Nuxz/Ex, -Nuyz/Ey, 1/Ez, Etayz_z/Gyz, Etaxz_z/Gxz, Etaxy_z/Gxy;
     Etax_yz/Ex, Etay_yz/Ey, Etaz_yz/Ez, 1/Gyz, Muxz_yz/Gxz, Muxy_yz/Gxy;
     Etax_xz/Ex, Etay_xz/Ey, Etaz_xz/Ez, Muyz_xz/Gyz, 1/Gxz, Muxy_xz/Gxy;
     Etax_xy/Ex, Etay_xy/Ey, Etaz_xy/Ez, Muyz_xy/Gyz, Muxz_xy/Gxz, 1/Gxy];