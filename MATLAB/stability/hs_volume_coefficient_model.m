function C_HT= hs_volume_coefficient_model(hs_cg_distance,hs_area, wing_chord, wing_area)
   C_HT = (hs_cg_distance)*hs_area/(wing_chord* wing_area);
    