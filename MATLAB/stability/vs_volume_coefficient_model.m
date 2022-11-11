function  C_VT = vs_volume_coefficient_model(vs_cg_distance, vs_area, wing_span, wing_area)
% x_avc
  C_VT = (vs_cg_distance)* vs_area/(wing_area*wing_span);