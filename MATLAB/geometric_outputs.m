function [AR,wing_area,wing_tip_thickness,wing_root_thickness] = geometric_outputs(taper_ratio,wing_chord,wing_span)
wing_area = wing_chord*wing_span;
AR = wing_span^2/wing_area;
wing_tip_thickness = 0.12*(2*wing_chord*taper_ratio)/(1+taper_ratio);
wing_root_thickness = wing_tip_thickness/taper_ratio;
end