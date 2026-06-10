function [x_out,y_out] = NACA_build(max_cam,max_cam_loc,thick,chord,panels)
%[x_out,y_out] = NACA_build(max_cam,max_cam_loc,thick,chord,panels)
%builds a naca airfoil for x/y values using a generalized equation with
%airfoil specfic paramters.

air_foil_shape_eq = @(x,t,c) (t*c)/.2 .* ( (.2969.*sqrt(x./c)) - (.1260.*(x./c)) - (.3516.*(x./c).^2) + (.2843.*(x./c).^3) - (.1036.*(x./c).^4) );
%air_foil_shape_eq = @(xOc,t) (t*c)/.2 * ( (.2969*sqrt(xOc)) - (.1260*(xOc)) - (.3516*(xOc).^2) + (.2843*(xOc).^3) - (.1036*(xOc).^4) );
%xOc is percentage of chord length, this version may be more useful



x = linspace((chord/(2*panels-1)),chord,(2*panels-1)); %x vector of points along chord diff from what we're calculating in fucntion

y_t = air_foil_shape_eq(x,thick,chord);

y_c_0pc_Eq = @(max_cam,max_cam_loc,x,chord) ( max_cam .* (x./max_cam_loc.^2) .* (2.*max_cam_loc - x./chord) ); %0 <= x =< pc
y_c_pcc_Eq = @(max_cam,max_cam_loc,x,chord) ( max_cam .* ( (chord-x)./ (1-max_cam_loc).^2 ) .* (1 + x./chord - 2.*max_cam_loc) ); %pc <= x <= c                                                )

if max_cam_loc*chord == 0
    y_c = y_c_pcc_Eq(max_cam,max_cam_loc,x,chord);
%this creates a single y_c vector from 0 -> c using the piecewise given
else
x_cut_left_truth = find(x<=max_cam_loc.*chord); %check for x <= pc
x_cut_left = x(1:x_cut_left_truth(end));
x_cut_right = x(x_cut_left_truth(end)+1:end);

y_c = [y_c_0pc_Eq(max_cam,max_cam_loc,x_cut_left,chord) y_c_pcc_Eq(max_cam,max_cam_loc,x_cut_right,chord)];
%this creates a single y_c vector from 0 -> c using the piecewise given
end

dy_c = [y_c(1) diff(y_c)]; %this creates a dy_c vector, incase it's not constant, if it is, well, it'll work just the same
dx = x(2) - x(1); %constant difference

zeta_vec = atan(dy_c./dx);

%just grinding through the final steps.
x_upper = x - y_t.*sin(zeta_vec);
x_lower = x + y_t.*sin(zeta_vec);

y_upper = y_c + y_t.*cos(zeta_vec);
y_lower = y_c - y_t.*cos(zeta_vec);

x_out = [flip(x_lower) x_upper(2:end)];
y_out = [flip(y_lower) y_upper(2:end)];




end