function [c_l,c_dw] = DiamondWing(M,alpha,epsilon1,epsilon2)
%[c_l,c_dw] = DiamondWing(M,alpha,epsilon1,epsilon2)
%alpha/episilon in degrees, both must be positive


gamma = 1.4;

%% Into the thick of it
%using epsilon 1
%let's do some basic testing. comparing alpha and epsilion matters a lot,
%as it depends what we get, where, we can't just assume we will always have
%two fans on the top side, though we will get shock->fan on the bottom
%always (who's gonna build a plane with a negative angle of attack????)
  
%% fan/fan shock/fan
if alpha > epsilon1
  %% top half
%expansion fan on top, shock on bottom

%basic stuff before we can get going

[~,~,P01oP1] = flowisentropic(gamma,M); %need our initial pressure ratio
theta_top = deg2rad(alpha - epsilon1); %find our angles for reg 1
theta_bottom = deg2rad(alpha + epsilon2); %angle for reg 4


%our path for top is:
%find m2 using pm function -> use that for pressure for the top side of the
%front half. We will deal with the back half, in a moment
%find m3 using oblique shock equation, use that to find p3. this the bottom
%side of the front half

[~,nu1_top,~] = flowprandtlmeyer(gamma,M,"mach"); %this gives us nu2. 
%the mach it gives us is not valid, we do have a nu2 that simply isn't nu1, 
% like if we had a bent flate plate 
nu2_top = nu1_top + rad2deg(theta_top); % we must solve this for a mach 2. good thing the function can do that for us
mach_2 = flowprandtlmeyer(gamma,nu2_top,"nu"); %this is the mach over the front top "plate"

theta_top_3 = epsilon2 + epsilon1; %this is what we will use going from R2 -> R3

nu3_top = nu2_top + theta_top_3;

mach_3 = flowprandtlmeyer(gamma,nu3_top,"nu"); %mach in region 3, or top back plate

[~,~,P03oP3] = flowisentropic(gamma,mach_3); %pressure ratio for region 3, or the top back plate;
[~,~,P02oP2] = flowisentropic(gamma,mach_2); %this will find our pressure ratios for the front top plate


  %% bottom half - this case, we have shock->fan

beta_bottom = ObliqueShockBeta(M,rad2deg(theta_bottom),gamma,"Weak"); %almost all shocks are weak shocks
%we now have beta
Mn_4 = M*sind(beta_bottom);

[Mn_4,~,P4oP1] = flownormalshock(gamma,Mn_4,"mach"); %pressure ratio for the front bottom plate + the normal mach after our oblique shock
%bottom front plate is "4"

mach_4 = Mn_4/sind(beta_bottom-rad2deg(theta_bottom));


[~,~,P04oP4] = flowisentropic(gamma,mach_4,"mach"); 

%going into reg 5, we have a fan for this case

nu1_bottom = flowprandtlmeyer(gamma,mach_4,"mach"); %find nu1 going into reg 5
nu2_bottom = nu1_bottom + (epsilon1+epsilon2);
mach_5 = flowprandtlmeyer(gamma,nu2_bottom,"nu"); %our mach for reg 5;

[~,~,P05oP5] = flowprandtlmeyer(gamma,mach_5,"mach"); %


%Pressure Ratios

P2oP1 = P01oP1/P02oP2;

P3oP1 = (P2oP1*P02oP2)/P03oP3;

P04oP1 = P04oP4*P4oP1;

P5oP1 = (P04oP1*P01oP1)/P05oP5;

%% shock/fan shock/fan

else %if the above is not true, we will have shock on both parts of the leading edge
%basic stuff before we can get going

[~,~,P01oP1] = flowisentropic(gamma,M); %need our initial pressure ratio
theta_top = deg2rad(alpha - epsilon1); %find our angles for reg 1
theta_bottom = deg2rad(alpha + epsilon2); %angle for reg 4




  %% top half
    %into region 2
    beta_top = ObliqueShockBeta(M,rad2deg(theta_top),gamma,"Weak"); %almost all shocks are weak shocks
    Mn_2 = M*sind(beta_top);
    [Mn_2,~,P2oP1] = flownormalshock(gamma,Mn_2,"mach");
    
    mach_2 = Mn_2/sind(beta_top-rad2deg(theta_top));

    [~,~,P02oP2] = flowisentropic(gamma,mach_2,"mach");

    %going into region 3 now

    [~,nu2_top,~] = flowprandtlmeyer(gamma,mach_2,"mach"); %this gives us nu2. 

    theta_top_3 = epsilon2 + epsilon1; %this is what we will use going from R2 -> R3

    nu3_top = nu2_top + theta_top_3;

    mach_3 = flowprandtlmeyer(gamma,nu3_top,"nu"); %mach in region 3, or top back plate

    [~,~,P03oP3] = flowisentropic(gamma,mach_3); %pressure ratio for region 3, or the top back plate;




  %% bottom half - this case, we have shock->fan
      %this part of the code stays the same as the previous case since we still have a shock->fan
beta_bottom = ObliqueShockBeta(Mn_lower,rad2deg(theta_bottom),gamma,"Weak"); %almost all shocks are weak shocks
%we now have beta
Mn_4 = M*sind(beta_bottom);

[Mn_4,~,P4oP1] = flownormalshock(gamma,Mn_4,"mach"); %pressure ratio for the front bottom plate + the normal mach after our oblique shock
%bottom front plate is "4"

mach_4 = Mn_4/sind(beta_bottom-rad2deg(theta_bottom));


[~,~,P04oP4] = flowisentropic(gamma,mach_4,"mach"); 

%going into reg 5, we have a fan for this case

nu1_bottom = flowprandtlmeyer(gamma,mach_4,"mach"); %find nu1 going into reg 5
nu2_bottom = nu1_bottom + (epsilon1+epsilon2);
mach_5 = flowprandtlmeyer(gamma,nu2_bottom,"nu"); %our mach for reg 5;

[~,~,P05oP5] = flowprandtlmeyer(gamma,mach_5,"mach"); %


%Pressure Ratios


P3oP1 = (P2oP1*P02oP2)*(1/P03oP3);

P04oP1 = P04oP4*P4oP1;

P5oP1 = (P04oP1*P01oP1)/P05oP5;


end


%% Force + c_l + c_dw calcs



% we must assume SOME chord length in order to actually do anything
%for simplicity's sake we're gonna assume c = 1

den = tand(epsilon2)/tand(epsilon1) + 1; %geometery
back_chord = 1/den;

front_chord = 1 - back_chord;

half_thick = front_chord*tand(epsilon1); %this eq comes from 
%deriving the stuff above actually


front_plates_L = half_thick/(cos(epsilon1)); %more geo, "mirrored" diamond
%so both front plates are the same length
back_plates_L = half_thick/cos(epsilon2); %as above, so below

%front sections
Fprime_2 = P2oP1*front_plates_L; %top, negative force 
Fprime_4 = P4oP1*front_plates_L; %bottom, positive force

%back sections
Fprime_3 = P3oP1*back_plates_L; %top, absb
Fprime_5 = P5oP1*back_plates_L; %bottom, absb

Lprime_front = Fprime_4*cos(theta_bottom) - Fprime_2*cos(theta_top);
Lprime_back = Fprime_5*cos(theta_top)  - Fprime_3*cos(theta_bottom);


front_plate_coeff = (2/(gamma*M^2)) * front_plates_L;
back_plate_coeff = (2/(gamma*M^2)) * back_plates_L;

c_l = (front_plate_coeff*Lprime_front) + (back_plate_coeff*Lprime_back);

Dprime_front = Fprime_4*sin(theta_bottom) - Fprime_2*sin(theta_top);
Dprime_back = Fprime_5*sin(theta_top)  - Fprime_3*sin(theta_bottom);


c_dw = (front_plate_coeff*Dprime_front) + (back_plate_coeff*Dprime_back);

end