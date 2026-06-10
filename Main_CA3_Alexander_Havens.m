%%ASEN 3111 CA3 - MAIN - Alexander Havens
%Collabs: John Davis,Tristan Seely,Professor Farnsworth
%Date: {4/2/2023} m/d/y
%P1: Compute lift of a thick AF
%P2: Study effect of thickness on lift
%P3: Study effect of camber on lift
%P4: Implement PLLT
%p5: Use P4 to estimate wing performance
clear
close all
clc
%tic


%% Problem 1 - Compute lift of a Thick Airfoil

%c_l = Vortex_Panel(x_b,y_b,V_inf,alpha)
%x_b,y_b, in ft
%v_inf in ft/s
%alpha in degrees
%go clockwise from TE
%make a naca airfoil construction function done

%% Find cl
% [xb yb] = NACA_build(0,0,.12,1,750); 
Vinf = 10;
alpha = 10;
% c_l = Vortex_Panel(xb,yb,Vinf,alpha);

%The correct answer, calculated with 5000 panels, is 1.2009. it took 1043
%seconds to calculate, pain

%% find panel count required

% Cl_error = 5; %big number just to initalize loop stuff
% cl_actual = 1.2009; %from above using 5000 panels
% panel_count = 10;
% while Cl_error > .01
% [xb,yb] = NACA_build(0,0,.12,1,panel_count);
% c_l = Vortex_Panel(xb,yb,10,10);
% Cl_error = abs((c_l - cl_actual)/cl_actual);
% panel_count = panel_count + 10; %groups of 10 seems like a good jump
% if panel_count == 750
%     break
%     %emergency exit
% end
% end
fprintf("PROBLEM 1----------------------\n")
panel_count = 60; %calculated with ^^, commenting out for performance enhancements
fprintf("The required number of panels to get 1%% relative error is: %.0f\n",panel_count)
fprintf("The c_l value is: %3.4f\n",1.2009)


%% Problem 2 - cl vs alpha plots
%naca 0006
%naca 0012
%naca 0024

[xb_6,yb_6] = NACA_build(0,0,.06,1,60);
[xb_12,yb_12] = NACA_build(0,0,.12,1,60);
[xb_24,yb_24] = NACA_build(0,0,.24,1,60);

cl_6 = zeros(1,length(-10:.1:10));
cl_12 = zeros(1,length(-10:.1:10));
cl_24 = zeros(1,length(-10:.1:10));
loop_count = 1;
for it0 = -10:.1:10
cl_6(loop_count) = Vortex_Panel(xb_6,yb_6,10,it0);
cl_12(loop_count) = Vortex_Panel(xb_12,yb_12,10,it0);
cl_24(loop_count) = Vortex_Panel(xb_24,yb_24,10,it0);
loop_count = loop_count + 1;
end
figure
hold on
plot([-10:.1:10],cl_6,"r")
plot([-10:.1:10],cl_12,"m")
plot([-10:.1:10],cl_24,"b")
legend("NACA 0006","NACA 0012","NACA 0024")
title("cl vs Angle of Attack")
xlabel("Angle of Attack [degrees]")
ylabel("cl [unitless]")
grid minor
hold off


poly_6_coeffs = polyfit([-10:.1:10],cl_6,1);
poly_12_coeffs = polyfit([-10:.1:10],cl_12,1);
poly_24_coeffs = polyfit([-10:.1:10],cl_24,1);

TAT_slope = 2*pi * deg2rad(alpha);

TAT_error_6 = abs( (TAT_slope - poly_6_coeffs(1))/poly_6_coeffs(1) );
TAT_error_12 = abs( (TAT_slope - poly_12_coeffs(1))/poly_12_coeffs(1) );
TAT_error_24 = abs( (TAT_slope - poly_24_coeffs(1))/poly_24_coeffs(1) );


fprintf("Lift slope for NACA 0006: %2.3f\n",poly_6_coeffs(1))
fprintf("Lift slope for NACA 0012: %2.3f\n",poly_12_coeffs(1))
fprintf("Lift slope for NACA 0024: %2.3f\n",poly_24_coeffs(1))
fprintf("----------------------\n")
%figured out by looking at graphs
fprintf("Zero Lift Angle of Attack for NACA 0006: 0\n")
fprintf("Zero Lift Angle of Attack for NACA 0012: 0\n")
fprintf("Zero Lift Angle of Attack for NACA 0024: 0\n")
fprintf("----------------------\n")
fprintf("Error for TAT for NACA 0006: %3.2f%%\n",TAT_error_6)
fprintf("Error for TAT for NACA 0012: %3.2f%%\n",TAT_error_12)
fprintf("Error for TAT for NACA 0024: %3.2f%%\n",TAT_error_24)



%this section ran in 4.25 seconds on my end, so im gonna leave it un
%commented

%% problem 3
%naca 0012
%naca 2412 naca 4412


[xb_012,yb_012]   = NACA_build(0,0,.12,1,panel_count);
[xb_2412,yb_2412] = NACA_build(.02,.4,.12,1,panel_count);
[xb_4412,yb_4412] = NACA_build(.04,.4,.12,1,panel_count);

cl_fat_12 = zeros(1,length([-10:.1:10]));
cl_fat_24 = zeros(1,length([-10:.1:10]));
cl_fat_44 = zeros(1,length([-10:.1:10]));
loop_count = 1;
for it1 = -10:.1:10
cl_fat_12(loop_count) = Vortex_Panel(xb_012,yb_012,10,it1);
cl_fat_24(loop_count) = Vortex_Panel(xb_2412,yb_2412,10,it1);
cl_fat_44(loop_count) = Vortex_Panel(xb_4412,yb_4412,10,it1);
loop_count = loop_count + 1;
end
figure
hold on
plot([-10:.1:10],cl_fat_12,"r")
plot([-10:.1:10],cl_fat_24,"m")
plot([-10:.1:10],cl_fat_44,"b")
legend("NACA 0012","NACA 2412","NACA 4412")
title("cl vs Angle of Attack")
xlabel("Angle of Attack [degrees]")
ylabel("cl [unitless]")
grid minor
hold off

poly_fat_12_coeffs = polyfit([-10:.1:10],cl_fat_12,1);
poly_fat_24_coeffs = polyfit([-10:.1:10],cl_fat_24,1);
poly_fat_44_coeffs = polyfit([-10:.1:10],cl_fat_44,1);

TAT_fat_12_error = abs( (TAT_slope - poly_fat_12_coeffs(1))/poly_fat_12_coeffs(1)  );
TAT_fat_24_error = abs( (TAT_slope - poly_fat_24_coeffs(1))/poly_fat_24_coeffs(1)  );
TAT_fat_44_error = abs( (TAT_slope - poly_fat_44_coeffs(1))/poly_fat_44_coeffs(1)  );

fprintf("PROBLEM 2----------------------\n")

fprintf("Lift slope for NACA 0012: %2.3f\n",poly_fat_12_coeffs(1))
fprintf("Lift slope for NACA 2412: %2.3f\n",poly_fat_24_coeffs(1))
fprintf("Lift slope for NACA 4412: %2.3f\n",poly_fat_44_coeffs(1))
fprintf("----------------------\n")
%looking at graphs
fprintf("Zero Lift AoA NACA 0012: 0 degrees\n")
fprintf("Zero Lift AoA NACA 2412: -2.2 degrees\n")
fprintf("Zero Lift AoA NACA 4412: -4.3 degrees\n")
fprintf("----------------------\n")
fprintf("Error for TAT for NACA 0012: %3.2f%%\n",TAT_fat_12_error)
fprintf("Error for TAT for NACA 2412: %3.2f%%\n",TAT_fat_24_error)
fprintf("Error for TAT for NACA 4412: %3.2f%%\n",TAT_fat_44_error)


%% Problem 4 - Implement PLLT & Verify
%PLLT Implement in function PLLT
% clc
% %Verify PLLT
AR = [4 6 8 10];
c_r = 1; %arbitrary, but set to 1 for ease of use
c_t = [0:.01:1];
a0_t = 2*pi;
a0_r = 2*pi;
aero_t = 0;
aero_r = 0;
geo_t = 5; %arbitrary but not = 0;
geo_r = 5; 
e = zeros(4,length(c_t));
c_L = zeros(4,length(c_t));
c_Di = zeros(4,length(c_t));
delta = zeros(4,length(c_t));
N = 50;
for it1 = 1:4
    b = .5*AR(it1)*(c_r+c_t);
    for it2 = 1:length(c_t)
     [e(it1,it2)] = PLLT(b(it2),a0_t,a0_r,c_t(it2),c_r,aero_t,aero_r,geo_t,geo_r,N);
     delta(it1,it2) = 1/e(it1,it2) - 1;
    end
end

figure
hold on
plot(c_t,delta);
xlabel("c_t/c_r")
ylabel("delta")
legend("AR = 4","AR = 6","AR = 8","AR = 10")
grid minor
title("Taper Ratio vs Delta")
hold off








%% Problem 5 - Cessna Performance
alpha = 3*pi/180;
b = 33.33333333333; %(ft)
c_r = 5.3333333333333; %(ft)
c_t = 3.708333333333; %(ft)
aero_r = (-2.2*pi)/180;% -2.2 degrees in radians
aero_t = 0;
a0_t = poly_fat_12_coeffs(1)*(180/pi);
a0_r = poly_fat_24_coeffs(1)*(180/pi);
geo_r = pi/180 + alpha; %1 degree in radians
geo_t = 0+alpha;

rho = (14.96)*10^(-4);
v_inf = 60*1.68781; %knots to ft/s

%[e_real,cl_real,cdi_real] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,25000);
%Uncomment this out (and wait like 10 mintes lol) to check it


e_real = .9806;
cl_real = .4314;
cdi_real = .0082;

cl_5_error = 5;
e_5_error = 5;
cdi_5_error = 5;
N5_10 = 1;
N5_100 = 1;
N5_1000 = 1;


% 10% error
loop_count = 1;

while cl_5_error > .1 && e_5_error > .1 && cdi_5_error > .1
[e5,cl5,cdi5] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N5_10);

cl_5_error = abs((cl5 - cl_real)/cl_real);
e_5_error = abs( (e5 - e_real)/e_real );
cdi_5_error = abs((cdi5 - cdi_real)/cdi_real);

N5_10 = N5_10 + 1; %step by 1 every time we check
loop_count = loop_count +1;
if loop_count == 200
    disp("N10 Emergency!")
    break %emergency exit
end
end


cl_5_error = 5;
e_5_error = 5;
cdi_5_error = 5;

loop_count = 1;
% 1% error
while cl_5_error > .01 && e_5_error > .01 && cdi_5_error > .01
[e5,cl5,cdi5] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N5_100);

cl_5_error = abs((cl5 - cl_real)/cl_real);
e_5_error = abs( (e5 - e_real)/e_real );
cdi_5_error = abs((cdi5 - cdi_real)/cdi_real);

N5_100 = N5_100 + 1; %step by 1 every time we check
loop_count = loop_count +1;
if loop_count == 200
    disp("N100 Emergency!")
    break %emergency exit
end
end


cl_5_error = 5;
e_5_error = 5;
cdi_5_error = 5;
loop_count = 1;
% .1% error
while cl_5_error > .001 && e_5_error > .001 && cdi_5_error > .001
[e5,cl5,cdi5] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N5_1000);

cl_5_error = abs((cl5 - cl_real)/cl_real);
e_5_error = abs( (e5 - e_real)/e_real );
cdi_5_error = abs((cdi5 - cdi_real)/cdi_real);

N5_1000 = N5_1000 + 1; %step by 1 every time we check
loop_count = loop_count +1;
if loop_count == 200
    disp("N1000 Emergency!")
    break %emergency exit
end
end

%I'm gonna leave this uncommented, it runs virtually instantly for me.
%I did it with c_l/cdi/e instead of lift/drag because the parts of the lift/drag equation
%won't change based on N, so there's point in recycling those calculations
%for no additional benefit. I included the e error for generalization's
%sake

fprintf("PROBLEM 5----------------------\n")

fprintf("Number of Odd terms for 10%% error:%2.0f\n",N5_10)
fprintf("Number of Odd terms for 1%% error:%2.0f\n",N5_100)
fprintf("Number of Odd terms for .1%% error:%2.0f\n",N5_1000)
fprintf("----------------------\n")
%We will use the real value for the actual lift/drag calculations, since
%those are hyper accurate(25,000 is a lot of terms!)

S = (b/2)*(c_t + c_r);
Lift_final = .5*rho*v_inf^2*cl_real*S;
Drag_final = .5*rho*v_inf^2*cdi_real*S;

fprintf("The final lift value is: %3.3f\n",Lift_final)
fprintf("The final drag value is: %3.3f\n",Drag_final)
fprintf("----------------------\n")
%toc
%ive run it like 10 times, takes less than 9 seconds everytime, so
%hopefully it also works on your computers






























