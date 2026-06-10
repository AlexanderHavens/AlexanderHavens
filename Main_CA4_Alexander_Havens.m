%%ASEN 3111 CA4 - MAIN - Alexander Havens
%Collabs: John Davis, Brad Lindworth
%Date: {4/2/2023} m/d/y
clear
clc
close all

%% Problem 1
theta_vec = [0:.1:58]; %x-axis
mach_vec = [1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.4 2.6 2.8 3 3.2 3.6 4.0 5 6 8 10 20];
gamma = 1.4;


for it = 1:length(mach_vec) %do the math for strong shocks

    for it1 = 1:length(theta_vec)



beta_strong_temp = ObliqueShockBeta(mach_vec(it),theta_vec(it1),gamma,"Strong");
%test to see if real, if not, we don't care
        if isreal(beta_strong_temp)
            beta_strong{1,it}(1,it1) = beta_strong_temp;
        else
            break
        end

    end


end




for it = 1:length(mach_vec) %math for weak shocks

    for it1 = 1:length(theta_vec)


beta_weak_temp = ObliqueShockBeta(mach_vec(it),theta_vec(it1),gamma,"Weak");
%test to see if real, if not, we don't care
        if isreal(beta_weak_temp)
        beta_weak{1,it}(1,it1) = beta_weak_temp;
        else
            break
        end


    end


end

%%plotting
fig_beta = figure;
hold on
format_vec = ["b" "g" "r" "c" "m" "k" "b:" "g:" "r:" "c:" "m:" "b--" "g--" "r--" "c--" "m--" "k--" "b-." "g-." "r-." "c-." "m-." "k-." "b" "g" "r" "c" "m"];
for it = 1:28

    theta_vec_temp_strong = theta_vec(1:length(beta_strong{it}) );
    theta_vec_temp_weak = theta_vec(1:length(beta_weak{it}));


   p_beta_strong(it) = plot(theta_vec_temp_strong,beta_strong{it},format_vec{it});
   p_beta_weak(it) = plot(theta_vec_temp_weak,beta_weak{it},format_vec{it});
end
lgd = legend([p_beta_strong],"1.05","1.1","1.15","1.2","1.25","1.3","1.35","1.4","1.45","1.5","1.6","1.7","1.8","1.9","2.0","2.2","2.4","2.6","2.8","3","3.2","3.6","4","5","6","8","10","20","Location","southeast");
lgd.NumColumns = 4;
lgd.Position(3:4) = lgd.Position(3:4)*.75;
lgd.Position(1:2) = lgd.Position(1:2)*(1/.68);
lgd.Position(2) = lgd.Position(2)*.89;
title(lgd,"Mach # =")
legend("boxoff")
grid minor
title("Theta & Mach vs Beta")
xlabel("Deflection Angle [Deg]")
ylabel("Beta [deg]")

%the legend will show up over some of the lines, nothing i can do it's a
%big legend im afraid, but i've done my best to make it legible as possible
%with the amount of information on the graph, expanding the window will
%make it more legible if you'd like to


%% Problem 2 %%
M = 3;
alpha = 10;
e1 = 7.5;
e2 = 5;

[c_l,c_dw] = DiamondWing(M,alpha,e1,e2);


fprintf("c_l:%2.3f \n",c_l)
fprintf("C_dw:%2.3f \n",c_dw)


%% Problem 3 %%
M_vec = [2:5];
alpha_3 = 10;
e1_3 = 7.5;
e2_3 = 5;

for itp3 = 1:4
  [c_l3(itp3),c_dw3(itp3)] = DiamondWing(M_vec(itp3),alpha_3,e1_3,e2_3);
end

den = tand(e2_3)/tand(e1_3) + 1; %geometery
back_chord = 1/den;

front_chord = 1 - back_chord;

half_thick = front_chord*tand(e1_3); %this eq comes from
%deriving the stuff above actually


front_plates_L = half_thick/(cos(e1_3)); %more geo, "mirrored" diamond
%so both front plates are the same length
back_plates_L = half_thick/cos(e2_3); %as above, so below

%gl and gu will be the same value, since the components that make up each
%will be the same

gl = front_plates_L + back_plates_L;

cl_flat = 4*deg2rad(alpha_3)./ sqrt(M_vec.^2 - 1);
cdw_flat = (2./sqrt(M_vec.^2 -1)).*(2.*deg2rad(alpha_3)^2+2.*gl^2);

hold on
figure
plot(M_vec,c_l3,M_vec,cl_flat);
title("c_l Analysis")
xlabel("Mach")
ylabel("c_l")
legend("S-E Theory","Linearized")
grid minor
hold off

hold on
figure
plot(M_vec,c_dw3,M_vec,cdw_flat);
title("c_w Analysis")
xlabel("Mach")
ylabel("c_w")
legend("S-E Theory","Linearized")
grid minor
hold off









