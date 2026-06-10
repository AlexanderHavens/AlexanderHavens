%%ASEN 3111 CA3 - MAIN - Alexander Havens
%Collabs: John Davis,Tristan Seely,Professor Farnsworth
%Date: {4/2/2023} m/d/y

theta_vec = [26:.1:58]; %x-axis
mach_vec = [2.2 2.4 2.6 2.8 3 3.2 3.6 4.0 5.5 6 8 10 20];
gamma = 1.4;
beta_strong = zeros(length(mach_vec),length(theta_vec));
beta_weak = zeros(length(mach_vec),length(theta_vec));
for it = 1:length(mach_vec)
    for it1 = 1:length(theta_vec)
beta_strong(it,it1) = ObliqueShockBeta(mach_vec(it),theta_vec(it1),gamma,"Strong");
beta_weak(it,it1) = ObliqueShockBeta(mach_vec(it),theta_vec(it1),gamma,"Strong");
    end
end



