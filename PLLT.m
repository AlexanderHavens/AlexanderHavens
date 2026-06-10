function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
%[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
%inputs
%b = span (ft)
%ao_t = cross-sectional lift slope at the tips (per radian
%a0_r = cross-sectional lift slope at the root (per radian)
%c_r = chord at tip (ft)
% c_r is the chord at the root (in feet)
%aero_t is the zero-lift angle of attack at the tips (in degrees)
%aero_r is the zero-lift angle of attack at the root (in degrees)
%geo_t is the geometric angle of attack at the tips (in degrees)
%geo_r is the geometric angle of attack at the root (in degrees),
%N = number of odd terms
%outputs
%e = span eff. factor
%c_L = coeff of lift
%c_Di = coeff of induced drag

%step 0 - prepare outside of the other stuff
%step 1: we must map all of our physical thingy ma bobbers to theta
%step 2: find theta guesses
%step 3: make 1 mathy operation per theta guess, so say we have 10 theta
%guesses, we need to apply the PLLT equation from class 10 times, each to a
%unique theta.
%step 4: funny square math after it's all put together



theta0_vec = ([1:N]*pi)./(2*N); %this creates all theta0 guesses at the same time
%we need a big ol eqn/mathy thing for EACH theta guess.
N_vec = [1:2*N];
n_odd_vec = N_vec(1:2:end); %we have odd Ns.


chord_theta = c_r - (c_r - c_t)*cos(theta0_vec); %here we do the mapping for a set theta guess
%for the first time round, it's for theta 1. 2nd, theta 3, and so on
alpha_0L_theta = aero_r - (aero_r - aero_t)*cos(theta0_vec);
alpha_geo_theta = geo_r - (geo_r - geo_t)*cos(theta0_vec);
a0_theta = a0_r - (a0_r - a0_t)*cos(theta0_vec);
%og eqn:
%alpha_geo(theta0) =
%[(2*b)/(pi*chord(theta0))]*SUM[An*sin(ntheta0)]+alpha_0L_theta(theta0) +
%SUM[(n*An*(sin(n*theta0)/sin(theta0)))]. we gotta work this into someting
%we can actually use.


%so we got alpha_geo - alpha0L = the other stuff
%we can stick that into a matrix of N length and do linear algebra to
%abracadbra numbers into existence because that's how math works
%the SUMs have to be done as vectors.

A_matrix = zeros(N,N);
left_matrix = (alpha_geo_theta - alpha_0L_theta);

for it1 = 1:length(theta0_vec)
    for it2 = 1:length(n_odd_vec)

    first_right_section = (4*b) ./ ( (a0_theta(it1) ).*(chord_theta(it1) ) );

    first_sum = sin(n_odd_vec(it2)*theta0_vec(it1));
    
    A_section_1 = first_right_section*first_sum;

    A_section_2 = n_odd_vec(it2).* (sin(n_odd_vec(it2).*theta0_vec(it1))./sin(theta0_vec(it1)));

    A_matrix(it1,it2) = A_section_1 + A_section_2;


    end
end

   
        
%ok, we've out our "B" and A matrices, and now we can do funny square math
%to solve them.

An_matrix = A_matrix\left_matrix'; %this matrix has every An term solved.
An_matrix = An_matrix';
%It goes A1 A3 A5 etc etc


S = (c_r + c_t).*(b/2);
AR = b^2/S;

%we can do an actual sum here because we want a scalar
delta = sum(n_odd_vec(2:end).*(An_matrix(2:end)./An_matrix(1)).^2); 

c_L = An_matrix(1)*pi*AR;

c_Di = ((c_L^2)/(pi*AR))*(1+delta);

e = 1/(1+delta);

%IT LIVES MWAHAHAHA
end



