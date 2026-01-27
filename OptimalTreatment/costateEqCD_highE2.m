%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CoState equation 
%
% Author: Tuğba Akman Date: Jan 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dp = costateEqCD_highE2(t,p,u1,u2,u3,Tu,x1,x2,x3,x4,xt,E2)

global Eq

dp = zeros(4,1);

x1 = interp1(xt,x1,t); 
x2 = interp1(xt,x2,t);
x3 = interp1(xt,x3,t); 
x4 = interp1(xt,x4,t); 
 
Tumor = x1;
ER = x2;
E2ER = x3;
CAF2 = x4;

u1 = interp1(Tu,u1,t); 
u2 = interp1(Tu,u2,t); 
u3 = interp1(Tu,u3,t); 

R1S = (1-Eq.m1*Tumor)*(Eq.k1_hat*(E2ER/(Eq.alpha1+E2ER)) + (1-u2)*Eq.k2*(CAF2/(Eq.alpha2+CAF2)));
R1R = 0;
R1I = Tumor*(1-Eq.m1*Tumor)*Eq.k1_hat*(Eq.alpha1/((Eq.alpha1+E2ER)^2));
R1D1 = Tumor*(1-Eq.m1*Tumor)*Eq.k2*(1-u2)*(Eq.alpha2/((Eq.alpha2+CAF2)^2));

R2S = 0;
R2R = -Eq.mu1 - Eq.db*0.5*(1-u3);
R2I = Eq.dub;
R2D1 = -Eq.beta*Eq.d3*(1/((1 + Eq.d3*CAF2)^2));

R3S = 0;
R3R = Eq.db*0.5*(1-u3);
R3I = -Eq.mu2 + Eq.dub;
R3D1 = 0;

R4S = (1-u1)*Eq.k1_hat*CAF2*(1-Eq.m2*CAF2)*(Eq.alpha3/((Eq.alpha3+Tumor)^2));
R4R = 0;
R4I = 0;
R4D1 = (1-u1)*(Eq.k1_hat/3)*(Tumor/((Eq.alpha3+Tumor)^2))*(1-2*Eq.m2*CAF2);

dp(1) = -p(1).*R1S - p(2).*R2S - p(3).*R3S - p(4).*R4S -1;

dp(2) = -p(1).*R1R - p(2).*R2R - p(3).*R3R - p(4).*R4R;

dp(3) = -p(1).*R1I - p(2).*R2I - p(3).*R3I - p(4).*R4I;
  
dp(4) = -p(1).*R1D1 - p(2).*R2D1 - p(3).*R3D1 - p(4).*R4D1;
 


end
