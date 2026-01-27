%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equation 
%
% Author: Tuğba Akman Date: Jan 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = stateEq_highE2(t,y,u1,u2,u3,Tu,E2)

global Eq

dy = zeros(4,1);

u1 = interp1(Tu,u1,t); 
u2 = interp1(Tu,u2,t); 
u3 = interp1(Tu,u3,t); 

T1 = y(1);
ER1 = y(2);
E2ER1 = y(3);
CAF21 = y(4);

if t<=60
  RHS2 = Eq.beta*(1/(1+Eq.d3*CAF21)) - 0.0042*ER1 - (1-u3)*Eq.db*0.5*ER1 + 0.0417*E2ER1; %high
  RHS3 = -0.0125*E2ER1 + (1-u3)*Eq.db*0.5*ER1 - 0.0417*E2ER1;
  % RHS6 = Eq.beta*(1/(1+Eq.d3*CAF22)) - 0.0042*ER2 - Eq.db*0.1*ER2 + 0.0417*E2ER2;%med
  % RHS7 = -0.0125*E2ER2 + Eq.db*0.1*ER2 - 0.0417*E2ER2;
  % RHS10 = Eq.beta*(1/(1+Eq.d3*CAF23)) - 0.0042*ER3 - Eq.db*0.025*ER3 + 0.0417*E2ER3;%low
  % RHS11 = -0.0125*E2ER3 + Eq.db*0.025*ER3 - 0.0417*E2ER3;
else
  RHS2 = Eq.beta*(1/(1+Eq.d3*CAF21)) - 0.0042*ER1 - Eq.db*0*ER1 + 0.0417*E2ER1;
  RHS3 = -0.0125*E2ER1 + Eq.db*0*ER1 - 0.0417*E2ER1;
end

dy(1) = T1*(1-(0.001*T1))*(Eq.k1_hat*(E2ER1/(Eq.alpha1 + E2ER1)) + (1-u2)*Eq.k2*(CAF21/(Eq.alpha2 + CAF21)));
dy(2) = RHS2;
dy(3) = RHS3;
dy(4) = (1-u1)*(Eq.k1_hat/3)*CAF21*(1-(0.001*CAF21))*(T1/(Eq.alpha3 + T1));

end
%%


