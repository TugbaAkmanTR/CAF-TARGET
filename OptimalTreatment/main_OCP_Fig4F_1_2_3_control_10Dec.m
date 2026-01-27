%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate optimal CAF-targeted treatment
%
% Author: Tuğba Akman Date: January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Reid, S. E., Pantaleo, J., Bolivar, P., Bocci, M., Sjölund, J., Morsing,
% M., ... & Pietras, K. (2024). %Cancer-associated fibroblasts rewire
% the estrogen receptor response in luminal breast cancer, enabling
%estrogen independence. Oncogene, 43(15), 1113-1126.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main_OCP_Fig4F_1_2_3_control_10Dec()

close all
clc
format long

set(0, 'defaultaxesfontsize',16)

global Eq

%With CAF2
y_data_with_CAF2_med_dose_mice1 = [0,0,0,0,0,0,0,6.6,12.3,14.8,8.1,13.4,18.4,23.6,32.6,35.11,43.8,50.6,66,117.5,137.9];
y_data_with_CAF2_med_dose_mice2 = [0,0,0,0,0,0,0,0,0,8.3,10.6,8.4,8.7,10.3,10.3,8.83,13.6,15.5,11.1,14.1,14.1];
y_data_with_CAF2_med_dose_mice3 = [0,0,0,0,0,0,0,4.3,5.8,6.6,5.4,6.4,7.2,27.2,12.6,29.04,35.1,43.6,56.9,83.1,100.7];
y_data_with_CAF2_med_dose_mice4 = [0,0,0,0,0,0,9,10.3,18.8,22.1,24.7,27.9,33.3,47.6,82.4,75.48,59.7,75.6,83.6,134.4,232.3];
y_data_with_CAF2_med_dose_mice5 = [0,0,0,0,0,5.1,4.8,6.3,6.9,8.4,8.1,12.7,14.1,29,33.9,35.11,51.7,62.1,73.9,84.8,120.4];

y_data_with_CAF2_med_dose = [y_data_with_CAF2_med_dose_mice1; y_data_with_CAF2_med_dose_mice2; y_data_with_CAF2_med_dose_mice3;...
    y_data_with_CAF2_med_dose_mice4; y_data_with_CAF2_med_dose_mice5];

%% Fig.4F
%Initial conditions
E2 = 0.1;
toll = 1e-5;
init_Tumor = 0.5;
init_ER = 1.300000e-04 ;
init_E2ER = 3.700000e-11;
init_CAF = 0.25*init_Tumor;

% Parameters - 5 mice
Eq.k1_hat = 4.193333e-02;
Eq.alpha1 = 4.446667e-05;
Eq.beta = 6.400000e-03;
Eq.db = 1.066000e-02;

% Weights in the cost functional
Eq.weight1 = 1; % omega_S for CD
Eq.weight2 = 1; % omega_R for CD
Eq.weight3 = 1; % omega_U for CD

Eq.m1=0.001;
Eq.m2=0.001;
Eq.mu1=0.0042;
Eq.mu2=0.0125;
Eq.dub = 0.0417;
Eq.k2 = 0.2046667;
Eq.alpha2 = 21.8413333;
Eq.d3 = 92.7513333;
Eq.alpha3 = 16.2820004;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);


% Initial conditions

initx = [init_Tumor, init_ER, init_E2ER, init_CAF];

% Time discretization
t0=0;                   % initial time
tf = 160;                % final time
Interval=tf*4;           % Number of subintervals
dt = (tf-t0)/Interval; % Temporal increment
Tu=t0:(1*dt):tf;        % ınterpolation points

Eq.time_E2_supply_ends = find(Tu==60);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uncontrolled case - No treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1; % one mouse
u1=zeros(size(Tu))';
u2=zeros(size(Tu))';
u3=zeros(size(Tu))';

% Solve the ODE
[~,X] = ode15s(@(t,x) stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu, initx, options);

MCF7_with_CAF2_med_dose(:,i) = X(:,1);
ER_with_CAF2_med_dose(:,i) = X(:,2);
E2ER_with_CAF2_med_dose(:,i) = X(:,3);
CAF2_with_CAF2_med_dose(:,i) = X(:,4);


%% Get ready for treatment
% Find the first nonzero volume in the data to start treatment
vec = y_data_with_CAF2_med_dose(i,:) ;
first_nonzero_idx = find(vec~= 0, 1, 'first'); % Returns the first index where vec is nonzero
first_nonzero_val = vec(first_nonzero_idx); % Get the value

Eq.threshold(i) = first_nonzero_val;

% Time discretization for intervention
%index_treatment(i) = min(find(MCF7_with_CAF2_med_dose(:,i) >= Eq.threshold(i) ))
index_treatment(i) = min(find(Tu == 45));
t0_treatment = Tu(index_treatment(i));

%index_treatment(i) = first_nonzero_idx;
%t0_treatment  = t_data_med_dose(first_nonzero_idx);
t0=t0_treatment;                   % initial time

tf = 160;                % final time
Interval=tf*4;           % Number of subintervals
dt = (tf-t0)/Interval; % Temporal increment
Tu_treatment(:,i)=t0:(1*dt):tf;        % ınterpolation points

initx_treatment = [MCF7_with_CAF2_med_dose(index_treatment(i),i),ER_with_CAF2_med_dose(index_treatment(i),i),...
    E2ER_with_CAF2_med_dose(index_treatment(i),i),CAF2_with_CAF2_med_dose(index_treatment(i),i)];

clear index_treatment t0_treatment


%%
%% Type I
%%

initx=initx_treatment'; % IC for state
initlambda=[0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)))+0.95;
u2=0*ones(size(Tu_treatment(:,i)));
u3=0*ones(size(Tu_treatment(:,i)));

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment(:,i)));
lambda2=zeros(size(Tu_treatment(:,i)));
lambda3=zeros(size(Tu_treatment(:,i)));
lambda4=zeros(size(Tu_treatment(:,i)));
x1=zeros(size(Tu_treatment(:,i)));
x2=zeros(size(Tu_treatment(:,i)));
x3=zeros(size(Tu_treatment(:,i)));
x4=zeros(size(Tu_treatment(:,i)));

% Control constraints
M1=0;    % Lower bound for control
M2=0.99; % Upper bound for control

nn=1; test = 11;
while(test > toll)

    oldu1 = u1;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_medE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
        flip(Tu_treatment(:,i)), initlambda, options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);

    % Optimality condition
    optCond1 = ((Eq.k1_hat/3).*lambda4.*x4.*(1-Eq.m2*x4).*x1.*(1./(Eq.alpha3 + x1)))/Eq.weight1;
    
    % Project the control
    u1 = min(M2, max(M1, optCond1));

    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
        JCD_for_inside(dd)=cost_u1(Tx,X,u1_convex_for_inside(:,dd));

    end

    [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);

    u1_convex_for_inside(:,pos_min_for_inside);
    u1 = u1_convex_for_inside(:,pos_min_for_inside);

    % Stopping criteria on relative error
    temp1 = norm((oldu1 - u1))/norm((u1));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);
    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp1, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, temp9))))))));

    XCD=X; TxCD=Tx;
    JCD_new(nn)=cost_u1(TxCD,XCD,u1);

    nn=nn+1;

end

MCF7_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,4);

u1_med_dose_control_TypeI(:,i) = u1;

[Ju1,P1,P2]=cost_u1(TxCD,MCF7_with_CAF2_med_dose_control_TypeI,u1,u1);
fprintf('Type I \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju1);


%%
%% Type II
%%

initx=initx_treatment'; % IC for state
initlambda=[0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)));
u2=0*ones(size(Tu_treatment(:,i)))+0.95;
u3=0*ones(size(Tu_treatment(:,i)));

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment(:,i)));
lambda2=zeros(size(Tu_treatment(:,i)));
lambda3=zeros(size(Tu_treatment(:,i)));
lambda4=zeros(size(Tu_treatment(:,i)));
x1=zeros(size(Tu_treatment(:,i)));
x2=zeros(size(Tu_treatment(:,i)));
x3=zeros(size(Tu_treatment(:,i)));
x4=zeros(size(Tu_treatment(:,i)));

% Control constraints
M1=0;    % Lower bound for control
M2=0.99; % Upper bound for control

nn=1; test = 11;
while(test > toll)

    oldu2 = u2;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_medE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
        flip(Tu_treatment(:,i)), initlambda, options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);

    % Optimality condition
    optCond2 = (lambda1.*Eq.k2.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha2+x4)))/Eq.weight2;

    % Project the control
    u2 = min(M2, max(M1, optCond2));
    
    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
        JCD_for_inside(dd)=cost_u2(Tx,X,u2_convex_for_inside(:,dd));

    end

    [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);

    u2_convex_for_inside(:,pos_min_for_inside);
    u2 = u2_convex_for_inside(:,pos_min_for_inside);
    
    % Stopping criteria on relative error
    temp10 = norm((oldu2 - u2))/norm((u2));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);
    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp10, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, temp9))))))));

    XCD=X; TxCD=Tx;
    JCD_new(nn)=cost_u2(TxCD,XCD,u2);

    nn=nn+1;

end

MCF7_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,4);

u2_med_dose_control_TypeII(:,i) = u2;

[Ju2,P1,P2]=cost_u2(TxCD,MCF7_with_CAF2_med_dose_control_TypeII,u2,u2);
fprintf('Type II \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju2);

%%
%% Type III
%%

initx=initx_treatment'; % IC for state
initlambda=[0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)));
u2=0*ones(size(Tu_treatment(:,i)));
u3=0*ones(size(Tu_treatment(:,i)))+0.95;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment(:,i)));
lambda2=zeros(size(Tu_treatment(:,i)));
lambda3=zeros(size(Tu_treatment(:,i)));
lambda4=zeros(size(Tu_treatment(:,i)));
x1=zeros(size(Tu_treatment(:,i)));
x2=zeros(size(Tu_treatment(:,i)));
x3=zeros(size(Tu_treatment(:,i)));
x4=zeros(size(Tu_treatment(:,i)));

% Control constraints
M1=0;    % Lower bound for control
M2=0.99; % Upper bound for control

nn=1; test = 11;
while(test > toll)

    oldu3 = u3;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_medE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
        flip(Tu_treatment(:,i)), initlambda, options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);

    % Optimality condition
    optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;

    % Project the control
    u3 = min(M2, max(M1, optCond3));

    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
        u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
        JCD_for_inside(dd)=cost_u3(Tx,X,u3_convex_for_inside(:,dd));

    end

    [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);

    u3_convex_for_inside(:,pos_min_for_inside);
    u3 = u3_convex_for_inside(:,pos_min_for_inside);

    % Stopping criteria on relative error
    temp11 = norm((oldu3 - u3))/norm((u3));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);

    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, max(temp9,temp11))))))));

    XCD=X; TxCD=Tx;
    JCD_new(nn)=cost_u3(TxCD,XCD,u3);

    nn=nn+1;

end

MCF7_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,4);

u3_med_dose_control_TypeIII(:,i) = u3;

[Ju3,P1,P2]=cost_u3(TxCD,MCF7_with_CAF2_med_dose_control_TypeIII,u3,u3);
fprintf('Type III \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju3);

%%
%% Type I and III
%%

initx=initx_treatment'; % IC for state
initlambda=[0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)))+0.95;
u2=0*ones(size(Tu_treatment(:,i)));
u3=0*ones(size(Tu_treatment(:,i)))+0.95;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment(:,i)));
lambda2=zeros(size(Tu_treatment(:,i)));
lambda3=zeros(size(Tu_treatment(:,i)));
lambda4=zeros(size(Tu_treatment(:,i)));
x1=zeros(size(Tu_treatment(:,i)));
x2=zeros(size(Tu_treatment(:,i)));
x3=zeros(size(Tu_treatment(:,i)));
x4=zeros(size(Tu_treatment(:,i)));

% Control constraints
M1=0;    % Lower bound for control
M2=0.99; % Upper bound for control

nn=1; test = 11;
while(test > toll)

    oldu1 = u1;
    oldu3 = u3;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_medE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
        flip(Tu_treatment(:,i)), initlambda, options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);

    % Optimality condition
    optCond1 = ((Eq.k1_hat/3)*lambda4.*x4.*(1-Eq.m2*x4).*x1.*(1./(Eq.alpha3 + x1)))/Eq.weight1;
    optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;

    % Project the control
    u1 = min(M2, max(M1, optCond1));
    u3 = min(M2, max(M1, optCond3));

    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
        u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
        JCD_for_inside(dd)=cost_u1_u3(Tx,X,u1_convex_for_inside(:,dd),u3_convex_for_inside(:,dd));

    end

    [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);

    u1_convex_for_inside(:,pos_min_for_inside);
    u1 = u1_convex_for_inside(:,pos_min_for_inside);
    u3_convex_for_inside(:,pos_min_for_inside);
    u3 = u3_convex_for_inside(:,pos_min_for_inside);

    % Stopping criteria on relative error
    temp1 = norm((oldu1 - u1))/norm((u1));
    temp11 = norm((oldu3 - u3))/norm((u3));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);
    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp1, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, max(temp9,temp11)))))))));

    XCD=X; TxCD=Tx;
    JCD_new(nn)=cost_u1_u3(TxCD,XCD,u1,u3);

    nn=nn+1;

end

MCF7_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,4);

u1_med_dose_control_TypeI_III(:,i) = u1;
u3_med_dose_control_TypeI_III(:,i) = u3;

[Ju1u3,P1,P2]=cost_u1_u3(TxCD,MCF7_with_CAF2_med_dose_control_TypeI_III,u1,u3);
fprintf('Type I+III \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju1u3);

%%
%% Type II and III
%%

initx=initx_treatment'; % IC for state
initlambda=[0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)));
u2=0*ones(size(Tu_treatment(:,i)))+0.95;
u3=0*ones(size(Tu_treatment(:,i)))+0.95;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment(:,i)));
lambda2=zeros(size(Tu_treatment(:,i)));
lambda3=zeros(size(Tu_treatment(:,i)));
lambda4=zeros(size(Tu_treatment(:,i)));
x1=zeros(size(Tu_treatment(:,i)));
x2=zeros(size(Tu_treatment(:,i)));
x3=zeros(size(Tu_treatment(:,i)));
x4=zeros(size(Tu_treatment(:,i)));

% Control constraints
M1=0;    % Lower bound for control
M2=0.99; % Upper bound for control

nn=1; test = 11;
while(test > toll)

    oldu2 = u2;
    oldu3 = u3;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_medE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
        flip(Tu_treatment(:,i)), initlambda, options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);

    % Optimality condition
    optCond2 = (lambda1.*Eq.k2.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha2+x4)))/Eq.weight2;
    optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;

    % Project the control
    u2 = min(M2, max(M1, optCond2));
    u3 = min(M2, max(M1, optCond3));

    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
        u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
        JCD_for_inside(dd)=cost_u2_u3(Tx,X,u2_convex_for_inside(:,dd),u3_convex_for_inside(:,dd));

    end

    [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);

    u2_convex_for_inside(:,pos_min_for_inside);
    u2 = u2_convex_for_inside(:,pos_min_for_inside);
    u3_convex_for_inside(:,pos_min_for_inside);
    u3 = u3_convex_for_inside(:,pos_min_for_inside);

    % Stopping criteria on relative error
    temp10 = norm((oldu2 - u2))/norm((u2));
    temp11 = norm((oldu3 - u3))/norm((u3));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);
    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp11, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, max(temp9,temp10)))))))));

    %test=temp1;

    XCD=X; TxCD=Tx;
    JCD_new(nn)=cost_u2_u3(TxCD,XCD,u2,u3);

    nn=nn+1;

end

MCF7_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,4);

u2_med_dose_control_TypeII_III(:,i) = u2;
u3_med_dose_control_TypeII_III(:,i) = u3;

[Ju2u3,P1,P2]=cost_u2_u3(TxCD,MCF7_with_CAF2_med_dose_control_TypeII_III,u2,u3);
fprintf('Type II+III \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju2u3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = zeros(1,7);
figure(930)
subplot(3,2,1)
h(1)=plot(Tu_treatment(:,i),u1_med_dose_control_TypeI_III(:,i),'k-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
h(2)=plot(Tu_treatment(:,i),u3_med_dose_control_TypeI_III(:,i),'r-','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
grid on
legend('u_1','u_3','Location','NorthWest')
xlabel('t')
ylabel('u_1(t), u_3(t)')
xlim([24,tf])
ylim([0,1])
sgtitle('Medium dose of E2 - Control functions - u_1(t) (circle), u_2(t) (diamond) , u_3(t) (square)');

subplot(3,2,2)
h(3)=plot(Tu_treatment(:,i),u2_med_dose_control_TypeII_III(:,i),'b-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
h(4)=plot(Tu_treatment(:,i),u3_med_dose_control_TypeII_III(:,i),'r-','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
grid on
legend('u_2','u_3','Location','NorthWest')
xlabel('t')
ylabel('u_2(t), u_3(t)')
xlim([24,tf])
ylim([0,1])

subplot(3,2,3)
h(5)=plot(Tu_treatment(:,i),u1_med_dose_control_TypeI(:,i),'k-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
grid on
xlabel('t')
ylabel('u_1(t)')
xlim([24,tf])
ylim([0,1])

subplot(3,2,4)
h(6)=plot(Tu_treatment(:,i),u2_med_dose_control_TypeII(:,i),'b-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
grid on
xlabel('t')
ylabel('u_2(t)')
xlim([24,tf])
ylim([0,1])

subplot(3,2,5)
h(7)=plot(Tu_treatment(:,i),u3_med_dose_control_TypeIII(:,i),'r-','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
grid on
xlabel('t')
ylabel('u_3(t)')
xlim([24,tf])
ylim([0,1])
sgtitle('Control functions - Medium dose E2 (0.5mg) (with CAF)','fontsize',18)

%%
Percentage_change_TypeI(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))* (MCF7_with_CAF2_med_dose(:,i)-MCF7_with_CAF2_med_dose_control_TypeI(:,i))
Percentage_change_TypeII(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))* (MCF7_with_CAF2_med_dose(:,i)-MCF7_with_CAF2_med_dose_control_TypeII(:,i))
Percentage_change_TypeIII(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))*(MCF7_with_CAF2_med_dose(:,i)- MCF7_with_CAF2_med_dose_control_TypeIII(:,i))
Percentage_change_TypeI_III(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))*(MCF7_with_CAF2_med_dose(:,i)- MCF7_with_CAF2_med_dose_control_TypeI_III(:,i))
Percentage_change_TypeII_III(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))* (MCF7_with_CAF2_med_dose(:,i)-MCF7_with_CAF2_med_dose_control_TypeII_III(:,i))

%%
%%Treatment
figure(412)
subplot(2,2,1)
plot(Tu,MCF7_with_CAF2_med_dose(:,1), '-b','LineWidth',5);
hold on
plot(Tu_treatment(:,1),MCF7_with_CAF2_med_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),MCF7_with_CAF2_med_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
hold on
plot(Tu_treatment(:,1),MCF7_with_CAF2_med_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
hold on
plot(Tu_treatment(:,1),MCF7_with_CAF2_med_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
hold on
plot(Tu_treatment(:,1),MCF7_with_CAF2_med_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
hold on
xline(Tu_treatment(1,1), '--r', 'Treatment starts', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'top', 'LineWidth', 2);
xline(60, '-.b', 'E2 supply ends', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'middle', 'LineWidth', 1);
legend('No treatment','Treatment - Type I','Treatment - Type II','Treatment - Type III','Treatment - Type I+III', 'Treatment - Type II+III','location', 'northwest',  'FontSize', 10)
sgtitle('Optimal treatment - Medium dose E2 (0.5mg)','fontsize',18)
xlabel('t','fontweight','normal','fontsize',18)
ylabel('T(t)','fontweight','normal','fontsize',18)
grid on
xlim([0,160])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

subplot(2,2,2)
plot(Tu,ER_with_CAF2_med_dose(:,1), '-b','LineWidth',5);
hold on
plot(Tu_treatment(:,1),ER_with_CAF2_med_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),ER_with_CAF2_med_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),ER_with_CAF2_med_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),ER_with_CAF2_med_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),ER_with_CAF2_med_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
xline(Tu_treatment(1,1), '--r', 'LineWidth', 2);
xline(60, '-.b', 'LineWidth', 1);
xlabel('t','fontweight','normal','fontsize',18)
ylabel('ER(t)','fontweight','normal','fontsize',18)
grid on
xlim([0,160])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

subplot(2,2,3)
plot(Tu,E2ER_with_CAF2_med_dose(:,1), '-b','LineWidth',5);
hold on
plot(Tu_treatment(:,1),E2ER_with_CAF2_med_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),E2ER_with_CAF2_med_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),E2ER_with_CAF2_med_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),E2ER_with_CAF2_med_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),E2ER_with_CAF2_med_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
xline(Tu_treatment(1,1), '--r', 'LineWidth', 2);
xline(60, '-.b', 'LineWidth', 1);
xlabel('t','fontweight','normal','fontsize',18)
ylabel('E2ER(t)','fontweight','normal','fontsize',18)
grid on
xlim([0,160])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

subplot(2,2,4)
plot(Tu,CAF2_with_CAF2_med_dose(:,1), '-b','LineWidth',5);
hold on
plot(Tu_treatment(:,1),CAF2_with_CAF2_med_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),CAF2_with_CAF2_med_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),CAF2_with_CAF2_med_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),CAF2_with_CAF2_med_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
hold on
plot(Tu_treatment(:,1),CAF2_with_CAF2_med_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
xline(Tu_treatment(1,1), '--r', 'LineWidth', 2);
xline(60, '-.b', 'LineWidth', 1);
xlabel('t','fontweight','normal','fontsize',18)
ylabel('CAF(t)','fontweight','normal','fontsize',18)
grid on
xlim([0,160])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

end

%Sub-functions to compute J

function [Ju1u3,P1,P2]=cost_u1_u3(Tx,X,u1,u3)
global Eq

P1 = trapz(Tx,X(:,1));
P2 = Eq.weight1*0.5*trapz(Tx,u1.^2)+ Eq.weight3*0.5*trapz(Tx,u3.^2);
Ju1u3=P1 + P2;
end

function [Ju2u3,P1,P2]=cost_u2_u3(Tx,X,u2,u3)
global Eq

P1 = trapz(Tx,X(:,1));
P2 = Eq.weight2*0.5*trapz(Tx,u2.^2)+ Eq.weight3*0.5*trapz(Tx,u3.^2);
Ju2u3=P1 + P2;
end

function [Ju1,P1,P2]=cost_u1(Tx,X,u1,u2)
global Eq

P1 = trapz(Tx,X(:,1));
P2 = Eq.weight1*0.5*trapz(Tx,u1.^2);
Ju1=P1 + P2;
end

function [Ju2,P1,P2]=cost_u2(Tx,X,u2,u1)
global Eq

P1 = trapz(Tx,X(:,1));
P2 = Eq.weight2*0.5*trapz(Tx,u2.^2);
Ju2=P1 + P2;
end

function [Ju3,P1,P2]=cost_u3(Tx,X,u3,u1)
global Eq

P1 = trapz(Tx,X(:,1));
P2 = Eq.weight3*0.5*trapz(Tx,u3.^2);
Ju3=P1 + P2;
end
