%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate optimal anti-hormonal treatment
%
% Author: Tuğba Akman Date: July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Obesity-Activated Adipose-Derived Stromal Cells Promote
% Breast Cancer Growth and Invasion
% Neoplasia (2018) 20, 1161-1174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_OCP()

close all
clear
clc
format long

set(0, 'defaultaxesfontsize',16)

global Eq

%Data
%t_data_high_dose = [17,24,32,39,45,52,62,67,73,80,87,95,102,108,115,119];
t_data_low_dose = [17,24,32,39,45,52,62,67,73,80,87,95,102,108,115,119,124];
%t_data_med_dose = [17,24,32,39,45,52,62,67,73,80,87,95,102,108,115,123,130,137,144,151,158];



% High dose E2 - Fig4D
%MCF7 alone
y_data_alone_high_dose_mice1 = [0,0,0,0,0,0,6.1,6.6,13.5,31.8,79.8,205.1,394.2,601.8,899.7,899.7];
y_data_alone_high_dose_mice2 = [0,0,0,0,0,0,2.9,6.4,6.3,9.9,18.8,42.6,98.5,174.4,229.3,264.3];
y_data_alone_high_dose_mice3 = [0,0,0,0,0,0,7.2,7.2,7.5,6.9,8.5,8.8,11.1,15.5,30.2,30.2];
y_data_alone_high_dose_mice4 = [0,0,0,0,0,5.8,21.2,23.1,25,64.2,123.2,99.4,144.5,197.8,204.5,204.5];
y_data_alone_high_dose_mice5 = [0,0,0,0,5.8,6.6,9.6,8.5,9.6,20.4,50.2,96.8,158.6,173.8,102.1,122.5];

%With CAF2
y_data_with_CAF2_high_dose_mice1 = [0,0,0,0,3.8,4.4,4.4,3.8,4,5.8,8.5,14.6,21.8,52.7,62.7,88.5];
y_data_with_CAF2_high_dose_mice2 = [0,0,0,0,0,0,0,4.4,5.1,8.7,9.5,11.9,12.8,14.1,8.8,8.8];
y_data_with_CAF2_high_dose_mice3 = [0,0,0,4,4.2,5.1,9.9,8.8,14.1,14.4,20.4,23.7,72.3,162,239.8,263];
y_data_with_CAF2_high_dose_mice4 = [0,4.2,4.4,9.6,9.9,9.9,13.2,16.6,19.4,27.2,39.3,60.4,177.8,443.5,663.5,692.4];
y_data_with_CAF2_high_dose_mice5 = [0,0,0,4.4,4.8,5.1,5.8,7.5,9.9,16,24.5,32.2,62.8,106.1,154.8,206.4];

y_data_with_CAF2_high_dose = [y_data_with_CAF2_high_dose_mice1; y_data_with_CAF2_high_dose_mice2; y_data_with_CAF2_high_dose_mice3;...
    y_data_with_CAF2_high_dose_mice4; y_data_with_CAF2_high_dose_mice5];

% % Low dose E2 - Fig4E
% %MCF7 alone
% y_data_alone_low_dose_mice1 = [0,0,0,0,8.2,9.9,15.2,15,9.9,19.6,21.1,16.8,22,26.9,25,43.5,54];
% y_data_alone_low_dose_mice2 = [0,0,0,0,0,0,0,0,0,4.6,6.1,6.6,9.2,17.7,33.9,68.9,79.6];
% y_data_alone_low_dose_mice3 = [0,0,0,0,0,0,6.6,10.3,13,21.4,29.6,61.5,74.8,100.5,141.3,187.9,196.8];
% y_data_alone_low_dose_mice4 = [0,0,0,0,0,0,4,5.1,5.6,7.2,8.8,11.1,20.4,25,47.7,83.9,91.8];
% y_data_alone_low_dose_mice5 = [0,0,0,0,0,0,5.3,8,7.2,8.5,9.9,9.9,14.1,16.5,25.8,38.5,40.8];
% 
% %With CAF2
y_data_with_CAF2_low_dose_mice1 = [0,0,0,0,0,0,0,6.3,7.2,15.1,21.2,31.5,88.7,107.6,116.2,116.2,123.9];
y_data_with_CAF2_low_dose_mice2 = [0,0,0,0,0,0,8,13.2,20.6,33,37,62.5,82.4,120.6,131.3,184.7,182.2];
y_data_with_CAF2_low_dose_mice3 = [0,0,0,0,0,0,0,8.7,14.2,24.1,66.7,96.8,136.1,294.4,294.9,487.5,491.7];
y_data_with_CAF2_low_dose_mice4 = [0,0,0,0,0,0,0,4.2,5.3,9.6,11.1,20.5,21.7,34.4,46.2,41.6,49.1];
y_data_with_CAF2_low_dose_mice5 = [0,0,0,0,0,0,0,0,0,13.6,19.4,26.5,49.8,96.6,169.4,183.9,200];

y_data_with_CAF2_low_dose = [y_data_with_CAF2_low_dose_mice1; y_data_with_CAF2_low_dose_mice2; y_data_with_CAF2_low_dose_mice3;...
    y_data_with_CAF2_low_dose_mice4; y_data_with_CAF2_low_dose_mice5];
% 
% % Med dose E2 - Fig4f
% %MCF7 alone
% y_data_alone_med_dose_mice1 = [0,0,0,0,0,4.8,5.1,5.1,8.8,10.1,7.8,8.5,9.2,11,12.7,14.08,14.8,16.8,18.2,35.8,21.7];
% y_data_alone_med_dose_mice2 = [0,0,0,0,0,0,0,0,0,5.3,7.5,7.8,10.1,12.3,13.9,12.38,11,14.9,19.8,24.4,31.5];
% y_data_alone_med_dose_mice3 = [0,0,0,0,0,0,0,0,0,12.7,11.7,12.6,14.8,16.7,24.1,24.12,25.7,25.7,30.2,39.3,67.4];
% y_data_alone_med_dose_mice4 = [0,0,0,0,0,0,0,0,0,4.2,7.8,9.8,13.9,14.4,9.2,10.61,12.6,12.6,12.6,30.1,17.6];
% y_data_alone_med_dose_mice5 = [0,0,0,0,0,0,0,0,0,0,4.4,5.1,5.8,6.1,6.9,9.14,6.9,9.2,10.6,23.1,26.6];
% 
% %With CAF2
% y_data_with_CAF2_med_dose_mice1 = [0,0,0,0,0,0,0,6.6,12.3,14.8,8.1,13.4,18.4,23.6,32.6,35.11,43.8,50.6,66,117.5,137.9];
% y_data_with_CAF2_med_dose_mice2 = [0,0,0,0,0,0,0,0,0,8.3,10.6,8.4,8.7,10.3,10.3,8.83,13.6,15.5,11.1,14.1,14.1];
% y_data_with_CAF2_med_dose_mice3 = [0,0,0,0,0,0,0,4.3,5.8,6.6,5.4,6.4,7.2,27.2,12.6,29.04,35.1,43.6,56.9,83.1,100.7];
% y_data_with_CAF2_med_dose_mice4 = [0,0,0,0,0,0,9,10.3,18.8,22.1,24.7,27.9,33.3,47.6,82.4,75.48,59.7,75.6,83.6,134.4,232.3];
% y_data_with_CAF2_med_dose_mice5 = [0,0,0,0,0,5.1,4.8,6.3,6.9,8.4,8.1,12.7,14.1,29,33.9,35.11,51.7,62.1,73.9,84.8,120.4];
% 
% y_data_with_CAF2_med_dose = [y_data_with_CAF2_med_dose_mice1; y_data_with_CAF2_med_dose_mice2; y_data_with_CAF2_med_dose_mice3;...
%     y_data_with_CAF2_med_dose_mice4; y_data_with_CAF2_med_dose_mice5];

%% Fig.4D - high estrogen
%Initial conditions
E2 = 0.025;

toll = 1e-5;
init_Tumor = 0.5;
init_ER = 2.2e-05;
init_E2ER = 5.4e-08;
init_CAF = 1.5;

% Parameters - 5 mice
Eq.k1_hat = 0.04413;
Eq.alpha1 = 0.001147;
Eq.beta = 0.0072;
Eq.db = 0.23;

% Weights in the cost functional
Eq.weight1 = 1; % omega_S for CD
Eq.weight2 = 1; % omega_R for CD
Eq.weight3 = 1; % omega_U for CD

Eq.m1=0.001;
Eq.m2=0.001;
Eq.mu1=0.0042;
Eq.mu2=0.0125;
Eq.dub = 0.0417;
Eq.k2 = 0.1651 ;
Eq.alpha2 = 13.02;
Eq.d3 = 69.77;
Eq.alpha3 = 2.778e+01;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

for i=1:1 % 5 mice
    
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
    
    u1=zeros(size(Tu))';
    u2=zeros(size(Tu))';
    u3=zeros(size(Tu))';
    
    % Solve the ODE
    [~,X] = ode15s(@(t,x) stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu, initx, options);
    
    MCF7_with_CAF2_low_dose(:,i) = X(:,1);
    ER_with_CAF2_low_dose(:,i) = X(:,2);
    E2ER_with_CAF2_low_dose(:,i) = X(:,3);
    CAF2_with_CAF2_low_dose(:,i) = X(:,4);
    
    
    %% Get ready for treatment
    % Find the first nonzero volume in the data to start treatment
    vec = y_data_with_CAF2_low_dose(2,:) ;
    first_nonzero_idx = find(vec~= 0, 1, 'first'); % Returns the first index where vec is nonzero
    first_nonzero_val = vec(first_nonzero_idx); % Get the value
    
    Eq.threshold(i) = first_nonzero_val;
    
    % Time discretization for intervention
    %index_treatment(i) = min(find(MCF7_with_CAF2_low_dose(:,i) >= Eq.threshold(i) ))
    index_treatment(i) = min(find(Tu == 45))
    t0_treatment = Tu(index_treatment(i))
    
    %index_treatment(i) = first_nonzero_idx;
    %t0_treatment  = t_data_low_dose(first_nonzero_idx);
    t0=t0_treatment;                   % initial time
    
    tf = 160;                % final time
    Interval=tf*4;           % Number of subintervals
    dt = (tf-t0)/Interval; % Temporal increment
    Tu_treatment(:,i)=t0:(1*dt):tf;        % ınterpolation points
    
    initx_treatment = [MCF7_with_CAF2_low_dose(index_treatment(i),i),ER_with_CAF2_low_dose(index_treatment(i),i),...
        E2ER_with_CAF2_low_dose(index_treatment(i),i),CAF2_with_CAF2_low_dose(index_treatment(i),i)];
    
    clear index_treatment t0_treatment
    
    
    %%
    %% Type I
    %%
    
    initx=initx_treatment'; % IC for state
    initlambda=[0,0,0,0]';                  % FC for adjoint
    size(Tu_treatment(:,i))
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Solve OCP for Type I + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nn=1; test = 11;
    while(test > toll)
        
        oldu1 = u1;
        %oldu2 = u2;
        %oldu3 = u3;
        oldx1=x1;
        oldx2=x2;
        oldx3=x3;
        oldx4=x4;
        oldlambda1=lambda1;
        oldlambda2=lambda2;
        oldlambda3=lambda3;
        oldlambda4=lambda4;
        
        % solve the state eqn forward
        [Tx,X] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);
        
        % solve the adjoint eqn backward
        x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
        [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_lowE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
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
        %optCond2 = (lambda1.*Eq.k3.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha3+x4)))/Eq.weight2;
        %optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;
        
        % Project the control
        u1 = min(M2, max(M1, optCond1));
        %u2 = min(M2, max(M1, optCond2));
        %u3 = min(M2, max(M1, optCond3));
        
        
        % Update the control
        for dd=1:1:length(theta_convex)
            
            thetaa=theta_convex(dd);
            u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
            %u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
            %u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
            JCD_for_inside(dd)=cost_u1(Tx,X,u1_convex_for_inside(:,dd));
            
        end
        
        [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);
        
        u1_convex_for_inside(:,pos_min_for_inside);
        u1 = u1_convex_for_inside(:,pos_min_for_inside);
        %u2_convex_for_inside(:,pos_min_for_inside);
        %u2 = u2_convex_for_inside(:,pos_min_for_inside);
        %u3_convex_for_inside(:,pos_min_for_inside);
        %u3 = u3_convex_for_inside(:,pos_min_for_inside);
        
        %Use same interpolation  points
        %Tu_treatment(:,i)=Tx;
        
        % Stopping criteria on relative error
        temp1 = norm((oldu1 - u1))/norm((u1));
        %temp10 = norm((oldu2 - u2))/norm((u2));
        %temp11 = norm((oldu3 - u3))/norm((u3));
        temp2 = norm((oldx1 - x1))/norm(x1);
        temp3 = norm((oldx2 - x2))/norm(x2);
        temp4 = norm((oldx3 - x3))/norm(x3);
        temp5 = norm((oldx4 - x4))/norm(x4);
        
        temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
        temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
        temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
        temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));
        
        test = max(temp1, max(temp2, max( temp3, max( temp4, max(temp5,...
            max(temp6, max(temp7, max(temp8, temp9))))))))
        
        %test=temp1;
        
        XCD=X; TxCD=Tx;
        JCD_new(nn)=cost_u1(TxCD,XCD,u1);
        
        nn=nn+1;
        
    end
    
    MCF7_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,4);
    
    u1_low_dose_control_TypeI(:,i) = u1;
    %u2_low_dose_control(:,i) = u2;
    %u3_low_dose_control_TypeI_III(:,i) = u3;
    
    
    %%
    %% Type II
    %%
    
    initx=initx_treatment'; % IC for state
    initlambda=[0,0,0,0]';                  % FC for adjoint
    size(Tu_treatment(:,i))
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Solve OCP for Type I + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nn=1; test = 11;
    while(test > toll)
        
        %oldu1 = u1;
        oldu2 = u2;
        %oldu3 = u3;
        oldx1=x1;
        oldx2=x2;
        oldx3=x3;
        oldx4=x4;
        oldlambda1=lambda1;
        oldlambda2=lambda2;
        oldlambda3=lambda3;
        oldlambda4=lambda4;
        
        % solve the state eqn forward
        [Tx,X] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);
        
        % solve the adjoint eqn backward
        x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
        [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_lowE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
            flip(Tu_treatment(:,i)), initlambda, options);
        
        lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
        lambda4 = Lambda(:,4);
        
        % Interpolation
        lambda1 = interp1(Tlambda,lambda1,Tx);
        lambda2 = interp1(Tlambda,lambda2,Tx);
        lambda3 = interp1(Tlambda,lambda3,Tx);
        lambda4 = interp1(Tlambda,lambda4,Tx);
        
        
        % Optimality condition
        %optCond1 = (lambda4.*x4.*(1-Eq.m2*x4).*x1.*(1./(Eq.alpha3 + x1)))/Eq.weight1;
        optCond2 = (lambda1.*Eq.k2.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha2+x4)))/Eq.weight2;
        %optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;
        
        
        % Project the control
        %u1 = min(M2, max(M1, optCond1));
        u2 = min(M2, max(M1, optCond2));
        %u3 = min(M2, max(M1, optCond3));
        
        
        % Update the control
        for dd=1:1:length(theta_convex)
            
            thetaa=theta_convex(dd);
            %u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
            u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
            %u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
            JCD_for_inside(dd)=cost_u2(Tx,X,u2_convex_for_inside(:,dd));
            
        end
        
        [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);
        
        %u1_convex_for_inside(:,pos_min_for_inside);
        %u1 = u1_convex_for_inside(:,pos_min_for_inside);
        u2_convex_for_inside(:,pos_min_for_inside);
        u2 = u2_convex_for_inside(:,pos_min_for_inside);
        %u3_convex_for_inside(:,pos_min_for_inside);
        %u3 = u3_convex_for_inside(:,pos_min_for_inside);
        
        %Use same interpolation  points
        %Tu_treatment(:,i)=Tx;
        
        % Stopping criteria on relative error
        %temp1 = norm((oldu1 - u1))/norm((u1));
        temp10 = norm((oldu2 - u2))/norm((u2));
        %temp11 = norm((oldu3 - u3))/norm((u3));
        temp2 = norm((oldx1 - x1))/norm(x1);
        temp3 = norm((oldx2 - x2))/norm(x2);
        temp4 = norm((oldx3 - x3))/norm(x3);
        temp5 = norm((oldx4 - x4))/norm(x4);
        
        temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
        temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
        temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
        temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));
        
        test = max(temp10, max(temp2, max( temp3, max( temp4, max(temp5,...
            max(temp6, max(temp7, max(temp8, temp9))))))))
        
        %test=temp1;
        
        XCD=X; TxCD=Tx;
        JCD_new(nn)=cost_u2(TxCD,XCD,u2);
        
        nn=nn+1;
        
    end
    
    MCF7_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,4);
    
    %u1_low_dose_control_TypeI_III(:,i) = u1;
    u2_low_dose_control_TypeII(:,i) = u2;
    % u3_low_dose_control_TypeII_III(:,i) = u3;
    
    
    %%
    %% Type III
    %%
    
    initx=initx_treatment'; % IC for state
    initlambda=[0,0,0,0]';                  % FC for adjoint
    size(Tu_treatment(:,i))
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Solve OCP for Type I + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nn=1; test = 11;
    while(test > toll)
        
        %oldu1 = u1;
        %oldu2 = u2;
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
        [Tx,X] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);
        
        % solve the adjoint eqn backward
        x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
        [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_lowE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
            flip(Tu_treatment(:,i)), initlambda, options);
        
        lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
        lambda4 = Lambda(:,4);
        
        % Interpolation
        lambda1 = interp1(Tlambda,lambda1,Tx);
        lambda2 = interp1(Tlambda,lambda2,Tx);
        lambda3 = interp1(Tlambda,lambda3,Tx);
        lambda4 = interp1(Tlambda,lambda4,Tx);
        
        
        % Optimality condition
        %optCond1 = (lambda4.*x4.*(1-Eq.m2*x4).*x1.*(1./(Eq.alpha3 + x1)))/Eq.weight1;
        %optCond2 = (lambda1.*Eq.k3.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha3+x4)))/Eq.weight2;
        optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;
        
        
        % Project the control
        %u1 = min(M2, max(M1, optCond1));
        %u2 = min(M2, max(M1, optCond2));
        u3 = min(M2, max(M1, optCond3));
        
        
        % Update the control
        for dd=1:1:length(theta_convex)
            
            thetaa=theta_convex(dd);
            %u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
            u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
            u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
            JCD_for_inside(dd)=cost_u3(Tx,X,u3_convex_for_inside(:,dd));
            
        end
        
        [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);
        
        %u1_convex_for_inside(:,pos_min_for_inside);
        %u1 = u1_convex_for_inside(:,pos_min_for_inside);
        %u2_convex_for_inside(:,pos_min_for_inside);
        %u2 = u2_convex_for_inside(:,pos_min_for_inside);
        u3_convex_for_inside(:,pos_min_for_inside);
        u3 = u3_convex_for_inside(:,pos_min_for_inside);
        
        %Use same interpolation  points
        %Tu_treatment(:,i)=Tx;
        
        % Stopping criteria on relative error
        %temp1 = norm((oldu1 - u1))/norm((u1));
        %temp10 = norm((oldu2 - u2))/norm((u2));
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
            max(temp6, max(temp7, max(temp8, max(temp9,temp11))))))))
        
        %test=temp1;
        
        XCD=X; TxCD=Tx;
        JCD_new(nn)=cost_u3(TxCD,XCD,u3);
        
        nn=nn+1;
        
    end
    
    MCF7_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,4);
    
    %u1_low_dose_control_TypeI_III(:,i) = u1;
    %u2_low_dose_control_TypeII_III(:,i) = u2;
    u3_low_dose_control_TypeIII(:,i) = u3;
    
    %%
    %% Type I and III
    %%
    
    initx=initx_treatment'; % IC for state
    initlambda=[0,0,0,0]';                  % FC for adjoint
    size(Tu_treatment(:,i))
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Solve OCP for Type I + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nn=1; test = 11;
    while(test > toll)
        
        oldu1 = u1;
        %oldu2 = u2;
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
        [Tx,X] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);
        
        % solve the adjoint eqn backward
        x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
        [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_lowE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
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
        %optCond2 = (lambda1.*Eq.k2.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha2+x4)))/Eq.weight2;
        optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;

        
        
        % Project the control
        u1 = min(M2, max(M1, optCond1));
        %u2 = min(M2, max(M1, optCond2));
        u3 = min(M2, max(M1, optCond3));
        
        
        % Update the control
        for dd=1:1:length(theta_convex)
            
            thetaa=theta_convex(dd);
            u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
            %u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
            u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
            JCD_for_inside(dd)=cost_u1_u3(Tx,X,u1_convex_for_inside(:,dd),u3_convex_for_inside(:,dd));
            
        end
        
        [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);
        
        u1_convex_for_inside(:,pos_min_for_inside);
        u1 = u1_convex_for_inside(:,pos_min_for_inside);
        %u2_convex_for_inside(:,pos_min_for_inside);
        %u2 = u2_convex_for_inside(:,pos_min_for_inside);
        u3_convex_for_inside(:,pos_min_for_inside);
        u3 = u3_convex_for_inside(:,pos_min_for_inside);
        
        %Use same interpolation  points
        %Tu_treatment(:,i)=Tx;
        
        % Stopping criteria on relative error
        temp1 = norm((oldu1 - u1))/norm((u1));
        %temp10 = norm((oldu2 - u2))/norm((u2));
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
            max(temp6, max(temp7, max(temp8, max(temp9,temp11)))))))))
        
        %test=temp1;
        
        XCD=X; TxCD=Tx;
        JCD_new(nn)=cost_u1_u3(TxCD,XCD,u1,u3);
        
        nn=nn+1;
        
    end
    
    MCF7_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,4);
    
    u1_low_dose_control_TypeI_III(:,i) = u1;
    %u2_low_dose_control(:,i) = u2;
    u3_low_dose_control_TypeI_III(:,i) = u3;
    
    
    %%
    %% Type II and III
    %%
    
    initx=initx_treatment'; % IC for state
    initlambda=[0,0,0,0]';                  % FC for adjoint
    size(Tu_treatment(:,i))
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Solve OCP for Type I + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nn=1; test = 11;
    while(test > toll)
        
        %oldu1 = u1;
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
        [Tx,X] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx, options);
        
        % solve the adjoint eqn backward
        x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
        [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD_lowE2(t,lambda,u1,u2,u3,Tu_treatment(:,i),x1,x2,x3,x4,Tx,E2), ...
            flip(Tu_treatment(:,i)), initlambda, options);
        
        lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
        lambda4 = Lambda(:,4);
        
        % Interpolation
        lambda1 = interp1(Tlambda,lambda1,Tx);
        lambda2 = interp1(Tlambda,lambda2,Tx);
        lambda3 = interp1(Tlambda,lambda3,Tx);
        lambda4 = interp1(Tlambda,lambda4,Tx);
        
        
        % Optimality condition
        %optCond1 = (k1_hat*lambda4.*x4.*(1-Eq.m2*x4).*x1.*(1./(Eq.alpha3 + x1)))/Eq.weight1;
        optCond2 = (lambda1.*Eq.k2.*x1.*(1-Eq.m1*x1).*x4.*(1./(Eq.alpha2+x4)))/Eq.weight2;
        optCond3 = (lambda3-lambda2).*Eq.db.*E2.*x2/Eq.weight3;

        
        
        % Project the control
        %u1 = min(M2, max(M1, optCond1));
        u2 = min(M2, max(M1, optCond2));
        u3 = min(M2, max(M1, optCond3));
        
        
        % Update the control
        for dd=1:1:length(theta_convex)
            
            thetaa=theta_convex(dd);
            %u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
            u2_convex_for_inside(:,dd) = (thetaa)*u2 + (1-thetaa)*oldu2;
            u3_convex_for_inside(:,dd) = (thetaa)*u3 + (1-thetaa)*oldu3;
            JCD_for_inside(dd)=cost_u2_u3(Tx,X,u2_convex_for_inside(:,dd),u3_convex_for_inside(:,dd));
            
        end
        
        [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);
        
        %u1_convex_for_inside(:,pos_min_for_inside);
        %u1 = u1_convex_for_inside(:,pos_min_for_inside);
        u2_convex_for_inside(:,pos_min_for_inside);
        u2 = u2_convex_for_inside(:,pos_min_for_inside);
        u3_convex_for_inside(:,pos_min_for_inside);
        u3 = u3_convex_for_inside(:,pos_min_for_inside);
        
        %Use same interpolation  points
        %Tu_treatment(:,i)=Tx;
        
        % Stopping criteria on relative error
        %temp1 = norm((oldu1 - u1))/norm((u1));
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
            max(temp6, max(temp7, max(temp8, max(temp9,temp10)))))))))
        
        %test=temp1;
        
        XCD=X; TxCD=Tx;
        JCD_new(nn)=cost_u2_u3(TxCD,XCD,u2,u3);
        
        nn=nn+1;
        
    end
    
    MCF7_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,4);
    
    %u1_low_dose_control_TypeI_III(:,i) = u1;
    u2_low_dose_control_TypeII_III(:,i) = u2;
    u3_low_dose_control_TypeII_III(:,i) = u3;
    
    %%
    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    h = zeros(1,7);
    figure(930)
    subplot(3,2,1)
    h(1)=plot(Tu_treatment(:,i),u1_low_dose_control_TypeI_III(:,i),'k-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    %hold on
    %plot(Tu_treatment(:,i),u2_low_dose_control(:,i),'r--','LineWidth',4)
    hold on
    h(2)=plot(Tu_treatment(:,i),u3_low_dose_control_TypeI_III(:,i),'r-','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    grid on
    %legend('I+III','II+III','I','II','III','Location','NorthWest')
    legend('u_1','u_3','Location','NorthWest')
    xlabel('t')
    ylabel('u_1(t), u_3(t)')
    xlim([24,tf])
    ylim([0,1])
    sgtitle('Low dose of E2 - Control functions - u_1(t) (circle), u_2(t) (diamond) , u_3(t) (square)');
    %title(sprintf('ID  %d', i));
    %xline(t0_treatment_CD,'--k','HandleVisibility','off')
    
    
    %figure(530+i)
    subplot(3,2,2)
    h(3)=plot(Tu_treatment(:,i),u2_low_dose_control_TypeII_III(:,i),'b-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    h(4)=plot(Tu_treatment(:,i),u3_low_dose_control_TypeII_III(:,i),'r-','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    %hold on
    %plot(Tu_treatment(:,i),u3_low_dose_control_TypeI_III(:,i),'k-','LineWidth',4)
    grid on
    legend('u_2','u_3','Location','NorthWest')
    xlabel('t')
    ylabel('u_2(t), u_3(t)')
    xlim([24,tf])
    ylim([0,1])
    %title(sprintf('ID  %d', i));
    %xline(t0_treatment_CD,'--k','HandleVisibility','off')
    
    %figure(630+i)
    subplot(3,2,3)
    h(5)=plot(Tu_treatment(:,i),u1_low_dose_control_TypeI(:,i),'k-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    %hold on
    %plot(Tu_treatment(:,i),u3_low_dose_control_TypeI_III(:,i),'k-','LineWidth',4)
    grid on
    %legend('u_1','Location','NorthWest')
    xlabel('t')
    ylabel('u_1(t)')
    xlim([24,tf])
    ylim([0,1])
    %title(sprintf('ID  %d', i));
    %xline(t0_treatment_CD,'--k','HandleVisibility','off')
    
    %figure(730+i)
    subplot(3,2,4)
    h(6)=plot(Tu_treatment(:,i),u2_low_dose_control_TypeII(:,i),'b-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    %hold on
    %plot(Tu_treatment(:,i),u3_low_dose_control_TypeI_III(:,i),'k-','LineWidth',4)
    grid on
    %legend('u_2','Location','NorthWest')
    xlabel('t')
    ylabel('u_2(t)')
    xlim([24,tf])
    ylim([0,1])
    %title(sprintf('ID  %d', i));
    %xline(t0_treatment_CD,'--k','HandleVisibility','off')
    
    %figure(830+i)
    subplot(3,2,5)
    h(7)=plot(Tu_treatment(:,i),u3_low_dose_control_TypeIII(:,i),'r-','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    %hold on
    %plot(Tu_treatment(:,i),u3_low_dose_control_TypeI_III(:,i),'k-','LineWidth',4)
    grid on
    %legend('u_3','Location','NorthWest')
    xlabel('t')
    ylabel('u_3(t)')
    xlim([24,tf])
    ylim([0,1])
    %legend(h, labels, 'Position', [0.35 0.05 0.3 0.05]); 
    %legend(h,'Treatment - Type I+III', 'Treatment - Type II+III','Treatment - Type I','Treatment - Type II','Treatment - Type III','location', 'northwest',  'FontSize', 8)
    sgtitle('Control functions - low dose E2 (0.5mg) (with CAF)','fontsize',18)
    
    %%
    Percentage_change_TypeI(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))* (MCF7_with_CAF2_low_dose(:,i)-MCF7_with_CAF2_low_dose_control_TypeI(:,i))
    Percentage_change_TypeII(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))* (MCF7_with_CAF2_low_dose(:,i)-MCF7_with_CAF2_low_dose_control_TypeII(:,i))
    Percentage_change_TypeIII(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))*(MCF7_with_CAF2_low_dose(:,i)- MCF7_with_CAF2_low_dose_control_TypeIII(:,i))
    Percentage_change_TypeI_III(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))*(MCF7_with_CAF2_low_dose(:,i)- MCF7_with_CAF2_low_dose_control_TypeI_III(:,i))
    Percentage_change_TypeII_III(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))* (MCF7_with_CAF2_low_dose(:,i)-MCF7_with_CAF2_low_dose_control_TypeII_III(:,i))
end

    
    %%
    %%Treatment
    figure(412)
    subplot(2,2,1)
    %plot(t_data_low_dose, y_data_with_CAF2_low_dose_mice1, 'sb','LineWidth',2);
    %hold on
    plot(Tu,MCF7_with_CAF2_low_dose(:,1), '-b','LineWidth',5);
    hold on
    plot(Tu_treatment(:,1),MCF7_with_CAF2_low_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),MCF7_with_CAF2_low_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),MCF7_with_CAF2_low_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),MCF7_with_CAF2_low_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),MCF7_with_CAF2_low_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6,'LineWidth',2);
    hold on
    xline(Tu_treatment(1,1), '--r', 'Treatment starts', 'LabelOrientation', 'horizontal', ...
      'LabelVerticalAlignment', 'top', 'LineWidth', 2);
    xline(60, '-.b', 'E2 supply ends', 'LabelOrientation', 'horizontal', ...
      'LabelVerticalAlignment', 'middle', 'LineWidth', 1);
    legend('No treatment','Treatment - Type I','Treatment - Type II','Treatment - Type III','Treatment - Type I+III', 'Treatment - Type II+III','location', 'northwest',  'FontSize', 10)
    %title('T(t)','fontweight','normal','fontsize',18)
    sgtitle('Optimal treatment - Low dose E2 (0.5mg)','fontsize',18)
    xlabel('t','fontweight','normal','fontsize',18)
    ylabel('T(t)','fontweight','normal','fontsize',18)
    grid on
    xlim([0,160])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    
    
    
    %figure(212)
    subplot(2,2,2)
    plot(Tu,ER_with_CAF2_low_dose(:,1), '-b','LineWidth',5);
    hold on
    plot(Tu_treatment(:,1),ER_with_CAF2_low_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),ER_with_CAF2_low_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),ER_with_CAF2_low_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),ER_with_CAF2_low_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),ER_with_CAF2_low_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    xline(Tu_treatment(1,1), '--r', 'LineWidth', 2);
    xline(60, '-.b', 'LineWidth', 1);
    %legend('ID 1 - No treatment','Treatment - Type I', 'location', 'northwest')
    %sgtitle('low dose E2 (0.5mg) (with CAF2)','fontsize',18)
    %xlabel('t (days)','fontweight','normal','fontsize',18)
    xlabel('t','fontweight','normal','fontsize',18)
    ylabel('ER(t)','fontweight','normal','fontsize',18)
    grid on
    xlim([0,160])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %saveas(fig,'ER_with_CAF2','eps')
    %saveas(fig,'ER_with_CAF2','fig')
    
   
    
    %figure(213)
    subplot(2,2,3)
    plot(Tu,E2ER_with_CAF2_low_dose(:,1), '-b','LineWidth',5);
    hold on
    plot(Tu_treatment(:,1),E2ER_with_CAF2_low_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),E2ER_with_CAF2_low_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),E2ER_with_CAF2_low_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),E2ER_with_CAF2_low_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),E2ER_with_CAF2_low_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    xline(Tu_treatment(1,1), '--r', 'LineWidth', 2);
    xline(60, '-.b', 'LineWidth', 1);
    %legend('ID 1 - No treatment','Treatment - Type I', 'location', 'northwest')
    %sgtitle('low dose E2 (0.5mg) (with CAF2)','fontsize',18)
    %xlabel('t (days)','fontweight','normal','fontsize',18)
    %title('E2ER(t)','fontweight','normal','fontsize',18)
    xlabel('t','fontweight','normal','fontsize',18)
    ylabel('E2ER(t)','fontweight','normal','fontsize',18)
    grid on
    xlim([0,160])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %saveas(fig,'E2ER_with_CAF2','eps')
    %saveas(fig,'E2ER_with_CAF2','fig')
    
   
    
    %figure(214)
    subplot(2,2,4)
    plot(Tu,CAF2_with_CAF2_low_dose(:,1), '-b','LineWidth',5);
    hold on
    plot(Tu_treatment(:,1),CAF2_with_CAF2_low_dose_control_TypeI(:,1), '-o','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),CAF2_with_CAF2_low_dose_control_TypeII(:,1), '-s','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),CAF2_with_CAF2_low_dose_control_TypeIII(:,1), '-^','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),CAF2_with_CAF2_low_dose_control_TypeI_III(:,1), '-d','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    hold on
    plot(Tu_treatment(:,1),CAF2_with_CAF2_low_dose_control_TypeII_III(:,1), '-v','MarkerIndices', 1:50:length(Tu), 'MarkerSize', 6, 'LineWidth',2);
    xline(Tu_treatment(1,1), '--r', 'LineWidth', 2);
    xline(60, '-.b', 'LineWidth', 1);
    %legend('ID 1 - No treatment','Treatment - Type I', 'location', 'northwest')
    %sgtitle('low dose E2 (0.5mg) (with CAF2)','fontsize',18)
    %ylabel('t (days)','fontweight','normal','fontsize',18)
    xlabel('t','fontweight','normal','fontsize',18)
    ylabel('CAF(t)','fontweight','normal','fontsize',18)
    grid on
    xlim([0,160])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %saveas(fig,'CAF2_with_CAF2','eps')
    %saveas(fig,'CAF2_with_CAF2','fig')
    
    


end

%Sub-function

function J=cost(Tx,X,u1,u2,u3)
global Eq

J=(trapz(Tx,X(:,1)) + Eq.weight1*0.5*trapz(Tx,u1.^2) + Eq.weight2*0.5*trapz(Tx,u2.^2) + Eq.weight3*0.5*trapz(Tx,u3.^2));
end

function J=cost_u1_u3(Tx,X,u1,u3)
global Eq

J=(trapz(Tx,X(:,1)) + Eq.weight1*0.5*trapz(Tx,u1.^2)+ Eq.weight3*0.5*trapz(Tx,u3.^2));
end

function J=cost_u2_u3(Tx,X,u2,u3)
global Eq

J=(trapz(Tx,X(:,1)) + Eq.weight2*0.5*trapz(Tx,u2.^2)+ Eq.weight3*0.5*trapz(Tx,u3.^2));
end

function J=cost_u1(Tx,X,u1)
global Eq

J=(trapz(Tx,X(:,1)) + Eq.weight1*0.5*trapz(Tx,u1.^2));
end

function J=cost_u2(Tx,X,u2)
global Eq

J=(trapz(Tx,X(:,1)) + Eq.weight2*0.5*trapz(Tx,u2.^2));
end

function J=cost_u3(Tx,X,u3)
global Eq

J=(trapz(Tx,X(:,1)) + Eq.weight3*0.5*trapz(Tx,u3.^2));
end


