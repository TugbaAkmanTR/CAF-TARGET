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

function main_OCP_Fig4F_1_2_3_constant_11Dec()

close all
clear
clc
format long

set(0, 'defaultaxesfontsize',16)

global Eq

%% Fig.4E
%Initial conditions
E2 = 0.1;
ss = 0.99;

init_Tumor = 0.5;
init_ER = 1.300000e-04 ;
init_E2ER = 3.700000e-11;
init_CAF = 0.3*init_Tumor;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uncontrolled case - No treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1; %one mouse
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

% Time discretization for intervention
%index_treatment(i) = min(find(MCF7_with_CAF2_med_dose(:,i) >= Eq.threshold(i) ))
index_treatment(i) = min(find(Tu == 45))
t0_treatment = Tu(index_treatment(i))

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the state eqn forward
u1=0*ones(size(Tu_treatment(:,i)))+ss;
u2=0*ones(size(Tu_treatment(:,i)));
u3=0*ones(size(Tu_treatment(:,i)));

[~,XCD] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

MCF7_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeI(:,i) = XCD(:,4);

[Ju1,P1,P2]=cost_u1(Tu_treatment(:,i),MCF7_with_CAF2_med_dose_control_TypeI,u1,u1);
fprintf('Type I \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)));
u2=0*ones(size(Tu_treatment(:,i)))+ss;
u3=0*ones(size(Tu_treatment(:,i)));

% solve the state eqn forward
[~,XCD] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

MCF7_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeII(:,i) = XCD(:,4);

[Ju2,P1,P2]=cost_u2(Tu_treatment(:,i),MCF7_with_CAF2_med_dose_control_TypeII,u2,u2);
fprintf('Type II \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)));
u2=0*ones(size(Tu_treatment(:,i)));
u3=0*ones(size(Tu_treatment(:,i)))+ss;

% solve the state eqn forward
[~,XCD] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

MCF7_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeIII(:,i) = XCD(:,4);

[Ju3,P1,P2]=cost_u3(Tu_treatment(:,i),MCF7_with_CAF2_med_dose_control_TypeIII,u3,u3);
fprintf('Type III \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type I + III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)))+ss;
u2=0*ones(size(Tu_treatment(:,i)));
u3=0*ones(size(Tu_treatment(:,i)))+ss;

% solve the state eqn forward
[~,XCD] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);


MCF7_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeI_III(:,i) = XCD(:,4);

[Ju1u3,P1,P2]=cost_u1_u3(Tu_treatment(:,i),MCF7_with_CAF2_med_dose_control_TypeI_III,u1,u3);
fprintf('Type I+III \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju1u3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type II + III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu_treatment(:,i)));
u2=0*ones(size(Tu_treatment(:,i)))+ss;
u3=0*ones(size(Tu_treatment(:,i)))+ss;

% solve the state eqn forward
[~,XCD] = ode15s(@(t,x)stateEq_medE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

MCF7_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,1);
ER_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,2);
E2ER_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,3);
CAF2_with_CAF2_med_dose_control_TypeII_III(:,i) = XCD(:,4);

[Ju2u3,P1,P2]=cost_u1_u3(Tu_treatment(:,i),MCF7_with_CAF2_med_dose_control_TypeII_III,u2,u3);
fprintf('Type II+III \n');
fprintf('Tumor size = %.3f, Control = %.3f, J = %.3f \n', P1,P2,Ju2u3);

%%
Percentage_change_TypeI(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))* (MCF7_with_CAF2_med_dose(:,i)-MCF7_with_CAF2_med_dose_control_TypeI(:,i))
Percentage_change_TypeII(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))* (MCF7_with_CAF2_med_dose(:,i)-MCF7_with_CAF2_med_dose_control_TypeII(:,i))
Percentage_change_TypeIII(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))*(MCF7_with_CAF2_med_dose(:,i)- MCF7_with_CAF2_med_dose_control_TypeIII(:,i))
Percentage_change_TypeI_III(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))*(MCF7_with_CAF2_med_dose(:,i)- MCF7_with_CAF2_med_dose_control_TypeI_III(:,i))
Percentage_change_TypeII_III(i) = 100*(1/MCF7_with_CAF2_med_dose(:,i))* (MCF7_with_CAF2_med_dose(:,i)-MCF7_with_CAF2_med_dose_control_TypeII_III(:,i))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%Treatment
figure(400)
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
sgtitle('Constant treatment - Medium dose E2 (0.1mg)','fontsize',16)
xlabel('t','fontweight','normal','fontsize',16)
ylabel('T(t)','fontweight','normal','fontsize',16)
grid on
xlim([0,160])
ylim([0,20])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

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
sgtitle('Constant treatment - Med dose E2 (0.1mg)','fontsize',18)
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
