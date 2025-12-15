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
t_data_high_dose = [17,24,32,39,45,52,62,67,73,80,87,95,102,108,115,119];
%t_data_low_dose = [17,24,32,39,45,52,62,67,73,80,87,95,102,108,115,119,124];
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
%
% y_data_with_CAF2_low_dose = [y_data_with_CAF2_low_dose_mice1; y_data_with_CAF2_low_dose_mice2; y_data_with_CAF2_low_dose_mice3;...
%     y_data_with_CAF2_low_dose_mice4; y_data_with_CAF2_low_dose_mice5];
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

for i=1:1 % 5 mice

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Type I
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve the state eqn forward
    u1=0*ones(size(Tu_treatment(:,i)))+0.5;
    u2=0*ones(size(Tu_treatment(:,i)));
    u3=0*ones(size(Tu_treatment(:,i)));

    [Tx,XCD] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

    MCF7_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeI(:,i) = XCD(:,4);


    %%
    %% Type II
    %
    %Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
    u1=0*ones(size(Tu_treatment(:,i)));
    u2=0*ones(size(Tu_treatment(:,i)))+0.5;
    u3=0*ones(size(Tu_treatment(:,i)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Type II
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve the state eqn forward
    [Tx,XCD] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

    MCF7_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeII(:,i) = XCD(:,4);

    %%
    %% Type III
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Type III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
    u1=0*ones(size(Tu_treatment(:,i)));
    u2=0*ones(size(Tu_treatment(:,i)));
    u3=0*ones(size(Tu_treatment(:,i)))+0.5;

    % solve the state eqn forward
    [Tx,XCD] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

    MCF7_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeIII(:,i) = XCD(:,4);

    %%
    %% Type I and III
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Solve OCP for Type I + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
    u1=0*ones(size(Tu_treatment(:,i)))+0.5;
    u2=0*ones(size(Tu_treatment(:,i)));
    u3=0*ones(size(Tu_treatment(:,i)))+0.5;

    % solve the state eqn forward
    [Tx,XCD] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);


    MCF7_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeI_III(:,i) = XCD(:,4);

    %%
    %% Type II and III
    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Type II + III
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
    u1=0*ones(size(Tu_treatment(:,i)));
    u2=0*ones(size(Tu_treatment(:,i)))+0.5;
    u3=0*ones(size(Tu_treatment(:,i)))+0.5;

    % solve the state eqn forward
    [Tx,XCD] = ode15s(@(t,x)stateEq_lowE2(t,x,u1,u2,u3,Tu,E2), Tu_treatment(:,i), initx_treatment, options);

    MCF7_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,1);
    ER_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,2);
    E2ER_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,3);
    CAF2_with_CAF2_low_dose_control_TypeII_III(:,i) = XCD(:,4);


    %%
    %%
    Percentage_change_TypeI(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))* (MCF7_with_CAF2_low_dose(:,i)-MCF7_with_CAF2_low_dose_control_TypeI(:,i))
    Percentage_change_TypeII(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))* (MCF7_with_CAF2_low_dose(:,i)-MCF7_with_CAF2_low_dose_control_TypeII(:,i))
    Percentage_change_TypeIII(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))*(MCF7_with_CAF2_low_dose(:,i)- MCF7_with_CAF2_low_dose_control_TypeIII(:,i))
    Percentage_change_TypeI_III(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))*(MCF7_with_CAF2_low_dose(:,i)- MCF7_with_CAF2_low_dose_control_TypeI_III(:,i))
    Percentage_change_TypeII_III(i) = 100*(1/MCF7_with_CAF2_low_dose(:,i))* (MCF7_with_CAF2_low_dose(:,i)-MCF7_with_CAF2_low_dose_control_TypeII_III(:,i))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
sgtitle('Constant treatment - Low dose E2 (0.025mg)','fontsize',18)
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


