%% 419 project Mark Choi, Davina Kang, Mary Chang
clear all; close all;
%% set working directory and read in table
traindata = readtable('training_data.csv', 'Delimiter', ',');
size(traindata)
summary(traindata)

% 1st column = patient ID = number
% 2nd column = response to treatment = 0 or 1 = number
% 3rd = protease sequence = string of C,G,T,A, etc
% 4th = rt sequence = string
% 5th = VL-t0 cell concentration = unit unclear
% 6th = CD4-t0 cell concentration = unit unclear

% make a subtable of patients who did not get better (Resp = 0)
train0 = traindata(traindata.Resp==0, :);
% make a subtable of patients who got better (Resp = 1)
train1 = traindata(traindata.Resp==1, :);



%% Visualize with histograms

% make subtable column of train0 and train1
vl0 = train0(:,5);
vl1 = train1(:,5);
cd0 = train0(:,6);
cd1 = train1(:,6);

% convert tables into homogenous arrays
vl0 =table2array(vl0);
vl1 =table2array(vl1);
cd0 =table2array(cd0);
cd1 =table2array(cd1);

figure;
histogram(vl0, 'FaceColor', 'r');
hold on;
histogram(vl1, 'FaceColor', 'b');
legend('Dead','Alive');
title('VL-t0 Counts of HIV Patients after Treatment');
xlabel('VL-t0 Cell Concentration (unit unclear)');
ylabel('# of Patients');

% need to edit and unify bins
figure;
histogram(cd0, 'FaceColor', 'r');
hold on;
histogram(cd1, 'FaceColor', 'b');
legend('Dead','Alive');
title('CD4-t0 Counts of HIV Patients after Treatment');
xlabel('CD4-t0 Cell Concentration (unit unclear)');
ylabel('# of Patients');


