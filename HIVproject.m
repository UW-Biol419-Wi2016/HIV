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

%% loading in test data
testdata = readtable('test_data.csv', 'Delimiter', ',');
size(testdata)
summary(testdata)

% 1st column = patient ID = number
% 2nd column = response to treatment = H(?) = string
% 3rd = protease sequence = string of C,G,T,A, etc
% 4th = rt sequence = string
% 5th = VL-t0 cell concentration = unit unclear
% 6th = CD4-t0 cell concentration = unit unclear

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
xlabel('VL-t0 Cell Concentration (per ml of blood)');
ylabel('# of Patients');

% need to edit and unify bins
figure;
histogram(cd0, 'FaceColor', 'r');
hold on;
histogram(cd1, 'FaceColor', 'b');
legend('Dead','Alive');
title('CD4-t0 Counts of HIV Patients after Treatment');
xlabel('CD4-t0 Cell Concentration (per ml of blood)');
ylabel('# of Patients');

%% bivariate (3D) histogram

figure;
hist3([cd0,vl0],[11,15],'FaceAlpha',0.65);
xlabel('CD4-t0 cell count');
ylabel('VL-t0 cell count');
title('3D Histogram of Patients who died');
figure;
hist3([cd1,vl1],[11,15]);
xlabel('CD4-t0 cell count');
ylabel('VL-t0 cell count');
title('3D Histogram of Patients who survived'); 
%% 

<<<<<<< HEAD
%% read in FASTA files using fastaread()

pr0 = fastaread('pr0 aligned');
pr1 = fastaread('pr1 aligned');
rt0 = fastaread('rt0 aligned');
rt1 = fastaread('rt1 aligned');

% convert structure into table
pr0 = struct2table(pr0);
pr1 = struct2table(pr1);
rt0 = struct2table(rt0);
rt1 = struct2table(rt1);
%% average nucleotides from pr1 aligned

a = pr1(1,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only

proteinbin = cell(187,length(c));

for i = 1:187, % i equals number of rows in pr1
a = pr1(i,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only
    for j = 1:length(c),
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = 0; % char(0) is a blank space
        end;
    end;
end; % now protein bin has a nucleotide in each cell
nucleotide_pr1 = proteinbin;
averageseq_pr1=[];

for i=1:297 % nucleotides 297 long
    temp = [proteinbin{:,i}];
    averageseq_pr1 = [averageseq_pr1, char(mode(double(temp)))];
end; % 297 nucleotides long, same as averageseq_pr1

% modesequence is an "average" sequence of nucleotides that appear most
% frequently at each site
%% average nucleotide sequence from pr0 aligned

a = pr0(1,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only

proteinbin = cell(733,length(c));

for i = 1:733, % i equals number of rows in pr0
a = pr0(i,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only
    for j = 1:length(c),
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = 0; % 
        end;
    end;
end; % now protein bin has a nucleotide in each cell
nucleotide_pr0 = proteinbin;
averageseq_pr0=[];

for i=1:length(c) % 297 nucleotides long
    temp = [proteinbin{:,i}];
    averageseq_pr0 = [averageseq_pr0, char(mode(double(temp)))]; 
end; %297 nucleotides long, same as averageseq_pr1

%% average nucleotide from rt0 aligned

a = rt0(1,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only

proteinbin = cell(794,length(c));

for i = 1:794, % i equals number of rows in rt0
a = rt0(i,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only
    for j = 1:length(c),
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = 0; % blank space
        end;
    end;
end; % now protein bin has a nucleotide in each cell
nucleotide_rt0 = proteinbin;
averageseq_rt0=[];

for i=1:length(c) % 1482 nucleotides long at most
    temp = [proteinbin{:,i}];
    averageseq_rt0 = [averageseq_rt0, char(mode(double(temp)))];
end;
%% average nucleotide from rt1 aligned

a = rt1(1,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only

proteinbin = cell(187,length(c));

for i = 1:187, % i equals number of rows in rt1
a = rt1(i,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only
    for j = 1:length(c),
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = 0; % char(0) is a blank space
        end;
    end;
end; % now protein bin has a nucleotide in each cell
nucleotide_rt1 = proteinbin;
averageseq_rt1=[];

for i=1:length(c) %1476
    temp = [proteinbin{:,i}];
    averageseq_rt1 = [averageseq_rt1, char(mode(double(temp)))];
end;
%%
mean(averageseq_pr0==averageseq_pr1)


% mean(averageseq_rt0==averageseq_rt1)

%% mutation table/cell

for i = 1:length(averageseq_rt0)
    avg_rt0(1,i)=cellstr(averageseq_rt0(i));
end;

iscellstr(averageseq_rt0(1,1)); % 0  no 
iscellstr(nucleotide_rt0(1,1)); % 1  yes
iscellstr(avg_rt0(1,1));        % 1  yes

x = size(nucleotide_rt0);
for i = 1:x(1)
    for j = 1:x(2)
        if strcmp(nucleotide_rt0(i,j), avg_rt0(1,j)) == 0
            mutation_rt0{i,j} = nucleotide_rt0{i,j};
        end;
    end;
end;

%% 
x = size(nucleotide_rt0);
for i = 1:x(1)
    for j = 1:x(2)
        if mutation_rt0{i,j} == 0
            mat_rt0(i,j) = 0;
        else
            mat_rt0(i,j) = 1;
        end;
=======
pr0 = struct2table(fastaread('pr0 protein aligned'));
pr1 = struct2table(fastaread('pr1 protein aligned'));
rt0 = struct2table(fastaread('rt0 protein aligned'));
rt1 = struct2table(fastaread('rt1 protein aligned'));

for i = 1:101
    for j = 1:733
        avgpr0 = mode(pr0{j,i})
>>>>>>> origin/master
    end;
end;

%%
for i = 1:101
    % run LDA using classify
    classify(test(X(:,i)),train training groups)
    strncmp()
    
% look for the position of the amino acid that is most associated with death

train0 = traindata(traindata.Resp==0, :);

train1 = traindata(traindata.Resp==1, :);
