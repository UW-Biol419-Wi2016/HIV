%% 
clear all; close all;
% pr = pr nucleotide aligned table,           920 by 297
% pr_protein = pr protein aligned table,      920 by 99
% rt = rt nucleotide aligned table trimmed to 920 by 1482
% rt_protein = rt protein nucleotide table    920 by 494

% nucleotide_pr = 920 by 297 cell,  each bin = pr nucleotide
% protein_pr = 920 by 99 cell,      each bin = pr protein
% nucleotide_rt = 920 by 1482 cell, each bin = rt nucleotide
% protein_rt = 920 by 494 cell,     each bin = rt protein

%% read in FASTA files using fastaread()
pr = fastaread('pr final'); % 920 by 297
pr_protein = fastaread('pr protein final'); % 920 by 101
rt = fastaread('rt final'); % 1000 rows
rt_protein = fastaread('rt protein final'); % 1000 by 

% convert structure into table
pr = struct2table(pr);
pr_protein = struct2table(pr_protein);
rt = struct2table(rt);
rt_protein = struct2table(rt_protein);
rt = rt(1:920,:); % trim the last 80 rows
rt_protein = rt_protein(1:920,:); % trim the last 80 rows

%% nucleotide_pr

proteinbin = cell(920,297); 

for i = 1:920, % i equals number of rows in pr final
a = pr(i,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only, length(c) = 297
    for j = 1:297,
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = char(0); % char(0) is a blank space
        end;
    end;
end; % now protein bin has a nucleotide in each cell
nucleotide_pr = proteinbin;

%% protein_pr
proteinbin=cell(920,99);
for i = 1:920, % i equals number of rows in pr final
a = pr_protein(i,2); % sequence cell with header
b = a{1,1}; % 'PQITL...'   sequence cell only
c = b{1}; % PQITL   protein sequence string only
    for j = 1:99, % 99 proteins long
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = char(0); % char(0) is a blank space
        end;
    end;
end; 
protein_pr = proteinbin;
%% nucleotide_pr
proteinbin = cell(920,1482);
for i = 1:920, % i equals number of rows in rt final
a = rt(i,2); % sequence cell with header
b = a{1,1}; % sequence cell only
c = b{1}; % sequence string only
    for j = 1:1482,
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = char(0); % blank space
        end;
    end;
end; % now protein bin has a nucleotide in each cell
nucleotide_rt = proteinbin;

%% protein_rt
proteinbin=cell(920,494);
for i = 1:920, % i equals number of rows in pr final
a = rt_protein(i,2); % sequence cell with header
b = a{1,1}; % 'PQITL...'   sequence cell only
c = b{1}; % PQITL   protein sequence string only
    for j = 1:494, % 494 proteins long
        proteinbin{i,j} = c(j);
        if proteinbin{i,j} == '-', % replace - into a space
            proteinbin{i,j} = char(0); % char(0) is a blank space
        end;
    end;
end; 
protein_rt = proteinbin;

%%
traindata = readtable('training_data.csv', 'Delimiter', ',');
testdata = readtable('test_data.csv', 'Delimiter',',');

train_response = traindata(1:920,2);
train_response = table2cell(train_response);
train_response = cell2mat(train_response);
VL = traindata(1:920,5); % VL cell count column of train csv 
VL_cell = table2cell(VL);
% VL_mat = cell2mat(VL_mat);
CD = traindata(1:920,6); % CD4 cell count column of train csv
CD_cell = table2cell(CD);
% CD_mat = cell2mat(CD_mat);

pr_nuc_cell = horzcat(VL_cell, CD_cell, nucleotide_pr); % 920 by 299 matrix
rt_nuc_cell = horzcat(VL_cell, CD_cell, nucleotide_rt); % 920 by 905 matrix

pr_prot_cell = horzcat(VL_cell, CD_cell, protein_pr); % 920 by 101 matrix
rt_prot_cell = horzcat(VL_cell, CD_cell, protein_rt); % 920 by 303 matrix

%%
[coeff, score, latent] = pca(pr_nuc_cell);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

pr_nuc_class = fitcdiscr(pr_nuc_cell(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

rt_predict = predict(pr_nuc_class,score(test,1:299));

cv_quad = mean(rt_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end
cv_1_acc =sum((cv_quad(:,1)==cv_quad(:,2)&cv_quad(:,2)==1))/sum(cv_quad(:,2)==1)
pr_nuc_cv_acc_quad = mean(cv_acc_quad) % 0.7884

%% 
[coeff, score, latent] = pca(rt_nuc_cell);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

rt_nuc_class = fitcdiscr(rt_nuc_cell(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

rt_predict = predict(rt_nuc_class,score(test,1:905));

cv_quad = mean(rt_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end
cv_1_acc =sum((cv_quad(:,1)==cv_quad(:,2)&cv_quad(:,2)==1))/sum(cv_quad(:,2)==1)
rt_nuc_cv_acc_quad = mean(cv_acc_quad) % 0.2074



