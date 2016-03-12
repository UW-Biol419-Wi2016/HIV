clear all; close all;
%% list of important variables

% pr = pr nucleotide aligned table,           920 by 297
% pr_protein = pr protein aligned table,      920 by 99
% rt = rt nucleotide aligned table trimmed to 920 by 1482
% rt_protein = rt protein nucleotide table    920 by 494

% nucleotide_pr = 920 by 297 cell,  each bin = pr nucleotide
% protein_pr = 920 by 99 cell,      each bin = pr protein
% nucleotide_rt = 920 by 1482 cell, each bin = rt nucleotide
% protein_rt = 920 by 494 cell,     each bin = rt protein

% avg_pr =         1 by 297 char
% avg_pr_protein = 1 by 99 char  
% avg_rt =         1 by 1482 char
% avg_rt_protein = 1 by 494 char

% avg_pr_cell =         1 by 297 cell
% avg_pr_protein_cell = 1 by 99 cell 
% avg_rt_cell =         1 by 1482 cell, 1 by 903, rest empty
% avg_rt_protein_cell = 1 by 494 cell, 1 by 301, rest empty

% pr_nuc_dev = 920 by 297 matrix of 0 and 1, 1 = deviation from avg
% pr_prot_dev = 920 by 99 matrix of deviation = 1
% rt_nuc_dev = 920 by 1482 matrix
% rt_prot_dev = 920 by 494 matrix

% pr_pca_mat = [VL_mat, CD_mat, pr_nuc_dev] % 920 by 299 matrix
% rt_pca_mat = [VL_mat, CD_mat, rt_nuc_dev] % 920 by 905 matrix
% pr_prot_pca_mat = [VL_mat, CD_mat, pr_prot_dev] % 920 by 101 matrix
% rt_prot_pca_mat = [VL_mat, CD_mat, rt_prot_dev] % 920 by 496 matrix

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

%% average nucleotides from 'pr final'

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

avg_pr=[];
for i=1:297 % 297 nucleotides long
    temp = [proteinbin{:,i}];
    avg_pr = [avg_pr, char(mode(double(temp)))];
end; % 297 nucleotides long, same as averageseq_pr1

% cell conversion of matrix
avg_pr_cell = cell(1,297);
for i=1:297 % 297 nucleotides long
    avg_pr_cell{1,i} = avg_pr(i);
end;

% nucleotide_pr = string cell where each bin = each pr nucleotdie
% 920 by 297

% avg_pr = string cell where each column = mode of each column
% 1by 297
%% mode protein sequence from 'pr protein final'
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

avg_pr_protein=[];
for i=1:99 % 99 proteins long
    temp = [proteinbin{:,i}];
    avg_pr_protein = [avg_pr_protein, char(mode(double(temp)))];
end; % 99 proteins long

% cell conversion of matrix
avg_pr__protein_cell = cell(1,99);
for i=1:99 % 99 proteins long
    avg_pr_protein_cell{1,i} = avg_pr_protein(i);
end;

%% average rt nucleotide from rt final

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
avg_rt=[];

for i=1:1482 % 1482 nucleotides long at most
    temp = [proteinbin{:,i}];
    avg_rt = [avg_rt, char(mode(double(temp)))];
end;

% cell conversion of matrix
avg_rt_cell = cell(1,1482);
for i=1:1482 % 1482 nucleotides long
    avg_rt_cell{1,i} = avg_rt(i);
end;

%% mode rt protein sequence from 'rt protein final'
proteinbin = cell(920,494);
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

avg_rt_protein=[];
for i=1:494 % 1482/3 = 494 proteins long
    temp = [proteinbin{:,i}];
    avg_rt_protein = [avg_rt_protein, char(mode(double(temp)))];
end; % 494 proteins long

% cell conversion of matrix
avg_rt_protein_cell = cell(1,494);
for i=1:494 % 494 proteins long
    avg_rt_protein_cell{1,i} = avg_rt_protein(i);
end;

%% a corresponding matrix of deviation = 1, same = 0

pr_nuc_dev = [];
for i = 1:920
    for j = 1:297
        if strcmp(nucleotide_pr{i,j},avg_pr_cell{1,j}) == 1 % if same
            pr_nuc_dev(i,j) = 0; % if same, no deviation, then 0
        elseif nucleotide_pr{i,j} == char(0)
            pr_nuc_dev(i,j) = 0;
        else
            pr_nuc_dev(i,j) = 1; % if different, yes deviation, 1
        end
    end
end

rt_nuc_dev = [];
for i = 1:920
    for j = 1:1482
        if strcmp(nucleotide_rt{i,j},avg_rt_cell{1,j}) == 1 % if same
            rt_nuc_dev(i,j) = 0; % if same, no deviation, then 0
        elseif nucleotide_rt{i,j} == char(0)
            rt_nuc_dev(i,j) = 0;
        else
            rt_nuc_dev(i,j) = 1; % if different, yes deviation, 1
        end
    end
end

rt_nuc_dev = rt_nuc_dev(:,1:903); % trim the rest of 903 columns

%%
traindata = readtable('training_data.csv', 'Delimiter', ',');
testdata = readtable('test_data.csv', 'Delimiter',',');

train_response = traindata(1:920,2);
train_response = table2cell(train_response);
train_response = cell2mat(train_response);
VL = traindata(1:920,5); % VL cell count column of train csv 
VL_mat = table2cell(VL);
VL_mat = cell2mat(VL_mat);
CD = traindata(1:920,6); % CD4 cell count column of train csv
CD_mat = table2cell(CD);
CD_mat = cell2mat(CD_mat);

pr_pca_mat = horzcat(VL_mat, CD_mat, pr_nuc_dev); % 920 by 299 matrix
rt_pca_mat = horzcat(VL_mat, CD_mat, rt_nuc_dev); % 920 by 905 matrix
%% pr nucleotide prediction
[coeff, score, latent] = pca(pr_pca_mat);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

pr_nuc_classify = fitcdiscr(score(train,1:10),train_response(train),'DiscrimType','pseudoQuadratic');

rt_predict = predict(pr_nuc_classify,score(test,1:10));

cv_quad = mean(rt_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end

pr_nuc_cv_acc_quad = mean(cv_acc_quad) % 0.7812

%% rt nucleotide prediction
[coeff, score, latent] = pca(rt_pca_mat);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

rt_nuc_classify = fitcdiscr(score(train,1:10),train_response(train),'DiscrimType','pseudoQuadratic');

rt_predict = predict(rt_nuc_classify,score(test,1:10));

cv_quad = mean(rt_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end

rt_nuc_cv_acc_quad = mean(cv_acc_quad) % 0.7840

%% 

pr_prot_dev = [];
for i = 1:920
    for j = 1:99
        if strcmp(protein_pr{i,j},avg_pr_protein_cell{1,j}) == 1 % if same
            pr_prot_dev(i,j) = 0; % if same, no deviation, then 0
        elseif protein_pr{i,j} == char(0)
            pr_prot_dev(i,j) = 0;
        else
            pr_prot_dev(i,j) = 1; % if different, yes deviation, 1
        end
    end
end

rt_prot_dev = [];
for i = 1:920
    for j = 1:494
        if strcmp(protein_rt{i,j},avg_rt_protein_cell{1,j}) == 1 % if same
            rt_prot_dev(i,j) = 0; % if same, no deviation, then 0
        elseif protein_rt{i,j} == char(0)
            rt_prot_dev(i,j) = 0;
        else
            rt_prot_dev(i,j) = 1; % if different, yes deviation, 1
        end
    end
end

rt_prot_dev = rt_prot_dev(:,1:301); % trim the rest of 903 columns

pr_prot_pca_mat = horzcat(VL_mat, CD_mat, pr_prot_dev); % 920 by 299 matrix
rt_prot_pca_mat = horzcat(VL_mat, CD_mat, rt_prot_dev); % 920 by 905 matrix

%% pr protein prediction
[coeff, score, latent] = pca(pr_prot_pca_mat);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

pr_prot_classify = fitcdiscr(score(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

rt_prot_predict = predict(pr_prot_classify,score(test,:));

cv_quad = mean(rt_prot_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end

pr_prot_cv_acc_quad = mean(cv_acc_quad) % 0.7836

%% rt protein prediction
[coeff, score, latent] = pca(rt_prot_pca_mat);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

rt_prot_classify = fitcdiscr(score(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

rt_prot_predict = predict(rt_prot_classify,score(test,:));

cv_quad = mean(rt_prot_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end

rt_prot_cv_acc_quad = mean(cv_acc_quad) % 0.7855 w/ 10 scores, 0.31 w/ all

%%  rt nucleotide deviation matrix only pca 
[coeff, score, latent] = pca(rt_nuc_dev);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

rt_nuc_classify2 = fitcdiscr(score(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

rt_nuc_predict2 = predict(rt_nuc_classify2,score(test,:));

cv_quad = horzcat(rt_nuc_predict2, train_response(test));
cv_acc_quad(i)= mean(cv_quad(:,1) == cv_quad(:,2));
end

rt_nuc_cv_acc_quad = mean(cv_acc_quad) % 0.2628 when all scores used


%%  rt protein deviation matrix only pca 
[coeff, score, latent] = pca(rt_prot_dev);

for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

rt_prot_classify2 = fitcdiscr(score(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

rt_prot_predict2 = predict(rt_prot_classify2,score(test,:));

cv_quad = horzcat(rt_prot_predict2, train_response(test));
cv_acc_quad(i)= mean(cv_quad(:,1) == cv_quad(:,2));
end

rt_prot_cv_acc_quad = mean(cv_acc_quad) % 0.3149
