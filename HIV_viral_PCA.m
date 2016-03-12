%%
clear all; close all;
%%
traindata = readtable('training_data.csv', 'Delimiter', ',');
testdata = readtable('test_data.csv', 'Delimiter',',');

train_response = traindata(1:920,2);
train_response = table2cell(train_response);
train_response = cell2mat(train_response);
VL = traindata(1:920,5); % VL cell count column of train csv 
VL_cell = table2cell(VL);
VL_mat = cell2mat(VL_cell); % conversion of table into matrix
CD = traindata(1:920,6); % CD4 cell count column of train csv
CD_cell = table2cell(CD);
CD_mat = cell2mat(CD_cell); % conversion of table into matrix

viral_mat = horzcat(VL_mat, CD_mat); % 920 by 2 matrix

%% PCA
[coeff, score, latent] = pca(viral_mat);

%% LDA
for i=1:100
test_frac = 0.2; % fraction of dataset to use for testing
permuted = randperm(920); 
test = permuted(1:floor(920*test_frac)); 
train = permuted(ceil((920*test_frac)):end);

viral_class = fitcdiscr(score(train,:),train_response(train),'DiscrimType','pseudoQuadratic');

viral_predict = predict(viral_class,score(test,:));

cv_quad = mean(viral_predict == train_response(test));
cv_acc_quad(i)= cv_quad;
end

viral_cv_acc_quad = mean(cv_acc_quad) % 0.7933
