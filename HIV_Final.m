clear all;

%%
test_data = readtable('test_data.csv','Delimiter',',');
% size(test_data)  % 692 by 6
test_ID = table2array(test_data(:,1));
test_pr_seq = table2array(test_data(:,3));
test_rt_seq = table2array(test_data(:,4));

%% create test_pr_seq and test_rt_seq fasta files, unaligned
% Do not run this section twice because it will add 1-692 rows on top of 
% the pre-existing file, respectively

% warnState = warning %Save the current warning state
% warning('off','Bioinfo:fastawrite:AppendToFile'); 
% for i = 1:692 
% fastawrite('test_pr_seq',num2str(i),test_pr_seq(i));
% fastawrite('test_rt_seq',num2str(i),test_rt_seq(i));
% end
% warning(warnState) %Reset warning state to previous settings

% Open the files in seaview  
% align -> align options -> muscle,  align -> align all
% props -> view as proteins,  file -> save prot alignment
%% 
% in the rt aligned file, row 437 is weird and it messes up alignment of all
% other sequences

%% histograms of all test patients 

% make subtable column of train0 and train1
test_VL = test_data(:,5);
test_CD = test_data(:,6);

% convert tables into homogenous arrays
test_VL = table2array(test_VL);
test_CD = table2array(test_CD);

figure;
histogram(test_VL, 'FaceColor', 'r');
title('VL-t0 Counts of test HIV Patients after Treatment');
xlabel('VL-t0 Cell Concentration (per mL of blood)');
ylabel('# of Patients');

figure;
histogram(test_CD, 'FaceColor', 'r');
title('CD4-t0 Counts of test HIV Patients after Treatment');
xlabel('CD4-t0 Cell Concentration (per mL of blood)');
ylabel('# of Patients');

%% bivariate (3D) histogram

figure;
hist3([test_CD,test_VL],[11,15],'FaceAlpha',0.65);
xlabel('CD4-t0 Cell Conc');
ylabel('VL-t0 Cell Conc');
title('Bivariate Histogram of Test HIV Patients');

%% 

test_pr = fastaread('test_pr_seq_aligned');
test_rt = fastaread('test_rt_seq_aligned');
test_pr_prot = fastaread('test_pr_prot_aligned');
test_rt_prot = fastaread('test_rt_prot_aligned');

% convert structure into table
test_pr = struct2table(test_pr);
test_rt = struct2table(test_rt);
test_pr_prot = struct2table(test_pr_prot);
test_rt_prot = struct2table(test_rt_prot);




