%% recitation 5: two way anova
[p, tbl, stats] = anova2(data,5) %5 entries per group(cell)
%format text
for i = 1:15
    if i <=5
        fprintf('A & %d & %d\\\\', data(i,1),data(i,2))
    elseif i <= 10
        fprintf('B & %d & %d\\\\', data(i,1),data(i,2))
    else
        fprintf('C & %d & %d\\\\', data(i,1),data(i,2))
    end
    fprintf('\n\\hline\n')
end


means = mean(A);
fprintf('A & %.2f & %.2f\\\\\n\\hline\n', means(1),means(2))
means = mean(B);
fprintf('B & %.2f & %.2f\\\\\n\\hline\n', means(1),means(2))
means = mean(C);
fprintf('C & %.2f & %.2f\\\\\n\\hline\n', means(1),means(2))

groupNames = ['A','B','C'];
%format mean
for i = 1:3
    fprintf('%s& %.2f & %.2f\\\\\n\\hline\n', groupNames(i), mean(data(i*5 - 4:i*5,1)),mean(data(i*5 - 4:i*5,2)))
end

%format std
for i = 1:3
    fprintf('%s& %.2f & %.2f\\\\\n\\hline\n', groupNames(i), std(data(i*5 - 4:i*5,1)),std(data(i*5 - 4:i*5,2)))
end


%test normality
hnorm=ones(3,2)*-1;
for i = 1:3
    for j = 1:2
        hnorm(i,j) = lillietest(data(i*5 - 4:i*5,j)); %0 dont reject, it's normal
    end
end

vars=ones(3,2)*-1;
vars2 = ones(3,2)*-1;
for i = 1:3
    for j = 1:2
        vars(i,j) = std(data(i*5 - 4:i*5,j)); %0 dont reject, it's normal
        temp = data(i*5 - 4:i*5,j);
        meantemp = mean(temp,'all');
        vars2(i,j) = sqrt(sum((temp - meantemp).^2)/(length(temp)-1));
    end
end


%test variance using lavene test
phomovar = vartestn([A,B,C],'TestType','LeveneAbsolute');%large p does not reject null of equal variance

%visualize variance using bar plot
x = 1:6;
dataInColumn = [A,B,C];
means = mean(dataInColumn);
error = std(dataInColumn);
bar(x,means)                
hold on
er = errorbar(x,means,error,error,'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlabel('Group Combinations')
ylabel('Group Average with SD')
title('Visualization of Variance of each Cell')
xticklabels({'A,Male','A,Female','B,Male','B,Female','C,Male','C,Female'})
set(gca,'FontSize',20)
hold off


%find group mean for treatment type
mean(data(1:5,:),'all')
mean(data(6:10,:),'all')
mean(data(11:15,:),'all')

%calculate MSAB
meanTable = [mean(A);mean(B);mean(C)];
meanTable = [meanTable; mean(meanTable)];
meanTable = [meanTable, mean(meanTable,2)]
overallMean = meanTable(4,3)
mserror = 4* (sum(std(A).^2) + sum(std(B).^2) + sum(std(C).^2)) / (30-2*3);
msTreatment = 0;
for i = 1:3
   d = data(i*5 - 4:i*5,:);
   msTreatment = msTreatment + (mean(d,'all') - overallMean)^2;
end
mstreatment = 5*2*msTreatment / (3-1);
Ftreatment = mstreatment/mserror;
ptreament = 1 - fcdf(Ftreatment, 2, 24);
fcrit_treatment = finv(0.95, 2,24);

msSex = 5*3*sum((mean(data) - overallMean).^2) / 1;
FSex = msSex / mserror;
pSex = 1 - fcdf(FSex, 1, 24);
fcrit_sex = finv(0.95, 1,24);

MSinter = 0;
for i = 1:3
    for j = 1:2
        MSinter = MSinter + 5 * (meanTable(i,j) - meanTable(i,end) - ...
            meanTable(end,j) + meanTable(end,end))^2;
    end
end
MSinter = MSinter / (2 * 1)
Finter = MSinter / mserror
pinter = 1 - fcdf(Finter, 2, 24)
fcrit_inter = finv(0.95, 2,24)

%% post-hoc
sex_post = multcompare(stats) %default tukey
treatment_post = multcompare(stats,'Estimate','row')

%% fake inclass data and test
% 2x3 expanded to (2x3) x 3
data = nan(a*b, b);
for i = 1:a
    for j = 1:b
        data(i*r - r + 1:i*r,j) = [means(i,j), means(i,j)-SDs(i,j),means(i,j)+SDs(i,j)];
    end
end

%%
r = 5;
[p, tbl, stats] = anova2(data,r) %5 entries per group(cell)
%column ABCD, row: pounds 
col_post = multcompare(stats) %default tukey, 
%col_post = post hoc of the variable at the column, should be of size: 
% num of combinations of columes in the data or # levels in the column variable x 6
% choose 2 out of the # of columns
row_post = multcompare(stats,'Estimate','row')
% size: row count / r x 6

%% Class example data
% SDs= repmat(20,2,3);
% means = [50 112 155; 94 200 43];
% % means = [50 112 155; 87 154 195];
% r = 3;
% a=2; %#of rows
% b = 3; % # of cols
% N = r*a*b;

%% Recitation example data
r = 5;
a=4; %#of rows
b = 4; % # of cols
N = r*a*b;

SDs = nan(a,b);
means = nan(a,b);
hnorm=nan(a,b);
reshaped_for_var_test = nan(r, a*b);
for row = 1:a
    for col = 1:b 
       means(row,col) = mean(data(row*r - r + 1:row*r,col));
       SDs(row,col) = std(data(row*r - r + 1:row*r,col));
       hnorm(row,col) = lillietest(data(row*r - r + 1:row*r,col));
       reshaped_for_var_test(:,(row-1)*b +col) = data(row*r - r + 1:row*r,col);
    end
end

%test variance using lavene test
phomovar = vartestn(reshaped_for_var_test,'TestType','LeveneAbsolute');%large p does not reject null of equal variance

%visualize variance using bar plot
x = 1:a*b;
means_reshaped = mean(reshaped_for_var_test);
error = std(reshaped_for_var_test);
bar(x,means_reshaped)                
hold on
er = errorbar(x,means_reshaped,error,error,'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlabel('Group Combinations')
ylabel('Group Average with SD')
title('Visualization of Variance of each Cell')
% xticklabels({'A,Male','A,Female','B,Male','B,Female','C,Male','C,Female'})
set(gca,'FontSize',20)
hold off

%% stats
MSerror = sum((r-1) * SDs.^2,'all') / (N - a*b)

means_row = mean(means,2);
means_col = mean(means,1);
grand_mean = mean(means,'all');

% disp('row is gender, col is age')
MS_row = r*b*sum((means_row - grand_mean).^2) / (a-1);
F_row = MS_row / MSerror;
prow = 1-fcdf(F_row, a-1, N-a*b);

MS_col = r*a*sum((means_col - grand_mean).^2) / (b-1);
F_col = MS_col / MSerror;
pcol = 1-fcdf(F_col, b-1, N-a*b);
f_crit_col = finv(0.95, b-1, N-a*b);

MS_interaction = r*sum((means - means_row - means_col + grand_mean) .^2,'all') / ((a-1) * (b-1));
F_interaction = MS_interaction / MSerror
p_interaction = 1 - fcdf(F_interaction, (a-1)*(b-1), N-a*b)

%% determine main effect: eg, age significance because it has 3 levels
% qcrit = qdist(0.05, b, N-b);
% hsd = qcrit * sqrt(MSerror / (a*r))
% age_diff = nan(1,3);
% age_diff(1) = abs(means_col(2) - means_col(1));
% age_diff(2) = abs(means_col(3) - means_col(1));
% age_diff(3) = abs(means_col(3) - means_col(2))
% 
% %
% qstats = age_diff / sqrt(MSerror / (a*r));
% for q = qstats
%     p = 1 - cdfTukey(q, N-b, b)
% end

%% with interaction
k = 11%from other calculation, k(k-1) /2 >= 48, k = 11
qcrit_interaction = qdist(0.05, k, N-a*b)
hsd_interaction = qcrit_interaction * sqrt(MSerror / r)

combs = nchoosek(1:a, 2); %generate all possible orders of variable A
combs_n = nchoosek(a,2);
diffs = nan(nchoosek(a,2)*b + nchoosek(b,2)*a,1);
idx = 1;
for col = 1:b
    for c = 1:combs_n
        diffs(idx) = (means(combs(c,1),col) - means(combs(c,2),col));
        idx = idx + 1;
    end
end
sigpairs = diffs > hsd_interaction

% all_diff = nan(a,3); %cross group of row, then 12, 13, 23
% % 
% all_diff(1,:) = abs(means(2,:) - means(1,:));
% diff = abs(bsxfun(@minus,means(1,:),means(1,:)')); 
% all_diff(2,1:2) = diff(1,2:3);
% all_diff(2,3) = diff(2,3);
% 
% diff = abs(bsxfun(@minus,means(2,:),means(2,:)')); 
% all_diff(3,1:2) = diff(1,2:3);
% all_diff(3,3) = diff(2,3);