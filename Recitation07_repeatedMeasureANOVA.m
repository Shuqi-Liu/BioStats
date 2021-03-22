%% maltab internal: doesn't seem to work
% t = table([1:10]', temp(:,1),temp(:,2),temp(:,3),temp(:,4),'VariableNames',{'subject','d1','d2','d3','d4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
rm = fitrm(t,'d4-d1~subject','WithinDesign',Meas)
ranovatbl = ranova(rm)
data = [t.d1, t.d2, t.d3, t.d4];

writematrix(data,'/Users/mac/Desktop/BioStats/Recitations/Recitation07/WedCytokinData.csv') 
save('/Users/mac/Desktop/BioStats/Recitations/Recitation07/WedCytokinData.mat','data')
%% hand calculation
[r,k] =size(data);
alpha = .05;
mean_subj = mean(data,2);
mean_groups = mean(data);
mean_grand = mean(data,'all');
sd_groups = std(data);
MSeffect = r * sum((mean_groups - mean_grand).^2) / (k-1);
SSwithin = (r-1) * sum(sd_groups.^2);
SSrandom = k* sum((mean_subj - mean_grand).^2);%diff in participants
SSerror = SSwithin - SSrandom;
MSerror = SSerror/ ((r-1) * (k-1));
F = MSeffect/ MSerror;
pvalue = 1 - fcdf(F,k-1, (r-1) * (k-1));
F_crit = finv(1 - alpha, k-1, (r-1) * (k-1));
%reject null that one of group is different
%% post hoc - bonferroni sequential
% bonferroni for repeated measure uses paired t-test, where as regular
% anova uses 2 sample t test
p_bonf = nan(1,k-1);
adjusted_alpha = alpha/(k-1); %numOfComparison = 3 for 4 groups in this case
for i = 1:k-1 
    [h,p] = ttest(data(:,i),data(:,i+1)); %default 2 tail
    p_bonf(i) = p;
    if (p <= adjusted_alpha) %reject null, this group is diff
        fprintf('\nsig diff pair: %d %d\n',i,i+1);
    end
end


%% post hoc - Tukey
qcrit = qdist(alpha, k, r*k-k);
hsd = qcrit * sqrt(MSerror / r);

combs = nchoosek(1:k, 2); %generate all possible orders of variable A
combs_n = nchoosek(k,2);
diffs = nan(1, combs_n);
for c = 1:combs_n
    diffs(c) = mean_groups(combs(c,1)) - mean_groups(combs(c,2));
end
diffs = abs(diffs);
sigpairs = diffs > hsd;

qstats = diffs / (sqrt(MSerror / r));
qstats = qstats/2;
pval_tukey = [];
for q = qstats
    pval_tukey(end+1) = 1-cdfTukey(q, (r-1)*(k-1), k);
end

%% posthoc -dunnet
%2nd arg index list of the experimental groups, 3rd arg the index of the control group
% pDunnett = dunnett(stats,[2,3],1);
tDunCrit = 3.1; %k=4, df = 6
tDunCirt = 2.49;
diffDun = tDunCrit * sqrt(MSerror / r);
diff = abs(mean_groups - mean_groups(1));
diff = diff(2:end);
resultsDunnett = diff >= diffDun
%seems to have no way to find p value by hand

%alternatively
tDunnett = diff / sqrt(MSerror / r);
% resultsDunnett2 = tDunnett >= tDunCrit;


%% if ignore repeat
[p,tbl,stats] = anova1(data);
%if ignoring the repeat, will get p > alpha, failed to reject null and fail
%to conclude there is a difference between the groups.
%hence, using the repeated measure increases power
