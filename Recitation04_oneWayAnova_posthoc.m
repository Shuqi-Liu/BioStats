%% recitation 4
close all
clc
clear all

a = normrnd(49.842105263157897, 9.337895992391575, [41,1]);
b = normrnd(37.022727272727273, 10.269595305455981, [41,1]);
c = normrnd(39.147058823529413, 9.773595699921593, [41,1]);
data = [a,b,c];

%test normality
pnorm=ones(1,3)*-1;
for i = 1:3
    pnorm(i) = lillietest(data(:,i)); %0 dont reject
end

%test variance using lavene test
phomovar = vartestn(data,'TestType','LeveneAbsolute');
%test variance of 2 groups
%hvarbc = vartest2(B,C)

%anova
[p,tbl,stats] = anova1(data); %each column is a group

%% Hand calculation
meanOverall = mean(data,'all');
k = 3;
N = numel(data);
nj = length(data);
msbtw = sum(nj* (mean(data) - meanOverall).^2) / (k-1)
mserror = sum((nj - 1)*std(data).^2)/(N - k)
SSerror = mserror * (N-k)

%% bonferroni, matlab
bonferroniRes = multcompare(stats,'CType','bonferroni');
%2 sample t by hand and check against corrected alpha -- by hand doesnt
%work
[~,P_bf_12,CI12,STATS12] = ttest2(a,b);
tstats12 = abs(mean(a) - mean(b)) / sqrt(mserror * 2 / nj);
p12_hand = 2*(1 - tcdf(tstats12, 2*nj - 2))
tcrit12 = tinv(1-0.05/3,80);
[H13,P_bf_13,CI13,STATS13] = ttest2(a,c);
tinv(1-0.05/3,64) %2.16 vs 2.22 reject
[H23,P_bf_23,CI23,STATS23] = ttest2(b,c);
tinv(1-0.05/3,70) %2.017 vs 1.6001 fail to reject
adjustedAlpha = 0.05/3

bft12 = abs(mean(A) - mean(B)) / sqrt(MSe / (na) + MSe / (nb));
p_bf_12 = 1 - tcdf(bft12, na+nb-2);
bft13 = abs(mean(A) - mean(C)) / sqrt(MSe * (1/ na + 1 / nc));
p_bf_13 = 1 - tcdf(bft13, na+nc-2);
bft23 = abs(mean(C) - mean(B)) / sqrt(MSe / (nc) + MSe / (nb));
p_bf_23 = 1 - tcdf(bft23, nc+nb-2);

%% tukey HSD, matlab, also generates a graph --- this works
hsdRes = multcompare(stats,'CType','hsd');

totalN = numel(data);
alpha = 0.05;
numGroup = 3;
qcritical = qdist(alpha,numGroup,totalN-numGroup);
MSe = tbl(3,4);
MSe = MSe{1,1};
MWbtw = tbl(2,4);
MWbtw = MWbtw{1,1};
ns = [length(a), length(b), length(c)];

%could compare q stats against q critical, or compare p value
%pvalue comparable to matlab hsd results output
q12 = abs(mean(a) - mean(b)) / sqrt(MSe / (2*ns(1)) + MSe / (2*ns(2)));
pval12 = 1 - cdfTukey(q12, totalN-groupsNum , groupsNum); 
% use the toal anova df
q13 = abs(mean(a) - mean(c)) / sqrt(MSe / (2*ns(1)) + MSe / (2*ns(3)));
pval13 = 1 - cdfTukey(q13, totalN-groupsNum , groupsNum);
q23 = abs(mean(b) - mean(c)) / sqrt(MSe / (2*ns(2)) + MSe / (2*ns(3)));
pval23 = 1 - cdfTukey(q23, totalN-groupsNum , groupsNum);
%cdfTukey is approximation

%%  Dunnett's test - only works for balanced data
%2nd arg index list of the experimental groups, 3rd arg the index of the control group
pDunnett = dunnett(stats,[2,3],1);
tDunCrit = 2.38;
diffDun = tDunCrit * sqrt(MSe / ns(1));
var = sum(sum((data - mean(data)).^2)) / (totalN - numGroup - 1);
minDiff = tDunCrit * sqrt(2* MSe / ns(1)); %reject 12,13
%16.428908031808959, 11.389537538859791, -5.039370492949168
diff = abs([mean(a) - mean(b), mean(a) - mean(c)]);
resultsDunnett = diff >= minDiff;
%seems to have no way to find p value by hand

%alternatively
tDunnett = diff / sqrt(2* MSe / ns(1));
resultsDunnett2 = tDunnett >= tDunCrit;

%% equivalently anova with unbalanced groups
A = [78,26,73,32,57,48,39,118,42,35,56,44,47,46,35,59,69,28,33,42,31,51,24,40,69,38,45,48,60,42,63,78,39,42,81,47,44,45];
B = [30,30,32,38,23,61,29,43,29,32,35,32,54,55,48,36,33,38,26,27,37,49,48,43,34,59,31,34,30,28,28,64,24,36,47,42,27,31,37,30,33,45,34,27];
C = [44,38,27,39,30,65,28,40,62,39,42,38,50,33,43,45,24,34,62,43,41,34,52,40,54,28,30,45,40,42,43,41,39,44];
data = [A, B, C];
groups = repelem([{'A'}, {'B'}, {'C'}], [length(A) length(B) length(C)]);

h = vartestn(data',groups','TestType','LeveneAbsolute');

n = max([numel(A), numel(B), numel(C)]);
A(end+1:n)=nan;
B(end+1:n)=nan;
C(end+1:n)=nan;
data2 = [A;B;C]';
[p2,tbl2,stats2] = anova1(data2);


