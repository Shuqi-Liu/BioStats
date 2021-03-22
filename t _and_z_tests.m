%% function lookup
% lookfor mean 
%look for any internal function that contains/relates to key words mean 
%(in function name or even in main doc)
whos %list the current variables in the workspace
pwd 
cd ..
cd 'Code' %step into directoy needs single quote
%% normal distribution
p_left = normcdf(2,0,1) %args: z-value, mu, sigma; returns: area to the left
z_val = norminv(p_left,0,1) %should give 2, args: area, mu and sigma

%% studnet's t-test
t_pdf_val = tpdf(-2.5706,5) %provides the pdf evaluated at 0 with df=50
p_left = tcdf(-2.5706,5) %provides the cdf evaluated from -inf to 0 with df=50, area to the left of 0
t_val = tinv(0.025,5) %corresponding t-value from -inf to t that gives the given area
%expected results: 0.0303, 0.0250, -2.5706

%t-stats and CI calculation needs to be done manulally (not calling matlab
%inheritent function)
%quiz question: 
n = [1623,1191];
sds = [17.5,20.1];
means = [128.2,126.5];
pooleds = sqrt(((n-1) * (sds.^2)') / (sum(n)-2));
t_crit = abs(tinv(0.025,sum(n)-2))
ci = [0,0];
ci(1) = means(1) - means(2) - t_crit * pooleds * sqrt(sum(1./n))
ci(2) = means(1) - means(2) + t_crit * pooleds * sqrt(sum(1./n))

ci_ver1 = [0,0];
ci_ver1(1) = means(1) - means(2) - t_crit * pooleds * sqrt(1/n(1) + 1/n(2))
ci_ver1(2) = means(1) - means(2) + t_crit * pooleds * sqrt(1/n(1) + 1/n(2))

%% sample size + power analysis
%p0: test-type, mean and var, p1 mean, power, 
%t:one sample t-test
n = sampsizepwr('t',[100 10],110,.8)
p=sampsizepwr('z',[0 8],5,[],40,'alpha',0.05,'tail','right')%power
%% recitation 03
se = std(difference)/sqrt(13) %no built func for this
[h,p,ci,stats] = ttest(difference) %this will give CI with default alpha or provided alpha, but needs the source data
[h,p,ci,stats] = ttest2(difference, test) %two sample t-test, this will give CI with default alpha or provided alpha, but needs the source data

%% sample size
nn = 1:100;
pwrout = sampsizepwr('z',[0 std(difference)],mean(difference),[],nn,'alpha',0.1,'tail','both');
p13 = sampsizepwr('z',[0 std(difference)],mean(difference),[],13,'alpha',0.1,'tail','both');
figure('Units','normalized','Position',[0.1,0.1,0.6,0.6])
plot(nn,pwrout,'b-',13,p13,'ro','LineWidth',3,'MarkerSize',20)
title('Power versus Sample Size')
xlabel('Sample Size')
ylabel('Power')
set(gca,'FontSize',20)

ef = mean(difference)-1:0.1:mean(difference)+1;
pwrout2 = sampsizepwr('z',[0 std(difference)],ef,[],50,'alpha',0.1,'tail','both');
figure('Units','normalized','Position',[0.1,0.1,0.6,0.6])
plot(ef,pwrout2,'b-',mean(difference),p13,'ro','LineWidth',3,'MarkerSize',20)
title('Power versus Alternatie Mean')
xlabel('\mu 1')
ylabel('Power')
set(gca,'FontSize',20)


%% examples of 2 normal distribution with non-normal sum
% https://stats.stackexchange.com/questions/30159/is-it-possible-to-have-a-pair-of-gaussian-random-variables-for-which-the-joint-d/30160#30160
X = normrnd(0,1,[100,1]);
B = binornd(1,0.5,[100,1]);
Y = X.*(2.*B-1);