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


