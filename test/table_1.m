
% chipPars to latex table
offset_fun = @(f,g,b1,b2) f.*(b1-b2)./(2*g-1);
% function for sigma from mean/variance relations
sigma_fun = @(b1,delta,f) sqrt(b1-1/12+delta./f);

% chipPars.slope is 1/f
% chipPars.offset b1
% b2 is  gainPars = cellfun(@(x) x(1),chipPars.gaincoeffs{2});
% g is  (gainPars./(chipPars.slope(perm)))/2;

deltaQ = zeros(4,3);
sigmaQ = zeros(4,3);
gainQ = zeros(4,3);
for ii=2:4
% now we get delta: 50
    gainPars = cellfun(@(x) x(2),chipPars.gaincoeffs{ii});
    Cgain = cellfun(@(x) x(1),chipPars.gaincoeffs{ii});
    perm = randperm(length(gainPars));
    g = (gainPars./(chipPars.slope(perm(1:length(gainPars)))))/2;

    offsetVals = offset_fun(1./chipPars.slope((1:length(gainPars))),g,chipPars.offset(perm(1:length(gainPars))),Cgain);

    deltaQ(ii,:) = quantile(offsetVals,3);
    sigmaVals = sigma_fun(chipPars.offset(perm(1:length(gainPars))),offsetVals,1./chipPars.slope(perm(1:length(gainPars))));
    sigmaQ(ii,:) = quantile(sigmaVals,3);
    gainQ(ii,:) = quantile(g,3);

end
Q = quantile(1./chipPars.slope,3);

% display('\begin{table}[ht]')
fprintf('\\begin{table}[ht]\n\\centering\n\\begin{tabular}{|l|l|l|l|l|}\n')
fprintf('\\hline\nGain & f (interquartile range) & $\\rm{g}$ & $\\rm{r}$ & $\\Delta$  \\\\ \n')
fprintf('\\hline\n')

fprintf('0&%4.2f ( %4.2f - %4.2f) & -&-&- \\\\ \n',Q([2 1 3]))
fprintf('\\hline\n')

% for i=2:4
fprintf('50& -& %4.2f ( %4.2f - %4.2f) & %4.2f ( %4.2f - %4.2f) & %4.2f ( %4.2f - %4.2f) \\\\ \n',gainQ(2,[2 1 3]),sigmaQ(2,[2 1 3]),deltaQ(2,[2 1 3]));
fprintf('\\hline\n')
fprintf('100& -& %4.2f ( %4.2f - %4.2f) & %4.2f ( %4.2f - %4.2f) & %4.2f ( %4.2f - %4.2f) \\\\ \n',gainQ(3,[2 1 3]),sigmaQ(3,[2 1 3]),deltaQ(3,[2 1 3]));
fprintf('\\hline\n')
fprintf('300& -& %4.2f ( %4.2f - %4.2f) & %4.2f ( %4.2f - %4.2f) & %4.2f ( %4.2f - %4.2f) \\\\ \n',gainQ(4,[2 1 3]),sigmaQ(4,[2 1 3]),deltaQ(4,[2 1 3]));
fprintf('\\hline\n')
% end
      
fprintf('\n\\end{tabular}\n\\caption{\\label{calibrationtable} Chip parameter estimates from the data shown in Fig.~\\ref{calibfig} calculated using Eqs. (\\ref{new_mean}) - (\\ref{nogainvar})}\n')
fprintf('\\end{table}\n')
%     G

% Condition & n & p \\
% \hline
% A & 5 & 0.1 \\
% \hline
% B & 10 & 0.01 \\
% \hline
% \end{tabular}
% \caption{\label{tab:example}Legend (350 words max). Example legend text.}
% \end{table}