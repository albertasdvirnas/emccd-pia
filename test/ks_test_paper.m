        
lamGuess = 30:1:44;
lamGuess = 10:30;

tv = zeros(1,length(lamGuess));
minks = zeros(1,length(lamGuess));

for i=1:length(lamGuess)
    i
    [tv(i),posMax,bins,minks(i)] = run_ktest(sortI,lamGuess(i),gain, adFactor, countOffset, roNoise);
end

figure,plot(lamGuess,tv)
xlabel('lambda bg value')
ylabel('thresh')

figure,plot(lamGuess,minks)
xlabel('lambda bg value')
ylabel('KS min')

%% 

tic
lamGuess = 38;
[thresh,bins,ks,minks,threshPos] = run_ktest(sortI,lamGuess,gain, adFactor, countOffset, roNoise);
toc

figure,plot(bins,ks)

%%
sortI = sortI(sortI<60);
randGrid = randi(length(sortI),1,1000);
sample = sortI(randGrid);
tic
lamGuess = 38;
[thresh,bins,ks,minks,threshPos] = run_ktest(sample,lamGuess,gain, adFactor, countOffset, roNoise);
toc

figure,plot(bins,ks)
