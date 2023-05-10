[nOutliers,hasOutliers,stats1, Nthresh] =  for_based_thresh(intensities, 35,gain,adFactor,countOffset,roNoise,qStar,U,stopIntensity);
[nOutliers,hasOutliers,stats2, Nthresh] =  for_based_thresh(intensities, 38,gain,adFactor,countOffset,roNoise,qStar,U,stopIntensity);
[nOutliers,hasOutliers,stats3, Nthresh] =  for_based_thresh(intensities, 40,gain,adFactor,countOffset,roNoise,qStar,U,stopIntensity);

figure,plot(stats1.xval,stats1.yval)
hold on
plot(stats2.xval,stats2.yval)
plot(stats3.xval,stats3.yval)


figure,plot(stats1.allVals,stats1.FOR)
hold on
plot(stats2.allVals,stats2.FOR)

plot(stats3.allVals,stats3.FOR)


%% -
pyenv('Version','python.exe')
pe = pyenv;
mod = py.importlib.import_module('sklearn');
pyrun('from sklearn.linear_model import TheilSenRegressor')

yval = py.numpy.array(stats1.yval);
Xval = py.numpy.array(stats1.xval);
reg = pyrun("reg = TheilSenRegressor(random_state=0).fit(X.reshape(-1, 1), y)","reg",X = Xval,y = yval);
coefV(j,:) = [double(reg.coef_) double(reg.intercept_)];

figure,plot(stats1.xval,stats1.yval)
