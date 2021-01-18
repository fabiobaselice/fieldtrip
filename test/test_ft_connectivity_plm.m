function test_ft_connectivity_plm

% Purpose: simulate linearly mixed data, and compute the PLM value

% create some instantaneously mixed data

% define some variables locally
nTrials  = 10;
nSamples = 5000;
fsample  = 250;

% mixing matrix
mixing   = [0.7 0.3 0;
              0 0.4 0.6];

data       = [];
data.trial = cell(1,nTrials);
data.time  = cell(1,nTrials);
for k = 1:nTrials
  dat = randn(3, nSamples);
  dat(2,:) = ft_preproc_bandpassfilter(dat(2,:), 1000, [7 15]);
  dat = 0.2.*(dat-repmat(mean(dat,2),[1 nSamples]))./repmat(std(dat,[],2),[1 nSamples]);
  data.trial{k} = mixing * dat;
  data.time{k}  = (0:nSamples-1)./fsample;
end
data.label = {'chan1' 'chan2' }';

figure;plot(dat'+repmat([0 1 2],[nSamples 1]));
title('original ''sources''');

figure;plot((mixing*dat)'+repmat([0 1],[nSamples 1])); 
axis([0 1000 -1 2]);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('mixed ''sources''');

% compute PLM connectivity
cfg = [];
cfg.method = 'plm';
cfg.fsample= fsample;
p = ft_connectivityanalysis(cfg, data);
disp(mean(p.plm,3))
