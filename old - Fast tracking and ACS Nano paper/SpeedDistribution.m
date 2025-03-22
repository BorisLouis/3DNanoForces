
step = 15;
alldata = [allSpeedR300;allSpeedR500;allSpeedR1000; allSpeedR2000];


edges = logspace(log10(min(alldata(:))),log10(max(alldata(:))));
edges = min(alldata(:)):step:max(alldata(:));

centers = edges(2:end);
centers = edges(1:end-1)+step;
[N300W,~] = histcounts(allSpeedR300,edges);

[N500W,~] = histcounts(allSpeedR500,edges);

[N1000W,~] = histcounts(allSpeedR1000,edges);

[N2000W,~] = histcounts(allSpeedR2000,edges);

Nf = {cat(1,centers,N300W)',cat(1,centers,N500W)',cat(1,centers,N1000W)',...
    cat(1,centers,N2000W)'};

%plot distribution
distributionPlot(Nf,'histOpt',0, 'showMM',0,'globalNorm',3);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  

ylim([-50 450]);
%set(gca,'YScale','log')