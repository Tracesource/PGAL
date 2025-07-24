clear;
clc;
warning off;
addpath(genpath('./'));

%% dataset
ds = {'ForestTypes'};
dsPath = './datasets/';
resPath = './res-lmd/';

for dsi =1:1:length(ds)
    dataName = ds{dsi}; disp(dataName);
    load(strcat(dsPath,dataName));
    k = length(unique(Y));
    n = length(Y);
    
    lambda = [1e-6 1e-2 10^2 10^6];
    anchor = [k 2*k 5*k];
   
    txtpathmean = strcat(resPath,strcat(dataName,'_mean.txt'));
    dlmwrite(txtpathmean, strcat('Dataset:',cellstr(dataName), '  Date:',datestr(now)),'-append','delimiter','','newline','pc');
    
    %%
    allresultmax = [];
    allresultmean = [];
    for ichor = 1:length(anchor)
        for id = 1:length(lambda)
            tic;
            [A,Z,U,iter,obj,obj_all] = PGAL(X,k,anchor(ichor),lambda(id));
            [resmean,std] = myNMIACCwithmean(U,Y,k);
            timer(ichor,id)  = toc;
            fprintf('Anchor:%d \t Lambda:%d\t Res:%12.6f %12.6f %12.6f %12.6f \tTime:%12.6f \n',[anchor(ichor) lambda(id) resmean(1) resmean(2) resmean(3) resmean(4) timer(ichor,id)]);
            dlmwrite(txtpathmean,  [anchor(ichor) lambda(id) resmean std timer(ichor,id)],'-append','delimiter','\t','newline','pc');
            allresultmean = [allresultmean;resmean std timer(ichor,id)];
        end
    end
    [c,d] = max(allresultmean(:,1));
    maxresultmean = allresultmean(d,:);
    dlmwrite('./totalResults/AllDatasetResultMean.txt',char(dataName),'-append','delimiter','\t','newline','pc');
    dlmwrite('./totalResults/AllDatasetResultMean.txt',maxresultmean,'-append','delimiter','\t','newline','pc');
end


