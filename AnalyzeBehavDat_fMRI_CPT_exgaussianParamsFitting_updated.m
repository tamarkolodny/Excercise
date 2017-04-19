function [params] = AnalyzeBehavDat_fMRI_CPT_exgaussianParamsFitting_updated(outputFolder,summaryFileName)
    
    %outputFolder = 'FittingOutput';
    %summaryFileName = 'FittingResults_ScannerData.txt';
    
    addpath('C:\Users\Tamaroy\Dropbox\Go-NoGo study\fMRI go-nogo\Behavioral Data\ExgaussianFitting_of_RTs\Raw_data_and_Analysis\Scripts_Zandbelt');

    inputDirectory = 'C:\Users\Tamaroy\Dropbox\Go-NoGo study\fMRI go-nogo\Behavioral Data\ExgaussianFitting_of_RTs\Raw_data_and_Analysis\scanner data';
    
    files = dir(fullfile(inputDirectory,'*.csv'));
    filenames = {files.name};
    for i=1:length(filenames)
        Subject_List{i} = strtok(filenames{i},'_');
    end
    Subject_List = unique(Subject_List);

    HeadlinesOfDataMat = {'ID','Task_Type','Block_Number','Trial_Number','StimID','IsTarget','Color','shape','time','TR','timeintoTR','Reaction_Time','Accuracy','AccuracyType'};
    n_blocks = 4; % Total number of blocks
    n = 160; % Num of trials in each block
    n_dummy = 4; % Num of dummy trials - trials to be excluded from the beginning of each block

    data = zeros(n_blocks*n,length(HeadlinesOfDataMat));
    params = zeros(length(Subject_List),11);
    
    for i=1:length(Subject_List)

        name = fullfile(inputDirectory,strcat(Subject_List{i},'*.csv'));
        file_list = dir(name); 

        % load data from csv files
        for f=1:length(file_list)
            inputfile = file_list(f);
            tmp = xlsread(fullfile(inputDirectory,inputfile.name));
            data(((f-1)*n)+1:f*n,:) = tmp(n_dummy+1:end,:); 
        end
        
        % organize data for analyses
        HF_Data = data(data(:,2)==1,:);
        LF_Data = data(data(:,2)==2,:);

        HF_Targets = HF_Data(HF_Data(:,6)==1,:);  
        HF_correctTargets = HF_Targets(HF_Targets(:,13)==1,:);
        RT_HF_targets = HF_correctTargets(:,12);
        
        LF_Targets = LF_Data(LF_Data(:,6)==1,:);
        LF_correctTargets = LF_Targets(LF_Targets(:,13)==1,:);
        RT_LF_targets = LF_correctTargets(:,12);

        % Exclude deviant RTs (4 SD from subjects' mean):
        meanRTtmp_HF = mean(RT_HF_targets); 
        stdRTtmp_HF = std(RT_HF_targets);
        upplim = meanRTtmp_HF + 4*stdRTtmp_HF;
        lowlim = meanRTtmp_HF - 4*stdRTtmp_HF;
        RT_HF_targets = RT_HF_targets(find((RT_HF_targets<upplim) & (RT_HF_targets>lowlim)),:); 
        
        meanRTtmp_LF = mean(RT_LF_targets); 
        stdRTtmp_LF = std(RT_LF_targets);
        upplim = meanRTtmp_LF + 4*stdRTtmp_LF;
        lowlim = meanRTtmp_LF - 4*stdRTtmp_LF;
        RT_LF_targets = RT_LF_targets(find((RT_LF_targets<upplim) & (RT_LF_targets>lowlim)),:); 
        
        csvwrite(strcat(Subject_List{i},'_scanner_HF_cleanRTs.csv'),RT_HF_targets);
        csvwrite(strcat(Subject_List{i},'_scanner_LF_cleanRTs.csv'),RT_LF_targets);
    end
    
        % Run exgaussian fitting and save resulting parameters and figures
        % of fit
        [X,fVal,exitFlag,solverOutput] = exgauss_fit(RT_HF_targets);
        params(i,3:5) = X;
        params(i,6) = fVal;
        figure;hold on
        exgauss_plot('both',RT_HF_targets,X,fullfile(outputFolder,strcat(Subject_List{i},'_HF')));
        close(gcf)
        
        clearvars X fVal exitFlag solverOutput
        
        [X,fVal,exitFlag,solverOutput] = exgauss_fit(RT_LF_targets);
        params(i,8:10) = X;
        params(i,11) = fVal;
        figure;hold on
        exgauss_plot('both',RT_LF_targets,X,fullfile(outputFolder,strcat(Subject_List{i},'_LF')));
        close(gcf)
              
        % Organize Output Mat
        params(i,1) = str2num(strtok(inputfile.name,'_'));
        params(i,2) = length(RT_HF_targets);
        params(i,7) = length(RT_LF_targets);
        
        clearvars   X fVal exitFlag solverOutput ...
                    name file_list inputfile tmp HF_Data HF_Targets HF_correctTargets RT_HF_targets...
                    f LF_Data LF_Targets LF_correctTargets RT_LF_targets;
    end

    % Print results to txt file
    fid=fopen(summaryFileName,'w');
    fprintf(fid,'subjectID,HF_N_trials,HF_Mu,HF_Sigma,HF_Tau,HF_minusLogLikelihood,LF_N_trials,LF_Mu,LF_Sigma,LF_Tau,LF_minusLogLikelihood\r\n');
    for j = 1:size(params,1)
        fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\r\n', params(j,:));
    end
    fclose(fid);
  
    
end
    
  






