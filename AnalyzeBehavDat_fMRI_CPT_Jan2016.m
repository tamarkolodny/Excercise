function [] = AnalyzeBehavDat_fMRI_CPT_Jan2016(DataFolder,outputFileName,format)
% format should be either 'short' (basic set of measures resported) 
%                      or 'long' (extensive report with manu additional measures)
%                      the default is 'short'.
%
% example: 
%          DataFolder = 'Y:\Experiments\fMRI_Go-nogo_Strooplike\BehavioralData\CPT_scanner_data\control';
%          Subject_List = {'1005275354','3215776500'}; 
%          outputFileName = 'tmp.txt'; 
%          format = 'long';
%          AnalyzeBehavDat_fMRI_CPT_Jan2016(DataFolder,Subject_List,outputFileName,format)
% 
%          DataFolder = 'C:\Users\Tamaroy\Dropbox\Go-NoGo study\fMRI go-nogo\Behavioral Data\ExgaussianFitting_of_RTs\Raw_data_and_Analysis\scanner data';
%          outputFileName = 'scanner_traditionalMeasures_4SDcleaning.txt'; 
%          format = 'short';

% Control_Final_Subject_List = {'1005275354','123123555','3215776500','3442650043','1941418687',...
%                                 '301111006','4161678883','3585355995','193594606','688082371',...
%                                 '3776145086','4180908251','1933321584','3991398969','2212915683',...
%                                 '3563617226','794844568','1305595106','3369278735','1671010669'};
% 
% Control_subjects_excluded_from_fMRI_analysis = ...
% {'80471670','1627945329','2380553191'};
%
% ADHD_Passive_Pretest = ...
% {'1104594723','1224900986','1234996303','1449532052','1489617418','1616783721','1780976643',...
%    '1968638271','2087417370','2319470116','2568859090','2725431208','2841447013','2945228906',...
%    '2946486023','3009342595','3055794468','3391965393','3560935931','3648076593','3804391337','4186984566'};
%    
% ADHD_Cpat_Pretest = ...
% {'81637888','130623847','214364796','1571335200','1806827974','2483101594','2707332970','2945228906',...
%   '3350238326','3350238326','3528489815','3714422859','3717015731','3781371592','4115549779'};
% 
% ADHD_Mindfulness_Pretest = ...
% {'77036516','191782068','529967520','886637171','1331405092','1505721852','1723541315','2313017966'...
% '2706350418','2889733261','3326878754','3950617229','4184373527','4192004997','4287150276'};
% 

    oldFolder = pwd;
    cd(DataFolder)
    HeadlinesOfDataMat = {'ID','Task_Type','Block_Number','Trial_Number','StimID','IsTarget','Color','shape','time','TR','timeintoTR','Reaction_Time','Accuracy','AccuracyType'};
    n_blocks = 4; % Total number of blocks
    n = 164; % Num of trials in each block
    n_dummy = 0; % Num of dummy trials - trials to be excluded from the beginning of each block
    
    files = dir(fullfile(DataFolder,'*.csv'));
    filenames = {files.name};
    for i=1:length(filenames)
        Subject_List{i} = strtok(filenames{i},'_');
    end
    Subject_List = unique(Subject_List);
    
    if strcmp(format,'long')
       output = zeros(length(Subject_List),75);
    else output = zeros(length(Subject_List),11);
    end
    
    for i=1:length(Subject_List)
        
        name = strcat(Subject_List{i},'*.csv')
        file_list = dir(name); 
        
        % load data from csv files
        data = zeros(n_blocks*n,length(HeadlinesOfDataMat));
        if length(file_list)<4 
            error('not enough input file for this subject');
        end
        for f=1:length(file_list)
            inputfile = file_list(f);
            tmp = xlsread(inputfile.name);
            data(((f-1)*n)+1:f*n,:) = tmp(n_dummy+1:end,:); 
        end
        
        HF_Data = data(data(:,2)==1,:);
        LF_Data = data(data(:,2)==2,:);

        HF_Targets = HF_Data(HF_Data(:,6)==1,:);  
        HF_NonTargets = HF_Data(HF_Data(:,6)==0,:);  
        LF_Targets = LF_Data(LF_Data(:,6)==1,:);  
        LF_NonTargets = LF_Data(LF_Data(:,6)==0,:); 
        
        HF_1stblocknum = HF_Data(1,3);
        HF_2ndblocknum = HF_Data(end,3);
        LF_1stblocknum = LF_Data(1,3);
        LF_2ndblocknum = LF_Data(end,3);

        % Compute Accuracy Measures
        HF_omissions = size(HF_Targets(HF_Targets(:,13)==0,:),1);
        HF_commissions = size(HF_NonTargets(HF_NonTargets(:,13)==0,:),1);
        LF_omissions = size(LF_Targets(LF_Targets(:,13)==0,:),1);
        LF_commissions = size(LF_NonTargets(LF_NonTargets(:,13)==0,:),1);

        % Compute RT Measures (after cleaning extreme trials)
        [HF_MeanRT,HF_stdRT,HF_Excluded] = computeRT(HF_Targets);
        [LF_MeanRT,LF_stdRT,LF_Excluded] = computeRT(LF_Targets);
        
        % Organize Output Mat
        output(i,1) = str2num(strtok(inputfile.name,'_'));
        output(i,2:5) = [HF_MeanRT HF_stdRT LF_MeanRT LF_stdRT];
        output(i,6:9) = [HF_omissions HF_commissions LF_omissions LF_commissions];
        output(i,10:11) = [HF_Excluded LF_Excluded];
        
        if strcmp(format,'long') % compute additional measures:
            
            % Compute RT of commissions
            HF_commission_trials = HF_NonTargets(HF_NonTargets(:,13)==0,:);
            if size(HF_commission_trials,1)>0
                HF_commissions_RT = mean(HF_commission_trials(:,12));
            else HF_commissions_RT = nan;
            end
            
            % Split commissions according to distractor type
            HF_commissions_SameColor = size(HF_commission_trials(HF_commission_trials(:,7)==1),1);
            HF_commissions_SameShape = size(HF_commission_trials(HF_commission_trials(:,8)==1),1);
            HF_commissions_Diff = size(HF_commission_trials(HF_commission_trials(:,7)~=1 & HF_commission_trials(:,8)~=1),1);            
            
            % Compute measures for each block seperately
            HF_omissions_1stBlock = size(HF_Targets(HF_Targets(:,13)==0 & HF_Targets(:,3)==HF_1stblocknum,:),1);
            HF_omissions_2ndBlock = size(HF_Targets(HF_Targets(:,13)==0 & HF_Targets(:,3)==HF_2ndblocknum,:),1);
            HF_commissions_1stBlock = size(HF_NonTargets(HF_NonTargets(:,13)==0 & HF_NonTargets(:,3)==HF_1stblocknum,:),1);
            HF_commissions_2ndBlock = size(HF_NonTargets(HF_NonTargets(:,13)==0 & HF_NonTargets(:,3)==HF_2ndblocknum,:),1);
            LF_omissions_1stBlock = size(LF_Targets(LF_Targets(:,13)==0 & LF_Targets(:,3)==LF_1stblocknum,:),1);
            LF_omissions_2ndBlock = size(LF_Targets(LF_Targets(:,13)==0 & LF_Targets(:,3)==LF_2ndblocknum,:),1);
            LF_commissions_1stBlock = size(LF_NonTargets(LF_NonTargets(:,13)==0 & LF_NonTargets(:,3)==LF_1stblocknum,:),1);
            LF_commissions_2ndBlock = size(LF_NonTargets(LF_NonTargets(:,13)==0 & LF_NonTargets(:,3)==LF_2ndblocknum,:),1);

            [HF_MeanRT_1stBlock,HF_stdRT_1stBlock,HF_Excluded_1stBlock] = computeRT(HF_Targets(HF_Targets(:,3)==HF_1stblocknum,:));
            [HF_MeanRT_2ndBlock,HF_stdRT_2ndBlock,HF_Excluded_2ndBlock] = computeRT(HF_Targets(HF_Targets(:,3)==HF_2ndblocknum,:));
            [LF_MeanRT_1stBlock,LF_stdRT_1stBlock,LF_Excluded_1stBlock] = computeRT(LF_Targets(LF_Targets(:,3)==LF_1stblocknum,:));
            [LF_MeanRT_2ndBlock,LF_stdRT_2ndBlock,LF_Excluded_2ndBlock] = computeRT(LF_Targets(LF_Targets(:,3)==LF_2ndblocknum,:));

            % Compute measures for first and second half of each block
            % (quartiles. could be used to examine time-on-task effects)
            HF_Data_1stQuarter = HF_Data(1:n/2,:);
            HF_Data_2ndQuarter = HF_Data(n/2+1:2*n/2,:);
            HF_Data_3rdQuarter = HF_Data(2*n/2+1:3*n/2,:);
            HF_Data_4thQuarter = HF_Data(3*n/2+1:end,:);
            LF_Data_1stQuarter = LF_Data(1:n/2,:);
            LF_Data_2ndQuarter = LF_Data(n/2+1:2*n/2,:);
            LF_Data_3rdQuarter = LF_Data(2*n/2+1:3*n/2,:);
            LF_Data_4thQuarter = LF_Data(3*n/2+1:end,:);
            
            HF_omissions_1stQuarter = size(HF_Data_1stQuarter(HF_Data_1stQuarter(:,13)==0 & HF_Data_1stQuarter(:,6)==1,:),1);
            HF_omissions_2ndQuarter = size(HF_Data_2ndQuarter(HF_Data_2ndQuarter(:,13)==0 & HF_Data_2ndQuarter(:,6)==1,:),1);
            HF_omissions_3rdQuarter = size(HF_Data_3rdQuarter(HF_Data_3rdQuarter(:,13)==0 & HF_Data_3rdQuarter(:,6)==1,:),1);
            HF_omissions_4thQuarter = size(HF_Data_4thQuarter(HF_Data_4thQuarter(:,13)==0 & HF_Data_4thQuarter(:,6)==1,:),1);
            HF_commissions_1stQuarter = size(HF_Data_1stQuarter(HF_Data_1stQuarter(:,13)==0 & HF_Data_1stQuarter(:,6)==0,:),1);
            HF_commissions_2ndQuarter = size(HF_Data_2ndQuarter(HF_Data_2ndQuarter(:,13)==0 & HF_Data_2ndQuarter(:,6)==0,:),1);
            HF_commissions_3rdQuarter = size(HF_Data_3rdQuarter(HF_Data_3rdQuarter(:,13)==0 & HF_Data_3rdQuarter(:,6)==0,:),1);
            HF_commissions_4thQuarter = size(HF_Data_4thQuarter(HF_Data_4thQuarter(:,13)==0 & HF_Data_4thQuarter(:,6)==0,:),1);
            LF_omissions_1stQuarter = size(LF_Data_1stQuarter(LF_Data_1stQuarter(:,13)==0 & LF_Data_1stQuarter(:,6)==1,:),1);
            LF_omissions_2ndQuarter = size(LF_Data_2ndQuarter(LF_Data_2ndQuarter(:,13)==0 & LF_Data_2ndQuarter(:,6)==1,:),1);
            LF_omissions_3rdQuarter = size(LF_Data_3rdQuarter(LF_Data_3rdQuarter(:,13)==0 & LF_Data_3rdQuarter(:,6)==1,:),1);
            LF_omissions_4thQuarter = size(LF_Data_4thQuarter(LF_Data_4thQuarter(:,13)==0 & LF_Data_4thQuarter(:,6)==1,:),1);
            LF_commissions_1stQuarter = size(LF_Data_1stQuarter(LF_Data_1stQuarter(:,13)==0 & LF_Data_1stQuarter(:,6)==0,:),1);
            LF_commissions_2ndQuarter = size(LF_Data_2ndQuarter(LF_Data_2ndQuarter(:,13)==0 & LF_Data_2ndQuarter(:,6)==0,:),1);
            LF_commissions_3rdQuarter = size(LF_Data_3rdQuarter(LF_Data_3rdQuarter(:,13)==0 & LF_Data_3rdQuarter(:,6)==0,:),1);
            LF_commissions_4thQuarter = size(LF_Data_4thQuarter(LF_Data_4thQuarter(:,13)==0 & LF_Data_4thQuarter(:,6)==0,:),1);

            [HF_MeanRT_1stQuarter,HF_stdRT_1stQuarter,HF_Excluded_1stQuarter] = computeRT(HF_Data_1stQuarter(HF_Data_1stQuarter(:,6)==1,:));
            [HF_MeanRT_2ndQuarter,HF_stdRT_2ndQuarter,HF_Excluded_2ndQuarter] = computeRT(HF_Data_2ndQuarter(HF_Data_2ndQuarter(:,6)==1,:));
            [HF_MeanRT_3rdQuarter,HF_stdRT_3rdQuarter,HF_Excluded_3rdQuarter] = computeRT(HF_Data_3rdQuarter(HF_Data_3rdQuarter(:,6)==1,:));
            [HF_MeanRT_4thQuarter,HF_stdRT_4thQuarter,HF_Excluded_4thQuarter] = computeRT(HF_Data_4thQuarter(HF_Data_4thQuarter(:,6)==1,:));
            [LF_MeanRT_1stQuarter,LF_stdRT_1stQuarter,LF_Excluded_1stQuarter] = computeRT(LF_Data_1stQuarter(LF_Data_1stQuarter(:,6)==1,:));
            [LF_MeanRT_2ndQuarter,LF_stdRT_2ndQuarter,LF_Excluded_2ndQuarter] = computeRT(LF_Data_2ndQuarter(LF_Data_2ndQuarter(:,6)==1,:));
            [LF_MeanRT_3rdQuarter,LF_stdRT_3rdQuarter,LF_Excluded_3rdQuarter] = computeRT(LF_Data_3rdQuarter(LF_Data_3rdQuarter(:,6)==1,:));
            [LF_MeanRT_4thQuarter,LF_stdRT_4thQuarter,LF_Excluded_4thQuarter] = computeRT(LF_Data_4thQuarter(LF_Data_4thQuarter(:,6)==1,:));

            % Organize Output Mat
            output(i,12:15) = [HF_commissions_RT HF_commissions_SameColor HF_commissions_SameShape HF_commissions_Diff];
            output(i,16:23) = [HF_omissions_1stBlock HF_omissions_2ndBlock HF_commissions_1stBlock HF_commissions_2ndBlock ...
                               LF_omissions_1stBlock LF_omissions_2ndBlock LF_commissions_1stBlock LF_commissions_2ndBlock];
            output(i,24:35) = [HF_MeanRT_1stBlock HF_stdRT_1stBlock HF_Excluded_1stBlock HF_MeanRT_2ndBlock HF_stdRT_2ndBlock HF_Excluded_2ndBlock ...
                               LF_MeanRT_1stBlock LF_stdRT_1stBlock LF_Excluded_1stBlock LF_MeanRT_2ndBlock LF_stdRT_2ndBlock LF_Excluded_2ndBlock];
            output(i,36:51) = [HF_omissions_1stQuarter HF_omissions_2ndQuarter HF_omissions_3rdQuarter HF_omissions_4thQuarter HF_commissions_1stQuarter HF_commissions_2ndQuarter HF_commissions_3rdQuarter HF_commissions_4thQuarter ...
                               LF_omissions_1stQuarter LF_omissions_2ndQuarter LF_omissions_3rdQuarter LF_omissions_4thQuarter LF_commissions_1stQuarter LF_commissions_2ndQuarter LF_commissions_3rdQuarter LF_commissions_4thQuarter];
            output(i,52:75) = [HF_MeanRT_1stQuarter HF_stdRT_1stQuarter HF_Excluded_1stQuarter HF_MeanRT_2ndQuarter HF_stdRT_2ndQuarter HF_Excluded_2ndQuarter HF_MeanRT_3rdQuarter HF_stdRT_3rdQuarter HF_Excluded_3rdQuarter HF_MeanRT_4thQuarter HF_stdRT_4thQuarter HF_Excluded_4thQuarter ...
                               LF_MeanRT_1stQuarter LF_stdRT_1stQuarter LF_Excluded_1stQuarter LF_MeanRT_2ndQuarter LF_stdRT_2ndQuarter LF_Excluded_2ndQuarter LF_MeanRT_3rdQuarter LF_stdRT_3rdQuarter LF_Excluded_3rdQuarter LF_MeanRT_4thQuarter LF_stdRT_4thQuarter LF_Excluded_4thQuarter];
                
            clearvars   HF_commission_trials HF_commissions_RT HF_commissions_SameColorHF_commissions_SameShape HF_commissions_Diff ...
                        HF_omissions_1stBlock HF_omissions_2ndBlock HF_commissions_1stBlock HF_commissions_2ndBlock ...
                        LF_omissions_1stBlock LF_omissions_2ndBlock LF_commissions_1stBlock LF_commissions_2ndBlock ...
                        HF_MeanRT_1stBlock HF_stdRT_1stBlock HF_Excluded_1stBlock HF_MeanRT_2ndBlock HF_stdRT_2ndBlock HF_Excluded_2ndBlock ...
                        LF_MeanRT_1stBlock LF_stdRT_1stBlock LF_Excluded_1stBlock LF_MeanRT_2ndBlock LF_stdRT_2ndBlock LF_Excluded_2ndBlock ...
                        HF_Data_1stQuarter HF_Data_2ndQuarter HF_Data_3rdQuarter HF_Data_4thQuarter LF_Data_1stQuarter LF_Data_2ndQuarter LF_Data_3rdQuarter LF_Data_4thQuarter ...
                        HF_omissions_1stQuarter HF_omissions_2ndQuarter HF_omissions_3rdQuarter HF_omissions_4thQuarter HF_commissions_1stQuarter HF_commissions_2ndQuarter HF_commissions_3rdQuarter HF_commissions_4thQuarter ...
                        LF_omissions_1stQuarter LF_omissions_2ndQuarter LF_omissions_3rdQuarter LF_omissions_4thQuarter LF_commissions_1stQuarter LF_commissions_2ndQuarter LF_commissions_3rdQuarter LF_commissions_4thQuarter ...          
                        HF_MeanRT_1stQuarter HF_stdRT_1stQuarter HF_Excluded_1stQuarter HF_MeanRT_2ndQuarter HF_stdRT_2ndQuarter HF_Excluded_2ndQuarter HF_MeanRT_3rdQuarter HF_stdRT_3rdQuarter HF_Excluded_3rdQuarter HF_MeanRT_4thQuarter HF_stdRT_4thQuarter HF_Excluded_4thQuarter ...
                        LF_MeanRT_1stQuarter LF_stdRT_1stQuarter LF_Excluded_1stQuarter LF_MeanRT_2ndQuarter LF_stdRT_2ndQuarter LF_Excluded_2ndQuarter LF_MeanRT_3rdQuarter LF_stdRT_3rdQuarter LF_Excluded_3rdQuarter LF_MeanRT_4thQuarter LF_stdRT_4thQuarter LF_Excluded_4thQuarter
        
        end

        clearvars   data name file_list inputfile tmp HF_commissions HF_Data HF_NonTargets HF_omissions...
                    HF_Targets f LF_commissions LF_Data HF_NonTargets LF_omissions LF_Targets;
        
    end
    
    fid=fopen(outputFileName,'w');
    if strcmp(format,'long')
        fprintf(fid,'subjectID,HF_meanRT,HF_stdRT,LF_meanRT,LF_stdRT,HF_omissions,HF_commissions,LF_omissions,LF_commissions,HF_Excluded,LF_Excluded,HF_commissions_RT,HF_commissions_SameColor,HF_commissions_SameShape,HF_commissions_Diff,');
        fprintf(fid,'HF_omissions_1stBlock,HF_omissions_2ndBlock,HF_commissions_1stBlock,HF_commissions_2ndBlock,LF_omissions_1stBlock,LF_omissions_2ndBlock,LF_commissions_1stBlock,LF_commissions_2ndBlock,');
        fprintf(fid,'HF_MeanRT_1stBlock,HF_stdRT_1stBlock,HF_Excluded_1stBlock,HF_MeanRT_2ndBlock,HF_stdRT_2ndBlock,HF_Excluded_2ndBlock,LF_MeanRT_1stBlock,LF_stdRT_1stBlock,LF_Excluded_1stBlock,LF_MeanRT_2ndBlock,LF_stdRT_2ndBlock,LF_Excluded_2ndBlock,');
        fprintf(fid,'HF_omissions_1stQuarter,HF_omissions_2ndQuarter,HF_omissions_3rdQuarter,HF_omissions_4thQuarter,HF_commissions_1stQuarter,HF_commissions_2ndQuarter,HF_commissions_3rdQuarter,HF_commissions_4thQuarter,');
        fprintf(fid,'LF_omissions_1stQuarter,LF_omissions_2ndQuarter,LF_omissions_3rdQuarter,LF_omissions_4thQuarter,LF_commissions_1stQuarter,LF_commissions_2ndQuarter,LF_commissions_3rdQuarter,LF_commissions_4thQuarter,');
        fprintf(fid,'HF_MeanRT_1stQuarter,HF_stdRT_1stQuarter,HF_Excluded_1stQuarter,HF_MeanRT_2ndQuarter,HF_stdRT_2ndQuarter,HF_Excluded_2ndQuarter,HF_MeanRT_3rdQuarter,HF_stdRT_3rdQuarter,HF_Excluded_3rdQuarter,HF_MeanRT_4thQuarter,HF_stdRT_4thQuarter,HF_Excluded_4thQuarter,');
        fprintf(fid,'LF_MeanRT_1stQuarter,LF_stdRT_1stQuarter,LF_Excluded_1stQuarter,LF_MeanRT_2ndQuarter,LF_stdRT_2ndQuarter,LF_Excluded_2ndQuarter,LF_MeanRT_3rdQuarter,LF_stdRT_3rdQuarter,LF_Excluded_3rdQuarter,LF_MeanRT_4thQuarter,LF_stdRT_4thQuarter,LF_Excluded_4thQuarter\r\n');
        for j = 1:size(output,1)
             fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,\r\n', output(j,:));
        end
    else fprintf(fid,'subjectID,HF_meanRT,HF_stdRT,LF_meanRT,LF_stdRT,HF_omissions,HF_commissions,LF_omissions,LF_commissions,HF_Excluded,LF_Excluded\r\n');
         for j = 1:size(output,1)
             fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\r\n', output(j,:));
         end
    end
    fclose(fid);
    cd(oldFolder)
end
   
  

function [meanRT,stdRT,excluded] = computeRT(datamat) 
    CorrectTrialsData = datamat(datamat(:,13)==1,:); %exclude errors
    numCorrectTrials = length(CorrectTrialsData);
    CorrectTrialsData = CorrectTrialsData(CorrectTrialsData(:,12)>100,:); %exclude very short RTs
    CorrectTrialsData = CorrectTrialsData(CorrectTrialsData(:,12)<3000,:); %exclude very long RTs
    meanRTtmp = mean(CorrectTrialsData(:,12)); 
    stdRTtmp = std(CorrectTrialsData(:,12));
    upplim = meanRTtmp + 4*stdRTtmp;
    lowlim = meanRTtmp - 4*stdRTtmp;
    CorrectTrialsData = CorrectTrialsData(find((CorrectTrialsData(:,12)<upplim) & (CorrectTrialsData(:,12)>lowlim)),:); %exclude trials with deviant RT, relative to one owns data
    excluded = numCorrectTrials - length(CorrectTrialsData);
    meanRT = mean(CorrectTrialsData(:,12));
    stdRT = std(CorrectTrialsData(:,12));
end







