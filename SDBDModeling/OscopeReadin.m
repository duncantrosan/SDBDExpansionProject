function [Time,Channel1,Channel2,Channel3,Channel4] = OscopeReadin(path)
%% This function is meant to read in data from the tektronix Oscilliscopes that we have in lab
% It deleates any Nans in the data as well as infinitys that sometimes
% appear in current measurments for reasons that I cannot explain. 
% path is a diretory not a file!


filelist = dir(path);
A = size(filelist);
filenamescell = struct2cell(filelist);
DSstore = 2;
[~,scheck]= size(filenamescell);
if scheck >= 3
logicDS = strcmp(filenamescell{1,3},'.DS_Store');
% Sometimes Importing files into Macs leave DSstore file that screws things
% up if you dont take it into account.

if logicDS == 1
    DSstore = 3;
end
end
Nfiles = A(1,1) - DSstore;

for i = DSstore+1:Nfiles + DSstore
    fullfiles{1,i-DSstore} = fullfile(filenamescell{2,i},filenamescell{1,i});
    filenames{i-DSstore,1} = filenamescell{1,i};
end
fullfiles = natsort(fullfiles);
clear path logicDS filelist A filenamescell i DSstore


%% Read in data


% Define Options for readin file
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [22, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Time", "Channel1", "Channel2","Channel3","Channel4"];
opts.VariableTypes = ["double", "double", "double","double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";



ignoredfiles = 0;
disp('Files Read In')
for i = 1:Nfiles
    % Import the data
    readinname = fullfiles{1,i};
    [path,name,ext] = fileparts(readinname);
    logicfiles = strcmp(ext,'.csv');
    if logicfiles == 1
    inputtable = readtable(readinname,opts);
    
    % Convert to output type
    tempmat = table2array(inputtable);
    templength = length(tempmat);
    tempsize = size(tempmat);
    
    % Finds 
    for j = 1:templength
        Time(i-ignoredfiles,j) = tempmat(j,1);
        Channel1(i-ignoredfiles,j) = tempmat(j,2);
        if tempsize(1,2) >= 3
        Channel2(i-ignoredfiles,j) = tempmat(j,3);
        end
         if tempsize(1,2) >= 4
        Channel3(i-ignoredfiles,j) = tempmat(j,4);
         end
        if tempsize(1,2) >= 5
        Channel4(i-ignoredfiles,j) = tempmat(j,5);
        end
    end
    end
    if logicfiles == 0
        ignoredfiles = ignoredfiles +1;
    end
    clear tempmat templength
    disp(readinname)
end
Nfiles = Nfiles - ignoredfiles;
clear readinname fullfiles logicfiles inputtable tempsize opts path name ext




dim = size(Time);
Npoints = dim(1,2);
Nfiles = dim(1,1);
%% Need to resutructure this list to insure theres no infs (oscopes interpret clipping as a inf)
for k = 1:Nfiles
    for n = 1:Npoints
        logicinf1 = isinf(Channel1(k,n));
        logicinf2 = isinf(Channel2(k,n));
        logicinf3 = isinf(Channel3(k,n));
        logicinf4 = isinf(Channel4(k,n));
%         logicinf3 = isinf(Channel3(k,n));
        
        j = 1;
        while logicinf1 == 1
            Channel1(k,n) = (Channel1(k,n-j)+Channel1(k,n+j))/2;
            j = j+1;
            logicinf1 = isinf(Channel1(k,n));
        end
         j = 1;
        while logicinf2 == 1
            Channel2(k,n) = (Channel2(k,n-j)+Channel2(k,n+j))/2;
            j = j+1;
            logicinf2 = isinf(Channel2(k,n));
        end  
         j = 1;
        while logicinf3 == 1
            Channel3(k,n) = (Channel3(k,n-j)+Channel3(k,n+j))/2;
            j = j+1;
            logicinf3 = isinf(Channel3(k,n));
        end 
         j = 1;
          while logicinf4 == 1
            Channel4(k,n) = (Channel4(k,n-j)+Channel4(k,n+j))/2;
            j = j+1;
            logicinf4 = isinf(Channel4(k,n));
        end  
        
        
    end
end
clear dim logicinf1 logicinf2 i j k n
%% Restructure to deleate all NANs
for k = 1:Nfiles
    for n = 1:Npoints
        if isnan(Channel1(k,n)) == 1 || isnan(Channel2(k,n)) == 1
            Channel1(:,n) = [];
            Channel2(:,n) = [];
            Channel3(:,n) = [];
            Channel4(:,n) = [];
            Time(:,n) = [];
            Npoints = Npoints-1;
    end
    end
end


