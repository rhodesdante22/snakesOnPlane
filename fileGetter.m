function [vids, fileListTables, tracingsListTables] = fileGetter()

f = fopen("training.txt");
trainingFilenames = fscanf(f,"%c");
namesArray = split(trainingFilenames);

% chose some pseudorandom files to run quick experiments on
C = cell(10,1);
for i=1:10
    C{i} = namesArray{i*500};
end


vids = cell(10,1);
vidPath = "/EchoNet-Dynamic/Videos/";
myDir = pwd;
filePath = fullfile(myDir, vidPath);
myFiles = dir(fullfile(filePath));
for i=1:10
    for k = 1:length(myFiles)
        baseFileName = myFiles(k).name;
        fullFileName = fullfile(filePath, baseFileName);
        contents = split(fullFileName, "s\");  % This line and the next work on MY file system. 
        localName = contents{3};
        if strcmp(localName, append(C{i}, ".avi"))
            disp(fullFileName)
            vids{i} = VideoReader(fullFileName);
            continue;
        end
    end
end
% the cell array "vids" now contains the video files 

%%
fileListTables = cell(10,1);
tracingsListTables = cell(10,1);
FileListCSV = readtable("EchoNet-Dynamic/FileList.csv");
VolumeTracingsCSV = readtable("EchoNet-Dynamic/VolumeTracings.csv");
FileListNames = (FileListCSV.FileName);

VolumeTracingNames = (VolumeTracingsCSV.FileName);
S = sprintf('%s*', VolumeTracingNames{:});
ar = split(S, ".avi*");
%% volumeTracings
for i=1:10
    currentFile = C{i};
    %convertedName = str2num(currentFile);
    rowNum = 1;
    while rowNum < height(VolumeTracingsCSV)
        if(strcmp(ar{rowNum},currentFile))
            init = rowNum;
            while(strcmp(ar{rowNum},currentFile))
               rowNum = rowNum + 1; 
            end
            fin = rowNum-1;
            tracingsListTables{i} = VolumeTracingsCSV(init:fin, :);
            continue;
        end
        rowNum = rowNum + 1;
    end
end
%% fileList
for i=1:10
    currentFile = C{i};
    convertedName = str2num(currentFile);
    rowNum = 1;
    while rowNum < height(FileListCSV)
        if(FileListNames(rowNum) == convertedName)
            fileListTables{i} = FileListCSV(rowNum, :);
        end
        rowNum = rowNum + 1;
    end
    
end

end
