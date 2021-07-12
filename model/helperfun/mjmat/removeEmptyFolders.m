function removeEmptyFolders

% thisDir = '/fmri2/PI/liu/quantum/liulab/experiments/';
allDir = dir(thisDir);
for i = 1:length(allDir)
    if strcmp(allDir(i).name(1),'.')
        continue
    end
    check = dir([thisDir,allDir(i).name]);
    if length(check)==2 && strcmp(check(1).name,'.') && strcmp(check(2).name,'..')
        rmdir([thisDir,allDir(i).name]);
    end
end
disp('check');