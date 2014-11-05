addpath(genpath(fullfile(pwd,'..','Projet_KI/')));

folderPath_LAB = '/Users/gong/Documents/projet_KI/train_set/fricative/friclab';

allFiles_LAB = dir( folderPath_LAB );
allNames_LAB = {allFiles_LAB.name};

ind2remove = cellfun(@(x)x(1) == '.',allNames_LAB);
allNames_LAB(:,ind2remove) =[];

for ii = 1 : length(allNames_LAB)
   writename = [allNames_LAB{ii}(1:end-4) 'w.lab'];
   audiosulptlab_convert( allNames_LAB{ii}, writename,2);    % type 2 
    
end

