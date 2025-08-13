close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Create general code of the project to run data based on only a
% folder, with little or no hard coding. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%taking in all data from the folder with Case* as the name
filenames = dir('Case*');


for i = 1:length(filenames)
    cd(filenames(i).name);

    cd ..
end

