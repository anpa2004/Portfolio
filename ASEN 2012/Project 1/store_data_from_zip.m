function [data] = store_data_from_zip(file)
    unzip(file,'Project_1');
    
    for i = 11:20
        data(i)=load('testrun',i,'.mat');



end