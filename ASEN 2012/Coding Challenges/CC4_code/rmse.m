function [rm] = rmse(d2,d1)
    
    N = length(d1);
%     fprintf('n');
%     disp(N);
%     num = sum((d2-d1).^2);
%     fprintf('numerator')
%     disp(num);
    rm = sqrt(sum((d2-d1).^2)/N);
end
