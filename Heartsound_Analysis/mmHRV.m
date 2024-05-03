function [IBI, MEAN, SDNN, r_MSSD] = mmHRV(heart, Tf)

%Get extreme points or Max points
extrMaxIndex = find(diff(sign(diff(heart)))==-2)+1;
extrMinIndex = find(diff(sign(diff(heart)))==+2)+1;

%To calculate IBI, the first line represents the result of the maximum value calculation, 
% and the second line represents the result of the minimum value calculation
[m, n] = size(extrMaxIndex);
for i = 1:n-1
    IBI(1,i) = (extrMaxIndex(i+1)-extrMaxIndex(i))*Tf;
end
[m, n] = size(extrMinIndex);
for i = 1:n-1
    IBI(2,i) = (extrMinIndex(i+1)-extrMinIndex(i))*Tf;
end

%Calculate MEAN
MEAN = mean(IBI, 2);

%Calculate SDNN
SDNN = std(IBI, 1, 2);

%Calculate r_MSSD
[m, n] = size(IBI);
for i = 1:n-1
    delta_RR(:,i) = IBI(:,i+1)-IBI(:,i);
end
r_MSSD = sqrt(sum(delta_RR,2).^2/(n-1));

end