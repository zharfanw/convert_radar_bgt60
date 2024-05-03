function [Output] = normalization(Input, lowMargin, highMargin)
%Input:enter
%lowMargin:normalized lower bound
%highMargin:normalized upper bound
%Output:output
maxValue = max(Input);
minValue = min(Input);
Output = (Input-minValue)/(maxValue-minValue);
Output = Output*(highMargin-lowMargin)+lowMargin;

end