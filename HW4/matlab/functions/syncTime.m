function [out1,out2,syncedTime] = syncTime(vec1,time1,vec2,time2)
%SYNCTIME Summary of this function goes here
%   Detailed explanation goes here
time_off = time2(1) - time1(1);
if time_off > 0
    n2 = 0;
    n1 = time_off;
    N = length(vec2) - n1;
    syncedTime = time1((1+n2):N-n1);
else
    n2 = abs(time_off);
    n1 = 0;
    N = length(vec1) - n2;
    syncedTime = time2((1+n1):N-n2);
end

out1 = vec1((1+n1):N-n2,:);
out2 = vec2((1+n2):N-n1,:);
end

