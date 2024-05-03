function t = calcLocalTime(lambda_i, GPStime)
t = 4.32e4*lambda_i + GPStime;
if t >= 86400
    t = t-86400;
elseif t<0
    t = t+86400;
end
end