% IFEN SX3:
% Intermediate Frequency = 5000445.88565834
% Sampling Frequency = 20000000
% Binary Data Type = 'int8'
% File Name = gpsBase_IFEN_IF.bin
% MATLAB Code to Read in Binary Signal Data
% dataSize = number of signal samples you would like to read-in

function [signalData, samplesRead] = ifen_parser(fileName)
    dataType = "int8";
    SF = 20e6; % 20000000 Hz (20 MHz)
    T = 0.001; % 1 ms of data
    dataSize = floor(T*SF);
    fileID = fopen(sprintf('%s', fileName));
    fseek(fileID, 0, 'bof');
    [signalData, samplesRead] = fread(fileID, dataSize, dataType);
end