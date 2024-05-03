function [Xs, svPrns] = sv_states(ephem, psr, time)
%sv_positions.m Provides SV states given ephemeris data and pseudoranges
%   Given ephemeris data, pseudoranges, and GPS time this function returns
%   SV position, clock bias, velocity, drift, and SV prns
% Inputs:
%    ephem  : a matrix of ephemeris data from all SVs
%    psr    : a vector of pseudoranges from all SVs [m]
%    time   : a vector of GPS time [s]
%
% Outputs:
%    Xs     : 3D matrix of all SV states [x,y,z,b,vx,vy,vz,d]
%             Matrix Structure: Xs(time,PRNs,state)
%       States:
%           [x,y,z]     : SV Positions [m]
%           [b]         : SV Clock Bias [m]
%           [vx,vy,vz]  : SV Velocities [m/s]
%           [d]         : SV Clock Drift [m/s]
%    svPrns : a vector of all visible SV PRNs
%
% Example(s):
%   Clock Drift for Time = 1 and PRN = 14
%       Xs(1,14,8)
%   XYZ Position at Time = 1 for all visible SVs
%       squeeze(Xs(1,svPrns(1,:),1:3))

    % CONSTANTS
    c = physconst('LightSpeed');

    % FUNCTION VARIABLES
    prnList = 1:min(size(psr));                     % List of all SV PRNs
    svPrns = false(length(psr), min(size(psr)));    % Matrix of seen SVs over time
    Xs = NaN(length(psr), min(size(psr)), 8);

    % CALCULATIONS
    for i = 1:length(psr)
        % Iterating over all visible SVs
        svPrns(i, ~isnan(psr(i,:))) = 1;        % List of currently visible SV PRNs
        prns = prnList(svPrns(i,:));            % List of currently visible SV PRNs
        N = sum(svPrns(i,:));                   % Number of currently visible SVs
        for j = 1:N
            tau = psr(i,prns(j)) / c;           % Message travel time [s]
            % Calculate SV state from ephemeris data
            [svPos, svB, svVel, svD] = ephem_to_sv_pos(ephem(prns(j)), time(i) - tau, tau);
            Xs(i,prns(j),:) = [svPos svB svVel svD];
        end
    end
end