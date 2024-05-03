function [smooth_psr] = carrier_smooth(psr, carr, M, static)
%CARRIER_SMOOTH Summary of this function goes here
%   Detailed explanation goes here
    smooth_psr = NaN(size(psr));
    smooth_psr(1,:) = psr(1,:);

    if ~exist('static', 'var')
        static = false;
    end

    if static
        for i = 2:length(psr)
            smooth_psr(i,:) = (1/M)*psr(i,:) + ((M-1)/M)*(smooth_psr(i-1,:) + ...
                (carr(i,:) - carr(i-1,:)));
        end
    else
        for i = 2:length(psr)
            for j = 1:length(psr(i,:))
                if isnan(carr(i-1,j)) && ~isnan(carr(i,j))
                    smooth_psr(i,j) = psr(i,j);
                elseif abs(carr(i,j) - carr(i-1,j)) > 1e3
                    smooth_psr(i,j) = psr(i,j);
                else
                    smooth_psr(i,j) = (1/M)*psr(i,j) + ...
                        ((M-1)/M)*(smooth_psr(i-1,j) + (carr(i,j) - carr(i-1,j)));
                end
            end
        end
    end
end

