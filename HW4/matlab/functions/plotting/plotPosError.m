function [] = plotPosError(title_text, Xu, mean, sigma, time)
%PLOTPOSERROR Summary of this function goes here
%   Detailed explanation goes here
    secs = seconds(time - time(1));
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    t = tiledlayout(3,1);
    title(t, title_text, 'FontSize', 20);
    nexttile(t);
    hold('on');
    grid('on');
    grid minor;
    plot(secs, Xu(1,:), 'LineWidth', 2);
    yline(mean(1) + sigma(1), '--r', 'LineWidth', 2);
    yline(mean(1) - sigma(1), '--r', 'LineWidth', 2);
    title('X-Position');
    xlabel('Time (s)');
    ylabel('Error (m)');
    ax = gca;
    ax.FontSize = 18;
    
    nexttile(t);
    hold('on');
    grid('on');
    grid minor;
    plot(secs, Xu(2,:), 'LineWidth', 2);
    yline(mean(2) + sigma(2), '--r', 'LineWidth', 2);
    yline(mean(2) - sigma(2), '--r', 'LineWidth', 2);
    title('Y-Position');
    xlabel('Time (s)');
    ylabel('Error (m)');
    ax = gca;
    ax.FontSize = 18;
    
    nexttile(t);
    hold('on');
    grid('on');
    grid minor;
    plot(secs, Xu(3,:), 'LineWidth', 2);
    yline(mean(3) + sigma(3), '--r', 'LineWidth', 2);
    yline(mean(3) - sigma(3), '--r', 'LineWidth', 2);
    title('Z-Position');
    xlabel('Time (s)');
    ylabel('Error (m)');
    ax = gca;
    ax.FontSize = 18;
    
    leg = legend('Error', '1\sigma', 'Orientation', 'horizontal');
    leg.Layout.Tile = 'north';
end

