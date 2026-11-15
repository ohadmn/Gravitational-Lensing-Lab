


% Reopen each saved figure invisibly
fig1 = openfig('part4_tau_fbl_contour_600_mes.fig', 'invisible');
fig2 = openfig('part4_umin_fbl_contour_600_mes.fig', 'invisible');
fig3 = openfig('part4_umin_tau_contour_600_mes.fig', 'invisible');

% Create a new figure for the combined plots
figure;

% Subplot 1: Copy the contents of the first figure
subplot(1, 3, 1); % 1 row, 3 columns, first subplot
ax1 = gca(fig1); % Get the axis handle of the opened figure
copyobj(allchild(ax1), gca); % Copy the contents into the subplot
title(ax1.Title.String); % Preserve the title
xlabel(ax1.XLabel.String); % Preserve the x-axis label
ylabel(ax1.YLabel.String); % Preserve the y-axis label
xlim(ax1.XLim); % Preserve x-axis limits
ylim(ax1.YLim); % Preserve y-axis limits
colormap(gca, colormap(ax1)); % Preserve the colormap
colorbar; % Add colorbar if needed
close(fig1); % Close the original figure

% Subplot 2: Copy the contents of the second figure
subplot(1, 3, 2); % Second subplot
ax2 = gca(fig2);
copyobj(allchild(ax2), gca);
title(ax2.Title.String);
xlabel(ax2.XLabel.String);
ylabel(ax2.YLabel.String);
xlim(ax2.XLim);
ylim(ax2.YLim);
colormap(gca, colormap(ax2));
colorbar;
close(fig2);

% Subplot 3: Copy the contents of the third figure
subplot(1, 3, 3); % Third subplot
ax3 = gca(fig3);
copyobj(allchild(ax3), gca);
title(ax3.Title.String);
xlabel(ax3.XLabel.String);
ylabel(ax3.YLabel.String);
xlim(ax3.XLim);
ylim(ax3.YLim);
colormap(gca, colormap(ax3));
colorbar;
close(fig3);

% Add a super title for the combined figure
sgtitle('Combined Contour Plots');
