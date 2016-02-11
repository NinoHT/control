% R2* corrected with its volume
clc, close all, clear
% load ~/Dropbox/ccsvi_swi/DGM_analysis/matlab/volumes.mat
% load ~/Dropbox/ccsvi_swi/DGM_analysis/matlab/R2starresults_10_02_2015_2px.mat
addpath ~/Dropbox/ccsvi_swi/DGM_analysis/matlab/support/

load('~/Dropbox/Matlab/CCSVI-1/R2s_in_DGM_ALL.mat')

% % Volumes of the structures and Volumes normalized
% The volumes of each structures were accounted (in voxels)
volB = DGM.Volumes.Brain;
Vdgm = [DGM.Volumes.Thal', DGM.Volumes.Caud', DGM.Volumes.Puta',...
    DGM.Volumes.Pall'];                         % Volumes not normalized

                                 
Vdgm_n = Vdgm./repmat(volB', 1, 4)*1000;    % Volumes normalized

indMS = find(DGM.vec_info(:, 1) ==1);
indSB = find(DGM.vec_info(:, 2) == 2 & DGM.vec_info(:, 1) ~= 1 );
indHC = find(DGM.vec_info(:, 2) == 3);
%
% =========================================================================
% % Statistics
% =========================================================================
% 
% %
set(0,'DefaultTextInterpreter', 'latex')
titles = {'$\it{Thalamus}$'; '$\it{Caudate}$'; '$\it{Putamen}$';'$\it{Pallidus}$'};
labels1 = [repmat('MS', length(indMS), 1);...
    repmat('HC', length(indHC), 1)];
% %
k2 = 0;
for k1=1:4
    k2 = k2+1;
    figure(1)
    clf
    data1 = [Vdgm_n(indMS, k1); Vdgm_n(indHC,k1)];
    bpwm = boxplot(data1, labels1);
    [yDb(:,:,k2), yDm(:,:,k2), yDw(:,:,k2), yDa(:,:,k2), yDoy.(['ol_' (num2str(k2))])] = boxplot_parts(bpwm);
    [p1, table1, stats1] = kruskalwallis(data1, labels1, 'off');
   StatS(k1) = stats1;
   Table.(['t' num2str(k2)]) = table1;
 
end
%
close
%% Plot the new box plots
close all
figure(2)
clf
% %
set(gcf, 'Position', [450    80   820   615])
set(gcf, 'Units', 'Points')
% % 
c = {[0.6 0.8 0.6]; [0.8 0.6 0.6]; [0.6 0.6 0.8]};
hold on
counter = 0;
for k1 =1:4
    for k2 = 1:2
        counter = counter + 1;
        xAxis(counter) = 0.65*k1+0.25*(k2-2);
        % Boxes
        REC(counter) = rectangle('Position', [xAxis(counter) yDb(1,k2,k1) 0.25 diff(yDb(1:2,k2,k1))], ...
            'FaceColor', c{k2});
        % Wiskers
        plot([1 1]*xAxis(counter)+0.25/2, [yDw(1,k2,k1) yDw(2,k2,k1)], '--k',...
             [1 1]*xAxis(counter)+0.25/2, [yDw(3,k2,k1) yDw(4,k2,k1)], '--k',...
             xAxis(counter)+0.25/2 + 0.25/5*[-1 1], yDw(1,k2,k1)*[1 1], '-k', ...
             xAxis(counter)+0.25/2 + 0.25/5*[-1 1], yDw(4,k2,k1)*[1 1], '-k')         
        % Mean
        plot([xAxis(counter) xAxis(counter)+0.25], yDm(:,k2, k1), '-r')
        % Outliers
        plot(xAxis(counter)+0.25/2, yDoy.(['ol_' (num2str(k1))]){k2}, '+r')        
    end
end
% Trick for labels
p1=plot(nan,nan,'s', 'MarkerSize', 15,'MarkerEdgeColor','k',...
    'MarkerFaceColor',get(REC(1),'facecolor'));
p2=plot(nan,nan,'s', 'MarkerSize', 15,'MarkerEdgeColor','k',...
    'MarkerFaceColor',get(REC(2),'facecolor'));
% XLimits
xlim([xAxis(1)-0.25 xAxis(end)+0.5])
ylim([1 10])
% Significance
% % Significance in Thalamus
plot(xAxis([1 1 2 2])+0.25/2, yDw(4,2,1)*[1.025 1.05 1.05 1.025],'-k', 'Color',...
    [0.3 0.3 0.5], 'LineWidth', 1.5)
plot(xAxis(2), (yDw(4,2,1)+0.6), 'dk', 'MarkerSize', 6,...
    'MarkerFaceColor', 'k')
% % Significance in Caudate
plot(xAxis([1 1 2 2]+2)+0.25/2, yDw(4,1,2)*[1.05 1.1 1.1 1.05],'-k', 'Color',...
    [0.3 0.3 0.5], 'LineWidth', 1.5)
plot(0.985*xAxis(4)*[1 1.03], (yDw(4,1,2)+0.6)*[1 1], '*k', 'MarkerSize', 6,...
    'MarkerFaceColor', 'k')
% % % Significance in Pallidus
plot(xAxis([1 1 2 2]+6)+0.25/2, yDw(4,2,4)*[1.05 1.1 1.1 1.05], 'Color',...
    [0.3 0.3 0.5], 'LineWidth', 1.5)
plot(xAxis(8), (yDw(4,2,4)+0.6), '*k', 'MarkerSize', 6,...
    'MarkerFaceColor', 'k')


% Edition
box on
set(gca, 'XTick', xAxis(2:2:end), 'XTickLabel', ' ', 'FontSize', 20, 'Units', 'Points')
legend([p1 p2], 'MS', 'HC', 'Location', 'Best')
ylabel('Volume [$\times10^{-3}\,a.u.$]', 'Interpreter', 'Latex', 'FontSize', 24)


Xpos = [0 113 208 315];
Ypos = [-47 -42 -44 -41];


for k=1:size(titles, 1)
    t1X(k) = text(Xpos(k), Ypos(k), titles{k}, 'FontSize', 20, 'Interpreter',...
        'Latex', 'Units', 'Points', 'EdgeColor', 'none', 'Rotation', 25);
end
%%
Xpos = [0 113 208 315];
Ypos = [-47 -42 -44 -41];

for k=1:4
    
    set(t1X(k), 'Pos', [Xpos(k) Ypos(k)],'Rotation', 25)
end
%
set(gca, 'Pos', [112 70 465 400])


%%
% This is a new comment
% =========================================================================
%% Stats and print in screen
% =========================================================================
clc
thr0 = 0.1;

fprintf('Significance Volume normalized\n')
a1 = sprintf('Thal\tCaud\tPut\tPall');
fprintf('%s\n%s\n%s\n', repmat('=', 1, size(a1,2)*2), a1, repmat('=', 1, size(a1,2)*2))
fprintf('%1.3f\t%1.3f\t%1.3f\t%1.3f\n', Table.t1{2,6}, Table.t2{2,6}, Table.t3{2,6}, Table.t4{2,6})

namStr = {'Thal'; 'Caud'; 'Put'; 'Pall'};
namGrp = {'MS'; 'HC'};

% %
for k1 =1:4
    c0 = multcompare(StatS(k1), 'Ctype', 'bonferroni', 'Display', 'off');
    inP = find(c0(:,6)< thr0, 6);
    if ~isempty(inP)
        fprintf('%s\n', namStr{k1})
        for k2 = 1:length(inP)
            fprintf('%s -> %s\t(%1.3f)\n', namGrp{c0(inP(k2),1)}, namGrp{c0(inP(k2),2)}, ...
                c0(inP(k2),6))
        end
    end
end

