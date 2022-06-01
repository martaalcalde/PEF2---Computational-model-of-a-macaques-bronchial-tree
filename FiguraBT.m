
load bronchis_model1_2.mat % Branches 
load lobes.mat          % Segments

% All branches in black + segments with jet colors

f = figure(1); % Create figure and set the background color to white
set(f,'Color','w');
hold on % Hold to represent the branches
axis equal % Same scale for the three axis

ax = gca; % Set the limit to the axes
ax.XLim = [-60 60];
ax.YLim = [-65 25];
ax.ZLim = [-98 50];
xlabel('Right - Left coordinate (mm)','Interpreter','latex')
ylabel('Anterior - Posterior coordinate (mm)','Interpreter','latex')
zlabel('Vertical coordinate (mm)','Interpreter','latex')

% Segments
colormap jet;
cmap = colormap;  % Colors for the plot
for k = 1:20
    if k < 12 % Right lung
        plcolor = cmap(k*floor(256/11),:);

    else  % Left lung
        plcolor = cmap((k-11)*floor(256/11),:);
    end

    patch(eval(['S' num2str(k)]),'FaceColor',plcolor,'FaceAlpha',0.1,'EdgeColor','none')
end

tic

for k = 1:size(Tubs,1)  
    branca(Tubs(k,5:7), Tubs(k,8:10), Tubs(k,3)/2, 100, 'k', .25)
    disp(['Diametre ' num2str(k) ' de 1433'])
end
timebranques = toc;

disp(['Les branques tarden ' num2str(timebranques) ' a ferse.'])

savefig(f,'./Figures/model1_2.fig')


f = openfig('./Figures/bronchisWithAngles2.fig');

print(f,'./Figures/model1_2_XY.png','-r1000','-dpng')
disp('Guardada primera figura')

view([180 0])

print(f,'./Figures/model1_2_XZ.png','-r1000','-dpng')
disp('Guardada segona figura')

view([90 0])

print(f,'./Figures/model1_2_YZ.png','-r1000','-dpng')
disp('Guardada tercera figura')






