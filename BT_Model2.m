clear all, close all

% First we load data from the lungs:
%   -  Xarxa: is a matrix obtained from the previous one and addition data.
%      The first three columns of the matrix contain cartesina coordinates
%      for each point of 2, this has been computed with addition data that
%      stablished the ditance in milimeteres between points in 1, in the x,
%      y and z directions. The fourth column of this matrix states numbers
%      1 or 2 if it correspond to the left or right lung, and the fifth
%      column has numbers from 1 to 20 for each lobe.    
load xarxa.mat    

% Separate the matrix Xarxa into a the column stablishing the segments and 
% the another matrix containing the cooridnates for each point and a fourth
% column to identify the points.
SE = Xarxa(:,5);      % Segments numbered (1 to 20)
Xarxa = Xarxa(:,1:4); % Forth colum 1 or 2 for right and left lung respectively

% General parameters
M = 3; % Proportional factor for diameter-length of the branload ("Segments_Macacs_SenseFiltres.mat",'segments');ch

% TRACHEA
dTraq = 8.1; % Tracha diameter

% Matrix contaning all branches: each row is a bronchial branch
% Meaning of each column:
%  1: How many branches there are between this one and the trachea
%  2: In which row is its mother branch
%  3: Diameter
%  4: Longitude
%  5: xi, initial x coordinate
%  6: yi, initial y coordinate
%  7: zi, initial z coordinate
%  8: xf, final x coordinate
%  9: yf, final y coordinate
% 10: zf, final z coordinate
% 11: birth, number of daughter branches emerging from it
% 12: vx, x coordinate of its mother's vector
% 13: vy, y coordinate of its mother's vector
% 14: vz, z coordinate of its mother's vector
% 15: angle, angle of bifrucation of its mother's branch 
% 16: lobe, integer refering to the lobe they irrigate, if they are leading
%     branches it is 0
Tubs = [0 0 dTraq 3*dTraq 0 0 3*dTraq 0 0 0 2 -1 0 0 0 0];

mothers = [    1,    2,   4,    5,  11,    9,  15,  16,     3,  21,    22,  25,    27,  31,  31];

section = { 1:11,  1:3,   1, 7:10,   4,  7:8,   7,   9, 12:13,  12, 14:15,  14, 17:18,  17,  19
           12:20, 4:11,   2,    6,   5, 9:10,   8,  10, 14:20,  13,    16,  15, 19:20,  18,  20
             nan,  nan,   3,  4:5, nan,  nan, nan, nan,   nan, nan, 17:20, nan,   nan, nan, nan
             nan,  nan, nan,   11, nan,  nan, nan, nan,   nan, nan,   nan, nan,   nan, nan, nan};

for k = 1:size(section,2)
    m = mothers(k);
    if all(isnan([section{3,k},section{4,k}])) % Bifrucation ------------
        [T1, T2] = tubes(SE,Xarxa,Tubs(m,:),m,section{1,k},section{2,k});
        Tubs = [Tubs; T1; T2];

    elseif isnan(section{4,k}) % Trifurcation ----------------------------
        [T1, T2, T3] = tubes(SE,Xarxa,Tubs(m,:),m,section{1,k},section{2,k},section{3,k});
        Tubs = [Tubs; T1; T2; T3];

    else % Quatrifurcation ------------------------------------------------
        [T1, T2, T3, T4] = tubes(SE,Xarxa,Tubs(m,:),m,section{1,k},section{2,k},section{3,k},section{4,k});
        Tubs = [Tubs; T1; T2; T3; T4];
    end
end

% Representation
load lobes.mat % Patch object for each lobe (to save computational time)
lob = Tubs(:,end) ~= 0; % Indexes of terminal branches
Terminal = Tubs(lob,:);
Others   = Tubs(~lob,:);
% Order in plotting the lobs depends on when the final branch has been
% computed:
seg = Tubs(lob,end);

figure(1);
hold on, axis equal
colormap jet;
cmap = colormap;  % Colors for the plot

for k = 1:20
    if k < 12
        plcolor = cmap(k*floor(256/11),:);
    else 
        plcolor = cmap((k-11)*floor(256/11),:);
    end
    patch(eval(['S' num2str(seg(k))]),'FaceColor',plcolor,'FaceAlpha',0.1,'EdgeColor','none')
    branca(Terminal(k,5:7),Terminal(k,8:10), Terminal(k,3)/2, 100, plcolor, .75);
end

for k = 1:size(Others,1)
    branca(Others(k,5:7),Others(k,8:10), Others(k,3)/2, 100, 'k', .25)
end
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')
save('BT_model2.mat','Tubs','-mat')
length(Tubs)