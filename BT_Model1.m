clear all, close all
load xarxa.mat
load ("Segments_Macacs_SenseFiltres.mat",'segments');

% Tracha and first bifuraction (left and right lungs) 
dTraq = 7.9; % Tracha diameter

% Matrix contaning all branches: each row is a bronchial branch
Tubs = [0 0 dTraq 3*dTraq 0 0 3*dTraq 0 0 0 2 -1 0 0 0];

SE = Xarxa(:,5); % Segments numbered (1 to 20)
Xarxa = Xarxa(:,1:4); % Forth colum 1 or 2 for right and left lung respectively

% General parameters
M = 3; % Proportional factor for diameter-length of the branch
n = 3; % Optimum factor for the function Angles

% First bifurcation 
N1 = sum(Xarxa(:,4)==1); % Volume right lung
N0 = sum(Xarxa(:,4)==1) + sum(Xarxa(:,4)==2); % Total volume

d1 = dTraq*(N1/N0)^(1/3); % Diameter for the branch irrigating the right lung
d2 = dTraq*(sum(Xarxa(:,4)==2)/size(Xarxa,1))^(1/3); % Diameter for the left one

l1 = M*d1; % Longitud of the right branch
l2 = M*d2; % Left branch

% fill1 and fill2 define the value in hte fourth colume of Xarxa that
% reference the points of the irrigated volume for each branch
fill1 = 1; 
fill2 = 2;
X     = Xarxa;
xf    = [0 0 0];  % End point of the mother branch
u     = [0 -1 0]; % Unitary orthogonal vector to the bifurcation plane
v     = [0 0 -1]; % Unitary vector in the directon of the mother's branch
vs    = [-1 0 0]; % Unitary orthogonal vector to the separation plane
[cos1,sin1,cos2,sin2,sin11,cos11,sin22,cos22] = anglesnew(fill1,fill2,X,xf,n,u,v,vs);
% cos and sin with one number reference the direction in the bifrucation
% plane, following the axes defined by the mother branch and the separation
% plane. The ones with two numbers define the direction in the separation
% plane, with respect to the axes defined by the mother branch and the
% bifurcation plane.

dir1 = v*(cos1+cos11)+vs*sin1+u*sin11; dir1 = dir1/norm(dir1); % direcction daugther 1
dir2 = v*(cos2+cos22)-vs*sin2+u*sin22; dir2 = dir2/norm(dir2); % direcction daugther 2

% The new branches are:    
Tubs = [Tubs;
    1 1 d1 l1 0 0 0 (v*cos1+vs*sin1)*l1 2 v acosd(cos1);
    1 1 d2 l2 0 0 0 (v*cos2-vs*sin2)*l2 2 v acosd(cos2)];

% RIGHT LUNG --------------------------------------------------------------

% Each colum of the cell array contains the referenced lobuls of each
% partition in each iteration. For example, in the first iteration will
% compute the branches that irrigate volume 1:3 and 4:11. The following
% iteration starts from the branch that irrigates the volume defined in the
% first row, in the case mentioned, it is 1:3 and is devided into 1:2 and 3
%sec = { 1:3, 1:2, 1,  4:5, 4, 6:10, 7:10,  7:8, 7, 9;
%       4:11,   3, 2, 6:11, 5,   11,    6, 9:10, 8, 10}; 
sec = { 1:3, 1:2,  1, 7:10,  4,  7:8,  7,  9;
       4:11,   3,  2,    6,  5, 9:10,  8, 10;
         [],  [], [],  4:5, [],   [], [], [];
         [],  [], [],   11, [],   [], [], []};

for k = 1:size(sec,2)
    Xn = Xarxa;
    Xn(ismember(SE,sec{1,k}),4) = 3; % Volume irrigated by daugther branch 1 is referenced with 3
    Xn(ismember(SE,sec{2,k}),4) = 4; % Volume irrigated by daugther branch 2 is referenced with 4

    g1 = Xn(:,4) == 3; 
    g2 = Xn(:,4) == 4; 

    % Trifurcacio (amb centre de massa)
    if k == 4
        Xn(ismember(SE,sec{3,k}),4) = 5; 
        Xn(ismember(SE,sec{4,k}),4) = 6; 

        g3 = Xn(:,4) == 5; g4 = Xn(:,4) == 6; 
        fill3 = 5;         fill4 = 6;
        in = 5;  % Row of the branch leading the trifurcation 4:5 - 11 - 6:10
        xf = Tubs(in,8:10);
        d1 = Tubs(in,3); % Mother's diameter

        cm1 = [0 0 0] + mean(Xn(Xn(:,4) == 3,1:3)); dir1 = (cm1-xf)/norm(cm1-xf);
        cm2 = [0 0 0] + mean(Xn(Xn(:,4) == 4,1:3)); dir2 = (cm2-xf)/norm(cm2-xf);
        cm3 = [0 0 0] + mean(Xn(Xn(:,4) == 5,1:3)); dir3 = (cm3-xf)/norm(cm3-xf);
        cm4 = [0 0 0] + mean(Xn(Xn(:,4) == 6,1:3)); dir4 = (cm4-xf)/norm(cm4-xf);
        v  = Tubs(in,8:10) - Tubs(in,5:7);    v = v/norm(v); % Unitary mother branch vector

        % Volums for each branch = sum number of points
        N1 = sum(Xn(:,4)==fill1);
        N2 = sum(Xn(:,4)==fill2);
        N3 = sum(Xn(:,4)==fill3);
        N4 = sum(Xn(:,4)==fill4);
        N0 = N1 + N2 + N3 + N4;

        % Diametres
        d2 = d1*(N1/N0)^(1/3);
        d3 = d1*(N2/N0)^(1/3);
        d4 = d1*(N3/N0)^(1/3);
        d5 = d1*(N4/N0)^(1/3);
        % Longituds
        l1 = M*d2;
        l2 = M*d3;
        l3 = M*d4;
        l4 = M*d5;
        Tubs = [Tubs;
            3 4 d2 l1 xf xf+dir1*l1 2 v acosd(dot(v,dir1));    
            3 4 d3 l2 xf xf+dir2*l2 2 v acosd(dot(v,dir2));    
            3 4 d4 l3 xf xf+dir3*l3 2 v acosd(dot(v,dir3));
            3 4 d5 l4 xf xf+dir4*l4 2 v acosd(dot(v,dir4))];
    else
        % Bifurcacio
        % Separation plane is find by Linear Discriminant Analisis method
        X = [Xn(g1,1:3); Xn(g2,1:3)];
        seg = cell(length(X),1);
        seg(1:sum(g1)) = {'s1'}; % Points corresponding to class 1 (volume 1)
        seg((1:sum(g2))+ sum(g1)) = {'s2'}; % Points corresponding to class 2 (volume 2)
        MdlLinear = fitcdiscr(X,seg); % Adjunt linear discriminant model
        L = MdlLinear.Coeffs(1,2).Linear; % Separation plane perpendicular vector
        vs = L'./norm(L); % Unitary normal separation plane vector
        % For the angles function:
        fill1 = 3;
        fill2 = 4;
        X     = Xn;
        if ismember(k, [6 8]) % In this columns of sec there is a change in the bronchial path
            switch k
%                 case 4
%                     in = 5; % Row of the branch leading to 4:11
%                     d1 = Tubs(in,3); % Mother's diameter
                case 6
                    in = 10; % Row of the branch leading to 6:10
                    d1 = Tubs(in,3);
                case 8
                    in = 17; % Row of the branch leading to 9:10
                    d1 = Tubs(in,3);
            end
            xf = Tubs(in,8:10);
            v  = Tubs(in,8:10) - Tubs(in,5:7);    v = v/norm(v); % Unitary mother branch vector
        else
            xf = Tubs(end-1,8:10);
            v  = Tubs(end-1,8:10) - Tubs(end-1,5:7);    v = v/norm(v); % Unitary mother branch vector
        end

        % Rule 1: The dot product between vs and cmp1 - xf has to be positive.
        % If not, vs = -vs.
        vs = sign(dot(vs,[0 0 0] + mean(X(X(:,4)==fill1,1:3))-xf))*vs;

        u     = cross(vs,v);                           u = u./norm(u); % Unitary normal bifurcation plane vector
        % Since separation vector does not necessarly be orhotgonal to mother
        % and bifrucation vectors, the cross product of these two is performed
        % to obtain the perpedindicular vector:
        us = cross(v,u); % Unitary normal separation vector for the axis
        [cos1,sin1,cos2,sin2,sin11,cos11,sin22,cos22,s1,s2] = angles(fill1,fill2,X,xf,n,u,v,us);
        dir1 = v*(cos1+cos11)+us*sin1+s1*u*sin11; dir1 = dir1/norm(dir1); % Direction daugther 1
        dir2 = v*(cos2+cos22)-us*sin2+s2*u*sin22; dir2 = dir2/norm(dir2); % Direction daugther 2
        if ismember(k,1:3)
            g = 1+k;
            m = 2*k;
        elseif k == 4
            g = 3;
            m = 5;
        elseif ismember(k,5:6)
            g = 4;
            m = k + 5;
        elseif k == 7
            g = 5;
            m = 14;
        elseif ismember(k,8:10)
            m = k + 9;
            if k == 8
                g = 6;
            else
                g = 7;
            end
        end

        % Volums for each branch = sum number of points
        N1 = sum(Xn(:,4)==3);
        N2 = sum(Xn(:,4)==4);
        N0 = N1 + N2;
        % Diametres
        d3 = d1*(N1/N0)^(1/3);
        d4 = d1*(N2/N0)^(1/3);
        % Longituds
        l1 = M*d3;
        l2 = M*d4;

        d1 = d3; % For the next iteration

        Tubs = [Tubs;
            g m d3 l1 xf xf+dir1*l1 2 v acosd(dot(v,dir1));
            g m d4 l2 xf xf+dir2*l2 2 v acosd(dot(v,dir2))];
    end
end

% Change natality of the branches irrigating the lobes:
Tubs([7:9 11 13:15 18:21],11) = 0; 

% LEFT LUNG ---------------------------------------------------------------
% Each colum of the cell array contains the referenced lobuls of each
% partition in each iteration. For example, in the first iteration will
% compute the branches that irrigate volume 1:3 and 4:11. The following
% iteration starts from the branch that irrigates the volume defined in the
% first row, in the case mentioned, it is 1:3 and is devided into 1:2 and 3
sec = {12:13, 12, 14:15, 14, 17:20, 17:18, 17, 19;
       14:20, 13, 16:20, 15, 16,    19:20, 18, 20}; 
st = length(Tubs(:,1));
for k = 1:size(sec,2)
    Xn = Xarxa;
    Xn(ismember(SE,sec{1,k}),4) = 3; % Volume irrigated by daugther branch 1 is referenced with 3
    Xn(ismember(SE,sec{2,k}),4) = 4; % Volume irrigated by daugther branch 2 is referenced with 4

    g1 = Xn(:,4) == 3; % Logical indexes for volume 1
    g2 = Xn(:,4) == 4; % Logical indexes for volume 2

    % Separation plane is find by Linear Discriminant Analisis method
    X = [Xn(g1,1:3); Xn(g2,1:3)];
    seg = cell(length(X),1);
    seg(1:sum(g1)) = {'s1'}; % Points corresponding to class 1 (volume 1)
    seg((1:sum(g2))+ sum(g1)) = {'s2'}; % Points corresponding to class 2 (volume 2)
    MdlLinear = fitcdiscr(X,seg); % Adjunt linear discriminant model
    L = MdlLinear.Coeffs(1,2).Linear; % Separation plane perpendicular vector
    vs = L'./norm(L); % Unitary normal separation plane vector   

    % For the angles function:
    fill1 = 3;
    fill2 = 4;
    X     = Xn;
    if ismember(k, [1,3,5,8]) % In this columns of sec there is a change in the bronchial path
        switch k
            case 1 
                in = 3; % Row of the branch leading to left lung
                d1 = Tubs(in,3); % Mother's diameter
            case 3
                in = 23; % Row of the branch leading to 14:20
                d1 = Tubs(in,3); % Mother's diameter
            case 5
                in = 27; % Row of the branch leading to 16:20
                d1 = Tubs(in,3);
            case 8
                in = 33; % Row of the branch leading to 19:20
                d1 = Tubs(in,3);
        end
        xf = Tubs(in,8:10);
        v  = Tubs(in,8:10) - Tubs(in,5:7);          v = v/norm(v); % Unitary mother branch vector 
    else
        xf = Tubs(end-1,8:10);
        v  = Tubs(end-1,8:10) - Tubs(end-1,5:7);    v = v/norm(v); % Unitary mother branch vector
    end
    if dot(vs,[0 0 0] + mean(X(X(:,4)==fill1,1:3))-xf) < 0
        vs = - vs;
    end
    u     = cross(vs,v);                            u = u./norm(u); % Unitary normal bifurcation plane vector
    % Since separation vector does not necessarly be orhotgonal to mother 
    % and bifrucation vectors, the cross product of these two is performed
    % to obtain the perpedindicular vector:
    us = cross(v,u); % Unitary normal separation vector for the axis
    
    [cos1,sin1,cos2,sin2,sin11,cos11,sin22,cos22,s1,s2] = angles(fill1,fill2,X,xf,n,u,v,us);
    
    dir1 = v*(cos1+cos11)+us*sin1+s1*u*sin11; dir1 = dir1/norm(dir1); % Direction daugther 1
    dir2 = v*(cos2+cos22)-us*sin2+s2*u*sin22; dir2 = dir2/norm(dir2); % Direction daugther 2
    
    % Generation and which is the mother
    if k == 1
        g = 2;
        m = 3;
    elseif k == 2
        g = 3;
        m = 24;
    elseif k == 3 
        g = k;
        m = 25;
    elseif k == 4 || k == 5
        g = 4;
        m = k + 24;
    elseif k == 6
        g = 5;
        m = 33;
    elseif k == 7 || k == 8
        g = 6;
        m = k + 27;
    end

    % Volums for each branch = sum number of points
    N1 = sum(Xn(:,4)==3);
    N2 = sum(Xn(:,4)==4);
    N0 = N1 + N2;
    % Diametres
    d3 = d1*(N1/N0)^(1/3);
    d4 = d1*(N2/N0)^(1/3);    
    % Longituds
    l1 = M*d3;
    l2 = M*d4;

    d1 = d3; % For the next iteration
    
    Tubs = [Tubs;
        g m d3 l1 xf xf+dir1*l1 2 v acosd(dot(v,dir1));
        g m d4 l2 xf xf+dir2*l2 2 v acosd(dot(v,dir2))];
end

% Change natality of the branches irrigating the lobes:
Tubs([24,25,28,29,31,34,35,36,37],11) = 0; 


%% REPRESENTATION ----------------------------------------------------------
pl = 1:11; 
pr = 12:20;
colormap jet;
cmap = colormap;  % Colors for the plot

figure(1);
hold on; axis equal
for k = 1:length(pl)
    b1 = bar3(nan,nan);
    b.FaceColor = cmap(k*floor(length(cmap)/length(pl)),:);
end
for k = 1:length(pr)
    b2 = bar3(nan,nan);
    b.FaceColor = cmap(k*floor(length(cmap)/length(pr)),:);
end
%title('Bronchial tree','FontSize',15,'Interpreter','latex')

i = 1; sec1 = [3,1,2,6,11,5,4,7,8,9,10,12,13,14,15,16,17,18,19,20];

for k = 1:size(Tubs,1)
    % All branches faded and black, except those irrigating a lobe. These
    % will be opaque and have the color of the lobe which they irrigate.
    if Tubs(k,11) == 0
        Tubs(k,11) = 1; 
        if k < 21
        plcolor = cmap(sec1(i)*floor(length(cmap)/length(pl)),:);
         R = double(ismember(segments,sec1(i))); % Get the desired volume of the matrix segments
        elseif k == 21
            plcolor = cmap(sec1(i)*floor(length(cmap)/length(pl)),:);
             R = double(ismember(segments,sec1(i))); % Get the desired volume of the matrix segments
            i = 0;
        else
            plcolor = cmap(i*floor(length(cmap)/length(pr)),:);
            R = double(ismember(segments,sec1(i+11))); % Get the desired volume of the matrix segments
        end
        branca(Tubs(k,5:7),Tubs(k,8:10), Tubs(k,3)/2, 100, plcolor, 1)
        fR = isosurface(R, 0.99);  % Create the patch object, isovalue = 0.99 since we want the contour
        fR.faces = fliplr(fR.faces); % Ensure normals point OUT
        coordenades_carina = [253 279 171];
        dX = 0.310547; % mida x pixel en mm
        dY = 0.310547; % mida y pixel en mm
        dZ = 0.625; % mida z pixel en mm
        fR = transformar_coordenades(fR,coordenades_carina,[dX dY dZ]); % Escalate to real values
        %figure(1);
        patch(fR,'FaceColor',plcolor,'FaceAlpha',0.1,'EdgeColor','none') % Represent patch object
        %display(sec1(i),'segment')
        i = i + 1;
        
    else
        branca(Tubs(k,5:7),Tubs(k,8:10), Tubs(k,3)/2, 100, 'k', .25);
    end

end
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

Tubs(7,16) = 3;  Tubs(8,16) = 1;  Tubs(9,16) = 2;   Tubs(18,16) = 7;
Tubs(19,16) = 8; Tubs(20,16) = 9; Tubs(21,16) = 10; Tubs(11,16) = 6;
Tubs(14,16) = 4; Tubs(15,16) = 5; Tubs(13,16) = 11;

Tubs([24,25,28,29,31,34,35,36,37],16) = [12:20];

save('BT_model1_2.mat','Tubs','-mat')

% leg1=legend(p,num2str(round(R(ind)')),'Location','NorthEast');set(leg1,'FontSize',9);
% ah1=axes('position',get(gca,'position'),'visible','off');
% leg2=legend(ah1,p1,magnitude{ind},'Location','NorthWest');set(leg2,'FontSize',9);
% legendTitle (leg1, 'R (km)' );
% legendTitle (leg2, 'Magnitude' );
% 
% l = legend('3','1','2','4','5','11','6','7','8','9','10','12','13','14','15','16','17','18','19','20');
% l.Location   = 'eastoutside';
% l.NumColumns = 2;
% axis equal




% Function to represent the volumes ---------------------------------------
function fv = transformar_coordenades(fv,xc,dx)
fv.vertices = (fv.vertices-xc).*dx;
end