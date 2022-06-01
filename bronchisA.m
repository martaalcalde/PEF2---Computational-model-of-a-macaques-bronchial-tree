clear all
%load BT_2CM3CM4CMdTraq81.mat % Tubs matrix contiaing the branches leading to the lobe    
load BT_model1_2.mat;
load xarxa.mat          % matrix Xarxa

% To generate the bronchis at each bifrucation it is followed the same
% process. 
% The direction of each bronchi is to follow the mass center of
% the irrigated volume.
% To distinguish terminal / already bifrucated / to be bifurcated branches
% the following totation for the natality column is used:
% Natality:
%   - 0: terminal branch
%   - 1: can bifrucate but has not done it yet
%   - 2: has been bifurcated
% There is an exception for those branches with only one daughter, this are
% stored and removed from the set of branches that will bifrucate next.
% Each time a new branch emerges, the volume that it irrigates is saved in
% the fourth column of the matrix Xn with the integer of its row in Tubs.

index = find(Tubs(:,11) == 1); % Branches that have to be bifurcated
oneD = []; % Matrix of branches that have only one daughter
Xn = Xarxa(:,1:4);
Xn(:,4) = Xarxa(:,5); % First iteration computs first lobe bifurcation
M = 3; n = 3;

tic
count = 1; % to keep track of the number of iterations
while ~isempty(index) % while there are branches to be bifurcated
    for k = 1:length(index)
        
        j = index(k); % Is the row of the branch being bifurcated, but also
        % it is the integer defining the volume that the branch irrigates

        xf = Tubs(j,8:10);                                 % final point
        d0 = Tubs(j,3);                                    % diameter
        v = Tubs(j,8:10) - Tubs(j,5:7);    v = v/norm(v);  % mother vector
        va = Tubs(j,12:14);                                % 'grandmother' vector
        vs = cross(v,va);    vs = vs/norm(vs);             % separation vector
        
        if count == 1
            vol = Xn(:,4) == Tubs(j,end);
        else 
            vol = Xn(:,4) == j; % logical indexes for the volume being irrigated
        end
        Xnn = Xn(vol,:);    % reduced matrix to lower the dimensions of the 
                            % calculations, thus the computational cost

        % Find the two subvolumes
        vS = repmat(vs,[sum(vol),1]);
        class = dot([Xnn(:,1)-xf(1), Xnn(:,2)-xf(2), Xnn(:,3)-xf(3)]', vS'); 
       
        % Make sure that both sub-volumes exist, if one of them is 0 it
        % means that that part is already irrigated so there is only
        % one duaghter branch. In that case the index of the mother branch
        % is save so it does not bifurcate again, even if its natality is
        % 1.
        right = class > 0;
        N1 = sum(right);
        N2 = sum(~right);

        if N1 == 0 || N2 == 0 % One volume is 0
            diam = d0;      % Diameter of the daughter = mother's diameter
            long = d0*M;    % Length of the daughter = mother's length
            
            % Update the integer defining the volume
            t = size(Tubs,1);         
            if N1 == 0
                Xnn(~right,4)  = t+1;
            else
                Xnn(right,4) = t+1;
            end

            % The direction changes:
            cm = [0 0 0] + mean(Xn(Xn(:,4)==j,1:3));
            dir = (cm-xf)/norm(cm-xf);

            % Save the integers in the whole matrix
            Xn(vol,4) = Xnn(:,4);

            Tubs = [Tubs;
                Tubs(j,1)+1 j diam long xf xf+dir*long 1 v acosd(dot(v,dir)) Tubs(j,16)];

            % Natality of mother branch stays 1, however this branch has to
            % be disntiguish from the non bifurcated yet, so the index j
            % will be stored in this matrix of one-daughter branches:
            oneD = [oneD j];

        else
            % The number asigned to the segment irrigated for each daughter branch
            % will be the row of the branch in the matrix Tubs:
            t = size(Tubs,1);
            Xnn(right,4)  = t+1;
            Xnn(~right,4) = t+2;

            % Save the integers in the whole matrix
            Xn(vol,4) = Xnn(:,4);

            % Diameters and lengths of the daughter branches
            d1 = d0 * (N1/(N1+N2))^(1/3);    l1 = M*d1;
            d2 = d0 * (N2/(N1+N2))^(1/3);    l2 = M*d2;
            
            u = cross(v,vs);
            [cos1,sin1,cos2,sin2] = Angle(t+1,t+2,Xnn,xf,n,u,v);

            dir1 = v*cos1+vs*sin1; dir1 = dir1/norm(dir1); % direcction daugther 1
            dir2 = v*cos2-vs*sin2; dir2 = dir2/norm(dir2); % direcction daugther 2

            Tubs(j,11) = 2;  % As j has bifuracted, its natality changes

            % A bronchi is termial if the diameter is smaller than 0.5mm
            if d1 > 0.5 
                nat1 = 1;
            else
                nat1 = 0;
            end
            if d2 > 0.5 
                nat2 = 1;
            else
                nat2 = 0;
            end

            Tubs = [Tubs;
                Tubs(j,1)+1 j d1 l1 xf xf+dir1*l1 nat1 v acosd(dot(v,dir1)) Tubs(j,16);
                Tubs(j,1)+1 j d2 l2 xf xf+dir2*l2 nat2 v acosd(dot(v,dir2)) Tubs(j,16)];

        end

    end
    index = find(Tubs(:,11) == 1);
    index(ismember(index,oneD)) = []; % Take out of the index list the branches that have only one daugther
    count = count + 1;
end
toc

save('bronchis_model1_2.mat','Tubs','-mat')