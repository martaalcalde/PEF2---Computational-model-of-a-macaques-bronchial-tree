close all
clear all

load('Segments_Macacs_SenseFiltres.mat','segments')
% parametres de les dades llegides
coordenades_carina = [253 279 171]; %carina = [252 280 163];
dX = 0.310547; % mida x pixel en mm
dY = 0.310547; % mida y pixel en mm
dZ = 0.625; % mida z pixel en mm
% Creem una superficie per cadascun dels pulmons
SR = 1:11;
SL = 12:20;
RL = double(ismember(segments,SR));
fRL = isosurface(RL, 0.99);  % Create the patch object, isovalue = 0.99 since we want the contour
fRL.faces = fliplr(fRL.faces); % Ensure normals point OUT
LL = double(ismember(segments,SL));
fLL = isosurface(LL, 0.99);  % Create the patch object
fLL.faces = fliplr(fLL.faces);% Ensure normals point OUT
% Transformem les coordenades
fRL = transformar_coordenades(fRL,coordenades_carina,[dX dY dZ]);
fLL = transformar_coordenades(fLL,coordenades_carina,[dX dY dZ]);
% Especifiquem la distancia de la Xarxa en mm
di = 0.2; %fem exemple amb di gran i despres ja ho farem amb un petit
x  = (min(min(fRL.vertices(:,1)),min(fLL.vertices(:,1)))-dX):di:(max(max(fRL.vertices(:,1)),max(fLL.vertices(:,1)))+dX);
y  = (min(min(fRL.vertices(:,2)),min(fLL.vertices(:,2)))-dX):di:(max(max(fRL.vertices(:,2)),max(fLL.vertices(:,2)))+dY);
z  = (min(min(fRL.vertices(:,3)),min(fLL.vertices(:,3)))-dX):di:(max(max(fRL.vertices(:,3)),max(fLL.vertices(:,3)))+dZ);
[x,y,z] = meshgrid(x,y,z);
x  = x(:);
y  = y(:);
z  = z(:);
Xarxa = [x y z 0.*x 0.*x];
%clear x y z LL RL
% Ens falta assignar cada punt de la xarxa a quina branca (pulmo) pertany
iR = inpolyhedron(fRL,Xarxa(:,1:3));
iL = inpolyhedron(fLL,Xarxa(:,1:3));
iB = iR & iL;
iR(iB) = 0; % Els punts que pertanyen als dos pulmons els fiquem com si no pertanyesin a cap
iL(iB) = 0;
clear iB
Xarxa(iR,4) = 1;
Xarxa(iL,4) = 2;
Xarxa(Xarxa(:,4)==0,:) = []; % Eliminem els punts de la xarxa que no formen part de cap pulmó
iR = Xarxa(:,4)==1;
iL = Xarxa(:,4)==2;
XX = []; se = [];
for k = SR
    SS = segments == k;
    fS = isosurface(SS, 0.99);  % Create the patch object
    fS.faces = fliplr(fS.faces);% Ensure normals point OUT
    fS = transformar_coordenades(fS,coordenades_carina,[dX dY dZ]);
    iS = inpolyhedron(fS,Xarxa(:,1:3));
    i0 = Xarxa(:,5)==0;
    Xarxa(iR & iS &  i0, 5) = k;
    Xarxa(iR & iS & ~i0, 5) = nan;
    XX = [XX; fS.vertices];
    se = [se; k*ones(size(fS.vertices,1),1)];
end
id = find((isnan(Xarxa(:,5)) | Xarxa(:,5)==0) & iR);
P  = Xarxa(id,1:3); % Mirem a quin segment ficar els indeterminats
for k = 1:size(P,1)
    d = sum((XX-P(k,:)).^2,2);
    [~,i] = min(d);
    Xarxa(id(k),5)  = se(i);
end
XX = []; se = [];
for k = SL
    SS = segments == k;
    fS = isosurface(SS, 0.99);  % Create the patch object
    fS.faces = fliplr(fS.faces);% Ensure normals point OUT
    fS = transformar_coordenades(fS,coordenades_carina,[dX dY dZ]);
    iS = inpolyhedron(fS,Xarxa(:,1:3));
    i0 = Xarxa(:,5)==0;
    Xarxa(iL & iS &  i0,5) = k;
    Xarxa(iL & iS & ~i0,5) = nan;
    XX = [XX; fS.vertices];
    se = [se; k*ones(size(fS.vertices,1),1)];
end
id = find((isnan(Xarxa(:,5)) | Xarxa(:,5)==0) & iL);
P  = Xarxa(id,1:3);
for k = 1:size(P,1)
    d = sum((XX-P(k,:)).^2,2);
    [~,i] = min(d);
    Xarxa(id(k),5)  = se(i);
end

save('xarxa02.mat','Xarxa','-mat')

%%
% Angle, punt_pla i branca els podem posar en fitxers separats
function fv = transformar_coordenades(fv,xc,dx)
fv.vertices = (fv.vertices-xc).*dx;
end
function [cos1,sin1,cos2,sin2,cmp1,cmp2] = Angle(fill1,fill2,X,xf,n,u,v)

% Calculem els cabals de cadascuna de les branques
q1 = sum(X(:,4)==fill1);
q2 = sum(X(:,4)==fill2);
qT = q1 + q2;
q1 = q1/qT;
q2 = q2/qT;

% valors de alpha (minmitzacio del treball)
cos1=(1+q1^(4/n)-(1-q1)^(4/n))/(2*q1^(2/n));
cos2=(1+q2^(4/n)-(1-q2)^(4/n))/(2*q2^(2/n));

% calculem el centre de masses
cm1 = [0 0 0] + mean(X(X(:,4)==fill1,1:3));
cm2 = [0 0 0] + mean(X(X(:,4)==fill2,1:3));
% calculem la projeccio del centre de masses sobre el pla de bifurcacio
A=u(1); B=u(2); C=u(3); D=-(A*xf(1)+B*xf(2)+C*xf(3));
if A*cm1(1)+B*cm1(2)+C*cm1(3)+D==0
    cmp1=cm1;
else
    cmp1=punt_pla(cm1,u,xf);
end
if A*cm2(1)+B*cm2(2)+C*cm2(3)+D==0
    cmp2=cm2;
else
    cmp2=punt_pla(cm2,u,xf);
end

% calculem l'angle mitjana
% branca fill1
angle_alpha = acos(cos1);
vect = cmp1-xf;
angle_beta = acos((vect*v')/(norm(vect)*norm(v)));
angle = min(mean([angle_alpha,angle_beta]),pi/2);
cos1 = cos(angle);
sin1 = sin(angle);
% branca fill2
angle_alpha = acos(cos2);
vect = cmp2-xf;
angle_beta = acos((vect*v')/(norm(vect)*norm(v)));
angle = min(0.5*angle_alpha+0.5*angle_beta,pi/2);
cos2 = cos(angle);
sin2 = sin(angle);

end
function [ projeccio ] = punt_pla( punt,normal,puntdelpla )

A=normal(1); B=normal(2); C=normal(3);
D=-(A*puntdelpla(1)+B*puntdelpla(2)+C*puntdelpla(3));

t=-(A*punt(1)+B*punt(2)+C*punt(3)+D)/(A*A+B*B+C*C);

projeccio(1)=punt(1)+A*t;
projeccio(2)=punt(2)+B*t;
projeccio(3)=punt(3)+C*t;

end
function [  ] = branca( inici, final, radi, punts, color, transparencia )

normal = inici-final;
longitud = norm(normal);
vx = normal(1)/longitud;
vy = normal(2)/longitud;
vz = normal(3)/longitud;
if vz<0
    vx=-vx;
    vy=-vy;
    vz=-vz;
end

[X, Y, Z] = cylinder(radi,punts);
Z=longitud*Z;
if vy==0 && vz==0
    phi=0;
else
    phi = asin(-vy/sqrt(vy^2+vz^2));
end
psi = asin(vx);
chi = 0;
a = [cos(psi)*cos(chi), -cos(psi)*sin(chi), sin(psi); ...
    cos(phi)*sin(chi)+sin(phi)*sin(psi)*cos(chi), ...
    cos(phi)*cos(chi)-sin(phi)*sin(psi)*sin(chi), ...
    -sin(phi)*cos(psi); ...
    sin(phi)*sin(chi)-cos(phi)*sin(psi)*cos(chi), ...
    sin(phi)*cos(chi)+cos(phi)*sin(psi)*sin(chi), ...
    cos(phi)*cos(psi)];

Xr = a(1,1)*X+a(1,2)*Y+a(1,3)*Z+inici(1);
Yr = a(2,1)*X+a(2,2)*Y+a(2,3)*Z+inici(2);
Zr = a(3,1)*X+a(3,2)*Y+a(3,3)*Z+inici(3);
MX = max(inici(1),final(1))+radi;
mX = min(inici(1),final(1))-radi;
XR = mean(Xr(:));
MY = max(inici(2),final(2))+radi;
mY = min(inici(2),final(2))-radi;
YR = mean(Yr(:));
MZ = max(inici(3),final(3))+radi;
mZ = min(inici(3),final(3))-radi;
ZR = mean(Zr(:));
condicio = XR>mX & XR<MX & YR>mY & YR<MY & ZR>mZ & ZR<MZ;
if ~condicio
    Xr = a(1,1)*X+a(1,2)*Y+a(1,3)*Z+final(1);
    Yr = a(2,1)*X+a(2,2)*Y+a(2,3)*Z+final(2);
    Zr = a(3,1)*X+a(3,2)*Y+a(3,3)*Z+final(3);
end

h1 = surf(Xr,Yr,Zr);
hold on
set(h1,'FaceColor',color);
set(h1,'EdgeColor','None');
color=get(h1,'FaceColor');
color_fosc=color/3;
h2=plot3(Xr(1,:),Yr(1,:),Zr(1,:),'linewidth',1);
set(h2,'Color',color_fosc);
h3=plot3(Xr(2,:),Yr(2,:),Zr(2,:),'linewidth',1);
set(h3,'Color',color_fosc);
alpha(h1,transparencia)
end