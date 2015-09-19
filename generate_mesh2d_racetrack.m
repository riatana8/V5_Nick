%Racetrack with prescribed peristalsis and pericardium
close all
clear all

L = 1;                              % length of computational domain (m)
N = 1024;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = L/(2*N);                       % space between boundary points in straight tube
ds2 = L/(4*N);                       % space between boundary points in straight tube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the racetrack

Let = 1.34*0.1;                          % Length of elastic tube (m)
Nend = 10;                           % Number of rigid points on each end of elastic section
Lt = Let+2*Nend*ds2;                 % Length of straight section with three rigid points on each end

diameter = 0.1;                     % diameter of the tube
R2 = 0.1;                           % radius of inner wall
R1 = R2+diameter;                   % radius of outer wall

Nstraight = 2*ceil(Lt/ds)          % number of points along each straight section
Ncurve = 2*ceil(pi*R1/ds);          % number of points along each curved section
Nrace = Nstraight+2*Ncurve;         % number of points making up the racetrack part
dtheta = pi/(Ncurve/2);             % angle increment for drawing curved edges

mesh_name = 'heart_';               % structure name

centery = 0;                        % y-position of center of curved sections
centerx1 = -0.5*Lt;                 % x-position of center of left curved section
centerx2 = 0.5*Lt;                  % x-position of center of right curved section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the pericardium
Dp = 2*diameter;                    %diameter of the pericardium
Nperi = 2*ceil((Dp-diameter)/ds);  % number of boundary points along the sides of the pericardium
Nperitot = Nperi + Nstraight;       % total number of pericardium points


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for markers
Nmarkersx = 11;                     %number of columns of markers
Nmarkersy = 11;                     %number of markers in each column
Nmarkers=Nmarkersx*Nmarkersy;       %total number of markers
dmx = Let/(Nmarkersx-1);            %space between markers in x-direction
dmy = diameter/(Nmarkersy-1);       %space between markers in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material parameters
kappa_spring = 1.0;               % spring constant (Newton)
kappa_beam = 0.1;                 % beam stiffness constant (Newton m^2) %2.5e-2 works for Wo>=5
kappa_target = kappa_spring;        % target point penalty spring constant (Newton)
Fmag = 4.0e0;                % this is my best guess at a reasonable applied force %4.0e0 works for Wo>=5
phase = 0;                      %initial phase of the oscillating force, where F=Fmag*phase and phase = (1+sin(2*pi*f*t-pi/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the elastic section of the tube
% Write out the vertex information

vertex_fid = fopen([mesh_name 'tube_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', 2*Nstraight);
figure(1)
hold on
%top part
for i=1:Nstraight,
    ytop = centery-R2;
    xtop = -Lt/2+(i-1)*ds2;
    plot(xtop,ytop,'k')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

%bottom part
for i=Nstraight+1:2*Nstraight,
    ybot = centery-R1;
    xbot = -Lt/2+(i-ceil(Nstraight/2)-1)*ds2;
    plot(xtop,ytop,'k')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end
fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make markers as vertices with no material properties

vertex_fid = fopen(['markers_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nmarkers);

%top part
for i=0:Nmarkersx-1,
    for j=0:Nmarkersy-1,
        y = centery-R2-j*dmy;
        x = -Let/2+i*dmx;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% race track part
% Write out the vertex information

vertex_fid = fopen([mesh_name 'race_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nrace);

%right curved part of racetrack
for i=1:ceil(Ncurve/2),
    theta = (i-1)*dtheta-pi/2;
    yin = centery+R2*sin(theta);
    xin = Lt/2+R2*cos(theta);
    plot(xin,yin,'b')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
end

for i=ceil(Ncurve/2)+1:Ncurve,
    theta=(i-Ncurve/2-1)*dtheta-pi/2;
    yout = centery+R1*sin(theta);
    xout = Lt/2+R1*cos(theta);
    plot(xout,yout,'b')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
end

%straight section on the top
for i = Ncurve+1:Ncurve+ceil(Nstraight/2),
    yin = centery+R2;
    xin = centerx2-(i-Ncurve-1)*ds;
    plot(xin,yin,'b')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
end

for i = Ncurve+ceil(Nstraight/2)+1:Ncurve+Nstraight,
    yout = centery+R1;
    xout = centerx2-(i-Ncurve-ceil(Nstraight/2)-1)*ds;
    plot(xout,yout,'b')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
end

%left curved part of racetrack
for i = Ncurve+Nstraight+1:Ncurve+Nstraight+ceil(Ncurve/2),
    theta = pi/2+(i-Ncurve-Nstraight-1)*dtheta;
    yin = centery+R2*sin(theta);
    xin = centerx1+R2*cos(theta);
    plot(xin,yin,'b')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
end

for i = Ncurve+Nstraight+ceil(Ncurve/2)+1:2*Ncurve+Nstraight,
    theta = pi/2+(i-Ncurve-Nstraight-ceil(Ncurve/2)-1)*dtheta;
    yout = centery+R1*sin(theta);
    xout = centerx1+R1*cos(theta);
    plot(xout,yout,'b')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
end
fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pericardium
% Write out the vertex information
%
vertex_fid = fopen([mesh_name 'peri_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nperitot);

% make the top and bottom of the pericardium
for i=1:ceil(Nstraight/2),
    ytop = centery-(R2-(Dp-diameter)/2);
    xtop = -Lt/2+i*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
    plot(xtop,ytop)
end

for i=ceil(Nstraight/2)+1:Nstraight,
    ybot = centery-R1-(Dp-diameter)/2;
    xbot = Lt/2-(i-ceil(Nstraight/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
    plot(xbot,ybot)
end

% make the four side pieces
for i=Nstraight+1:Nstraight+ceil(Nperi/4),
    y = centery-(R1+(Dp-diameter)/2)+(i-Nstraight-1)*ds;
    x = -Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y,'r')
end

for i=Nstraight+ceil(Nperi/4)+1:Nstraight+ceil(Nperi/2),
    y = centery-R2+(i-Nstraight-ceil(Nperi/4)-1)*ds;
    x = -Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y,'b')
end

for i=Nstraight+ceil(Nperi/2)+1:Nstraight+ceil(3*Nperi/4),
    y = centery-(R1+(Dp-diameter)/2)+(i-Nstraight-ceil(Nperi/2)-1)*ds;
    x = Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y,'g')
end

for i=Nstraight+ceil(3*Nperi/4)+1:Nperitot,
    y = centery-R2+(i-Nstraight-ceil(3*Nperi/4)-1)*ds;
    x = Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y,'k')
end
fclose(vertex_fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the spring information for the elastic section

spring_fid = fopen([mesh_name 'tube_' num2str(N) '.spring'], 'w');
fprintf(spring_fid, '%d\n', Nstraight-2);

%elastic part of tube
for i = 0:ceil(Nstraight/2)-2,
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds/(ds^2), ds);
end

for i = ceil(Nstraight/2):Nstraight-2,
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds/(ds^2), ds);
end

fclose(spring_fid);
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Write out the beam information for the elastic section

beam_fid = fopen([mesh_name 'tube_' num2str(N) '.beam'], 'w');
fprintf(beam_fid, '%d\n', Nstraight-4);

%elastic part of tube
for i = 0:ceil(Nstraight/2)-3,
    fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds/(ds^4));
end

for i = ceil(Nstraight/2):Nstraight-3,
    fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds/(ds^4));
end
fclose(beam_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the target point information for the ends of the elastic tube
target_fid = fopen([mesh_name 'tube_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', 2*Nstraight); 


2*Nstraight
for i = 0:(Nstraight)-1,   
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds2/(ds2^2));
end

for i = (Nstraight):(2*Nstraight)-1,   
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds2/(ds2^2));
end


% 
% for i = 0:Nend-1,   %OLD
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% for i = ceil(Nstraight/2)-Nend:ceil(Nstraight/2)-1,  %OLD
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% for i = ceil(Nstraight/2):ceil(Nstraight/2)+Nend-1,   %OLD
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% for i = Nstraight-Nend:Nstraight-1,    %OLD
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end

fclose(target_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the racetrack
target_fid = fopen([mesh_name 'race_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nrace);

for i = 0:Nrace-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write out the target point information for the pericardium
target_fid = fopen([mesh_name 'peri_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nperitot);

for i = 0:Nperitot-1,
   fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integrating nick's code for pumping

 HeartTube(diameter, Lt, L, R1,R2, ds2, centery, kappa_spring, kappa_beam, kappa_target)

