%% V5.2

%% 3/29 TS notes: 

close all
clear all
close all


%% In update_target_point_positions, changed R2 to 0.177 (From graph).
% Same in update_target_point_positions_peri


%% 
L = 1;                              % length of computational domain (m)
N = 1024;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = L/(2*N);                       % space between boundary points in straight tube

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the racetrack

Let = 1.34*0.1;                          % Length of elastic tube (m) scaled wrt to diameter of heart
% V5 TS. Nick's Heart_Tube.m code uses a fraction of the length of the tube to determine
% the number of rigid points on each end of elastic section. Made changes
% to use the same number he's using (can output these numbers from his code). 
Nend = 8
% Nend = 10;                           % Number of rigid points on each end of elastic section
Lt = Let+2*Nend*ds;                 % Length of main tube straight section with Nend rigid points on each end



% main tube
diameter = 0.1;                     % diameter of the main tube (scaling factor for everything else)
R2 = 0.1;                           % radius of inner wall
R1 = R2+diameter;                   % radius of outer wall

%Uppermost branch
diamTop = diameter*(0.0206/0.0442);		%diameter of the uppermost branch
LtTop = Lt+2*(diameter-diamTop); %length of the straight parts for the upper section
R1Top = R2 + diamTop;		%radius of outer wall for upper section

%bottom branch
diamBot = diameter*(0.0158/0.0442);
R2Bot = 1/2*(diamBot);


Ls = 2*R2+3*diamBot;			%length of the straight outer side walls


Nstraight = 2*ceil(Lt/ds)          % number of points along each shorter straight section (top+bottom)
NstraightBranch = 2*ceil(LtTop/ds);	% number of points along longer straight sections
NstraightSide = 2*ceil(Ls/ds);	%number of points along each straight (vertical)  outer side section
Ncurve = 2*ceil(pi*R1/ds);          % number of points along a circle of diameter R1. 
                                       %Will use fractions of this number
                                       %to draw points for curved sections
%%NcurveSmall = 2*ceil(pi*R2Bot/ds)    %number of points along a circle of diameter R2Bot (small)


Nrace = (1/2)*Nstraight+4*Ncurve + (3/2)*NstraightBranch + NstraightSide; %+NcurveSmall % number of points making up the racetrack part. More added later in code

%For Nrace: (1/2)*Nstraight = one straight section of length Lt,
%Nstraight/2 each
    % 4*Ncurve = 6 semi circles and 4 quarter circles = 4 circles, Ncurve
    % each
    %(3/2)*NstraightBranch = 3 straight sections of length Ls,
    %NstraightBranch/2 each.
    % NstraightSide = 2 straight (vertical) sections, NstraightSide/2 each.
    %Ncurve small = points for 1 cirle of radius R2Bot
    

dthetaHalf = pi/(Ncurve/2);             % angle increment for drawing curved edges (for a semi circle)
dthetaQuart = (pi/2)/(Ncurve/4);             % angle increment for drawing curved edges (for a quarter of a circle)


mesh_name = 'heart_';               % structure name


centery = -0.1;                        % y-position of center of curved sections in the lower area
centeryMid = centery+R2+(3)*R2Bot; %y-position of the center of the curved sections in the middle area
centeryTop = centeryMid+3*R2Bot +R2;	%y-position of center of curved sections in the upper area
centerx1 = -0.5*Lt;                 % x-position of center of left curved section of main tube
centerx2 = 0.5*Lt;                  % x-position of center of right curved section of main tube
centerx1Top = -0.5*LtTop;   %shifting of upper branch center of left cuved section
centerx2Top = 0.5*LtTop ; %shifting of uppder branch center of right curved section

%Note: more parameters for middle section added later in the code. Includes
%adding more points to Nrace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the pericardium
Dp = 2*diameter;                    %diameter of the pericardium
Nperi = 2*ceil((Dp-diameter)/ds);  % number of boundary points along the sides of the pericardium
Nperitot = Nperi + Nstraight;       % total number of pericardium points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the actuator
La = 0.04;                          % length of the actuator section
NLa = ceil(La/ds);                  % number of points along each actuator
Ca = 0.25*Lt;                       % center of the actuator section
NCa = ceil(ceil(Nstraight/2)*Ca/Lt); % index of the center point
Na1 = NCa - ceil(NLa/2);            % index of the starting point with respect to the elastic section
Na2 = Na1+NLa-1;                    % index of the ending point with respect to the elastic section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for peristalsis
Lnperi = 0.005;                          % length of the end of the flexible tube without peristalsis (modified from 0.05)
Lperi = Lt-2*Lnperi;                      % length of the peristalsis section of tube
Nperist = ceil(Lperi/ds);                 % number of peristaltic points on top or bottom of tube
NCent = ceil(Nstraight/4);                 %index of center of elastic tube
Nperi1 = NCent-ceil(0.5*Nperist);           % index of the starting point with respect to the elastic section
Nperi2 = Nperi1+Nperist-1;                  %index of the ending point with respect to the elastic section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the prescribed peristalsis
Lap = Let*0.75;                          % length of the actuator section
NLap = ceil(Lap/ds);                  % number of points along each actuator
Cap = 0.5*Lt;                       % center of the actuator section
NCap = ceil(ceil(Nstraight/2)*Cap/Lt); % index of the center point
Na1p = NCap - ceil(NLap/2);            % index of the starting point with respect to the elastic section
Na2p = Na1p+NLap-1;                    % index of the ending point with respect to the elastic section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for markers
Nmarkersx = 11;                     %number of columns of markers
Nmarkersy = 11;                     %number of markers in each column
Nmarkers=Nmarkersx*Nmarkersy       %total number of markers
dmx = Let/(Nmarkersx-1);            %space between markers in x-direction
dmy = diameter/(Nmarkersy-1);       %space between markers in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% material parameters
kappa_spring = 1.0e-1;               % spring constant (Newton)
kappa_beam = 1.0e-1;                 % beam stiffness constant (Newton m^2)
kappa_target = 2.0e0*15;        % target point penalty spring constant (Newton)
Fmag = kappa_spring;                % this is my best guess at a reasonable applied force
phase = 0;                      %initial phase of the oscillating force, where F=Fmag*phase and phase = (1+sin(2*pi*f*t-pi/2));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the elastic section of the tube
% Write out the vertex information

vertex_fid = fopen([mesh_name 'tube_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nstraight);

figure
hold on
n = 0;
p = 0;

%bot part
for i=1:ceil(Nstraight/2),
    ybot = centery-R1;
    xbot = -Lt/2+(i-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
    plot(xbot,ybot,'b*')
     n = n+1; 
end

%bottom part
for i=ceil(Nstraight/2)+1:Nstraight,
    ytop = centery-R2;
    xtop = -Lt/2+(i-ceil(Nstraight/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
    plot(xtop,ytop,'k*')
     p = p+1;
end
p
n


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
    %plot(x,y)
    end
end
fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% actuator part
% Write out the vertex information

%top part
% vertex_fid = fopen(['actuator_top_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', NLa);
% 
% for i=Na1:Na2,
%     ytop = centery-R2;
%     xtop = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
%     plot(xtop,ytop);
% end
% fclose(vertex_fid);
% 
% %bottom part
% vertex_fid = fopen(['actuator_bot_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', NLa);
% 
% for i=Na1:Na2,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%      plot(xbot,ybot);
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% valveless pumping applied force
% Use either the target point actuator or the applied force, but not both.
% Write out the vertex information

% vertex_fid = fopen(['vp_aforce_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', 2*NLa);
% 
% %top points
% for i=Na1:Na2,
%     ytop = centery-R2;
%     xtop = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
%      plot(xtop,ytop);
% end
% 
% %bottom points
% for i=Na1:Na2,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%      plot(xbot,ybot);
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peristalsis applied force
% Use either the target point actuator or the applied force, but not both.
% Write out the vertex information

% vertex_fid = fopen(['peri_aforce_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', 2*Nperist);
% 
% %top points
% for i=Nperi1:Nperi2,
%     ytop = centery-R2;
%     xtop = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
%      plot(xtop,ytop);
% end
% 
% %bottom points
% for i=Nperi1:Nperi2,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%      plot(xbot,ybot);
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%
% prescribed peristalsis (motion) part
% Write out the vertex information
% 
% %top part
% vertex_fid = fopen(['pperi_top_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', NLap);
% 
% for i=Na1p:Na2p,
%     ytop = centery-R2;
%     xtop = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
%     plot(xtop,ytop,'r*')
% end
% fclose(vertex_fid);
% 
% %bottom part
% vertex_fid = fopen(['pperi_bot_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', NLap);
% 
% for i=Na1p:Na2p,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%     plot(xbot,ybot,'r*')
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% race track part
% %%%calculations necessary before writing out vertex info:
%outer quarter circle
for i=ceil(Ncurve/2)+1:ceil(3*Ncurve/4),
    theta=(i-Ncurve/2-1)*dthetaQuart-(pi/2);
    yout = centery+R1*sin(theta);
    xout = Lt/2+R1*cos(theta);
end

%starting points for the right outer straight section
yrStart = yout+ds;
xrStart = xout+ds;

%outer quarter circle
for i = 1:ceil(Ncurve/4)
    theta = (-pi/2)-(i)*dthetaQuart;
    yout = centery+R1*sin(theta);
    xout = centerx1+R1*cos(theta);
end
%starting points for the left outer stright section
ylStart = yout+ds;
xlStart = xout+ds;

%additional parameters for middle branches:
centerx2Mid = xrStart-diameter;
centerx1Mid = xlStart+diameter;
Lm = abs (centerx2Mid-centerx1Mid);%length of middle branches (to ensure width of one diameter on either side, vertically)
NstraightMid = 2*ceil(Lm/ds);

Nrace = Nrace+NstraightMid;


% Write out the vertex information

vertex_fid = fopen([mesh_name 'race_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nrace);

%%%right curved part of lower racetrack section

%inner semicircle
for i=1:ceil(Ncurve/2)
    theta = (i-1)*dthetaHalf- pi/2;
    yin = centery+R2*sin(theta);
    xin = Lt/2+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
    % plot(xin,yin)
end

%outer quarter circle
for i=ceil(Ncurve/2)+1:ceil(3*Ncurve/4),
    theta=(i-Ncurve/2-1)*dthetaQuart-(pi/2);
    yout = centery+R1*sin(theta);
    xout = Lt/2+R1*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
    %plot(xout, yout)
end

%starting points for the right outer straight section
yrStart = yout+ds;
xrStart = xout+ds;

%right outer straight section
for i = 1:ceil(NstraightSide/2)
    xRside = xrStart;
    yRside = yrStart +i*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xRside, yRside);
    %plot(xRside,yRside)
end



%straight section for the upper part of the lower section (bottom wall of
%bottom branch)
for i = Ncurve+1:Ncurve+ceil(Nstraight/2),
    yin = centery+R2;
    xin = centerx2-(i-Ncurve-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)

end


%%%left curved part of racetrack
%inner semicircle
for i = 1:ceil(Ncurve/2),
    theta = pi/2+i*dthetaHalf;
    yin = centery+R2*sin(theta);
    xin = centerx1+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
    %plot(xin, yin)
end

%outer quarter circle
for i = 1:ceil(Ncurve/4)
    theta = (-pi/2)-(i)*dthetaQuart;
    yout = centery+R1*sin(theta);
    xout = centerx1+R1*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
   % plot(xout, yout)
end

ylStart = yout+ds;
xlStart = xout+ds;

%left outer straight section
for i = 1:ceil(NstraightSide/2)
     xLside = xlStart;
     yLside = ylStart + (i)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xLside, yLside);
   % plot(xLside,yLside)
end

%%Middle Branches


%%top wall of lower skinny branch
for i = 1:ceil(NstraightMid/2)
    yin = centeryMid-R2Bot;
    xin = centerx2Mid-(i)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)
end

%%bottom wall of upper skinny branch
for i = 1:ceil(NstraightMid/2)
    yin = centeryMid+R2Bot;
    xin = centerx2Mid-(i)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)
end

%%Left inner semi circle

for i = 1:ceil(Ncurve/2),
    theta = pi/2+i*dthetaHalf;
    yin = centeryMid+R2Bot*sin(theta);
    xin = centerx1Mid+R2Bot*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin, yin)
end
%% right inner semi circle
for i=1:ceil(Ncurve/2)
    theta = (i-1)*dthetaHalf- pi/2;
    yin = centeryMid+R2Bot*sin(theta);
    xin = centerx2Mid+R2Bot*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)
end


%%%%%%%%%% upper branches

%right curved part of racetrack
for i=1:ceil(Ncurve/2)
    theta = (i-1)*dthetaHalf- (pi/2);
    yin = centeryTop+R2*sin(theta);
    xin = centerx2Top+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)
end


for i=1:ceil(Ncurve/4),
    theta=0+i*dthetaQuart;
    yout = centeryTop+R1Top*sin(theta);
    xout = centerx2Top+R1Top*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
   % plot(xout, yout)
end



%straight (horizontal) sections on the upper chamber 
for i = Ncurve+ceil(NstraightBranch/2)+1:Ncurve+NstraightBranch,
    yout = centeryTop-R2;
    xout = centerx2Top-(i-Ncurve-ceil(NstraightBranch/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
   % plot(xout, yout)
    
end
%make the walls of the upper
%branch
for i = Ncurve+1:Ncurve+ceil(NstraightBranch/2),
    yin = centeryTop+R2;
    xin = centerx2Top-(i-Ncurve-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)
    
end

for i = Ncurve+ceil(NstraightBranch/2)+1:Ncurve+NstraightBranch,
    yout = centeryTop+R1Top;
    xout = centerx2Top-(i-Ncurve-ceil(NstraightBranch/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
   % plot(xout, yout)
end



%%%left curved part of upper racetrack

%inner semicircle
for i = Ncurve+NstraightBranch+1:Ncurve+NstraightBranch+ceil(Ncurve/2),
    theta = pi/2+(i-Ncurve-NstraightBranch-1)*dthetaHalf;
    yin = centeryTop+R2*sin(theta);
    xin = centerx1Top+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin, yin);
   % plot(xin,yin)
end

%outer quartercircle
for i = 1:Ncurve/4
    theta = -pi-(i)*dthetaQuart;
    yout = centeryTop+R1Top*sin(theta);
    xout = centerx1Top+R1Top*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout, yout);
   % plot(xout, yout)
end




fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pericardium
% Write out the vertex information

vertex_fid = fopen([mesh_name 'peri_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nperitot);

% make the top and bottom of the pericardium
for i=1:ceil(Nstraight/2),
    ytop = centery-(R2-(Dp-diameter)/2);
    xtop = -Lt/2+(i-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
    plot(xtop,ytop)
end

for i=ceil(Nstraight/2)+1:Nstraight,
    ybot = centery-R1-(Dp-diameter)/2;
    xbot = -Lt/2+(i-ceil(Nstraight/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
    plot(xbot,ybot)
end

% make the four side pieces
for i=Nstraight+1:Nstraight+ceil(Nperi/4),
    y = centery-(R1+(Dp-diameter)/2)+(i-Nstraight-1)*ds;
    x = -Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y)
end

for i=Nstraight+ceil(Nperi/4)+1:Nstraight+ceil(Nperi/2),
    y = centery-R2+(i-Nstraight-ceil(Nperi/4)-1)*ds;
    x = -Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y)
end

for i=Nstraight+ceil(Nperi/2)+1:Nstraight+ceil(3*Nperi/4),
    y = centery-(R1+(Dp-diameter)/2)+(i-Nstraight-ceil(Nperi/2)-1)*ds;
    x = Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y)
end

for i=Nstraight+ceil(3*Nperi/4)+1:Nperitot,
    y = centery-R2+(i-Nstraight-ceil(3*Nperi/4)-1)*ds;
    x = Lt/2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y)
end
fclose(vertex_fid);

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the spring information for the valveless pumping applied force

% spring_fid = fopen(['vp_aforce_' num2str(N) '.spring'], 'w');
% fprintf(spring_fid, '%d\n', NLa);
% 
% %elastic part of tube
% for i = 0:(NLa-1),
%     fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', i, i+NLa, Fmag, phase, 1);
% end
% 
% fclose(spring_fid);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the spring information for the peristalsis applied force
% 
% spring_fid = fopen(['peri_aforce_' num2str(N) '.spring'], 'w');
% fprintf(spring_fid, '%d\n', Nperist);
% 
% %elastic part of tube
% for i = 0:(Nperist-1),
%     fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', i, i+Nperist, Fmag, phase, 2);
% end
% 
% fclose(spring_fid);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write out the beam information for the elastic section

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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the target point information for the ends of the elastic tube
target_fid = fopen([mesh_name 'tube_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nstraight-4*Nend); 


Nstraight
for i = Nend:(Nstraight/2-Nend)-1,   
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = (Nstraight/2+Nend):(Nstraight-Nend)-1,   
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
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
% Write out the target point information for the actuator

%top actuator
% target_fid = fopen(['actuator_top_' num2str(N) '.target'], 'w');
% fprintf(target_fid, '%d\n', NLa);
% 
% for i = 0:NLa-1,
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% fclose(target_fid);
% 
% %bottom actuator
% target_fid = fopen(['actuator_bot_' num2str(N) '.target'], 'w');
% fprintf(target_fid, '%d\n', NLa);
% 
% for i = 0:NLa-1,
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the actuator

%top prescribed peristalsis
% target_fid = fopen(['pperi_top_' num2str(N) '.target'], 'w');
% fprintf(target_fid, '%d\n', NLap);
% 
% for i = 0:NLap-1,
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% fclose(target_fid);
% 
% %bottom prescribed peristalsis
% target_fid = fopen(['pperi_bot_' num2str(N) '.target'], 'w');
% fprintf(target_fid, '%d\n', NLap);
% 
% for i = 0:NLap-1,
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the racetrack
target_fid = fopen([mesh_name 'race_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nrace);

for i = 0:Nrace-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the pericardium
target_fid = fopen([mesh_name 'peri_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nperitot);

for i = 0:Nperitot-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integrating nick's code for pumping

 HeartTube(diameter, Lt, L, R1,R2, ds, centery, kappa_spring, kappa_beam, kappa_target)

