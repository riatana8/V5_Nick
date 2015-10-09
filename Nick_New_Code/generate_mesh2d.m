%% V5.2

%% 3/29 TS notes: 


clear all
%close all


%% In update_target_point_positions, changed R2 to 0.177 (From graph).
% Same in update_target_point_positions_peri


%% 
L = 1;                              % length of computational domain (m)
N = 1024;                           % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = L/(2*N);                       % space between boundary points in straight tube
ds2 = L/(4*N);
ppm = 1/ds2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the racetrack

Let = 1.34*0.1;                          % Length of elastic tube (m) scaled wrt to diameter of heart
% V5 TS. Nick's Heart_Tube.m code uses a fraction of the length of the tube to determine
% the number of rigid points on each end of elastic section. Made changes
% to use the same number he's using (can output these numbers from his code). 
%Nend = 25;
Nend = 40;                           % Number of rigid points on each end of elastic section
Lt = Let+2*Nend*ds2;                 % Length of main tube straight section with Nend rigid points on each end



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


Nstraight = 2*ceil(Lt/ds2)          % number of points along each shorter straight section (top+bottom)
NstraightBranch = 2*ceil(LtTop/ds);	% number of points along longer straight sections
NstraightSide = 2*ceil(Ls/ds);	%number of points along each straight (vertical)  outer side section
Ncurve = 2*ceil(pi*R1/ds);          % number of points along a circle of diameter R1. 
                                       %Will use fractions of this number
                                       %to draw points for curved sections
NcurveSmall = 2*ceil(pi*R2Bot/ds)    %number of points along a circle of diameter R2Bot (small)


Nrace = ceil((1/2)*0.5*Nstraight+4*Ncurve + (3/2)*NstraightBranch + NstraightSide) %+NcurveSmall % number of points making up the racetrack part. More added later in code

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
Nperitot = Nperi + Nstraight/2;       % total number of pericardium points

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
kappa_spring = 1.0e-1;          % spring constant (Newton) FOR TUBE WHERE PERI OCCURS
kappa_beam = 1.0e-1;            % beam stiffness constant (Newton m^2) for TUBE WHERE PERI OCCURS
kappa_target = 2000;            % target point penalty spring constant (Newton)
Fmag = kappa_spring;            % this is my best guess at a reasonable applied force
kappa_beamR = 0.004;            % beam coefficient for RACETRACK
bForce = 1e9;          %THIS JUST SAYS WHAT BEAM FORCE SHOULD BE MUST DO CALC. TO FIND COEFF., k_beam
kappa_springR = 25000.0;         % spring coefficient for RACETRACK
phase = 0;                      % initial phase of the oscillating force, where F=Fmag*phase and phase = (1+sin(2*pi*f*t-pi/2));
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
    xbot = -Lt/2+(i-1)*ds2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
 %   %plot(xbot,ybot,'b*')
     n = n+1; 
end

%bottom part
for i=ceil(Nstraight/2)+1:Nstraight,
    ytop = centery-R2;
    xtop = -Lt/2+(i-ceil(Nstraight/2)-1)*ds2;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
  %  %plot(xtop,ytop,'k*')
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
  %  %plot(x,y)
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
%     %plot(xtop,ytop);
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
%      %plot(xbot,ybot);
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
%      %plot(xtop,ytop);
% end
% 
% %bottom points
% for i=Na1:Na2,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%      %plot(xbot,ybot);
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
%      %plot(xtop,ytop);
% end
% 
% %bottom points
% for i=Nperi1:Nperi2,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%      %plot(xbot,ybot);
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
%     %plot(xtop,ytop,'r*')
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
%     %plot(xbot,ybot,'r*')
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%***************************** race track part **************************
% 
%******** calculations necessary before writing out vertex info:********
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%starting points for the left outer straight section
ylStart = yout+ds;
xlStart = xout+ds;

%additional parameters for middle branches:
centerx2Mid = xrStart-diameter;
centerx1Mid = xlStart+diameter;
Lm = abs (centerx2Mid-centerx1Mid);%length of middle branches (to ensure width of one diameter on either side, vertically)
NstraightMid = 2*ceil(Lm/ds);

Nrace = Nrace+NstraightMid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write out the vertex information for RACE TRACK!!!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen([mesh_name 'race_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nrace+1);

%%%right curved part of lower racetrack section

%inner semicircle
for i=1:ceil(Ncurve/2)
    theta = (i-1)*dthetaHalf- pi/2;
    yin(i) = centery+R2*sin(theta);
    xin(i) = Lt/2+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i),yin(i),'*'); hold on;
    %pause(0.01);
end
Ninfo(1) = ceil(Ncurve/2);
clear xin yin; 
%pause();

%outer quarter circle
for i=ceil(Ncurve/2)+1:ceil(3*Ncurve/4),
    jj=i-ceil(Ncurve/2);
    theta=(i-Ncurve/2-1)*dthetaQuart-(pi/2);
    yout(jj) = centery+R1*sin(theta);
    xout(jj) = Lt/2+R1*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout(jj), yout(jj) );
    %plot(xout(jj), yout(jj),'*'); hold on;
    %pause(0.001);
end
Ninfo(2) = length(xout)+Ninfo(1);
%pause();

%starting points for the right outer straight section
yrStart = yout(end)+ds;
xrStart = xout(end)+ds;
clear xout yout;


%right outer straight section
for i = 1:ceil(NstraightSide/2)
    xRside(i) = xrStart;
    yRside(i) = yrStart +i*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xRside(i), yRside(i) );
    %plot(xRside(i),yRside(i) ,'r*'); hold on;
    %pause(0.01);
end
Ninfo(3) = length(xRside) + Ninfo(2);
clear xRside yRside;
%pause();



%straight section for the upper part of the lower section (bottom wall of
%bottom branch)
for i = Ncurve+1:Ncurve+ceil(Nstraight/4),
    jj=i-Ncurve;
    yin(jj) = centery+R2;
    xin(jj) = centerx2-(i-Ncurve-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(jj), yin(jj) );
    %plot(xin(jj),yin(jj),'b*'); hold on;
    %pause(0.01);
%%%%
end
Ninfo(4) = length(xin) + Ninfo(3);
clear xin yin;
%pause();

%%%left curved part of racetrack
%inner semicircle
for i = 1:ceil(Ncurve/2),
    theta = pi/2+i*dthetaHalf;
    yin(i) = centery+R2*sin(theta);
    xin(i) = centerx1+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i), yin(i), 'g*'); hold on;
    %pause(0.001);
end
Ninfo(5) = length(xin)+Ninfo(4);
clear xin yin;
%pause();

%outer quarter circle
for i = 1:ceil(Ncurve/4)
    theta = (-pi/2)-(i)*dthetaQuart;
    yout(i) = centery+R1*sin(theta);
    xout(i) = centerx1+R1*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout(i), yout(i) );
    %plot(xout(i), yout(i),'k*'); hold on;
    %pause(0.001);
end
Ninfo(6) = length(xout) + Ninfo(5);
%pause();

ylStart = yout(end)+ds;
xlStart = xout(end)+ds;
clear xout yout;

%left outer straight section
for i = 1:ceil(NstraightSide/2)
     xLside(i) = xlStart;
     yLside(i) = ylStart + (i)*ds;
     fprintf(vertex_fid, '%1.16e %1.16e\n', xLside(i), yLside(i) );
     %plot(xLside(i),yLside(i),'r*'); hold on;
     %pause(0.001);
end
Ninfo(7) = length(xLside) + Ninfo(6);
clear xLside yLside;
%pause();


%%Middle Branches


%%top wall of lower skinny branch
for i = 1:ceil(NstraightMid/2)
    yin(i) = centeryMid-R2Bot;
    xin(i) = centerx2Mid-(i)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i),yin(i),'m*'); hold on;
    %pause(0.001);
end
Ninfo(8) = length(xin) + Ninfo(7);
clear xin yin;
%pause();

%%bottom wall of upper skinny branch
for i = 1:ceil(NstraightMid/2)
    yin(i) = centeryMid+R2Bot;
    xin(i) = centerx2Mid-(i)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i),yin(i),'r*'); hold on;
    %pause(0.001);
end
Ninfo(9) = length(xin) + Ninfo(8);
clear xin yin;
%pause();

%%Left inner semi circle

for i = 1:ceil(Ncurve/2),
    theta = pi/2+i*dthetaHalf;
    yin(i) = centeryMid+R2Bot*sin(theta);
    xin(i) = centerx1Mid+R2Bot*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i), yin(i),'k*'); hold on;
    %pause(0.001); 
end
Ninfo(10) = length(xin) + Ninfo(9);
clear xin yin;
%pause();

%% right inner semi circle
for i=1:ceil(Ncurve/2)
    theta = (i-1)*dthetaHalf- pi/2;
    yin(i) = centeryMid+R2Bot*sin(theta);
    xin(i) = centerx2Mid+R2Bot*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i),yin(i),'*'); hold on;
    %pause(0.001);
end
Ninfo(11) = length(xin) + Ninfo(10);
clear xin yin;
%pause();


%%%%%%%%%% upper branches

%right curved part of racetrack
for i=1:ceil(Ncurve/2)
    theta = (i-1)*dthetaHalf- (pi/2);
    yin(i) = centeryTop+R2*sin(theta);
    xin(i) = centerx2Top+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(i), yin(i) );
    %plot(xin(i),yin(i),'k*'); hold on;
    %pause(0.001);
end
Ninfo(12) = length(xin) + Ninfo(11);
clear xin yin;
%pause();


for i=1:ceil(Ncurve/4),
    theta=0+i*dthetaQuart;
    yout(i) = centeryTop+R1Top*sin(theta);
    xout(i) = centerx2Top+R1Top*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout(i), yout(i) );
    %plot(xout(i), yout(i),'m*'); hold on;
    %pause(0.001);
end
Ninfo(13) = length(xout) + Ninfo(12);
clear xout yout;
%pause();


%straight (horizontal) sections on the upper chamber 
for i = Ncurve+ceil(NstraightBranch/2)+1:Ncurve+NstraightBranch,
    jj = i - Ncurve - ceil(NstraightBranch/2);
    yout(jj) = centeryTop-R2;
    xout(jj) = centerx2Top-(i-Ncurve-ceil(NstraightBranch/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout(jj), yout(jj) );
    %plot(xout(jj), yout(jj),'g*'); hold on;
    %pause(0.001);
end
Ninfo(14) = length(xout) + Ninfo(13);
clear xout yout;
%pause();

%make the walls of the upper
%branch
for i = Ncurve+1:Ncurve+ceil(NstraightBranch/2),
    jj = i - Ncurve;
    yin(jj) = centeryTop+R2;
    xin(jj) = centerx2Top-(i-Ncurve-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(jj), yin(jj) );
    %plot(xin(jj),yin(jj),'b*'); hold on;
    %pause(0.001);
end
Ninfo(15) = length(xin) + Ninfo(14);
clear xin yin;
%pause();

for i = Ncurve+ceil(NstraightBranch/2)+1:Ncurve+NstraightBranch,
    jj = i - Ncurve - ceil(NstraightBranch/2);
    yout(jj) = centeryTop+R1Top;
    xout(jj) = centerx2Top-(i-Ncurve-ceil(NstraightBranch/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout(jj), yout(jj) );
    %plot(xout(jj), yout(jj),'g*'); hold on;
    %pause(0.001);
end
Ninfo(16) = length(xout) + Ninfo(15);
clear xout yout;
%pause();


%%%left curved part of upper racetrack

%inner semicircle
for i = Ncurve+NstraightBranch+1:Ncurve+NstraightBranch+ceil(Ncurve/2),
    theta = pi/2+(i-Ncurve-NstraightBranch-1)*dthetaHalf;
    jj = i - Ncurve - NstraightBranch;
    yin(jj) = centeryTop+R2*sin(theta);
    xin(jj) = centerx1Top+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xin(jj), yin(jj) );
    %plot(xin(jj),yin(jj),'k*'); hold on;
    %pause(0.001);
end
Ninfo(17) = length(xin) + Ninfo(16);
clear yin xin;
%pause();

%outer quartercircle
for i = 1:Ncurve/4
    theta = -pi-(i)*dthetaQuart;
    yout(i) = centeryTop+R1Top*sin(theta);
    xout(i) = centerx1Top+R1Top*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xout(i), yout(i) );
    %plot(xout(i), yout(i),'r*'); hold on;
    %pause(0.001);
end
Ninfo(18) = length(xout) + Ninfo(17);
clear xout yout;
%pause();

fprintf('\n\nBelow is Information about 18 pieces of geometry that make racetrack:\n');
fprintf('Ninfo(i) = how many Lag. Pts are in Racetrack AFTER piece "i" is added\n'); 
Ninfo'

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% **** Write out the SPRING information for the elastic section of RACETRACK ***
%
% ************************ ALL STRAIGHT PORTIONS!!!!!!! ****************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note: Need Springs On Piece {4,8,9,14,15,16}
%Note: Will Need Beams There As Well

spring_fid = fopen([mesh_name 'race_' num2str(N) '.spring'], 'w');

N4 = Ninfo(4) - Ninfo(3);
N8 = Ninfo(8) - Ninfo(7);
N9 = Ninfo(9) - Ninfo(8);
N14= Ninfo(14)- Ninfo(13);
N15= Ninfo(15)- Ninfo(14);
N16= Ninfo(16)- Ninfo(15);

Ntot = N4+N8+N9+N14+N15+N16+6;
%NOTE: +6 is to connect FIRST PT in piece to PREVIOUS POINT BEHIND IT

fprintf(spring_fid, '%d\n', Ntot);

%SPRINGS ALONG PIECE 4 (note connects flush w/ piece 5 -> no special case needed at end)
%1st: connect FIRST pt of PIECE 1 to FIRST pt in PIECE 4 (cpp:Ninfo(1)-1 -> Ninfo(3) )
%Last: connect LAST pt of PIECE 4 to FIRST pt in PIECE 5 (cpp:Ninfo(4)-1 -> Ninfo(4) )
for i = 0:N4
    if i==0
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(1)-1, Ninfo(3), kappa_springR*ds/(ds^2), ds);
    else
        ii=i-1+Ninfo(3); %FOR C++ IMPLEMENTATION
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, ii+1, kappa_springR*ds/(ds^2), ds);
    end
end

%SPRINGS ALONG PIECE 8 
%1st: connect FIRST pt of PIECE 11 to FIRST pt in PIECE 8 (cpp:Ninfo(10) -> Ninfo(7) )
%Last: connect LAST pt of PIECE 8 to LAST pt in PIECE 10 (cpp:Ninfo(8)-1 -> Ninfo(10)-1 )
for i = 0:N8
    if i==0
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(10), Ninfo(7), kappa_springR*ds/(ds^2), ds);
    elseif i<N8
        ii=(i-1)+Ninfo(7);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, ii+1, kappa_springR*ds/(ds^2), ds);
    else
        ii=(i-1)+Ninfo(7);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, Ninfo(10)-1, kappa_springR*ds/(ds^2), ds);
    end
end

%SPRINGS ALONG PIECE 9 
%1st: connect LAST pt of PIECE 11 to FIRST pt in PIECE 9 (cpp:Ninfo(11)-1 -> Ninfo(8) )
%Last: connect LAST pt of PIECE 9 to FIRST pt in PIECE 10 (cpp:Ninfo(9)-1 -> Ninfo(9) )
for i = 0:N9
    if i==0
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(11)-1, Ninfo(8), kappa_springR*ds/(ds^2), ds);
    elseif i<N9
        ii=(i-1)+Ninfo(8);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, ii+1, kappa_springR*ds/(ds^2), ds);
    else
        ii=(i-1)+Ninfo(8);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(9)-1, Ninfo(9), kappa_springR*ds/(ds^2), ds);
    end
end

%SPRINGS ALONG PIECE 14 
%1st: connect first pt of piece 12 to BEGIN piece 14 (cpp: Ninfo(11) -> Ninfo(13)
%Last: connect last pt of piece 14 to END of piece 17 (cpp: Ninfo(14)-1 -> Ninfo(17)-1
for i = 0:N14
    if i==0
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(11), Ninfo(13), kappa_springR*ds/(ds^2), ds);
    elseif i<N14
        ii=(i-1)+Ninfo(13);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, ii+1, kappa_springR*ds/(ds^2), ds);
    else
        ii=(i-1)+Ninfo(13);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(14)-1, Ninfo(17)-1, kappa_springR*ds/(ds^2), ds);
    end
end

%SPRINGS ALONG PIECE 15 
%1st: connect LAST pt of piece 12 to FIRST pt in piece 15 (cpp: Ninfo(12)-1 -> Ninfo(14)
%Last: connect last pt of piece 15 to FIRST pt of piece 17 (cpp: Ninfo(15)-1 -> Ninfo(16)
for i = 0:N15
    if i==0
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(12)-1, Ninfo(14), kappa_springR*ds/(ds^2), ds);
    elseif i<N15
        ii=(i-1)+Ninfo(14);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, ii+1, kappa_springR*ds/(ds^2), ds);
    else
        ii=(i-1)+Ninfo(14);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(15)-1, Ninfo(16), kappa_springR*ds/(ds^2), ds);
    end
end


%SPRINGS ALONG PIECE 16 
%1st: connect LAST pt of piece 13 to FIRST pt in piece 16 (cpp: Ninfo(13)-1 -> Ninfo(15)
%Last: connect last pt of piece 16 to LAST pt of piece 18 (cpp: Ninfo(16)-1 -> Ninfo(18)-1 
for i = 0:N16
    if i==0
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(13)-1, Ninfo(15), kappa_springR*ds/(ds^2), ds);
    elseif i<N16
        ii=(i-1)+Ninfo(15);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', ii, ii+1, kappa_springR*ds/(ds^2), ds);
    else
        ii=(i-1)+Ninfo(15);
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Ninfo(16)-1, Ninfo(18)-1, kappa_springR*ds/(ds^2), ds);
    end
end


fclose(spring_fid);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% **** Write out the BEAM information for the elastic section of RACETRACK ***
%
% ************************ ALL STRAIGHT PORTIONS!!!!!!! ****************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam_fid = fopen([mesh_name 'race_' num2str(N) '.beam'], 'w');

N4 = Ninfo(4) - Ninfo(3);
N8 = Ninfo(8) - Ninfo(7);
N9 = Ninfo(9) - Ninfo(8);
N14= Ninfo(14)- Ninfo(13);
N15= Ninfo(15)- Ninfo(14);
N16= Ninfo(16)- Ninfo(15);

Ntot = N4+N8+N9+N14+N15+N16;

fprintf(beam_fid, '%d\n', Ntot);

fprintf('The BEAM FORCE for RACETRACK is: %d\n',bForce);

%BEAMS ALONG PIECE 4 (note connects flush w/ piece 5 -> no special case needed at end)
%1st: connect FIRST pt of PIECE 1 to FIRST pt in PIECE 4 (cpp:Ninfo(1)-1 -> Ninfo(3) )
%Last: connect LAST pt of PIECE 4 to FIRST pt in PIECE 5 (cpp:Ninfo(4)-1 -> Ninfo(4) )
for i = 1:N4
    if i==1
        ii=i-1+Ninfo(3); %FOR C++ IMPLEMENTATION
        %fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(1)-1, ii, ii+1, kappa_beamR*ds2/(ds2^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(1)-1, ii, ii+1, bForce);
    else
        ii=i-1+Ninfo(3); %FOR C++ IMPLEMENTATION
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, kappa_beamR*ds2/(ds2^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, bForce);
    end
end

%BEAMS ALONG PIECE 8 
%1st: connect FIRST pt of PIECE 11 to FIRST pt in PIECE 8 (cpp:Ninfo(10) -> Ninfo(7) )
%Last: connect LAST pt of PIECE 8 to LAST pt in PIECE 10 (cpp:Ninfo(8)-1 -> Ninfo(10)-1 )
for i = 1:N8
    if i==1
        ii=(i-1)+Ninfo(7);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(10), ii, ii+1, kappa_beamR*ds2/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(10), ii, ii+1, bForce);
    elseif i<N8
        ii=(i-1)+Ninfo(7);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, bForce);
    else
        ii=(i-1)+Ninfo(7);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(10)-1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(10)-1, bForce);
    end
end

%BEAMS ALONG PIECE 9 
%1st: connect LAST pt of PIECE 11 to FIRST pt in PIECE 9 (cpp:Ninfo(11)-1 -> Ninfo(8) )
%Last: connect LAST pt of PIECE 9 to FIRST pt in PIECE 10 (cpp:Ninfo(9)-1 -> Ninfo(9) )
for i = 1:N9
    if i==1
        ii=(i-1)+Ninfo(8);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(11)-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(11)-1, ii, ii+1, bForce);
    elseif i<N9
        ii=(i-1)+Ninfo(8);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, bForce);
    else
        ii=(i-1)+Ninfo(8);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(9), kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(9), bForce);
    end
end

%BEAMS ALONG PIECE 14 
%1st: connect first pt of piece 12 to BEGIN piece 14 (cpp: Ninfo(11) -> Ninfo(13)
%Last: connect last pt of piece 14 to END of piece 17 (cpp: Ninfo(14)-1 -> Ninfo(17)-1
for i = 1:N14
    if i==1
        ii=(i-1)+Ninfo(13);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(11), ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(11), ii, ii+1, bForce);
    elseif i<N14
        ii=(i-1)+Ninfo(13);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, bForce);       
    else
        ii=(i-1)+Ninfo(13);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(17)-1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(17)-1, bForce);
    end
end


%BEAMS ALONG PIECE 15 
%1st: connect LAST pt of piece 12 to FIRST pt in piece 15 (cpp: Ninfo(12)-1 -> Ninfo(14)
%Last: connect last pt of piece 15 to FIRST pt of piece 17 (cpp: Ninfo(15)-1 -> Ninfo(16)
for i = 1:N15
    if i==1
        ii=(i-1)+Ninfo(14);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(12)-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(12)-1, ii, ii+1, bForce);
    elseif i<N15
        ii=(i-1)+Ninfo(14);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, bForce);
    else
        ii=(i-1)+Ninfo(14);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(16), kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(16), bForce);
    end
end


%BEAMS ALONG PIECE 16 
%1st: connect LAST pt of piece 13 to FIRST pt in piece 16 (cpp: Ninfo(13)-1 -> Ninfo(15)
%Last: connect last pt of piece 16 to LAST pt of piece 18 (cpp: Ninfo(16)-1 -> Ninfo(18)-1 
for i = 1:N16
    if i==1
        ii=(i-1)+Ninfo(15);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(13)-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', Ninfo(13)-1, ii, ii+1, bForce );
    elseif i<N16
        ii=(i-1)+Ninfo(15);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, ii+1, bForce );        
    else
        ii=(i-1)+Ninfo(15);
        %fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(18)-1, kappa_beamR*ds/(ds^4));
        fprintf(beam_fid, '%d %d %d %1.16e\n', ii-1, ii, Ninfo(18)-1,  bForce);
    end
end


fclose(beam_fid);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ****************************** Pericardium ******************************
% 
% ******************** Write out the vertex information *******************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vertex_fid = fopen([mesh_name 'peri_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', Nperitot);
% 
% % make the top and bottom of the pericardium
% for i=1:ceil(Nstraight/4),
%     ytop = centery-(R2-(Dp-diameter)/2);
%     xtop = -Lt/2+i*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
%    % %plot(xtop,ytop)
% end
% 
% for i=ceil(Nstraight/4)+1:ceil(Nstraight/2),
%     ybot = centery-R1-(Dp-diameter)/2;
%     xbot = Lt/2+(i-ceil(Nstraight/2)-1)*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
%     %plot(xbot,ybot)
% end
% 
% % make the four side pieces
% for i=ceil(Nstraight/2)+1:ceil(Nstraight/2)+ceil(Nperi/4),
%     y = centery-(R1+(Dp-diameter)/2)+(i-ceil(Nstraight/2)-1)*ds;
%     x = -Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
%     %plot(x,y,'r')
% end
% 
% for i=ceil(Nstraight/2)+ceil(Nperi/4)+1:ceil(Nstraight/2)+ceil(Nperi/2),
%     y = centery-R2+(i-ceil(Nstraight/2)-ceil(Nperi/4)-1)*ds;
%     x = -Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
%    % %plot(x,y,'b')
% end
% 
% for i=ceil(Nstraight/2)+ceil(Nperi/2)+1:ceil(Nstraight/2)+ceil(3*Nperi/4),
%     y = centery-(R1+(Dp-diameter)/2)+(i-ceil(Nstraight/2)-ceil(Nperi/2)-1)*ds;
%     x = Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
%    % %plot(x,y,'g')
% end
% 
% for i=ceil(Nstraight/2)+ceil(3*Nperi/4)+1:Nperitot,
%     y = centery-R2+(i-ceil(Nstraight/2)-ceil(3*Nperi/4)-1)*ds;
%     x = Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
%    % %plot(x,y,'k')
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% **** Write out the spring information for the elastic section of TUBE *****
%
% ************ WHERE DYNAMIC SUCTION PUMPING (PERISTALSIS) WILL OCCUR!!!! ****************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%spring_fid = fopen([mesh_name 'tube_' num2str(N) '.spring'], 'w');
%fprintf(spring_fid, '%d\n', Nstraight-2);

%elastic part of tube
%for i = 0:ceil(Nstraight/2)-2,
%    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds2/(ds2^2), ds2);
%end

%for i = ceil(Nstraight/2):Nstraight-2,
%    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds2/(ds2^2), ds2);
%end

%fclose(spring_fid);



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
% 
% beam_fid = fopen([mesh_name 'tube_' num2str(N) '.beam'], 'w');
% fprintf(beam_fid, '%d\n', Nstraight-4);
% 
% %elastic part of tube
% for i = 0:ceil(Nstraight/2)-3,
%     fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds2/(ds2^4));
% end
% 
% for i = ceil(Nstraight/2):Nstraight-3,
%     fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds2/(ds2^4));
% end
% fclose(beam_fid);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the target point information for the ends of the elastic tube
target_fid = fopen([mesh_name 'tube_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nstraight); 


Nstraight
for i = 0:(Nstraight/2)-1,   
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds2/(ds2^2));
end

for i = (Nstraight/2):(Nstraight)-1,   
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
%
% ************ Write out the TARGET point information for the RACETRACK **********
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_fid = fopen([mesh_name 'race_' num2str(N) '.target'], 'w');

N1 = Ninfo(1);
N2 = Ninfo(2) - Ninfo(1);
N3 = Ninfo(3) - Ninfo(2);
N5= Ninfo(5)- Ninfo(4);
N6= Ninfo(6)- Ninfo(5);
N7= Ninfo(7)- Ninfo(6);
N10 = Ninfo(10) - Ninfo(9);
N11= Ninfo(11)- Ninfo(10);
N12= Ninfo(12)- Ninfo(11);
N13= Ninfo(13)- Ninfo(12);
N17= Ninfo(17)- Ninfo(16);
N18= Ninfo(18)- Ninfo(17);

Ntot = N1+N2+N3+N5+N6+N7+N10+N11+N12+N13+N17+N18;

fprintf(target_fid, '%d\n', Ntot);

% Target Pts along Piece 1
for i = 0:Ninfo(1)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 2 + 3
for i = Ninfo(1):Ninfo(3)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 5
for i = Ninfo(4):Ninfo(5)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 6 + 7
for i = Ninfo(5):Ninfo(7)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 10
for i = Ninfo(9):Ninfo(10)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 11
for i = Ninfo(10):Ninfo(11)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 12
for i = Ninfo(11):Ninfo(12)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 13
for i = Ninfo(12):Ninfo(13)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 17
for i = Ninfo(16):Ninfo(17)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


% Target Pts along Piece 18
for i = Ninfo(17):Ninfo(18)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end


fclose(target_fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *************** TARGET POINT information for the PERICARDIUM *******************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% target_fid = fopen([mesh_name 'peri_' num2str(N) '.target'], 'w');
% 
% fprintf(target_fid, '%d\n', Nperitot);
% 
% for i = 0:Nperitot-1,
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Integrating nick's code for pumping

 HeartTube(diameter, Lt, L, R1,R2, ds2, centery, kappa_spring, kappa_beam, kappa_target)

