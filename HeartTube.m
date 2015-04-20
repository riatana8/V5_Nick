%  function HeartTube()
 function HeartTube(d, L, L_d, R_o,R_i, ds, centery, spring_force, beam_force, kappa_target)

%This script creates the input files for HeartTube simulations using
%Dynamic Suction Pumping w/ Biconcave Disc "Red Blood Cells" following the
%equation in Crowl and Fogelson; however, it 'reflects' the disc, making 
%a more butterfly-ish shape.

%NOTE: THE blood_cells.spring file needs to be updated because of the
%construction of new RBC shape!!!!

%It has equally spaced grid points on the straight tube on the bottom
%This function lets you specify the number of red blood cells and the
%volume fraction. 
% %%V5 TS Parameters%%%%%%%%%%%%%%%%%
% d = .1
% L = 0.134
% L_d = 1
% R_o = .1
% R_i = .2
% ds = L_d/(2*512)
% centery = -0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = 0.2;               %Diameter of tube
% L = 1.5;               %Length of Straight-Leg in Tube
% L_d= 5.0;              %Length of computational domain [-L_d/2, L_d/2]
% R_o = 0.75;            %Radius of the Outer Bndry Circle (NOTE: R_o > d)        
% R_i = R_o - d;         %Radius of the Inner Bndry Circle 
% ds= L_d / (2*1024);    %Lagrangian Point spacing

%Constructs the Outer and Inner Heart Tube Boundaries
%NOTE: num_o/i tells number of points on bottom flat section (used for target point
%input textfile)
[X_tube,Y_tube,NTvec,Info_Target] = Make_Boundary(ds,L,R_o,R_i,d,centery);



% %Parameters needed for constructing the blood cells in the heart tube
% %Also finds center of each blood cell and then places a blood cell (circle)
% %around it with npts.
% Vol = Volume(d,L,R_o,R_i);                       %Function computes volume of the heart tube
% Vf = 0.15;                                       %Volume Fraction of Red Blood Cells
% Nc = 50;                                         %Number of Red Blood Cells
% [x_m, y_m ] = Put_RBCs(Nc,L,d,R_i,R_o);          %Find pts along axis where Red Blood Cells go
% 
% %Creates blood cell -> finds points around center of each RBC
% %Choice of what kind of blood cell desired (circular, elliptical,'biconcave' projection)
% %
% %[X_c,Y_c,Npts] = Make_Circle_Cells(Nc,Vol,Vf,x_m,y_m,ds);  
% %
% e = 0.9;    %Ellipticity of Blood cells
% [X_c,Y_c,NCvec] = Make_Elliptical_Cells(Nc,Vol,Vf,L,e,x_m,y_m,ds);
% %
% %pts = 110;    %Number of Points for each biconcave red blood cell 
% %[X_c,Y_c] = Make_Biconcave_Cells(Nc,Vol,Vf,L,x_m,y_m,pts);
% 
% 
% 
% %For HeartTube Inputs
% %kappa_spring = 500.0; %up from 5.0
% spring_force = 10e7; %makes a uniform spring force around tube geometry
% %kappa_beam = 0.00000001; %up from 0.0001
% beam_force = 0.00014; %htc_31 was 1.4
% kappa_target = 10000.0;
% %k_btwn = 0.05;

% %For Cell Inputs
% k_spring_cell = 5000.0;
% k_beam_cell = 0.0000001;
% 
%% 03/29/15 V5 TS note:  Commenting write_Input_Files out to use only generate_mesh2d to 
% create vertex, spring, and beam files for heart tube and attached
% morphology.

% write_Input_Files(NTvec,NCvec,X_tube,Y_tube,X_c,Y_c)



% for i=1:length(X_tube)
%    plot(X_tube(i),Y_tube(i),'*'); hold on;
%    axis([-R_o-3*L/4 R_o+3*L/4 -R_o-d R_o+d]);
%    pause(0.0005);
% end


%PLOTS IT
%plot_it_all(X_tube,Y_tube,x_m,y_m,X_c,Y_c,R_o,L,Nc,Vf,d);    
% figure
% plot(X_tube,Y_tube,'o')
% 
% fprintf(' \n\nINFO FOR UPDATE-TARGET-PT-FILE: \n\n');
% fprintf('w (wave-width) = %d\n',Info_Target(1));
% fprintf('R_o (outer radius) = %d\n',Info_Target(2));
% fprintf('R_i (inner radius) = %d\n',Info_Target(3));
% fprintf('A (amplitude) = %d\n',Info_Target(4));
% fprintf('A_tilde (revised amplitude) = %d\n',Info_Target(5));
% fprintf('Outer Bottom: 1st Pt = %d, 2nd Pt = %d\n',Info_Target(6),Info_Target(7));
% fprintf('Inner Bottom: 1st Pt = %d, 2nd Pt = %d\n',Info_Target(8),Info_Target(9));
% fprintf('d (diameter of tube) = %d\n\n\n',d);
% 
% 


    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of entire Heart Tube Boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,NVec,Info_Target] = Make_Boundary(ds,L,R_o,R_i,d,centery)
    X= R_o
    Y=R_i
    NVec=[1,1];
    %L_tube = 2*L + 2*pi*R;
    %Constructs vectors containing pts for straight portion (including ends x=[-L/2,L/2] )
    [xS yS] = make_Flat_Portion(ds,L,centery);
    xS_T = xS(end:-1:1); %x-Values for top of tube
    yS_T_o = yS - R_o;   %y-Values for OUTER top of tube
    yS_T_i = yS - R_i;   %y-Values for INNER top of tube
   % figure
    plot(xS, yS_T_o,'r-')
    hold on
    plot(xS, yS_T_i,'r-')
    
    [yS_B yS_B2 Nbot Info_Target] = make_Gaussian_Wave(xS,d,L,R_o,R_i); %Gives Final Position for Wave
    yP_o_2 = centery - yS_B - R_o;     %y-Values for OUTER bottom of tube for PHASE 2 (right when peri. starts)
    yP_i_2 = centery +yS_B - R_i;    %y-Values for INNER bottom of tube for PHASE 2
   % figure
    plot(xS, yP_o_2)
    hold on
    plot (xS, yP_i_2)
    
    yP_o_3 = centery - yS_B2 - R_o;    %y-Values for OUTER bottom of tube for PHASE 3 (right when peri. ends)
    yP_i_3 = centery + yS_B2 - R_i;   %y-Values for INNER bottom of tube for PHASE 3
   % figure
    plot(xS, yP_o_3)
    hold on
    plot(xS, yP_i_3)
    %%%% V5 TS made all entries in yB -0.1 instead of 0 (our model has a center
% y coord of -0.1 not 0
    yB = zeros(1,length(xS));
    yB(1,1:end) = centery;
    yB_o = yB - R_o;         %y-Values for OUTER bottom of tube for PHASE 1 (all flat)
    yB_i = yB - R_i;         %y-Values for INNER bottom of tube for PHASE 1
    xP_B = xS;               %x-Values for bottom of tube
                 %Nbot(1);   % # pts. before peri. region
                 %Nbot(2);   % # pts. in peri. region
                 %Nbot(3)    % # pts. total for bottom region
 
%     figure(2)
%     plot(xS,yB_o,'r*'); hold on;
%     plot(xS,yB_i,'r*'); hold on;
% 
%     figure(3)
%     plot(xS,yP_o_2,'k*'); hold on;
%     plot(xS,yP_i_2,'m*'); hold on;
%     
%     figure(4)
%     plot(xS,yP_o_3,'c*'); hold on;
%     plot(xS,yP_i_3,'g*'); hold on;
    
    length(xS)
    
    %PRINTS BOTTOM-TUBE PTS. TO .txt FILES%
    print_Them_Vertices(xS,'X.txt'); %Prints x-Values (all same regardless top, bottom)
    print_Them_Vertices(yS_T_o,'yOut_1.txt'); % changed to yS_T_o from yB_o 4/20
    print_Them_Vertices(yP_o_2,'yOut_2.txt');
    print_Them_Vertices(yP_o_3,'yOut_3.txt');
    print_Them_Vertices(yS_T_i,'yIn_1.txt');  % changed to yS_T_i from yB_i 4/20
    print_Them_Vertices(yP_i_2,'yIn_2.txt');
    print_Them_Vertices(yP_i_3,'yIn_3.txt');
    
    
    %Constructs vector of pts. for bottom peristalsis wave part
    %[Xbot,Ybot,Nbot,Info_Target] = make_Peristalsis_Wave(ds,d,L,R_o,R_i);
    %xP_B = Xbot;       %pts. for bottom straight-peri-straight part
    %yP_o = -Ybot - R_o;
    %yP_i = Ybot - R_i;
                        %Nbot(1);   % # pts. before peri. region
                        %Nbot(2);   % # pts. in peri. region
                        %Nbot(3)    % # pts. total for bottom region
    
%     %Constructs vectors containing pts for round OUTER portion (not-including ends)
%     half = 1; a=R_o;  b=R_o; ang = pi;
%     [X_o,Y_o] = make_Elliptical_Geometry(a,b,ds,half);
%     xR_o = X_o+L/2;
%     yR_o = Y_o;
%     [xL_o,yL_o] = rotate_Geometry(X_o,Y_o,ang);
%     xL_o = xL_o - L/2;
% 
%     
%     %Constructs vectors containing pts for round INNER portion (not-including ends)
%     half = 1;  a=R_i;   b=R_i; ang = pi;
%     [X_i,Y_i] = make_Elliptical_Geometry(a,b,ds,half);
%     xR_i = X_i+L/2;
%     yR_i = Y_i;
%     [xL_i,yL_i] = rotate_Geometry(X_i,Y_i,ang);
%     xL_i = xL_i - L/2;
% 
%     
%     %X = [xS_B xR_o xS_T xL_o xS_B xR_i xS_T xL_i];
%     %Y = [yS_B_o yR_o yS_T_o yL_o yS_B_i yR_i yS_T_i yL_i];
%     
%     X = [xP_B xR_o xS_T   xL_o xP_B xR_i xS_T xL_i];
%     Y = [yB_o yR_o yS_T_o yL_o yB_i yR_i yS_T_i yL_i];
%       
%     Nf = length(xS);     %Number of points in flat region of tube on top
%     NR_o = length(xR_o); %Number of points in round side of tube on outside
%     NR_i = length(xR_i); %Number of points in round side of tube on inside
%     
%     N1 = Nbot(1);                  % # of pts. before wave region starts on bottom (still flat portion)
%     N2 = Nbot(2);                  % # of pts. in peristalsis region on bottom (WAVE PART ONLY)
%     N3 = Nbot(3);                  % # of pts. in bottom part of tube (flat+peri+flat) up to rounded part
%     N4 = N3+NR_o;                  % # of pts leading up to straight-top on OUTER of tube
%     N5 = Nf;                       % # of pts on straight section of tube on TOP
%     N6 = N3+ Nf + 2*NR_o;          % # of pts leading up to straight-bottom portion on INSIDE of tube
%     N7 = 2*N3 + Nf + 2*NR_o + NR_i;% # of pts leading up to straight-top portion on INSIDE of tube
%     N8 = N7 + NR_i + Nf;           % # of pts on entire tube (OUTER + INNER)
%     
%     NVec = [N1 N2 N3 N4 N5 N6 N7 N8]
%     
%     %Gives Lagrangian Pt. #'s for Peristaltic Regions (for C++ indexing
%     %starting at 0...)
%     Info_Target(6) = N1;        % FIRST Pt. of Peristaltic Region on OUTER bottom
%     Info_Target(7) = N1+N2-1;   % LAST Pt. of Peristalstic Region on OUTER bottom
%     Info_Target(8) = N6+N1;     % FIRST Pt. of Peristaltic Region on INNER bottom
%     Info_Target(9) = N6+N1+N2-1;% LAST Pt. of Peristaltic Region on INNER bottom
%     
%     
%     %Prints Lagrangian-Pts. for Beginning of Peristaltic Wave
%     %fileID = fopen('Y_Peri_Outer.txt','w');
%     %fprintf(fileID,'%1.16e\n',yP_o_2);
%     %fclose(fileID);
%     %
%     %fileID2 = fopen('Y_Peri_Inner.txt','w');
%     %fprintf(fileID2,'%1.16e\n',yP_i_2);
%     %fclose(fileID2);
%     
%     return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of a Boundary: Makes beginning of sine wave for peristalsis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

function [yS_B yS_B2 NVec Info_Target] = make_Gaussian_Wave(xS,d,L,R_o,R_i)
%% v5 TS Having issues resolving the number of points in heart tube. Thought: 
% Nick's  A and w make sense with a d of 0.2 and an L of 1.5. The width of
% our wave will be almost the length of our tube.
% Any adjustments needed for a geometery that is not centered at a y = 0 is
% made where this code is returned.

frac = 0.1;     %frac of straight-portion not doing peristalsis
Am = 0.3*d      %amplitude of peristalsis wave
w = 0.1*d          %width of wave

Lend = -L/2 + frac*L/2; %Last pt. going R->L on straight part on LHS of peri. region
Rend = L/2 -  frac*L/2; %First pt. going R-> on straight part of RHS of peri. region


yS_B= zeros(1,length(xS));
yS_B2 = yS_B;

%Find pt BEFORE peristaltic region
n=1; x=xS(1); flag = 0;
while flag == 0
    if x > Lend
        N_beg_P = n-1;
        flag = 1;
    end
    n = n+1; x = xS(n);
end

%Find pt AFTER peristaltic region
n=1; x=xS(1); flag = 0;
while flag == 0
    if x > Rend
        N_end_P = n;
        flag = 1;
    end
    n = n+1; x = xS(n);
end

%Find pt. AFTER width, w.
n=1; x=xS(1); flag = 0;
while flag == 0
    if x > ( xS(N_beg_P+1)+w )
        Nw = n;
        flag = 1;
    end
    n = n+1; x = xS(n);
end

%Find pt. AFTER width, w. BEFORE end of Peristalsis Region
n=N_end_P+1; x=xS(n); flag = 0;
while flag == 0
    if x < ( xS(N_end_P-1)-w )
        Nw_2 = n;
        flag = 1;
    end
    n = n-1; x = xS(n);
end

%Get model amplitude from desired amplitude, A:
A_tilde = get_Amplitude_For_Wave(Am,w)

x0 = xS(N_beg_P+1);  x1 = xS(Nw); xC = (x0+x1)/2;
for i=N_beg_P+1:Nw
   yS_B(i) = A_tilde*(xS(i)-x0)^2*(x1-xS(i))^2*exp( -(xS(i)-xC)^2 / (w/2)^2 ); 
end

x0 = xS(Nw_2);  x1 = xS(N_end_P-1); xC = (x0+x1)/2;
for i=Nw_2:N_end_P
   yS_B2(i) = A_tilde*(xS(i)-x0)^2*(x1-xS(i))^2*exp( -(xS(i)-xC)^2 / (w/2)^2 ); 
end

    Nbef =  N_beg_P;      % # pts. before peri. region
    Nperi = Nw - N_beg_P; % # pts. in peri. region
    Ntot  = length(xS) ;  % # pts. total for bottom region
        
    NVec = [Nbef Nperi Ntot]

    %Store Info for update_target_points.C file
    Info_Target(1) = w;      %Stores width
    Info_Target(2) = R_o;    %Stores Outer Radius
    Info_Target(3) = R_i;    %Stores Inner Radius
    Info_Target(4) = Am;      %Stores Amplitude of Sine Wave
    Info_Target(5) = A_tilde;%Stores Revised Amplitude for Wave-Function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: get Amplitude for Peristaltic Wave
% Note: -Based on what wave-form you use. Need to get data from regression
%        and then invert to find necessary height from height desired.
%       -i.e., A is desired height, A_tilde is height needed for use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A_tilde = get_Amplitude_For_Wave(A,w)

%A: desired wave height

%Note: A = m*A_Tilde + b  from simulations in test_Peri_Shape.m (needs to
%know width, w.)

m = 1.0e-0;
b = 0.0;

%A_tilde = (A - b)/m;  %for linear regression data <--- Note depends on lots of parameters

A_tilde = 16/w^4*A; % <-- Comes from evaluating wave at x_C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of a Boundary: Makes beginning of sine wave for peristalsis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% function [Xnew,Ynew,Nvec,Info_Target] = make_Peristalsis_Wave(ds,d,L,R_o,R_i)
% 
% frac = 0.9;     %frac of straight-portion not doing peristalsis
% A = 0.4*d;      %amplitude of peristalsis wave
% 
% Lend = -L/2 + frac*L/2; %Last pt. going R->L on straight part on LHS of peri. region
% Rend = L/2 -  frac*L/2; %First pt. going R-> on straight part of RHS of peri. region
% Lnew = Rend;
% 
% xL = [-L/2:ds:Lend Lend];
% yL = zeros(1,length(xL));
% xR = [Rend:ds:L/2 L/2];
% yR = zeros(1,length(xR));
% 
% n = 1;          %counting for array
% xp = Lend;      %first x-value at (Lend,0) [Note: doesn't store]
% yp = 0;         %first y-value at (Lend,0) [NOTE: doesn't store]
% tol = ds / 10;  %error tolerance for root-finding algorithm
% 
%     while xp < Rend; 
%             
%         xn = xp + ds;      %Guess for next x-value
%         
%         xfar = xp + 2.0*ds;  %Far guess for next x-value
%             
%         %initiating guess for bisection-algorithm
%         yn = A*sin(pi*xn/Lnew);
%         errSign = ( ds - sqrt( (xn-xp)^2 + (yn-yp)^2 ) );
%         err = abs(errSign);
%     
%         %Bisection algorithm to make points equally spaced
%         ct = 0;
%         while ( err > tol )
%          
%             if errSign < 0
%                 xfar = xn;
%                 xn = (xn+xp)/2;
%             elseif errSign > 0
%                 xp = xn;
%                 xn = (xn+xfar)/2;
%             end
%             
%             if xn<xp
%                 fprintf('NOT CONVERGING AT x = %d\n',tprev)
%                 break;  
%             end
%             
%             yn = A*sin(pi*xn/Lnew);
%             errSign = ( ds - sqrt( (xn-xp)^2 + (yn-yp)^2 ) );
%             err = abs(errSign);
%             
%             if ct>3
%                 fprintf('ct = %d\n',ct);
%                 fprintf('xp = %d\n',xp);
%                 fprintf('xn = %d\n',xn);
%                 fprintf('xfar = %d\n',xfar);
%                 fprintf('sqrt() = %d\n',sqrt( (xn-xp)^2 + (yn-yp)^2));
%                 
%                 fprintf('yn = %d\n',yn);
%                 fprintf('yp = %d\n',yp);
%                 fprintf('ds = %d\n',ds);
%                 fprintf('err = %d\n\n',errSign); 
%                 pause();
%             end
%             ct=ct+1;
%         end
%         
%         X(n) = xn;     %Store X-value
%         Y(n) = yn;     %Store Y-value
%         
%         xp= xn;     %Save prev. x-value of bisection
%         yp= yn;     %Save prev. y-value of bisection
%         
%         n = n+1;       %Update counter
%     
%         %plot(xn,yn,'*'); hold on;
%         %pause();
%         %axis([-L/2 L/2 -A A]);
%       
%     end
%         
%     %Length of vector before taking out extra point 
%     N = length(X);
% 
%     %Get rid of extra point
%     Xperi = X(1:N-1);
%     Yperi = Y(1:N-1);
% 
%     Xnew = [xL Xperi xR];
%     Ynew = [yL Yperi yR];
%     
%     %plot(Xnew,Ynew,'*'); hold on;
%     %axis([-L/2 L/2 -A A]);
%     %pause();
%     
%     Nbef = length(xL);                 % # pts. before peri. region
%     Nperi = length(Xperi);             % # pts. in peri. region
%     Ntot  = Nbef + Nperi + length(xR); % # pts. total for bottom region
%     
%     Nvec = [Nbef Nperi Ntot];
% 
%     %Store Info for update_target_points.C file
%     Info_Target(1) = Lnew;  %Stores period
%     Info_Target(2) = R_o;   %Stores Outer Radius
%     Info_Target(3) = R_i;   %Stores Inner Radius
%     Info_Target(4) = A;     %Stores Amplitude of Sine Wave
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of a Boundary: rotation function
% Note: rotates the (x,y)-pts by an angle, ang.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 
% function [Xnew,Ynew] = rotate_Geometry(x,y,ang)
% 
% for i=1:length(x)
%    Xnew(i) = x(i)*cos(ang) - y(i)*sin(ang);
%    Ynew(i) = x(i)*sin(ang) + y(i)*cos(ang);
% end
% 
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of a Boundary: makes flat portion of tube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x y] = make_Flat_Portion(ds,L,centery)

n=1;         %Counter
xVal = -L/2; %Left-most point of straight-tube
x(n) = xVal; %Saves that point.

%Stores x-values for straight-portion
xVal = xVal+ds;
while xVal < L/2
    n = n+1;
    x(n) = xVal;
    xVal = xVal+ds;
end

%Makes sure we store the right most end-pt
len = length(x);
x(len+1) = L/2;

%Stores y-values at 0.
len2 = length(x);
%% V5 TS 3/29/15 changed y(lens) = 0 to = centery
y(1:len2) = centery; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Construct on boundary: Makes Right/Left Arcs.
% Note that the counting starts at the point (0,-b) or the top point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [X_C,Y_C] = make_Elliptical_Geometry(a,b,ds,half)
% 
% %a = semi-major (horizontal) axis of ellipse 
% %b = semi-minor (vertical) axis of ellipse
% %ds = lag. point spacing
% 
% n = 1;           %counting for array
% Xprev = 0;       %first x-value at (0,-b) [Note: doesn't store]
% Yprev = -b;      %first y-value at (0,-b) [NOTE: doesn't store]
% tprev = -pi/2;   %initial angle 
% tol = ds / 10;   %error tolerance for root-finding algorithm
% maxy = max(a,b); %To be used for guessing next pt.
% 
% if half == 1
%     ang_end = pi/2;   %For half a circle
% else
%     ang_end = 3*pi/2; %For full circle
% end
% 
%     while tprev < ang_end 
%     
%         tnext = tprev + ds/maxy; %Guess for next angle value
%         
%         if tprev < -pi/4
%             tfar =  0.0;%abs(2*tprev);      %Far guess
%         else
%             tfar = 3*pi/2;
%         end
%             
%         %initiating guess for bisection-algorithm
%         xn = a*cos(tnext);
%         yn = b*sin(tnext);
%         errSign = ( ds - sqrt( (xn-Xprev)^2 + (yn-Yprev)^2 ) );
%         err = abs(errSign);
%     
%         %Bisection algorithm to make points equally spaced
%         while ( err > tol )
%         
%             if errSign < 0
%                 tfar = tnext;
%                 tnext = (tnext+tprev)/2;
%             elseif errSign > 0
%                 tprev = tnext;
%                 tnext = (tnext+tfar)/2;
%             end
%             
%             if tnext<tprev
%                 fprintf('NOT CONVERGING AT ANGLE %d\n',tprev)
%                 break;  
%             end
%             
%             xn = a*cos(tnext);
%             yn = b*sin(tnext);
%             errSign = ( ds - sqrt( (xn-Xprev)^2 + (yn-Yprev)^2 ) );
%             err = abs(errSign);
%         end
%         
%         X(n) = xn;     %Store X-value
%         Y(n) = yn;     %Store Y-value
%         
%         Xprev= xn;     %Save prev. x-value of bisection
%         Yprev= yn;     %Save prev. y-value of bisection
%         tprev = tnext; %Update previous angle
%         
%         n = n+1;       %Update counter
%     
%         %plot(xn,yn,'*'); hold on;
%       
%     end
%         
%     %Length of vector before taking out extra point 
%     N = length(X);
% 
%     if half == 1
%         %Get rid of extra point
%         X_C = X(1:N-1);
%         Y_C = Y(1:N-1);
%     else
%         X_C = X;
%         Y_C = Y;
%     end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Computes the Volume of the Heart Tube
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function Vol = Volume(d,L,R_o,R_i)
% 
% Vol = 2*d*L + pi*(R_o^2 - R_i^2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Puts Red Blood Cells in
% % Note: Finds (x,y)-centers for blood cells
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [x,y] = Put_RBCs(N,L,d,R_i,R_o)
% 
% %Assume at first r < d/2
% 
% R = (R_i + R_o ) / 2;
% 
% L_tube = 2*L + 2*pi*R;
% 
% ds = L_tube /(N);          %Distance between points
% 
% x = zeros(N,1);
% y = x;
% 
% for i=1:N
%     
%    z = (i-1)*ds;
%    
%    if mod(i,2) == 0
%        w = d/3;
%    else
%        w = 2*d/3; 
%    end
%    
%    if z < L
%        
%        x(i) = z;
%        y(i) = -R_i - w;
%        
%    elseif ( ( z >= L ) && ( z < (L + pi*R)) )
%        
%        zz = z - L;
%        x(i) = (R_i+w)*cos( zz/R + 3*pi/2 ) + L;
%        y(i) = (R_i+w)*sin( zz/R + 3*pi/2 );
%        
%    elseif ( ( z >= (L + pi*R) ) && ( z < (2*L + pi*R) )  )
% 
%        zz = z - ( L + pi*R );
%        x(i) = L - zz; 
%        y(i) = R_i + w ;
%        
%    else
%        
%        zz = z - (2*L +pi*R);
%        x(i) = (R_i+w)*cos( zz/R + pi/2 );
%        y(i) = (R_i+w)*sin( zz/R + pi/2 );
%        
%    end
%         
% end
% 
%    %Shift so symmetric about y-axis
%     for j=1:length(x)
%        x(j)= x(j) - L/2;
%     end
% 
%     %Place blood cells out of region where tube is contracted for
%     %peristalsis (need to move RBCs 2,3,4,5
%     
%     %UNCOMMENT FOR SINE WAVE PERISTALSIS
% %     dely = abs( y(6) - y(5) );
% %     
% %     y(2) = y(2) - 0.25*dely;
% %     
% %     x(3) = x(9);
% %     y(3) = y(8) + dely;
% %     
% %     x(4) = x(8);
% %     y(4) = y(8) - 2*dely;
% %     
% %     x(5) = ( x(3) + x(10) )/2;
% %     y(5) = y(10) - 1.5*dely;
% %     
% %     x(9) = ( x(9) + x(8) )/2;
% %     y(9) = y(9) + 0.5*dely;
% %     
% %     x(8) = ( x(8) + x(7) )/2;
% %     y(8) = y(8) + 0.5*dely;
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Makes the [Circular] Red Blood Cells
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [X_c,Y_c,Npts] = Make_Circle_Cells(Nc,Vol,Vf,x_m,y_m,ds)
% 
% Vc = Vol*Vf;  
% r = sqrt( Vc / (Nc*pi) ); 
% 
% theta = 0:ds:2*pi;
% 
% for i=1:Nc
%    
%     x = x_m(i,1);
%     y = y_m(i,1);
%     
%     for j=1:length(theta)
%         
%         X_c(i,j) = x + r*cos(theta(j));
%         Y_c(i,j) = y + r*sin(theta(j));
%         
%     end
% end
% 
% %Stores number of pts for each cell as vector
% Npts = length(theta)*ones(1,Nc); 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Makes the [Elliptical] Red Blood Cells
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [X_c,Y_c,Npts] = Make_Elliptical_Cells(Nc,V,Vf,L,e,x_m,y_m,ds)
% 
% 
% %Vf = volume fraction of cells
% %V  = volume of full tube
% %Vc = volume of each blood cell
% %Nc = specified # of blood cells
% %[x_m, y_m] = (x,y) coordinates of center-point
% %pts = number of discretized pts on boundary
% 
% Vc = (Vf*V)/Nc;
% 
% r2 = ( sqrt(1-e^2)*Vc / pi )^(1/2);
% r1 = r2 / sqrt(1-e^2);
% 
% %X_c = zeros(Nc,pts);
% %Y_c = X_c;
% 
% X_HT = L/2;
% Y_HT = 0;
% 
% %theta = 0:ds:2*pi;
% 
% for i=1:Nc
%    
%     x = x_m(i,1);
%     y = y_m(i,1);
%     
%     %To add rotation into cells if lying in rounding part of HT
%     if ( x > L/2 )    
%         m = -(x-X_HT)/(y-Y_HT);
%         angle =  atan(m);          
%     elseif (x < -L/2)                
%         m = -(x+X_HT)/(y-Y_HT);  %X_HT -> (-1)*X_HT since half-circle on left side
%         angle =  atan(m);                   
%     else
%         angle = 0.0;
%     end
%     
%     %ACTUALLY MAKE ELLIPTICAL BLOOD CELL GEOMETRY BELOW%
% 
%     %Make ellipse as though around origin
%     half = 0;
%     [xE,yE] = make_Elliptical_Geometry(r1,r2,ds,half);
%     
%     for j=1:length(xE) %pts %(# of pts on blood cell)
%         
%         X_o = xE(j); %X_o = r1*cos(theta(j));
%         Y_o = yE(j); %Y_o = r2*sin(theta(j));
%         
%         %Rotate as though about origin
%         X_R = X_o*cos(angle)-Y_o*sin(angle);
%         Y_R = X_o*sin(angle)+Y_o*cos(angle);
%         
%         %Not translate to new position
%         X_T = X_R + x;
%         Y_T = Y_R + y;
%         
%         %Actual points on boundary of ellipse
%         X_c(i,j)=X_T;
%         Y_c(i,j)=Y_T;
%     
%     end    
% 
%     %Saves # of pts for each cell in vector
%     Npts(i) = length(xE);
%     
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Makes the Biconcave_Cells Red Blood Cells
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [X_c,Y_c] = Make_Biconcave_Cells(Nc,V,Vf,L,x_m,y_m,pts)
% 
% 
% %Vf = volume fraction of cells
% %V  = volume of full tube
% %Vc = volume of each blood cell
% %Nc = specified # of blood cells
% %[x_m, y_m] = (x,y) coordinates of center-point
% %pts = number of discretized pts on boundary
% 
% Vc = (Vf*V)/Nc;
% 
% a = sqrt(Vc/0.322837);
% 
% %X_c = zeros(Nc,pts);
% %Y_c = X_c;
% 
% X_HT = L/2;
% Y_HT = 0;
% 
% theta = 0:2*pi/(0.5*pts):2*pi;
% 
% for i=1:Nc
%    
%     x = x_m(i,1);
%     y = y_m(i,1);
%     
%     %To add rotation into cells if lying in rounding part of HT
%     if ( x > L/2 )    
%         m = -(x-X_HT)/(y-Y_HT);
%         angle =  atan(m);          
%     elseif (x < -L/2)                
%         m = -(x+X_HT)/(y-Y_HT);  %X_HT -> (-1)*X_HT since half-circle on left side
%         angle =  atan(m);                   
%     else
%         angle = 0.0;
%     end
%     
%     %Actually make boundary of blood cells on RHS
%     for j=1:pts/2 
%         %Make biconcave cell as though around origin
%         %ang = 0;
%         ang = theta(j);   %Theta for polar coordinates
%         r = a/2*( 0.207 + 2.003*sin(ang)^2 - 1.123*sin(ang)^4 )*cos(ang);
%         X_o = r*cos(ang);
%         Y_o = r*sin(ang);
%         %pause();
%         
%         %[xBi, yBi] = make_Biconcave_Geometry(a,ds);
%         
%         %len = length(xBi);
%         %for j=1:len
%         
%         %X_o = xBi(j);
%         %Y_o = yBi(j);
%         
%         X_n = -X_o;
%         
%         %Rotate as though about origin
%         X_R = X_o*cos(angle)-Y_o*sin(angle);
%         Y_R = X_o*sin(angle)+Y_o*cos(angle);
%         
%         X_R2= X_n*cos(angle)-Y_o*sin(angle);
%         Y_R2= X_n*sin(angle)+Y_o*cos(angle);
%         
%         %Now translate to new position
%         X_T = X_R + x;
%         Y_T = Y_R + y;
%         X_T2= X_R2+ x;
%         Y_T2= Y_R2+ y;
%         
%         %Actual points on boundary of ellipse
%         X_c(i,j)=X_T;
%         Y_c(i,j)=Y_T;
%         X_c(i,j+pts/2) = X_T2;
%         Y_c(i,j+pts/2) = Y_T2;
%         
%     end
%  
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Construct on boundary: Makes Right/Left Arcs.
% % Note that the counting starts at the point (0,-b) or the top point.
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [X_C,Y_C] = make_Biconcave_Geometry(a,ds)
% 
% %a = semi-major (horizontal) axis of ellipse 
% %b = semi-minor (vertical) axis of ellipse
% %ds = lag. point spacing
% 
% n = 1;                     %counting for array
% Xprev = 0.011875878869741; %first x-value
% Yprev = 0.0;               %first y-values
% tprev = 0.0;               %initial angle 
% tol = ds / 10;             %error tolerance for root-finding algorithm
% maxy = max(a,a);           %To be used for guessing next pt.
% 
% X(n) = Xprev;  Y(n) = Yprev; %Store first-value
% 
% ang_end = pi; %Ending angle for bisection
% 
%     while tprev < ang_end 
% 
%         n = n+1    %Update counter
% 
%         tnext = tprev + ds/maxy; %Guess for next angle value
%         
%         if tprev < -pi/4
%             tfar =  0.0;%abs(2*tprev);      %Far guess
%         else
%             tfar = pi;
%         end
%             
%         %initiating guess for bisection-algorithm
%         r = a/2*( 0.207 + 2.003*sin(tnext)^2 - 1.123*sin(tnext)^4 )*cos(tnext);
%         %X_o = r*cos(tnext)
%         %Y_o = r*sin(tnext)
%         
%         xn = r*cos(tnext);
%         yn = r*sin(tnext);
%         errSign = ( ds - sqrt( (xn-Xprev)^2 + (yn-Yprev)^2 ) );
%         err = abs(errSign);
%         
%         %Bisection algorithm to make points equally spaced
%         while ( err > tol )
%         
%             if errSign < 0
%                 tfar = tnext;
%                 tnext = (tnext+tprev)/2;
%             elseif errSign > 0
%                 tprev = tnext;
%                 tnext = (tnext+tfar)/2;
%             end
%             
%             if tnext<tprev
%                 fprintf('NOT CONVERGING AT ANGLE %d\n',tprev)
%                 break;  
%             end
%             
%             r = a/2*( 0.207 + 2.003*sin(tnext)^2 - 1.123*sin(tnext)^4 )*cos(tnext);
%             xn = r*cos(tnext);
%             yn = r*sin(tnext);
%             errSign = ( ds - sqrt( (xn-Xprev)^2 + (yn-Yprev)^2 ) );
%             err = abs(errSign);
%             
% %             for j=1:length(X)
% %                xtest = X(j);
% %                ytest = Y(j);
% %                
% %                diffy = sqrt( (xn-xtest)^2 + (yn-ytest)^2 );
% %                if diffy < 1e-5
% %                    err = 0;
% %                end
% %                
% %             end
%         end
%         
%         X(n) = xn;     %Store X-value
%         Y(n) = yn;     %Store Y-value
%         
%         Xprev= xn;     %Save prev. x-value of bisection
%         Yprev= yn;     %Save prev. y-value of bisection
%         tprev = tnext; %Update previous angle
%             
%         tprev
%         
%         %plot(xn,yn,'*'); hold on;
%         %pause();
%       
%     end
%         
%     %Length of vector before taking out extra point 
%     N = length(X);
% 
%     if half == 1
%         %Get rid of extra point
%         X_C = X(1:N-1);
%         Y_C = Y(1:N-1);
%     else
%         X_C = X;
%         Y_C = Y;
%     end
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FUNCTION TO PLOT GEOMETRY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function plot_it_all(X_tube,Y_tube,x_m,y_m,X_c,Y_c,R_o,L,Nc,Vf,d)
% 
% plot(X_tube,Y_tube,'o'); hold on;
% plot(x_m,y_m,'g*'); hold on;
% for i=1:Nc
%    plot(X_c(i,:),Y_c(i,:),'r.'); hold on; 
% end
% axis([-R_o-3*L/4 R_o+3*L/4 -R_o-d R_o+d]);
% t1 = strcat('Volume Fraction = ',num2str(Vf));
% t2 = strcat('Tube diameter = ',num2str(d),{'  /  '});
% t3 = strcat('Num. RBCs = ',num2str(Nc),{'  /  '});
% str_title = strcat(t3,t2,t1);
% 
% title(str_title);
% xlabel('x');
% ylabel('y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to Plot Pts. to Text File for Bottom of Tube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Them_Vertices(X,str)

fileID = fopen(str,'w');
fprintf(fileID,'%1.16e\n',X);
fclose(fileID);