function write_Input_Files(NTvec,NCvec,Xtube,Ytube,X_c,Y_c)

%NOTE: Didn't write beams to file, yet <-- Not sure if needed
%NOTE: Only setup for elliptical/circular RBCs -> NOT biconcave yet!
% MW: 03/29 commented out RBC

%NTvec Entries
N1 = NTvec(1);  % # of pts. before wave region starts on bottom (still flat portion)
N2 = NTvec(2);  % # of pts. in peristalsis region on bottom (WAVE PART ONLY)
N3 = NTvec(3);  % # of pts. in bottom part of tube (flat+peri+flat) up to rounded part
N4 = NTvec(4);  % # of pts leading up to straight-top on OUTER of tube
N5 = NTvec(5);  % # of pts on straight section of tube on TOP
N6 = NTvec(6);  % # of pts leading up to peristaltic-bottom portion on INSIDE of tube
N7 = NTvec(7);  % # of pts leading up to straight-top portion on INSIDE of tube
N8 = NTvec(8);  % # of pts on entire tube (OUTER + INNER)

NTvec'

%NCvec Entries
%NCvec(i) = # of lag-pts on ith RBC


%Heart-Tube Vertices
print_HT_Vertices(NTvec,Xtube,Ytube);

%Heart-Tube Springs%
%Springs_Side_Force  /  Springs_On_Top  /   Springs_Btwn_Outer_&_Inner
sF_connect = 1e8;       sF_top = 2.5e7;         sF_btwn = 5e3;
print_HT_Springs(NTvec,Xtube,Ytube,sF_connect,sF_top,sF_btwn);

%Heart-Tube Beams%
bF_top = 1e8;
print_HT_Beams(NTvec,Xtube,Ytube,bF_top)

%Heart-Tube Target Pts%
tForce = 1e6; t_Force_Rnd = 1e6;
print_HT_Target_Pts(NTvec,tForce,t_Force_Rnd);

% %RBC Vertices
% print_RBC_Vertices(NCvec,X_c,Y_c);
% 
% %RBC Springs
% %RBC Spring-Force
% sF_RBC = 1e7;   sF_inside = 1e7;
% write_RBC_springs(NCvec,X_c,Y_c,sF_RBC,sF_inside);
% 
% %RBC Beams
% %RBC Beam-Force
% %bF_RBC = 1e6;
% %write_RBC_beams(NCvec,X_c,Y_c,bF_RBC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Write out the VERTEX information for Heart Tube %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_HT_Vertices(NTvec,Xtube,Ytube)

N6 = NTvec(6);  % # of pts leading up to straight-bottom portion on INSIDE of tube
N8 = NTvec(8);  % # of pts on entire tube (OUTER + INNER)
Ni = N8-N6;


    vertex_fid = fopen(['heart_tube' '.vertex'], 'w');

    fprintf(vertex_fid, '%d\n', N8 );

    %Prints Outer Tube Vertices
    for s = 0:N6-1
    
        X(1) = Xtube(s+1);
        X(2) = Ytube(s+1);
        fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
    end

    %Prints Inner Tube
    for s = 0:Ni-1
       
        X(1) = Xtube(N6+s+1);
        X(2) = Ytube(N6+s+1);     
        fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
    end

    fclose(vertex_fid); %%%ENDS HEART_TUBE VERTICES
    
    
    %plot(Xtube(1:N6),Ytube(1:N6),'m*'); hold on;
    %plot(Xtube(N6+1:N8),Ytube(N6+1:N8),'k*'); hold on;
    %pause()

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Write out the SPRING information for Heart Tube  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_HT_Springs(NTvec,Xtube,Ytube,sF_connect,sF_top,sF_topbottom)

    %NTvec Entries
    N1 = NTvec(1);  % # of pts. before wave region starts on bottom (still flat portion)
    N2 = NTvec(2);  % # of pts. in peristalsis region on bottom (WAVE PART ONLY)
    N3 = NTvec(3);  % # of pts. in bottom part of tube (flat+peri+flat) up to rounded part
    N4 = NTvec(4);  % # of pts leading up to straight-top on OUTER of tube
    N5 = NTvec(5);  % # of pts on straight section of tube on TOP
    N6 = NTvec(6);  % # of pts leading up to flat-peristalstic-flat-bottom portion on INSIDE of tube
    N7 = NTvec(7);  % # of pts leading up to straight-top portion on INSIDE of tube
    N8 = NTvec(8);  % # of pts on entire tube (OUTER + INNER)

    N9 = N3-N2-N1;   % # of pts. after wave-region on bottom before rounded part
    
    %Nc = 2 + 1 + 2 + 1; % # of connecting springs (2 extra on bottom, 1 on top for both INNER + OUTER)
    %N_Springs = 2*N1 + 2*N9 + 2*N5 + Nc;
    
    Nc = 2*(1 + 1);        %Extra connecting springs on top (2 is to add one extra spring on left)
    N_Springs = 3*N5 + Nc; %Along Outer Top, Along Inner Top, Between Inner/Outer Top 
    
    spring_fid = fopen(['heart_tube' '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N_Springs ); 
    
    sF_1 = sF_connect;
    sF_2 = sF_top;
    sF_3 = sF_topbottom;
    
    %d_eq = equilibrium distance between pts

    %SPRINGS BETWEEN OUTER TUBE VERTICES ON BOTTOM-LEFT
%     for s = 0:N1
%         if s == 0 
%             Xc = Xtube(N6);     %bc indexing starts at 1
%             Yc = Ytube(N6);
%             Xnext = Xtube(s+1); %bc indexing starts at 1
%             Ynext = Ytube(s+1);
%             ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%             d_eq = ds;
%             %spring_force = kappa_spring*ds/(ds^2)
%             fprintf(spring_fid, '%d %d %1.16e %1.16e\n', N6-1, s, sF_1, d_eq);
%         else
%             Xc = Xtube(s);
%             Yc = Ytube(s);
%             Xnext = Xtube(s+1);
%             Ynext = Ytube(s+1);
%             ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%             d_eq = ds;
%             %spring_force = kappa_spring*ds/(ds^2)
%             fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s-1, s, sF_1, d_eq);  
%         end      
%     end
    
    %SPRINGS BETWEEN OUTER TUBE VERTICES ON BOTTOM-RIGHT
%     for s = 0:N9
%             Xc = Xtube(N1+N2+s);      %bc indexing starts at 1
%             Yc = Ytube(N1+N2+s);
%             Xnext = Xtube(N1+N2+s+1); %bc indexing starts at 1 
%             Ynext = Ytube(N1+N2+s+1);
%             ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%             d_eq = ds;
%             %spring_force = kappa_spring*ds/(ds^2)
%             fprintf(spring_fid, '%d %d %1.16e %1.16e\n', N1+N2+s-1,N1+N2+s, sF_1, d_eq); 
%     end
    
    %Springs between OUTER vertices on TOP
    for s = 0:N5+1
            Xc = Xtube(N4+s);      %bc indexing starts at 1
            Yc = Ytube(N4+s);
            Xnext = Xtube(N4+s+1); %bc indexing starts at 1
            Ynext = Ytube(N4+s+1);
            ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
            d_eq = ds;
            %spring_force = kappa_spring*ds/(ds^2)
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n',N4+s-1,N4+s, sF_2, d_eq);
    end
    

    
%     %SPRINGS BETWEEN INNER TUBE VERTICES ON BOTTOM
%     for s = 0:N1
%         if s == 0 
%             Xc = Xtube(N8);       %bc indexing starts at 1
%             Yc = Ytube(N8);
%             Xnext = Xtube(N6+s+1);%bc indexing starts at 1
%             Ynext = Ytube(N6+s+1);
%             ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%             d_eq = ds;
%             %spring_force = kappa_spring*ds/(ds^2)
%             fprintf(spring_fid, '%d %d %1.16e %1.16e\n', N8-1,N6+s, sF_1, d_eq);
%         else
%             Xc = Xtube(N6+s);     %bc indexing starts at 1
%             Yc = Ytube(N6+s);
%             Xnext = Xtube(N6+s+1);%bc indexing starts at 1
%             Ynext = Ytube(N6+s+1);
%             ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%             d_eq = ds;
%             %spring_force = kappa_spring*ds/(ds^2)
%             fprintf(spring_fid, '%d %d %1.16e %1.16e\n',N6+s-1,N6+s, sF_1, d_eq);  
%         end      
%     end
%     
%     %SPRINGS BETWEEN INNER TUBE VERTICES ON BOTTOM-RIGHT
%     for s = 0:N9
%             Xc = Xtube(N6+N1+N2+s);     %bc indexing starts at 1
%             Yc = Ytube(N6+N1+N2+s);
%             Xnext = Xtube(N6+N1+N2+s+1);%bc indexing starts at 1
%             Ynext = Ytube(N6+N1+N2+s+1);
%             ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%             d_eq = ds;
%             %spring_force = kappa_spring*ds/(ds^2)
%             fprintf(spring_fid, '%d %d %1.16e %1.16e\n', N6+N1+N2+s-1, N6+N1+N2+s, sF_1, d_eq); 
%     end
    
    %Springs between INNER vertices on TOP
    for s = 0:N5+1
            Xc = Xtube(N7+s);      %bc indexing starts at 1
            Yc = Ytube(N7+s);
            Xnext = Xtube(N7+s+1); %bc indexing starts at 1 
            Ynext = Ytube(N7+s+1);
            ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
            d_eq = ds;
            %spring_force = kappa_spring*ds/(ds^2)
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n',N7+s-1,N7+s, sF_2, d_eq);
    end
    
    
    %Springs between OUTER AND INNER vertices on TOP
    for s = 0:N5-1
            Xout = Xtube(N4+s+1);   %bc indexing starts at 1
            Yout = Ytube(N4+s+1);
            Xin = Xtube(N7+s+1);    %bc indexing starts at 1
            Yin = Ytube(N7+s+1);
            ds = sqrt( (Xout-Xin)^2 + (Yout-Yin)^2 );
            d_eq = ds;
            %spring_force = kappa_spring*ds/(ds^2)
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n',N4+s,N7+s, sF_3, d_eq);
    end

    
    fclose(spring_fid);  %%%%END SPRING PORTION FOR HEART TUBE

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Write out the BEAM information for HEARTTUBE %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_HT_Beams(NTvec,Xtube,Ytube,bF_top)
    
    %NTvec Entries
    N1 = NTvec(1);  % # of pts. before wave region starts on bottom (still flat portion)
    N2 = NTvec(2);  % # of pts. in peristalsis region on bottom (WAVE PART ONLY)
    N3 = NTvec(3);  % # of pts. in bottom part of tube (flat+peri+flat) up to rounded part
    N4 = NTvec(4);  % # of pts leading up to straight-top on OUTER of tube
    N5 = NTvec(5);  % # of pts on straight section of tube on TOP
    N6 = NTvec(6);  % # of pts leading up to flat-peristalstic-flat-bottom portion on INSIDE of tube
    N7 = NTvec(7);  % # of pts leading up to straight-top portion on INSIDE of tube
    N8 = NTvec(8);  % # of pts on entire tube (OUTER + INNER)

    N9 = N3-N2-N1;   % # of pts. after wave-region on bottom before rounded part
    
    Nc = 2 + 2; %Extra connecting beams on top
    
    %N_Springs = 2*N1 + 2*N9 + 2*N5 + Nc;
    N_beams = 2*N5 + Nc;
    
    beam_fid = fopen(['heart_tube' '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N_beams ); 

    beam_force = bF_top;
    
    %BEAMS BETWEEN OUTER TUBE TOP VERTICES
    for s = 0:N5+1     %Loops over top vertices in outer boundary
        
        C1=0.0; %?
        C2=0.0; %?
        
            %Xc = x_out(s+1);
            %Yc = y_out(s+1);
            %Xnext = x_out(s+2);
            %Ynext = y_out(s+2);
            %ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
            %beam_force = kappa_beam*ds/(ds^4)
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n',N7+s-2,N7+s-1,N7+s, beam_force , C1, C2);
    end
    
    %BEAMS BETWEEN INNER TUBE TOP VERTICES
    for s = 0:N5+1     %Loops over top vertices in inner boundary
        
        C1=0.0; %?
        C2=0.0; %?
        
            %Xc = x_out(s+1);
            %Yc = y_out(s+1);
            %Xnext = x_out(s+2);
            %Ynext = y_out(s+2);
            %ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
            %beam_force = kappa_beam*ds/(ds^4)
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n',N4+s-2,N4+s-1,N4+s, beam_force , C1, C2);
    end
    
    
    
    fclose(beam_fid);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Write out the TARGET information for HEARTTUBE %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_HT_Target_Pts(NTvec,tForce,tForceRnd)

    %NTvec Entries
    N1 = NTvec(1);  % # of pts. BEFORE wave region starts on bottom (still flat portion)
    N2 = NTvec(2);  % # of pts. in peristalsis region on bottom (WAVE PART ONLY)
    N3 = NTvec(3);  % # of pts. in bottom part of tube (flat+peri+flat) up to rounded part
    N4 = NTvec(4);  % # of pts leading up to straight-top on OUTER of tube
    N5 = NTvec(5);  % # of pts on straight section of tube on TOP
    N6 = NTvec(6);  % # of pts leading up to straight-bottom portion on INSIDE of tube
    N7 = NTvec(7);  % # of pts leading up to straight-top portion on INSIDE of tube
    N8 = NTvec(8);  % # of pts on entire tube (OUTER + INNER)
    
    NR_o = N4-N3;     % # of pts. in OUTER rounded region (one side)
    NR_i = N7-(N6+N3);% # of pts. in INNER rounded region (one side)
    
    target_force = tForce;
    target_force2= tForceRnd;
    
 
    %Write out the target point information
    target_fid = fopen(['heart_tube' '.target'], 'w');
    
    fprintf(target_fid, '%d\n', 3*N3 + 2*NR_i + 2*NR_o );  %NOTE: num_o should equal num_i

    %Target Pts. For Bottom OUTER and OUTER Right side Curve
    for s = 0:N4-1
            %fprintf(target_fid, '%d %1.16e\n', s, kappa_target*ds/(ds^2));
            fprintf(target_fid, '%d %1.16e\n', s, target_force);
    end
    
    %Target Pts. for OUTER Left Side Curve, Bottom INNER, and Inner Right Side
    %for s = 0:N7-(N4+N5)-1
    for s = 0:NR_o+(N8-N6)-1
            %fprintf(target_fid, '%d %1.16e\n', s, kappa_target*ds/(ds^2));
            fprintf(target_fid, '%d %1.16e\n', N4+N5 + s, target_force);
    end
    
    %Target Pts. for INNER Left Side of Tube
    %for s = 0:NR_i-1
            %fprintf(target_fid, '%d %1.16e\n', s, kappa_target*ds/(ds^2));
    %        fprintf(target_fid, '%d %1.16e\n', N7+N5+s, target_force2);
    %end
    

    fclose(target_fid);
    
    %N3s = N2 + NR_o + NL_o;
    
    %N(1) = 0;     % START Location of Target Pts on OUTER on BOTTOM for peristalsis
    %N(2) = N2-1;  % END Location of Target Pts on OUTER on BOTTOM for peristalsis
    %N(3) = N3s;   % START Location of Target Pts on INNER on BOTTOM for peristalsis
    %N(4) = N3s+N2;% END Location of Target Pts on INNER on BOTTOM for peristalsis
    
    %N';
    
    return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% Write out the VERTEX information for RBCs  %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function print_RBC_Vertices(NCvec,X_c,Y_c)
% 
%     vertex_fid = fopen(['blood_cells' '.vertex'], 'w');
% 
%     Nc = length(X_c(:,1));
%     
%     %Counts total number of pts for RBCs (Total Vertices Needed...in case RBC aren't same size)
%     pts_tot = 0; 
%     for i=1:Nc
%         pts_tot = pts_tot + NCvec(i);
%     end
%     
%     fprintf(vertex_fid, '%d\n', pts_tot );
% 
%     for s = 0:Nc-1   %Loops over each RBC 
%         
%         ss = s+1;
%         pts = NCvec(ss); %Gives number of pts. for specific RBC
%         
%         for i=1:pts  %Loops over # pts. in each RBC
%             X(1) = X_c(ss,i);
%             X(2) = Y_c(ss,i);
%             fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
%         end 
%     end
% 
%     fclose(vertex_fid); %%ENDS BLOOD CELL VERTICES
%         
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% Write out the SPRING information for RBCs %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function write_RBC_springs(NCvec,X_c,Y_c,sF_RBC,sF_inside)
% 
%     spring_fid = fopen(['blood_cells' '.spring'], 'w');
%     
%     Nc = length(X_c(:,1));
%     
%     %Counts total number of pts for RBCs (Total Vertices Needed...in case RBC aren't same size)
%     pts_tot = 0; 
%     for i=1:Nc
%         pts_tot = pts_tot + NCvec(i);
%     end
%     
%     %Counts total number of springs INSIDE for RBCs (connecting 'side-to-side')
%     count = 0;
%     for i = 1:Nc                       %LOOPS OVER # OF BLOOD CELLS
%         pts = NCvec(i);                %Gives # of pts for specific RBC
%         count = count + floor(pts/2);
%     end
%     fprintf(spring_fid, '%d\n', pts_tot+count );
%     
%     spring_force = sF_RBC;
%     spring_force_inside = sF_inside;
%     
%     pts_so_far = 0;                    %Initialize                   
%     for i = 1:Nc                       %LOOPS OVER # OF BLOOD CELLS
%         pts = NCvec(i);                %Gives # of pts for specific RBC
%         for s=0:pts-1                  %LOOPS OVER PTS. AROUND CELLS
%             ss = s+1;
%             if s < pts-1 
%                 Xc = X_c(i,ss);
%                 Yc = Y_c(i,ss);
%                 Xnext = X_c(i,ss+1);
%                 Ynext = Y_c(i,ss+1);
%                 ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%                 s_label = s + pts_so_far;       %ORDERS THE SPRINGS
%                 s_label_next = s_label + 1;
%                 %fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s_label, s_label_next, k_spring_cell*ds/(ds^2), ds);
%                 fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s_label, s_label_next, spring_force, ds);
% 
%             else
%                 Xc = X_c(i,pts);
%                 Yc = Y_c(i,pts);
%                 Xnext = X_c(i,1);
%                 Ynext = Y_c(i,1);
%                 ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );   
%                 s_label = (pts-1) + pts_so_far;  %ORDERS THE SPRINGS
%                 s_label_next = pts_so_far;
%                 %fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s_label, s_label_next, k_spring_cell*ds/(ds^2), ds);
%                 fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s_label, s_label_next, spring_force, ds);  
%             end 
%         end
%         
%         for s=0:floor(pts/2)-1               %LOOPS OVER PTS. AROUND CELLS
%                 s1 = s + pts_so_far;         %ORDERS THE SPRINGS
%                 s2 = s1 + floor(pts/2) + 1;  %SPRING TO CONNECT '"ON OPPOSITE SIDE'"
%                 ss = s+1;
%                 Xc = X_c(i,ss+1);                %+1 bc indexing starts at 1 in Matlab
%                 Yc = Y_c(i,ss+1);                %+1 bc indexing starts at 1 in Matlab
%                 Xnext = X_c(i,ss+1+floor(pts/2));%Connects to spring on other side
%                 Ynext = Y_c(i,ss+1+floor(pts/2));%Connects to spring on other side
%                 ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%                 fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, spring_force_inside, ds);
%         end
%         
%         pts_so_far = pts_so_far + pts; %Number of RBC pts adding up
% 
%     end
% 
%     fclose(spring_fid); %%ENDS BLOOD CELL SPRINGS
%     
%     
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% Write out the BEAM information for BLOOD CELLS %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% function write_RBC_beams(NCvec,X_c,Y_c,sF_RBC)
%     
%     beam_fid = fopen(['blood_cells' '.beam'], 'w');
%     
%     Nc = length(X_c(:,1));
%     
%     %Counts total number of pts for RBCs (Total Vertices Needed...in case RBC aren't same size)
%     pts_tot = 0; 
%     for i=1:Nc
%         pts_tot = pts_tot + NCvec(i);
%     end
% 
%     fprintf(beam_fid, '%d\n', pts_tot );
%     
%     beam_force = sF_RBC;
%     
%     C1=0;
%     C2=0;
% 
%     pts_so_far = 0;                    %Initialize                   
%     for i = 1:Nc                       %LOOPS OVER # OF BLOOD CELLS
%         pts = NCvec(i);                %Gives # of pts for specific RBC
%         for s=0:pts-1                  %LOOPS OVER PTS. AROUND CELLS
%             ss = s+1;
%             if s < pts-2 
%                 Xc = X_c(i,ss);
%                 Yc = Y_c(i,ss);
%                 Xnext = X_c(i,ss+1);
%                 Ynext = Y_c(i,ss+1);
%                 ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%                 s_lab = s + pts_so_far;       %ORDERS THE SPRINGS
%                 s_lab_p1 = s_lab + 1;
%                 s_lab_p2 = s_lab + 2;
%                 %fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s_l, s_l+1, s_l+2, k_beam_cell*ds/(ds^4), C1, C2);
%                 fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s_lab, s_lab_p1, s_lab_p2, beam_force, C1, C2);
% 
%             elseif s == pts-2
%                     Xc = X_c(i,ss);
%                     Yc = Y_c(i,ss);
%                     Xnext = X_c(i,ss+1);
%                     Ynext = Y_c(i,ss+1);
%                     ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%                     s_lab = s + pts_so_far;       %ORDERS THE SPRINGS
%                     s_lab_p1 = s_lab + 1;
%                     s_lab_p2 = pts_so_far;
%                     %fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s_l, s_l+1, s_first, k_beam_cell*ds/(ds^4), C1, C2);
%                     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s_lab, s_lab_p1, s_lab_p2, beam_force, C1, C2);
% 
%             else
%                     %( s == pts-1 )
%                     %Xc = X_c(i,ss);
%                     %Yc = Y_c(i,ss);
%                     %Xnext = X_c(i,ss+1);
%                     %Ynext = Y_c(i,ss+1);
%                     %ds = sqrt( (Xc-Xnext)^2 + (Yc-Ynext)^2 );
%                     s_lab = s + pts_so_far;       %ORDERS THE SPRINGS
%                     s_lab_p1 = pts_so_far;
%                     s_lab_p2 = pts_so_far+1;
%                     %fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s_l, s_l+1, s_first, k_beam_cell*ds/(ds^4), C1, C2);
%                     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s_lab, s_lab_p1, s_lab_p2, beam_force, C1, C2);
%             end   
%         end
%         
%         pts_so_far = pts_so_far + pts; %Number of RBC pts adding up
%     
%     end
%     
%     fclose(beam_fid);
    