#include "update_target_point_positions.h"
#include <ibamr/IBTargetPointForceSpec.h>
#include <stdio.h>
#include <stdlib.h>

void update_target_point_positions(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt)
{
    
    const int finest_ln = hierarchy->getFinestLevelNumber();

    static const double pi = 4*atan(1);	
	
	double spring_hard = 25.0;	  // spring stiffness for the gaussian wave
	double spring_soft = 0.001;	  // spring stiffness for the rest of the tube
	
	//Info for Peristaltic Region stuff
	static const int numPts = 569;//Number of Points in Bottom of Tube, changed from 155 4/20
	// is this the number of points period or points that move for peristalsis?
	static double X[numPts];	  //Stores X-Values for Bottom of Tube
	static double yOut_1[numPts]; //Stores Y-Values for OUTER Region in PHASE 1
	static double yOut_2[numPts]; //Stores Y-Values for OUTER Region in PHASE 2
	static double yOut_3[numPts]; //Stores Y-Values for OUTER Region in PHASE 3
	static double yIn_1[numPts];  //Stores Y-Values for INNER Region in PHASE 1
	static double yIn_2[numPts];  //Stores Y-Values for INNER Region in PHASE 2
	static double yIn_3[numPts];  //Stores Y-Values for INNER Region in PHASE 3

    //
    // Storing all the Y-Values for Peristaltic Phase
    //
	FILE *fileID1, *fileID2,*fileID3, *fileID4,*fileID5, *fileID6, *fileID7;
    fileID1 = fopen("X.txt", "r");
    fileID2 = fopen("yOut_1.txt", "r");
	fileID3 = fopen("yOut_2.txt", "r");
    fileID4 = fopen("yOut_3.txt", "r");
    fileID5 = fopen("yIn_1.txt", "r");
    fileID6 = fopen("yIn_2.txt", "r");
    fileID7 = fopen("yIn_3.txt", "r");
	for(int i=0; i<numPts; i++){
		fscanf(fileID1,"%lf\n",&X[i]);
		fscanf(fileID2,"%lf\n",&yOut_1[i]);
		fscanf(fileID3,"%lf\n",&yOut_2[i]);
		fscanf(fileID4,"%lf\n",&yOut_3[i]);
		fscanf(fileID5,"%lf\n",&yIn_1[i]);
		fscanf(fileID6,"%lf\n",&yIn_2[i]);
		fscanf(fileID7,"%lf\n",&yIn_3[i]);
		
		
		//cout << yIn_2[i] << "\n";
	}
	
	cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
	
    fclose(fileID1);
	fclose(fileID2);
	fclose(fileID3);
	fclose(fileID4);    
	fclose(fileID5);
	fclose(fileID6);
	fclose(fileID7);
		
	double period = 0.5;          //period of peristalsis
	double t1 = 0.05*period;	  //time to make gaussian wave
	double t2= period - 2*t1;	  //time for wave to move
	double t3 = t1;	              //time to reduce gaussian wave
	
	
	int nO_1st = 0;      //First Pt. of Peristaltic Region on OUTER
	int nO_2nd = 568;  //Last Pt. of Peristaltic Region on OUTER
	int nI_1st = 569;   //First Pt. of Peristaltic Region on INNER
	int nI_2nd = 1137; //Last Pt. of Peristaltic Region on INNER
	int ppm = 4096;		//Number of points per meter (1/ds2 in generate mesh)
	//Need to check these. Note: V5 heart tube vertex files write inner (top) first, then outer (top). Guessing that there is a region
	// before and after the peristalsis region of each tube that is 8 points long. length of peristaltic region is 104.
	//NOTE: -1's are bc counting starts at 0 in arrays in C++ 
	
	int p0 = 31;	// Left PT of wave at initial position
	int p1 = 121;	// Right Pt of wave at initial position
	int pC = round((p0+p1)/2);	// Pt closest to center at initial position
			
	double x0 = X[p0];		      //Left- position of Wave at INITIAL POSITION, changed from 31 4/20
	double x1 = X[p1];		      //Right- position of Wave at INITIAL POSITION
	// change these to scale? with X = 616, x0 = X[31], x1 = X[113]
	double xC = (x0+x1)/2;		  //Center-Pt of Wave at INITIAL POSITION
	double tt = fmod(current_time,period); //Current time in simulation (remainder of time/period for phases)
	double c = (-2*xC)/t2;        //Wave-Speed (minus sign bc xC is negative)
	double tp = c*(tt-t1);        // Wave-Speed * time
	double xCn = xC + tp;		  //Center-Pt of Wave at time, t
	double x0n = x0 + tp;		  //Left-Pt of Wave at time, t
	double x1n = x1 + tp;		  //Right-Pt of Wave at time, t
	
	int pCn = pC + round(tp*ppm);		// Position of center of wave at time t
	int p0n = p0 + round(tp*ppm);		// Position of left end of wave at time t
	int p1n = p1 + round(tp*ppm);		// Position of right of wave at time t
		
	
	double g1;		//Interpolation Function between Phase 1 and 2
	double g3;		//Interpolation Function between Phase 3 and 4
	
	double x;				  //x-Pt specified (rolls over each lag-pt)
	double w = 0.025;           //Width of Gaussian Wave
	double A_tilde = 1638400.0*2;   //"Fixed" Amplitude for Gaussian Wave
	double R_o = 0.2;	      //OUTER Radius
	double R_i = 0.1;        //INNER Radius

    //
    // Find out the Lagrangian index ranges.
    //
    //const std::pair<int,int>& lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);


    //
    // Get the LMesh (which we assume to be associated with the finest level of
    // the patch hierarchy).  Note that we currently need to update both "local"
    // and "ghost" node data.
    //
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());


    //
    // Update the target point positions in their associated target point force specs.
    //
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
 
 
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        LNode* node_idx = *it;
        IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec == NULL) continue;  // skip to next node


        // Here we update the position of the target point.
        //
        // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
        //        X_target     is the target position of the target point
        //        X_target[0]  is the x component of the target position
        //        X_target[1]  is the y component of the target position
        //        X_target[2]  is the z component of the target position (for a 3D simulation)
        //


        const int lag_idx = node_idx->getLagrangianIndex();
		

        //
        //Depending on the version of IBAMR, you need to select one of the ways of accessing target point positions
        //
        //FOR KD MODULE:
        Point& X_target = force_spec->getTargetPointPosition();
        double X_spring = force_spec->getStiffness();
        //FOR NEMOs / KD (NOT MODULE)
        //TinyVector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
        //
        //OLD:
        //IBTK::Vector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
		

		x = X_target[0]; //x-pt for specific target point associated with lag_idx
        
		//cout << "tt = " << tt << "\n";
		//cout << "t1 = " << t1 << "\n";
		
		
		if (tt<=t1) {  //TIME FOR PHASE 1	-> Make Gaussian Wave from straight tube
			
				//g1 = (1.0/t1)*current_time + 0.0;
				
				if ((lag_idx>=nO_1st)&&(lag_idx<=nO_2nd)) {
					
					
					g1 = (1.0/t1)*tt + 0.0;
					
					X_target[1] = yOut_1[lag_idx] + g1*( yOut_2[lag_idx] - yOut_1[lag_idx] );
					
					//cout << "tt = " << tt << "\n";
					//cout << "g1 = " << g1 << "\n";
					//cout << "lag_idx =" << lag_idx << "\n";
					//cout << "yOut_1 = " << yOut_1[lag_idx] << "\n";
					//cout << "yOut_2 = " << yOut_2[lag_idx] << "\n";
					//cout << "X_target[1] = " << X_target[1] << "\n\n\n\n";
					
					if ((lag_idx>=p0) && (lag_idx<=p1)) {
		
						X_spring = spring_hard;
		
					} else {
		
						X_spring = spring_soft;
		
					} 
			
				} else if ((lag_idx>=nI_1st)&&(lag_idx<=nI_2nd)) {
			
					g1 = (1.0/t1)*tt + 0.0;
					
					X_target[1] = yIn_1[lag_idx-nI_1st] + g1*( yIn_2[lag_idx-nI_1st] - yIn_1[lag_idx-nI_1st] );
					
					//cout << "tt = " << tt << "\n";
					//cout << "g1 = " << g1 << "\n";
					//cout << "lag_idx =" << lag_idx << "\n";
					//cout << "yIn_1 = " << yIn_1[lag_idx-nI_1st] << "\n";
					//cout << "yIn_2 = " << yIn_2[lag_idx-nI_1st] << "\n";
					//cout << "X_target[1] = " << X_target[1] << "\n\n\n\n";

					if ((lag_idx>=(p0+nI_1st)) && (lag_idx<=(p1+nI_1st))) {
					
						X_spring = spring_hard;
					
					} else {
						
						X_spring = spring_soft;
		
					} 
					
				}
			
				
			        //END TIME FOR PHASE 1 -> Gaussian Wave is Made 
			
		} else if ( (tt>t1) && (tt<=(t1+t2)) ) {
					
					//TIME FOR PHASE 2 -> Translate Gaussian Wave
			
				
				if (x < x0n) {
				
					if ((lag_idx>=nO_1st)&&(lag_idx<=nO_2nd)) {
						X_target[1] = -R_o;
						X_spring = spring_soft;
					} else if ((lag_idx>=nI_1st)&&(lag_idx<=nI_2nd)) {
						X_target[1] = -R_i;
						X_spring = spring_soft;
					}
				
				} else if (x > x1n) {
				
					if ((lag_idx>=nO_1st)&&(lag_idx<=nO_2nd)) {
						X_target[1] = -R_o;
						X_spring = spring_soft;
					} else if ((lag_idx>=nI_1st)&&(lag_idx<=nI_2nd)) {
						X_target[1] = -R_i;
						X_spring = spring_soft;
					}
				
				} else {    
				
					if ((lag_idx>=nO_1st)&&(lag_idx<=nO_2nd)) {
					
						X_target[1] =    A_tilde*pow(x-x0n,2)*pow(x1n-x,2)*exp( -1*pow(x-xCn,2) / pow(w/2.0,2)  ) - R_o;
						X_spring = spring_hard;
					
					} else if ((lag_idx>=nI_1st)&&(lag_idx<=nI_2nd)) {
					
						X_target[1] = -1*A_tilde*pow(x-x0n,2)*pow(x1n-x,2)*exp( -1*pow(x-xCn,2) / pow(w/2.0,2)  ) - R_i;
						X_spring = spring_hard;
					
					}
				
				}
			
					//END TIME FOR PHASE 2 -> Gaussian Wave Finished Moving
		
		} else if ( (tt > (t1+t2)) && ( tt <= (t1+t2+t3)) ) {
					
					//TIME FOR PHASE 3 -> Crunch Gaussian Wave to straight tube 
						
				if ((lag_idx>=nO_1st)&&(lag_idx<=nO_2nd)) {
				
					g3 = (1.0/t3)*tt - (t1+t2)/t3;
					X_target[1] = yOut_3[lag_idx] + g3*( yOut_1[lag_idx] - yOut_3[lag_idx] );
					
					if ((lag_idx>=p0) && (lag_idx<=p1)) {
		
						X_spring = spring_hard;
		
					} else {
		
						X_spring = spring_soft;
		
					}
				
				} else if ((lag_idx>=nI_1st)&&(lag_idx<=nI_2nd)) {
				
					g3 = (1.0/t3)*tt - (t1+t2)/t3;
					X_target[1] = yIn_3[lag_idx-nI_1st] +  g3*( yIn_1[lag_idx-nI_1st] - yIn_3[lag_idx-nI_1st] );
				
					if ((lag_idx>=p0) && (lag_idx<=p1)) {
		
						X_spring = spring_hard;
		
					} else {
		
						X_spring = spring_soft;
		
					}
				}
					
					//END TIME FOR PHASE 3 -> Now back to straight tube
			
		}

		
		//cout << lag_idx << "\n"

        
    
    } //Ends for-loop over all lag_pts
    
    return;

}// update_target_point_positions
