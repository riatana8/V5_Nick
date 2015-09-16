#include "update_springs.h"
#include <ibamr/IBSpringForceSpec.h>

void
update_springs(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt)

{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    static const double pi = 4*atan(1);
    static const double L1 = 1; // length of computational domain (meters)
    static const int N1 = 1024; // number of cartesian grid meshwidths at the finest level of the AMR grid
	
	static const int numPts = 599;//Number of Points in Bottom of Tube, changed from 155 4/20
	// is this the number of points period or points that move for peristalsis?
	static double X[numPts];	  //Stores X-Values for Bottom of Tube

	FILE *fileID1;
	fileID1 = fopen("X.txt", "r");
	for(int i=0; i<numPts; i++){
		fscanf(fileID1,"%lf\n",&X[i]);
	}
	
	cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
	
    fclose(fileID1);


	double period = 1.0;          //period of peristalsis
	double t1 = 0.05*period;	  //time to make gaussian wave
	double t2= period - 2*t1;	  //time for wave to move
	double t3 = t1;	              //time to reduce gaussian wave
	double spring_hard = 0.1;	  // spring stiffness for the gaussian wave
	double spring_soft = 0.001;	  // spring stiffness for the rest of the tube
	
	int nO_1st = 0;      //First Pt. of Peristaltic Region on OUTER
	int nO_2nd = 598;  //Last Pt. of Peristaltic Region on OUTER
	int nI_1st = 599;   //First Pt. of Peristaltic Region on INNER
	int nI_2nd = 1198; //Last Pt. of Peristaltic Region on INNER
	//Need to check these. Note: V5 heart tube vertex files write inner (top) first, then outer (top). Guessing that there is a region
	// before and after the peristalsis region of each tube that is 8 points long. length of peristaltic region is 104.
	//NOTE: -1's are bc counting starts at 0 in arrays in C++ 
	
	double x0 = X[31];		      //Left-Pt of Wave at INITIAL POSITION, changed from 31 4/20
	double x1 = X[91];		      //Right-Pt of Wave at INITIAL POSITION
	// change these to scale? with X = 616, x0 = X[31], x1 = X[113]
	double xC = (x0+x1)/2;		  //Center-Pt of Wave at INITIAL POSITION
	double tt = fmod(current_time,period); //Current time in simulation (remainder of time/period for phases)
	double c = (-2*xC)/t2;        //Wave-Speed (minus sign bc xC is negative)
	double tp = c*(tt-t1);        // Wave-Speed * time
	double xCn = xC + tp;		  //Center-Pt of Wave at time, t
	double x0n = x0 + tp;		  //Left-Pt of Wave at time, t
	double x1n = x1 + tp;		  //Right-Pt of Wave at time, t
	
	double g1;		//Interpolation Function between Phase 1 and 2
	double g3;		//Interpolation Function between Phase 3 and 4
	
	double x;				  //x-Pt specified (rolls over each lag-pt)


    // Find out the Lagrangian index ranges.
    const std::pair<int,int>& plate2d_left_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);
    //const std::pair<int,int>& plate2d_rght_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(1, finest_ln);

    // Get the LMesh (which we assume to be associated with the finest level of
    // the patch hierarchy).  Note that we currently need to update both "local"
    // and "ghost" node data.
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

    // Update the spring lengths in their associated spring specs.
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        LNode* node_idx = *it;
        IBSpringForceSpec* spring_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
		
        if (spring_spec == NULL) continue;  // skip to next node

        // Here we update the resting length of the spring
        //
        // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
        //        resting_length    	is the resting length of the current Lagrangian point that is the "master index'
		//		  						Since there may be more than one spring associated with each point, it's a vector.
        //        resting_length[0]  	is the resting length of the first spring
        //        resting_length[1]  	would be the resting length of the second spring associated with that Lagrangian point.
        //
        // In this example, the resting length is increased by 0.01*dt each time step.

        const int lag_idx = node_idx->getLagrangianIndex();
	//there are two ways to get resting lenghts depending on the IBAMR version. Both are copied here.
	//std::vector<double>& resting_length = spring_spec->getRestingLengths();
	//double resting_length = spring_spec->getParameters()[0][1];
    
	//Note that you can also getStiffnesses
	double spring_stiffness = spring_spec->getParameters()[0][0];
	
	if (tt<=t1) {
	
		if ((lag_idx>=x0) && (lag_idx<=x1)) {
		
			spring_stiffness[lag_idx] = spring_spec->spring_hard;
		
		} else {
			
			spring_stiffness[lag_idx] = spring_spec->spring_soft;
			
		}
	} else if ( (tt>t1) && (tt<=(t1+t2)) ) {
	
		if ((lag_idx>=x0n) && (lag_idx<=x1n)) {
		
			spring_stiffness[lag_idx] = spring_hard;
		
		} else {
			
			spring_stiffness[lag_idx] = spring_soft;
			
		}
	} else if ( (tt > (t1+t2)) && ( tt <= (t1+t2+t3)) ) {
	
		if ((lag_idx>=x0n) && (lag_idx<=x1n)) {
		
			spring_stiffness[lag_idx] = spring_hard;
		
		} else {
			
			spring_stiffness[lag_idx] = spring_soft;
			
		}
	}
		//Note that you can also getStiffnesses
	  	// if (plate2d_left_lag_idxs.first <= lag_idx && lag_idx < plate2d_left_lag_idxs.second)
		// 	{
		// 	   //resting_length[0]+=0.01*dt;
		// 	   //I've commented this out, but you could make the left plate grow and the right shrink
		// 	 }
		// if (plate2d_rght_lag_idxs.first <= lag_idx && lag_idx < plate2d_rght_lag_idxs.second)
	    {
			//resting_length[0]-=0.01*dt;
			//I've commented this out, but you could make the left plate grow and the right shrink
	    }
		

    }
    return;
}// update_springs
