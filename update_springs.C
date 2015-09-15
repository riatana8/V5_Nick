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
    static const int N1 = 512; // number of cartesian grid meshwidths at the finest level of the AMR grid

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
	double resting_length = spring_spec->getParameters()[0][1];
    
	//Note that you can also getStiffnesses
	double spring_stiffness = spring_spec->getParameters()[0][0];
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
