#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

__global__ void pairwise_accels(vector3* d_hPos, double* mass, vector3* d_hAccels, int n){
    //set which entities to iterate on
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    //ensure blocks stay in scope
    if (i >= n || j >= n){
        return;
    }

    //no gravity on self
    if (i == j){
        FILL_VECTOR(d_hAccels[i * n + j], 0,0,0)
    }
    else {
        vector3 d_dist;

        //create distance vector
        double dx = d_hPos[i][0] - d_hPos[j][0];
        double dy = d_hPos[i][1] - d_hPos[j][1];
        double dz = d_hPos[i][2] - d_hPos[j][2];
        FILL_VECTOR(d_dist, dx,dy,dz);

        //magnitude
        double magnitude_sq = d_dist[0] * d_dist[0] + d_dist[1] * d_dist[1] + d_dist[2] * d_dist[2];
        double magnitude = sqrt(magnitude_sq);

        //acceleration magnitude
        double accel_mag = -GRAV_CONSTANT * mass[j] / magnitude_sq;

        //create acceleration vector
        vector3 accel;
        FILL_VECTOR(accel, accel_mag * d_dist[0] / magnitude,
                           accel_mag * d_dist[1] / magnitude,
                           accel_mag * d_dist[2] / magnitude);
        
        //store acceleration vector
        d_hAccels[i * n + j][0] = accel[0];
        d_hAccels[i * n + j][1] = accel[1];
        d_hAccels[i * n + j][2] = accel[2];
    }

}

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
void compute(){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	int i,j,k;
	vector3* values=(vector3*)malloc(sizeof(vector3)*NUMENTITIES*NUMENTITIES);
	vector3** accels=(vector3**)malloc(sizeof(vector3*)*NUMENTITIES);
	for (i=0;i<NUMENTITIES;i++)
		accels[i]=&values[i*NUMENTITIES];
	//first compute the pairwise accelerations.  Effect is on the first argument.
	for (i=0;i<NUMENTITIES;i++){
		for (j=0;j<NUMENTITIES;j++){
			if (i==j) {
				FILL_VECTOR(accels[i][j],0,0,0);
			}
			else{
				vector3 distance;
				//calculate distances between bodies
				for (k=0;k<3;k++) distance[k]=hPos[i][k]-hPos[j][k];

				//calculate magnitude between bodies
				double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
				double magnitude=sqrt(magnitude_sq);

				//calculate acceleration magnitude
				double accelmag=-1*GRAV_CONSTANT*mass[j]/magnitude_sq;
				
				//store acceleration magnitude of j on i in accels[i][j]
				FILL_VECTOR(accels[i][j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
			}
		}
	}
	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	for (i=0;i<NUMENTITIES;i++){
		vector3 accel_sum={0,0,0};
		for (j=0;j<NUMENTITIES;j++){
			for (k=0;k<3;k++)
				//running sum of all accelerations based on every other entity
				accel_sum[k]+=accels[i][j][k];
		}
		//compute the new position based on the velocity and time interval
		for (k=0;k<3;k++){
			//compute the new velocity based on the acceleration and time interval
			hVel[i][k]+=accel_sum[k]*INTERVAL;
			//compute the new position based on the velocity and time interval
			hPos[i][k]+=hVel[i][k]*INTERVAL;
		}
	}
	free(accels);
	free(values);
}
