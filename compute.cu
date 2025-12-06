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

__global__ void sum_accels(vector3* d_hAccels, vector3* accel_sum, int n){
    //set which entity to iterate on
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    //ensure in scope
    if (i > n){
        return;
    }

    //sum for entity
    vector3 sum = {0,0,0};
    for (int j = 0; j < n; j++){
        sum[0] += d_hAccels[i * n + j][0];
        sum[1] += d_hAccels[i * n + j][1];
        sum[2] += d_hAccels[i * n + j][2];
    }

    //store accel for entity
    accel_sum[i][0] = sum[0];
    accel_sum[i][1] = sum[1];
    accel_sum[i][2] = sum[2];
}

__global__ void update_positions(vector3* d_hPos, vector3* d_hVel, vector3* accel_sum, int n){
    //set which entity to iterate on
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    //ensure in scope
    if (i > n){
        return;
    }

    //create vectors for updated values
    vector3 new_v, new_p;
    FILL_VECTOR(new_v, d_hVel[i][0], d_hVel[i][1], d_hVel[i][2]);
    FILL_VECTOR(new_v, d_hPos[i][0], d_hPos[i][1], d_hPos[i][2]);

    //update values
    new_v[0] += accel_sum[i][0] * INTERVAL;
    new_v[1] += accel_sum[i][1] * INTERVAL;
    new_v[2] += accel_sum[i][2] * INTERVAL;

    new_p[0] += new_v[0] * INTERVAL;
    new_p[0] += new_v[1] * INTERVAL;
    new_p[0] += new_v[2] * INTERVAL;

    //store updated values
    d_hVel[i][0] = new_v[0];
    d_hVel[i][1] = new_v[1];
    d_hVel[i][2] = new_v[2];

    d_hPos[i][0] = new_p[0];
    d_hPos[i][1] = new_p[1];
    d_hPos[i][2] = new_p[2];
}

void compute() {
    //device pointers
    vector3* d_hPos;
    vector3* d_hVel;
    double* d_mass;
    vector3* d_hAccels;
    vector3* accel_sum;

    size_t vecSize = sizeof(vector3) * NUMENTITIES;
    size_t matSize = sizeof(vector3) * NUMENTITIES * NUMENTITIES;

    //allocate device memory
    cudaMalloc(&d_hPos, vecSize);
    cudaMalloc(&d_hVel, vecSize);
    cudaMalloc(&d_mass, sizeof(double) * NUMENTITIES);
    cudaMalloc(&d_hAccels, matSize);
    cudaMalloc(&accel_sum, vecSize);

    //copy to device
    cudaMemcpy(d_hPos, hPos, vecSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_hVel, hVel, vecSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mass, mass, sizeof(double) * NUMENTITIES, cudaMemcpyHostToDevice);

    //Pairwise acceleration calculation
    dim3 blockDim(16,16);
    dim3 gridDim((NUMENTITIES + blockDim.x - 1) / blockDim.x,
                 (NUMENTITIES + blockDim.y - 1) / blockDim.y);
    pairwise_accels<<<gridDim, blockDim>>>(d_hPos, d_mass, d_hAccels, NUMENTITIES);

    //Sum accelerations calculation
    int block1D = 256;
    int grid1D = (NUMENTITIES + block1D - 1) / block1D;
    sum_accels<<<grid1D, block1D>>>(d_hAccels, accel_sum, NUMENTITIES);

    //Update positions and velocities
    update_positions<<<grid1D, block1D>>>(d_hPos, d_hVel, accel_sum, NUMENTITIES);

    //copy back to host
    cudaMemcpy(hPos, d_hPos, vecSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(hVel, d_hVel, vecSize, cudaMemcpyDeviceToHost);

    //free device memory
    cudaFree(d_hPos);
    cudaFree(d_hVel);
    cudaFree(d_mass);
    cudaFree(d_hAccels);
    cudaFree(accel_sum);
}
