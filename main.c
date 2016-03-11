#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001
#define TAG 0

int main (int argc, char* argv[]){
	struct node *nodehead;
	int nodecount;
	int *num_in_links, *num_out_links;
	double *r, *r_pre, *r_recv;
	int i, j, k;
	double damp_const;
	int iterationcount = 0;
	int collected_nodecount;
	double cst_addapted_threshold;
	FILE *fp;

	int numprocs;
	int myid;
	MPI_Status stat;

	// Load the data and simple verification
	if ((fp = fopen("data_output", "r")) == NULL ) {
		printf("Error loading the data_output.\n");
		return 253;
	}
	if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;

	fclose(fp);

	// Adjust the threshold according to the problem size
	cst_addapted_threshold = THRESHOLD;
	
	// Calculate the result
	if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;
	
	r = malloc(nodecount * sizeof(double));
	r_pre = malloc(nodecount * sizeof(double));
	for ( i = 0; i < nodecount; ++i)
		r[i] = 1.0 / nodecount;
	damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
	// CORE CALCULATION
	do{
		++iterationcount;
		vec_cp(r, r_pre, nodecount);
		for ( i = 0; i < nodecount; ++i){
			r[i] = 0;

			MPI_Init(NULL, NULL);
			MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
			MPI_Comm_rank(MPI_COMM_WORLD, &myid);

			if(myid == 0) {
				printf("%d: We have %d processors\n", myid, numprocs);
				MPI_Bcast(r, nodecount, MPI_DOUBLE, k, TAG, MPI_COMM_WORLD);

				for(k=1;k<numprocs;k++) {
					MPI_Recv(r_recv, nodecount, MPI_DOUBLE, k, TAG, MPI_COMM_WORLD, &stat);
				}
			}
			else {
				MPI_Recv(r, nodecount, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);
				
				for ( j = 0; j < nodehead[i].num_in_links; ++j)
					r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
				r[i] *= DAMPING_FACTOR;
				r[i] += damp_const;
				
				MPI_Send(r, nodecount, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
			}

			MPI_Finalize();
		}

	} while(rel_error(r, r_pre, nodecount) >= EPSILON);

	// post processing
	node_destroy(nodehead, nodecount);
	//free(num_in_links); free(num_out_links);
	
	Lab4_saveoutput(r, nodecount, 0.0);

	free(r); free(r_pre);
}
