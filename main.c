#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Lab4_IO.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001

int main (int argc, char* argv[]){
	struct node *nodehead;
	int nodecount;
	int *num_in_links, *num_out_links;
	double *r, *r_pre, *r_local;
	double start, end;
	int i, j;
	double damp_const, relative_error;
	int collected_nodecount;
	double cst_addapted_threshold;
	FILE *fp;
	int my_rank;
	int numprocs;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	// Load the data and simple verification
	if ((fp = fopen("data_output", "r")) == NULL ){
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
	if (my_rank == 0)
		GET_TIME(start);
	do{
		vec_cp(r, r_pre, nodecount);
		//Operate on subsection of r based on rank
		for ( i = my_rank*nodecount/numprocs; i < (my_rank+1)*nodecount/numprocs; ++i){
			r[i] = 0;
			for(j = 0; j < nodehead[i].num_in_links; j++)
				r[i] += r_pre[nodehead[i].inlinks[j]]/num_out_links[nodehead[i].inlinks[j]];
			r[i] *= DAMPING_FACTOR;
			r[i] += damp_const;
		}
		MPI_Allgather(MPI_IN_PLACE, nodecount/numprocs, MPI_DOUBLE, r, nodecount/numprocs, MPI_DOUBLE, MPI_COMM_WORLD);
	}while(rel_error(r, r_pre, nodecount) >= EPSILON);


	
	if (my_rank == 0){
		GET_TIME(end);
		Lab4_saveoutput(r, nodecount, end - start);
	}
	node_destroy(nodehead, nodecount);
	free(r); 
	free(r_pre); 
	MPI_Finalize();
	return 0;
}
