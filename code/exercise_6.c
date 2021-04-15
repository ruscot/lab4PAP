/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication with non-blocking
//       messages.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - MPI type for non contiguous cells
// NEW:
//     - Non-blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex6(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5
	lbm_comm_init_ex5(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_release_ex6(lbm_comm_t * comm)
{
	//we use the same implementation than ext 5
	lbm_comm_release_ex5(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex6(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - non-blocking MPI functions
	//         - use MPI type for non contiguous side 
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	// TIP: create a function to get the target rank from x,y task coordinate.
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/write communication
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.
	// TIP: The previous trick require to make two batch of non-blocking communications.

	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications
	double *cell;
	int rank;
	int coords[2];
	//Debut du test
	//printf("\nCoords %d, %d, comm x : %d, comm y : %d\n", comm->rank_x,comm->rank_y, comm->x, comm->y);
	//X axis
	if(comm->rank_x == 0){
        // droite
		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Request request;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
        MPI_Isend(cell,(9 * (comm->height)),MPI_DOUBLE, rank, 0,comm->communicator, &request);

		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		MPI_Request request2;
        cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
        MPI_Irecv(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator,&request2);
    }
    else if(comm->rank_x == comm->nb_x - 1){
		//gauche
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
		MPI_Request request3;
        MPI_Irecv(cell,(9 * (comm->height)),MPI_DOUBLE,rank, 0,comm->communicator,&request3);
        
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, 1, 0);
		MPI_Request request4;
        MPI_Isend(cell,(9 * (comm->height)),MPI_DOUBLE,rank, 0,comm->communicator, &request4);
    }
    else{ 
		//gauche
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
		MPI_Request request5;
        MPI_Irecv(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator,&request5);
        
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, 1, 0);
		MPI_Request request6;
        MPI_Isend(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator, &request6);
		
		// droite
		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
		MPI_Request request7;
        MPI_Isend(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator, &request7);

		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
		MPI_Request request8;
        MPI_Irecv(cell,(9 * (comm->height)),MPI_DOUBLE, rank, 0,comm->communicator,&request8);

	}
	
	//For the y axis
	if(comm->rank_y == 0){
		//Bas
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		MPI_Request yrequest;
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 2);
        MPI_Isend(cell, 1,comm->type, rank, 0, comm->communicator, &yrequest);

		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
		MPI_Request yrequest1;
        MPI_Irecv(cell, 1,comm->type, rank, 0, comm->communicator,&yrequest1);
    }
    else if(comm->rank_y == comm->nb_y - 1){
		//Haut
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
		MPI_Request yrequest2;
		MPI_Irecv(cell, 1,comm->type,rank, 0, comm->communicator,&yrequest2);

		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y -1 ;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 1);
		MPI_Request yrequest3;
        MPI_Isend(cell, 1,comm->type, rank, 0, comm->communicator, &yrequest3);
	}
    else{ 
		//Haut
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
		MPI_Request yrequest4;
		MPI_Irecv(cell, 1,comm->type, rank, 0, comm->communicator,&yrequest4);
		
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 1);
		MPI_Request yrequest5;
        MPI_Isend(cell, 1,comm->type, rank, 0, comm->communicator, &yrequest5);
		
		//Bas
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 2);
		MPI_Request yrequest6;
        MPI_Isend(cell, 1,comm->type, rank, 0, comm->communicator, &yrequest6);

		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
		MPI_Request yrequest7;
        MPI_Irecv(cell, 1,comm->type, rank, 0, comm->communicator,&yrequest7);	
	}

	//MPI_Barrier(comm->communicator);

	//For the diag
	//todo
	if(comm->rank_x == 0){
		if(comm->rank_y == 0){
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Request drequest;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Isend(cell, 9,MPI_DOUBLE,rank, 0,comm->communicator, &drequest);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Request drequest1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &drequest1);
		} else if(comm->rank_y == comm->nb_y - 1){
			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest2;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &drequest2);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest3;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &drequest3);
		}
    	else{ 
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest4;
			MPI_Isend(cell, 9,MPI_DOUBLE,rank, 0,comm->communicator, &drequest4);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest5;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &drequest5);

			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest6;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &drequest6);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest7;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &drequest7);
		}
	} else if(comm->rank_x == comm->nb_x - 1){
		if(comm->rank_y == 0){
			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest8;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &drequest8);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request drequest9;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &drequest9);
		} else if(comm->rank_y == comm->nb_y - 1){
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d2request;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d2request);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d2request1;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d2request1);
		} else{ 
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d2request2;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d2request2);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d2request3;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d2request3);

			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d2request4;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d2request4);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d2request5;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d2request5);
		}
	} else {
		if(comm->rank_y == 0){
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request;
			MPI_Isend(cell, 9,MPI_DOUBLE, rank, 0,comm->communicator, &d3request);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request1;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d3request1);

			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request2;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d3request2);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request3;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d3request3);
		} else if(comm->rank_y == comm->nb_y - 1){
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request4;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d3request4);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request5;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d3request5);

			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request6;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d3request6);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request7;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d3request7);
		} else{ 
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request8;
			MPI_Isend(cell, 9,MPI_DOUBLE, rank, 0,comm->communicator, &d3request8);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d3request9;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d3request9);
			
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d4request;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d4request);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d4request1;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d4request1);
			
			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d4request2;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d4request2);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d4request3;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d4request3);

			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d4request4;
			MPI_Irecv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, &d4request4);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Request d4request5;
			MPI_Isend(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator, &d4request5);
		}
	}
	//MPI_Barrier(comm->communicator);
	//MPI_Status status;
	//MPI_Wait(&request, &status);
}
