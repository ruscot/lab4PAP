/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication scheme with
//      8 neighbors using MPI types for non contiguous
//      side.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
// NEW:
//     - >>> MPI type for non contiguous cells <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex5(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5 execpt for type creation
	lbm_comm_init_ex4(comm, total_width, total_height);

	//TODO: create MPI type for non contiguous side in comm->type.
	MPI_Type_vector( comm->width, 9, comm->height*9, MPI_DOUBLE, &(comm->type));
	MPI_Type_commit(&(comm->type));
}

/****************************************************/
void lbm_comm_release_ex5(lbm_comm_t * comm)
{
	//we use the same implementation than ex5 except for type destroy
	lbm_comm_release_ex4(comm);

	//TODO: release MPI type created in init.
	MPI_Type_free(&(comm->type));
}

/****************************************************/
void lbm_comm_ghost_exchange_ex5(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
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
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
        MPI_Send(cell,(9 * (comm->height)),MPI_DOUBLE, rank, 0,comm->communicator);

		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);

        cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
        MPI_Recv(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator,MPI_STATUS_IGNORE);
    }
    else if(comm->rank_x == comm->nb_x - 1){
		//gauche
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
        MPI_Recv(cell,(9 * (comm->height)),MPI_DOUBLE,rank, 0,comm->communicator,MPI_STATUS_IGNORE);
        
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, 1, 0);
        MPI_Send(cell,(9 * (comm->height)),MPI_DOUBLE,rank, 0,comm->communicator);
    }
    else{ 
		//gauche
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
        MPI_Recv(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator,MPI_STATUS_IGNORE);
        
		coords[0] = comm->rank_x - 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, 1, 0);
        MPI_Send(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator);
		
		// droite
		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
        MPI_Send(cell,(9 * (comm->height)),MPI_DOUBLE,rank,0,comm->communicator);

		coords[0] = comm->rank_x + 1;
		coords[1] = comm->rank_y;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
        MPI_Recv(cell,(9 * (comm->height)),MPI_DOUBLE, rank, 0,comm->communicator,MPI_STATUS_IGNORE);

	}


	//printf("\nAvant Y\n");
	
	//For the y axis
	if(comm->rank_y == 0){
		//Bas
		//int i,j;
		/*for(i = 0; i < comm->width; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_down[i*9+j] = lbm_mesh_get_cell(mesh, i, comm->height - 2)[j];
            }
        }*/
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		//Buffer full
		
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 2);
        MPI_Send(cell, 1,comm->type, rank, 0, comm->communicator);

		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
        MPI_Recv(cell, 1,comm->type, rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		/*for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i, comm->height - 1)[j] = comm->buffer_recv_down[i*9+j];
            }
        }*/
    }
    else if(comm->rank_y == comm->nb_y - 1){
		//Haut
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
		MPI_Recv(cell, 1,comm->type,rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		//int i,j;
		/*for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i, 0)[j] = comm->buffer_recv_up[i*9+j];
            }
        }*/

        /*for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_up[i*9+j] = lbm_mesh_get_cell(mesh, i , 1)[j];
            }
        }*/
		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y -1 ;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 1);
        MPI_Send(cell, 1,comm->type, rank, 0, comm->communicator);
	}
    else{ 
		//int i,j;
		//Haut
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 0);
		MPI_Recv(cell, 1,comm->type, rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		/*for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i, 0)[j] = comm->buffer_recv_up[i*9+j];
            }
        }*/

        /*for(i = 0; i < comm->width; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_up[i*9+j] = lbm_mesh_get_cell(mesh, i , 1)[j];
            }
        }*/
		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, 1);
        MPI_Send(cell, 1,comm->type, rank, 0, comm->communicator);
		
		//Bas
		/*for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_down[i*9+j] = lbm_mesh_get_cell(mesh, i , comm->height - 2)[j];
            }
        }*/
		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 2);
        MPI_Send(cell, 1,comm->type, rank, 0, comm->communicator);

		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
        MPI_Recv(cell, 1,comm->type, rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		/*for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i , comm->height - 1)[j] = comm->buffer_recv_down[i*9+j];
            }
        }*/
		
	}

	//printf("\nAvant diagonal\n");

	//For the diag
	//todo
	if(comm->rank_x == 0){
		if(comm->rank_y == 0){
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9,MPI_DOUBLE,rank, 0,comm->communicator);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
		} else if(comm->rank_y == comm->nb_y - 1){
			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
		}
    	else{ 
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9,MPI_DOUBLE,rank, 0,comm->communicator);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
		}
	} else if(comm->rank_x == comm->nb_x - 1){
		if(comm->rank_y == 0){
			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
		} else if(comm->rank_y == comm->nb_y - 1){
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
		} else{ 
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);

			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
		}
	} else {
		if(comm->rank_y == 0){
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9,MPI_DOUBLE, rank, 0,comm->communicator);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
		} else if(comm->rank_y == comm->nb_y - 1){
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);

			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
		} else{ 
			//Bas droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width-2, comm->height-2);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9,MPI_DOUBLE, rank, 0,comm->communicator);
			//Bas droite receive
			cell = lbm_mesh_get_cell(mesh, comm->width-1, comm->height-1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
			
			//haut gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Haut gauche send 
			cell = lbm_mesh_get_cell(mesh, 1, 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
			
			//Haut droite send 
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 1);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
			//Haut droite recv
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			coords[0] = comm->rank_x + 1;
			coords[1] = comm->rank_y - 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);

			//Bas gauche recv
			cell = lbm_mesh_get_cell(mesh, 0, comm->height - 1);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Recv(cell, 9, MPI_DOUBLE, rank, 0, comm->communicator, MPI_STATUS_IGNORE);
			//Bas gauche send
			cell = lbm_mesh_get_cell(mesh, 1, comm->height - 2);
			coords[0] = comm->rank_x - 1;
			coords[1] = comm->rank_y + 1;
			MPI_Cart_rank(comm->communicator, coords, &rank);
			MPI_Send(cell, 9, MPI_DOUBLE, rank, 0,comm->communicator);
		}
	}

	//printf("\nEnd function\n");
	//free(coords);
}
