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
//       8 neighbors using manual copy for non
//       contiguous side and blocking communications
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
//     - Manual copy for non continguous cells
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/

int calculateDimY(int size){
	int y = 2;
	int i;
	for( i = y ; i < size / 2; i++ ){
		if( size % i == 0){
			y = i;
		}
	}

	return y;
}

void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
	//
	// TODO: calculate the splitting parameters for the current task.
	//
	//int rank; 
	//printf("\nStart of Init exercise 4\n");
	int comm_size;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	int dim_y = calculateDimY(comm_size) ;
	//todo a changer
	if(comm_size % dim_y != 0 || (total_width % (comm_size/dim_y)) != 0 || total_height % dim_y != 0){
		exit(1);
	}
	
	//Creation de la MPI_cart
	
	
	int dims[2] = {comm_size/dim_y, dim_y};
	MPI_Dims_create(comm_size, 2, dims);

	int periods[2] = {0, 0};

	int reorder = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &(comm->communicator));
    
	// TODO: calculate the number of tasks along X axis and Y axis.
	comm->nb_x = comm_size / dim_y;
	comm->nb_y = dim_y;

	//Get the current ranl
    int my_rank = 0;
    MPI_Comm_rank((comm->communicator), &my_rank);
	
	//Get the coord
	int my_coords[2] = {0, 0};
    MPI_Cart_coords((comm->communicator), my_rank, 2, my_coords);
	//printf("\nThere before exit\n");
	// TODO: calculate the current task position in the splitting
	comm->rank_x = my_coords[0];
	comm->rank_y = my_coords[1];

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = ( total_width ) / (comm_size / dim_y ) + 2 ;
	comm->height = ((total_height)/dim_y) + 2;

	// TODO : calculate the absolute position  (in cell number) in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = (comm->width - 2) * my_coords[0] ;
	comm->y = (comm->height - 2) * my_coords[1] ;

	//OPTIONAL : if you want to avoid allocating temporary copy buffer
	//           for every step :
	//comm->buffer_recv_down, comm->buffer_recv_up, comm->buffer_send_down, comm->buffer_send_up
	//allocating bufffer
	comm->buffer_send_up = malloc(sizeof(double) * 9 * (comm->width));
	comm->buffer_send_down = malloc(sizeof(double) * 9 * (comm->width));
	comm->buffer_recv_up = malloc(sizeof(double) * 9 * (comm->width));
	comm->buffer_recv_down = malloc(sizeof(double) * 9 * (comm->width));
	//if debug print comm
	//lbm_comm_print(comm);
	#ifndef NDEBUG
	lbm_comm_print( comm );
	#endif
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	//free allocated ressources
	/*free(comm->buffer_send_up);
	free(comm->buffer_send_down);
	free(comm->buffer_recv_up);
	free(comm->buffer_recv_down);*/
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - manual copy in temp buffer for non contiguous side 
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
		int i,j;
		for(i = 0; i < comm->width; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_down[i*9+j] = lbm_mesh_get_cell(mesh, i, comm->height - 2)[j];
            }
        }
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		//Buffer full
        MPI_Send(comm->buffer_send_down, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator);

		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        MPI_Recv(comm->buffer_recv_down, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i, comm->height - 1)[j] = comm->buffer_recv_down[i*9+j];
            }
        }
    }
    else if(comm->rank_y == comm->nb_y - 1){
		//Haut
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		MPI_Recv(comm->buffer_recv_up, (9 * (comm->width)),MPI_DOUBLE,rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		int i,j;
		for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i, 0)[j] = comm->buffer_recv_up[i*9+j];
            }
        }

        for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_up[i*9+j] = lbm_mesh_get_cell(mesh, i , 1)[j];
            }
        }
		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y -1 ;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        MPI_Send(comm->buffer_send_up, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator);
	}
    else{ 
		int i,j;
		//Haut
		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		MPI_Recv(comm->buffer_recv_up, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i, 0)[j] = comm->buffer_recv_up[i*9+j];
            }
        }

        for(i = 0; i < comm->width; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_up[i*9+j] = lbm_mesh_get_cell(mesh, i , 1)[j];
            }
        }
		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        MPI_Send(comm->buffer_send_up, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator);
		
		//Bas
		for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                comm->buffer_send_down[i*9+j] = lbm_mesh_get_cell(mesh, i , comm->height - 2)[j];
            }
        }
		//Buffer full
		coords[0] = comm->rank_x ;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        MPI_Send(comm->buffer_send_down, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator);

		coords[0] = comm->rank_x;
		coords[1] = comm->rank_y + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
        MPI_Recv(comm->buffer_recv_down, (9 * (comm->width)),MPI_DOUBLE, rank, 0, comm->communicator,MPI_STATUS_IGNORE);
		for(i = 0; i < comm->width ; i++){
            for(j = 0 ; j < 9 ; j++){
                lbm_mesh_get_cell(mesh, i , comm->height - 1)[j] = comm->buffer_recv_down[i*9+j];
            }
        }
		
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
