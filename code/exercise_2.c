/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement odd/even 1D blocking communication scheme 
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
// NEW:
//     - >>> Odd/even communication ordering <<<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex2(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex2(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with blocking MPI functions using
	//       odd/even communications.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
	double * cell;
	if(comm->rank_x == 0){
		//gauche
		cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
        MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x + 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
		// droite
        cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
        MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x + 1),0,MPI_COMM_WORLD);
    }
    else if(comm->rank_x == comm->nb_x - 1){
		if(comm->rank_x % 2 == 0){
			//gauche
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			cell = lbm_mesh_get_cell(mesh, 1, 0);
			MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD);
		} else {
			//gauche
			cell = lbm_mesh_get_cell(mesh, 1, 0);
			MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD);

			cell = lbm_mesh_get_cell(mesh, 0, 0);
			MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
    }
    else{
		if(comm->rank_x % 2 == 0){
			//gauche	
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			//droite
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x + 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			//gauche
			cell = lbm_mesh_get_cell(mesh, 1, 0);
			MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD);
			
			// droite
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
			MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x + 1),0,MPI_COMM_WORLD);

			
		} else {
			//gauche
			cell = lbm_mesh_get_cell(mesh, 1, 0);
			MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD);
			
			// droite
			cell = lbm_mesh_get_cell(mesh, comm->width - 2, 0);
			MPI_Send(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x + 1),0,MPI_COMM_WORLD);

			//gauche	
			cell = lbm_mesh_get_cell(mesh, 0, 0);
			MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x - 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			//droite 
			cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
			MPI_Recv(cell,(9 * comm->height),MPI_DOUBLE,(comm->rank_x + 1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
		}
    }
}
