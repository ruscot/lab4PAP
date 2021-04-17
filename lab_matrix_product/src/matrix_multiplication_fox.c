#include <stdio.h>
#include <mpi.h>

#include <stdlib.h>

#include "utils.h"

/* a matrix multiplication without locality (column-first)*/

void parallelMatrixMultiplication(int dimension, double *A, double *B, double *C, int rank, int w_size)
{
    int sizePerRank = dimension*dimension / (w_size - 1);
    int i = 0;
    int j = 0;
    int k = 0;
    //Nombre dans C
    for (i = sizePerRank * (rank - 1); i < sizePerRank; i++) {
        for (j = 0; j < dimension; j++) {
            C[i + (rank - 1) * dimension] += A[(rank - 1) + j * dimension] * B[j + (rank - 1) * dimension];
        }
    }
}

void parallelFoxMatrixMultiplication(int dimension, double *A, double *B, double *C, int rank, int w_size)
{
    //Avec moins de coeur que de dimension
    int nbCasePerRank = (dimension*dimension) / (w_size - 1);
    int i = 0;
    int j = 0;
    int k = 0;
    int colonne,ligne;
    double nbCourantA = 0;
    double saveFirst;
    //printf("\nDebut du calcul local\n");
    //On suppose que l'on a un coeur pour chaque case de notre matrice
    for(k = 0; k < dimension; k++){
        //Pour chacune des cases que l'on doit calculer
        for(i = nbCasePerRank * (rank - 1); i < nbCasePerRank * rank; i++){
            //La case que l'on souhaite changer est la case i
            //Il faut récupérer le la diagonale courante 
            ligne = (nbCasePerRank * (rank - 1)) % dimension;
            nbCourantA = A[ligne + dimension * k];
            //Local computation
            C[i] = C[i] + nbCourantA * B[i];
        }
        
        //Vertical shift of B
        for(i = 0; i < dimension ; i++){
            saveFirst = B[i * dimension];
            for(j = 0; j < dimension - 2; j++){
                B[j + i * dimension] = B[j + i * dimension + 1];
            }
            B[(i + 1) * dimension - 1] = saveFirst;
        }
    } 
    //printf("\nFin du calcul local\n");
    
}

int main(int argc, char *argv[])
{
    unsigned int exp ;
    double *A, *B ,*C;
    double *A_check, *B_check ,*C_check;

    unsigned int mat_size=0;

    int my_rank;
    int w_size;

    double start=0, av=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &w_size);

    if(argc != 2){
        printf("usage: %s matrix_size\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    else{
        mat_size = atoi(argv[1]);
        if(w_size <= 1 || mat_size % (w_size-1) != 0 ) {
            printf("\nmat size can't be divide by number of rank - 1 : %d matrix_size, %d number of rank - 1\n",mat_size, w_size - 1);
            MPI_Finalize();
            return 0;
        }
    }
    
    if(my_rank == 0){

        printf("test with a matrix of size %u x %u\n",mat_size, mat_size);

        
        A = allocMatrix(mat_size);
        B = allocMatrix(mat_size);
        C = allocMatrix(mat_size);

    }

#ifdef PERF_EVAL
    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        if(my_rank == 0){
            
            initMatrix(mat_size, A);
            initMatrix(mat_size, B);
            initMatrixZero(mat_size, C);
            
            start = MPI_Wtime();
            sequentialMatrixMultiplication_REF(mat_size, A, B , C);
            
            experiments [exp] = MPI_Wtime() - start;
        } 
    }

    if(my_rank == 0){
        av = average_time() ;  
        printMatrix(mat_size, A);
        printMatrix(mat_size, B);
        printMatrix(mat_size, C);
        printf ("\n REF sequential time \t\t\t %.3lf seconds\n\n", av) ;
    }
    
    
    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        if(my_rank == 0){
            initMatrix(mat_size, A);
            initMatrix(mat_size, B);
            initMatrixZero(mat_size, C);
            
            start = MPI_Wtime();
            int numberOfRank;
            //sequentialMatrixMultiplication(mat_size, A, B , C);
            //On envoie en Broadcast les données de A et B
            for(numberOfRank = 1; numberOfRank < w_size; numberOfRank++){
                //On récupère le calcul de C dans des différents coeurs
                //On l'ajoute à la matrice C
                MPI_Send(A, mat_size*mat_size, MPI_DOUBLE, numberOfRank, 0, MPI_COMM_WORLD);
                MPI_Send(B, mat_size*mat_size, MPI_DOUBLE, numberOfRank, 0, MPI_COMM_WORLD);

            } 
            
            double tmp[mat_size*mat_size];
            for(numberOfRank = 1; numberOfRank < w_size; numberOfRank++){
                //On récupère le calcul de C dans des différents coeurs
                //printf("\nAvant de tous mettre dans C\n");
                MPI_Recv(tmp, mat_size*mat_size, MPI_DOUBLE, numberOfRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //On l'ajoute à la matrice C
                int i = 0;
                int nbCasePerRank = (mat_size*mat_size) / (w_size - 1);
                for(i = nbCasePerRank * (numberOfRank - 1); i < nbCasePerRank * numberOfRank; i++){
                    C[i] = tmp[i];
                }
                //printf("\nApres avoir tous mit dans C\n");
            } 
            //free(tmp);
            experiments [exp] = MPI_Wtime() - start;
        } else {
            //On réceptionne les matrices A et B de 0
            double tmpA[mat_size*mat_size];
            double tmpB[mat_size*mat_size];
            MPI_Recv(tmpA, mat_size*mat_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(tmpB, mat_size*mat_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            //On initialise la matrice locale C
            double tmp1[mat_size*mat_size];
            int i;
            for(i = 0; i < mat_size*mat_size; i++){
                tmp1[i] = 0;
            }
            
            //On commence l'algorithme de fox
            parallelFoxMatrixMultiplication(mat_size, tmpA, tmpB, tmp1, my_rank, w_size);
            //printf("\nAvant d'envoyer à 0 \n");
            //Il faut envoyer au rank 0 la matrice C que l'on a calculé
            MPI_Send(tmp1, mat_size*mat_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            //printf("\nApres avoir envoyé à 0 \n");
        
        }
        
    }

    if(my_rank == 0){
        av = average_time() ;  
        printMatrix(mat_size, A);
        printMatrix(mat_size, B);
        printMatrix(mat_size, C);
        
        printf ("\n my mat_mult \t\t\t %.3lf seconds\n\n", av) ;
    }
    

#endif /*PERF_EVAL*/

#ifdef CHECK_CORRECTNESS
    /* running my sequential implementation of the matrix
       multiplication */
    if(my_rank == 0){
        initMatrix(mat_size, A);
        initMatrix(mat_size, B);
        initMatrixZero(mat_size, C);

        A_check = createMatrixCopy(mat_size, A);
        B_check = createMatrixCopy(mat_size, B);
        C_check = allocMatrix(mat_size);

        initMatrixZero(mat_size, C_check);
    }

    /* check for correctness */
    if(my_rank == 0){
        sequentialMatrixMultiplication(mat_size, A, B , C);
        
        sequentialMatrixMultiplication_REF(mat_size, A_check, B_check , C_check);


        if(checkMatricesEquality(mat_size, C, C_check)){
            printf("\t CORRECT matrix multiplication result \n");
        }
        else{
            printf("\t FAILED matrix multiplication !!! \n");
        }

        /* printMatrix(mat_size, C); */
        /* printMatrix(mat_size, C_check); */
        
        free(A_check);
        free(B_check);
        free(C_check);
    }

#endif /* CHECK_CORRECTNESS */

    if(my_rank == 0){
        free(A);
        free(B);
        free(C);
    }

    MPI_Finalize();

    return 0;
}
