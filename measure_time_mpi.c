#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdbool.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>


#define START_TIMEMEASUREMENT(name)   struct timeval __FILE__##__func__##name##actualtime; \
gettimeofday(&__FILE__##__func__##name##actualtime, NULL); \
double __FILE__##__func__##name##s_time = (double)__FILE__##__func__##name##actualtime.tv_sec+((double)__FILE__##__func__##name##actualtime.tv_usec/1000000.0)


#define END_TIMEMEASUREMENT(name, res) gettimeofday(&__FILE__##__func__##name##actualtime, NULL); \
res = (double)__FILE__##__func__##name##actualtime.tv_sec+((double)__FILE__##__func__##name##actualtime.tv_usec/1000000.0) -__FILE__##__func__##name##s_time


int main (int c, char **v) {

  int rank, size, i;
  MPI_Status status;

  MPI_Init(&c, &v);

  // get rank and number of processes and print it out
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("size: %d \n",size );
  printf("rank: %d \n",rank );
  // run task funktion and measure duration and calculate average and max times
  //start measure time
  double* elapsed_time_send = malloc(sizeof(double)*4);
  double* elapsed_time_recv = malloc(sizeof(double)*4);
  if(rank == 0){
  for(int i = 0; i<4; i++){
    elapsed_time_send[i] = 0;
    elapsed_time_recv[i] = 0;
    }
  }
  START_TIMEMEASUREMENT(measure_game_time);
  int sleeptime = 0;
  while (sleeptime < 1000000) {
      sleeptime++;
  }
  usleep(sleeptime);
  //end measure time
  END_TIMEMEASUREMENT(measure_game_time, elapsed_time_recv[rank]);
  elapsed_time_send[rank] = elapsed_time_recv[rank];
  printf("rank: %d time elapsed: %lf sec\n", rank, elapsed_time_recv[rank]);
  int comm = 99;
    if (rank != 0 && rank < size-1) {
      for(int i = 0; i <4; i ++){
        printf("%d sent %lf !!\n",rank, elapsed_time_send[i]);
        printf("%d recv %lf !!\n",rank, elapsed_time_recv[i]);
      }
      // Receive from left worker
      MPI_Recv(elapsed_time_recv, 4, MPI_INT,(rank-1), comm, MPI_COMM_WORLD, &status);
      elapsed_time_send[rank-1] = elapsed_time_recv[rank-1];
      // Send to right
      MPI_Send(elapsed_time_send, 4, MPI_INT,(rank+1), comm, MPI_COMM_WORLD);
      for(int i = 0; i <4; i ++){
        printf("%d sent %lf !!\n",rank, elapsed_time_send[i]);
        printf("%d recv %lf !!\n",rank, elapsed_time_recv[i]);
      }
    }
    else if (rank == 0) {
      // Send to right
      MPI_Send(elapsed_time_send, 4, MPI_INT,(rank+1), comm, MPI_COMM_WORLD);
      for(int i = 0; i <4; i ++){
        printf("%d sent %lf !! \n",rank, elapsed_time_send[i]);
        printf("%d recv %lf !! \n",rank, elapsed_time_recv[i]);
      }
    }
    else{
      MPI_Recv(elapsed_time_recv, 4, MPI_INT,(rank-1), comm, MPI_COMM_WORLD, &status);
      elapsed_time_send[rank-1] = elapsed_time_recv[rank-1];
      for(int i = 0; i < 4; i++){
        printf("%lf \n", elapsed_time_recv[i]);
      }
    }
  MPI_Finalize();
  return 0;
}
