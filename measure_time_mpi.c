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

  double elapsed_time[4];

  // get rank and number of processes and print it out
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("size: %d \n",size );
  printf("rank: %d \n",rank );
  // run task funktion and measure duration and calculate average and max times
  //start measure time
  START_TIMEMEASUREMENT(measure_game_time);
  int sleeptime = 0;
  while (sleeptime < 1000000) {
      sleeptime++;
  }
  usleep(sleeptime);
  //end measure time
  END_TIMEMEASUREMENT(measure_game_time, elapsed_time[rank]);
  printf("rank: %d time elapsed: %lf sec\n", rank, elapsed_time[rank]);
  if(rank != 0){
    MPI_Send(elapsed_time, 10, MPI_INT, 1, 123, MPI_COMM_WORLD);
  }
  else{
    printf("%s \n",elapsed_time[0]);
  }
  MPI_Finalize();
  return 0;
}
