#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdbool.h>
#include <mpi.h>

//OPTIONAL: comment this out for console output
//#define CONSOLE_OUTPUT
#define START_TIMEMEASUREMENT(name)   struct timeval __FILE__##__func__##name##actualtime; \
gettimeofday(&__FILE__##__func__##name##actualtime, NULL); \
double __FILE__##__func__##name##s_time = (double)__FILE__##__func__##name##actualtime.tv_sec+((double)__FILE__##__func__##name##actualtime.tv_usec/1000000.0)


#define END_TIMEMEASUREMENT(name, res) gettimeofday(&__FILE__##__func__##name##actualtime, NULL); \
res = (double)__FILE__##__func__##name##actualtime.tv_sec+((double)__FILE__##__func__##name##actualtime.tv_usec/1000000.0) -__FILE__##__func__##name##s_time


#define calcIndex(width, x,y)  ((y)*(width) + (x))
#define ALIVE 1
#define DEAD 0

#define X 0                        //
#define Y 1                        //

int           rank = 0;            // The current MPI rank in the global communicator.
int           rank_cart = 0;       // The current MPI rank in the cart communicator.
int           num_tasks;           // The number of processes



MPI_Datatype  filetype;            //
MPI_Comm      cart_comm;           // Communicator for the cartesian grid
MPI_File      file;                // A shared file pointer
MPI_Datatype  memtype;             // A new created type for the inner and outer data including the ghost layer



void myexit (const char * s, ...) {
  va_list args;
  va_start(args, s);
  if(rank == 0){
    vprintf (s, args);
    printf ("\n");
  }
  va_end(args);
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

char vtk_header[2048];
void create_vtk_header (char * header, int width, int height, int timestep) {
  char buffer[1024];
  header[0] = '\0';
  strcat (header, "# vtk DataFile Version 3.0\n");
  snprintf (buffer, sizeof(buffer), "Gameoflife timestep %d \n", timestep);
  strcat (header, buffer);
  strcat (header, "BINARY\n");
  strcat (header, "DATASET STRUCTURED_POINTS\n");
  snprintf (buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat (header, buffer);
  strcat (header, "SPACING 1.0 1.0 1.0\n");
  strcat (header, "ORIGIN 0 0 0\n");
  snprintf (buffer, sizeof(buffer), "POINT_DATA %ld\n", width * height);
  strcat (header, buffer);
  strcat (header, "SCALARS data char 1\n");
  strcat (header, "LOOKUP_TABLE default\n");
}


void write_field (char* currentfield, int width, int height, int timestep) {
  if (timestep == 0){
    if(rank_cart == 0) {
      mkdir("./gol/", 0777);
    }
    create_vtk_header (vtk_header, width, height, timestep);
  }

  //printf ("writing timestep %d\n", timestep);
  char filename[1024];
  snprintf (filename, 1024, "./gol/gol-%05d.vtk", timestep);
  MPI_Offset header_offset = (MPI_Offset)strlen(vtk_header);

    /* TODO Create a new file handle for collective I/O
     *      Use the global 'file' variable.
   */
  /* TODO Set the file view for the file handle using collective I/O
   *
   */
  // rc = ...

  /* TODO Write the data using collective I/O
   *
   */


  /* TODO Close the file handle.
   *
   */
}

int countLifingsPeriodic(char* currentfield, int x, int y, int w, int h){
  int n = 0;
  for(int y1 = y - 1; y1 <= y+1; y1++){
    for (int x1 = x - 1; x1 <= x + 1; x1++){
      if(currentfield[calcIndex(w, (x1 +w ) % w, (y1 +h) %h)]){
        n++;
      }
    }
  }
  return n;
}

void evolve (char* currentfield, char* newfield, int width, int height) {
  // TODO traverse through each voxel and implement game of live logic
  // HINT: avoid boundaries
  for (int y = 0; y < height; y++){
    for (int x = 0; x < width; x++){
      int n = countLifingsPeriodic(currentfield, x , y , width, height);
      if (currentfield[calcIndex(width, x, y)]) n--;
      newfield[calcIndex(width, x, y)] = (n == 3 || (n == 2 && currentfield[calcIndex(width,x,y)]));
    }
  }
}

void filling_random (char * currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand () < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner (char * currentfield, int width, int height) {
  currentfield[calcIndex(width, width/2+0, height/2+1)] = ALIVE;
  currentfield[calcIndex(width, width/2+1, height/2+2)] = ALIVE;
  currentfield[calcIndex(width, width/2+2, height/2+0)] = ALIVE;
  currentfield[calcIndex(width, width/2+2, height/2+1)] = ALIVE;
  currentfield[calcIndex(width, width/2+2, height/2+2)] = ALIVE;
}

void apply_periodic_boundaries(char * field, int width, int height){
  //TODO: implement periodic boundary copies
  int i, j, k, l;
  for (int y = 0; y < height - 1; y++) {
      i = calcIndex(width, width - 1, y);
      j = calcIndex(width, 1, y);
      l = calcIndex(width, 0, y);
      k = calcIndex(width, width - 2, y);
      field[i] = field[j];
      field[l] = field[k];
  }
  int a, b, c, d;
  for (int x = 1; x < width - 1; x++) {
    a = calcIndex(width, x, height - 1);
    b = calcIndex(width, x, 1);
    d = calcIndex(width, x, 0);
    c = calcIndex(width, x, height - 2);
    field[a] = field[b];
    field[d] = field[c];
  }
}

void game (int width, int height, int num_timesteps, int gsizes[2]) {
  char *currentfield = calloc (width * height, sizeof(char));
  char *newfield = calloc (width * height, sizeof(char));

  // TODO 1: use your favorite filling
  //filling_random (currentfield, width, height);
  filling_runner (currentfield, width, height);

  int time = 0;
  write_field (currentfield, gsizes[X], gsizes[Y], time);

  for (time = 1; time <= num_timesteps; time++) {
    // TODO 2: implement evolve function (see above)
    evolve (currentfield, newfield, width, height);
    write_field (newfield, gsizes[X], gsizes[Y], time);
    apply_periodic_boundaries(newfield,width,height);
    // TODO 3: implement SWAP of the fields
    char *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }

  free (currentfield);
  free (newfield);
}

int main (int c, char **v) {

  MPI_Init(&c, &v);

  int width, height, num_timesteps;
  int process_numX;
  int process_numY;

  if (c == 6) {
    width = atoi (v[1]); ///< read width + 2 boundary cells (low x, high x)
    height = atoi (v[2]); ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi (v[3]); ///< read timesteps

    if (width <= 0) {
      width = 32; ///< default width
    }
    if (height <= 0) {
      height = 32; ///< default height
    }
    process_numX = atoi (v[4]); ///< read number of processes in X
    process_numY = atoi (v[5]); ///< read number of processes in Y

  }
  else {
    printf("%d\n",c );
   myexit("Too less arguments");
  }



  // TODO get the global rank of the process and save it to rank_global
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);        //-
  // TODO get the number of processes and save it to num_tasks variable
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);   //-

  /* Abort if the number of processes does not match with the given configuration.
   */
  if (num_tasks != (process_numX*process_numY)) {
    myexit("ERROR: %d MPI processes needed.\n", process_numX*process_numY);
  }


  /* TODO Create a new cartesian communicator of the worker communicator and get the information.
  */
  int dims[2] = {process_numX, process_numY};
  int periods[2] = {1,1};
  int coords[2];
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
  MPI_Comm_rank(cart_comm, &rank_cart);
  MPI_Cart_coords(cart_comm, rank_cart, 2, coords);
  int gsizes[2] = {width, height};  // global size of the domain without boundaries
  int lsizes[2] = {width/process_numX, height/process_numY};

  MPI_Status status;
  double* elapsed_time_send = malloc(sizeof(double)*4);
  double* elapsed_time_recv = malloc(sizeof(double)*4);
  if(rank == 0){
  for(int i = 0; i < size; i++){
    elapsed_time_send[i] = 0;
    elapsed_time_recv[i] = 0;
    }
  }
  START_TIMEMEASUREMENT(measure_game_time);
  /* TODO create and commit a subarray as a new filetype to describe the local
   *      worker field as a part of the global field.
   *      Use the global variable 'filetype'.
   * HINT: use MPI_Type_create_subarray and MPI_Type_commit functions
  */


  /* TODO Create a derived datatype that describes the layout of the inner local field
   *      in the memory buffer that includes the ghost layer (local field).
   *      This is another subarray datatype!
   *      Use the global variable 'memtype'.
  */

  game(lsizes[X], lsizes[Y], num_timesteps, gsizes);
  END_TIMEMEASUREMENT(measure_game_time, elapsed_time_recv[rank]);
  elapsed_time_send[rank] = elapsed_time_recv[rank];
  printf("rank: %d time elapsed: %lf sec\n", rank, elapsed_time_recv[rank]);
  int comm = 99;
    if (rank != 0 && rank < num_tasks-1) {
      // Receive from left worker
      MPI_Recv(elapsed_time_recv, num_tasks, MPI_DOUBLE,(rank-1), comm, MPI_COMM_WORLD, &status);
      for(int i = 0; i <rank; i ++){
        elapsed_time_send[i] = elapsed_time_recv[i];
      }
      // Send to right
      MPI_Send(elapsed_time_send, num_tasks, MPI_DOUBLE,(rank+1), comm, MPI_COMM_WORLD);
    }
    else if (rank == 0) {
      // Send to right
      MPI_Send(elapsed_time_send, num_tasks, MPI_DOUBLE,(rank+1), comm, MPI_COMM_WORLD);
    }
    else{
      MPI_Recv(elapsed_time_recv, num_tasks, MPI_DOUBLE,(rank-1), comm, MPI_COMM_WORLD, &status);
      for(int i = 0; i <rank; i ++){
        elapsed_time_send[i] = elapsed_time_recv[i];
      }
      double max = 0;
      double average = 0;
      for(int i = 0; i < num_tasks; i++){
        if(elapsed_time_send[i] > max){
          max = elapsed_time_send[i];
        }
          average = average + elapsed_time_send[i];
      }
      printf("Max is %lf. Average is %lf \n", max, average/num_tasks);
    }
  MPI_Finalize();
}
