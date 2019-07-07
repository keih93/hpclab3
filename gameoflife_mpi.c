#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdbool.h>
#include <mpi.h>
#include <sys/time.h>

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
int           coords[2];
int           gsizes[2];           // global size of the domain without boundaries
int           lsizes[2];           // local size without boundaries
int           dims[2];

MPI_Datatype  filetype;            //
MPI_Comm      cart_comm;           // Communicator for the cartesian grid
MPI_File      file;                // A shared file pointer
MPI_Datatype  memtype;             // A new created type for the inner and outer data including the ghost layer
MPI_Status status;


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

  char filename[1024];
  snprintf (filename, 1024, "./gol/gol-%05d.vtk", timestep);
  MPI_Offset header_offset = (MPI_Offset)strlen(vtk_header);

  MPI_File_delete(filename,MPI_INFO_NULL);
  MPI_File_open(cart_comm,  filename ,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

  if(rank_cart == 0){
  MPI_File_write(file, vtk_header, strlen(vtk_header), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  MPI_File_set_view(file, header_offset, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(file, currentfield,1, memtype, MPI_STATUS_IGNORE);

  MPI_File_close(&file);
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

void filling_rank (char * currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] = rank_cart;
    }
  }
}

void filling_runner1 (char * currentfield, int width, int height) {
  int offset_x = width/3;
  int offset_y = height/2;
  currentfield[calcIndex(width, offset_x+0, offset_y+1)] = ALIVE;
  currentfield[calcIndex(width, offset_x+1, offset_y+2)] = ALIVE;
  currentfield[calcIndex(width, offset_x+2, offset_y+0)] = ALIVE;
  currentfield[calcIndex(width, offset_x+2, offset_y+1)] = ALIVE;
  currentfield[calcIndex(width, offset_x+2, offset_y+2)] = ALIVE;
}

void filling_runner (char * currentfield, int width, int height) {
  if( ((gsizes[0]/2+0) >= (coords[0]*lsizes[0])) && ((gsizes[1]/2+1) >= (coords[1]*lsizes[1])) ){
    if( ((gsizes[0]/2+0) < (coords[0]*lsizes[0] +lsizes[0])) && ((gsizes[1]/2+1) < (coords[1]*lsizes[1] + lsizes[1])) ){
        currentfield[calcIndex(width, (gsizes[0]/2+0) - (coords[0]*lsizes[0])+1, (gsizes[1]/2+1)-(coords[1]*lsizes[1])+1)] = ALIVE;
        printf("rank %d: first cell gindex %d %d lindex %d %d \n",rank_cart,gsizes[0]/2+0,gsizes[1]/2+1,coords[0]*lsizes[0],coords[1]*lsizes[1]);
    }
  }
  if( ((gsizes[0]/2+1) >= (coords[0]*lsizes[0])) && ((gsizes[1]/2+2) >= (coords[1]*lsizes[1])) ){
    if( ((gsizes[0]/2+1) < (coords[0]*lsizes[0] +lsizes[0])) && ((gsizes[1]/2+2) < (coords[1]*lsizes[1] + lsizes[1])) ){
        currentfield[calcIndex(width, (gsizes[0]/2+1) - (coords[0]*lsizes[0])+1, (gsizes[1]/2+2)-(coords[1]*lsizes[1])+1)] = ALIVE;
        printf("rank %d: second cell gindex %d %d lindex %d %d \n",rank_cart,gsizes[0]/2+1,gsizes[1]/2+2,coords[0]*lsizes[0],coords[1]*lsizes[1]);
    }
  }
  if( ((gsizes[0]/2+2) >= (coords[0]*lsizes[0])) && ((gsizes[1]/2+0) >= (coords[1]*lsizes[1])) ){
    if( ((gsizes[0]/2+2) < (coords[0]*lsizes[0] +lsizes[0])) && ((gsizes[1]/2+0) < (coords[1]*lsizes[1] + lsizes[1])) ){
        currentfield[calcIndex(width, (gsizes[0]/2+2) - (coords[0]*lsizes[0])+1, (gsizes[1]/2+0)-(coords[1]*lsizes[1])+1)] = ALIVE;
        printf("rank %d: third cell gindex %d %d lindex %d %d \n",rank_cart,gsizes[0]/2+2,gsizes[1]/2+0,coords[0]*lsizes[0],coords[1]*lsizes[1]);
    }
  }
  if( ((gsizes[0]/2+2) >= (coords[0]*lsizes[0])) && ((gsizes[1]/2+1) >= (coords[1]*lsizes[1])) ){
    if( ((gsizes[0]/2+2) < (coords[0]*lsizes[0] +lsizes[0])) && ((gsizes[1]/2+1) < (coords[1]*lsizes[1] + lsizes[1])) ){
        currentfield[calcIndex(width, (gsizes[0]/2+2) - (coords[0]*lsizes[0])+1, (gsizes[1]/2+1)-(coords[1]*lsizes[1])+1)] = ALIVE;
        printf("rank %d: fourth cell gindex %d %d lindex %d %d \n",rank_cart,gsizes[0]/2+2,gsizes[1]/2+1,coords[0]*lsizes[0],coords[1]*lsizes[1]);
    }
  }
  if( ((gsizes[0]/2+2) >= (coords[0]*lsizes[0])) && ((gsizes[1]/2+2) >= (coords[1]*lsizes[1])) ){
    if( ((gsizes[0]/2+2) < (coords[0]*lsizes[0] +lsizes[0])) && ((gsizes[1]/2+2) < (coords[1]*lsizes[1] + lsizes[1])) ){
        currentfield[calcIndex(width, (gsizes[0]/2+2) - (coords[0]*lsizes[0])+1, (gsizes[1]/2+2)-(coords[1]*lsizes[1])+1)] = ALIVE;
        printf("rank %d: fifth cell gindex %d %d lindex %d %d \n",rank_cart,gsizes[0]/2+2, gsizes[1]/2+2,coords[0]*lsizes[0],coords[1]*lsizes[1]);
    }
  }
}

void apply_periodic_boundaries(char * field, int width, int height){
  //TODO: implement periodic boundary copies
  printf("switching boundaries\n");
    char *sendcellstb[2];
    char *recvcellstb[2];
    char *sendcellslr[2];
    char *recvcellslr[2];
    for (int a = 0; a < 2; a++) {
        recvcellstb[a] = calloc(width + 1, sizeof(char));
        sendcellstb[a] = calloc(width + 1, sizeof(char));
    }
    for (int a = 0; a < 2; a++) {
        recvcellslr[a] = calloc(height + 1, sizeof(char));
        sendcellslr[a] = calloc(height + 1, sizeof(char));
    }
  int toprank, botrank, leftrank, rightrank;
  int siderank[4] = {num_tasks,num_tasks,num_tasks,num_tasks};
  int topcoords[2], botcoords[2], leftcoords[2], rightcoords[2];
  int sidecoords[4][2];
  int maxcoords[2];
  // side cells
  char sidecells[4]={DEAD,DEAD,DEAD,DEAD};//sidedownleft, sideupleft, sidedownright, sideupright
  char recvsidecells[4]={DEAD,DEAD,DEAD,DEAD};
  MPI_Cart_coords(cart_comm, num_tasks-1, 2, maxcoords);
  if(coords[1] == maxcoords[1]){
    topcoords[0]=coords[0];
    topcoords[1]=0;
    MPI_Cart_rank(cart_comm, topcoords,&toprank);
  }else{
    topcoords[0]=coords[0];
    topcoords[1]=coords[1]+1;
    MPI_Cart_rank(cart_comm, topcoords,&toprank);
  }
  if(coords[1] == 0){
    botcoords[0]=coords[0];
    botcoords[1]=maxcoords[0];
    MPI_Cart_rank(cart_comm, botcoords,&botrank);
  }else{
    botcoords[0]=coords[0];
    botcoords[1]=coords[1]-1;
    MPI_Cart_rank(cart_comm, botcoords,&botrank);
  }
  if(coords[0] == maxcoords[0]){
    rightcoords[0]=0;
    rightcoords[1]=coords[1];
    MPI_Cart_rank(cart_comm, rightcoords,&rightrank);
  }else{
    rightcoords[0]=coords[0]+1;
    rightcoords[1]=coords[1];
    MPI_Cart_rank(cart_comm, rightcoords,&rightrank);
  }
  if(coords[0] == 0){
    leftcoords[0]=maxcoords[0];
    leftcoords[1]=coords[1];
    MPI_Cart_rank(cart_comm, leftcoords,&leftrank);
  }else{
    leftcoords[0]=coords[0]-1;
    leftcoords[1]=coords[1];
    MPI_Cart_rank(cart_comm, leftcoords,&leftrank);
  }
  //siderank
  int s;
  if((coords[0]-1) >= 0){
    if((coords[1]-1) >= 0){
      s = 0;
      sidecoords[s][0]=coords[0]-1;
      sidecoords[s][1]=coords[1]-1;
      MPI_Cart_rank(cart_comm, sidecoords[s],&siderank[s]);
      int dl = calcIndex(width, 1, 1);
      sidecells[s] = field[dl];
    }
    if((coords[1]+1) <= maxcoords[1]){
      s = 1;
      sidecoords[s][0]=coords[0]-1;
      sidecoords[s][1]=coords[1]+1;
      MPI_Cart_rank(cart_comm, sidecoords[s],&siderank[s]);
      int ul = calcIndex(width, 1, height-2);
      sidecells[s] = field[ul];
    }
  }
  if((coords[0]+1) <= maxcoords[0]){
    if((coords[1]-1) >= 0){
      s = 2;
      sidecoords[s][0]=coords[0]+1;
      sidecoords[s][1]=coords[1]-1;
      MPI_Cart_rank(cart_comm, sidecoords[s],&siderank[s]);
      int dr = calcIndex(width, width-2, 1);
      sidecells[s] = field[dr];
    }
    if((coords[1]+1) <= maxcoords[1]){
      s = 3;
      sidecoords[s][0]=coords[0]+1;
      sidecoords[s][1]=coords[1]+1;
      MPI_Cart_rank(cart_comm, sidecoords[s],&siderank[s]);
      int ur = calcIndex(width, width-2, height-2);
      sidecells[s] = field[ur];
    }
  }
  //count number of sides
    int countside = 0;
    for(int h = 0; h < 4; h++){
      if(siderank[h] != num_tasks)
      countside++;
    }
  //send siderank
  MPI_Request request1[2*countside];
  MPI_Status status1[2*countside];
  int numrequest = 0;
  for(int h = 0; h < 4; h++){
    //printf("%d siderank %d num_tasks %d h %d \n",rank_cart,siderank[h], num_tasks,h);
    if(siderank[h] != num_tasks){
      MPI_Isend(&sidecells[h], 1, MPI_CHAR, siderank[h], 1, cart_comm, &(request1[numrequest]));
      MPI_Irecv(&recvsidecells[h], 1, MPI_CHAR, siderank[h], 1, cart_comm, &(request1[numrequest+1]));
      numrequest = numrequest +2;
      //printf("%d h %d numrequest %d \n",rank_cart,h, numrequest );
    }
  }
  //printf("%d before MPI_Waitall countside %d request1 %d\n",rank_cart,countside,(sizeof(request1)/sizeof(MPI_Request)));
  MPI_Waitall(countside, request1, status1);
  //printf("%d out\n",rank_cart );
  // put side cells in place
  int a1, a2, a3, a4;
  for(int h = 0; h < 4; h++){
   if(siderank[h] < num_tasks){
      switch (h) {
        case 0:
        a1 = calcIndex(width, 0, 0);
        field[a1] = recvsidecells[h];
        break;
        case 1:
        a2 = calcIndex(width, 0, height - 1);
        field[a2] = recvsidecells[h];
        break;
        case 2:
        a3 = calcIndex(width, width-1, 0);
        field[a3] = recvsidecells[h];
        break;
        case 3:
        a4 = calcIndex(width, width-1, height - 1);
        field[a4] = recvsidecells[h];
        break;
      }
    }
  }

  //prepare sendcells
  for (int x = 0; x < width - 1; x++) {
    int b = calcIndex(width, x, 1);
    int c = calcIndex(width, x, height - 2);
    sendcellstb[0][x] = field[c];
    sendcellstb[1][x] = field[b];
  }
  //mark the buffer top or bottom
  sendcellstb[0][width] = 't';
  sendcellstb[1][width] = 'b';
  for (int y = 0; y < height - 1; y++) {
      int j = calcIndex(width, 1, y);
      int k = calcIndex(width, width - 2, y);
      sendcellslr[0][y] = field[j];
      sendcellslr[1][y] = field[k];
  }
  //mark the buffer left or right
  sendcellslr[0][height] = 'l';
  sendcellslr[1][height] = 'r';
    printf("send %c empty receive %c \n",sendcellslr[0][height],recvcellslr[0][height]);
    MPI_Request request[8];
    MPI_Status status[8];
    MPI_Isend(sendcellstb[0], width+1, MPI_CHAR, toprank, 1, cart_comm, &(request[0]));
    MPI_Isend(sendcellstb[1], width+1, MPI_CHAR, botrank, 1, cart_comm, &(request[1]));
    MPI_Isend(sendcellslr[0], height+1, MPI_CHAR, leftrank, 1, cart_comm, &(request[2]));
    MPI_Isend(sendcellslr[1], height+1, MPI_CHAR, rightrank, 1, cart_comm, &(request[3]));

    MPI_Irecv(recvcellstb[0], width+1, MPI_CHAR, toprank, 1, cart_comm, &(request[4]));
    MPI_Irecv(recvcellstb[1], width+1, MPI_CHAR, botrank, 1, cart_comm, &(request[5]));
    MPI_Irecv(recvcellslr[0], height+1, MPI_CHAR, leftrank, 1, cart_comm, &(request[6]));
    MPI_Irecv(recvcellslr[1], height+1, MPI_CHAR, rightrank, 1, cart_comm, &(request[7]));
    MPI_Waitall(8, request, status);
    printf("send %c received %c \n",sendcellslr[0][height],recvcellslr[0][height]);
    printf("send %c received %c \n",sendcellslr[1][height],recvcellslr[1][height]);
  for(int i = 0; i < 2; i++){
    if(recvcellstb[i][width] == 'b'){
      for (int x = 0; x < width - 1; x++) {
        int a = calcIndex(width, x, height - 1);
        field[a] = recvcellstb[i][x];
      }
    }
    if(recvcellstb[i][width] == 't'){
      for (int x = 0; x < width - 1; x++) {
        int d = calcIndex(width, x, 0);
        field[d] = recvcellstb[i][x];
      }
    }
    if(recvcellslr[i][height] == 'r'){
      printf("%d right Received\n",rank_cart);
      for (int y = 0; y < height - 1; y++) {
          int l = calcIndex(width, 0, y);
          field[l] = recvcellslr[i][y];
      }
    }
    if(recvcellslr[i][height] == 'l'){
      printf("%d left Received\n", rank_cart);
      for (int y = 0; y < height - 1; y++) {
          int u = calcIndex(width, width - 1, y);
          field[u] = recvcellslr[i][y];
      }
    }
  }
}

void game (int width, int height, int num_timesteps, int gsizes[2]) {
  char *currentfield = calloc (width * height, sizeof(char));
  char *newfield = calloc (width * height, sizeof(char));
  //variables for switching boudaries
  //filling_random (currentfield, width, height);
  if(rank_cart == 0)
  filling_runner1 (currentfield, width, height);
  //filling_rank (currentfield, width, height);
  apply_periodic_boundaries(currentfield,width,height);
  int time = 0;
  write_field (currentfield, gsizes[X], gsizes[Y], time);
  for (time = 1; time <= num_timesteps; time++) {
    evolve (currentfield, newfield, width, height);
    write_field (newfield, gsizes[X], gsizes[Y], time);
    apply_periodic_boundaries(newfield,width,height);
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
  dims[0] = process_numX;
  dims[1] = process_numY;
  int periods[2] = {1,1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
  MPI_Comm_rank(cart_comm, &rank_cart);
  MPI_Cart_coords(cart_comm, rank_cart, 2, coords);
  gsizes[0] = width;  // global size of the domain without boundaries
  gsizes[1] = height;
  lsizes[0] = width/process_numX;
  lsizes[1] = height/process_numY;
  int starts[2] = {coords[0]*lsizes[0],coords[1]*lsizes[1]};
  int startindices[2] = {1,1};
  int memsize[2] = {lsizes[0]+2,lsizes[1]+2};
  double* elapsed_time_send = malloc(sizeof(double)*4);
  double* elapsed_time_recv = malloc(sizeof(double)*4);
  if(rank == 0){
  for(int i = 0; i < num_tasks; i++){
    elapsed_time_send[i] = 0;
    elapsed_time_recv[i] = 0;
    }
  }

  START_TIMEMEASUREMENT(measure_game_time);

  MPI_Type_create_subarray(2,gsizes,lsizes,starts, MPI_ORDER_FORTRAN,MPI_CHAR,&filetype);
  MPI_Type_commit(&filetype);

  MPI_Type_create_subarray(2,memsize,lsizes,startindices, MPI_ORDER_FORTRAN,MPI_CHAR,&memtype);
  MPI_Type_commit(&memtype);

  game(memsize[X], memsize[Y], num_timesteps, gsizes);

  END_TIMEMEASUREMENT(measure_game_time, elapsed_time_recv[rank]);

  //time calculate
  elapsed_time_send[rank] = elapsed_time_recv[rank];
  printf("rank: %d time elapsed: %lf sec\n", rank, elapsed_time_recv[rank]);
  int comm = 99;
    if (rank != 0 && rank < num_tasks-1) {
      // Receive from left worker
      MPI_Recv(elapsed_time_recv, num_tasks, MPI_DOUBLE,(rank-1), comm, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
      MPI_Recv(elapsed_time_recv, num_tasks, MPI_DOUBLE,(rank-1), comm, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
