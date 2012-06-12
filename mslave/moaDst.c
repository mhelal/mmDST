#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <errno.h>
#include "moaDst.h"
#include "globals.h"
#include "utils.h"

void checkMutexErrorCode (int caller, int errCode) {
  
  switch (errCode) {
  case EINVAL:
    mprintf(1, "the mutex has not been properly initialized.\n", caller);
    break;
  case EDEADLK:
    mprintf(1, "the  mutex  is  already  locked  by  the  calling  thread (``error checking'' mutexes only).\n", caller);
    break;
  case EBUSY:
    mprintf(1, "the mutex could not be acquired because it was  currently locked.\n", caller);
    break;
  case EPERM:
    mprintf(1, "the calling thread does not own the mutex (``error checking'' mutexes only).\n", caller);
    break;
  }
  
}

void checkMPIErrorCode (int caller, int errmsg) {
  switch(errmsg) {
  case MPI_SUCCESS:
    //mprintf(1, "Proc [%d]success\n", 0, myProcid);
    break;
  case MPI_ERR_COMM:
    mprintf(1, "Invalid communicator\n", caller);
    break;
  case MPI_ERR_TYPE:
    mprintf(1, "Invalid datatype argument\n", caller);
    break;
  case MPI_ERR_COUNT:
    mprintf(1, "Invalid count argument\n", caller);
    break;
  case MPI_ERR_TAG:
    mprintf(1, "Invalid  tag  argument\n", caller);
    break;
  case MPI_ERR_RANK:
    mprintf(1, "Invalid source or destination rank\n", caller);
    break;
  }
}

int NBReceive (void * buffer, int count, MPI_Datatype type, int source, int tag, MPI_Request * request, MPI_Status * status) {
	struct timeval timeout;
	int flag, MPI_return;
	char msg[SHORT_MESSAGE_SIZE];

		
	MPI_return = MPI_Recv (buffer, count, type, source, tag, MOAMSA_COMM_WORLD, status);
	checkMPIErrorCode (2, MPI_return);
	/*
	MPI_return = MPI_Irecv (buffer, count, type, source, tag, MOAMSA_COMM_WORLD, request);
	checkMPIErrorCode (MPI_return);
	MPI_Test(request, &flag, status);
	while (flag==0) {
		timeout.tv_sec = 0;
		timeout.tv_usec = 0;
		select (10, NULL, NULL, NULL, &timeout);
		MPI_Test(request, &flag, status);
	}
*/
	return MPI_return;
}

int NBSend (void * buffer, int count, MPI_Datatype type, int destination, int tag, MPI_Request * request) {
 
	int MPI_return;

  MPI_return = MPI_Send (buffer, count, type, destination, tag, MOAMSA_COMM_WORLD);
  checkMPIErrorCode (3, MPI_return);
  return MPI_return;
}
