#pragma once
#include <iostream>
#include <vector>
#include <mpi.h>
using namespace std;

class boundary_communication
{
public:
	// Constructor
	boundary_communication(bool**& grid, int imax_local, int jmax_local);

	// Constructs MPI datatypes
	void make_types();

	int imax;
	int jmax;

	bool** grid;

	// Declare different MPI types to be created and activated
	MPI_Datatype left_send_type;
	MPI_Datatype left_recv_type;
	MPI_Datatype right_send_type;
	MPI_Datatype right_recv_type;
	MPI_Datatype top_send_type;
	MPI_Datatype top_recv_type;
	MPI_Datatype bottom_send_type;
	MPI_Datatype bottom_recv_type;
	MPI_Datatype top_left_send_type;
	MPI_Datatype top_left_recv_type;
	MPI_Datatype top_right_send_type;
	MPI_Datatype top_right_recv_type;
	MPI_Datatype bottom_left_send_type;
	MPI_Datatype bottom_left_recv_type;
	MPI_Datatype bottom_right_send_type;
	MPI_Datatype bottom_right_recv_type;

	~boundary_communication();
};
