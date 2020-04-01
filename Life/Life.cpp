#include <chrono>
#include <locale>
#include <sstream>
#include <fstream>
#include <string>
#include <mpi.h>

#include "boundaries.h"

using namespace std;

/*======================================================================================*/
// Set environment parameters
int imax_global = 108; int jmax_global = 113;		// Domain size in rows and columns
bool input_picture = true;						// Whether we run the game of life on a picture
bool test = false;								// Whether we run the game of life on a glider as a test
// Defaults to generating a random grid if input_picture and test are false
double init_alive = 40;							// Percentage of cells originally alive in random grid

bool periodic = true;							// Boundary conditions
int max_steps = 150;							// # generations to run game for
/*======================================================================================*/


/*======================================================================================*/
// Set parameters that the program requires (non-user defined)
int id, p;										// Stores processor ID and number of processors used to run program
int tag_num;									// Transaction ID for each send/recv pair
int proc_i, proc_j;								// Number of processors along i and j directions 
bool** grid;									// Declare double boolean pointer array that's going to act as our base grid (refreshed on each update)
bool** new_grid;								// Declare double boolean pointer array that's going to act as our update storage grid
int i_id, j_id;									// i and j position of processor in domain
int imax_local, jmax_local;						// Will store how many rows and columns each processor will be responsible for
bool* receive_list, * send_list;				// details which processor sends to/receive from each processor 
int gen;										// Stores current generation of game
/*======================================================================================*/


void domain_decomposition()
{
	// Peer-to-Peer communications:
	// Work out the optimal dimensions of the grid based on p
	// Do distribution calculations on zeroth processor
	// Change imax and jmax in here
	if (p > 1) {
		int gap = p;
		for (int n = 1; n < p; n++) {
			if (p % n == 0) {
				if (gap > abs(n - p / n)) {
					gap = abs(n - p / n);
					proc_i = n;
					proc_j = p / n;
				}
			}
		}
	}
	else {
		proc_i = 1;
		proc_j = 1;
	}

	i_id = id / proc_j;					// Col position of node
	j_id = id % proc_j;					// Row position of node

	// calculate the how many rows I am responsible for
	int rows_rem = imax_global;
	for (int i = 0; i <= i_id; i++)
	{
		// Calculate the number of rows I have
		imax_local = rows_rem / (proc_i - i);
		// Decrement the remaining rows
		rows_rem -= imax_local;
	}

	// calculate how many columns I am responsible for
	int cols_rem = jmax_global;
	for (int j = 0; j <= j_id; j++)
	{
		// Calculate the number of cols I have
		jmax_local = cols_rem / (proc_j - j);
		// Decrement the remaining cols
		cols_rem -= jmax_local;
	}
}

void create_communications_list(int* send_list, int* recv_list, bool& periodic, int& n_cnt, MPI_Request* requests)
{
	// Do send and receives for 8 neighbours
	n_cnt = 0;
	for (int ii = -1; ii <= 1; ii++) {
		for (int jj = -1; jj <= 1; jj++) {
			// Go back to original coordinates
			int i = ii + i_id;			// Convert to global coordinates
			int j = jj + j_id;			// Convert to global coordinates
			int neighbour_id;
			if (periodic) {
				// Allow it to flip over domain walls
				if (i != i_id || j != j_id) {
					// Then calculate the neighbour id using absolute values
					if (i < 0) {			// if neighbour coordinates negative in global coordinates
						i += proc_i;		// relocate to other side
					}
					if (i >= proc_i) {		// if  neighbour coordinates beyond scope of global coordinates
						i -= proc_i;		// relocate to other side
					}
					if (j < 0) {			// if neighbour coordinates negative in global coordinates
						j += proc_j;		// relocate to other side
					}
					if (j >= proc_j) {		// if neighbour coordinates beyond scope of global coordinates
						j -= proc_j;		// relocate to other side
					}
					neighbour_id = i * proc_j + j;
					send_list[n_cnt * 3] = id;
					send_list[n_cnt * 3 + 1] = -ii;		// Use relative (unflipped) values to infer direction of communication
					send_list[n_cnt * 3 + 2] = -jj;		// Use relative (unflipped) values to infer direction of communication
					// wrap up this information and send to every neighbour we have, and then receive same information
					// store [id|ii|jj]
					// start pointer ------¬
					MPI_Isend(&send_list[n_cnt * 3], 3, MPI_INT, neighbour_id, tag_num, MPI_COMM_WORLD, &requests[n_cnt * 2]);		// count is number of total sets of information
					MPI_Irecv(&recv_list[n_cnt * 3], 3, MPI_INT, neighbour_id, tag_num, MPI_COMM_WORLD, &requests[n_cnt * 2 + 1]);	// count only goes up once for every two communications (i.e. number of neighbours)
					n_cnt++;
				}
			}
			else if (!periodic) {
				// Don't allow it to cross over domain walls
				if ((i >= 0 && i < proc_i) && (j >= 0 && j < proc_j) && (i != i_id || j != j_id)) {
					neighbour_id = i * proc_j + j;		// translate to neighbour id
					send_list[n_cnt * 3] = id;			// add current id to send list
					send_list[n_cnt * 3 + 1] = -ii;		// Use relative (unflipped) values to infer direction of communication
					send_list[n_cnt * 3 + 2] = -jj;		// Use relative (unflipped) values to infer direction of communication
					// store [id| i| j |ii|jj]
					// start pointer ------¬
					MPI_Isend(&send_list[n_cnt * 3], 3, MPI_INT, neighbour_id, tag_num, MPI_COMM_WORLD, &requests[n_cnt * 2]);		// count is number of total sets of information
					MPI_Irecv(&recv_list[n_cnt * 3], 3, MPI_INT, neighbour_id, tag_num, MPI_COMM_WORLD, &requests[n_cnt * 2 + 1]);	// count only goes up once for every two communications (i.e. number of neighbours)
					n_cnt++;
				}
			}
		}
	}
	//MPI_Waitall(n_cnt * 2, requests, MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
}

void make_boundaries(bool**& grid, int* recv_list, int n_cnt, int tag_num)
{
	// Declare a class that extracts send/recv sections from the local grid
	boundary_communication section(grid, imax_local, jmax_local);
	// Generate MPI datatypes in make_types() method
	section.make_types();
	// Make requests list to bind our non-blocking communications
	MPI_Request* section_requests = new MPI_Request[n_cnt * 2 + 1];		// 2 sends and receives per neighbour

	// Communications count
	int comm_count = 0;

	// TAG NUMBER is written in such a way that the communication becomes unique for the generation AND the neighbour it is communicating with.
	for (int it = 0; it < n_cnt; it++) {							// Iterate over the number of neighbours we're assigned to communicate with

		// unpack information in communications list
		int partner = recv_list[it * 3];							// neighbouring processor id
		int i_rel_pos = recv_list[it * 3 + 1];						// +1 => neighbour is below of me
		int j_rel_pos = recv_list[it * 3 + 2];						// +1 => neighbour is to the right me

		// Left side
		if (i_rel_pos == 0 && j_rel_pos < 0) {
			// Simultaneously send left side and receive into left halo
			MPI_Irecv(MPI_BOTTOM, 1, section.left_recv_type, partner, gen+5, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.left_send_type, partner, gen+4, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Right side
		else if (i_rel_pos == 0 && j_rel_pos > 0) {
			// Simultaneously send right side and receive into right halo
			MPI_Irecv(MPI_BOTTOM, 1, section.right_recv_type, partner, gen+4, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.right_send_type, partner, gen+5, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Top side
		else if (i_rel_pos < 0 && j_rel_pos == 0) {
			// Simultaneously send top side and receive into top halo
			MPI_Irecv(MPI_BOTTOM, 1, section.top_recv_type, partner, gen+7, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.top_send_type, partner, gen+2, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Bottom side
		else if (i_rel_pos > 0 && j_rel_pos == 0) {
			// Simultaneously send bottom side and receive into bottom halo
			MPI_Irecv(MPI_BOTTOM, 1, section.bottom_recv_type, partner, gen+2, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.bottom_send_type, partner, gen+7, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Top left corner
		else if (i_rel_pos < 0 && j_rel_pos < 0) {
			// Simultaneously send top left corner and receive into top left corner of halo
			MPI_Irecv(MPI_BOTTOM, 1, section.top_left_recv_type, partner, gen+8, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.top_left_send_type, partner, gen+1, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Top right corner
		else if (i_rel_pos < 0 && j_rel_pos > 0) {
			// Simultaneously send top-rigt corner and receive into top-right corner of halo
			MPI_Irecv(MPI_BOTTOM, 1, section.top_right_recv_type, partner, gen+6, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.top_right_send_type, partner, gen+3, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Bottom left corner
		else if (i_rel_pos > 0 && j_rel_pos < 0) {
			// Simultaneously send bottom-left corner and receive into bottom-left corner of halo
			MPI_Irecv(MPI_BOTTOM, 1, section.bottom_left_recv_type, partner, gen+3, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.bottom_left_send_type, partner, gen+6, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}

		// Bottom right corner
		else if (i_rel_pos > 0 && j_rel_pos > 0) {
			// Simultaneously send bottom right corner and receive into bottom right corner of halo
			MPI_Irecv(MPI_BOTTOM, 1, section.bottom_right_recv_type, partner, gen+1, MPI_COMM_WORLD, &section_requests[comm_count * 2]);
			MPI_Isend(MPI_BOTTOM, 1, section.bottom_right_send_type, partner, gen+8, MPI_COMM_WORLD, &section_requests[comm_count * 2 + 1]);
			comm_count++;
		}
	}

	MPI_Waitall(2 * n_cnt, section_requests, MPI_STATUSES_IGNORE) == MPI_SUCCESS;
	delete[] section_requests;		// Relinquish memory associated with section_requests
}

int cell_neighbours(bool**& grid, int x, int y)
{
	// Calculate neighbours surrounding current cell at x (i) and y (j) position
	// This taps into the halo, so we needn't worry about boundary conditions
	// Initialize neighbour counts with zero
	int neighbour_count = 0;
	for (int i = -1; i <= +1; i++) {
		for (int j = -1; j <= +1; j++) {
			// Look round neighbourhood and add up how many alive neighbours there are
			if (i == 0 && j == 0) {}									// Don't count self
			else { neighbour_count += grid[x + i][y + j]; }				// Append to sum	
		}
	}
	// Return how many cells that are currently alive near me
	return neighbour_count;
}

int game_rules(int n_neighbours, bool state)
{
	/* Rules:
		1 - if a living cell has fewer than 2 neighbours, it dies from inability to breed;
		2 - if a living cell has 2-3 neighbours, it survives;
		3 - if a living cell has 4 or more neighbours, it dies from overpopulation, and;
		4 - if a dead cell has exactly 3 neighbours, it respawns.
	*/
	if (state == 1) {									// already alive
		if (n_neighbours < 2)							// Rule 1: If already alive with less than two neighbours => dies
			return 0;
		if (n_neighbours == 2 || n_neighbours == 3)		// Rule 2: If already alive with two/three neighbours => survives
			return 1;
		if (n_neighbours >= 4)							// Rule 3: If already alive with at least four neighbours => dies
			return 0;
	}
	else {												// already dead
		if (n_neighbours == 3)							// Rule 4: If dead with exactly three neighbours => revives
			return 1;
	}
	// Return whether the current cell is alive (1) or dead (0)
	return state;
}

void game_update(bool**& grid, bool periodic, int* recv_list, int n_cnt, int tag_num)
{
	// Add in the pads
	make_boundaries(grid, recv_list, n_cnt, tag_num);

	// Do calculations on the current grid and paint them onto the new canvas grid
	for (int i = 1; i < imax_local + 1; i++) {
		for (int j = 1; j < jmax_local + 1; j++) {
			new_grid[i][j] = game_rules(cell_neighbours(grid, i, j), grid[i][j]);
		}
	}
	tag_num++;			// Increment the tag number by one to prevent race conditions

	// Swap the old and new grids to complete the update
	bool** temp = grid;
	grid = new_grid;
	new_grid = temp;
}

void generate_random_grid(int imax_local, int jmax_local, double init_alive, bool**& grid)
{
	// Resize the grid without having to use loop
	grid = new bool* [imax_local + 2];		// Include pads at boundaries
	new_grid = new bool* [imax_local + 2];

	// Loop over all elements and assign random bool
	for (int i = 0; i < imax_local + 2; i++) {
		grid[i] = new bool[jmax_local + 2];
		new_grid[i] = new bool[jmax_local + 2];

		for (int j = 0; j < jmax_local + 2; j++) {
			grid[i][j] = (bool)0;
			new_grid[i][j] = (bool)0;
		}
	}
	// Loop over rows and columns and set to random bool
	for (int i = 1; i <= imax_local; i++) {
		for (int j = 1; j <= jmax_local; j++) {
			grid[i][j] = (bool)((rand() % 100) < init_alive);
		}
	}
}

void write_in(int imax_local, int jmax_local, bool**& grid)
{
	// Make arrow shape to run test
	grid = new bool* [imax_local + 2];		// Include pads at boundaries
	new_grid = new bool* [imax_local + 2];

	// Initially fill in with zeros
	for (int i = 0; i < imax_local + 2; i++) {
		grid[i] = new bool[jmax_local + 2];
		new_grid[i] = new bool[jmax_local + 2];
		for (int j = 0; j < jmax_local + 2; j++) {
			grid[i][j] = 0;
			new_grid[i][j] = 0;
		}
	}

	// Then manually fill in values for a glider
	for (int j = 3; j < 6; j++) {
		grid[3][j] = 1;
	}
	grid[1][4] = 1;
	grid[2][5] = 1;
}

void get_picture(bool**& grid)
{
	// Preallocate memory for grids
	grid = new bool* [imax_local + 2];		// Include pads at boundaries
	new_grid = new bool* [imax_local + 2];	

	// Initially fill in with zeros
	for (int i = 0; i < imax_local + 2; i++) {
		grid[i] = new bool[jmax_local + 2];
		new_grid[i] = new bool[jmax_local + 2];
		for (int j = 0; j < jmax_local + 2; j++) {
			grid[i][j] = 0;
			new_grid[i][j] = 0;
		}
	}
	// Read in file specified during pre-processing
	ifstream infile;
	infile.open("./indata/grid_" + to_string(id) + '_' + to_string(i_id) + '_' + to_string(j_id) + ".txt");
	// Check that file exists
	if (!infile) {
		cerr << "Unable to open file datafile";
		exit(1);   // call system to stop
	}

	string x;
	int cnt = 0;
	// Act only if file exists
	while (infile.good()) {
		getline(infile, x, ',');
		bool b = x == "1";									// Convert string to bool
		grid[cnt / jmax_local][cnt % jmax_local] = b;		// Fill in the code
		cnt++;												// Increment count
	}
	infile.close();
}

void write_out_meta(double time)
{
	// Print out information to be used by user and post-processor
	/*Information in each file:
	- Filename: [i_id (processor row position in domain)]_[j_id (processor col position in domain)]_[generation]_[# cores]_info.txt
	- imax_local
	- jmax local
	- runtime. */
	stringstream filename;		// Declare stream
	filename << "./meta/" << i_id << '_' << j_id << '_' << gen - 1 << '_' << p << "_info.txt";
	fstream f;					// Declare file
	f.open(filename.str().c_str(), ios_base::out);	// Open to be written into
	f << imax_local << " ";							// rows that the processor is responsible for
	f << jmax_local << " ";							// cols that the processor is responsible for
	f << time << " ";								// runtime
	f.close();
}

void write_out(bool**& grid, int it)
{
	// Write out each generation state for each timestep on each processor:
	// Format of filename: "./data/`i coordinates`_`j coordinates`_`generation`.txt"
	stringstream filename;
	filename << "./data/" << i_id << '_' << j_id << '_' << it << ".txt";
	fstream f;
	f.open(filename.str().c_str(), ios_base::out);
	// Written in row-major ravelled format
	// 1 2 3
	// 4 5 6  =>  1 2 3 4 5 6 7 8 9
	// 7 8 9

	// such that
	// 1 0 1
	// 1 1 1  =>  1 1 1 1 1 1 0 1 1
	// 0 1 1

	// Write out into file
	for (int i = 1; i <= imax_local; i++) {
		for (int j = 1; j <= jmax_local; j++) {
			f << grid[i][j] << " ";
		}
	}
	f.close();
}

int main(int argc, char* argv[])
{
	// Initialize environment
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 10);
	int tag_num = 0;

	// Peer-to-Peer communications:
	/*Outcomes:
	- i_id, j_id calculated;
	- proc_i, proc_j calculated;
	- imax_local, jmax_local calculated.*/
	domain_decomposition();

	// Inform user of regional affairs
	if (id == 0) {
		cout << "Dividing domain into (" << proc_i << " x " << proc_j << ") blocks.\n\n";
		cout.flush();
	}

	MPI_Request* requests = new MPI_Request[16];	// Maximum number of requests for sending and receiving from 8 surrounding elements
	int* send_list = new int[8 * 3];				// 8 potential communications with 3 actual sends
	int* recv_list = new int[8 * 3];				// 8 potential communications with 3 actual sends
	int n_cnt;										// Neighbour count on this processor

	// Before we begin, we must equip every processor with information regarding
	// [ id | i | j |rel_i|rel_j]
	// id - identification of what processor we are communicating with;
	// rel_i - relative (to us) i index of processor we are communicating with
	// rel_j - relative (to us) j index of processor we are communicating with
	create_communications_list(send_list, recv_list, periodic, n_cnt, requests);
	tag_num++;

	// Inform user of who is communicated with whom
	cout << "Process " << id << " received data from " << n_cnt << " processors : \n";
	for (int n = 0; n < n_cnt; n++) {
		cout << recv_list[n * 3] << " in relative direction i, j = (" << recv_list[n * 3 + 1] << "," << recv_list[n * 3 + 2] << "). \n";
	}

	// Generate a grid which each processor needs to solve
	if (test) write_in(imax_local, jmax_local, grid);
	else if (input_picture) get_picture(grid);
	else (generate_random_grid(imax_local, jmax_local, init_alive, grid));		// Returns ((1 + len_x + 1) * (1 + len_y + 1)))
	// Write out first iteration to file
	write_out(grid, 0);

	// For time tests, we will keep the barrier in place to make sure all processors have aligned before executing time measurements
	MPI_Barrier(MPI_COMM_WORLD);
	// Declare start clock
	chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

	// Start timestepping through game generations
	for (gen = 1; gen <= max_steps; gen++) {
		// Update the grid for each iteration of the game
		game_update(grid, periodic, recv_list, n_cnt, tag_num);
		// Write out each grid's contents to file to be used by the post-processor
		write_out(grid, gen);
	}

	// Impose another barrier to make sure all processors are at the same point such that we are now measuring the slowest processor
	MPI_Barrier(MPI_COMM_WORLD);
	// Declare end clock
	chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
	// Calculate duration
	chrono::duration<double> time_it = chrono::duration_cast<chrono::duration<double>>(end - start);
	// Print out information to be used by user and post-processor
	/*Information in each file:
	- Filename: [i_id (processor row position in domain)]_[j_id (processor col position in domain)]_[generation]_[# cores]_info.txt
	- imax_local
	- jmax local
	- runtime.*/
	write_out_meta(time_it.count());

	for (int i = 0; i < imax_local + 2; i++) {
		delete[] grid[i];
		delete[] new_grid[i];
	}
	delete[] grid;
	delete[] new_grid;

	delete[] requests;
	delete[] send_list;
	delete[] recv_list;

	MPI_Finalize();
}
