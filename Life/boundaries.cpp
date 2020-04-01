#include "boundaries.h"

using namespace std;

// Constructor
boundary_communication::boundary_communication(bool**& grid, int imax_local, int jmax_local) : grid(grid), imax(imax_local), jmax(jmax_local)
{}

void boundary_communication::make_types()
{
	// Declare storage vectors
	vector<int> block_length;
	vector<MPI_Aint> addresses;
	vector<MPI_Datatype> typelist;

	/*========================================================================================*/
	// LEFT
	/*========================================================================================*/
	// Data is discontiguous in memory space, so we must extract the entire column iteratively
	//Left side send type
	for (int i = 1; i <= imax; i++)
	{
		block_length.push_back(1);
		MPI_Aint ls_temp;
		MPI_Get_address(&grid[i][1], &ls_temp);
		addresses.push_back(ls_temp);
		typelist.push_back(MPI_C_BOOL);
	}
	// Create datatype for sending left strip of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->left_send_type);
	MPI_Type_commit(&this->left_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();
	// Data is discontiguous in memory space, so we must extract the entire column iteratively
	//Left side recv type
	for (int i = 1; i <= imax; i++)
	{
		block_length.push_back(1);
		MPI_Aint lr_temp;
		MPI_Get_address(&grid[i][0], &lr_temp);
		addresses.push_back(lr_temp);
		typelist.push_back(MPI_C_BOOL);
	}
	// Create datatype for sending left strip of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->left_recv_type);
	MPI_Type_commit(&this->left_recv_type);

	/*========================================================================================*/
	// RIGHT
	/*========================================================================================*/
	
	block_length.clear();
	addresses.clear();
	typelist.clear();
	// Data is discontiguous in memory space, so we must extract the entire column iteratively
	// Right side send type
	for (int i = 1; i <= imax; i++)
	{
		block_length.push_back(1);
		MPI_Aint rs_temp;
		MPI_Get_address(&grid[i][jmax], &rs_temp);
		addresses.push_back(rs_temp);
		typelist.push_back(MPI_C_BOOL);
	}
	// Create datatype for sending right strip of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->right_send_type);
	MPI_Type_commit(&this->right_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();
	// Data is discontiguous in memory space, so we must extract the entire column iteratively
	// Right side recv type
	for (int i = 1; i <= imax; i++)
	{
		block_length.push_back(1);
		MPI_Aint rr_temp;
		MPI_Get_address(&grid[i][jmax + 1], &rr_temp);
		addresses.push_back(rr_temp);
		typelist.push_back(MPI_C_BOOL);
	}
	// Create datatype for sending right strip of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->right_recv_type);
	MPI_Type_commit(&this->right_recv_type);

	/*========================================================================================*/
	// TOP
	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();
	MPI_Aint ts_temp;
	// Top side send type
	// Data is contiguous in space, so we can just grab the first pointer in that array
	block_length.push_back(jmax);
	MPI_Get_address(&grid[1][1], &ts_temp);
	addresses.push_back(ts_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending top strip of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->top_send_type);
	MPI_Type_commit(&this->top_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();
	MPI_Aint tr_temp;
	// Top side recv type
	// Data is contiguous in space, so we can just grab the first pointer in that array
	block_length.push_back(jmax);
	MPI_Get_address(&grid[0][1], &tr_temp);
	addresses.push_back(tr_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending top strip of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &top_recv_type);
	MPI_Type_commit(&this->top_recv_type);

	/*========================================================================================*/
	// BOTTOM
	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	// bottom side send type
	// Data is contiguous in space, so we can just grab the first pointer in that array
	block_length.push_back(jmax);
	MPI_Aint bs_temp;
	MPI_Get_address(&grid[imax][1], &bs_temp);
	addresses.push_back(bs_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending bottom strip of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->bottom_send_type);
	MPI_Type_commit(&this->bottom_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//top side recv type
	block_length.push_back(jmax);
	MPI_Aint br_temp;
	MPI_Get_address(&grid[imax + 1][1], &br_temp);
	addresses.push_back(br_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending bottom strip of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->bottom_recv_type);
	MPI_Type_commit(&this->bottom_recv_type);


	/*========================================================================================*/
	// TOP - LEFT
	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//top-left corner send type
	block_length.push_back(1);
	MPI_Aint tls_temp;
	MPI_Get_address(&grid[1][1], &tls_temp);
	addresses.push_back(tls_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending top-left corner of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->top_left_send_type);
	MPI_Type_commit(&this->top_left_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//top-left corner recv type
	block_length.push_back(1);
	MPI_Aint tlr_temp;
	MPI_Get_address(&grid[0][0], &tlr_temp);
	addresses.push_back(tlr_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending top - left corner of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->top_left_recv_type);
	MPI_Type_commit(&this->top_left_recv_type);

	/*========================================================================================*/
	// TOP - RIGHT
	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//top-right corner send type
	block_length.push_back(1);
	MPI_Aint trs_temp;
	MPI_Get_address(&grid[1][jmax], &trs_temp);
	addresses.push_back(trs_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending top - right corner of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->top_right_send_type);
	MPI_Type_commit(&this->top_right_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//top-right corner recv type
	block_length.push_back(1);
	MPI_Aint trr_temp;
	MPI_Get_address(&grid[0][jmax + 1], &trr_temp);
	addresses.push_back(trr_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending top - right corner of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->top_right_recv_type);
	MPI_Type_commit(&this->top_right_recv_type);

	/*========================================================================================*/
	// BOTTOM - LEFT
	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	// bottom-left corner send type
	block_length.push_back(1);
	MPI_Aint bls_temp;
	MPI_Get_address(&grid[imax][1], &bls_temp);
	addresses.push_back(bls_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending bottom - left corner of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->bottom_left_send_type);
	MPI_Type_commit(&this->bottom_left_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	// bottom-left corner recv type
	block_length.push_back(1);
	MPI_Aint blr_temp;
	MPI_Get_address(&grid[imax + 1][0], &blr_temp);
	addresses.push_back(blr_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending bottom - left corner of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->bottom_left_recv_type);
	MPI_Type_commit(&this->bottom_left_recv_type);

	/*========================================================================================*/
	// BOTTOM - RIGHT
	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//bottom-right corner send type
	block_length.push_back(1);
	MPI_Aint brs_temp;
	MPI_Get_address(&grid[imax][jmax], &brs_temp);
	addresses.push_back(brs_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending bottom - right corner of inner grid
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->bottom_right_send_type);
	MPI_Type_commit(&this->bottom_right_send_type);

	/*========================================================================================*/

	block_length.clear();
	addresses.clear();
	typelist.clear();

	//Bottom-right corner recv type
	block_length.push_back(1);
	MPI_Aint brr_temp;
	MPI_Get_address(&grid[imax + 1][jmax + 1], &brr_temp);
	addresses.push_back(brr_temp);
	typelist.push_back(MPI_C_BOOL);
	// Create datatype for sending bottom - right corner of halo
	MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->bottom_right_recv_type);
	MPI_Type_commit(&this->bottom_right_recv_type);

	block_length.clear();
	addresses.clear();
	typelist.clear();
}


boundary_communication::~boundary_communication() 
{
	// Free the newly created datatypes when they're no longer needed
	MPI_Type_free(&this->left_send_type);
	MPI_Type_free(&this->left_recv_type);
	MPI_Type_free(&this->right_send_type);
	MPI_Type_free(&this->right_recv_type);
	MPI_Type_free(&this->top_send_type);
	MPI_Type_free(&this->top_recv_type);
	MPI_Type_free(&this->bottom_send_type);
	MPI_Type_free(&this->bottom_recv_type);
	MPI_Type_free(&this->top_left_send_type);
	MPI_Type_free(&this->top_left_recv_type);
	MPI_Type_free(&this->top_right_send_type);
	MPI_Type_free(&this->top_right_recv_type);
	MPI_Type_free(&this->bottom_left_send_type);
	MPI_Type_free(&this->bottom_left_recv_type);
	MPI_Type_free(&this->bottom_right_send_type);
	MPI_Type_free(&this->bottom_right_recv_type);
};

