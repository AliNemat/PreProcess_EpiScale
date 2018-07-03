//main.cpp
//===================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "coord.h"
#include "algorithm.h"

using namespace std;


//=============================================
// Declare Helper Functions
//=============================================

Algorithm* which_Algorithm();

//=============================================
// Main Function
//=============================================
int main(int argc, char** argv) {

	if ( argc != 2) {
		cout << "Error: incorrect number of arguments" << endl;
		cout << "Usage: Executable   outputfilename" << endl;
		return 1;
	}

	string filename = argv[1];

	//Ask user about their structure
	Algorithm* alg = which_Algorithm();

	if (alg == NULL) {
		return 0;
	}
	cout << "created alg" << endl;
	//Find Cell Centers
	alg->find_Cell_Centers();
	cout << "found cell centers" << endl;
	//Output Cell Centers to Outpufile
	alg->output_Cell_Centers(filename);
	

	return 0;
}



//===============================================
// Define Helper Functions
//===============================================

Algorithm* which_Algorithm() {

	Algorithm* alg;
	string type;
	double len, width, circ_rad, cell_rad;

	cout << "Cell Radius: ";
	cin >> cell_rad;
	cout << endl;

	cout << "What shape would you like your cell layer?" << endl;
	cout << "Options: 'Rectangle' or 'Circle'" << endl;

	cin >> type;

	if (type == "Rectangle") {
		cout << "Length of rectangle: ";
		cin >> len;
		cout << endl;
		cout << "Width of rectangle: ";
		cin >> width;
		cout << endl;
		alg = new Rect_Alg(cell_rad, len, width);
	}
	else if (type == "Circle") {
		cout << "Radius of Circle: ";
		cin >> circ_rad;
		cout << endl;
		alg = new Circle_Alg(cell_rad, circ_rad);
	}
	else {
		cout << "Error: No shape of that type." << endl;
		return NULL;
	}

	return alg;
}







//=============================================
//end of file

