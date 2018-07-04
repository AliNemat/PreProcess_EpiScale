//Algorithm.cpp
//==================================================
// Include Dependencies
#include <iostream>
#include <fstream>
#include <vector>

#include "phys.h"
#include "coord.h"
#include "algorithm.h"
//=======================


//========================================================
// Class Algorithm

Algorithm::Algorithm(double cell_r) {
	cell_radius = cell_r;
}

void Algorithm::output_Cell_Centers(string filename) {

	ofstream ofs(filename.c_str());

	ofs << cell_centers.size() << endl;
	Coord curr;

	for (unsigned int i = 0; i < cell_centers.size(); i++) {
		curr = cell_centers.at(i);
		ofs << curr.get_X() << ' ' << curr.get_Y() << ' ' << 0.0 << endl;
	}

	ofs.close();

	return;
}

//========================================================
//Class Rectangle Algorithm

Rect_Alg::Rect_Alg(double cell_r, double l, double w) : Algorithm(cell_r + eps) {
	length = l;
	width = w;

	Coord temp;
	// Vector for going up a level
	temp.set_X(2 * cell_radius * cos(pi / 3));
	temp.set_Y(2 * cell_radius * sin(pi / 3));
	up_vect = temp;
	// vector for moving to left
	temp.set_X((-2) * cell_radius);
	temp.set_Y(0.0);
	left_vect = temp;
	// vector for moving to right
	temp.set_X(2 * cell_radius);
	right_vect = temp;

}

void Rect_Alg::find_Cell_Centers() {

	Coord orig;
	Coord curr;

	while (within_Boundary(orig)) {
		
		cell_centers.push_back(orig);

		curr = orig + left_vect;
		//iterate to the left
		while (within_Boundary(curr)) {
			cell_centers.push_back(curr);
			curr += left_vect;
		}

		curr = orig + right_vect;
		//iterate to the right
		while (within_Boundary(curr)) {
			cell_centers.push_back(curr);
			curr += right_vect;
		}

		//move up a level
		orig += up_vect;
	}

	cout << "Number of cells added: " << cell_centers.size() << endl;

	return;
}

bool Rect_Alg::within_Boundary(Coord curr) {

	if ( (curr.get_X() < width) && (curr.get_X() >= 0.0) && 
		 (curr.get_Y() < length) && (curr.get_Y() >= 0.0) ) {
		return true;
	}
	else {
		return false;
	}
}

//=========================================================
// Class Circle Algorithm

Circle_Alg::Circle_Alg(double cell_r, double circle_r) : Algorithm (cell_r + eps) {
	
	circ_radius = circle_r;
	circ_center = Coord();

	Coord temp;
	// Vector for going up 
	temp.set_Y(2 * cell_radius);
	up_vect = temp;
	// vector for moving down
	temp.set_Y((-2) * cell_radius);
	down_vect = temp;
	// vectors for moving to right
	temp.set_X(2 * cell_radius * cos(11 * pi / 6));
	temp.set_Y(2 * cell_radius * sin(11 * pi / 6));
	right_down = temp;
	temp.set_X(2 * cell_radius * cos(pi / 6));
	temp.set_Y(2 * cell_radius * sin(pi / 6));
	right_up = temp;

	// vectors for moving left
	temp.set_X(2 * cell_radius * cos(7 * pi / 6));
	temp.set_Y(2 * cell_radius * sin(7 * pi / 6));
	left_down = temp;
	temp.set_X(2 * cell_radius * cos(5 * pi / 6));
	temp.set_Y(2 * cell_radius * sin(5 * pi / 6));
	left_up = temp;
	
}

void Circle_Alg::find_Cell_Centers() {

	Coord orig;
	Coord curr;

	move_Vertically(orig);
	
	//move right
	curr = orig + right_down;
	while (move_Vertically(curr)) {
		if (curr.get_Y() < 0.0) {
			curr += right_up;
		}
		else {
			curr += right_down;
		}
	}

	//move left
	curr = orig + left_down;
	while (move_Vertically(curr)) {
		if (curr.get_Y() < 0.0) {
			curr += left_up;
		}
		else {
			curr += left_down;
		}
	}

	cout << "Number of cells added: " << cell_centers.size() << endl;

	return;
}

bool Circle_Alg::move_Vertically(Coord orig) {
	
	if (! within_Boundary(orig)) {
		return false;
	}

	cell_centers.push_back(orig);

	//move up
	Coord curr = orig + up_vect;
	while ( within_Boundary(curr) ) {
		cell_centers.push_back(curr);
		curr += up_vect;
	}

	//move down
	curr = orig + down_vect;
	while ( within_Boundary(curr) ) {
		cell_centers.push_back(curr);
		curr += down_vect;
	}

	return true;
}

bool Circle_Alg::within_Boundary(Coord curr) {

	if ( curr.length() < circ_radius ) {
		return true;
	}
	else {
		return false;
	}
}




//===========================================================
//end of file

