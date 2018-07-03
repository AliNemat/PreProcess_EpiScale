// Algorithm.h
//=============================================
// Include Guards
#ifndef _ALG_H_INCLUDED_
#define _ALG_H_INCLUDED_
//=====================
// Forward Declarations

//=====================
// Include Dependencies
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

#include "coord.h"
//=====================



//===============================================
// Class Algorithm

class Algorithm {

	protected:
		vector<Coord> cell_centers;
		double cell_radius;

	public:
		Algorithm(double cell_r);

		virtual void find_Cell_Centers() = 0;
		virtual void output_Cell_Centers(string filename);
		virtual bool within_Boundary(Coord curr) = 0;

};


//=================================================
// Class Square Algorithm

class Rect_Alg : public Algorithm {

	private:
		double length;
		double width;
		Coord up_vect;
		Coord left_vect;
		Coord right_vect;

	public:
		Rect_Alg(double cell_r, double l, double w);
	
		virtual void find_Cell_Centers();
		bool within_Boundary(Coord curr);

};

//=================================================
// Class Circle Algorithm



class Circle_Alg : public Algorithm {

	private:
		double circ_radius;
		Coord circ_center;

		Coord up_vect;
		Coord down_vect;
		Coord left_up;
		Coord left_down;
		Coord right_up;
		Coord right_down;


	public:
		Circle_Alg(double cell_r, double circle_r);

		virtual void find_Cell_Centers();
		bool move_Vertically(Coord orig);
		bool within_Boundary(Coord curr);

};




//=====================================================
//end of file

#endif



