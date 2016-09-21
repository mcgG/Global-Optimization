#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cmath>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Circle_2<K> Circle;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef Delaunay::Vertex_handle Vertec_handle;
typedef K::Point_2 Point;

std::vector<Point> vertices; // vector stores all vertices of delaunay triangles
Delaunay dt;
Point Vn_Point0, Vn_Point1, Vn_Point2, Beta_Point0, Beta_Point1, Beta_Point2, C;

int i;	//initial number of Delaunay triangles 
int j;	//initial number of iterations
//double minMn = ; //smallest function value
double preVn; //Vn-1


//This is the quadratic function 
double func(double x1, double x2) {
	double result = pow(x1+x2-1, 2.0) + 10*pow(x1-x2, 2.0);
	return result;
}

//  Branin function
// double func(double x1, double x2) {
// 	x1 = 15.0*x1-5.0; //rescale to [-5,10]x[0,15]
//   x2 = 15.0*x2;
// 	double result = pow(x2-5.1/4.0/pow(M_PI, 2.0)*pow(x1, 2.0)+5.0*x1/M_PI-6.0, 2.0) + 10*(1.0-1.0/8.0/M_PI)*cos(x1) + 10.0;
// 	return result;
// }

//  Branin function
// double func(double x1, double x2) {
//     x1 = 15.0*x1-5.0; //rescale to [-5,10]x[0,15]
//     x2 = 15.0*x2;
//     double z1 = x2-(5.1/(4.0*pow(M_PI, 2.0)))*pow(x1, 2.0)+5.0*x1/M_PI-6.0;
//     double z2 = 10.0*(1.0-1.0/(8.0*M_PI))*cos(x1);
//     double result = pow(z1, 2.0) + z2 +10.0;
//     return result;
// }

// Rastrigin function
// double func(double x1, double x2) {
// 	return x1*x1 + x2*x2 - cos(18*x1) - cos(18*x2);
// }

//Initial the delaunay triangle
void init_Triangle() {
	i = 2;
	j = 0;
	dt.insert(vertices.begin(), vertices.end());
}

//The average of the function values at the three vertices
double getFi(Point p0, Point p1, Point p2) {
	return ( func( p0.hx(), p0.hy() )
				 + func( p1.hx(), p1.hy() )
				 + func( p2.hx(), p2.hy() ) )
				 / 3.0;
}

// Get the smallest function value among all vertices
double getMn() {
	double Mn = std::numeric_limits<double>::max();
	for (Point p : vertices) {
		double currFunc = func(p.hx(), p.hy());
		if (currFunc < Mn)	Mn = currFunc;
	}
	return Mn;
}

// Get the smallest area among all delaunay facets
double getVn() {
	double Vn = std::numeric_limits<double>::max();
	Delaunay::Finite_faces_iterator fit;
	for (fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++) {
		Point p0 = Point(fit->vertex(0)->point().hx(), fit->vertex(0)->point().hy());  
  	Point p1 = Point(fit->vertex(1)->point().hx(), fit->vertex(1)->point().hy());  
  	Point p2 = Point(fit->vertex(2)->point().hx(), fit->vertex(2)->point().hy()); 
   	double currArea = CGAL::area(p0, p1, p2);
  	if (currArea < Vn) {
  		Vn = currArea;
  		Vn_Point0 = p0;
  		Vn_Point1 = p1;
  		Vn_Point2 = p2;
  	}
	}
	return Vn;
}

//Get g(Vn) value
double getGVn(double Vn) {
	double q = Circle(Vn_Point0, Vn_Point1, Vn_Point2).squared_radius() 
						/ CGAL::area(Vn_Point0, Vn_Point1, Vn_Point2);
	if (Vn > 0 && Vn <= 0.5)	return q * Vn * log(1.0 / Vn);
	return q;
}

/**
*	Get the best value(largest) of Phi
*	Beta_Point0-2: vertex of the best triangle's centroid point
* C is the centroid point
*/
double getBestPhi(double Mn, double GVn) {
	double BestPhi = std::numeric_limits<double>::min();
	Delaunay::Finite_faces_iterator fit;
	for (fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++) {
			Point p0 = Point(fit->vertex(0)->point().hx(), fit->vertex(0)->point().hy());  
	  	Point p1 = Point(fit->vertex(1)->point().hx(), fit->vertex(1)->point().hy());  
	  	Point p2 = Point(fit->vertex(2)->point().hx(), fit->vertex(2)->point().hy()); 
	   	double T = CGAL::area(p0, p1, p2);
	   	//double currPhi = T / (getFi(p0, p1, p2) - Mn + GVn);
	   	double currPhi = T / (getFi(p0, p1, p2) - Mn);
	   	if (currPhi > BestPhi)	{
	   		BestPhi = currPhi;
	   		Beta_Point0 = p0;	
	   		Beta_Point1 = p1;
	   		Beta_Point2 = p2;
	   		C = CGAL::centroid(p0, p1, p2);
	   	}
		}
	return BestPhi;
}

// Rebuild the delaunay triangle, insert Point C into the Delaunay mesh
void rebuildDelaunay() {
	vertices.push_back(C);
	dt.insert(C);
}


//main function
int main(int argc, char *argv[]) {
	//get the iteration number
	int iteration = std::stoi(argv[1]);
	//initial points
	// Point p0 = Point(-2, -2);
	// Point p1 = Point(2, -2);
	// Point p2 = Point(2, 2);
	// Point p3 = Point(-2, 2);

	Point p0 = Point(0, 0);
	Point p1 = Point(0, 1);
	Point p2 = Point(1, 1);
	Point p3 = Point(1, 0);
	vertices.push_back(p0);
	vertices.push_back(p1);
	vertices.push_back(p2);
	vertices.push_back(p3);

	//Initial delaunay tirangle
	init_Triangle();


	double Mn = getMn();
	double Vn = getVn();
	double GVn = getGVn(Vn);

	preVn = Vn;


	while (j < iteration) {
	 	double BestPhi = getBestPhi(Mn, GVn); // get the best Phi value
		rebuildDelaunay();										// rebuild the delaunay triangle 
		double C_Value = func(C.hx(), C.hy());// evaluate f(c)
		if (C_Value < Mn) Mn = C_Value;				// if f(c) is smaller than Mn(the smallest value among all vertices than Mn=f(c))
		
		Vn = getVn();		
		GVn = getGVn(Vn);											// Among all triangle facets calculate the smallest area
		// if (Vn < preVn) {											// If the smallest area Vn is smaller than the previous Vn, then update g(Vn)
		// 	GVn = getGVn(Vn);
		// 	preVn = Vn;
		// } else {
		// 	Vn = preVn;													// If the smallest area Vn is not smaller than the previous Vn, do nothing
		// }
		j++;																	//	increase iteratoin number
	}

	//This is the output part
	std::ofstream myfile;
	myfile.open("Vertices.txt");
	for (Point p : vertices) {
		myfile << p << ";";
		std::cout << "value: " << func(p.hx(), p.hy()) << std::endl;
	}	
	myfile.close();

	return 0;
}