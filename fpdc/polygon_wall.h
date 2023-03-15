#ifndef __POLYGON_WALL__H__
#define __POLYGON_WALL__H__

//#include "mps.h"
#include <set>
#include <vector>
#include "stl_reader.h"


class cellGrid;
class mesh {
public:

	/// Initializes the mesh from the STL-file specified through path
	/// Constructor
  	//mesh(const std::string& path, const int numMesh);
  	mesh(void);
  	/// Destructor
  	virtual ~mesh(void);
  	/// Vector operations
  	static void subtraction(double *a, double *b, double *c);
	static double squaredDist(double *a, double *b);
  	/// Get variables outside polygon_wall.cpp
    long getnElements(void);
    long getnNodes(long ielem);
    double getNodePositionX(long ielem, long inode);
    double getNodePositionY(long ielem, long inode);
    double getNodePositionZ(long ielem, long inode);
	double getNormalX(long ielem);
	double getNormalY(long ielem);
	double getNormalZ(long ielem);

  	/// Read mesh file (STL)
  	void readMeshFile(const char * meshfilename);
  	/// Initialize the component values of PND and number of neighbors

	void findCandidates(int cx, int cy, int cz, int ncx, int ncy, int ncz, cellGrid *grid);
	void fillCandidates(cellGrid *grid);

	long calcDist(double *p, double reL, int cpx, int cpy, int cpz, double *closestPoint, double &cqd);
    /// Find closest point on the mesh from a particle
  	/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
    void closestPointMesh(int nMesh, int *type,  double *quaddev, double *pndS, double pndS0, double reL, int nPart, long *index, double *dMesh,
        long *idNearMesh, long *idNearElement, double *x, double *y, double *z, double *mx, double *my, double *mz, cellGrid *grid);


	///
	double calcNearestPoint(double *a, double *b, double *c, double *p, double *nearest);


private:
	struct position
	{
		double x, y, z;
	};
	/// Mesh Normals
	struct normals
	{
		position *pos;
	};
	/// Mesh vertexes
	struct nodes
	{
		position *pos;
//		long *ID;							/// specifies the node index
//		double *dx, *dy, *dz;				/// specifies the node displacement
//		double *ux, *uy, *uz;				/// specifies the node velocity
//		double *forceX, *forceY, *forceZ;	/// specifies the node force
	};
	/// Mesh elements
	struct element
	{
		nodes *node;
		long nNodes;		/// number of nodes belonging to the element
		normals normal;   /// Mesh transformation matrix
//		long ID;			/// specifies the element index
//		double press;		/// specifies the element pressure
//		long *nodeID;		/// specifies the nodes index
//		double *nodeDX, *nodeDY, *nodeDZ;	/// specifies the node displacement
//		double *nodeUX, *nodeUY, *nodeUZ;	/// specifies the node velocity

//		double* partN1;		/// specifies the shape function
//		double* partN2;
//		double* partN3;
	};
	long nElements;		                                /// number of finite elements
	element *elem;
 					/// Particles near the mesh
	std::vector<std::vector<std::vector<std::set<long> > > > candidates_set;                 
};

class cellGrid{
		public:
	/// Constructor
	cellGrid(void);

  	/// Destructor
  	virtual ~cellGrid(void);

	/// Initialization and dividing domain to get nCell && dCell
	void initCell(double dCx, double lSxMin, double lSyMin ,double lSzMin, double lSxMax, double lSyMax ,double lSzMax);

	void fillBuffer(mesh *mesh);

	void getBounds(int a[][2], double *b, double *c, double *d);
	
	///Boolean AABB - Triangle intersection Function
	///ref: https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox_tam.pdf	
	bool triangleAABB_intersection(double *cen, double normX, double normY, double normZ, double v0x, double v0y, double v0z, 
	double v1x, double v1y, double v1z, double v2x, double v2y, double v2z);

	/// Function that return the cell position
	int getCellPosX(double cx);
	int getCellPosY(double cy);
	int getCellPosZ(double cz);

	/// Gets coordinates and return the cell center

	double getCellCenterX(double cx);
	double getCellCenterY(double cy);
	double getCellCenterZ(double cz);
	
	/// Return the (d-th - 1) element id on the cell [a][b][c] if it exists, return -1 otherwise
	long getCellBuffer(int a, int b, int c, long d);
	int getNumCells(int dir);

	private:

	
	int numCells[3];
	double eCells[3];
	double minSpc[3];
	double maxSpc[3];

	///For each cell gives the triangles intersecting it
	std::vector<std::vector<std::vector<std::vector<long> > > > cellBuffer;	
	//		//		//		//

	long nEl;
};

#endif