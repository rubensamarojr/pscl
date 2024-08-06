#ifndef __POLYGON_WALL__H__
#define __POLYGON_WALL__H__

//#include "mps.h"
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <array>
#include "stl_reader.h"


class Point2 {
	public:
		Point2():x(0),y(0){}
		Point2(float a, float b):x(a),y(b){}
		float x, y;
};
class Point: public Point2 {
	public:
		Point(): Point2(), z(0) {} 
    	Point(float a, float b, float c): Point2(a, b), z(c) {}
		float z;
};

const std::array<std::array<int, 3>, 27> neighPos = {{
    {{-1, -1, -1}}, // 0
    {{-1, -1, 0}},  // 1
    {{-1, -1, 1}},  // 2
    {{-1, 0, -1}},  // 3
    {{-1, 0, 0}},   // 4
    {{-1, 0, 1}},   // 5
    {{-1, 1, -1}},  // 6
    {{-1, 1, 0}},   // 7
    {{-1, 1, 1}},   // 8
    {{0, -1, -1}},  // 9
    {{0, -1, 0}},   // 10
    {{0, -1, 1}},   // 11
    {{0, 0, -1}},   // 12
    {{0, 0, 0}},    // 13
    {{0, 0, 1}},    // 14
    {{0, 1, -1}},   // 15
    {{0, 1, 0}},    // 16
    {{0, 1, 1}},    // 17
    {{1, -1, -1}},  // 18
    {{1, -1, 0}},   // 19
    {{1, -1, 1}},   // 20
    {{1, 0, -1}},   // 21
    {{1, 0, 0}},    // 22
    {{1, 0, 1}},    // 23
    {{1, 1, -1}},   // 24
    {{1, 1, 0}},    // 25
    {{1, 1, 1}}     // 26
}};
struct hash_tuple { 
  
    template <class T1, class T2, class T3> 
  
    std::size_t operator()(const std::tuple<T1, T2, T3>& x) const
    { 
        return std::get<0>(x) 
               ^ std::get<1>(x) 
               ^ std::get<2>(x); 
    } 
}; 
  
class cellGrid{
		public:
	/// Constructor
	cellGrid(void);

  	/// Destructor
  	virtual ~cellGrid(void);

	/// Initialization and dividing domain to get nCell && dCell
	void initCell(float dCx, float lSxMin, float lSyMin ,float lSzMin, float lSxMax, float lSyMax ,float lSzMax, int nel);


	void getBounds(int a[][2], float *b, float *c, float *d);
	
	std::vector<std::array<int,3>> getNeighbors(int x, int y, int z);
	void precomputeNeighborIndices();
	int flatIndex(int x, int y, int z);
	void tdimIndex(int flatIndex, int& x, int& y, int& z) const ;
	///Boolean AABB - Triangle intersection Function
	///ref: https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox_tam.pdf	
	bool triangleAABB_intersection(float *cen, float normX, float normY, float normZ, float v0x, float v0y, float v0z, 
	float v1x, float v1y, float v1z, float v2x, float v2y, float v2z);

	/// Function that return the cell position
	int getCellPosX(float cx);
	int getCellPosY(float cy);
	int getCellPosZ(float cz);

	/// Gets coordinates and return the cell center

	float getCellCenterX(float cx);
	float getCellCenterY(float cy);
	float getCellCenterZ(float cz);

	
	/// Return the (d-th - 1) element id on the cell [a][b][c] if it exists, return -1 otherwise
	int getCellBuffer(int a, int b, int c, int d);
	int getNumCells(int dir);

	protected:
		
	int res;
	int numCells[3];
	float eCells[3];
	float minSpc[3];
	float maxSpc[3];

	bool enoughmemory;	

	///For each cell gives the triangles intersecting it
	std::vector<std::vector<int>> neighborIndices;


	std::vector<std::vector<unsigned short int> > cellBufferFlattened;
	std::vector<std::vector<std::vector<std::vector<unsigned short int> > > > cellBuffer;
	std::vector<std::unordered_set<unsigned short int>> candidatesSet_flat;
	std::vector<std::vector<std::vector<std::unordered_set<unsigned short int> > > > candidates_set;

	//		//		//		//
};
class mesh: public cellGrid {
public:

	/// Initializes the mesh from the STL-file specified through path
	/// Constructor
  	//mesh(const std::string& path, const int numMesh);
  	mesh(void);
  	/// Destructor
  	virtual ~mesh(void);
  	/// Vector operations
  	static void subtraction(const float *a, const float *b, float *c);
	static double squaredDist(const float *a, const float *b);
  	/// Get variables outside polygon_wall.cpp
    int getnElements(void);
    int getnNodes(int ielem);
    float getNodePositionX(int ielem, int inode);
    float getNodePositionY(int ielem, int inode);
    float getNodePositionZ(int ielem, int inode);
	float getNormalX(int ielem);
	float getNormalY(int ielem);
	float getNormalZ(int ielem);
	int getCellPosX(float cx);
	int getCellPosY(float cy);
	int getCellPosZ(float cz);
	void set_partdist(float a){partdist2 = double(a)*double(a);};

	void set_reserve(bool flat);
	bool do_reserve(){return enoughmemory;};
  	/// Read mesh file (STL)
  	void readMeshFile(const char * meshfilename,const float dcell,const float  VminX,const float  VminY,const float  VminZ,const float  VmaxX,const float  VmaxY,const float  VmaxZ, int opt, int simR);
  	/// Initialize the component values of PND and number of neighbors

	static Point2 intersectionPointZ(const Point& p1, const Point& p2, const float p);
	bool intersectTrianglePlaneZ(const float z, std::vector<Point2> &intersectionPoints, const int t);
	void planeInterpolationXY(const Point2 &a, const Point2 &b, const float e,std::vector<std::set<int> > &cellcoords);
	void planeInterpolationXZ(const Point2 &a, const Point2 &b, const float e,std::vector<std::set<int> > &cellcoords);
	void planeInterpolationYZ(const Point2 &a, const Point2 &b, const float e,std::vector<std::set<int> > &cellcoords);
	void clearCellBuffer();
	void clearCellBuffer_flat();
	void fillBuffer();

	// void fillbx(const std::vector<std::vector<int> >&cellcd, const int x, const int t);
	// void fillby(const std::vector<std::vector<int> >&cellcd, const int y, const int t);
	// void fillb(const std::vector<std::vector<int> >&cellcd, const int z, const int t);
	void fillBopt();
	void fillBopt_flat();
	void fillBopt_reserve();
	void fillBopt_reserveflat();
	// void findCandidates(int cx, int cy, int cz, int ncx, int ncy, int ncz);
	void clearCandidates();
	void clearCandidates_flat();

	void fillCandidates();
	void fillCandidates2();
	void fillCandidates3();
	void fillCandidates4();
	void fillCandidates5();
	void fillCandidates_seq();
	void fillCandidates_flat();

	int calcDist(const float *p, const double reL2, const int cpx, const int cpy, const int cpz, float *closestPoint, double &cqd);
	int calcDist_flat(const float *p, const double reL2, const int fIndex, float *closestPoint, double &cqd);
    /// Find closest point on the mesh from a particle
  	/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
    void closestPointMesh(int nMesh, int *type,  float *quaddev, float *pndS, float pndS0, float reL, int nPart, int *index, float *dMesh,
        int *idNearMesh, int *idNearElement, float *x, float *y, float *z, float *mx, float *my, float *mz);


	///
	double calcNearestPoint(const float *a, const float *b, const float *c, const float *p, float *nearest);
	void writeCellBufferToFile(const std::string& filename);
private:
	struct position
	{
		float x, y, z;
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
//		int *ID;							/// specifies the node index
//		float *dx, *dy, *dz;				/// specifies the node displacement
//		float *ux, *uy, *uz;				/// specifies the node velocity
//		float *forceX, *forceY, *forceZ;	/// specifies the node force
	};
	/// Mesh elements
	struct element
	{
		nodes *node;
		int nNodes;		/// number of nodes beinting to the element
		normals normal;   /// Mesh transformation matrix
//		int ID;			/// specifies the element index
//		float press;		/// specifies the element pressure
//		int *nodeID;		/// specifies the nodes index
//		float *nodeDX, *nodeDY, *nodeDZ;	/// specifies the node displacement
//		float *nodeUX, *nodeUY, *nodeUZ;	/// specifies the node velocity

//		float* partN1;		/// specifies the shape function
//		float* partN2;
//		float* partN3;
	};
	int nElements;		                                /// number of finite elements
	int simResolution;
	element *elem;
 					/// Particles near the mesh
	double partdist2;
};


#endif