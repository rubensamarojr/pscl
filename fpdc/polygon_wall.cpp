 /** MPS - Moving Particle Semi-
 ***(c) 2019 USP/TPN
 *** Rubens Amaro, Lucas Pereira
 **/

#include "polygon_wall.h"
//#include "mps.h"
#include <immintrin.h> // For SIMD intrinsics
#include <algorithm>   // For std::min, std::max
#include <cfloat>      // For FLT_MAX
#include <iostream>
#include <cmath>
#include <omp.h>
#define eps static_cast<float>(10e-6)

/*
mesh::mesh(const std::string& path, const int numMesh) {

	//solid = new solid_fem[numSolids];
	solidWall = new mesh[numMesh];

}
*/

///
/// CELL GRID CLASS METHODS
///
cellGrid::cellGrid()
{
}

cellGrid::~cellGrid(void)
{
}




void cellGrid::initCell(float dCx, float lSxMin, float lSyMin ,float lSzMin, float lSxMax, float lSyMax ,float lSzMax, int nel){

	///Need to be called from init.cpp with the mesh
    ///Memory allocation


	///Assigning values to the cell class atributes	
	///// Same as the parameters from MPS particles
	
    eCells[0] = dCx;		//dCell[0]
    eCells[1] = dCx;		//dCell[1]
    eCells[2] = dCx;		//dCell[2]

	//fprintf(stderr, "\n* QUANTIDADE DE CELULAS [x][y][z]:::[%d][%d][%d] ", nCx,nCy,nCz);
	fprintf(stderr, "\n* CELLSIZE [x][y][z]:::[%f][%f][%f] ", dCx,dCx,dCx);

	int nCx = static_cast<int>((lSxMax - lSxMin)/dCx) +1;
	int nCy = static_cast<int>((lSyMax - lSyMin)/dCx) +1;
	int nCz = static_cast<int>((lSzMax - lSzMin)/dCx) +1;


    numCells[0] = nCx; 	//nCell[0]
    numCells[1] = nCy; 	//nCell[1]
    numCells[2] = nCz;	//nCell[2]
	printf("\nnumcells = %d,%d,%d\n",numCells[0],numCells[1],numCells[2]);

	/// Simulation boundaries(max and min coordinates on each direction)

	minSpc[0] = lSxMin;				//limSpace[0][1]
	minSpc[1] = lSyMin;				//limSpace[1][1]
	minSpc[2] = lSzMin;				//limSpace[2][1]
	

	maxSpc[0] = lSxMax;				//limSpace[0][1]
	maxSpc[1] = lSyMax;				//limSpace[1][1]
	maxSpc[2] = lSzMax;				//limSpace[2][1]

	// double multi_factor = 1;
	// if (nCx*nCy*nCz > 10000){
	// 	multi_factor = 10;
	// 	if (double(std::min(nCx,std::min(nCy,nCz)))/double(std::max(nCx,std::max(nCy,nCz))) < 0.1) multi_factor = 1;
	// }
	// if (nCx*nCy*nCz > 100000){
	// 	multi_factor = 5;
	// 	if (double(std::min(nCx,std::min(nCy,nCz)))/double(std::max(nCx,std::max(nCy,nCz))) < 0.1) multi_factor = 1;
	// }
	// if(nCx*nCy*nCz > 1000000){
	// 	multi_factor = 2;
	// 	if (double(std::min(nCx,std::min(nCy,nCz)))/double(std::max(nCx,std::max(nCy,nCz))) < 0.1) multi_factor = 1;
	// }
	// if(nCx*nCy*nCz > 5000000){
	// 	multi_factor = 0.5;
	// 	if (double(std::min(nCx,std::min(nCy,nCz)))/double(std::max(nCx,std::max(nCy,nCz))) < 0.1) multi_factor = 1;
	// }
	// if(nCx*nCy*nCz > 10000000){
	// 	multi_factor = 0.1;
	// 	if (double(std::min(nCx,std::min(nCy,nCz)))/double(std::max(nCx,std::max(nCy,nCz))) < 0.1) multi_factor = 1;
	// }
	// printf("numElements: %d\n",nel);
	// printf("totalcells = %d\n", nCx*nCy*nCz);
	// multi_factor = 0.25;
	// if(nCx*nCy*nCz < 10000) res = nel;
	// if(res > nel) res = nel;

	int totalCells = numCells[0] * numCells[1] * numCells[2];
	
	candidatesSet_flat.resize(totalCells);
	cellBufferFlattened.resize(totalCells);
	precomputeNeighborIndices();
	

	cellBuffer.resize(numCells[0]);
	for(int i=0; i<numCells[0]; i++)
	{
		cellBuffer[i].resize(numCells[1]);
		for(int j=0; j<numCells[1]; j++){
			cellBuffer[i][j].resize(numCells[2]);
			// 	for(int k=0; k<numCells[2]; k++){
	// 		// 		candidates_set[i][j][k].reserve(res);
				// }
		}
	
	}

}

void cellGrid::getBounds(int a[][2], float *b, float *c, float *d){
	// float e = 0.55*eCells[0];
	float e = 0.0*eCells[0];
	/// Minimum
	///(X_Vertex0, X_Vertex1, X_Vertex2)
	int x0 = getCellPosX(b[0]-e);
	int x1 = getCellPosX(b[1]-e);
	int x2 = getCellPosX(b[2]-e);

	int aux = std::min(x0, x1);

	a[0][0] = std::min(aux, x2);

	///(Y_Vertex1, Y_Vertex2, Y_Vertex3)
	int y0 = getCellPosY(c[0]-e);
	int y1 = getCellPosY(c[1]-e);
	int y2 = getCellPosY(c[2]-e);

	aux = std::min(y0, y1);
	a[1][0] = std::min(aux, y2);
	
	///(Z_Vertex1, Z_Vertex2, Z_Vertex3)
	int z0 = getCellPosZ(d[0]-e);
	int z1 = getCellPosZ(d[1]-e);
	int z2 = getCellPosZ(d[2]-e);
	
	aux = static_cast<int>(std::min(z0, z1));
	a[2][0] = std::min(aux, z2);
	
	/// Maximum
	///(X_Vertex1, X_Vertex2, X_Vertex3)
	x0 = getCellPosX(b[0]+e);
	x1 = getCellPosX(b[1]+e);
	x2 = getCellPosX(b[2]+e);
	a[0][1] = std::max(std::max(x0, x1), x2);

	///(Y_Vertex1, Y_Vertex2, Y_Vertex3)
	y0 = getCellPosY(c[0]+e);
	y1 = getCellPosY(c[1]+e);
	y2 = getCellPosY(c[2]+e);
	a[1][1] = std::max(std::max(y0, y1), y2);

	///(Z_Vertex1, Z_Vertex2, Z_Vertex3)
	z0 = getCellPosZ(d[0]+e);
	z1 = getCellPosZ(d[1]+e);
	z2 = getCellPosZ(d[2]+e);
	a[2][1] = std::max(std::max(z0, z1), z2);
	/*
	fprintf(stderr, "\n* EXTREMOS CELULAS A SEREM TESTADAS");
	fprintf(stderr, "\n* MIN[x][y][z]:::[%d][%d][%d] ", a[0][0],a[1][0],a[2][0]);
	fprintf(stderr, "\n* MAX[x][y][z]:::[%d][%d][%d] ", a[0][1],a[1][1],a[2][1]);
	*/
}


bool cellGrid::triangleAABB_intersection(float *cen, float normX, float normY, float normZ, float v0x, float v0y, float v0z, 
float v1x, float v1y, float v1z, float v2x, float v2y, float v2z){
	/*
	/ ref: https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox_tam.pdf
	/ cen[2] = center of aabb coordinates
	/ v[3][3] = vertices of mesh element(triangle)
	/ cell->ecell->spc->i = extent of aabb on i direction; i=0,1,2
	*/

	int X,Y,Z;
	float v0[3] = {0.0,0.0,0.0};
	float v1[3] = {0.0,0.0,0.0};
	float v2[3] = {0.0,0.0,0.0};
	
	float normalAux[3] = {normX,normY,normZ};
	float v[3][3]={v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z};

	X = 0;
	Y = 1;
	Z = 2;
	
	
	//Moving AABB and vertices to origin
	mesh::subtraction(v[0],cen,v0);
	mesh::subtraction(v[1],cen,v1);
	mesh::subtraction(v[2],cen,v2); /*
	if (f) fprintf(stderr, "\n cen = [%f][%f][%f]*", cen[X], cen[Y], cen[Z]);
	if (f) fprintf(stderr, "\n v0x = [%f][%f][%f]*", v0x, v0y, v0z);
	if (f) fprintf(stderr, "\n v0[X] = [%f][%f][%f]*", v0[X], v0[Y], v0[Z]);
	if (f) fprintf(stderr, "\n 1*");
	*/

	if (std::max(std::max(v0[X], v1[X]), v2[X]) < -(eCells[0])||std::min(std::min(v0[X], v1[X]), v2[X]) > eCells[0])
		return false;
		
	//if (f) fprintf(stderr, "\n 2*");
	if (std::max(std::max(v0[Y], v1[Y]), v2[Y]) < -(eCells[1])||std::min(std::min(v0[Y], v1[Y]), v2[Y]) > eCells[1])
		return false;
		
	//if (f) fprintf(stderr, "\n 3*");
	if (std::max(std::max(v0[Z], v1[Z]), v2[Z]) < -(eCells[2])||std::min(std::min(v0[Z], v1[Z]), v2[Z]) > eCells[2])
		return false;

	
	//9 Axis tests considering the cross product between each aabb normal and triangle edge
	//a_ij = ei ^ fj ; i,j = 0,1,2

	//Edge vectors
	// f0 = v1 - v0
	// f1 = v2 - v1
	// f2 = v0 - v2

	//Projection array(projecting triangle vertices onto aij vectors)
	//p[0]-> value that repeats among the three projections; p[1]-> unique value among the three projections
	float p[2], r;

	//
	//a0:
	//

	//a00:
	//a[0][0] = [0,-f0[Z],f0[Y]] = [0,v0[Z]- v1[Z], v1[Y] - v0[Y]] -> array A unnecessary

	// p = dot product <a,v>; 
	// p0 = p1 = p[0]
	// p2 = p[1]
	p[0] = v0[Z]*v1[Y] - v0[Y]*v1[Z];
	p[1] = v2[Y]*(v0[Z]-v1[Z]) + v2[Z]*(v1[Y]-v0[Y]);

	r = eCells[1] * fabs(v1[Z] - v0[Z]) + eCells[2] * fabs(v1[Y] - v0[Y]);
	
	//if (f) fprintf(stderr, "\n 4*");
	if (std::max(p[0],p[1])<-r||std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	//a01:
	//a[0][1] = [0,-f1[Z],f1[Y]] = [0,v1[Z]- v2[Z], v2[Y] - v1[Y]] ; 

	//P_i = dot product <a,v_i>; 
	//p1 = p2 = p[0]
	//p0 = p[1]
	p[0] = v1[Z]*v2[Y] - v1[Y]*v2[Z];
	p[1] = v0[Y]*(v1[Z]-v2[Z]) + v0[Z]*(v2[Y]-v1[Y]);

	r = eCells[1] * fabs(v2[Z] - v1[Z]) + eCells[2] * fabs(v2[Y] - v1[Y]);
	
	//if (f) fprintf(stderr, "\n 5*");
	if (std::max(p[0],p[1])<-r||std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	//a02:
	//a[0][2] = [0,-f2[Z],f2[Y]] = [0,v2[Z]- v0[Z], v0[Y] - v2[Y]] 

	//P_i = dot product <a,v_i>; 
	//p0 = p2 = p[0]
	//p1 = p[1]
	p[0] = v2[Z]*v0[Y] - v2[Y]*v0[Z];
	p[1] = v1[Y]*(v2[Z]-v0[Z]) + v1[Z]*(v0[Y]-v2[Y]);

	r = eCells[1] * fabs(v0[Z] - v2[Z]) + eCells[2] * fabs(v0[Y] - v2[Y]);
	//if (f) fprintf(stderr, "\n 6*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	
	//
	//a1:
	//

	//a10:
	//a[1][0] = [f0[Z],0,-f0[X]] = [v1[Z]- v0[Z], 0, v0[X] - v1[X]] 
	// p = dot product <a,v>; 
	// p0 = p1 = p[0]
	// p2 = p[1]
	p[0] = v0[X]*v1[Z] - v0[Z]*v1[X];
	p[1] = v2[X]*(v1[Z]-v0[Z]) + v2[Z]*(v0[X]-v1[X]);

	r = eCells[0] * fabs(v1[Z] - v0[Z]) + eCells[2] * fabs(v1[X] - v0[X]);
	//if (f) fprintf(stderr, "\n 7*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	//a11:
	//a[1][1] = [f1[Z],0,-f1[X]] = [v2[Z]- v1[Z], 0, v1[X] - v2[X]] 

	//p = dot product <a,v>; 
	//p1 = p2 = p[0]
	//p0 = p[1]

	p[0] = v1[X]*v2[Z] - v1[Z]*v2[X];
	p[1] = v0[X]*(v2[Z]-v1[Z]) + v0[Z]*(v1[X]-v2[X]);

	r = eCells[0] * fabs(v2[Z] - v1[Z]) + eCells[2] * fabs(v2[X] - v1[X]);
	//if (f) fprintf(stderr, "\n 8*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	//a12:
	//a[1][2] = [f2[Z],0,-f2[X]] = [v0[Z]- v2[Z], 0, v2[X] - v0[X]] 

	//p = dot product <a,v>; 
	//p0 = p2 = p[0]
	//p1 = p[1]

	p[0] = v2[X]*v0[Z] - v2[Z]*v0[X];
	p[1] = v1[X]*(v0[Z]-v2[Z]) + v1[Z]*(v2[X]-v0[X]);

	r = eCells[0] * fabs(v0[Z] - v2[Z]) + eCells[2] * fabs(v0[X] - v2[X]);
	//if (f) fprintf(stderr, "\n 9*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;	
		
		
	//
	//a2:
	//

	//a20:
	//a[2][0] = [-f0[Y],f0[X],0] = [v0[Y]- v1[Y], v1[X] - v0[X],0] 
	// p = dot product <a,v>; 
	// p0 = p1 = p[0]
	// p2 = p[1]
	p[0] = v0[Y]*v1[X] - v0[X]*v1[Y];
	p[1] = v2[X]*(v0[Y]-v1[Y]) + v2[Y]*(v1[X]-v0[X]);

	r = eCells[0] * fabs(v1[Y] - v0[Y]) + eCells[1] * fabs(v1[X] - v0[X]);
	//if (f) fprintf(stderr, "\n 10*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	//a21:
	//a[2][1] = [-f1[Y],f1[X],0] = [v1[Y]- v2[Y], v2[X] - v1[X],0] 

	// p = dot product <a,v>; 
	// p1 = p2 = p[0]
	// p0 = p[1]
	p[0] = v1[Y]*v2[X] - v1[X]*v2[Y];
	p[1] = v0[X]*(v1[Y]-v2[Y]) + v0[Y]*(v2[X]-v1[X]);

	r = eCells[0] * fabs(v2[Y] - v1[Y]) + eCells[1] * fabs(v2[X] - v1[X]);
	//if (f) fprintf(stderr, "\n 11*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;


	//a22:
	//a[2][0] = [-f2[Y],f2[X],0] = [v2[Y]- v0[Y], v0[X] - v2[X],0] 

	// p = dot product <a,v>; 
	// p0 = p2 = p[0]
	// p1 = p[1]
	p[0] = v2[Y]*v0[X] - v2[X]*v0[Y];
	p[1] = v1[X]*(v2[Y]-v0[Y]) + v1[Y]*(v0[X]-v2[X]);

	r = eCells[0] * fabs(v0[Y] - v2[Y]) + eCells[1] * fabs(v0[X] - v2[X]);
	//if (f) fprintf(stderr, "\n 12*");
	if (std::max(p[0],p[1])<-r || std::min(p[0],p[1])>r) 
	// if the greater projection cannot be greater than negative radius or the lower projection cannot be lower than the radius => separating axis found
		return false;

	// Testing separating axis from the triangle normal
	//Note that calculating the plane distance on this case requires v0 before translation to origin
	const float plane_dist = normalAux[X]*v0x + normalAux[Y]*v0y + normalAux[Z]*v0z;

	const float s = normalAux[X]*cen[X] + normalAux[Y]*cen[Y] + normalAux[Z]*cen[Z] - plane_dist;
	r = eCells[0] * fabs(normalAux[X]) + eCells[1] * fabs(normalAux[Y]) + eCells[2] * fabs(normalAux[Z]);
	return (fabs(s)<=r);
}


/// Calculate the Cell position on direction by the given direction coordinate

int cellGrid::getCellPosX(float cx){
	
	//fprintf(stderr, "\n*cell coordinate x - %f", cx);
	//fprintf(stderr, "\n*cell position x - %d", static_cast<int>((cx - minSpc[0])/eCells[0]));
	return static_cast<int>((cx - minSpc[0])/eCells[0]);
}
int cellGrid::getCellPosY(float cy){
	
	//fprintf(stderr, "\n*cell coordinate y - %f", cy);
	//fprintf(stderr, "\n*cell position y - %d", static_cast<int>((cy - minSpc[1])/eCells[1]));
	return static_cast<int>((cy - minSpc[1])/eCells[1]);
}
int cellGrid::getCellPosZ(float cz){
	
	return static_cast<int>((cz - minSpc[2])/eCells[2]);
}


/// Calculate the Cell center coordinate on direction by the given direction position
float cellGrid::getCellCenterX(float cx){

	return (cx+0.5)*eCells[0] + minSpc[0];
}
float cellGrid::getCellCenterY(float cy){

	return (cy+0.5)*eCells[1] + minSpc[1];
}
float cellGrid::getCellCenterZ(float cz){

	return (cz+0.5)*eCells[2] + minSpc[2];
}



int cellGrid::getNumCells(int dir){
	return numCells[dir];
}

/// Return the (d-th - 1) element id on the cell [a][b][c] if it exists, return -1 otherwise
int cellGrid::getCellBuffer(int a, int b, int c, int d){
	if (cellBuffer[a][b][c].size() <= d) return -1; ///Returning -1 means the cell does not have any more elements
	return cellBuffer[a][b][c][d];

}

///
///MESH CLASS METHODS
///
mesh::mesh(void)
{
}

mesh::~mesh(void)
{
}


/// Function to find
/// addition of two vector array.

/// Function to find
/// subtraction of two vector array.
void mesh::subtraction(const float *a, const float *b, float *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

double mesh::squaredDist(const float *a,const float *b)
{
	return (a[0]-b[0]) * (a[0]-b[0]) + (a[1]-b[1]) * (a[1]-b[1]) + (a[2]-b[2]) * (a[2]-b[2]);
}

int mesh::getnElements(void)
{
   return nElements;
}

int mesh::getnNodes(int ielem)
{
   return 3;
}

float mesh::getNodePositionX(int ielem, int inode)
{
	if (inode == 0){
		return meshData.x0[ielem];
	}
	if (inode == 1){
		return meshData.x1[ielem];
	}
	if (inode == 2){
		return meshData.x2[ielem];
	}
}

float mesh::getNodePositionY(int ielem, int inode)
{
	if (inode == 0){
		return meshData.y0[ielem];
	}
	if (inode == 1){
		return meshData.y1[ielem];
	}
	if (inode == 2){
		return meshData.y2[ielem];
	}
}

float mesh::getNodePositionZ(int ielem, int inode)
{
	if (inode == 0){
		return meshData.z0[ielem];
	}
	if (inode == 1){
		return meshData.z1[ielem];
	}
	if (inode == 2){
		return meshData.z2[ielem];
	}
}

float mesh::getNormalX(int ielem)
{
   return elem[ielem].normal.pos->x;
}
float mesh::getNormalY(int ielem)
{
   return elem[ielem].normal.pos->y;
}
float mesh::getNormalZ(int ielem)
{
   return elem[ielem].normal.pos->z;
}

/// Calculate the Cell position on direction by the given direction coordinate

int mesh::getCellPosX(float cx){
	
	//fprintf(stderr, "\n*cell coordinate x - %f", cx);
	//fprintf(stderr, "\n*cell position x - %d", static_cast<int>((cx - minSpc[0])/eCells[0]));
	return static_cast<int>((cx - minSpc[0])/eCells[0]);
}
int mesh::getCellPosY(float cy){
	
	//fprintf(stderr, "\n*cell coordinate y - %f", cy);
	//fprintf(stderr, "\n*cell position y - %d", static_cast<int>((cy - minSpc[1])/eCells[1]));
	return static_cast<int>((cy - minSpc[1])/eCells[1]);
}
int mesh::getCellPosZ(float cz){
	
	return static_cast<int>((cz - minSpc[2])/eCells[2]);
}

void mesh::set_reserve(bool flat){
	enoughmemory = false;
	res = static_cast<int>(nElements/std::min(numCells[0],std::min(numCells[1],numCells[2])));
	
	printf("Estimate Reserve for each cell = %d\n",res);
	printf("Estimate usage of memory: %d Mb\n", int(res*numCells[0]*numCells[1]*numCells[2]*10e-6*2));
	if ((res*numCells[0]*numCells[1]*numCells[2]*10e-6*2) < 1500){
		enoughmemory = true;
		if (flat){
			for(int i=0; i<numCells[0]*numCells[1]*numCells[2]; i++) {cellBufferFlattened[i].reserve(res);}//candidates_set[i][j][k].reserve(res); }
		
		}else{
			for(int i=0; i<numCells[0]; i++) for(int j=0; j<numCells[1]; j++) for(int k=0; k<numCells[2]; k++) {cellBuffer[i][j][k].reserve(res);}//candidates_set[i][j][k].reserve(res); }
		}
	}
}
/// https://github.com/sreiter/stl_reader
void mesh::readMeshFile(const char * meshfilename,const float dcell,const float  VminX,const float  VminY,const float  VminZ,const float  VmaxX,const float  VmaxY,const float  VmaxZ, int opt, int simR)
{	
	std::vector<float> coords, normal_vecs;
	std::vector<unsigned int> elems, solids;
	try {

		// Load a mesh
		bool success = stl_reader::ReadStlFile(meshfilename, coords, normal_vecs, elems, solids);
		
		nElements = elems.size() / 3;
		// printf("\nnElements = %d\n", nElements );
		elem = (struct element*) malloc(nElements*sizeof(struct element));
		
		
		meshData.x0.resize(nElements);
		meshData.y0.resize(nElements);
		meshData.z0.resize(nElements);
		meshData.x1.resize(nElements);
		meshData.y1.resize(nElements);
		meshData.z1.resize(nElements);
		meshData.x2.resize(nElements);
		meshData.y2.resize(nElements);
		meshData.z2.resize(nElements);

		for(size_t ielem = 0; ielem < nElements; ++ielem) {
			
			float* c0 = &coords[3 * elems[3 * ielem]];
			meshData.x0[ielem] = c0[0];
			meshData.y0[ielem] = c0[1];
			meshData.z0[ielem] = c0[2];
	
			float* c1 = &coords[3 * elems[3 * ielem + 1]];
			meshData.x1[ielem] = c1[0];
			meshData.y1[ielem] = c1[1];
			meshData.z1[ielem] = c1[2];

			float* c2 = &coords[3 * elems[3 * ielem + 2]];
			meshData.x2[ielem] = c2[0];
			meshData.y2[ielem] = c2[1];
			meshData.z2[ielem] = c2[2];	
			
			//std::cout << std::endl;
		
			float* n = &normal_vecs[3 * ielem];
			//std::cout   << "normal of triangle " << ielem << ": "
			//			<< "(" << n[0] << ", " << n[1] << ", " << n[2] << ")\n";
						
			elem[ielem].normal.pos = (struct position *) malloc(1*sizeof(struct position));
			elem[ielem].normal.pos->x = n[0];
			elem[ielem].normal.pos->y = n[1];
			elem[ielem].normal.pos->z = n[2];
		}
	}
	catch (std::exception& e) {
	std::cout << e.what() << std::endl;
	}

	simResolution = simR;
	initCell(dcell, VminX, VminY, VminZ, VmaxX, VmaxY, VmaxZ, nElements);
	
    set_reserve(opt == 2);		
	candidates_set.resize(numCells[0]);
		for(int i=0; i<numCells[0]; i++)
		{
			candidates_set[i].resize(numCells[1]);
			for(int j=0; j<numCells[1]; j++){
				candidates_set[i][j].resize(numCells[2]);
				// for (int k = 0; k < numCells[2];k++){
				// 	candidates_set[i][j][k].reserve(3*res);
				// }
			}

		}

}

// CALCULA INTERSECÇÃO - 0.5 -> RETORNA PONTOS 
// USADO EM 1
Point2 mesh::intersectionPointZ(const Point& p1, const Point& p2, const float p) {
    Point edgeVector = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    float t = -((p1.z - p) / edgeVector.z);
    return {p1.x + t * edgeVector.x, p1.y + t * edgeVector.y};
}

// CORTAPLANO - 1 -> RETORNA LINHA
bool mesh::intersectTrianglePlaneZ(const float z, std::vector<Point2> &intersectionPoints, const int t) {
    int c = 0;

    // Define triangle vertices
    Point p0(meshData.x0[t], meshData.y0[t], meshData.z0[t]);
    Point p1(meshData.x1[t], meshData.y1[t], meshData.z1[t]);
    Point p2(meshData.x2[t], meshData.y2[t], meshData.z2[t]);

    Point triangle[3] = {p0, p1, p2};

    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        Point p1 = triangle[i];
        Point p2 = triangle[j];

        if (fabs(p1.z - p2.z) < eps && fabs(p1.z - z) < eps) {
            intersectionPoints[0] = Point2(p1.x, p1.y);
            intersectionPoints[1] = Point2(p2.x, p2.y);
            return true;
        }

        if ((p1.z - z) * (p2.z - z) <= 0) {
            intersectionPoints[c] = intersectionPointZ(p1, p2, z);
            c++;
        }
    }

    if (c == 1) intersectionPoints[1] = intersectionPoints[0];
    return c != 0;
}
// LINHA PRA CELULA - 2 - RETORNA CELULA
void mesh::planeInterpolationXY(const Point2 &a, const Point2 &b, const float e, std::vector<std::set<int>>& cellcoords) {
    const float einv = 1/e;
    // Ajustando as coordenadas dos pontos a e b para considerar o deslocamento do ponto de origem do domínio
    const int p0x = static_cast<int>((a.x - minSpc[0]) *einv);
    const int p1x = static_cast<int>((b.x - minSpc[0]) *einv);
    const int p0y = static_cast<int>((a.y - minSpc[1]) *einv);
    const int p1y = static_cast<int>((b.y - minSpc[1]) *einv);

    // Adicionando as coordenadas do ponto inicial ao vetor de coordenadas da célula
    cellcoords[p0x].insert(p0y);

    // Calculando os deltas de interpolação
    const float dx = (b.x - a.x) / (b.y - a.y);
    const float dy = (b.y - a.y) / (b.x - a.x);
    
    if (p1x > p0x) {
        for (int i = p0x + 1; i <= p1x; i++) {
            // float inter = a.y - minSpc[1] + (i * e + minSpc[0] - a.x ) * dy;
            cellcoords[i].insert(static_cast<int>((a.y - minSpc[1] + (i * e + minSpc[0] - a.x ) * dy) *einv));
            // cellcoords[i].insert(static_cast<int>((a.y - minSpc[1] + (i * e + minSpc[0] - a.x +eps) * dy) *einv));
        }
    } else {
        for (int i = p0x; i > p1x; i--) {
            // float inter = a.y - minSpc[1] + (i * e - a.x + minSpc[0]) * dy;
            cellcoords[i-1].insert(static_cast<int>((a.y - minSpc[1] + (i * e - a.x + minSpc[0]) * dy) *einv));
            // cellcoords[i-1].insert(static_cast<int>((a.y - minSpc[1] + (i * e - a.x + minSpc[0]-eps) * dy) *einv));
        }
    }

    if (p1y > p0y) {
        for (int j = p0y + 1; j <= p1y; j++) {
            // float inter = a.x - minSpc[0] + (j * e - a.y + minSpc[1]) * dx;
            cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[1]) * dx) *einv)].insert(j);
            // cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[1]+eps) * dx) *einv)].insert(j);
        }
    } else {
        for (int j = p0y; j > p1y; j--) {
            // float inter = a.x - minSpc[0] + (j * e - a.y + minSpc[1]) * dx;
            cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[1]) * dx) *einv)].insert(j-1);
            // cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[1]-eps) * dx) *einv)].insert(j-1);
        }
    }

}

// LINHA PRA CELULA - 2 - RETORNA CELULA
void mesh::planeInterpolationXZ(const Point2 &a, const Point2 &b, const float e, std::vector<std::set<int>>& cellcoords) {
    const float einv = 1/e;
    
    // Ajustando as coordenadas dos pontos a e b para considerar o deslocamento do ponto de origem do domínio
    const int p0x = static_cast<int>((a.x - minSpc[0]) *einv);
    const int p1x = static_cast<int>((b.x - minSpc[0]) *einv);
    const int p0z = static_cast<int>((a.y - minSpc[2]) *einv);
    const int p1z = static_cast<int>((b.y - minSpc[2]) *einv);

    // Adicionando as coordenadas do ponto inicial ao vetor de coordenadas da célula
    cellcoords[p0x].insert(p0z);

    // Calculando os deltas de interpolação
    const float dx = (b.x - a.x) / (b.y - a.y);
    const float dz = (b.y - a.y) / (b.x - a.x);
    
    if (p1x > p0x) {
        for (int i = p0x + 1; i <= p1x; i++) {
            // float inter = a.y - minSpc[2] + (i * e + minSpc[0] - a.x ) * dz;
            cellcoords[i].insert(static_cast<int>((a.y - minSpc[2] + (i * e + minSpc[0] - a.x ) * dz) *einv));
            // cellcoords[i].insert(static_cast<int>((a.y - minSpc[2] + (i * e + minSpc[0] - a.x +eps) * dz) *einv));
        }
    } else {
        for (int i = p0x; i > p1x; i--) {
            // float inter = a.y - minSpc[2] + (i * e - a.x + minSpc[0]) * dz;
            cellcoords[i-1].insert(static_cast<int>((a.y - minSpc[2] + (i * e - a.x + minSpc[0]) * dz) *einv));
            // cellcoords[i-1].insert(static_cast<int>((a.y - minSpc[2] + (i * e - a.x + minSpc[0]-eps) * dz) *einv));
        }
    }

    if (p1z > p0z) {
        for (int j = p0z + 1; j <= p1z; j++) {
            // float inter = a.x - minSpc[0] + (j * e - a.y + minSpc[2]) * dx;
            cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[2]) * dx) *einv)].insert(j);
            // cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[2]) * dx+eps) *einv)].insert(j);
        }
    } else {
        for (int j = p0z; j > p1z; j--) {
            // float inter = a.x - minSpc[0] + (j * e - a.y + minSpc[2]) * dx;
            cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[2]) * dx) *einv)].insert(j-1);
            // cellcoords[static_cast<int>((a.x - minSpc[0] + (j * e - a.y + minSpc[2]) * dx-eps) *einv)].insert(j-1);
        }
    }

}

// LINHA PRA CELULA - 2 - RETORNA CELULA
void mesh::planeInterpolationYZ(const Point2 &a, const Point2 &b, const float e, std::vector<std::set<int>>& cellcoords) {
    const float einv = 1/e;
    
    // Ajustando as coordenadas dos pontos a e b para considerar o deslocamento do ponto de origem do domínio
    const int p0y = static_cast<int>((a.x - minSpc[1]) *einv);
    const int p0z = static_cast<int>((a.y - minSpc[2]) *einv);
    const int p1y = static_cast<int>((b.x - minSpc[1]) *einv);
    const int p1z = static_cast<int>((b.y - minSpc[2]) *einv);

    // Adicionando as coordenadas do ponto inicial ao vetor de coordenadas da célula
    cellcoords[p0y].insert(p0z);

    // Calculando os deltas de interpolação
    const float dy = (b.x - a.x) / (b.y - a.y);
    const float dz = (b.y - a.y) / (b.x - a.x);
    
    if (p1y > p0y) {
        for (int i = p0y + 1; i <= p1y; i++) {
            // float inter = a.y - minSpc[2] + (i * e + minSpc[1] - a.x ) * dz;
            cellcoords[i].insert(static_cast<int>((a.y - minSpc[2] + (i * e + minSpc[1] - a.x ) * dz) *einv));
            // cellcoords[i].insert(static_cast<int>((a.y - minSpc[2] + (i * e + minSpc[1] - a.x +eps) * dz) *einv));
        }
    } else {
        for (int i = p0y; i > p1y; i--) {
            // float inter = a.y - minSpc[2] + (i * e - a.x + minSpc[1]) * dz;
            cellcoords[i-1].insert(static_cast<int>((a.y - minSpc[2] + (i * e - a.x + minSpc[1]) * dz) *einv));
            // cellcoords[i-1].insert(static_cast<int>((a.y - minSpc[2] + (i * e - a.x + minSpc[1]-eps) * dz) *einv));
        }
    }

    if (p1z > p0z) {
        for (int j = p0z + 1; j <= p1z; j++) {
            // float inter = a.x - minSpc[1] + (j * e - a.y + minSpc[2]) * dy;
            cellcoords[static_cast<int>((a.x - minSpc[1] + (j * e - a.y + minSpc[2]) * dy) *einv)].insert(j);
            // cellcoords[static_cast<int>((a.x - minSpc[1] + (j * e - a.y + minSpc[2]+eps) * dy) *einv)].insert(j);
        }
    } else {
        for (int j = p0z; j > p1z; j--) {
            // float inter = a.x - minSpc[1] + (j * e - a.y + minSpc[2]) * dy;
            cellcoords[static_cast<int>((a.x - minSpc[1] + (j * e - a.y + minSpc[2]) * dy) *einv)].insert(j-1);
            // cellcoords[static_cast<int>((a.x - minSpc[1] + (j * e - a.y + minSpc[2]-eps) * dy) *einv)].insert(j-1);
        }
    }

}

void mesh::clearCellBuffer(){
	int  i, j, k;
	for(i=0; i<(numCells[0]); i++)
		for(j=0; j<(numCells[1]); j++)
			for(k=0;k<numCells[2]; k++)
				cellBuffer[i][j][k].clear();
}
void mesh::clearCellBuffer_flat(){
	for(int i=0; i<numCells[0]*numCells[1]*numCells[2]; i++)
				cellBufferFlattened[i].clear();
}
void mesh::clearCandidates_flat(){
	for(int i=0; i<numCells[0]*numCells[1]*numCells[2]; i++)
				candidatesSet_flat[i].clear();
}

void mesh::clearCandidates() {
    for (auto& plane : candidates_set) {
        for (auto& row : plane) {
            for (auto& cell : row) {
                cell.clear();
            }
        }
    }
}
int cellGrid::flatIndex(const int x, const int y, const int z) {
    return x * numCells[1] * numCells[2] + y * numCells[2] + z;
}
void cellGrid::tdimIndex(int flatIndex, int& x, int& y, int& z) const {
        z = flatIndex % numCells[2];
        y = (flatIndex / numCells[2]) % numCells[1];
        x = flatIndex / (numCells[1] * numCells[2]);
    }


std::vector<std::array<int,3>> cellGrid::getNeighbors(int x, int y, int z) {
	std::vector<std::array<int,3>> neighbors;
    // Define the possible neighbor offsets

    for (int c = 0; c < 27; c++) {
        int nx = neighPos[c][0] + x;
		int ny = neighPos[c][1] + y;
		int nz = neighPos[c][2] + z;
        // Check if the neighbor is within bounds
        if (nx >= 0 && nx < numCells[0] && ny >= 0 && ny < numCells[1] && nz >= 0 && nz < numCells[2]) {
            neighbors.push_back({nx, ny, nz});
        }
    }
	return neighbors;
}

void cellGrid::precomputeNeighborIndices() {
    neighborIndices.resize(numCells[0] * numCells[1] * numCells[2]);

    for (int x = 0; x < numCells[0]; x++) {
        for (int y = 0; y < numCells[1]; y++) {
            for (int z = 0; z < numCells[2]; z++) {
                int fIndex = flatIndex(x, y, z);
                std::vector<std::array<int,3>> neighbors = getNeighbors(x, y, z);
                for (const auto& neighbor : neighbors) {
                    int neighborFlatIndex = flatIndex(neighbor[0], neighbor[1], neighbor[2]);
					// printf("cell %d,%d,%d, flatindex = %d, flatneigh = %d\n",x,y,z,fIndex,neighborFlatIndex);
                    neighborIndices[fIndex].push_back(neighborFlatIndex);
                }
            }
        }
    }

}


void mesh::fillBopt(){    
	if (do_reserve()){ return fillBopt_reserve();}


	const float invEc0 = 1/eCells[0];
	const float invEc1 = 1/eCells[1];
	const float invEc2 = 1/eCells[2];

	#pragma omp parallel
    {
		
    	std::vector<std::array<int, 3>> localInserts;

        #pragma omp for nowait
		for (int t = 0; t<nElements;t++){

			bool generalcase = true;
			std::vector<std::set<int> > cellcd;
				
			// test if in plane yz
			if (static_cast<int>((meshData.x0[t]-minSpc[0])*invEc0) == static_cast<int>((meshData.x1[t]-minSpc[0])*invEc0) && static_cast<int>((meshData.x0[t]-minSpc[0])*invEc0) == static_cast<int>((meshData.x2[t]-minSpc[0])*invEc0)){
				cellcd.resize(numCells[1]);
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				planeInterpolationYZ(Point2(meshData.y2[t],meshData.z2[t]),Point2(meshData.y0[t],meshData.z0[t]),eCells[0],cellcd);
				planeInterpolationYZ(Point2(meshData.y1[t],meshData.z1[t]),Point2(meshData.y0[t],meshData.z0[t]),eCells[0],cellcd);
				planeInterpolationYZ(Point2(meshData.y2[t],meshData.z2[t]),Point2(meshData.y1[t],meshData.z1[t]),eCells[0],cellcd);
				const int x = ((meshData.x0[t]-minSpc[0])*invEc0);
				for (int i = 0; i< cellcd.size();i++){
					if(cellcd[i].size() == 0) continue;
					
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[x][i][j].push_back(t);}
						localInserts.push_back({x, i, j});
					}
				} 
				generalcase = false;
					
			}

			// test if in plane xz
			if (generalcase == true)
			if (static_cast<int>((meshData.y0[t]-minSpc[1])*invEc1) == static_cast<int>((meshData.y1[t]-minSpc[1])*invEc1) && static_cast<int>((meshData.y0[t]-minSpc[1])*invEc1) == static_cast<int>((meshData.y2[t]-minSpc[1])*invEc1)){
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				planeInterpolationXZ(Point2(meshData.x2[t],meshData.z2[t]),Point2(meshData.x0[t],meshData.z0[t]),eCells[0],cellcd);
				planeInterpolationXZ(Point2(meshData.x1[t],meshData.z1[t]),Point2(meshData.x0[t],meshData.z0[t]),eCells[0],cellcd);
				planeInterpolationXZ(Point2(meshData.x2[t],meshData.z2[t]),Point2(meshData.x1[t],meshData.z1[t]),eCells[0],cellcd);
				const int y = ((meshData.y0[t]-minSpc[1])*invEc1);
				for (int i = 0; i< cellcd.size();i++){
					if(cellcd[i].size() == 0) continue;
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][y][j].push_back(t);}
						localInserts.push_back({i, y, j});
					}
				} 
				generalcase = false;
			}
			
			// test if in plane xy
			if (generalcase == true)
			if (static_cast<int>((meshData.z0[t]-minSpc[2])*invEc2) == static_cast<int>((meshData.z1[t]-minSpc[2])*invEc2) && static_cast<int>((meshData.z0[t]-minSpc[2])*invEc2) == static_cast<int>((meshData.z2[t]-minSpc[2])*invEc2)){
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				planeInterpolationXY(Point2(meshData.x2[t],meshData.y2[t]),Point2(meshData.x0[t],meshData.y0[t]),eCells[0],cellcd);
				planeInterpolationXY(Point2(meshData.x1[t],meshData.y1[t]),Point2(meshData.x0[t],meshData.y0[t]),eCells[0],cellcd);
				planeInterpolationXY(Point2(meshData.x2[t],meshData.y2[t]),Point2(meshData.x1[t],meshData.y1[t]),eCells[0],cellcd);
				const int z = ((meshData.z0[t]-minSpc[2])*invEc2);
				for (int i = 0; i< cellcd.size();i++){
					if(cellcd[i].size() == 0) continue;
					
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][j][z].push_back(t);}
						localInserts.push_back({i, j, z});
					}
				} 
				generalcase = false;
			}

			if(generalcase){
				float zvert[] = {meshData.z0[t], meshData.z1[t], meshData.z2[t]};

				auto [min, max] = std::minmax_element(zvert, zvert + 3);
				int mid = 3 - (min - zvert) - (max - zvert);

				// Extracting corresponding x and y values
				float midx, midy, midz;
				float minx, miny, minz;
				float maxx, maxy, maxz;

				if ((min - zvert) == 0) {
					minx = meshData.x0[t]; miny = meshData.y0[t]; minz = meshData.z0[t];
				} else if ((min - zvert) == 1) {
					minx = meshData.x1[t]; miny = meshData.y1[t]; minz = meshData.z1[t];
				} else {
					minx = meshData.x2[t]; miny = meshData.y2[t]; minz = meshData.z2[t];
				}

				if ((max - zvert) == 0) {
					maxx = meshData.x0[t]; maxy = meshData.y0[t]; maxz = meshData.z0[t];
				} else if ((max - zvert) == 1) {
					maxx = meshData.x1[t]; maxy = meshData.y1[t]; maxz = meshData.z1[t];
				} else {
					maxx = meshData.x2[t]; maxy = meshData.y2[t]; maxz = meshData.z2[t];
				}

				if (mid == 0) {
					midx = meshData.x0[t]; midy = meshData.y0[t]; midz = meshData.z0[t];
				} else if (mid == 1) {
					midx = meshData.x1[t]; midy = meshData.y1[t]; midz = meshData.z1[t];
				} else {
					midx = meshData.x2[t]; midy = meshData.y2[t]; midz = meshData.z2[t];
				}

				
				// maxz = *max ; minz = *min

				//min-zvert = min index
				//max-zvert = max index
				std::vector<Point2> intersections(2), intersectold(2);
				int k = static_cast<int>((minz-minSpc[2])*invEc2);
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				// checking first vertex
				{
				intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t);
				// if min,max:{0,1}->mid=2 / min,max:{0,2}->mid=1 / min,max:{1,2}-> mid = 0
				
				planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd); 
				planeInterpolationXY(Point2(minx,miny),intersections[0],eCells[0],cellcd); 
				planeInterpolationXY(Point2(minx,miny),intersections[1],eCells[0],cellcd); 
				
				//fill buffer from cellcd 
				for (int i = 0; i< cellcd.size();i++){
					if(cellcd[i].size() == 0) continue;
					
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = minz; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][j][k].push_back(t);}
						localInserts.push_back({i, j, k});
					}
				}
				
				intersectold[0] = intersections[0]; //Salva a interseccao de z-1
				intersectold[1] = intersections[1];
				}
				// from bot+1 to mid-1
				for(k = static_cast<int>((minz-minSpc[2])*invEc2)+1; k<static_cast<int>((midz-minSpc[2])*invEc2)-1;k++){
				// for(k = static_cast<int>((minz-minSpc[2])*invEc2)+1; k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

					if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
					
						for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
						planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
						planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					
				//fill buffer from cellcd 
						for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
								
								localInserts.push_back({i, j, k});
							}
						}
						intersectold[0] = intersections[0]; //Salva a interseccao de z-1
						intersectold[1] = intersections[1];
						}
				}		
				// checking mid vertex 
				{
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				if (intersectTrianglePlaneZ(midz-eps, intersections,t)){
					planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					
					intersectold[0] = intersections[0]; //Salva a interseccao de z-1
					intersectold[1] = intersections[1];

					k = static_cast<int>((midz-minSpc[2])*invEc2);
					for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){

								localInserts.push_back({i, j, k});
							}
						}
				}
				}
				// from mid to top 
				for(k = static_cast<int>((midz-minSpc[2])*invEc2); k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

					if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
					
						for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
						planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
						planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
						//fill buffer from cellcd 
						for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
							// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
							// for(int_fast32_t j = minz; j<=*max;++j){
								// #pragma omp critical		
								// {cellBuffer[i][j][k].push_back(t);}
								localInserts.push_back({i, j, k});
							}
						}
						intersectold[0] = intersections[0]; //Salva a interseccao de z-1
						intersectold[1] = intersections[1];
						}
				}
				// checking last vertex
				{
					for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
					k =static_cast<int>((maxz-minSpc[2])*invEc2);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd); 
					planeInterpolationXY(Point2(maxx,maxy),intersectold[0],eCells[0],cellcd); 
					planeInterpolationXY(Point2(maxx,maxy),intersectold[1],eCells[0],cellcd); 
					
					//fill buffer from cellcd 
					for (int i = 0; i< cellcd.size();i++){
						if(cellcd[i].size() == 0) continue;
						
						for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
						// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
						// for(int_fast32_t j = minz; j<=*max;++j){
							// #pragma omp critical
							// {cellBuffer[i][j][k].push_back(t);}
							localInserts.push_back({i, j, k});
						}
					}
				}
			}
			#pragma omp critical
			{
				for (const auto& insert : localInserts) {
					int i = insert[0];
					int j = insert[1];
					int k = insert[2];
					cellBuffer[i][j][k].push_back(t);
				}
				localInserts.clear();
			}	
		}
	}
}

void mesh::fillBopt_flat(){    
	
	if (do_reserve()){ return fillBopt_reserveflat();}

	const float invEc0 = 1/eCells[0];
	const float invEc1 = 1/eCells[1];
	const float invEc2 = 1/eCells[2];

	#pragma omp parallel
    {
		
    	std::vector<unsigned int> localInserts;

        #pragma omp for 
		for (int t = 0; t<nElements;t++){

			bool generalcase = true;
			std::vector<std::set<int> > cellcd;

			float zvert[] = {meshData.z0[t], meshData.z1[t], meshData.z2[t]};

			auto [min, max] = std::minmax_element(zvert, zvert + 3);
			int mid = 3 - (min - zvert) - (max - zvert);

			// Extracting corresponding x and y values
			float midx, midy, midz;
			float minx, miny, minz;
			float maxx, maxy, maxz;

			if ((min - zvert) == 0) {
				minx = meshData.x0[t]; miny = meshData.y0[t]; minz = meshData.z0[t];
			} else if ((min - zvert) == 1) {
				minx = meshData.x1[t]; miny = meshData.y1[t]; minz = meshData.z1[t];
			} else {
				minx = meshData.x2[t]; miny = meshData.y2[t]; minz = meshData.z2[t];
			}

			if ((max - zvert) == 0) {
				maxx = meshData.x0[t]; maxy = meshData.y0[t]; maxz = meshData.z0[t];
			} else if ((max - zvert) == 1) {
				maxx = meshData.x1[t]; maxy = meshData.y1[t]; maxz = meshData.z1[t];
			} else {
				maxx = meshData.x2[t]; maxy = meshData.y2[t]; maxz = meshData.z2[t];
			}

			if (mid == 0) {
				midx = meshData.x0[t]; midy = meshData.y0[t]; midz = meshData.z0[t];
			} else if (mid == 1) {
				midx = meshData.x1[t]; midy = meshData.y1[t]; midz = meshData.z1[t];
			} else {
				midx = meshData.x2[t]; midy = meshData.y2[t]; midz = meshData.z2[t];
			}

			// test if in plane yz
			if (static_cast<int>((minx-minSpc[0])*invEc0) == static_cast<int>((midx-minSpc[0])*invEc0) && static_cast<int>((minx-minSpc[0])*invEc0) == static_cast<int>((maxx-minSpc[0])*invEc0)){
				cellcd.resize(numCells[1]);
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				planeInterpolationYZ(Point2(maxy,maxz),Point2(miny,minz),eCells[0],cellcd);
				planeInterpolationYZ(Point2(midy,midz),Point2(miny,minz),eCells[0],cellcd);
				planeInterpolationYZ(Point2(maxy,maxz),Point2(midy,midz),eCells[0],cellcd);
				const int x = ((minx-minSpc[0])*invEc0);
				const auto [ymin,ymax] = std::minmax({static_cast<int>((miny-minSpc[1])*invEc1),static_cast<int>((midy-minSpc[1])*invEc1),static_cast<int>((maxy-minSpc[1])*invEc1)});
				for (int i = ymin; i<= ymax;i++){					
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[x][i][j].push_back(t);}
						localInserts.push_back(flatIndex(x, i, j));
					}
				} 
				generalcase = false;
					
			}


			// test if in plane xz
			if (generalcase == true)
			if (static_cast<int>((miny-minSpc[1])*invEc1) == static_cast<int>((midy-minSpc[1])*invEc1) && static_cast<int>((miny-minSpc[1])*invEc1) == static_cast<int>((maxy-minSpc[1])*invEc1)){
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				planeInterpolationXZ(Point2(maxx,maxz),Point2(minx,minz),eCells[0],cellcd);
				planeInterpolationXZ(Point2(midx,midz),Point2(minx,minz),eCells[0],cellcd);
				planeInterpolationXZ(Point2(maxx,maxz),Point2(midx,midz),eCells[0],cellcd);
				const int y = ((miny-minSpc[1])*invEc1);
				const auto [xmin,xmax] = std::minmax({static_cast<int>((minx-minSpc[0])*invEc0),static_cast<int>((midx-minSpc[0])*invEc0),static_cast<int>((maxx-minSpc[0])*invEc0)});
				for (int i = xmin; i<= xmax;i++){		
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][y][j].push_back(t);}
						localInserts.push_back(flatIndex(i, y, j));
					}
				} 
				generalcase = false;
			}
			
			// test if in plane xy
			if (generalcase == true)
			if (static_cast<int>((minz-minSpc[2])*invEc2) == static_cast<int>((maxz-minSpc[2])*invEc2)){
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				planeInterpolationXY(Point2(maxx,maxy),Point2(minx,miny),eCells[0],cellcd);
				planeInterpolationXY(Point2(midx,midy),Point2(minx,miny),eCells[0],cellcd);
				planeInterpolationXY(Point2(maxx,maxy),Point2(midx,midy),eCells[0],cellcd);
				const int z = ((minz-minSpc[2])*invEc2);
				const auto [xmin,xmax] = std::minmax({static_cast<int>((minx-minSpc[0])*invEc0),static_cast<int>((midx-minSpc[0])*invEc0),static_cast<int>((maxx-minSpc[0])*invEc0)});
				for (int i = xmin; i<= xmax;i++){				
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][j][z].push_back(t);}
						localInserts.push_back(flatIndex(i, j, z));
					}
				} 
				generalcase = false;
			}

			if(generalcase){
				cellcd.resize(numCells[0]);
				// maxz = *max ; minz = *min

				//min-zvert = min index
				//max-zvert = max index
				std::vector<Point2> intersections(2), intersectold(2);
				int k = static_cast<int>((minz-minSpc[2])*invEc2);
				intersectTrianglePlaneZ(minz, intersections,t);
				intersectold[0] = intersections[0]; //Salva a interseccao de z-1
				intersectold[1] = intersections[1];
				// from bot+1 to mid-1
				for(k; k<static_cast<int>((midz-minSpc[2])*invEc2);k++){
				// for(k = static_cast<int>((minz-minSpc[2])*invEc2)+1; k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

					if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
					
						for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
						planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
						planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					
				//fill buffer from cellcd 
						for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
								
								localInserts.push_back(flatIndex(i, j, k));
							}
						}
						intersectold[0] = intersections[0]; //Salva a interseccao de z-1
						intersectold[1] = intersections[1];
						}
				}		
				// checking mid vertex 
				{
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				if (intersectTrianglePlaneZ(midz-eps, intersections,t)){
					planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					
					intersectold[0] = intersections[0]; //Salva a interseccao de z-1
					intersectold[1] = intersections[1];

					k = static_cast<int>((midz-minSpc[2])*invEc2);
					for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){

								localInserts.push_back(flatIndex(i, j, k));
							}
						}
				}
				}
				// from mid to top 
				for(k = static_cast<int>((midz-minSpc[2])*invEc2); k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

					if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
					
						for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
						planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
						planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
						//fill buffer from cellcd 
						for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
							// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
							// for(int_fast32_t j = minz; j<=*max;++j){
								// #pragma omp critical		
								// {cellBuffer[i][j][k].push_back(t);}
								localInserts.push_back(flatIndex(i, j, k));
							}
						}
						intersectold[0] = intersections[0]; //Salva a interseccao de z-1
						intersectold[1] = intersections[1];
						}
				}
				
				// checking last vertex
				{
					for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
					k =static_cast<int>((maxz-minSpc[2])*invEc2);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd); 
					planeInterpolationXY(Point2(maxx,maxy),intersectold[0],eCells[0],cellcd); 
					planeInterpolationXY(Point2(maxx,maxy),intersectold[1],eCells[0],cellcd); 
					
					//fill buffer from cellcd 
					for (int i = 0; i< cellcd.size();i++){
						if(cellcd[i].size() == 0) continue;
						
						for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
						// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
						// for(int_fast32_t j = minz; j<=*max;++j){
							// #pragma omp critical
							// {cellBuffer[i][j][k].push_back(t);}
							localInserts.push_back(flatIndex(i, j, k));
						}
					}
				}
			}
			#pragma omp critical
			{			
			for (const auto& insert : localInserts) {
					cellBufferFlattened[insert].push_back(t);
			}
			localInserts.clear();
			}
		}

	}

		
}

void mesh::fillBopt_reserveflat(){    
	

	const float invEc0 = 1/eCells[0];
	const float invEc1 = 1/eCells[1];
	const float invEc2 = 1/eCells[2];

	#pragma omp parallel for 
		for (int t = 0; t<nElements;t++){

			bool generalcase = true;
			std::vector<std::set<int> > cellcd;

			float zvert[] = {meshData.z0[t], meshData.z1[t], meshData.z2[t]};

			auto [min, max] = std::minmax_element(zvert, zvert + 3);
			int mid = 3 - (min - zvert) - (max - zvert);

			// Extracting corresponding x and y values
			float midx, midy, midz;
			float minx, miny, minz;
			float maxx, maxy, maxz;

			if ((min - zvert) == 0) {
				minx = meshData.x0[t]; miny = meshData.y0[t]; minz = meshData.z0[t];
			} else if ((min - zvert) == 1) {
				minx = meshData.x1[t]; miny = meshData.y1[t]; minz = meshData.z1[t];
			} else {
				minx = meshData.x2[t]; miny = meshData.y2[t]; minz = meshData.z2[t];
			}

			if ((max - zvert) == 0) {
				maxx = meshData.x0[t]; maxy = meshData.y0[t]; maxz = meshData.z0[t];
			} else if ((max - zvert) == 1) {
				maxx = meshData.x1[t]; maxy = meshData.y1[t]; maxz = meshData.z1[t];
			} else {
				maxx = meshData.x2[t]; maxy = meshData.y2[t]; maxz = meshData.z2[t];
			}

			if (mid == 0) {
				midx = meshData.x0[t]; midy = meshData.y0[t]; midz = meshData.z0[t];
			} else if (mid == 1) {
				midx = meshData.x1[t]; midy = meshData.y1[t]; midz = meshData.z1[t];
			} else {
				midx = meshData.x2[t]; midy = meshData.y2[t]; midz = meshData.z2[t];
			}

			// test if in plane yz
			if (static_cast<int>((minx-minSpc[0])*invEc0) == static_cast<int>((midx-minSpc[0])*invEc0) && static_cast<int>((minx-minSpc[0])*invEc0) == static_cast<int>((maxx-minSpc[0])*invEc0)){
				cellcd.resize(numCells[1]);
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				planeInterpolationYZ(Point2(maxy,maxz),Point2(miny,minz),eCells[0],cellcd);
				planeInterpolationYZ(Point2(midy,midz),Point2(miny,minz),eCells[0],cellcd);
				planeInterpolationYZ(Point2(maxy,maxz),Point2(midy,midz),eCells[0],cellcd);
				const int x = ((minx-minSpc[0])*invEc0);
				const auto [ymin,ymax] = std::minmax({static_cast<int>((miny-minSpc[1])*invEc1),static_cast<int>((midy-minSpc[1])*invEc1),static_cast<int>((maxy-minSpc[1])*invEc1)});
				for (int i = ymin; i<= ymax;i++){					
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[x][i][j].push_back(t);}
						cellBufferFlattened[flatIndex(x, i, j)].push_back(t);
					}
				} 
				generalcase = false;
					
			}


			// test if in plane xz
			if (generalcase == true)
			if (static_cast<int>((miny-minSpc[1])*invEc1) == static_cast<int>((midy-minSpc[1])*invEc1) && static_cast<int>((miny-minSpc[1])*invEc1) == static_cast<int>((maxy-minSpc[1])*invEc1)){
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				planeInterpolationXZ(Point2(maxx,maxz),Point2(minx,minz),eCells[0],cellcd);
				planeInterpolationXZ(Point2(midx,midz),Point2(minx,minz),eCells[0],cellcd);
				planeInterpolationXZ(Point2(maxx,maxz),Point2(midx,midz),eCells[0],cellcd);
				const int y = ((miny-minSpc[1])*invEc1);
				const auto [xmin,xmax] = std::minmax({static_cast<int>((minx-minSpc[0])*invEc0),static_cast<int>((midx-minSpc[0])*invEc0),static_cast<int>((maxx-minSpc[0])*invEc0)});
				for (int i = xmin; i<= xmax;i++){		
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][y][j].push_back(t);}
						cellBufferFlattened[flatIndex(i, y, j)].push_back(t);
					}
				} 
				generalcase = false;
			}
			
			// test if in plane xy
			if (generalcase == true)
			if (static_cast<int>((minz-minSpc[2])*invEc2) == static_cast<int>((maxz-minSpc[2])*invEc2)){
				
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				cellcd.resize(numCells[0]);
				planeInterpolationXY(Point2(maxx,maxy),Point2(minx,miny),eCells[0],cellcd);
				planeInterpolationXY(Point2(midx,midy),Point2(minx,miny),eCells[0],cellcd);
				planeInterpolationXY(Point2(maxx,maxy),Point2(midx,midy),eCells[0],cellcd);
				const int z = ((minz-minSpc[2])*invEc2);
				const auto [xmin,xmax] = std::minmax({static_cast<int>((minx-minSpc[0])*invEc0),static_cast<int>((midx-minSpc[0])*invEc0),static_cast<int>((maxx-minSpc[0])*invEc0)});
				for (int i = xmin; i<= xmax;i++){				
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = *min; j<=*max;++j){
						// #pragma omp critical
						// {cellBuffer[i][j][z].push_back(t);}
						cellBufferFlattened[flatIndex(i, j, z)].push_back(t);
					}
				} 
				generalcase = false;
			}

			if(generalcase){
				cellcd.resize(numCells[0]);
				// maxz = *max ; minz = *min

				//min-zvert = min index
				//max-zvert = max index
				std::vector<Point2> intersections(2), intersectold(2);
				int k = static_cast<int>((minz-minSpc[2])*invEc2);
				intersectTrianglePlaneZ(minz, intersections,t);
				intersectold[0] = intersections[0]; //Salva a interseccao de z-1
				intersectold[1] = intersections[1];
				// from bot+1 to mid-1
				for(k; k<static_cast<int>((midz-minSpc[2])*invEc2);k++){
				// for(k = static_cast<int>((minz-minSpc[2])*invEc2)+1; k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

					if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
					
						for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
						planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
						planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					
				//fill buffer from cellcd 
						for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
								
								cellBufferFlattened[flatIndex(i, j, k)].push_back(t);
							}
						}
						intersectold[0] = intersections[0]; //Salva a interseccao de z-1
						intersectold[1] = intersections[1];
						}
				}		
				// checking mid vertex 
				{
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				if (intersectTrianglePlaneZ(midz-eps, intersections,t)){
					planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					
					intersectold[0] = intersections[0]; //Salva a interseccao de z-1
					intersectold[1] = intersections[1];

					k = static_cast<int>((midz-minSpc[2])*invEc2);
					for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){

								cellBufferFlattened[flatIndex(i, j, k)].push_back(t);
							}
						}
				}
				}
				// from mid to top 
				for(k = static_cast<int>((midz-minSpc[2])*invEc2); k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

					if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
					
						for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
						planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
						planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
						planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
						//fill buffer from cellcd 
						for (int i = 0; i< cellcd.size();i++){
							if(cellcd[i].size() == 0) continue;
							
							for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
							// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
							// for(int_fast32_t j = minz; j<=*max;++j){
								// #pragma omp critical		
								// {cellBuffer[i][j][k].push_back(t);}
								cellBufferFlattened[flatIndex(i, j, k)].push_back(t);
							}
						}
						intersectold[0] = intersections[0]; //Salva a interseccao de z-1
						intersectold[1] = intersections[1];
						}
				}
				
				// checking last vertex
				{
					for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
					k =static_cast<int>((maxz-minSpc[2])*invEc2);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd); 
					planeInterpolationXY(Point2(maxx,maxy),intersectold[0],eCells[0],cellcd); 
					planeInterpolationXY(Point2(maxx,maxy),intersectold[1],eCells[0],cellcd); 
					
					//fill buffer from cellcd 
					for (int i = 0; i< cellcd.size();i++){
						if(cellcd[i].size() == 0) continue;
						
						for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
						// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
						// for(int_fast32_t j = minz; j<=*max;++j){
							// #pragma omp critical
							// {cellBuffer[i][j][k].push_back(t);}
							cellBufferFlattened[flatIndex(i, j, k)].push_back(t);
						}
					}
				}
			}

	}

		
}

void mesh::fillCandidates_flat() {
    #pragma omp parallel for
    for (int i = 0; i < numCells[0] * numCells[1] * numCells[2]; ++i) {
        auto& local = candidatesSet_flat[i];

        for (const int n : neighborIndices[i]) {
            const auto& cb = cellBufferFlattened[n];

            if (!cb.empty()) {
                local.insert(cb.begin(), cb.end());
            }
        }
    }
}


void mesh::fillBopt_reserve(){    
	const float invEc0 = 1/eCells[0];
	const float invEc1 = 1/eCells[1];
	const float invEc2 = 1/eCells[2];

	#pragma omp parallel for
	for (int t = 0; t<nElements;t++){

		bool generalcase = true;
		std::vector<std::set<int> > cellcd;
			
		// test if in plane yz
		if (static_cast<int>((meshData.x0[t]-minSpc[0])*invEc0) == static_cast<int>((meshData.x1[t]-minSpc[0])*invEc0) && static_cast<int>((meshData.x0[t]-minSpc[0])*invEc0) == static_cast<int>((meshData.x2[t]-minSpc[0])*invEc0)){
			cellcd.resize(numCells[1]);
			
			for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
			planeInterpolationYZ(Point2(meshData.y2[t],meshData.z2[t]),Point2(meshData.y0[t],meshData.z0[t]),eCells[0],cellcd);
			planeInterpolationYZ(Point2(meshData.y1[t],meshData.z1[t]),Point2(meshData.y0[t],meshData.z0[t]),eCells[0],cellcd);
			planeInterpolationYZ(Point2(meshData.y2[t],meshData.z2[t]),Point2(meshData.y1[t],meshData.z1[t]),eCells[0],cellcd);
			const int x = ((meshData.x0[t]-minSpc[0])*invEc0);
			for (int i = 0; i< cellcd.size();i++){
				if(cellcd[i].size() == 0) continue;
				
				for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
				// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
				// for(int_fast32_t j = *min; j<=*max;++j){
					// #pragma omp critical
					{cellBuffer[x][i][j].push_back(t);}
				}
			} 
			generalcase = false;
				
		}

		// test if in plane xz
		if (generalcase == true)
		if (static_cast<int>((meshData.y0[t]-minSpc[1])*invEc1) == static_cast<int>((meshData.y1[t]-minSpc[1])*invEc1) && static_cast<int>((meshData.y0[t]-minSpc[1])*invEc1) == static_cast<int>((meshData.y2[t]-minSpc[1])*invEc1)){
			
			for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
			cellcd.resize(numCells[0]);
			planeInterpolationXZ(Point2(meshData.x2[t],meshData.z2[t]),Point2(meshData.x0[t],meshData.z0[t]),eCells[0],cellcd);
			planeInterpolationXZ(Point2(meshData.x1[t],meshData.z1[t]),Point2(meshData.x0[t],meshData.z0[t]),eCells[0],cellcd);
			planeInterpolationXZ(Point2(meshData.x2[t],meshData.z2[t]),Point2(meshData.x1[t],meshData.z1[t]),eCells[0],cellcd);
			const int y = ((meshData.y0[t]-minSpc[1])*invEc1);
			for (int i = 0; i< cellcd.size();i++){
				if(cellcd[i].size() == 0) continue;
				for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
				// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
				// for(int_fast32_t j = *min; j<=*max;++j){
					// #pragma omp critical
					{cellBuffer[i][y][j].push_back(t);}
				}
			} 
			generalcase = false;
		}
		
		// test if in plane xy
		if (generalcase == true)
		if (static_cast<int>((meshData.z0[t]-minSpc[2])*invEc2) == static_cast<int>((meshData.z1[t]-minSpc[2])*invEc2) && static_cast<int>((meshData.z0[t]-minSpc[2])*invEc2) == static_cast<int>((meshData.z2[t]-minSpc[2])*invEc2)){
			
			for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
			cellcd.resize(numCells[0]);
			planeInterpolationXY(Point2(meshData.x2[t],meshData.y2[t]),Point2(meshData.x0[t],meshData.y0[t]),eCells[0],cellcd);
			planeInterpolationXY(Point2(meshData.x1[t],meshData.y1[t]),Point2(meshData.x0[t],meshData.y0[t]),eCells[0],cellcd);
			planeInterpolationXY(Point2(meshData.x2[t],meshData.y2[t]),Point2(meshData.x1[t],meshData.y1[t]),eCells[0],cellcd);
			const int z = ((meshData.z0[t]-minSpc[2])*invEc2);
			for (int i = 0; i< cellcd.size();i++){
				if(cellcd[i].size() == 0) continue;
				
				for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
				// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
				// for(int_fast32_t j = *min; j<=*max;++j){
					// #pragma omp critical
					{cellBuffer[i][j][z].push_back(t);}
				}
			} 
			generalcase = false;
		}

		if(generalcase){
			float zvert[] = {meshData.z0[t], meshData.z1[t], meshData.z2[t]};

			auto [min, max] = std::minmax_element(zvert,zvert+3);
			
			int mid = 3 - (min - zvert) - (max - zvert);

			// Extracting corresponding x and y values
			float midx, midy, midz;
			float minx, miny, minz;
			float maxx, maxy, maxz;

			if ((min - zvert) == 0) {
				minx = meshData.x0[t]; miny = meshData.y0[t]; minz = meshData.z0[t];
			} else if ((min - zvert) == 1) {
				minx = meshData.x1[t]; miny = meshData.y1[t]; minz = meshData.z1[t];
			} else {
				minx = meshData.x2[t]; miny = meshData.y2[t]; minz = meshData.z2[t];
			}

			if ((max - zvert) == 0) {
				maxx = meshData.x0[t]; maxy = meshData.y0[t]; maxz = meshData.z0[t];
			} else if ((max - zvert) == 1) {
				maxx = meshData.x1[t]; maxy = meshData.y1[t]; maxz = meshData.z1[t];
			} else {
				maxx = meshData.x2[t]; maxy = meshData.y2[t]; maxz = meshData.z2[t];
			}

			if (mid == 0) {
				midx = meshData.x0[t]; midy = meshData.y0[t]; midz = meshData.z0[t];
			} else if (mid == 1) {
				midx = meshData.x1[t]; midy = meshData.y1[t]; midz = meshData.z1[t];
			} else {
				midx = meshData.x2[t]; midy = meshData.y2[t]; midz = meshData.z2[t];
			}

			// maxz = *max ; minz = *min

			//min-zvert = min index
			//max-zvert = max index
			std::vector<Point2> intersections(2), intersectold(2);
			int k = static_cast<int>((minz-minSpc[2])*invEc2);
			
			for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
			cellcd.resize(numCells[0]);
			// checking first vertex
			{
			intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t);
			// if min,max:{0,1}->mid=2 / min,max:{0,2}->mid=1 / min,max:{1,2}-> mid = 0
			
			planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd); 
			planeInterpolationXY(Point2(minx,miny),intersections[0],eCells[0],cellcd); 
			planeInterpolationXY(Point2(minx,miny),intersections[1],eCells[0],cellcd); 
			
			if (static_cast<int>((midz-minSpc[2])*invEc2) == k){
				planeInterpolationXY(Point2(minx,miny),Point2(midx,midy),eCells[0],cellcd);
				planeInterpolationXY(Point2(midx,midy),intersections[0],eCells[0],cellcd); 
				planeInterpolationXY(Point2(midx,midy),intersections[1],eCells[0],cellcd);
			}

			//fill buffer from cellcd 
			for (int i = 0; i< cellcd.size();i++){
				if(cellcd[i].size() == 0) continue;
				
				for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
				// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
				// for(int_fast32_t j = minz; j<=*max;++j){
					// #pragma omp critical
					{cellBuffer[i][j][k].push_back(t);}
				}
			}
			
			intersectold[0] = intersections[0]; //Salva a interseccao de z-1
			intersectold[1] = intersections[1];
			}
			// from bot+1 to mid-1
			for(k = static_cast<int>((minz-minSpc[2])*invEc2)+1; k<static_cast<int>((midz-minSpc[2])*invEc2)-1;k++){
			// for(k = static_cast<int>((minz-minSpc[2])*invEc2)+1; k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

				if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
				
					for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
					planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
				
			//fill buffer from cellcd 
					for (int i = 0; i< cellcd.size();i++){
						if(cellcd[i].size() == 0) continue;
						
						for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
							cellBuffer[i][j][k].push_back(t);
						}
					}
					intersectold[0] = intersections[0]; //Salva a interseccao de z-1
					intersectold[1] = intersections[1];
					}
			}		
			// checking mid vertex 
			{
			for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
			if (intersectTrianglePlaneZ(midz-eps, intersections,t)){
				planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
				planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
				planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
				planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
				planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
				planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
				
				intersectold[0] = intersections[0]; //Salva a interseccao de z-1
				intersectold[1] = intersections[1];

				k = static_cast<int>((midz-minSpc[2])*invEc2);
				for (int i = 0; i< cellcd.size();i++){
						if(cellcd[i].size() == 0) continue;
						
						for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){

							cellBuffer[i][j][k].push_back(t);
						}
					}
			}
			}
			// from mid to top 
			for(k = static_cast<int>((midz-minSpc[2])*invEc2); k<static_cast<int>((maxz-minSpc[2])*invEc2);k++){

				if (intersectTrianglePlaneZ(static_cast<float>(k+1)*eCells[2]+minSpc[2], intersections,t)){ // if intersect triang/plane z
				
					for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
					planeInterpolationXY(intersections[0],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd);
					planeInterpolationXY(intersectold[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[0],intersectold[1],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[0],eCells[0],cellcd);
					planeInterpolationXY(intersections[1],intersectold[1],eCells[0],cellcd);
					//fill buffer from cellcd 
					for (int i = 0; i< cellcd.size();i++){
						if(cellcd[i].size() == 0) continue;
						
						for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
						// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
						// for(int_fast32_t j = minz; j<=*max;++j){
							// #pragma omp critical		
							{cellBuffer[i][j][k].push_back(t);}
						}
					}
					intersectold[0] = intersections[0]; //Salva a interseccao de z-1
					intersectold[1] = intersections[1];
					}
			}
			// checking last vertex
			{
				for (auto &p : cellcd) p.clear();	// clear cellcd previous iteration
				k =static_cast<int>((maxz-minSpc[2])*invEc2);
				intersectTrianglePlaneZ(static_cast<float>(k)*eCells[2]+minSpc[2], intersections,t);
				planeInterpolationXY(intersections[0],intersections[1],eCells[0],cellcd); 
				planeInterpolationXY(Point2(maxx,maxy),intersections[0],eCells[0],cellcd); 
				planeInterpolationXY(Point2(maxx,maxy),intersections[1],eCells[0],cellcd); 
				
				if (static_cast<int>((midz-minSpc[2])*invEc2) == k){
					planeInterpolationXY(Point2(maxx,maxy),Point2(midx,midy),eCells[0],cellcd);
					planeInterpolationXY(Point2(midx,midy),intersections[0],eCells[0],cellcd); 
					planeInterpolationXY(Point2(midx,midy),intersections[1],eCells[0],cellcd);
				}
				//fill buffer from cellcd 
				for (int i = 0; i< cellcd.size();i++){
					if(cellcd[i].size() == 0) continue;
					
					for(int j = *cellcd[i].begin(); j<=*--cellcd[i].end();++j){
					// const auto [min, max] = std::minmax_element(begin(cellcd[i]), end(cellcd[i]));
					// for(int_fast32_t j = minz; j<=*max;++j){
						// #pragma omp critical
							cellBuffer[i][j][k].push_back(t);
					}
				}
			}
		}
	
	}
}

void mesh::fillBuffer() {
    #pragma omp parallel
    {
        std::vector<std::array<int, 3>> localInserts;

        #pragma omp for
        for (int t = 0; t < nElements; t++) {
            float cen[3];
            float normalx, normaly, normalz;
            float vx[3], vy[3], vz[3];
            int bound[3][2] = {0, 0, 0, 0, 0, 0};

            for (int n = 0; n < 3; n++) {
                vx[n] = getNodePositionX(t, n);
                vy[n] = getNodePositionY(t, n);
                vz[n] = getNodePositionZ(t, n);
            }

            getBounds(bound, vx, vy, vz);
            normalx = getNormalX(t);
            normaly = getNormalY(t);
            normalz = getNormalZ(t);

            for (int i = bound[0][0]; i <= bound[0][1]; i++) {
                if (i < 0 || (i >= numCells[0])) continue;
                for (int j = bound[1][0]; j <= bound[1][1]; j++) {
                    if (j < 0 || (j >= numCells[1])) continue;
                    for (int k = bound[2][0]; k <= bound[2][1]; k++) {
                        if (k < 0 || (k >= numCells[2])) continue;

                        cen[0] = getCellCenterX(i);
                        cen[1] = getCellCenterY(j);
                        cen[2] = getCellCenterZ(k);

                        if (triangleAABB_intersection(cen, normalx, normaly, normalz, vx[0], vy[0], vz[0], vx[1], vy[1], vz[1], vx[2], vy[2], vz[2])) {
                            localInserts.push_back({i, j, k});
                        }
                    }
                }
            }

            #pragma omp critical
            {
                for (const auto& insert : localInserts) {
                    int i = insert[0];
                    int j = insert[1];
                    int k = insert[2];
                    cellBuffer[i][j][k].push_back(t);
                }
                localInserts.clear();
            }
        }
    }
}
void mesh::fillCandidates() {
	int x = std::max_element(numCells,numCells+3) - numCells;
	int xmin = std::min_element(numCells,numCells+3) - numCells;
	int xmid = (3 - x - xmin);
	if (numCells[xmin] > 12 && numCells[xmid] > 12)
	// 3 2 2
    #pragma omp parallel num_threads(12)
    {	
        int numThreads = omp_get_num_threads();
        int tid = omp_get_thread_num(); // Thread ID
		// printf("\nthread = %d",tid);
		
		int threaddomain[3];

		threaddomain[x] = static_cast<int>(tid/4);
		threaddomain[(x+1)%3] = static_cast<int>(tid/2)%2;
		threaddomain[(x+2)%3] = static_cast<int>(tid%2);
		
		int chunkSize[3];

        chunkSize[x] = (numCells[x] + 2) / 3; // Ceiling division

        chunkSize[(x+1)%3] = (numCells[(x+1)%3] + 1) /2; // Ceiling division
		chunkSize[(x+2)%3] = (numCells[(x+2)%3] + 1) /2; // Ceiling division

		// printf("thread: %d , domain: %d,%d,%d\n",tid, threaddomain[0],threaddomain[1],threaddomain[2]);

        int startIdxI = threaddomain[0] * chunkSize[0];
        int startIdxJ = threaddomain[1] * chunkSize[1];
        int startIdxK = threaddomain[2] * chunkSize[2];

        int endIdxI = std::min((threaddomain[0] + 1) * chunkSize[0], numCells[0]);
        int endIdxJ = std::min((threaddomain[1] + 1) * chunkSize[1], numCells[1]);
        int endIdxK = std::min((threaddomain[2] + 1) * chunkSize[2], numCells[2]);
		
        for (int i = startIdxI; i < endIdxI; i++) {
			
            const int startX = std::max(i - 1, 0);
            const int endX = std::min(i + 1, numCells[0] - 1);

            for (int j = startIdxJ; j < endIdxJ; j++) {

                    const int startY = std::max(j - 1, 0);
                    const int endY = std::min(j + 1, numCells[1] - 1);

                for (int k = startIdxK; k < endIdxK; k++) {
					
                    const int startZ = std::max(k - 1, 0);
                    const int endZ = std::min(k + 1, numCells[2] - 1);
					
                    auto& cell = cellBuffer[i][j][k];

                    for (int ii = startX; ii <= endX; ++ii) {
                        for (int jj = startY; jj <= endY; ++jj) {
                            for (int kk = startZ; kk <= endZ; ++kk) {
								//#pragma omp critical
                                
                                {candidates_set[i][j][k].insert(cellBuffer[ii][jj][kk].begin(), cellBuffer[ii][jj][kk].end());}
                            }
                        }
                    }
					
					// cellBuffer[i][j][k].clear();
                }
            }
        }
    }
	else if (numCells[xmin] < 12 && numCells[xmid] > 12)
	// 4 3 1
	#pragma omp parallel num_threads(12)
    {
        int numThreads = omp_get_num_threads();
        int tid = omp_get_thread_num(); // Thread ID
		// printf("\nthread = %d",tid);
		
		int threaddomain[3];
		threaddomain[x] = static_cast<int>(tid %4);
		threaddomain[xmid] = static_cast<int>(tid %3);
		threaddomain[xmin] = 0;

		int chunkSize[3];

		chunkSize[x] = (numCells[x] + 3) / 4; 
		chunkSize[xmid] = (numCells[xmid] + 2 ) / 3;
		chunkSize[xmin] = numCells[xmin];

		// printf("thread: %d , domain: %d,%d,%d\n",tid, threaddomain[0],threaddomain[1],threaddomain[2]);

        int startIdxI = threaddomain[0] * chunkSize[0];
        int startIdxJ = threaddomain[1] * chunkSize[1];
        int startIdxK = threaddomain[2] * chunkSize[2];

        int endIdxI = std::min((threaddomain[0] + 1) * chunkSize[0], numCells[0]);
        int endIdxJ = std::min((threaddomain[1] + 1) * chunkSize[1], numCells[1]);
        int endIdxK = std::min((threaddomain[2] + 1) * chunkSize[2], numCells[2]);
		
        for (int i = startIdxI; i < endIdxI; i++) {
			
            const int startX = std::max(i - 1, 0);
            const int endX = std::min(i + 1, numCells[0] - 1);

            for (int j = startIdxJ; j < endIdxJ; j++) {

                    const int startY = std::max(j - 1, 0);
                    const int endY = std::min(j + 1, numCells[1] - 1);

                for (int k = startIdxK; k < endIdxK; k++) {
					
                    const int startZ = std::max(k - 1, 0);
                    const int endZ = std::min(k + 1, numCells[2] - 1);
					
                    auto& cell = cellBuffer[i][j][k];

                    for (int ii = startX; ii <= endX; ++ii) {
                        for (int jj = startY; jj <= endY; ++jj) {
                            for (int kk = startZ; kk <= endZ; ++kk) {
								//#pragma omp critical
                                
                                {candidates_set[i][j][k].insert(cellBuffer[ii][jj][kk].begin(), cellBuffer[ii][jj][kk].end());}
                            }
                        }
                    }
					
					// cellBuffer[i][j][k].clear();
                }
            }
        }
    }
	else if (numCells[xmid] < 12 && numCells[x]>24)
	// 6 2 1
	#pragma omp parallel num_threads(12)
    {
        int numThreads = omp_get_num_threads();
        int tid = omp_get_thread_num(); // Thread ID
		// printf("\nthread = %d",tid);
		
		int threaddomain[3];
		threaddomain[x] = static_cast<int>(tid % 6);
		threaddomain[xmid] = static_cast<int>(tid/ 6);
		threaddomain[xmin] = 0;

		int chunkSize[3];

		chunkSize[x] = (numCells[x] + 5) / 6; 
		chunkSize[xmid] = (numCells[xmid] + 1 ) / 2;
		chunkSize[xmin] = numCells[xmin];

		// printf("thread: %d , domain: %d,%d,%d\n",tid, threaddomain[0],threaddomain[1],threaddomain[2]);

        int startIdxI = threaddomain[0] * chunkSize[0];
        int startIdxJ = threaddomain[1] * chunkSize[1];
        int startIdxK = threaddomain[2] * chunkSize[2];

        int endIdxI = std::min((threaddomain[0] + 1) * chunkSize[0], numCells[0]);
        int endIdxJ = std::min((threaddomain[1] + 1) * chunkSize[1], numCells[1]);
        int endIdxK = std::min((threaddomain[2] + 1) * chunkSize[2], numCells[2]);
		
        for (int i = startIdxI; i < endIdxI; i++) {
			
            const int startX = std::max(i - 1, 0);
            const int endX = std::min(i + 1, numCells[0] - 1);

            for (int j = startIdxJ; j < endIdxJ; j++) {

                    const int startY = std::max(j - 1, 0);
                    const int endY = std::min(j + 1, numCells[1] - 1);

                for (int k = startIdxK; k < endIdxK; k++) {
					
                    const int startZ = std::max(k - 1, 0);
                    const int endZ = std::min(k + 1, numCells[2] - 1);
					
                    auto& cell = cellBuffer[i][j][k];

                    for (int ii = startX; ii <= endX; ++ii) {
                        for (int jj = startY; jj <= endY; ++jj) {
                            for (int kk = startZ; kk <= endZ; ++kk) {
								//#pragma omp critical
                                {candidates_set[i][j][k].insert(cellBuffer[ii][jj][kk].begin(), cellBuffer[ii][jj][kk].end());}
                            }
                        }
                    }
					
					// cellBuffer[i][j][k].clear();
                }
            }
        }
    }
	else fillCandidates_seq();
}

void mesh::fillCandidates2() {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells[0]; i++) {
        for (int j = 0; j < numCells[1]; j++) {
            for (int k = 0; k < numCells[2]; k++) {
                const int startX = std::max(i - 1, 0);
                const int startY = std::max(j - 1, 0);
                const int startZ = std::max(k - 1, 0);
                const int endX = std::min(i + 1, numCells[0] - 1);
                const int endY = std::min(j + 1, numCells[1] - 1);
                const int endZ = std::min(k + 1, numCells[2] - 1);
				
                auto& cell = cellBuffer[i][j][k];
                for (int ii = startX; ii <= endX; ++ii) {
                    for (int jj = startY; jj <= endY; ++jj) {
                        for (int kk = startZ; kk <= endZ; ++kk) {
                            candidates_set[ii][jj][kk].insert(cell.begin(), cell.end());
                        }
                    }
                }
				
            	cell.clear();
            }
        }
    }
}

void mesh::fillCandidates3() {
	
    #pragma omp parallel
    {
		
		std::unordered_map<std::tuple<int, int, int>, std::set<int>, hash_tuple> localcandidatesSet_flat;


        #pragma omp for collapse(3)
        for (int i = 0; i < numCells[0]; i++) {
            for (int j = 0; j < numCells[1]; j++) {
                for (int k = 0; k < numCells[2]; k++) {
                    const int startX = std::max(i - 1, 0);
                    const int startY = std::max(j - 1, 0);
                    const int startZ = std::max(k - 1, 0);
                    const int endX = std::min(i + 1, numCells[0] - 1);
                    const int endY = std::min(j + 1, numCells[1] - 1);
                    const int endZ = std::min(k + 1, numCells[2] - 1);

                    auto& cell = cellBuffer[i][j][k];

                    for (int ii = startX; ii <= endX; ++ii) {
                        for (int jj = startY; jj <= endY; ++jj) {
                            for (int kk = startZ; kk <= endZ; ++kk) {
                                auto key = std::make_tuple(ii, jj, kk);
                                localcandidatesSet_flat[key].insert(cell.begin(), cell.end());
                            }
                        }
                    }
                }
            }
        }

        // Synchronize insertion into candidates_set
        
        {
            for (auto& entry : localcandidatesSet_flat) {
                auto& key = entry.first;
                auto& localSet = entry.second;
				 // Separate tuple key into components i, j, k
                int i = std::get<0>(key);
                int j = std::get<1>(key);
                int k = std::get<2>(key);

                // Access and insert into candidates_set
                candidates_set[i][j][k].insert(localSet.begin(), localSet.end());
            }
        }
    }
}

// void mesh::fillCandidates4(){
	
// 	#pragma omp parallel for collapse(3) 
// 		for(int i =0; i < numCells[0]; i++){
// 			for(int j = 0; j < numCells[1]; j++){
// 				for(int k = 0; k < numCells[2]; k++){
// 					for (int c = 0; c<27; c++){
// 						if(neighPos[c][0]+i < 0 or neighPos[c][0]+i >= numCells[0]) continue;
// 						if(neighPos[c][1]+j < 0 or neighPos[c][1]+j >= numCells[1]) continue;
// 						if(neighPos[c][2]+k < 0 or neighPos[c][2]+k >= numCells[2]) continue;
// 						candidates_set[i][j][k].insert(cellBuffer[i+neighPos[c][0]][j+neighPos[c][1]][k+neighPos[c][2]].begin(), cellBuffer[i+neighPos[c][0]][j+neighPos[c][1]][k+neighPos[c][2]].end());
						
// 					}
// 				}
// 			}
			
// 		}
	
// }
void mesh::fillCandidates4() {
	#pragma omp parallel for collapse(3)
	for (int i = 0; i < numCells[0]; ++i) {
		for (int j = 0; j < numCells[1]; ++j) {
			for (int k = 0; k < numCells[2]; ++k) {

				auto& cellCandidates = candidates_set[i][j][k];


				for (int c = 0; c < 27; ++c) {
					const int ni = i + neighPos[c][0];
					const int nj = j + neighPos[c][1];
					const int nk = k + neighPos[c][2];

					if ((unsigned)ni >= (unsigned)numCells[0] ||
						(unsigned)nj >= (unsigned)numCells[1] ||
						(unsigned)nk >= (unsigned)numCells[2])
						continue;

					const auto& neighborCell = cellBuffer[ni][nj][nk];

					if (!neighborCell.empty()) {
						cellCandidates.insert(neighborCell.begin(), neighborCell.end());
					}
					
				}
			}
		}
	}
}

void mesh::fillCandidates5(){
	
	#pragma omp parallel for collapse(3)
		for(int i =0; i < numCells[0]; i++){
			for(int j = 0; j < numCells[1]; j++){
				for(int k = 0; k < numCells[2]; k++){
                    auto& cell = cellBuffer[i][j][k];
						for (int c = 0; c<27; c++){
							if(neighPos[c][0]+i < 0 or neighPos[c][0]+i >= numCells[0]) continue;
							if(neighPos[c][1]+j < 0 or neighPos[c][1]+j >= numCells[1]) continue;
							if(neighPos[c][2]+k < 0 or neighPos[c][2]+k >= numCells[2]) continue;
							
							candidates_set[i+neighPos[c][0]][j+neighPos[c][1]][k+neighPos[c][2]].insert(cell.begin(), cell.end());
							
						}
					}
					// cellBuffer[i][j][k].clear();
			}
			
		}
	
}
void mesh::fillCandidates_seq(){
	
		for(int i =0; i < numCells[0]; i++){
			for(int j = 0; j < numCells[1]; j++){
				for(int k = 0; k < numCells[2]; k++){
					const int startX = std::max(i - 1, 0);
					const int startY = std::max(j - 1, 0);
					const int startZ = std::max(k - 1, 0);
					const int endX = std::min(i + 1, numCells[0] - 1);
					const int endY = std::min(j + 1, numCells[1] - 1);
					const int endZ = std::min(k + 1, numCells[2] - 1);
					for (int ii = startX; ii <= endX; ++ii) {
							for (int jj = startY; jj <= endY; ++jj) {
								for (int kk = startZ; kk <= endZ; ++kk) {
									const int bufferSize = cellBuffer[ii][jj][kk].size(); // Calculate size once
									for (int w = 0; w < bufferSize; ++w) {
										candidates_set[i][j][k].insert(cellBuffer[ii][jj][kk][w]);
									}
								}
							}
						}
					// cellBuffer[i][j][k].clear();
				}
			}
			
		}
	
}

// int mesh::calcDist(const float *p, const double range, const int cpx, const int cpy, const int cpz, float *closestPoint, double &cqd) {
//     if (candidates_set[cpx][cpy][cpz].empty()) {
//         return -1;
//     }

//     int cid = -1;
//     double threshold = 0.2 * partdist2;
//     float point[3] = {-1.0, -1.0, -1.0};

//     // Loop through candidates (use iterator since unordered_set doesn't support indexing)
//     int i = 0;
//     for (auto it = candidates_set[cpx][cpy][cpz].begin(); it != candidates_set[cpx][cpy][cpz].end(); ++it) {
//         int id = *it;

//         const float a[3] = { meshData.x0[id], meshData.y0[id], meshData.z0[id] };
//         const float b[3] = { meshData.x1[id], meshData.y1[id], meshData.z1[id] };
//         const float c[3] = { meshData.x2[id], meshData.y2[id], meshData.z2[id] };

//         // Compute the closest point (normal calculation)
//         double quadDist = calcNearestPoint(a, b, c, p, point);

//         if (quadDist < cqd) {
//             cqd = quadDist;
//             cid = id;
			
// 			closestPoint[0] = point[0];
// 			closestPoint[1] = point[1];
// 			closestPoint[2] = point[2]; 
//             // memcpy(closestPoint, point, 3 * sizeof(float));

//             if (cqd < threshold) return cid; // Early exit if the distance is below threshold
//         }
//         i++;
//     }

//     return cid;
// }


int mesh::calcDist(const float * __restrict p, const double range, const int cpx, const int cpy, const int cpz, float * __restrict closestPoint, double &cqd) {
    const auto& candidates = candidates_set[cpx][cpy][cpz];
    if (candidates.empty()) return -1;

    int cid = -1;
    const double threshold = 0.2 * partdist2;

    float tmpClosest[3]; 

    for (int id : candidates) {
		
        const float a[3] = { meshData.x0[id], meshData.y0[id], meshData.z0[id] };
        const float b[3] = { meshData.x1[id], meshData.y1[id], meshData.z1[id] };
        const float c[3] = { meshData.x2[id], meshData.y2[id], meshData.z2[id] };

        const double quadDist = calcNearestPoint(a, b, c, p, tmpClosest);

        if (quadDist < cqd) {
			
            cqd = quadDist;
            cid = id;

            closestPoint[0] = tmpClosest[0];
            closestPoint[1] = tmpClosest[1];
            closestPoint[2] = tmpClosest[2];

            if (quadDist < threshold) break; 
        }
    }

    return cid;
}

// std::array<int, 2> elemId2 = mesh.calc2Dist(testPoint, range, grid.getCellPosX(testPoint[0]), grid.getCellPosY(testPoint[1]), grid.getCellPosZ(testPoint[2]), closestPoint, sclosestPoint, closestPointDist2, sclosestPointDist2);
      
// std::array<int, 2> mesh::calc2Dist(const double *p, double range, int cpx, int cpy, int cpz, double *closestPoint, double *sclosestPoint, double &cqd, double &scqd){

// 	cqd = range;
// 	scqd = range/3;
// 	std::array<int, 2> cid = {-1,-1};
// 	double quadDist;
// 	double point[3] = {-1.0,-1.0,-1.0};
// 	//fprintf(stderr, "i = [%f], j = [%f], k = [%f]", p[X], p[Y], p[Z]);
// 	//fprintf(stderr, "Pi = [%d], Pj = [%d], Pk = [%d]", grid->getCellPosX(p[X]), grid->getCellPosY(p[Y]), grid->getCellPosZ(p[Z]));

// 	double a[3], b[3], c[3];
// 	long id;
	
// 	for (std::set<long>::iterator it=candidates_set[cpx][cpy][cpz].begin(); it!=candidates_set[cpx][cpy][cpz].end(); ++it){
// 		id = *it;
// 		a[0] = elem[id].node[0].pos->x;
// 		a[1] = elem[id].node[0].pos->y;
// 		a[2] = elem[id].node[0].pos->z;

// 		b[0] = elem[id].node[1].pos->x;
// 		b[1] = elem[id].node[1].pos->y;
// 		b[2] = elem[id].node[1].pos->z;

// 		c[0] = elem[id].node[2].pos->x;
// 		c[1] = elem[id].node[2].pos->y;
// 		c[2] = elem[id].node[2].pos->z;
// 		quadDist = calcNearestPoint(a,b,c,p,point);
		
// 		if(quadDist < cqd){
// 			cqd = quadDist;
// 			cid[0] = id;
// 			closestPoint[0] = point[0];
// 			closestPoint[1] = point[1];
// 			closestPoint[2] = point[2]; 
// 		}
// 		else if(quadDist < scqd){
// 			scqd = quadDist;
// 			cid[1] = id;
// 			sclosestPoint[0] = point[0];
// 			sclosestPoint[1] = point[1];
// 			sclosestPoint[2] = point[2]; 
// 		}

// 	}
// 	return  cid;

// }
int mesh::calcDist_flat(const float* __restrict p, const double range, const int fIndex, float* __restrict closestPoint, double& cqd) {
    const auto& candidates = candidatesSet_flat[fIndex];
    if (candidates.empty()) return -1;

    int cid = -1;
    const double threshold = 0.2 * partdist2;

    float tmpClosest[3]; // Scratch space

    for (const int id : candidates) {
        const float a[3] = { meshData.x0[id], meshData.y0[id], meshData.z0[id] };
        const float b[3] = { meshData.x1[id], meshData.y1[id], meshData.z1[id] };
        const float c[3] = { meshData.x2[id], meshData.y2[id], meshData.z2[id] };

        const double quadDist = calcNearestPoint(a, b, c, p, tmpClosest);

        if (quadDist < cqd) {
            cqd = quadDist;
            cid = id;

            closestPoint[0] = tmpClosest[0];
            closestPoint[1] = tmpClosest[1];
            closestPoint[2] = tmpClosest[2];

            if (cqd < threshold) return cid; // Early exit
        }
    }

    return cid;
}


/*
/// Find closest point on the mesh from a particle i
/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
void mesh::closestPointMesh(int nMesh, int *type, float *quaddev, float *pndS, float pndS0, float reL, int nPart, int *index, float *dMesh,
    int *idNearMesh, int *idNearElement, float *x, float *y, float *z, float *mx, float *my, float *mz, cellGrid *grid)
{
	float reL2 = reL*reL;
	//fprintf(stderr, "\n * RAIO LARG = %f \n",reL);
	struct timeval t1, t2;
	float beta = 0.90;
	float elapsedTime;
	int idElement[nPart];
	for(int i=0; i<nPart; i++)
	{
		idElement[i] = 0;
	}
	
	for(int i = 0; i<=2; i++) old[i] = -1;
	/// start timer
	gettimeofday(&t1, NULL);
	#pragma omp parallel for
	for(int i=0; i<nPart; i++)
	{
		/// type: flow=50, fluid=51, solid=53, solidMesh=61, boundCond=62
		/// boundcary condition: noPress=40, noFreeSurf=41, freeSurf=42,
		//fprintf(stderr, "\n * RAIO LARGE = %f \n",reL);
		int m = index[i];
		
		///fprintf(stderr,"%f", quaddev[i]);
		if(type[m]==51)
		{
			
			float testPoint[3] = {x[i], y[i], z[i]};
			int cpx, cpy, cpz;
			cpx = grid->getCellPosX(testPoint[0]);
			cpy = grid->getCellPosY(testPoint[1]);
			cpz = grid->getCellPosZ(testPoint[2]);

			//fprintf(stderr, "\n * x,y,z = %f - %f - %f \n",testPoint[0], testPoint[1], testPoint[2]);
			///Find the closest points to TestPoint
			float closestPoint[3] = {0.0,0.0,0.0};     ///the coordinates of the closest point will be returned here
			float closestPointDist2 = 0;   ///the squared distance to the closest point will be returned here
			int elemId;           ///the cell id of the cell containing the closest point will be returned here
			//fprintf(stderr, "\n * closestPointDist2 = %f \n",closestPointDist2);
			/// Get the coordinates of the closest point
			//fprintf(stderr, "\n * x,y,z = %f - %f - %f \n", closestPoint[0], closestPoint[1], closestPoint[2]);
			//if(pndS[i] < pndS0 * beta)
			elemId = calcDist(testPoint, reL, cpx, cpy, cpz, closestPoint, closestPointDist2);
			if (elemId != -1)
//			vtkStaticCellLocator::FindClosestPointWithinRadius
			{
				
				//fprintf(stderr, "\n * calcDistTrue ");
				//fprintf(stderr, "\n * closestPointDist2 = %f ", closestPointDist2);
				//fprintf(stderr, "\n * Part ID = %ld ", i);
				//fprintf(stderr, "\n * x,y,z = %f - %f - %f \n", closestPoint[0], closestPoint[1], closestPoint[2]);
				
				/// Mirror particle position Xm = Xi + 2*(Xw - Xi)
				mx[i] = x[i] + 2*(closestPoint[0] - x[i]);
				my[i] = y[i] + 2*(closestPoint[1] - y[i]);
				mz[i] = z[i] + 2*(closestPoint[2] - z[i]);
				idElement[i] = elemId;
			}
			else
			{
				
				//fprintf(stderr, "\n * calcDistFalse ");
				//fprintf(stderr, "\n * closestPointDist2 = %f ", closestPointDist2);
				//fprintf(stderr, "\n * x,y,z = %f - %f - %f \n", closestPoint[0], closestPoint[1], closestPoint[2]);
				/// Mirror particle position Xm -> 00 i.e. Xm =  Xi + 10*reL
				mx[i] = x[i] + 10*reL;
				my[i] = y[i] + 10*reL;
				mz[i] = z[i] + 10*reL;
				idElement[i] = -1;
			}
		}
	}

	/// stop timer
	gettimeofday(&t2, NULL);

	#pragma omp parallel
	{
		for(int i=0; i<nPart; i++)
		{
  			/// flow=50, fluid=51, solid=53, solidMesh=61, boundCond=62
  			int m = index[i];
			if(type[m]==51)
			{
				float p[3]; /// Position of particle i
				p[0] = x[i];
				p[1] = y[i];
				p[2] = z[i];
				float mp[3], pmp[3];
				mp[0] = mx[i];
				mp[1] = my[i];
				mp[2] = mz[i];
				subtraction(p,mp,pmp); /// pmp = p - mp
				/// Square of distance between particle i and wall particle
				float dpw2 = squaredNorm(pmp)/4.0;
				/// Calculate the distance between mesh and particle
				if (dpw2 < reL2)
				{
					float dpw = sqrt(dpw2);
					if(dMesh[i] == 0.0)
					{
						dMesh[i] = dpw;
						idNearElement[i] = idElement[i];
						idNearMesh[i] = nMesh;
					}
					else if(dMesh[i]>dpw)
					{
						dMesh[i] = dpw;
						idNearElement[i] = idElement[i];
						idNearMesh[i] = nMesh;
					}						
				}
			}
		}
	}
}
*/
	
double mesh::calcNearestPoint(const float* __restrict a, const float* __restrict b, const float* __restrict c, const float* __restrict p, float* __restrict nearest) {
	/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
	/// Returns the nearest point of the triangle a,b,c to particle p on pointer *nearest and squaredDistance(nearest, p)

		float abx = b[0] - a[0], aby = b[1] - a[1], abz = b[2] - a[2];
		float acx = c[0] - a[0], acy = c[1] - a[1], acz = c[2] - a[2];
		float apx = p[0] - a[0], apy = p[1] - a[1], apz = p[2] - a[2];
	
		float d1 = abx * apx + aby * apy + abz * apz;
		float d2 = acx * apx + acy * apy + acz * apz;
	
		// Region 1: Vertex A
		if (d1 <= 0.0f && d2 <= 0.0f) {
			nearest[0] = a[0]; nearest[1] = a[1]; nearest[2] = a[2];
			return squaredDist(nearest, p);
		}
	
		float bpx = p[0] - b[0], bpy = p[1] - b[1], bpz = p[2] - b[2];
		float d3 = abx * bpx + aby * bpy + abz * bpz;
		float d4 = acx * bpx + acy * bpy + acz * bpz;
	
		// Region 2: Vertex B
		if (d3 >= 0.0f && d4 <= d3) {
			nearest[0] = b[0]; nearest[1] = b[1]; nearest[2] = b[2];
			return squaredDist(nearest, p);
		}
	
		float vc = d1 * d4 - d3 * d2;
		// Region 3: Edge AB
		if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
			float v = d1 / (d1 - d3);
			nearest[0] = a[0] + v * abx;
			nearest[1] = a[1] + v * aby;
			nearest[2] = a[2] + v * abz;
			return squaredDist(nearest, p);
		}
	
		float cpx = p[0] - c[0], cpy = p[1] - c[1], cpz = p[2] - c[2];
		float d5 = abx * cpx + aby * cpy + abz * cpz;
		float d6 = acx * cpx + acy * cpy + acz * cpz;
	
		// Region 4: Vertex C
		if (d6 >= 0.0f && d5 <= d6) {
			nearest[0] = c[0]; nearest[1] = c[1]; nearest[2] = c[2];
			return squaredDist(nearest, p);
		}
	
		float vb = d5 * d2 - d1 * d6;
		// Region 5: Edge AC
		if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
			float w = d2 / (d2 - d6);
			nearest[0] = a[0] + w * acx;
			nearest[1] = a[1] + w * acy;
			nearest[2] = a[2] + w * acz;
			return squaredDist(nearest, p);
		}
	
		float va = d3 * d6 - d5 * d4;
		// Region 6: Edge BC
		if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
			float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
			nearest[0] = b[0] + w * (c[0] - b[0]);
			nearest[1] = b[1] + w * (c[1] - b[1]);
			nearest[2] = b[2] + w * (c[2] - b[2]);
			return squaredDist(nearest, p);
		}
	
		// Region 7: Inside face
		float denom = 1.0f / (va + vb + vc);
		float v = vb * denom;
		float w = vc * denom;
		nearest[0] = a[0] + abx * v + acx * w;
		nearest[1] = a[1] + aby * v + acy * w;
		nearest[2] = a[2] + abz * v + acz * w;
		return squaredDist(nearest, p);
	}
	

// double mesh::calcNearestPoint(const float  *a, const float  *b, const float  *c, const float *p, float *nearest){
// 		/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
// 		/// Returns the nearest point of the triangle a,b,c to particle p on pointer *nearest and squaredDistance(nearest, p)
		
		
// 		const float ab[3] = {b[0]-a[0],b[1]-a[1],b[2]-a[2]};
// 		const float ac[3] = {c[0]-a[0],c[1]-a[1],c[2]-a[2]};
// 		const float ap[3] = {p[0]-a[0],p[1]-a[1],p[2]-a[2]};
	
	
// 		const float  d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
// 		const float  d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];
	
	
// 		//First Region Check: nearest -> a vertex
// 		if (d1 <= 0.0 && d2 <= 0.0){
			
// 			for(int i = 0; i<=2; i++) nearest[i]=a[i];
// 			return squaredDist(nearest,p);
// 		}
	
// 		const float bp[3] = {p[0]-b[0],p[1]-b[1],p[2]-b[2]};
	
// 		const float  d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
// 		const float  d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];
	
	
// 		//Second Region Check: nearest -> b vertex
// 		if (d3 >= 0.0 && d4 <= d3){
			
// 			for(int i = 0; i<=2; i++) nearest[i]=b[i];
// 			return squaredDist(nearest,p);
// 		}
	
// 		const float  vc = d1*d4 - d3*d2;
// 		//Third Region Check: nearest -> point in ab edge
// 		if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
			
// 			const float  v = d1 / (d1 - d3);
// 			for(int i = 0; i<=2; i++) nearest[i]= a[i] + v * ab[i];
// 			return squaredDist(nearest,p);
// 	}
		
// 		const float cp[3] = {p[0]-c[0],p[1]-c[1],p[2]-c[2]};
// 		const float  d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
// 		const float  d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];	
	
	
// 		//Fourth Region Check: nearest -> c vertex
// 		if  (d6 >= 0.0 && d5 <= d6){
			
// 			for(int i = 0; i<=2; i++) nearest[i]=c[i];
// 			return squaredDist(nearest,p);
// 		} 
	
// 		const float  vb = d5*d2 - d1*d6;
		
// 		//Fifth Region Check: nearest -> point in ac edge 
// 		if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
			
// 			 const float  w = d2 / (d2 - d6);
// 			for( int  i = 0; i<=2; i++) nearest[i] = a[i] + w * ac[i];
// 			return squaredDist(nearest,p);
// 		}
	
// 		const float  va = d3*d6 - d5*d4;
		
// 		//Sixth Region Check: nearest -> point in bc edge 
// 		if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
			
// 			const float  w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
// 			for( int  i = 0; i<=2; i++) nearest[i] = b[i] + w * (c[i] - b[i]); 
// 			return squaredDist(nearest,p);
// 		}
// 		 const float  denom = 1.0 / (va + vb + vc);
// 		 const float  v = vb * denom;
// 		 const float  w = vc * denom;
		
// 		//Seventh Region Check: nearest -> point inside triangle 
		
// 		for( int  i = 0; i<=2; i++) nearest[i] = a[i] + ab[i] * v + ac[i] * w;
// 		return squaredDist(nearest,p);
	
// 	}



void mesh::writeCellBufferToFile(const std::string& filename) {
    std::ofstream outFile(filename);

    // Verifica se o arquivo foi aberto corretamente
    if (!outFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    // Percorre o cellBuffer
    for (int x = 0; x < cellBuffer.size(); ++x) {
        for (int y = 0; y < cellBuffer[x].size(); ++y) {
            for (int z = 0; z < cellBuffer[x][y].size(); ++z) {
				for (const auto& it : cellBuffer[x][y][z]) 
				outFile << " " << x << " " << y << " " << z << " " << it << std::endl;
                // for (int w = 0; w < cellBuffer[x][y][z].size(); ++w) {
                //     // Escreve o ID no arquivo
                //     // outFile << "Celula (" << x << "," << y << "," << z << "): ID = " << cellBuffer[x][y][z][w] << std::endl;
                //     outFile << " " << x << " " << y << " " << z << " " << cellBuffer[x][y][z][w] << std::endl;
                // }
            }
        }
    }

    // Fecha o arquivo
    outFile.close();
}