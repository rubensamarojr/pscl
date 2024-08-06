 /** MPS - Moving Particle Semi-
 ***(c) 2019 USP/TPN
 *** Rubens Amaro, Lucas Pereira
 **/

#include "polygon_wall.h"
//#include "mps.h"
#include <iostream>
#include <cmath>
#include <omp.h>

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




void cellGrid::initCell(double dCx, double lSxMin, double lSyMin ,double lSzMin, double lSxMax, double lSyMax ,double lSzMax){

	///Need to be called from init.cpp with the mesh
    ///Memory allocation


	///Assigning values to the cell class atributes	
	///// Same as the parameters from MPS particles
	
    eCells[0] = dCx;		//dCell[0]
    eCells[1] = dCx;		//dCell[1]
    eCells[2] = dCx;		//dCell[2]

	//fprintf(stderr, "\n* QUANTIDADE DE CELULAS [x][y][z]:::[%d][%d][%d] ", nCx,nCy,nCz);
	//fprintf(stderr, "\n* TAMANHO DAS CELULAS [x][y][z]:::[%f][%f][%f] ", dCx,dCy,dCz);

	int nCx = int((lSxMax - lSxMin)/dCx) +1;
	int nCy = int((lSyMax - lSyMin)/dCx) +1;
	int nCz = int((lSzMax - lSzMin)/dCx) +1;


    numCells[0] = nCx; 	//nCell[0]
    numCells[1] = nCy; 	//nCell[1]
    numCells[2] = nCz;	//nCell[2]
	//printf("\nnumcells = %d,%d,%d",numCells[0],numCells[1],numCells[2]);

	/// Simulation boundaries(max and min coordinates on each direction)

	minSpc[0] = lSxMin;				//limSpace[0][1]
	minSpc[1] = lSyMin;				//limSpace[1][1]
	minSpc[2] = lSzMin;				//limSpace[2][1]

	maxSpc[0] = lSxMax;				//limSpace[0][1]
	maxSpc[1] = lSyMax;				//limSpace[1][1]
	maxSpc[2] = lSzMax;				//limSpace[2][1]

	cellBuffer.resize(numCells[0]);
	for(int i=0; i<numCells[0]; i++)
	{
		cellBuffer[i].resize(numCells[1]);
		for(int j=0; j<numCells[1]; j++){
			cellBuffer[i][j].resize(numCells[2]);
		}

	}
	
}

void cellGrid::getBounds(int a[][2], double *b, double *c, double *d){
	// double e = 0.55*eCells[0];
	double e = 0.0*eCells[0];
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
	
	aux = int(std::min(z0, z1));
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

// CALCULA INTERSECÇÃO - 0.5 -> RETORNA PONTOS / USADO EM 1
Point2 cellGrid::intersectionPointZ(const node& p1, const Point& p2, const double p) {
    Point edgeVector = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    double t = -((p1.z - p) / edgeVector.z);
    return {p1.x + t * edgeVector.x, p1.y + t * edgeVector.y};
}
void cellGrid::fillBuffer(mesh *mesh){
	
	std::vector<std::vector<long> > cellcd;	
    cellcd.resize(numCells[0]);
    std::vector<Point2> intersections, intersectold;
    intersections.reserve(2);
    intersectold.reserve(2);

	
	int a, b, c;
	int n, t, i, j, k;
	for(i=0; i<(numCells[0]); i++)
		for(j=0; j<(numCells[1]); j++)
			for(k=0;(k< numCells[2]); k++)
				cellBuffer[i][j][k].clear();
			
	nEl = mesh->getnElements(); 
	double cen[3];
	double normalx, normaly, normalz;
	double vx[3], vy[3], vz[3];
	int bound[3][2] = {0,0,0,0,0,0};
	for (t=0; t<nEl; t++){
			
		for(n = 0; n<3; n++){
			vx[n] = mesh->getNodePositionX(t,n);
			vy[n] = mesh->getNodePositionY(t,n);
			vz[n] = mesh->getNodePositionZ(t,n);
		}
	
		///store the cells of the vertices as cells containing the triangle defined by its vertices
		//for(n= 0; n < 3; n++){
		//	a = getCellPosX(vx[n]);
		//	b = getCellPosY(vy[n]);
		//	c = getCellPosZ(vz[n]);
		//	cellBuffer[a][b][c].push_back(t);			
		//	}
		///define the domain of cells to test the triangle

		
		getBounds(bound, vx, vy, vz);
		normalx = mesh->getNormalX(t);
		normaly = mesh->getNormalY(t);
		normalz = mesh->getNormalZ(t);
		
		for (i = bound[0][0]; i<=bound[0][1]; i++){
			for (j = bound[1][0]; j<=bound[1][1]; j++){
				for (k = bound[2][0]; k<=bound[2][1]; k++){	
					if((i<0) || (j<0) || (k<0))
						continue;
					if((i>=numCells[0]) || (j>=numCells[1]) || (k>=numCells[2]))
						continue;
					//cellBuffer[i][j][k].size() == 0 means theres no element in it yet
					//if(cellBuffer[i][j][k].size()==0 || cellBuffer[i][j][k].back()!=t){
						
						cen[0] = getCellCenterX(i);
						cen[1] = getCellCenterY(j);
						cen[2] = getCellCenterZ(k);
						
						if(triangleAABB_intersection(cen, normalx, normaly, normalz, vx[0], vy[0], vz[0], vx[1], vy[1], vz[1], vx[2], vy[2], vz[2])){
						
							cellBuffer[i][j][k].push_back(t);
						
						}	
					//}	
					//cellBuffer[i][j][k].back() = means triangle vertex in cell->avoid double test
					/*else if (){						
						cen[0] = getCellCenterX(i);
						cen[1] = getCellCenterY(j);
						cen[2] = getCellCenterZ(k);
						if(triangleAABB_intersection(cen, normalx, normaly, normalz, vx[0], vy[0], vz[0], vx[1], vy[1], vz[1], vx[2], vy[2], vz[2])){
							
							cellBuffer[i][j][k].push_back(t);
						}						
					}
					*/
				}
			}
		}
	}	
	
/*	
	int ncellx = cellBuffer.size();
	int ncelly = cellBuffer[0].size();
	int ncellz = cellBuffer[0][0].size();
	
	fprintf(stderr, "\n*cell nx - %d", ncellx);
	fprintf(stderr, "\n*cell ny - %d", ncelly);
	fprintf(stderr, "\n*cell nz - %d", ncellz);
	
	for (i = 0; i < ncellx*0.5 ; i++){
			for (j = 0 ; j < 2; j++){
				for (k = 0; k < ncellz; k++){	
					long nelem = cellBuffer[i][j][k].size();
					if (nelem == 0) fprintf(stderr, "\n*cellBuffer[%d][%d][%d]-> [/] ",i,j,k);
					else{
						fprintf(stderr, "\n*cellBuffer[%d][%d][%d]",i,j,k);
						for(int l = 0; l<nelem;l++){
							fprintf(stderr, "\n*cellBuffer[%d][%d][%d][%ld] = %ld",i,j,k,l,cellBuffer[i][j][k][l]);
						}
					}
				}
			}
		}
	
	fprintf(stderr, "\n*cell nx - %d", ncellx);
	fprintf(stderr, "\n*cell ny - %d", ncelly);
	fprintf(stderr, "\n*cell nz - %d", ncellz);
*/
	
	


}

bool cellGrid::triangleAABB_intersection(double *cen, double normX, double normY, double normZ, double v0x, double v0y, double v0z, 
double v1x, double v1y, double v1z, double v2x, double v2y, double v2z){
	/*
	/ ref: https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox_tam.pdf
	/ cen[2] = center of aabb coordinates
	/ v[3][3] = vertices of mesh element(triangle)
	/ cell->ecell->spc->i = extent of aabb on i direction; i=0,1,2
	*/

	int X,Y,Z;
	double v0[3] = {0.0,0.0,0.0};
	double v1[3] = {0.0,0.0,0.0};
	double v2[3] = {0.0,0.0,0.0};
	
	double normalAux[3] = {normX,normY,normZ};
	double v[3][3]={v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z};

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
	double p[2], r;

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
	double plane_dist = normalAux[X]*v0x + normalAux[Y]*v0y + normalAux[Z]*v0z;

	double s = normalAux[X]*cen[X] + normalAux[Y]*cen[Y] + normalAux[Z]*cen[Z] - plane_dist;
	r = eCells[0] * fabs(normalAux[X]) + eCells[1] * fabs(normalAux[Y]) + eCells[2] * fabs(normalAux[Z]);
	return (fabs(s)<=r);
}


/* UNUTILIZED
/// Return the cell position on the cellIndex list from the coordinates given
long cellGrid::getPositionByCoords(double cx, double cy, double cz){

	///cx = coordinate in x direction
	///cy = coordinate in y direction
	///cz = coordinate in z direction

	long position;
	double dx, dy, dz;

	dx = cx - minSpc[0];
	dy = cy - minSpc[1];
	dz = cz - minSpc[2];

	position = long(dz/eCells[2])*getCellsPerDepth(); //adding the cells per unit of depth times depth(Z)
	position += long(dy/eCells[1])*numCells[0];
	position += long(dx/eCells[0]);

	return position;
}

/// Update the coords array with the coordinates of the cell center from the cell position on the cellIndex list
void cellGrid::getCoordByPosition(long pos, double *coords){
	int posY, posZ;

	posZ = int(pos / getCellsPerDepth()); 		//Z coordinate considering the Z extent as the unit
	pos = pos % getCellsPerDepth();

	posY = int(pos / numCells[0]);		//Y coordinate considering the Y extent as the unit
	pos = pos % numCells[0];		//X coordinate considering the X extent as the unit



	//the (position on the list in each direction plus half) times the extent of the cell resulting in its center
	coords[0] = (pos+0.5)*eCells[0] + minSpc[0];
	coords[1] = (posY+0.5)*eCells[1] + minSpc[1];
	coords[2] = (posZ+0.5)*eCells[2] + minSpc[2];

}
*/

/// Calculate the Cell position on direction by the given direction coordinate

int cellGrid::getCellPosX(double cx){
	
	//fprintf(stderr, "\n*cell coordinate x - %f", cx);
	//fprintf(stderr, "\n*cell position x - %d", int((cx - minSpc[0])/eCells[0]));
	return int((cx - minSpc[0])/eCells[0]);
}
int cellGrid::getCellPosY(double cy){
	
	//fprintf(stderr, "\n*cell coordinate y - %f", cy);
	//fprintf(stderr, "\n*cell position y - %d", int((cy - minSpc[1])/eCells[1]));
	return int((cy - minSpc[1])/eCells[1]);
}
int cellGrid::getCellPosZ(double cz){
	
	return int((cz - minSpc[2])/eCells[2]);
}


/// Calculate the Cell center coordinate on direction by the given direction position
double cellGrid::getCellCenterX(double cx){

	return (cx+0.5)*eCells[0] + minSpc[0];
}
double cellGrid::getCellCenterY(double cy){

	return (cy+0.5)*eCells[1] + minSpc[1];
}
double cellGrid::getCellCenterZ(double cz){

	return (cz+0.5)*eCells[2] + minSpc[2];
}



int cellGrid::getNumCells(int dir){
	return numCells[dir];
}

/// Return the (d-th - 1) element id on the cell [a][b][c] if it exists, return -1 otherwise
long cellGrid::getCellBuffer(int a, int b, int c, long d){
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
void mesh::subtraction(double *a, double *b, double *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

double mesh::squaredDist(double *a, double *b)
{
	return (a[0]-b[0]) * (a[0]-b[0]) + (a[1]-b[1]) * (a[1]-b[1]) + (a[2]-b[2]) * (a[2]-b[2]);
}

long mesh::getnElements(void)
{
   return nElements;
}

long mesh::getnNodes(long ielem)
{
   return elem[ielem].nNodes;
}

double mesh::getNodePositionX(long ielem, long inode)
{
   return elem[ielem].node[inode].pos->x;
}

double mesh::getNodePositionY(long ielem, long inode)
{
   return elem[ielem].node[inode].pos->y;
}

double mesh::getNodePositionZ(long ielem, long inode)
{
   return elem[ielem].node[inode].pos->z;
}

double mesh::getNormalX(long ielem)
{
   return elem[ielem].normal.pos->x;
}
double mesh::getNormalY(long ielem)
{
   return elem[ielem].normal.pos->y;
}
double mesh::getNormalZ(long ielem)
{
   return elem[ielem].normal.pos->z;
}

/// https://github.com/sreiter/stl_reader
void mesh::readMeshFile(const char * meshfilename)
{	
	std::vector<double> coords, normal_vecs;
	std::vector<unsigned int> elems, solids;
	try {

		// Load a mesh
		bool success = stl_reader::ReadStlFile(meshfilename, coords, normal_vecs, elems, solids);
		
		nElements = elems.size() / 3;
		// printf("\nnElements = %d\n", nElements );
		elem = (struct element*) malloc(nElements*sizeof(struct element));
		
		for(size_t ielem = 0; ielem < nElements; ++ielem) {
			//std::cout << "\ncoordinates of triangle " << ielem << ": ";
			elem[ielem].nNodes = 3; /// number of nodes per element
			elem[ielem].node = (struct nodes*) malloc(3*sizeof(struct nodes));
			for(int inode = 0; inode < 3; ++inode) {
				double* c = &coords[3 * elems [3 * ielem + inode]];
          		//std::cout << "(" << c[0] << ", " << c[1] << ", " << c[2] << ") \n";
				elem[ielem].node[inode].pos = (struct position*) malloc(1*sizeof(struct position));
				elem[ielem].node[inode].pos->x = c[0];
				elem[ielem].node[inode].pos->y = c[1];
				elem[ielem].node[inode].pos->z = c[2];
          		//std::cout << "(" << elem[ielem].node[inode].pos->x << ", " << elem[ielem].node[inode].pos->y << ", " << elem[ielem].node[inode].pos->z << ") \n";
			}
			
			//std::cout << std::endl;
		
			double* n = &normal_vecs[3 * ielem];
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


	/// Print to check
 	/*for(size_t ielem = 0; ielem < nElements; ++ielem)
 	{
 		std::cerr << "coordinates of element " << ielem << std::endl;
		for(size_t inode = 0; inode < elem[ielem].nNodes; ++inode)
		{
			std::cerr << "" << elem[ielem].node[inode].pos->x << ", " << elem[ielem].node[inode].pos->y << ", " << elem[ielem].node[inode].pos->z << " ; ";
		}
		std::cerr	<< "normal of element " << ielem << ": "
 		  				<< "" << elem[ielem].normal.pos->x << ", " << elem[ielem].normal.pos->y << ", " << elem[ielem].normal.pos->z << std::endl;
        std::cerr << "transformation matrix of element " << ielem << std::endl;
        std::cerr << elem[ielem].Rref[0] << ", " << elem[ielem].Rref[1] << ", " << elem[ielem].Rref[2] << ", "
                  << elem[ielem].Rref[3] << ", " << elem[ielem].Rref[4] << ", " << elem[ielem].Rref[5] << ", "
                  << elem[ielem].Rref[6] << ", " << elem[ielem].Rref[7] << ", " << elem[ielem].Rref[8] << ", " << std::endl;
	}*/
}


void mesh::fillCandidates(cellGrid *grid){

	int nc0 = grid->getNumCells(0);
	int nc1 = grid->getNumCells(1);
	int nc2 = grid->getNumCells(2);

	candidates_set.resize(nc0);
		for(int i=0; i<nc0; i++)
		{
			candidates_set[i].resize(nc1);
			for(int j=0; j<nc1; j++){
				candidates_set[i][j].resize(nc2);
			}

		}
	#pragma omp parallel for 
	for(int i =0; i < nc0; i++){
		for(int j = 0; j < nc1; j++){
			for(int k = 0; k < nc2; k++){
				findCandidates(i,j,k,nc0,nc1,nc2,grid);
				//if (candidates_set[i][j][k].size() == 0){
				//	 findCandidates(i,j,k,nc0,nc1,nc2,grid,2);
				//}
			}
		}
		
	}
}
void mesh::findCandidates(int cx, int cy, int cz, int ncx, int ncy, int ncz, cellGrid *grid){
	long id;
	//if (cx==0 && cy == 15 && cz == 3) fprintf(stderr, "\n * ID = %ld \n",id);
	for(int i=cx-1 ; i<= cx+1; i++){
		for(int j =cy-1 ; j<= cy+1; j++){
			for(int k =cz-1 ; k<= cz+1; k++){

				if((i<0) || (j<0) || (k<0))
					continue;
				if((i>=ncx) || (j>=ncy) || (k>=ncz))
					continue;
				
				for (long w = 0; w < nElements; w++){
					id = grid->getCellBuffer(i,j,k,w);
					//if (cx==0 && cy == 15 && cz == 3) fprintf(stderr, "\n * cell[%ld][%ld][%ld] \n",i,j,k);
					//if (cx==0 && cy == 15 && cz == 3) fprintf(stderr, "\n * w = %ld \n",w);
					//if (cx==0 && cy == 15 && cz == 3) fprintf(stderr, "\n * id = %ld \n",id);
					if (id == -1) break;
					candidates_set[cx][cy][cz].insert(id);
					//fprintf(stderr, "\n * ID = %ld \n",id);
				}
				
			}
		}
	}
			
	//fprintf(stderr, "\n * Cell[%ld][%ld][%ld] \n",cx, cy, cz);
	//fprintf(stderr, "\n * Cell[%ld][%ld][%ld] \n",old[0], old[1], old[2]);
	//if(candidates_set.empty()) fprintf(stderr, "\n * Cell[%ld][%ld][%ld]  VAZIA!!!!!!!!!!!!!!!!!!!! \n",cx, cy, cz);
	/*
	if(cy <=1){
		
		fprintf(stderr, "\n * Cell[%ld][%ld][%ld] \n",cx, cy, cz);
		for (std::set<long>::iterator it=candidates_set.begin(); it!=candidates_set.end(); ++it){
			fprintf(stderr, "\n * ID = %ld \n",*it);
		}
	}*/

}
long mesh::calcDist(double *p, double range, int cpx, int cpy, int cpz, double *closestPoint, double &cqd){
	cqd = range*range;
	long cid = -1;
	double quadDist;
	double point[3] = {-1.0,-1.0,-1.0};
	//fprintf(stderr, "i = [%f], j = [%f], k = [%f]", p[X], p[Y], p[Z]);
	//fprintf(stderr, "Pi = [%d], Pj = [%d], Pk = [%d]", grid->getCellPosX(p[X]), grid->getCellPosY(p[Y]), grid->getCellPosZ(p[Z]));
	if(candidates_set[cpx][cpy][cpz].empty()){
		
		return cid;
	}
	double a[3], b[3], c[3];
	long id;
	
	for (std::set<long>::iterator it=candidates_set[cpx][cpy][cpz].begin(); it!=candidates_set[cpx][cpy][cpz].end(); ++it){
		id = *it;
		a[0] = elem[id].node[0].pos->x;
		a[1] = elem[id].node[0].pos->y;
		a[2] = elem[id].node[0].pos->z;

		b[0] = elem[id].node[1].pos->x;
		b[1] = elem[id].node[1].pos->y;
		b[2] = elem[id].node[1].pos->z;

		c[0] = elem[id].node[2].pos->x;
		c[1] = elem[id].node[2].pos->y;
		c[2] = elem[id].node[2].pos->z;
		quadDist = calcNearestPoint(a,b,c,p,point);
		
		if(quadDist < cqd){
			cqd = quadDist;
			cid = id;
			for(int i =0; i<=2; i++) closestPoint[i] = point[i]; 
		}

	}
	return cid;

}

/*
/// Find closest point on the mesh from a particle i
/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
void mesh::closestPointMesh(int nMesh, int *type, double *quaddev, double *pndS, double pndS0, double reL, int nPart, long *index, double *dMesh,
    long *idNearMesh, long *idNearElement, double *x, double *y, double *z, double *mx, double *my, double *mz, cellGrid *grid)
{
	double reL2 = reL*reL;
	//fprintf(stderr, "\n * RAIO LARG = %f \n",reL);
	struct timeval t1, t2;
	double beta = 0.90;
	double elapsedTime;
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
			
			double testPoint[3] = {x[i], y[i], z[i]};
			int cpx, cpy, cpz;
			cpx = grid->getCellPosX(testPoint[0]);
			cpy = grid->getCellPosY(testPoint[1]);
			cpz = grid->getCellPosZ(testPoint[2]);

			//fprintf(stderr, "\n * x,y,z = %f - %f - %f \n",testPoint[0], testPoint[1], testPoint[2]);
			///Find the closest points to TestPoint
			double closestPoint[3] = {0.0,0.0,0.0};     ///the coordinates of the closest point will be returned here
			double closestPointDist2 = 0;   ///the squared distance to the closest point will be returned here
			long elemId;           ///the cell id of the cell containing the closest point will be returned here
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
				double p[3]; /// Position of particle i
				p[0] = x[i];
				p[1] = y[i];
				p[2] = z[i];
				double mp[3], pmp[3];
				mp[0] = mx[i];
				mp[1] = my[i];
				mp[2] = mz[i];
				subtraction(p,mp,pmp); /// pmp = p - mp
				/// Square of distance between particle i and wall particle
				double dpw2 = squaredNorm(pmp)/4.0;
				/// Calculate the distance between mesh and particle
				if (dpw2 < reL2)
				{
					double dpw = sqrt(dpw2);
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

double mesh::calcNearestPoint(double *a, double *b, double *c, double *p, double *nearest){
	/// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
	/// Returns the nearest point of the triangle a,b,c to particle p on pointer *nearest and squaredDistance(nearest, p)
	
	int X,Y,Z;
	X = 0; Y = 1; Z = 2;
	double ab[3], ac[3], ap[3];
	double d1, d2;
	subtraction(b,a,ab);
	subtraction(c,a,ac);
	subtraction(p,a,ap);

	/*
	int flag = 0;
	if(p[X]> 0.009 && p[X]<0.011 && p[Y] < -0.29 && p[Y]>-0.291 && p[Z] == 0.14) flag = 1;
	if(flag == 1){
	fprintf(stderr, "\n *a = [%f][%f][%f]", a[X], a[Y], a[Z]);
	fprintf(stderr, "\n *b = [%f][%f][%f]", b[X], b[Y], b[Z]);
	fprintf(stderr, "\n *c = [%f][%f][%f]", c[X], c[Y], c[Z]);
	fprintf(stderr, "\n *p = [%f][%f][%f]", p[X], p[Y], p[Z]);
	fprintf(stderr, "\n *ab = [%f][%f][%f]", ab[X], ab[Y], ab[Z]);
	fprintf(stderr, "\n *ac = [%f][%f][%f]", ac[X], ac[Y], ac[Z]);
	fprintf(stderr, "\n *ap = [%f][%f][%f]", ap[X], ap[Y], ap[Z]);
	}
	*/

	d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
	d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];

	/*
	if(flag == 1){
		fprintf(stderr, "\n *d1 = [%f]", d1);
		fprintf(stderr, "\n *d2 = [%f]", d2);
	}
	*/

	//First Region Check: nearest -> a vertex
	if (d1 <= 0.0 && d2 <= 0.0){
		//if (flag == 1) fprintf(stderr, "\n *1 - a vertex");
		for(int i = 0; i<=2; i++) nearest[i]=a[i];
		return squaredDist(nearest,p);
	} 

	double bp[3];
	double d3, d4;
	subtraction(p,b,bp);

	d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
	d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];
	/*
	if(flag == 1){
		fprintf(stderr, "\n *d3 = [%f]", d3);
		fprintf(stderr, "\n *d4 = [%f]", d4);
	}
	*/

	//Second Region Check: nearest -> b vertex
	if (d3 >= 0.0 && d4 <= d3){
		//if (flag == 1) fprintf(stderr, "\n *2 - b vertex");
		for(int i = 0; i<=2; i++) nearest[i]=b[i];
		return squaredDist(nearest,p);
	} 

	double vc = d1*d4 - d3*d2;
	//Third Region Check: nearest -> point in ab edge
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
		//if (flag == 1) fprintf(stderr, "\n *3 - ab edge");
		double v = d1 / (d1 - d3);
		for(int i = 0; i<=2; i++) nearest[i]= a[i] + v * ab[i];
		return squaredDist(nearest,p);
	}

	double cp[3];
	double d5, d6;
	subtraction(p,c,cp);
	d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
	d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];	

	/*if(flag == 1){
		fprintf(stderr, "\n *d5 = [%f]", d5);
		fprintf(stderr, "\n *d6 = [%f]", d6);
	}*/

	//Fourth Region Check: nearest -> c vertex
	if  (d6 >= 0.0 && d5 <= d6){
		//if (flag == 1) fprintf(stderr, "\n *4 - c vertex");
		for(int i = 0; i<=2; i++) nearest[i]=c[i];
		return squaredDist(nearest,p);
	} 

	double vb = d5*d2 - d1*d6;
	
	//Fifth Region Check: nearest -> point in ac edge 
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
		//if (flag == 1) fprintf(stderr, "\n *5 - ac edge");
		double w = d2 / (d2 - d6);
		for(int i = 0; i<=2; i++) nearest[i] = a[i] + w * ac[i];
		return squaredDist(nearest,p);
	}

	double va = d3*d6 - d5*d4;
	
	//Sixth Region Check: nearest -> point in bc edge 
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
		//if (flag == 1) fprintf(stderr, "\n *6 - bc edge");
		double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		for(int i = 0; i<=2; i++) nearest[i] = b[i] + w * (c[i] - b[i]); 
		return squaredDist(nearest,p);
	}
	double denom = 1.0 / (va + vb + vc);
	double v = vb * denom;
	double w = vc * denom;
	
	//Seventh Region Check: nearest -> point inside triangle 
	//if (flag == 1) fprintf(stderr, "\n *7 - triangle");
	for(int i = 0; i<=2; i++) nearest[i] = a[i] + ab[i] * v + ac[i] * w;
	return squaredDist(nearest,p);

}