#include <fcpw/fcpw.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeDMAT.h>
#include <igl/readDMAT.h>
#include <igl/parallel_for.h>
#include <igl/get_seconds.h>
#include <igl/AABB.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <cstdio>
#include "polygon_wall.h"
#include <fstream>
#include <iostream>
#include <string>

// COMMENT LINE TO DISABLE TEST
//#define LIBIGL_TEST
// #define FCPW_TEST
#define FPDC_TEST
// #define BRUTE_TEST

//#define WRITE_DMAT_FILES

int main(int argc, char * argv[])
{
	const auto & tictoc = []()
	{
		static double t_start = igl::get_seconds();
		double diff = igl::get_seconds()-t_start;
		t_start += diff;
		return diff;
	};


	///////////////// SET VALUES /////////////////
	// Number of tests
	int numTests = 3;
	
	// Declare the variables
	int numPoints, division;

	// Input the integer
	std::cout << "Enter the number of points to be generated: aa";
	std::cin >> numPoints;

	std::cout << "Enter an integer [20-60] to divide the minimum mesh dimension (used to define particle distance): ";
	std::cin >> division;

	std::string pathfilename = argv[1]; // the path and name of the file
	std::string outfilename = pathfilename.substr(5, pathfilename.length() - 9); // the file name
	outfilename += "_nPoints" + std::to_string(numPoints) + "_div" + std::to_string(division) + ".txt"; // adds ".txt" to the end of the file name

	std::cout << "Output file name: " << outfilename << '\n';
	std::ofstream outfile;
	outfile.open(outfilename); // opens the file with the given name

	typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> MatrixdX3R;
	MatrixdX3R V;
	// Use RowMajor so we can use direclty in setObjectTriangle
	Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> F;

	igl::read_triangle_mesh(argv[1],V,F);

	const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
	const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();

	double partdist = std::min(Vmax(0)-Vmin(0), Vmax(1)-Vmin(1));
	partdist = std::min(partdist, Vmax(2)-Vmin(2));

	partdist /= division; // Minimum dimension divided by 20
	double reL = 2.1*partdist;  // 3D MPS effective radius
	double reL2 = reL*reL;
	double dcell = 1.02*reL;
	printf("Particle distance: %.3e\n", partdist);
	outfile << "Number of Particles: " << numPoints << std::endl;
	outfile << "Particle distance: " << partdist << " m" << std::endl;
	outfile << "Effective radius: " << reL << " m" << std::endl;
	outfile << "Min mesh dimensions: " << Vmin(0) << " m " << Vmin(1) << " m " << Vmin(2) << " m " << std::endl;
	outfile << "Max mesh dimensions: " << Vmax(0) << " m " << Vmax(1) << " m " << Vmax(2) << " m " << std::endl;

	// Generate a list of random query points in the bounding box
	MatrixdX3R Q;
	if(argc<=2 || !igl::readDMAT(argv[2],Q))
	{
		printf("generating random points\n");
		Q = MatrixdX3R::Random(numPoints,3);
		// const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
		// const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
		const Eigen::RowVector3d Vdiag = Vmax-Vmin;
		for(int q = 0;q<Q.rows();q++)
		{
			Q.row(q) = (Q.row(q).array()*0.5+0.5)*Vdiag.array() + Vmin.array();
		}
	}

	///////////////// SET VALUES END /////////////////


	///////////////// LIBIGL /////////////////
#ifdef LIBIGL_TEST

	Eigen::MatrixXi SVI, SVJ, NF;
	Eigen::MatrixXd NV; 
	igl::remove_duplicate_vertices(V, F, 2.0e-15, NV, SVI, SVJ, NF);
	igl::AABB<Eigen::MatrixXd,3> treeMesh;
	
	// Build the acceleration structure
	tictoc();
	treeMesh.init(NV, NF);
	double tPreLibigl = tictoc();
	printf("Pre process LIBIGL %g \n",tPreLibigl);
	outfile << std::endl << "Pre process LIBIGL: " << tPreLibigl << std::endl;
	
	printf("Starting query LIBIGL\n");
	///QUERY CLOSESTPOINTMESH

	Eigen::MatrixXd closestPoint = Eigen::MatrixXd::Zero(Q.rows(), 3);
	Eigen::VectorXd closestPointDist2 = Eigen::VectorXd::Zero(Q.rows());
	Eigen::VectorXi elemId_igl = Eigen::VectorXi::Zero(Q.rows());
	Eigen::VectorXi I_igl(Q.rows());
	Eigen::VectorXd CQD_igl(Q.rows());
	
	double tQueryLibiglMean = 0.0;

	for (int ii = 0; ii < numTests; ++ii)
	{
		tictoc();
		treeMesh.squared_distance(NV, NF, Q, closestPointDist2, elemId_igl, closestPoint);
		double tQueryLibigl = tictoc();
		tQueryLibiglMean += tQueryLibigl;
		printf("Query LIBIGL (%d): %g\n",ii, tQueryLibigl);
		printf("Pre process + Query LIBIGL %g \n",tPreLibigl+tQueryLibigl);
		outfile << "Query LIBIGL (" << ii << "): " << tQueryLibigl << " s" << std::endl;
		outfile << "Pre process + Query LIBIGL: " << tPreLibigl+tQueryLibigl << " s" << std::endl;
	}
	
	tQueryLibiglMean /= numTests;
	printf("Pre process + Query LIBIGL Mean %g \n\n",tPreLibigl+tQueryLibiglMean);
	outfile << "Pre process + Query LIBIGL Mean: " << tPreLibigl+tQueryLibiglMean << " s" << std::endl;

	//#pragma omp parallel for
	for(int q = 0;q<Q.rows();q++)
	{
		I_igl(q) = elemId_igl[q];
		CQD_igl(q) = closestPointDist2[q];
	}
	
#ifdef WRITE_DMAT_FILES
	igl::writeDMAT("Q_igl.dmat",Q);
	igl::writeDMAT("CQD_igl.dmat",CQD_igl);
	igl::writeDMAT("I_igl.dmat",I_igl);
#endif

#endif
	/////////////// LIBIGL END /////////////////



	///////////////// FCPW /////////////////
#ifdef FCPW_TEST

	const int nVertices = V.rows();
	const int nTriangles = F.rows();
	using namespace fcpw;
	// initialize a 3d scene
	Scene<3> scene;
	
	// set the types of primitives the objects in the scene contain;
	// in this case, we have a single object consisting of only triangles
	scene.setObjectTypes({{PrimitiveType::Triangle}});
	
	// set the vertex and triangle count of the (0th) object
	scene.setObjectVertexCount(nVertices, 0);
	scene.setObjectTriangleCount(nTriangles, 0);
	
	// specify the vertex positions
	for (int i = 0; i < nVertices; i++) {
		scene.setObjectVertex(Vector<3>(V(i,0),V(i,1),V(i,2)), i, 0);
	}
	
	// specify the triangle indices
	for (int i = 0; i < nTriangles; i++) {
		scene.setObjectTriangle(&F(i,0), i, 0);
	}
	
	tictoc();
	// now that the geometry has been specified, build the acceleration structure
	scene.build(AggregateType::Bvh_SurfaceArea, true); // the second boolean argument enables vectorization
	double tPreFcpw = tictoc();
	printf("Pre process FCPW: %g s\n",tPreFcpw);
	outfile << std::endl << "Pre process FCPW: " << tPreFcpw << std::endl;
	
	printf("Starting query FCPW\n");
	
	Eigen::Matrix<double,Eigen::Dynamic,2,Eigen::RowMajor> UV(Q.rows(),2);
	Eigen::VectorXi I(Q.rows());

	double tQueryFcpwMean = 0.0;
	for (int ii = 0; ii < numTests; ++ii)
	{
		// printf("Query FCPW %d\n", ii);
		tictoc();
		igl::parallel_for(Q.rows(),[&](const int q)//(int q = 0;q<Q.rows();q++)
		{
			// perform a closest point query
			Interaction<3> interaction;
			Vector<3> queryPoint(Q(q,0),Q(q,1),Q(q,2));
			// scene.findClosestPoint(queryPoint, interaction);
			scene.findClosestPoint(queryPoint, interaction, reL2);
			I(q) = interaction.primitiveIndex;
			UV(q,0) = interaction.uv[0];
			UV(q,1) = interaction.uv[1];
		});
		double tQueryFcpw = tictoc();
		tQueryFcpwMean += tQueryFcpw;
		printf("Query FCPW (%d): %g s\n",ii, tQueryFcpw);
		printf("Pre process + Query FCPW %g s\n",tPreFcpw+tQueryFcpw);
		outfile << "Query FCPW (" << ii << "): " << tQueryFcpw << " s" << std::endl;
		outfile << "Pre process + Query FCPW: " << tPreFcpw+tQueryFcpw << " s" << std::endl;
	}
	
	tQueryFcpwMean /= numTests;
	printf("Pre process + Query FCPW Mean %g s\n\n",tPreFcpw+tQueryFcpwMean);
	outfile << "Pre process + Query FCPW Mean: " << tPreFcpw+tQueryFcpwMean << " s" << std::endl;
	
#ifdef WRITE_DMAT_FILES 
	// igl::writeDMAT("Q.dmat",Q);
	igl::writeDMAT("I.dmat",I);
	igl::writeDMAT("UV.dmat",UV);
#endif

#endif
///////////////// FCPW END /////////////////



#ifdef FPDC_TEST

	///////////////// FPDC /////////////////
	double a;
	printf("\n1:\n");
	std::cin>>a;
	mesh mesh1;
	int initOption = 1;
	mesh1.readMeshFile(argv[1],dcell, Vmin(0), Vmin(1), Vmin(2), Vmax(0), Vmax(1), Vmax(2),initOption);

	tictoc();
	mesh1.fillCandidates4();
	double tPreFpdc = tictoc();
	printf("Pre process FPDC: %g s\n",tPreFpdc);
	outfile << std::endl << "Pre process FPDC: " << tPreFpdc << std::endl;
	
	printf("Starting query FPDC\n");
	
	Eigen::VectorXi I_fpdc(Q.rows());
	Eigen::VectorXd CQD_fpdc(Q.rows());

	double tQueryFpdcMean = 0.0;
	for (int ii = 0; ii < numTests; ++ii)
	{
		tictoc();
		#pragma omp parallel for
		for(int q = 0;q<Q.rows();q++)
		{
			long elemId_fpdc;
			double closestPointDist2;
			double testPoint[3] = {Q(q,0),Q(q,1),Q(q,2)};
			double closestPoint[3] = {0.0,0.0,0.0};
			closestPointDist2 = 0; 
			elemId_fpdc = mesh1.calcDist(testPoint, reL, mesh1.getCellPosX(testPoint[0]), mesh1.getCellPosY(testPoint[1]), mesh1.getCellPosZ(testPoint[2]), closestPoint, closestPointDist2);
			I_fpdc(q) = elemId_fpdc;
			CQD_fpdc(q) = closestPointDist2;
		}
		double tQueryFpdc = tictoc();
		tQueryFpdcMean += tQueryFpdc;
		printf("Query FPDC (%d): %g s\n",ii, tQueryFpdc);
		printf("Pre process + Query FPDC %g s\n",tPreFpdc+tQueryFpdc);
		outfile << "Query FPDC (" << ii << "): " << tQueryFpdc << " s" << std::endl;
		outfile << "Pre process + Query FPDC: " << tPreFpdc+tQueryFpdc << " s" << std::endl;
	}

	tQueryFpdcMean /= numTests;
	printf("Pre process + Query FPDC Mean %g s\n\n",tPreFpdc+tQueryFpdcMean);
	outfile << "Pre process + Query FPDC Mean: " << tPreFpdc+tQueryFpdcMean << " s" << std::endl;
	mesh1.writeCellBufferToFile("id_mode0"+outfilename);

#ifdef WRITE_DMAT_FILES
	// igl::writeDMAT("Q_fpdc.dmat",Q);
	igl::writeDMAT("I_fpdc.dmat",I_fpdc);
	igl::writeDMAT("UV_fpdc.dmat",CQD_fpdc);
#endif

#endif
	///////////////// FPDC TPN /////////////////



	///////////////// FPDC (Brute force) /////////////////
#ifdef BRUTE_TEST
	
	mesh mesh2;
	mesh2.readMeshFile(argv[1],dcell, Vmin(0), Vmin(1), Vmin(2), Vmax(0), Vmax(1), Vmax(2));

	tictoc();
	mesh2.fillCandidates();

	double tPreBrute = tictoc();
	printf("Pre process Brute: %g s\n",tPreBrute);
	outfile << std::endl << "Pre process Brute: " << tPreBrute << std::endl;
	
	printf("Starting query Brute\n");

	Eigen::VectorXi I_brute(Q.rows());
	Eigen::VectorXd CQD_brute(Q.rows());
	
	double tQueryBruteMean = 0.0;
	for (int ii = 0; ii < numTests; ++ii)
	{
		tictoc();
		#pragma omp parallel for
		for(int q = 0;q<Q.rows();q++)
		{
			long elemId_brute;
			double closestPointDist2;
			double testPoint[3] = {Q(q,0),Q(q,1),Q(q,2)};
			double closestPoint[3] = {0.0,0.0,0.0};
			closestPointDist2 = 0; 
			elemId_brute = mesh2.calcDist(testPoint, reL, mesh2.getCellPosX(testPoint[0]), mesh2.getCellPosY(testPoint[1]), mesh2.getCellPosZ(testPoint[2]), closestPoint, closestPointDist2);
			I_brute(q) = elemId_brute;
			CQD_brute(q) = closestPointDist2;
		}
		double tQueryBrute = tictoc();
		tQueryBruteMean += tQueryBrute;
		printf("Query Brute (%d): %g s\n",ii, tQueryBrute);
		printf("Pre process + Query Brute %g s\n",tPreBrute+tQueryBrute);
		outfile << "Query Brute (" << ii << "): " << tQueryBrute << " s" << std::endl;
		outfile << "Pre process + Query Brute: " << tPreBrute+tQueryBrute << " s" << std::endl;
	}

	tQueryBruteMean /= numTests;
	printf("Pre process + Query Brute Mean %g s\n\n",tPreBrute+tQueryBruteMean);
	outfile << "Pre process + Query Brute Mean: " << tPreBrute+tQueryBruteMean << " s" << std::endl;

#ifdef WRITE_DMAT_FILES
	// igl::writeDMAT("Q_brute.dmat",Q);
	igl::writeDMAT("I_brute.dmat",I_brute);
	igl::writeDMAT("UV_brute.dmat",CQD_brute);
#endif

#endif
	///////////////// FPDC (Brute force) TPN /////////////////

	outfile.close(); // close the output file

}