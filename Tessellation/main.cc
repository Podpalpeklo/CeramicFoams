/// Author: Roman Papšík
/// Email: roman.papsik@email.cz
/// Date:   2016-10-23

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>

#include "./lib_voro/voro++.hh"

// Definition of keypoint
struct keypoint {
	double x;
	double y;
	double z;

	bool operator==(const keypoint& other) const
	{
		return (std::abs(x-other.x) < 1e-12 && std::abs(y-other.y) < 1e-12 && std::abs(z-other.z) < 1e-12);
	}
	double inSurface(const double& maxX, const double& maxY, const double& maxZ) const
	{
		return (x < 1e-9 || y < 1e-9 || z < 1e-9 || std::abs(x - maxX) < 1e-9 || std::abs(y - maxY) < 1e-9 || std::abs(z - maxZ) < 1e-9);
	}
};

// Definition of line
struct line
{
	keypoint a;
	keypoint b;
	bool operator==(const line& other) const
	{
		return ((a == other.a && b == other.b) || (a == other.b && b == other.a));
	}
	double length() const
	{
		double dx = b.x - a.x;
		double dy = b.y - a.y;
		double dz = b.z - a.z;
		return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
	}
};

// Declaration of storages
std::vector<keypoint> keypoints;
std::vector<line> lines;
std::vector<std::vector<keypoint>> areas;
std::vector<std::vector<unsigned int>> volumes;

// Return index based on coordinates
unsigned int FindKeypointIndex(keypoint &k)
{
	for (unsigned int i = keypoints.size() - 1; i >= 0 ; i--) {
		if (k == keypoints[i]) return i + 1;
	}
	return -1;
}

// Sum lengths of all lines
double SumLinesLength()
{
	double totalLength = 0;
	for (unsigned int l = 0; l < lines.size(); l++) {
		totalLength = totalLength + lines[l].length();
	}
	return totalLength;
}

double CalculateStandardDeviation()
{
	double lines_total_length = SumLinesLength();
	double lines_mean_length = lines_total_length / lines.size();
	double lines_variance = 0.0;
	double lines_deviation = 0.0;

	for (unsigned int l = 0; l < lines.size(); l++)
	{
		lines_variance += pow(lines[l].length() - lines_mean_length, 2);
	}
	lines_variance = lines_variance / lines.size();

	lines_deviation = sqrt(lines_variance);
	return lines_deviation;
}

// Write APDL code into file
void WriteAPDL(FILE *fp, const double &x, const double &y, const double &z, const char * p)
{
	fputs("/PREP7\n", fp);

	fputs("\n!----- FOAM PROPERTIES --------------------------------\n", fp);
	std::fprintf(fp, "sizeX = %f\nsizeY = %f\nsizeZ = %f\n", x, y, z);
	std::fprintf(fp, "lines_total_length = %f\n", SumLinesLength());
	std::fprintf(fp, "lines_mean_length = %f\n", SumLinesLength()/lines.size());
	std::fprintf(fp, "lines_length_deviation = %f\n", CalculateStandardDeviation());

	// Write information about number of elements
	if (p[0] == 'k' || p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		std::fprintf(fp, "Kcount = %zd\n", keypoints.size());
	}
	if (p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		std::fprintf(fp, "Lcount = %zd\n", lines.size());
	}
	if (p[0] == 'a' || p[0] == 'v') {
		std::fprintf(fp, "Acount = %zd\n", areas.size());
	}
	if (p[0] == 'v') {
		std::fprintf(fp, "Vcount = %zd\n", volumes.size());
	}

	// print each type
	if (p[0] == 'k' || p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		fputs("\n!----- KEYPOINTS -------------------------------------------\n", fp);
		for (unsigned int k = 0; k < keypoints.size(); k++) {
			std::fprintf(fp, "K,%d,%.10f,%.10f,%.10f\n", k + 1, keypoints[k].x, keypoints[k].y, keypoints[k].z);
		}
	}

	if (p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		fputs("\n!----- LINES -----------------------------------------------\n", fp);
		for (unsigned int l = 0; l < lines.size(); l++) {
			std::fprintf(fp, "L,%d,%d\n", FindKeypointIndex(lines[l].a), FindKeypointIndex(lines[l].b));
		}
	}

	if (p[0] == 'a' || p[0] == 'v') {
		fputs("\n!----- AREAS -----------------------------------------------\n", fp);
		for (unsigned int a = 0; a < areas.size(); a++) {
			std::fprintf(fp, "A");
			for (unsigned int b = 0; b < areas[a].size(); b++) {
				std::fprintf(fp, ",%d", FindKeypointIndex(areas[a][b]));
			}
			std::fprintf(fp, "\n");
		}
	}

	if (p[0] == 'v') {
		fputs("\n!----- VOLUMES ---------------------------------------------\n", fp);
		for (unsigned int v = 0; v < volumes.size(); v++) {
			std::fprintf(fp, "ASEL,S,AREA,,%d\n", volumes[v][0]);
			for (unsigned int w = 1; w < volumes[v].size(); w++) {
				std::fprintf(fp, "ASEL,A,AREA,,%d\n", volumes[v][w]);
			}
			std::fprintf(fp, "VA,ALL\n");
		}
		std::fprintf(fp, "ALLSEL,ALL\n");
	}

	fputs("\n/EOF", fp);
}

void RemoveDuplicates()
{
	// keypoints
	for (unsigned int i = 1; i < keypoints.size(); i++) {
		for (unsigned int j = 0; j < i; j++) {
			if (keypoints[j] == keypoints[i]) {
				keypoints.erase(keypoints.begin() + i);
				i--;
				break;
			}
		}
	}

	// lines
	for (unsigned int i = 1; i < lines.size(); i++) {
		for (unsigned int j = 0; j < i; j++) {
			if (lines[j] == lines[i]) {
				lines.erase(lines.begin() + i);
				i--;
				break;
			}
		}
	}
}

void RemoveBorderLines(const double& maxX, const double& maxY, const double& maxZ)
{
	for (unsigned int i = 0; i < lines.size(); i++) {
		if (lines[i].a.inSurface(maxX, maxY, maxZ) && lines[i].b.inSurface(maxX, maxY, maxZ)) {
			lines.erase(lines.begin() + i);
			i--;
		}
	}
	
}

// Store faces
void storePolygon(std::vector<int> &faceVertices, std::vector<double> &v, int j)
{
	int k, l;
	int n = faceVertices[j]; // Number of vertices

	for (k = 0; k < n; k++) {
		keypoint kp;
		l = 3 * faceVertices[j + k + 1];
		kp.x = v[l];
		kp.y = v[l + 1];
		kp.z = v[l + 2];
		keypoints.push_back(kp);
	}

	for (k = 0; k < n; k++) {
		line ln;
		l = (k + 1) % n; // ensures connection of last with first
		ln.a = keypoints[keypoints.size() - n + k];
		ln.b = keypoints[keypoints.size() - n + l];
		lines.push_back(ln);
	}

	std::vector<keypoint> area;
	for (k = 0; k < n; k++) {
		area.push_back(keypoints[lines.size() - n + k]);
	}
	areas.push_back(area);
}

int main(int argc, char** argv)
{
	if (argc != 7) {
		std::cerr << "Incorrect number of received parameters (6 expected) but " << (argc - 1) << "given." << std::endl;
		return 1;
	}

	const char * exportType = argv[1];
	const double maxX = atof(argv[2]);
	const double maxY = atof(argv[3]);
	const double maxZ = atof(argv[4]); // Coordinates of bounding box
	const char * coordinates = argv[5]; // name of file containing particles coordinates
	const char * apdl = argv[6];

	char last_char_extension = coordinates[strlen(coordinates) - 1];
	if (last_char_extension == 'r')
	{
		const int initial_memory = 8;
		int nx, ny, nz; // number of blocks
		const double epsilon = 1e-10;

		static voro::pre_container_poly foam_pre_container(0, maxX + epsilon, 0, maxY + epsilon, 0, maxZ + epsilon, false, false, false); // declaration
		foam_pre_container.import(coordinates); // read file with particle coordinates
		foam_pre_container.guess_optimal(nx, ny, nz); // Guess optimal number of blocks

		static voro::container_poly foam_container(0, maxX + epsilon, 0, maxY + epsilon, 0, maxZ + epsilon, nx, ny, nz, false, false, false, initial_memory); // declaration
		foam_pre_container.setup(foam_container); // transfer of particles into proper container

		int current_cell_id;
		double x, y, z;
		voro::c_loop_all foam_loop(foam_container);
		voro::voronoicell_neighbor computedCell;
		std::vector<int> cellsNeighbours;
		std::vector<int> faceVertices;
		std::vector<double> cellsVertices; // vector of the vertex vectors

		if (foam_loop.start())
		{
			do {
				if (foam_container.compute_cell(computedCell, foam_loop)) // do if there are any particles to consider and we can compute its voronoi cells
				{
					current_cell_id = foam_loop.pid(); // id of currently considered particle
					foam_loop.pos(x, y, z); // coordinates of currently considered particle

					computedCell.neighbors(cellsNeighbours); // store id of current cell's neighbours
					computedCell.face_vertices(faceVertices); // store ids of vertices that make up a face of k vertices
					computedCell.vertices(x, y, z, cellsVertices); // store vector of vertices in global coordinate system
																   // without x,y,z would return relative position of vertices to particle

					std::vector<unsigned int> volume;
					for (unsigned int i = 0, j = 0; i < cellsNeighbours.size(); i++) // proces all faces i
					{
						storePolygon(faceVertices, cellsVertices, j);
						j += faceVertices[j] + 1;
						volume.push_back(areas.size());

					}
					volumes.push_back(volume);
				}
			} while (foam_loop.inc());
		}
	}
	else {
		const int initial_memory = 8;
		int nx, ny, nz; // number of blocks
		const double epsilon = 1e-10;

		static voro::pre_container foam_pre_container(0, maxX + epsilon, 0, maxY + epsilon, 0, maxZ + epsilon, false, false, false); // declaration
		foam_pre_container.import(coordinates); // read file with particle coordinates
		foam_pre_container.guess_optimal(nx, ny, nz); // Guess optimal number of blocks

		static voro::container foam_container(0, maxX + epsilon, 0, maxY + epsilon, 0, maxZ + epsilon, nx, ny, nz, false, false, false, initial_memory); // declaration
		foam_pre_container.setup(foam_container); // transfer of particles into proper container

		int current_cell_id;
		double x, y, z;
		voro::c_loop_all foam_loop(foam_container);
		voro::voronoicell_neighbor computedCell;
		std::vector<int> cellsNeighbours;
		std::vector<int> faceVertices;
		std::vector<double> cellsVertices; // vector of the vertex vectors

		if (foam_loop.start())
		{
			do {
				if (foam_container.compute_cell(computedCell, foam_loop)) // do if there are any particles to consider and we can compute its voronoi cells
				{
					current_cell_id = foam_loop.pid(); // id of currently considered particle
					foam_loop.pos(x, y, z); // coordinates of currently considered particle

					computedCell.neighbors(cellsNeighbours); // store id of current cell's neighbours
					computedCell.face_vertices(faceVertices); // store ids of vertices that make up a face of k vertices
					computedCell.vertices(x, y, z, cellsVertices); // store vector of vertices in global coordinate system
																   // without x,y,z would return relative position of vertices to particle

					std::vector<unsigned int> volume;
					for (unsigned int i = 0, j = 0; i < cellsNeighbours.size(); i++) // proces all faces i
					{
						storePolygon(faceVertices, cellsVertices, j);
						j += faceVertices[j] + 1;
						volume.push_back(areas.size());

					}
					volumes.push_back(volume);
				}
			} while (foam_loop.inc());
		}
	}

	// Merging keypoints and lines makes it singificantly faster in Ansys
	if (exportType[0] == 'k' || exportType[0] == 'l' || exportType[0] == 'a')
	{
		RemoveDuplicates();
	}

	// Remove border lines if no areas or volumes should be created
	if (exportType[0] == 'k' || exportType[0] == 'l')
	{
		RemoveBorderLines(maxX, maxY, maxZ);
	}

	FILE *output = voro::safe_fopen(apdl, "w");
	WriteAPDL(output, maxX, maxY, maxZ, exportType);
	_fcloseall();

	return 0;
}
