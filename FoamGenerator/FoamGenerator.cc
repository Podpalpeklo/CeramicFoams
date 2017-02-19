/// Author: Roman Papšík <roman.papsik@vutbr.cz>
/// Date:   2016-10-23

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>

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
};

// definition of line
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
	for (unsigned int i = 0; i < keypoints.size(); i++) {
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

// Write APDL code into file
void WriteAPDL(FILE *fp, const double &x, const double &y, const double &z, const char * p)
{
	fputs("/PREP7\n", fp);

	fputs("\n!----- STRUCTURE PROPERTIES --------------------------------\n", fp);
	std::fprintf(fp, "SIZEx = %f\nSIZEy = %f\nSIZEz = %f\n", x, y, z);
	std::fprintf(fp, "LLENtot = %f\n", SumLinesLength());

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
			std::fprintf(fp, "K,%d,%.18f,%.18f,%.18f\n", k + 1, keypoints[k].x, keypoints[k].y, keypoints[k].z);
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
			std::fprintf(fp, "VA");
			for (unsigned int w = 0; w < volumes[v].size(); w++) {
				std::fprintf(fp, ",%d", volumes[v][w]);
			}
			std::fprintf(fp, "\n");
		}
	}

	fputs("/EOF", fp);
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
		std::cerr << "Incorrect number of received parameters (6 expected)" << std::endl;
		return 1;
	}

	const char * exportType = argv[1];
	const double maxX = atof(argv[2]);
	const double maxY = atof(argv[3]);
	const double maxZ = atof(argv[4]); // Coordinates of bounding box
	const char * coordinates = argv[5]; // name of file containing particles coordinates
	const char * apdl = argv[6];

	const int initial_memory = 8;
	int nx, ny, nz; // number of blocks

	voro::pre_container foam_pre_container(0, maxX + 0.0001, 0, maxY + 0.0001, 0, maxZ + 0.0001, false, false, false); // declaration
	foam_pre_container.import(coordinates); // read file with particle coordinates
	foam_pre_container.guess_optimal(nx, ny, nz); // Guess optimal number of blocks

	voro::container foam_container(0, maxX + 0.0001, 0, maxY + 0.0001, 0, maxZ + 0.0001, nx, ny, nz, false, false, false, initial_memory); // declaration
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

	RemoveDuplicates();

	FILE *output = voro::safe_fopen(apdl, "w");
	WriteAPDL(output, maxX, maxY, maxZ, exportType);
	_fcloseall();

	return 0;
}
