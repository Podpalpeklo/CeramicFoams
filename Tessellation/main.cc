/// Author: Roman Papšík
/// Email:  roman.papsik@email.cz
/// Date:   2019-03-22

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

#include "./lib_voro/voro++.hh"

// Definition of keypoint
struct keypoint
{
	double x;
	double y;
	double z;
	size_t index;

	bool operator==(const keypoint& other) const
	{
		return (std::abs(x-other.x) < 1e-12 && std::abs(y-other.y) < 1e-12 && std::abs(z-other.z) < 1e-12);
	}

	double inSurface(const double& maxX, const double& maxY, const double& maxZ) const
	{
		return (x < 1e-9 || y < 1e-9 || z < 1e-9 || std::abs(x - maxX) < 1e-9 || std::abs(y - maxY) < 1e-9 || std::abs(z - maxZ) < 1e-9);
	}
};

// Declaration of keypoints storage
std::vector<keypoint> keypoints;

// Definition of line
struct line
{
	size_t first_keypoint_index;
	size_t second_keypoint_index;

	bool operator==(const line &other) const
	{
		keypoint first_keypoint = keypoints.at(first_keypoint_index);
		keypoint second_keypoint = keypoints.at(second_keypoint_index);
		keypoint other_first_keypoint = keypoints.at(other.first_keypoint_index);
		keypoint other_second_keypoint = keypoints.at(other.second_keypoint_index);

		return ((first_keypoint == other_first_keypoint && second_keypoint == other_second_keypoint) || (first_keypoint == other_second_keypoint && second_keypoint == other_first_keypoint));
	}

	double length() const
	{
		keypoint first_keypoint = keypoints.at(first_keypoint_index);
		keypoint second_keypoint = keypoints.at(second_keypoint_index);

		double dx = second_keypoint.x - first_keypoint.x;
		double dy = second_keypoint.y - first_keypoint.y;
		double dz = second_keypoint.z - first_keypoint.z;

		return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
	}

	double inSurface(const double& maxX, const double& maxY, const double& maxZ) const
	{
		keypoint first_keypoint = keypoints.at(first_keypoint_index);
		keypoint second_keypoint = keypoints.at(second_keypoint_index);

		return (first_keypoint.inSurface(maxX, maxY, maxZ) && second_keypoint.inSurface(maxX, maxY, maxZ));
	}
};

// Declaration of lines storage
std::vector<line> lines;

// Definition of area
struct area
{
	std::vector<size_t> keypoints_indices;
	size_t index;

	bool operator==(const area &other) const
	{
		std::vector<size_t> this_keypoints_indices = keypoints_indices;
		std::vector<size_t> other_keypoints_indices = other.keypoints_indices;

		std::sort(this_keypoints_indices.begin(), this_keypoints_indices.end());
		std::sort(other_keypoints_indices.begin(), other_keypoints_indices.end());

		return (this_keypoints_indices == other_keypoints_indices);
	}
};

// Declaration of areas storage
std::vector<area> areas;

// Definition of volume
struct volume
{
	std::vector<size_t> areas_indices;
};

// Declaration of volumes storage
std::vector<volume> volumes;

/*/// Return index of keypoint
size_t FindKeypointIndexByCoordinates(keypoint &sought_keypoint)
{
	for (size_t k = 0; k < keypoints.size(); k++) {
		if (keypoints.at(k) == sought_keypoint) {
			return keypoints.at(k).index;
		}
	}

	return -1;
}*/

/// Sum lengths of all lines
double SumLinesLength()
{
	double lines_total_length = 0.0;

	for (size_t l = 0; l < lines.size(); l++) {
		lines_total_length += lines.at(l).length();
	}

	return lines_total_length;
}

/// Calculate standard deviation of lines' length [mm]
double CalculateStandardDeviation()
{
	double lines_total_length = SumLinesLength();
	double lines_mean_length = lines_total_length / lines.size();
	double lines_variance = 0.0;
	double lines_deviation = 0.0;

	for (size_t l = 0; l < lines.size(); l++) {
		lines_variance += pow(lines.at(l).length() - lines_mean_length, 2);
	}

	lines_variance = lines_variance / lines.size();
	lines_deviation = sqrt(lines_variance);

	return lines_deviation;
}

/// Removes duplicates of lines and keypoints from storage
void RemoveDuplicates()
{
	// Merge coincident lines by shifting keypoins indices
	std::cout << "Coincident lines are being merged." << std::endl;
	for (size_t line = 0; line < lines.size(); line++) {
		keypoint current_first_keypoint = keypoints.at(lines.at(line).first_keypoint_index);
		for (size_t kp = 0; kp < keypoints.size(); kp++) {
			keypoint tested_keypoint = keypoints.at(kp);
			if (tested_keypoint == current_first_keypoint) {
				lines.at(line).first_keypoint_index = tested_keypoint.index;
				break;
			}
		}
		keypoint current_second_keypoint = keypoints.at(lines.at(line).second_keypoint_index);
		for (size_t kp = 0; kp < keypoints.size(); kp++) {
			keypoint tested_keypoint = keypoints.at(kp);
			if (tested_keypoint == current_second_keypoint) {
				lines.at(line).second_keypoint_index = tested_keypoint.index;
				break;
			}
		}
	}

	// Merge coincident areas by shifting keypoints indices
	std::cout << "Coincident areas are being merged." << std::endl;
	for (size_t area = 0; area < areas.size(); area++)
	{
		for (size_t area_keypoint = 0; area_keypoint < areas.at(area).keypoints_indices.size(); area_keypoint++)
		{
			for (size_t tested_keypoint = 0; tested_keypoint < keypoints.size(); tested_keypoint++)
			{
				if (keypoints.at(tested_keypoint) == keypoints.at(areas.at(area).keypoints_indices.at(area_keypoint)))
				{
					areas.at(area).keypoints_indices.at(area_keypoint) = keypoints.at(tested_keypoint).index;
					break;
				}
			}
		}
	}

	// Adjust volumes to indices of merged areas
	std::cout << "Adjusting volumes to merged areas." << std::endl;
	for (size_t volume = 0; volume < volumes.size(); volume++)
	{
		for (size_t volume_area = 0; volume_area < volumes.at(volume).areas_indices.size(); volume_area++)
		{
			for (size_t tested_area = 0; tested_area < areas.size(); tested_area++)
			{
				if (areas.at(tested_area) == areas.at(volumes.at(volume).areas_indices.at(volume_area)))
				{
					volumes.at(volume).areas_indices.at(volume_area) = areas.at(tested_area).index;
					break;
				}
			}
		}
	}

	// Delete duplicate areas
	std::cout << "Duplicate areas are being deleted." << std::endl;
	for (size_t area = 0; area < areas.size(); area++)
	{
		for(size_t tested_area = area + 1; tested_area < areas.size(); tested_area++)
		{
			if(areas.at(tested_area) == areas.at(area))
			{
				areas.erase(areas.begin() + tested_area);
				tested_area--;
			}
		}
	}

	// Delete duplicate lines
      std::cout << "Duplicate lines are being deleted." << std::endl;
	for (size_t line = 0; line < lines.size(); line++)
	{
		for (size_t tested_line = line + 1; tested_line < lines.size(); tested_line++)
		{
			if (lines.at(tested_line) == lines.at(line))
			{
				lines.erase(lines.begin() + tested_line);
				tested_line--;
			}
		}
	}

	// Delete unused coincident keypoints
	std::cout << "Unused coincident keypoints are being deleted." << std::endl;
	for (size_t kp = 0; kp < keypoints.size(); kp++)
	{
		for (size_t tested_kp = kp + 1; tested_kp < keypoints.size(); tested_kp++) {
			if (keypoints.at(tested_kp) == keypoints.at(kp)) {
				keypoints.erase(keypoints.begin() + tested_kp);
				tested_kp--;
			}
		}
	}
}

/// Removes lines, which are on the boundary
void RemoveBorderLines(const double& maxX, const double& maxY, const double& maxZ)
{
	std::cout << "Lines on boundary are being deleted." << std::endl;

	for (size_t i = 0; i < lines.size(); i++) {
		if (lines.at(i).inSurface(maxX, maxY, maxZ)) {
			lines.erase(lines.begin() + i);
			i--;
		}
	}
	
}

/// This function is called for each generated cell and it processes it into geometrical entities
void storeCell(std::vector<int> &cellsNeighbours, std::vector<int> &faceVertices, std::vector<double> &coordinates_of_vertices_of_cell)
{
	// Initialize new volume of cell
	volume new_volume;
	// ?
	int j = 0;
	// Do for each face of the cell
	for (size_t face = 0; face < cellsNeighbours.size(); face++)
	{
		// Number of vertices
		int n = faceVertices.at(j); 

		for (int k = 0; k < n; k++) {
			keypoint new_keypoint;
			int l = 3 * faceVertices.at(j + k + 1);
			new_keypoint.x = coordinates_of_vertices_of_cell.at(l);
			new_keypoint.y = coordinates_of_vertices_of_cell.at(l + 1);
			new_keypoint.z = coordinates_of_vertices_of_cell.at(l + 2);
			new_keypoint.index = keypoints.size();
			keypoints.push_back(new_keypoint);
		}

		for (int k = 0; k < n; k++) {
			line new_line;
			int l = (k + 1) % n; // ensures connection of last with first
			new_line.first_keypoint_index = keypoints.size() - n + k;
			new_line.second_keypoint_index = keypoints.size() - n + l;
			lines.push_back(new_line);
		}

		area new_area;
		for (int k = 0; k < n; k++) {
			new_area.keypoints_indices.push_back(lines.size() - n + k);
		}
		new_area.index = areas.size();
		areas.push_back(new_area);

		new_volume.areas_indices.push_back(areas.back().index);

		j += faceVertices.at(j) + 1;
	}
	volumes.push_back(new_volume);
}

/// Write APDL code into file
void WriteAPDL(FILE *fp, const double &x, const double &y, const double &z, const char * p)
{
	std::cout << "Input file for APDL is being written." << std::endl;
	fputs("/PREP7\n", fp);

	fputs("\n!----- FOAM PROPERTIES --------------------------------\n", fp);
	std::fprintf(fp, "sizeX = %f\nsizeY = %f\nsizeZ = %f\n", x, y, z);
	// std::fprintf(fp, "lines_total_length = %f\n", SumLinesLength());
	// std::fprintf(fp, "lines_mean_length = %f\n", SumLinesLength()/lines.size());
	// std::fprintf(fp, "lines_length_deviation = %f\n", CalculateStandardDeviation());

	// Write information about number of elements
	if (p[0] == 'k' || p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		std::fprintf(fp, "Kcount = %zu\n", keypoints.size());
	}
	if (p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		std::fprintf(fp, "Lcount = %zu\n", lines.size());
	}
	if (p[0] == 'a' || p[0] == 'v') {
		std::fprintf(fp, "Acount = %zu\n", areas.size());
	}
	if (p[0] == 'v') {
		std::fprintf(fp, "Vcount = %zu\n", volumes.size());
	}

	// print each type
	if (p[0] == 'k' || p[0] == 'l' || p[0] == 'a' || p[0] == 'v') {
		fputs("\n!----- KEYPOINTS -------------------------------------------\n", fp);
		for (size_t k = 0; k < keypoints.size(); k++) {
			std::fprintf(fp, "K,%zu,%.10f,%.10f,%.10f\n", keypoints.at(k).index + 1, keypoints.at(k).x, keypoints.at(k).y, keypoints.at(k).z);
		}
	}

	if (p[0] == 'l') {
		fputs("\n!----- LINES -----------------------------------------------\n", fp);
		for (size_t l = 0; l < lines.size(); l++) {
			std::fprintf(fp, "L,%zu,%zu\n", lines.at(l).first_keypoint_index + 1, lines.at(l).second_keypoint_index + 1);
		}
	}

	if (p[0] == 'a' || p[0] == 'v') {
		fputs("\n!----- AREAS -----------------------------------------------\n", fp);
		for (size_t a = 0; a < areas.size(); a++)
		{
			// Areas jsou seøazené APDLkem jinak, než jsou jejich indexy - to musíš øešit pøi výpisu jejich hledáním
			std::fprintf(fp, "NUMSTR,AREA,%zu $ A", areas.at(a).index + 1);
			for (size_t k = 0; k < areas.at(a).keypoints_indices.size(); k++) {
				std::fprintf(fp, ",%zu", areas.at(a).keypoints_indices.at(k) + 1);
			}
			std::fprintf(fp, "\n");
		}
	}

	if (p[0] == 'v') {
		fputs("\n!----- VOLUMES ---------------------------------------------\n", fp);
		for (size_t v = 0; v < volumes.size(); v++) {
			std::fprintf(fp, "ASEL,S,AREA,,%zu\n", volumes.at(v).areas_indices.at(0) + 1);
			for (size_t a = 1; a < volumes.at(v).areas_indices.size(); a++) {
				std::fprintf(fp, "ASEL,A,AREA,,%zu\n", volumes.at(v).areas_indices.at(a) + 1);
			}
			//std::fprintf(fp, "LSLA,S $ KSLL,S $ NUMMRG,KP $ VA,ALL\n");
			std::fprintf(fp, "VA,ALL\n");
		}
	}

	//std::fprintf(fp, "ALLSEL,ALL $ NUMMRG,ALL\n");
	std::fprintf(fp, "ALLSEL,ALL\n");
	fputs("\n/EOF", fp);
}

int main(int argc, char** argv)
{
	if (argc != 7) {
		std::cerr << "Incorrect number of received parameters. Six expected, but only " << (argc - 1) << "given." << std::endl;
		return 1;
	}

	const char* exportType = argv[1];
	const double maxX = atof(argv[2]);
	const double maxY = atof(argv[3]);
	const double maxZ = atof(argv[4]); // Coordinates of bounding box
	const char* coordinates = argv[5]; // name of file containing particles coordinates
	const char* apdl = argv[6];

	std::cout << "Size of specimen is " << maxX << "x" << maxY << "x" << maxZ << " mm." << std::endl;

	char last_char_extension = coordinates[strlen(coordinates) - 1];
	if (last_char_extension == 'r')
	{
		std::cout << "Laguerre tesselation from nuclei is being created." << std::endl;

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
		std::vector<double> coordinates_of_vertices_of_cell; // vector of the vertex vectors

		if (foam_loop.start())
		{
			do {
				if (foam_container.compute_cell(computedCell, foam_loop)) // do if there are any particles to consider and we can compute its voronoi cells
				{
					current_cell_id = foam_loop.pid(); // id of currently considered particle
					foam_loop.pos(x, y, z); // coordinates of currently considered particle

					computedCell.neighbors(cellsNeighbours); // store id of current cell's neighbours
					computedCell.face_vertices(faceVertices); // store ids of vertices that make up a face of k vertices
					computedCell.vertices(x, y, z, coordinates_of_vertices_of_cell); // store vector of vertices in global coordinate system (without x,y,z would return relative position of vertices to particle)

					storeCell(cellsNeighbours, faceVertices, coordinates_of_vertices_of_cell);
				}
			} while (foam_loop.inc());
		}
	}
	else
	{
		std::cout << "Voronoi tessellation from nuclei is being created." << std::endl;

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
		std::vector<double> coordinates_of_vertices_of_cell; // vector of the vertex vectors

		if (foam_loop.start())
		{
			do {
				if (foam_container.compute_cell(computedCell, foam_loop)) // do if there are any particles to consider and we can compute its voronoi cells
				{
					current_cell_id = foam_loop.pid(); // id of currently considered particle
					foam_loop.pos(x, y, z); // coordinates of currently considered particle

					computedCell.neighbors(cellsNeighbours); // store id of current cell's neighbours
					computedCell.face_vertices(faceVertices); // store ids of vertices that make up a face of k vertices
					computedCell.vertices(x, y, z, coordinates_of_vertices_of_cell); // store vector of vertices in global coordinate system (without x,y,z would return relative position of vertices to particle)

					storeCell(cellsNeighbours, faceVertices, coordinates_of_vertices_of_cell);
				}
			} while (foam_loop.inc());
		}
	}

	// Merging solid geometry singificantly impoves speed in Ansys
	RemoveDuplicates();

	/*// Remove border lines if no areas or volumes should be created
	if (exportType[0] == 'l')
	{
		RemoveBorderLines(maxX, maxY, maxZ);
	}*/

	FILE *output = voro::safe_fopen(apdl, "w");
	WriteAPDL(output, maxX, maxY, maxZ, exportType);
	_fcloseall();

	return 0;
}
