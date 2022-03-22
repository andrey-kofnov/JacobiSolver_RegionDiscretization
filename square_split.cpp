#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <assert.h> 

bool is_prime(int n) {

	bool res = true;

	if (!n || n == 1) {
		res = false;
	}

	for (auto i = 2; i <= n / 2; ++i) {
		if (!(n % i)) {
			res = false;
			break;
		}
	}
	return res;
}


void optimal_set(int* out, int n) {

	int res[2] = { -1, -1 };

	int min_diff = 32000;
	int a = -1;
	int b = -1;

	for (auto i = 2; i <= n / 2; ++i) {
		if (!(n % i)) {
			if (abs((n / i) - i) < min_diff) {
				min_diff = abs((n / i) - i);
				a = (n / i);
				b = i;
			}
		}
	}
	res[0] = a; // axis Y
	res[1] = b; // axis X
	std::cout << res[0] << ' ' << res[1] << std::endl;
	*out = res[0];
	*(++out) = res[1];
}

class bounds {
	double start;
	double stop;
public:
	bounds(double _start = 0.0, double _stop = 0.0) {
		this->start = _start;
		this->stop = _stop;
	}

	bounds(const bounds& old_obj) {
		this->start = old_obj.start;
		this->stop = old_obj.stop;
	}

	double _start() { return this->start; }
	double _stop() { return this->stop; }
	~bounds() {}
};


typedef struct Point {
	double x = 0;
	double y = 0;
	Point* west = nullptr;
	Point* east = nullptr;
	Point* south = nullptr;
	Point* north = nullptr;
	int owner = -1;
	int ghost_for = -1;
} Point;

class Point_ptr {
	Point* A;
public:
	Point_ptr(double x, double y, int _owner = -1, int _ghost_for = -1) {
		this->A = new Point;
		this->A->x = x;
		this->A->y = y;
		this->A->owner = _owner;
		this->A->ghost_for = _ghost_for;
	}

	double _x() { return this->A->x; }
	double _y() { return this->A->y; }
	int _owner() { return this->A->owner; }
	int _ghost_for() { return this->A->ghost_for; }

	Point& operator* () {
		return *(this->A);
	}

	Point* operator-> () {
		return this->A;
	}

	Point* ptr() {
		return this->A;
	}

	~Point_ptr() {
		this->A->east = nullptr;
		this->A->west = nullptr;
		this->A->north = nullptr;
		this->A->south = nullptr;
		// delete[] (this->A);
	}
};

class Process_params {
	int proc_id;
	std::vector<Point_ptr> grid;
	std::vector<Point_ptr> ghost_layer;
	int south_neighbour;
	int north_neighbour;
	int west_neighbour;
	int east_neighbour;
public:
	Process_params(std::vector<Point_ptr> _grid,
		std::vector<Point_ptr> _ghost_layer,
		int _proc_id = -1,
		int _south_neighbour = -1, 
		int _north_neighbour = -1, 
		int _west_neighbour = -1,
		int _east_neighbour = -1) {

		this->proc_id= _proc_id;
		this->grid = _grid;
		this->ghost_layer = _ghost_layer;
		this->south_neighbour = _south_neighbour;
		this->north_neighbour = _north_neighbour;
		this->west_neighbour = _west_neighbour;
		this->east_neighbour = _east_neighbour;
	}
	Process_params(int _proc_id = -1,
		int _south_neighbour = -1,
		int _north_neighbour = -1,
		int _west_neighbour = -1,
		int _east_neighbour = -1) {
		this->proc_id = _proc_id;
		this->south_neighbour = _south_neighbour;
		this->north_neighbour = _north_neighbour;
		this->west_neighbour = _west_neighbour;
		this->east_neighbour = _east_neighbour;
	}

	int get_id() { return this->proc_id; }
	int get_south_neighbour() { return this->south_neighbour; }
	int get_north_neighbour() { return this->north_neighbour; }
	int get_east_neighbour() { return this->east_neighbour; }
	int get_west_neighbour() { return this->west_neighbour; }
	void set_myself(int x) { this->proc_id = x; }
	void set_south_neighbour(int x) { this->south_neighbour = x; }
	void set_north_neighbour(int x) { this->north_neighbour = x; }
	void set_east_neighbour(int x) { this->east_neighbour = x; }
	void set_west_neighbour(int x) { this->west_neighbour = x; }

	std::vector<Point_ptr> get_grid() { return this->grid; }
	std::vector<Point_ptr> get_ghost_layer() { return this->ghost_layer; }

	void add_grid_elem(Point_ptr x) { this->grid.push_back(x); }
	void add_ghost_layer_elem(Point_ptr x) { this->ghost_layer.push_back(x); }


	~Process_params() {
		this->grid.clear();
		this->ghost_layer.clear();
	}
};

class Output_features {
	int init_num_proc;
	int num_proc;
	int opt_set[2];
	std::string mod;
	std::vector<Process_params> params;

	std::vector < Point_ptr> points;
public:
	Output_features(std::vector<Process_params> _params, 
		int _init_num_proc = 1, 
		int _num_proc = 1, 
		std::string _mod = "1D") {

		this->params = _params;
		this->init_num_proc = _init_num_proc;
		this->num_proc = _num_proc;
		this->mod = _mod;
		if (this->mod == "1D") {
			this->opt_set[0] = num_proc;
			this->opt_set[1] = 1;
		}
		else {
			optimal_set(&((this->opt_set)[0]), num_proc);
		}
	}

	Output_features(int _init_num_proc = 1, int _num_proc = 1) {
		this->init_num_proc = _init_num_proc;
		this->num_proc = _num_proc;
		this->mod = "1D";
		this->opt_set[0] = num_proc;
		this->opt_set[1] = 1;
	}

	int get_init_num_proc() { return this->init_num_proc; }
	int get_num_proc() { return this->num_proc; }
	int* get_opt_set() { return this->opt_set; }
	std::string get_mod() { return this->mod; }

	std::vector<Process_params>
		get_params() { return this->params; }

	void add_params(Process_params x) { this->params.push_back(x); }
	void set_points(std::vector<Point_ptr> _points) { this->points = _points; }

	std::vector<Point_ptr> get_points() { 
		return this->points; 
	}

	void add_Point(int proc_id, Point_ptr x) {
		this->params[proc_id].add_grid_elem(x);
	}

	void add_Ghost_Point(int proc_id, Point_ptr x) {
		this->params[proc_id].add_ghost_layer_elem(x);
	}

	void replace_proc_params(Process_params x, int id) {
		Process_params tmp = this->params[id];
		this->params[id] = x;
	}

	std::vector<Point_ptr>
		get_proc_grid(int id) {
		return this->params[id].get_grid();
	}

	std::vector<Point_ptr>
		get_proc_ghost_layer(int id) {
		return this->params[id].get_ghost_layer();
	}

	void set_opt_set(int* a) {
		this->opt_set[0] = a[0];
		this->opt_set[1] = a[1];
	}

	void set_opt_set(int a, int b) {
		this->opt_set[0] = a;
		this->opt_set[1] = b;
	}

	void set_mod(std::string s) { this->mod = s; }


	~Output_features() {
		this->params.clear();
	}
};



int main(void) {


	std::cout << "Dimension: (1D or 2D) ";
	std::string s;
	std::cin >> s;
	assert(s == "1D" || s == "2D");
	std::string mod = s;

	double H = 1.0;
	double h;
	int resolution;

	std::cout << "Resolution is ";
	std::cin >> resolution;
	h = H / (resolution - 1);

	int init_num_proc, num_proc;

	std::cout << "Number of processes: ";
	std::cin >> init_num_proc;
	num_proc = init_num_proc;
	int os[2] = { -1, -1 };

	if (mod == "1D") {
		if (init_num_proc >= resolution - 2) {
			num_proc = resolution - 2;
		}
		os[0] = num_proc;
		os[1] = 1;
	}
	else {
		if (init_num_proc >= (resolution - 2) * (resolution - 2)) {
			num_proc = (resolution - 2) * (resolution - 2);
		}
		if (is_prime(num_proc)) {
			mod = "1D";
			os[0] = num_proc;
			os[1] = 1;
		}
		else {
			optimal_set(&os[0], num_proc);
		}
	}

	// Check
	std::cout << "Current mod: " << mod << "; Number of processes: " << num_proc << std::endl;
	std::cout << "Optimal set: " << os[0] << ", " << os[1] << std::endl;
	

	std::vector<double> PointsX;
	std::vector<double> PointsY;

	double p = 0.0;

	for (auto i = 0; i < resolution; ++i) {
		PointsX.push_back(p);
		PointsY.push_back(p);
		p += h;
	}

	std::vector<Point_ptr> points;
	for (std::vector<double>::iterator itX = PointsX.begin(); itX != PointsX.end(); ++itX) {
		for (std::vector<double>::iterator itY = PointsY.begin(); itY != PointsY.end(); ++itY) {
			Point_ptr tmp(*itX, *itY);
			points.push_back(tmp);
		}
	}

	// Initially here

	


	std::vector<double> Points;
	p = 0.0;
	for (auto i = 0; i < resolution; ++i) {
		Points.push_back(p);
		p += h;
	}

	// Check
	/*
	std::cout << "Points: ";
	for (std::vector<double>::iterator it = Points.begin(); it != Points.end(); ++it) {
		std::cout << *it << ' ';
	}
	std::cout << std::endl;
	*/

	///////////////////////////////////////////////////////////////////////////////////

	std::vector<bounds> proc_pointsY; // os[0]
	std::vector<bounds> proc_pointsX; // os[1]

	// Points for processes along the axis Y   : os[0] 

	// Each Process operates over k or k+1 points
	int c_num_proc = os[0];
	int k = int((resolution - 2) / c_num_proc);
	int resid = (resolution - 2) % c_num_proc;
	std::vector<int> num_pointsY;
	for (auto i = 0; i < c_num_proc; ++i) {
		if (resid > 0) {
			num_pointsY.push_back(k + 1);
			--resid;
		}
		else {
			num_pointsY.push_back(k);
		}
	}

	int processed_points = 0;
	for (auto i = 0; i < c_num_proc; ++i) {
		auto cur_point = Points.begin() + processed_points;
		bounds tmp(*cur_point, *(cur_point + num_pointsY[i] + 1));
		processed_points += num_pointsY[i];
		proc_pointsY.push_back(tmp);
	}

	// Points for processes along the axis X  : os[1] 

	// Each Process operates over k or k+1 points
	c_num_proc = os[1];
	k = int((resolution - 2) / c_num_proc);
	resid = (resolution - 2) % c_num_proc;
	std::vector<int> num_pointsX;
	for (auto i = 0; i < c_num_proc; ++i) {
		if (resid > 0) {
			num_pointsX.push_back(k + 1);
			--resid;
		}
		else {
			num_pointsX.push_back(k);
		}
	}

	processed_points = 0;
	for (auto i = 0; i < c_num_proc; ++i) {
		auto cur_point = Points.begin() + processed_points;
		bounds tmp(*cur_point, *(cur_point + num_pointsX[i] + 1));
		processed_points += num_pointsX[i];
		proc_pointsX.push_back(tmp);
	}

	num_pointsX.clear();
	num_pointsY.clear();


	Output_features output(init_num_proc, num_proc);
	output.set_mod(mod);
	output.set_opt_set(os);
	output.set_points(points);
	for (int i = 0; i < os[1]; i++) {
		for (int j = 0; j < os[0]; j++) {
			Process_params tmp;
			tmp.set_myself(i * os[0] + j);
			if (i < os[1] - 1){ tmp.set_east_neighbour((i + 1) * os[0] + j); }	
			if (i > 0) { tmp.set_west_neighbour((i - 1) * os[0] + j); }
			if (j < os[0] - 1) { tmp.set_north_neighbour(i * os[0] + j + 1); }
			if (j > 0) { tmp.set_south_neighbour(i * os[0] + j - 1); }
			output.add_params(tmp);
		}
	}

	for (auto i = 0; i < PointsX.size(); i++) {
		for (auto j = 0; j < PointsY.size(); j++) {
			for (int l = 0; l < os[1]; l++) {
				for (int k = 0; k < os[0]; k++) {
					if (proc_pointsY[k]._start() < PointsY[j] &&
						proc_pointsY[k]._stop() > PointsY[j] &&
						proc_pointsX[l]._start() < PointsX[i] &&
						proc_pointsX[l]._stop() > PointsX[i]) {

						points[i * PointsY.size() + j]->owner = l * os[0] + k;
						output.add_Point(l* os[0] + k, points[i * PointsY.size() + j]);

					}
					if (proc_pointsY[k]._start() < PointsY[j] &&
						proc_pointsY[k]._stop() > PointsY[j] ) {
						if (proc_pointsX[l]._start() == PointsX[i] || 
							proc_pointsX[l]._stop() == PointsX[i]) {

							points[i * PointsY.size() + j]->ghost_for = l * os[0] + k;
							output.add_Ghost_Point(l * os[0] + k, points[i * PointsY.size() + j]);
						}
					}
					if (proc_pointsX[l]._start() < PointsX[i] &&
						proc_pointsX[l]._stop() > PointsX[i]) {
						if (proc_pointsY[k]._start() == PointsY[j] ||
							proc_pointsY[k]._stop() == PointsY[j]) {

							points[i * PointsY.size() + j]->ghost_for = l * os[0] + k;
							output.add_Ghost_Point(l* os[0] + k, points[i * PointsY.size() + j]);
						}
					}
				}
			}
			if (j > 0) {
				points[i * PointsY.size() + j]->south = points[i * PointsY.size() + j - 1].ptr();
			}
			if (j < PointsY.size() - 1) {
				points[i * PointsY.size() + j]->north = points[i * PointsY.size() + j + 1].ptr();
			}
			if (i > 0) {
				points[i * PointsY.size() + j]->west = points[(i - 1) * PointsY.size() + j].ptr();
			}
			if (i < PointsX.size() - 1) {
				points[i * PointsY.size() + j]->east = points[(i + 1) * PointsY.size() + j].ptr();
			}
		}
	}
	
	// Check
	/*
	std::cout << "Points: " << std::endl;
	std::vector<Point_ptr> pt = output.get_points();
	for (std::vector<Point_ptr>::iterator it = pt.begin(); it != pt.end(); ++it) {
		std::cout << (*it)->x << ' ' << (*it)->y << std::endl;
		std::cout << "OWNER: " << (*it)->owner << std::endl;
		std::cout << "GHOST FOR: " << (*it)->ghost_for << std::endl;
		if ((*it)->east != nullptr) {
			std::cout << "EAST: " << (*it)->east->x << ' ' << (*it)->east->y << std::endl;
		}
		if ((*it)->west != nullptr) {
			std::cout << "WEST: " << (*it)->west->x << ' ' << (*it)->west->y << std::endl;
		}
		if ((*it)->north != nullptr) {
			std::cout << "NORTH: " << (*it)->north->x << ' ' << (*it)->north->y << std::endl;
		}
		if ((*it)->south != nullptr) {
			std::cout << "SOUTH: " << (*it)->south->x << ' ' << (*it)->south->y << std::endl;
		}

		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/


	std::cout << "-------------- Output params -----------------" << std::endl;
	std::cout << "Output: Initial number of processes: " << output.get_init_num_proc() << std::endl;
	std::cout << "Output: Number of processes: " << output.get_num_proc() << std::endl;
	std::vector<Process_params> par = output.get_params();
	for (std::vector<Process_params>::iterator it = par.begin(); it != par.end(); ++it) {
		std::cout << "Output: Process id: " << (*it).get_id() << std::endl;
		std::cout << "Output: Ghost layer: " << (*it).get_ghost_layer().size() << std::endl;
		std::cout << "Output: Grid: ";
		std::vector<Point_ptr> cur_grid = (*it).get_grid();
		for (std::vector<Point_ptr>::iterator itr = cur_grid.begin(); itr != cur_grid.end(); ++itr) {
			std::cout << ' ' << ( *itr)->x << ';' << (*itr)->y;
		}
		std::cout << std::endl;
		std::cout << "Output: South neighbour: " << (*it).get_south_neighbour() << "; North neighbour: " << (*it).get_north_neighbour() << std::endl;
		std::cout << "Output: West neighbour: " << (*it).get_west_neighbour() << "; East neighbour: " << (*it).get_east_neighbour() << std::endl;
	}
	std::cout << "Output: Optimal set: " << output.get_opt_set()[0] << ' ' << output.get_opt_set()[1] << std::endl;
	std::cout << "Output: Current mod: " << output.get_mod() << std::endl;
	return 0;
}