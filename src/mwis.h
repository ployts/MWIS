#ifndef _MWIS_H_
#define _MWIS_H_

#include "others.h"

struct list_node
{
	size_t L, R, U, D, label;
};

class GRAPH
{
private:
	enum MS
	{
		REMOVED,
		ADDED,
		WEIGHT_MODIFIED,
		TOKEN
	};
	enum IS
	{
		INCLUDED,
		EXCLUDED
	};

#ifdef DEBUG
	enum VS
	{
		VS_EXCLUDED, //v=1
		VS_INCLUDED, //v=0
		VS_IF,		 //v=1 if u
		VS_IFNOT,	 //v=1 if !u
		VS_NOTIF,	 //v=0 if u
		VS_NOTIFNOT, //v=0 if !u
		VS_TOKEN
	};

	struct vertex_status
	{
		vertex v, u;
		int vs;
		vertex_status(vertex v_, int vs_, vertex u_)
		{
			v = v_;
			u = u_;
			vs = vs_;
		}
	};

	stack<vertex_status> local_solution;
	stack<vertex_status> best_solution;
#endif

	struct modified_node
	{
		vertex idx;
		int status;
		weight_node weight;
		size_t lbcnt;

		modified_node(vertex idx_, int status_, weight_node weight_, size_t list_buffer_cnt_)
		{
			idx = idx_;
			status = status_;
			weight = weight_;
			lbcnt = list_buffer_cnt_;
		}

		modified_node(vertex idx_, int status_, weight_node weight_)
		{
			idx = idx_;
			status = status_;
			weight = weight_;
		}
	};
	struct branching_node
	{
		vertex idx;
		weight_node local_weight;
		int status;
		size_t pre_n;
		branching_node(vertex idx_, weight_node local_weight_, int status_, size_t pre_)
		{
			idx = idx_;
			local_weight = local_weight_;
			status = status_;
			pre_n = pre_;
		}
	};

	list_node *list_buffer;
	size_t *degree;
	size_t list_buffer_cnt;
	size_t cnt_n;
	size_t max_n;
	size_t global_n, global_m;
	size_t local_n;
	size_t pre_n;
	dynamic_bitset<> is_available;

	weight_node *weights;
	weight_node *all_neighbors_weight;
	weight_node is_weight = 0;

	stack<modified_node> modified_stack;
	stack<branching_node> branching_stack;

	dynamic_bitset<> is_in_reduction_queue[10]; //associate with reduction queue
	queue<vertex> reduction_queue[10];
#ifdef DEBUG
	dynamic_bitset<> IS_STATUS;
#endif
	bool is_clique_neighborhood;
	bool is_independent_neighborhood;

	void insert_vertex(vertex vertex_idx);
	void remove_lr(vertex &removed_vertex);
	void remove_ud(vertex &removed_vertex);
	void remove(vertex removed_vertex);
	void recover_lr(vertex &removed_vertex);
	void recover_ud(vertex &removed_vertex);
	void recover(vertex removed_vertex);

	void branch(weight_node local_weight);
	weight_node apply_all_reductions();
	void recover_reductions();
	void get_available_vertices(vector<vertex> &vertex_order);
	weight_node computing_upper_bound(vector<vertex> &vertex_order, weight_node target_weight);
	void add_to_reduction_queue(vertex idx, size_t queue_idx);
	void add_neighbors_to_reduction_queue(vertex idx, size_t queue_idx);

	weight_node (GRAPH::*reductions[10])() = {
		&GRAPH::neighborhood_reduction,
		&GRAPH::degree2_reduction,
		&GRAPH::twin_reduction,
		&GRAPH::heavy_set_reduction,
		&GRAPH::unconfined_reduction,
		&GRAPH::critical_set_reduction,
		&GRAPH::double_heavy_set_reduction};
	bool is_closed(vertex u, vertex v);
	void modify_weight(vertex v, weight_node obj_weight);
	weight_node include_vertex(vertex v);
	void exclude_vertex(vertex v);
	weight_node degree01_process(vertex v);
	weight_node fold_vertices(vertex v);
	weight_node degree2_process(vertex v);
	weight_node degree2_reduction();
	weight_node neighborhood_reduction();

	weight_node general_degree2_cases(vertex v);
	void merge_simultaneous_set(vector<vertex> &SimS);
	void get_dis2_neighbors(vertex v, vector<vertex> &neighbors);
	bool get_unconfined_set_exactly(vertex &v, vector<vertex> &S_v);
	weight_node unconfined_reduction();
	weight_node twin_reduction();
	weight_node critical_set_reduction();
	weight_node heavy_set_reduction();
	weight_node double_heavy_set_reduction();

	void split_components(vector<vector<vertex>> &components);
	weight_node branch_in_component(GRAPH *componnent_alg, vector<vertex> &component, weight_node lower_weight, bool called_flag);

	enum reductions_type
	{
		neighborhood,
		degree2,
		twin,
		heavy_set,
		unconfined,
		critical_set,
		double_heavy_set,
		independent
	};

	void update_solution(weight_node local_weight);
	void local_search(weight_node local_weigh);
	void cout_best_is_now();
	bool is_running_time_out();

	vertex find_max_available_vertex();
	weight_node local_max_weight(vector<vertex> &vertices, weight_node lower_weight);

public:
	std::chrono::duration<float> reduction_time;
	size_t remaining_vertices;
	size_t reduction_offset_all;
	bool is_called_from_split_direct = false;
	bool is_fisrt_runed = false;
	weight_node ls_lower_weight = 0;
	size_t depth;

	void build(size_t n, size_t m, pair<vertex, vertex> edges[], weight_node weights_[]);
	void build(size_t n, size_t m);
	void insert_edge(vertex &edge_from, vertex &edge_to);
	void assign_weight(vertex v, weight_node weight_v);
	weight_node MAX_IS();
	void run_branch(weight_node lower_weight);
	void reset(list_node *list_buffer_start, weight_node *weight_buffer_start, weight_node *neighbors_weight_buffer_start, size_t *degree_buffer_start);
	void reduce();
	void output_reduced_graph(string filepath);
#ifdef DEBUG
	dynamic_bitset<> export_best_is();
#endif
};

class GLOBAL_ARGS
{
private:
	bool is_inited = false;

public:
	float all_reduction_cost = 0;
	unsigned long long reductions_cnt[10] = {0};
	float reductions_tims_cnt[10] = {0};
	weight_node reductions_weight_cnt[10] = {0};
	unsigned long long first_reductions_cnt[10] = {0};
	float first_reductions_tims_cnt[10] = {0};
	weight_node first_reductions_weight_cnt[10] = {0};

	size_t N;
	size_t M;
	std::chrono::_V2::system_clock::time_point START_TIME;
	SMALL_GRAPH SG;		// MWIS solver for small graph, brute-force
	SET SET_BUFFERS[5]; //sets to be used temporarily
	pair<size_t, size_t> *EDGES;
	weight_node *WEIGHT;
	weight_node *CAP_WEIGHT; // buffers for speeding up unconfined reduction.
	size_t *REMAP;
	ISAP *FLOW_GRAPH; // MAXFLOW solver for critical independent set.
	list_node *LARGE_BUFFER;
	weight_node *WEIGHT_BUFFER;
	weight_node *NEIGHBORS_WEIGHT_BUFFER;
	size_t *DEGREE_BUFFER;
	GRAPH *GRAPH_BUFFER;

	void read(string filepath);

	~GLOBAL_ARGS();
};

void gogogogo(int argc, char *argv[]);
#endif