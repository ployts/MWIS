#include "mwis.h"
#include <unistd.h>
#include <getopt.h>

//Global arguments to be used.

opt_arg OA;

//Statistics
float all_reduction_cost = 0;
unsigned long long reductions_cnt[10] = {0};
float reductions_tims_cnt[10] = {0};
weight_node reductions_weight_cnt[10] = {0};
unsigned long long first_reductions_cnt[10] = {0};
float first_reductions_tims_cnt[10] = {0};
weight_node first_reductions_weight_cnt[10] = {0};

SMALL_GRAPH small_graph; // MWIS solver for small graph, brute-force
SET set_buffer0;		 //sets to be used temporarily
SET set_buffer1;		 //sets to be used temporarily
SET set_buffer2;		 //sets to be used temporarily
SET set_buffer3;		 //sets to be used temporarily
SET set_buffer4;		 //sets to be used temporarily
pair<size_t, size_t> *init_edges;
weight_node *init_weights;
weight_node *cap_weight; // buffers for speeding up unconfined reduction.
std::chrono::_V2::system_clock::time_point start_time;
size_t *REMAP;
ISAP *flow_graph; // MAXFLOW solver for critical independent set.
list_node *LARGE_BUFFER;
weight_node *WEIGHT_BUFFER;
weight_node *NEIGHBORS_WEIGHT_BUFFER; //
size_t *DEGREE_BUFFER;
GRAPH *GRAPH_BUFFER;

void GRAPH::build(size_t n, size_t m, pair<vertex, vertex> edges[], weight_node weights_[])
{
	list_buffer_cnt = 4 * n + 1;
	cnt_n = n + 1;
	max_n = 2 * n;
	global_n = n;
	global_m = m;
	max_n = 2 * n;
	local_n = n;
	is_available.resize(max_n, false);
	for (auto i = 0; i < 10; i++)
		is_in_reduction_queue[i].resize(max_n, false);
#ifdef DEBUG
	IS_STATUS.resize(max_n, false);
#endif

	for (size_t i = 1; i <= global_n; i++)
	{
		LBL(i) = i - 1;
		LBR(i) = i + 1;
		LBU(i) = LBD(i) = i;
		LBLB(i) = 0;

		LBL(i + max_n) = LBR(i + max_n) = i + max_n;
		LBU(i + max_n) = i + max_n - 1;
		LBD(i + max_n) = i + max_n + 1;
		LBLB(i + max_n) = i;
		is_available[i] = true;
		all_neighbors_weight[i] = 0;
		weights[i] = weights_[i];
		degree[i] = 0;
	}

	LBL(0) = global_n;
	LBR(0) = 1;
	LBU(0) = max_n + global_n;
	LBD(0) = 1 + max_n;

	LBLB(0) = 0;
	size_t head = 0;
	recover_lr(head);
	recover_ud(head);

	for (size_t i = 0; i < global_m; i++)
	{
		insert_edge(edges[i].first, edges[i].second);
	}
}

void GRAPH::build(size_t n, size_t m)
{
	list_buffer_cnt = 4 * n + 1;
	cnt_n = n + 1;
	max_n = 2 * n;
	global_n = n;
	global_m = m;
	max_n = 2 * n;
	local_n = n;

	is_available.resize(max_n, false);
	for (auto i = 0; i < 10; i++)
		is_in_reduction_queue[i].resize(max_n, false);
#ifdef DEBUG
	IS_STATUS.resize(max_n, false);
#endif

	for (size_t i = 1; i <= global_n; i++)
	{
		LBL(i) = i - 1;
		LBR(i) = i + 1;
		LBU(i) = LBD(i) = i;
		LBLB(i) = 0;

		LBL(i + max_n) = LBR(i + max_n) = i + max_n;
		LBU(i + max_n) = i + max_n - 1;
		LBD(i + max_n) = i + max_n + 1;
		LBLB(i + max_n) = i;
		is_available[i] = true;
		all_neighbors_weight[i] = 0;
		degree[i] = 0;
	}

	LBL(0) = global_n;
	LBR(0) = 1;
	LBU(0) = max_n + global_n;
	LBD(0) = 1 + max_n;

	LBLB(0) = 0;
	size_t head = 0;
	recover_lr(head);
	recover_ud(head);
}

void GRAPH::assign_weight(vertex v, weight_node weight_v)
{
	weights[v] = weight_v;
}

inline void GRAPH::insert_vertex(vertex vertex_idx)
{
	add_to_reduction_queue(vertex_idx, 0);
	all_neighbors_weight[vertex_idx] = 0;
	local_n++;
	degree[vertex_idx] = 0;
	is_available[vertex_idx] = true;
	LBL(vertex_idx) = LBL(0);
	LBR(vertex_idx) = 0;
	LBU(vertex_idx) = LBD(vertex_idx) = vertex_idx;
	LBLB(vertex_idx) = 0;
	recover_ud(vertex_idx);
	recover_lr(vertex_idx);
	vertex_idx += max_n;
	LBU(vertex_idx) = LBU(0);
	LBD(vertex_idx) = 0;
	LBL(vertex_idx) = LBR(vertex_idx) = vertex_idx;
	LBLB(vertex_idx) = vertex_idx - max_n;
	recover_ud(vertex_idx);
	recover_lr(vertex_idx);
}

inline void GRAPH::insert_edge(vertex &edge_from, vertex &edge_to)
{
	all_neighbors_weight[edge_from] += weights[edge_to];
	LBD(list_buffer_cnt) = edge_from;
	LBU(list_buffer_cnt) = LBU(edge_from);
	LBR(list_buffer_cnt) = edge_to + max_n;
	LBL(list_buffer_cnt) = LBL(edge_to + max_n);
	LBLB(list_buffer_cnt) = edge_to;
	recover_ud(list_buffer_cnt);
	recover_lr(list_buffer_cnt);
	list_buffer_cnt++;
	degree[edge_from]++;
}
inline void GRAPH::remove_lr(vertex &removed_vertex)
{
	LBL(LBR(removed_vertex)) = LBL(removed_vertex);
	LBR(LBL(removed_vertex)) = LBR(removed_vertex);
}
inline void GRAPH::remove_ud(vertex &removed_vertex)
{
	LBD(LBU(removed_vertex)) = LBD(removed_vertex);
	LBU(LBD(removed_vertex)) = LBU(removed_vertex);
}
void GRAPH::remove(vertex removed_vertex)
{
	is_available[removed_vertex] = false;
	local_n--;
	size_t now = removed_vertex;
	remove_lr(now);
	now = LBD(now);
	while (now != removed_vertex)
	{
		remove_lr(now);
		degree[LBLB(now)]--;
		all_neighbors_weight[LBLB(now)] -= weights[removed_vertex];
		add_to_reduction_queue(LBLB(now), 0);
		now = LBD(now);
	}

	removed_vertex += max_n;
	now = removed_vertex;
	remove_ud(now);
	now = LBR(now);
	while (now != removed_vertex)
	{
		remove_ud(now);
		now = LBR(now);
	}
}
inline void GRAPH::recover_lr(vertex &removed_vertex)
{
	LBL(LBR(removed_vertex)) = removed_vertex;
	LBR(LBL(removed_vertex)) = removed_vertex;
}
inline void GRAPH::recover_ud(vertex &removed_vertex)
{
	LBD(LBU(removed_vertex)) = removed_vertex;
	LBU(LBD(removed_vertex)) = removed_vertex;
}
void GRAPH::recover(vertex removed_vertex)
{
	is_available[removed_vertex] = true;
	//all_neighbors_weight[removed_vertex] = 0;
	local_n++;
	size_t now = removed_vertex;
	recover_lr(now);
	now = LBD(now);
	while (now != removed_vertex)
	{
		recover_lr(now);
		degree[LBLB(now)]++;
		all_neighbors_weight[LBLB(now)] += weights[removed_vertex];
		//all_neighbors_weight[removed_vertex] += weights[LBLB(now)];
		now = LBD(now);
	}
	removed_vertex += max_n;
	now = removed_vertex;
	recover_ud(now);
	now = LBR(now);
	while (now != removed_vertex)
	{
		recover_ud(now);
		now = LBR(now);
	}
}

void GRAPH::run_branch(weight_node lower_weight)
{
	weight_node local_weight = 0;

	is_weight = max(lower_weight - 1, 0);

	pre_n = local_n;
	size_t old_local_weight;

	if (!is_called_from_split_direct)
	{
		while (true)
		{
			old_local_weight = local_weight;

			for (size_t i = LBR(0); i != 0; i = LBR(i))
			{
				reduction_queue[0].push(i);
			}

			local_weight += apply_all_reductions();

			if (old_local_weight == local_weight)
			{
				break;
			}
			if (!is_fisrt_runed)
			{
				break;
			}
		}
	}

	branch(local_weight);
}

void GRAPH::branch(weight_node local_weight)
{
	reduction_time = std::chrono::system_clock::now() - start_time;
	remaining_vertices = local_n;
	reduction_offset_all = local_weight;

	if (is_fisrt_runed)
	{
		cout << reduction_time.count() << " " << local_n << " " << local_weight << endl;
	}

	update_solution(local_weight);

	if (local_n == 0)
	{
		return;
	}

	GRAPH *component_alg = &GRAPH_BUFFER[depth + 1];

	if (!is_called_from_split_direct && local_n >= SPLIT_LIMIT)
	{
		vector<vector<vertex>> components;
		split_components(components);

		if (components.size() > 1)
		{
			vector<weight_node> component_upperbound(components.size(), 0);
			vector<size_t> component_idx(components.size(), 0);
			vector<double> avg_degree(components.size(), 0);
			weight_node total_clique_cover_weight = 0;
			for (size_t i = 0; i < components.size(); i++)
			{
				component_upperbound[i] = computing_upper_bound(components[i], INF);
				total_clique_cover_weight += component_upperbound[i];
				component_idx[i] = i;
				for (auto u : components[i])
				{
					avg_degree[i] += degree[u];
				}
				avg_degree[i] /= components[i].size();
			}

			//Here are some strategies to select connected components to run first.
			//sort(component_idx.begin(), component_idx.end(), [&components, &component_upperbound](const size_t lhs, const size_t rhs) { return component_upperbound[lhs] < component_upperbound[rhs] || (component_upperbound[lhs] == component_upperbound[rhs] && components[lhs].size() < components[rhs].size()); });
			//sort(component_idx.begin(), component_idx.end(), [&components, &component_upperbound](const size_t lhs, const size_t rhs) { return components[lhs].size() < components[rhs].size() || (components[lhs].size() == components[rhs].size() && component_upperbound[lhs] < component_upperbound[rhs]); });
			//sort(component_idx.begin(), component_idx.end(), [&components, &component_upperbound](const size_t lhs, const size_t rhs) { return components[lhs].size() > components[rhs].size() || (components[lhs].size() == components[rhs].size() && component_upperbound[lhs] > component_upperbound[rhs]); });
			sort(component_idx.begin(), component_idx.end(), [&](const size_t lhs, const size_t rhs) { return avg_degree[lhs] < avg_degree[rhs]; });

			bool is_better = true;
#ifdef DEBUG
			IS_STATUS.reset();
#endif
			size_t remaining_n = local_n;

			for (size_t i = 0; i < components.size(); i++)
			{
				if (local_weight + total_clique_cover_weight <= is_weight)
				{
					is_better = false;
					break;
				}

				component_alg->reset(list_buffer + list_buffer_cnt, weights + cnt_n, all_neighbors_weight + cnt_n, degree + cnt_n);
				local_weight += branch_in_component(component_alg, components[component_idx[i]], is_weight - (local_weight + total_clique_cover_weight - component_upperbound[component_idx[i]]), true);
				total_clique_cover_weight -= component_upperbound[component_idx[i]];
			}

			if (!is_better)
			{
				return;
			}
			if (local_weight > is_weight)
			{
				update_solution(local_weight);
#ifdef DEBUG
				for (auto component : components)
				{
					for (vertex vc : component)
					{
						if (IS_STATUS[vc])
						{
							best_solution.push(vertex_status(vc, VS::VS_INCLUDED, 0));
						}
					}
				}
#endif
			}
			return;
		}
		else if (local_n <= global_n / 2)
		{
			component_alg->reset(list_buffer + list_buffer_cnt, weights + cnt_n, all_neighbors_weight + cnt_n, degree + cnt_n);
			local_weight += branch_in_component(component_alg, components[0], is_weight - local_weight, true);
			if (local_weight > is_weight)
			{
				update_solution(local_weight);
#ifdef DEBUG
				for (vertex vc : components[0])
				{
					if (IS_STATUS[vc])
					{
						best_solution.push(vertex_status(vc, VS::VS_INCLUDED, 0));
					}
				}
#endif
			}
			return;
		}
	}

	if (local_n > LOCAL_SEARCH_LIMIT || is_running_time_out())
		local_search(local_weight);
	if (is_running_time_out())
		return;

	vertex v = find_max_available_vertex();
	vertex first_branch = v;
	pre_n = local_n;
	while (true)
	{
		if (is_running_time_out())
		{
			break;
		}

		if (v == 0)
		{
			update_solution(local_weight);

			v = branching_stack.top().idx;
			continue;
		}

		if (branching_stack.empty() || branching_stack.top().idx != v)
		{
			vector<vertex> S_v;
			if (get_unconfined_set_exactly(v, S_v))
			{
				exclude_vertex(v);
				branching_stack.push(branching_node(v, local_weight, IS::EXCLUDED, pre_n));
			}
			else
			{
				branching_stack.push(branching_node(v, local_weight, IS::INCLUDED, pre_n));
				for (auto u : S_v)
				{
					local_weight += include_vertex(u);
				}
			}
		}
		else
		{
			auto top_branch = branching_stack.top();
			branching_stack.pop();
			local_weight = top_branch.local_weight;
			recover_reductions();
			if (top_branch.status == IS::INCLUDED)
			{
				top_branch.status = IS::EXCLUDED;
				exclude_vertex(v);

				branching_stack.push(top_branch);
				pre_n = top_branch.pre_n;
			}
			else
			{
				if (v == first_branch)
				{
					break;
				}
				v = branching_stack.top().idx;
				continue;
			}
		}

		local_weight += apply_all_reductions();

		if (local_n >= SPLIT_LIMIT && 9 * pre_n >= 10 * local_n) //When the vertices decrease 10%, run global reductions and branch in components independently.
		{
			pre_n = local_n;
			vector<vector<vertex>> components;
			split_components(components);

			if (components.size() > 1)
			{
				vector<weight_node> component_upperbound(components.size(), 0);
				vector<size_t> component_idx(components.size(), 0);
				weight_node total_clique_cover_weight = 0;
				for (size_t i = 0; i < components.size(); i++)
				{
					component_upperbound[i] = computing_upper_bound(components[i], INF);
					total_clique_cover_weight += component_upperbound[i];
					component_idx[i] = i;
				}

				sort(component_idx.begin(), component_idx.end(), [&components, &component_upperbound](const size_t lhs, const size_t rhs) { return components[lhs].size() < components[rhs].size() || (components[lhs].size() == components[rhs].size() && component_upperbound[lhs] < component_upperbound[rhs]); });

				bool is_better = true;
#ifdef DEBUG
				IS_STATUS.reset();
#endif
				for (size_t i = 0; i < components.size(); i++)
				{
					if (local_weight + total_clique_cover_weight <= is_weight)
					{
						is_better = false;
						break;
					}

					component_alg->reset(list_buffer + list_buffer_cnt, weights + cnt_n, all_neighbors_weight + cnt_n, degree + cnt_n);
					local_weight += branch_in_component(component_alg, components[component_idx[i]], is_weight - (local_weight + total_clique_cover_weight - component_upperbound[component_idx[i]]), false);
					total_clique_cover_weight -= component_upperbound[component_idx[i]];
				}
				modified_stack.push(modified_node(0, MS::TOKEN, 0));
#ifdef DEBUG
				local_solution.push(vertex_status(0, VS::VS_TOKEN, 0));
#endif
				if (!is_better)
				{
					recover_reductions();
					continue;
				}

				if (local_weight > is_weight)
				{
					update_solution(local_weight);
#ifdef DEBUG
					best_solution = local_solution;

					for (auto component : components)
					{
						for (vertex vc : component)
						{
							if (IS_STATUS[vc])
							{
								best_solution.push(vertex_status(vc, VS::VS_INCLUDED, 0));
							}
						}
					}
#endif
				}
				recover_reductions();
				continue;
			}
			else if (local_n <= global_n / 2)
			{
				component_alg->reset(list_buffer + list_buffer_cnt, weights + cnt_n, all_neighbors_weight + cnt_n, degree + cnt_n);
				local_weight += branch_in_component(component_alg, components[0], is_weight - local_weight, true);

				modified_stack.push(modified_node(0, MS::TOKEN, 0));
#ifdef DEBUG
				local_solution.push(vertex_status(0, VS::VS_TOKEN, 0));
#endif

				if (local_weight > is_weight)
				{
#ifdef DEBUG
					for (vertex vc : components[0])
					{
						if (IS_STATUS[vc])
						{
							local_solution.push(vertex_status(vc, VS::VS_INCLUDED, 0));
						}
					}
#endif
					update_solution(local_weight);
				}
				recover_reductions();
				continue;
			}
		}

		if (9 * pre_n < 10 * local_n)
		{
			v = find_max_available_vertex();
		}
		else //We use the clique cover to prune
		{
			vector<vertex> clique_cover;
			get_available_vertices(clique_cover);
			if (local_weight + computing_upper_bound(clique_cover, is_weight - local_weight) > is_weight)
			{
				v = find_max_available_vertex();
			}
		}
	}

	return;
}
void GRAPH::reset(list_node *list_buffer_start, weight_node *weight_buffer_start, weight_node *neighbors_weight_buffer_start, size_t *degree_buffer_start)
{
	list_buffer = list_buffer_start;
	weights = weight_buffer_start;
	all_neighbors_weight = neighbors_weight_buffer_start;
	degree = degree_buffer_start;
	list_buffer_cnt = 0;
	is_available.reset();
	is_weight = 0;
	modified_stack = stack<modified_node>();
	branching_stack = stack<branching_node>();
#ifdef DEBUG
	IS_STATUS.reset();
#endif
	ls_lower_weight = 0;
}

vertex GRAPH::find_max_available_vertex()
{
	vertex max_vertex = LBR(0);

	for (size_t i = LBR(0); i != 0; i = LBR(i))
	{
		if (degree[i] > degree[max_vertex] || (degree[i] == degree[max_vertex] && weights[i] > weights[max_vertex]))
		{
			max_vertex = i;
		}
		// if (all_neighbors_weight[i] - weights[i] > all_neighbors_weight[max_vertex] - weights[max_vertex])
		// {
		// 	max_vertex = i;
		// }
	}

	return max_vertex;
}

weight_node GRAPH::apply_all_reductions()
{
	auto s_tttt = std::chrono::system_clock::now();

	weight_node reduction_offset = 0;
	bool reduction_restart_flag = true;

	while (reduction_restart_flag)
	{
		reduction_restart_flag = false;
		for (size_t i = 0; i < 7; i++)
		{
			size_t old_n = local_n;
			weight_node old_offset = reduction_offset;
			auto re_t = std::chrono::system_clock::now();
			reduction_offset += (this->*(reductions[i]))();
			std::chrono::duration<float> cost_time = std::chrono::system_clock::now() - re_t;
			reductions_tims_cnt[i] += cost_time.count();
			reductions_cnt[i] += old_n - local_n;
			reductions_weight_cnt[i] += reduction_offset - old_offset;

			if (old_n != local_n)
			{
				reduction_restart_flag = true;
				break;
			}
		}
	}

	std::chrono::duration<float> cost_time = std::chrono::system_clock::now() - s_tttt;
	all_reduction_cost += cost_time.count();

	modified_stack.push(modified_node(0, MS::TOKEN, 0));
#ifdef DEBUG
	local_solution.push(vertex_status(0, VS::VS_TOKEN, 0));
#endif
	return reduction_offset;
}

void GRAPH::recover_reductions()
{
	if (!modified_stack.empty() && modified_stack.top().status == MS::TOKEN)
		modified_stack.pop();

	while (!modified_stack.empty() && modified_stack.top().status != MS::TOKEN)
	{
		auto top_mod = modified_stack.top();
		modified_stack.pop();
		if (top_mod.status == MS::REMOVED)
		{
			recover(top_mod.idx);
		}
		else if (top_mod.status == MS::WEIGHT_MODIFIED)
		{
			For(i, top_mod.idx)
			{
				all_neighbors_weight[LBLB(i)] += top_mod.weight - weights[top_mod.idx];
			}
			weights[top_mod.idx] = top_mod.weight;
			continue;
		}
		else if (top_mod.status == MS::ADDED)
		{
			remove(top_mod.idx);
			cnt_n = top_mod.idx;
			list_buffer_cnt = top_mod.lbcnt;
			continue;
		}
	}

#ifdef DEBUG
	if (!local_solution.empty() && local_solution.top().vs == VS::VS_TOKEN)
		local_solution.pop();

	while (!local_solution.empty() && local_solution.top().vs != VS::VS_TOKEN)
	{
		local_solution.pop();
	}
#endif
}

void GRAPH::get_available_vertices(vector<vertex> &vertex_order)
{
	for (size_t i = LBR(0); i != 0; i = LBR(i))
	{
		vertex_order.push_back(i);
	}
}

weight_node GRAPH::computing_upper_bound(vector<vertex> &vertex_order, weight_node target_weight)
{
	if (vertex_order.size() == 0)
	{
		return 0;
	}
	weight_node upper_weight = 0;
	SET &is_colored = set_buffer0;
	is_colored.clear();
	sort(vertex_order.begin(), vertex_order.end(), [&](const vertex lhs, const vertex rhs) {
		return weights[lhs] > weights[rhs] || (weights[lhs] == weights[rhs] && degree[lhs] > degree[rhs]);
		//return degree[lhs] > degree[rhs] || (degree[lhs] > degree[rhs] && weights[lhs] > weights[rhs]);
	});

	vector<size_t> clique_cnt(vertex_order.size() + 1, 0);
	vector<weight_node> clique_weight(vertex_order.size() + 1, 0);
	size_t clique_colors_cnt = 0;
	vector<size_t> neighbors_cnt;
	bool is_alone;
	for (size_t k = 0; k < vertex_order.size(); k++)
	{
		auto v = vertex_order[k];
		neighbors_cnt.assign(clique_colors_cnt + 1, 0);

		For(i, v)
		{
			if (is_colored.get(LBLB(i)))
				neighbors_cnt[REMAP[LBLB(i)]]++;
		}
		is_alone = true;
		for (size_t i = 1; i <= clique_colors_cnt; i++)
		{
			if (neighbors_cnt[i] == clique_cnt[i])
			{
				REMAP[v] = i;
				is_colored.add(v);
				clique_cnt[i]++;
				if (weights[v] > clique_weight[i])
					upper_weight += weights[v] - clique_weight[i];
				is_alone = false;
				break;
			}
		}
		if (is_alone)
		{
			clique_colors_cnt++;
			clique_cnt[clique_colors_cnt] = 1;
			REMAP[v] = clique_colors_cnt;
			is_colored.add(v);
			upper_weight += weights[v];
			clique_weight[clique_colors_cnt] = weights[v];

			if (upper_weight > target_weight)
			{
				return upper_weight;
			}
		}
	}

	if (clique_colors_cnt == 1)
		is_clique_neighborhood = true;
	if (clique_colors_cnt == vertex_order.size())
		is_independent_neighborhood = true;
	return upper_weight;
}

inline void GRAPH::add_to_reduction_queue(vertex idx, size_t queue_idx)
{
	if (!is_in_reduction_queue[queue_idx][idx])
	{
		reduction_queue[queue_idx].push(idx);
		is_in_reduction_queue[queue_idx][idx] = true;
	}
}
bool GRAPH::is_closed(vertex u, vertex v)
{
	if (degree[u] > degree[v])
	{
		swap(u, v);
	}

	For(i, u)
	{
		if (LBLB(i) == v)
		{
			return true;
		}
	}
	return false;
}

void GRAPH::modify_weight(vertex v, weight_node obj_weight)
{
	modified_stack.push(modified_node(v, MS::WEIGHT_MODIFIED, weights[v]));
	add_to_reduction_queue(v, reductions_type::neighborhood);

	For(i, v)
	{
		all_neighbors_weight[LBLB(i)] += obj_weight - weights[v];
		add_to_reduction_queue(LBLB(i), reductions_type::neighborhood);
	}
	weights[v] = obj_weight;
}

weight_node GRAPH::include_vertex(vertex v) //include v into the independent set
{
	modified_stack.push(modified_node(v, MS::REMOVED, weights[v]));
#ifdef DEBUG
	local_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
	remove(v);

	For(i, v)
	{
		exclude_vertex(LBLB(i));
	}

	return weights[v];
}

inline void GRAPH::exclude_vertex(vertex v) //not exclude vertex, just remove
{
	modified_stack.push(modified_node(v, MS::REMOVED, weights[v]));
	remove(v);
}

weight_node GRAPH::degree01_process(vertex v)
{
	if (weights[v] == 0)
	{
		exclude_vertex(v);
	}
	else if (degree[v] == 1)
	{
		auto u = LBLB(LBD(v));
		if (weights[v] >= weights[u])
		{
			return include_vertex(v);
		}
		else
		{
			exclude_vertex(v);
			modify_weight(u, weights[u] - weights[v]);

#ifdef DEBUG
			local_solution.push(vertex_status(v, VS::VS_NOTIF, u));
			local_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
			return weights[v];
		}
	}
	return 0;
}

weight_node GRAPH::fold_vertices(vertex v)
{
	exclude_vertex(v);
	SET &neighbors = set_buffer0;
	neighbors.clear();
	neighbors.add(v);

#ifdef DEBUG
	local_solution.push(vertex_status(v, VS::VS_IFNOT, cnt_n));
#endif

	For(i, v)
	{
		vertex u = LBLB(i);
		neighbors.add(u);
		exclude_vertex(u);

#ifdef DEBUG
		local_solution.push(vertex_status(u, VS::VS_IF, cnt_n));
#endif
	}

	insert_vertex(cnt_n);
	weights[cnt_n] = all_neighbors_weight[v] - weights[v];
	modified_stack.push(modified_node(cnt_n, MS::ADDED, weights[cnt_n], list_buffer_cnt));

	For(i, v)
	{
		vertex j = LBLB(i);

		For(k, j)
		{
			vertex u = LBLB(k);

			if (!neighbors.get(u))
			{
				add_to_reduction_queue(u, reductions_type::neighborhood);
				insert_edge(cnt_n, u);
				insert_edge(u, cnt_n);
				neighbors.add(u);
			}
		}
	}
	add_to_reduction_queue(cnt_n, reductions_type::neighborhood);
	cnt_n++;
	return weights[v];
}
weight_node GRAPH::degree2_process(vertex v)
{
	auto offset = general_degree2_cases(v); //some genearal cases of degree 2
	if (!is_available[v])
		return offset;

	vertex v1, v2, v3, v4, v5;
	v3 = v;
	v2 = LBLB(LBU(v3));
	v4 = LBLB(LBD(v3));

	size_t cnt_flag = 0;
	offset = 0;

	while (cnt_flag++ < 2) // a special case of degree 2, newly added; when combined with folding, it has a powerful performence
	{
		if (weights[v2] >= weights[v3] && weights[v3] > weights[v4])
		{
			offset = weights[v3] - weights[v4];
			modify_weight(v2, weights[v2] - offset);
			modify_weight(v3, weights[v4]);
#ifdef DEBUG
			// exactly one of v2 and v3 is in the mwis, assume that v3 is in mwis by default, then v2->!v3 and v3->!v4
			local_solution.push(vertex_status(v4, VS::VS_NOTIF, v3));
			local_solution.push(vertex_status(v3, VS::VS_NOTIF, v2));
			local_solution.push(vertex_status(v3, VS::VS_INCLUDED, 0));
#endif
			break;
		}
		swap(v2, v4);
	}

	cnt_flag = 0;
	while (cnt_flag++ < 2) //4-cycle containing 2 degree-2 vertices.
	{
		if (degree[v4] == 2 && weights[v2] >= weights[v3] && weights[v3] >= weights[v4])
		{
			For(i, v4)
			{
				if (LBLB(i) != v3)
				{
					v5 = LBLB(i);
					break;
				}
			}

			if (!is_closed(v2, v5))
			{
				continue;
			}
			return offset + fold_vertices(v3);
		}
		swap(v2, v4);
	}

	return offset;
}

weight_node GRAPH::degree2_reduction()
{
	weight_node reduction_offset = 0;
	size_t old_n = local_n;

	while (!reduction_queue[reductions_type::degree2].empty())
	{
		auto v = reduction_queue[reductions_type::degree2].front();
		reduction_queue[reductions_type::degree2].pop();
		is_in_reduction_queue[reductions_type::degree2][v] = false;
		if (!is_available[v])
		{
			continue;
		}

		if (weights[v] >= all_neighbors_weight[v])
		{
			reduction_offset += include_vertex(v);
			continue;
		}

		if (degree[v] < 2 || weights[v] == 0)
		{
			reduction_offset += degree01_process(v);
			continue;
		}
		else if (degree[v] == 2)
		{
			reduction_offset += degree2_process(v);
			continue;
		}
		add_to_reduction_queue(v, reductions_type::twin);
	}
	return reduction_offset;
}

weight_node GRAPH::neighborhood_reduction()
{
	weight_node reduction_offset = 0;
	while (!reduction_queue[reductions_type::neighborhood].empty())
	{
		auto v = reduction_queue[reductions_type::neighborhood].front();
		reduction_queue[reductions_type::neighborhood].pop();
		is_in_reduction_queue[reductions_type::neighborhood][v] = false;
		if (!is_available[v])
		{
			continue;
		}

		weight_node neighbor_weight = 0;

		if (weights[v] >= all_neighbors_weight[v])
		{
			reduction_offset += include_vertex(v);
			continue;
		}
		if (degree[v] < 2 || weights[v] == 0)
		{
			reduction_offset += degree01_process(v);
			continue;
		}

		add_to_reduction_queue(v, reductions_type::degree2);
	}
	return reduction_offset;
}

void GRAPH::merge_simultaneous_set(vector<vertex> &SimS)
{
	weight_node total_weight = 0;
	for (auto v : SimS)
	{
		total_weight += weights[v];
		exclude_vertex(v);
	}
	insert_vertex(cnt_n);
	weights[cnt_n] = total_weight;
	modified_stack.push(modified_node(cnt_n, MS::ADDED, weights[cnt_n], list_buffer_cnt));
	SET &neighbors = set_buffer0;
	neighbors.clear();

	for (auto v : SimS)
	{
#ifdef DEBUG
		local_solution.push(vertex_status(v, VS::VS_IF, cnt_n));
#endif
		For(k, v)
		{
			vertex u = LBLB(k);

			if (!neighbors.get(u))
			{
				add_to_reduction_queue(u, reductions_type::neighborhood);
				insert_edge(cnt_n, u);
				insert_edge(u, cnt_n);
				neighbors.add(u);
			}
		}
	}
	add_to_reduction_queue(cnt_n, reductions_type::neighborhood);
	cnt_n++;
}

bool GRAPH::get_unconfined_set_exactly(vertex &v, vector<vertex> &S_v) // Check v is unconfined and get confined set S_v
{
	SET &is_visited = set_buffer3;
	is_visited.clear();
	queue<vertex> Q;
	SET &is_in_Q = set_buffer2;
	is_in_Q.clear();
	SET &set_S_v = set_buffer0;
	set_S_v.clear();
	set_S_v.add(v);
	S_v.push_back(v);
	SET &set_neighbos_S_v = set_buffer1;
	set_neighbos_S_v.clear();

	bool is_unconfined = false;

	For(i, v)
	{
		Q.push(LBLB(i));
		is_in_Q.add(LBLB(i));
		set_neighbos_S_v.add(LBLB(i));
		is_visited.add(LBLB(i));
		cap_weight[LBLB(i)] = weights[v];
	}

	while (!Q.empty())
	{
		vertex u = Q.front();
		Q.pop();
		is_in_Q.remove(u);

		if (weights[u] < cap_weight[u] || degree[u] > 16 + set_S_v.size() + set_neighbos_S_v.size()) //Here we introduce a parameter 16 to reduce the search tree of unconfined set, it is not the best one
		{
			continue;
		}

		if (weights[u] >= all_neighbors_weight[u])
		{
			is_unconfined = true;
			break;
		}

		vector<vertex> uncovered_vertices;
		weight_node max_neighbor_weight = 0;

		weight_node all_weight = 0;
		weight_node upper_weight = all_neighbors_weight[u];
		bool is_large = false;
		bool is_possible_heavy = true;
		For(i, u)
		{
			if (!set_S_v.get(LBLB(i)))
			{
				if (!set_neighbos_S_v.get(LBLB(i)))
				{
					uncovered_vertices.push_back(LBLB(i));
					max_neighbor_weight = max(max_neighbor_weight, weights[LBLB(i)]);
					all_weight += weights[LBLB(i)];
					if (uncovered_vertices.size() > 1 && max_neighbor_weight + cap_weight[u] > weights[u])
					{
						break;
					}
					if (uncovered_vertices.size() > 16)
					{
						is_large = true;
						break;
					}
					if (uncovered_vertices.size() > 1 && max_neighbor_weight + cap_weight[u] > weights[u])
					{
						is_possible_heavy = false;
						break;
					}
				}
				else
				{
					upper_weight -= weights[LBLB(i)];
					if (upper_weight <= weights[u])
					{
						break;
					}
				}
			}
		}

		if (upper_weight <= weights[u])
		{
			is_unconfined = true;
			break;
		}

		if (!is_possible_heavy)
		{
			continue;
		}

		if (is_large)
		{
			continue;
		}

		if (uncovered_vertices.size() == 1)
		{
			S_v.push_back(uncovered_vertices[0]);
			set_S_v.add(uncovered_vertices[0]);

			For(i, uncovered_vertices[0])
			{
				set_neighbos_S_v.add(LBLB(i));
				if (!is_in_Q.get(LBLB(i)))
				{
					Q.push(LBLB(i));
					is_in_Q.add(LBLB(i));
				}

				if (!is_visited.get(LBLB(i)))
				{
					cap_weight[LBLB(i)] = 0;
					is_visited.add(LBLB(i));
				}
				cap_weight[LBLB(i)] = cap_weight[LBLB(i)] + weights[uncovered_vertices[0]];
			}
			continue;
		}

		if (max_neighbor_weight + cap_weight[u] > weights[u])
		{
			continue;
		}

		if (uncovered_vertices.size() <= MAX_SMALL_GRAPH_SIZE)
		{
			vector<vertex> &vertices = uncovered_vertices;
			SET &is_mapped = set_buffer4;
			is_mapped.clear();
			size_t idx_cnt = 0;
			sort(vertices.begin(), vertices.end(), [&](const vertex &lhs, const vertex &rhs) { return weights[lhs] < weights[rhs]; });
			for (; idx_cnt < vertices.size(); idx_cnt++)
			{
				REMAP[vertices[idx_cnt]] = idx_cnt;
				is_mapped.add(vertices[idx_cnt]);
			}
			small_graph.build(idx_cnt, weights[u] - cap_weight[u]);

			for (vertex u : vertices)
			{
				small_graph.assign_weight(REMAP[u], weights[u]);
				For(i, u)
				{
					if (is_mapped.get(LBLB(i)))
					{
						small_graph.add_edge(REMAP[u], REMAP[LBLB(i)]);
					}
				}
			}
			if (small_graph.is_independent_neighborhood())
			{
				if (cap_weight[u] + all_weight - weights[vertices[0]] <= weights[u])
				{
					for (vertex uc : uncovered_vertices)
					{
						S_v.push_back(uc);
						set_S_v.add(uc);

						For(i, uc)
						{
							set_neighbos_S_v.add(LBLB(i));
							if (!is_in_Q.get(LBLB(i)))
							{
								Q.push(LBLB(i));
								is_in_Q.add(LBLB(i));
							}

							if (!is_visited.get(LBLB(i)))
							{
								cap_weight[LBLB(i)] = 0;
								is_visited.add(LBLB(i));
							}
							cap_weight[LBLB(i)] = cap_weight[LBLB(i)] + weights[uc];
						}
					}
				}
				continue;
			}
			else if (small_graph.is_clique_neighborhood() || small_graph.max_is() + cap_weight[u] <= weights[u])
			{
				is_unconfined = true;
				break;
			}
		}
	}

	return is_unconfined;
}
weight_node GRAPH::unconfined_reduction()
{

	weight_node reduction_offset = 0;
	size_t old_n = local_n;

	while (!reduction_queue[reductions_type::unconfined].empty())
	{
		auto v = reduction_queue[reductions_type::unconfined].front();
		reduction_queue[reductions_type::unconfined].pop();

		if (!is_available[v] || !is_in_reduction_queue[reductions_type::unconfined][v])
		{
			is_in_reduction_queue[reductions_type::unconfined][v] = false;
			continue;
		}

		is_in_reduction_queue[reductions_type::unconfined][v] = false;

		if (weights[v] >= all_neighbors_weight[v])
		{
			reduction_offset += include_vertex(v);
			continue;
		}

		if (degree[v] < 2 || weights[v] == 0)
		{
			reduction_offset += degree01_process(v);
			continue;
		}

		vector<vertex> S_v;
		if (get_unconfined_set_exactly(v, S_v))
		{
			exclude_vertex(v);

			for (auto u : S_v)
			{
				add_to_reduction_queue(u, reductions_type::unconfined);
			}

			continue;
		}
		else
		{
			for (auto u : S_v)
			{
				is_in_reduction_queue[reductions_type::unconfined][u] = false;
			}

			//Here we can reduce further vertices by merging simultaneous set that is getting from confining set.

			// 	if (S_v.size() > 1)
			// 	{

			// 		vector<vertex> simu_set;
			// 		SET &is_checked = set_buffer0;
			// 		is_checked.clear();

			// 		simu_set.push_back(v);
			// 		is_checked.add(v);
			// 		for (auto u : S_v)
			// 		{
			// 			if (!is_checked.get(v))
			// 			{
			// 				vector<vertex> S_u;
			// 				get_unconfined_set_exactly(u, S_u);
			// 				if (S_u.size() == S_v.size())
			// 				{
			// 					simu_set.push_back(u);
			// 				}
			// 				else
			// 				{
			// 					for (auto w : S_u)
			// 						is_checked.add(w);
			// 				}
			// 				is_checked.add(u);
			// 			}
			// 			is_in_reduction_queue[reductions_type::unconfined][u] = false;
			// 		}

			// 		if (simu_set.size() > 1)
			// 		{
			// 			merge_simultaneous_set(simu_set);

			// 		}
			// 	}
		}
	}
	return reduction_offset;
}

weight_node GRAPH::general_degree2_cases(vertex v)
{
	vector<vertex> neighbors;
	weight_node neighbor_weight = all_neighbors_weight[v];
	weight_node min_neighbor_weight = INF;
	For(i, v)
	{
		neighbors.push_back(LBLB(i));
		min_neighbor_weight = min(min_neighbor_weight, weights[LBLB(i)]);
	}

	is_clique_neighborhood = false;
	weight_node clique_cover_weight = weights[neighbors[0]] + weights[neighbors[1]];

	if (is_closed(neighbors[0], neighbors[1]))
	{
		clique_cover_weight = max(weights[neighbors[0]], weights[neighbors[1]]);
		is_clique_neighborhood = true;
	}

	if (weights[v] >= clique_cover_weight)
	{
		return include_vertex(v);
	}
	else if (!is_clique_neighborhood && weights[v] >= neighbor_weight - min_neighbor_weight)
	{
		return fold_vertices(v);
	}
	else if (is_clique_neighborhood)
	{
		vertex max_isolated = v, max_neighbor = v;
		for (vertex u : neighbors)
		{
			if (degree[u] == degree[v])
			{
				if (weights[u] > weights[max_isolated])
				{
					max_isolated = u;
				}
			}
			else if (weights[u] > weights[max_neighbor])
			{
				max_neighbor = u;
			}
		}

		v = max_isolated;

		if (weights[v] >= weights[max_neighbor])
		{
			return include_vertex(v);
		}

		exclude_vertex(v);

		For(i, v)
		{
			auto u = LBLB(i);

			if (weights[v] >= weights[u])
			{
				exclude_vertex(u);
			}
			else
			{
				modify_weight(u, weights[u] - weights[v]);
#ifdef DEBUG
				local_solution.push(vertex_status(v, VS::VS_NOTIF, u));
#endif
			}
		}
#ifdef DEBUG
		local_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
		return weights[v];
	}
	add_to_reduction_queue(v, reductions_type::unconfined);
	return 0;
}

weight_node GRAPH::twin_reduction()
{
	weight_node reduction_offset = 0;
	size_t old_n = local_n;

	while (!reduction_queue[reductions_type::twin].empty())
	{
		auto v = reduction_queue[reductions_type::twin].front();
		reduction_queue[reductions_type::twin].pop();
		is_in_reduction_queue[reductions_type::twin][v] = false;

		if (!is_available[v])
			continue;
		add_to_reduction_queue(v, reductions_type::heavy_set);

		if (weights[v] >= all_neighbors_weight[v])
		{
			reduction_offset += include_vertex(v);
			continue;
		}

		if (degree[v] < 2 || weights[v] == 0)
		{
			reduction_offset += degree01_process(v);
			continue;
		}

		if (degree[v])
		{
			vector<vertex> twins;
			twins.push_back(v);
			weight_node twins_weight = weights[v];

			SET &is_checked_twin = set_buffer0;
			SET &neighbors = set_buffer1;
			is_checked_twin.clear();
			neighbors.clear();
			is_checked_twin.add(v);
			vertex min_neighbor = LBLB(LBD(v));
			For(i, v)
			{
				is_checked_twin.add(LBLB(i));
				neighbors.add(LBLB(i));

				if (degree[min_neighbor] > degree[LBLB(i)])
					min_neighbor = LBLB(i);
			}

			For(i, min_neighbor)
			{
				auto u = LBLB(i);
				if (!is_checked_twin.get(u) && degree[u] == degree[v])
				{
					bool is_twin = true;

					For(l, u)
					{
						if (!neighbors.get(LBLB(l)))
						{
							is_twin = false;
							break;
						}
					}
					if (is_twin)
					{
						twins.push_back(u);
						twins_weight += weights[u];
					}
				}
				is_checked_twin.add(u);
			}

			if (twins.size() == 1)
			{
				continue;
			}

			if (twins_weight >= all_neighbors_weight[v])
			{
				reduction_offset += twins_weight;

				for (auto u : twins)
				{
#ifdef DEBUG
					local_solution.push(vertex_status(u, VS::VS_INCLUDED, 0));
#endif
					exclude_vertex(u);
				}
				For(i, v)
				{
					auto u = LBLB(i);
					exclude_vertex(u);
				}
				continue;
			}
			else
			{
				for (size_t i = 1; i < twins.size(); i++)
				{
					exclude_vertex(twins[i]);
#ifdef DEBUG
					local_solution.push(vertex_status(twins[i], VS::VS_IF, v));
#endif
				}
				modify_weight(v, twins_weight);
				continue;
			}
		}
	}
	return reduction_offset;
}

weight_node GRAPH::critical_set_reduction()
{
	if ((branching_stack.size() != 0) && (local_n < 100 || 9 * pre_n < 10 * local_n)) // && (branching_stack.size() % 3 != 0)
	{
		return 0;
	}

	weight_node reduction_offset = 0;
	vector<vertex> vertices;
	get_available_vertices(vertices);

	if (vertices.size() <= 2)
	{
		return 0;
	}

	auto num = vertices.size();
	vertex source = 0;
	vertex sink = 2 * num + 1;

	for (size_t i = 0; i < num; i++)
	{
		REMAP[vertices[i]] = i + 1;
	}

	flow_graph->init(num * 2 + 2);

	for (size_t i = 1; i <= num; i++)
	{
		vertex cor = vertices[i - 1];
		flow_graph->addedge(source, i, weights[cor]);
	}
	for (size_t i = 1; i <= num; i++)
	{
		vertex cor = vertices[i - 1];
		For(j, cor)
		{
			flow_graph->addedge(i, REMAP[LBLB(j)] + num, weights[cor]);
		}

		flow_graph->addedge(i + num, sink, weights[cor]);
	}

	flow_graph->Max_flow();
	vector<vertex> cst;
	SET &cst_flag = set_buffer0;
	cst_flag.clear();
	for (size_t i = 0; i < num; i++)
	{
		if (flow_graph->flow[i * 2] > 0)
		{
			cst.push_back(vertices[i]);
			cst_flag.add(vertices[i]);
		}
	}

	for (auto v : cst)
	{
		bool is_isolated = true;

		For(i, v)
		{
			if (cst_flag.get(LBLB(i)))
			{
				is_isolated = false;
				break;
			}
		}

		if (is_isolated)
		{
			reduction_offset += include_vertex(v);
		}
	}

	return reduction_offset;
}

weight_node GRAPH::heavy_set_reduction()
{
	weight_node reduction_offset = 0;
	while (!reduction_queue[reductions_type::heavy_set].empty())
	{
		auto v = reduction_queue[reductions_type::heavy_set].front();
		reduction_queue[reductions_type::heavy_set].pop();
		is_in_reduction_queue[reductions_type::heavy_set][v] = false;

		if (!is_available[v])
			continue;

		add_to_reduction_queue(v, reductions_type::unconfined);

		if (weights[v] >= all_neighbors_weight[v])
		{
			reduction_offset += include_vertex(v);
			continue;
		}

		if (degree[v] < 2 || weights[v] == 0)
		{
			reduction_offset += degree01_process(v);
			continue;
		}

		if (degree[v] >= 3 && degree[v] <= MAX_SMALL_GRAPH_SIZE)
		{
			vector<vertex> vertices;
			SET &is_mapped = set_buffer4;
			is_mapped.clear();

			size_t idx_cnt = 0;
			For(i, v)
			{
				vertices.push_back(LBLB(i));
			}
			sort(vertices.begin(), vertices.end(), [&](const vertex &lhs, const vertex &rhs) { return weights[lhs] < weights[rhs]; });
			for (; idx_cnt < vertices.size(); idx_cnt++)
			{
				REMAP[vertices[idx_cnt]] = idx_cnt;
				is_mapped.add(vertices[idx_cnt]);
			}

			small_graph.build(idx_cnt, weights[v]);
			for (vertex u : vertices)
			{
				small_graph.assign_weight(REMAP[u], weights[u]);
				For(i, u)
				{
					if (is_mapped.get(LBLB(i)))
					{
						small_graph.add_edge(REMAP[u], REMAP[LBLB(i)]);
					}
				}
			}

			if (small_graph.is_independent_neighborhood() && weights[v] >= all_neighbors_weight[v] - weights[vertices[0]])
			{
				reduction_offset += fold_vertices(v);
				continue;
			}
			else if (small_graph.is_clique_neighborhood())
			{
				vertex max_isolated = v, max_neighbor = v;
				for (vertex u : vertices)
				{
					if (degree[u] == degree[v])
					{
						if (weights[u] > weights[max_isolated])
						{
							max_isolated = u;
						}
					}
					else if (weights[u] > weights[max_neighbor])
					{
						max_neighbor = u;
					}
				}

				v = max_isolated;

				if (weights[v] >= weights[max_neighbor])
				{
					reduction_offset += include_vertex(v);
					continue;
				}

				exclude_vertex(v);

				For(i, v)
				{
					auto u = LBLB(i);

					if (weights[v] >= weights[u])
					{
						exclude_vertex(u);
					}
					else
					{
						modify_weight(u, weights[u] - weights[v]);
#ifdef DEBUG
						local_solution.push(vertex_status(v, VS::VS_NOTIF, u));
#endif
					}
				}
#ifdef DEBUG
				local_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
				reduction_offset += weights[v];
				continue;
			}
			else if (weights[v] >= small_graph.max_is())
			{
				reduction_offset += include_vertex(v);
				continue;
			}
		}
	}

	return reduction_offset;
}

void GRAPH::get_dis2_neighbors(vertex v, vector<vertex> &neighbors)
{
	SET &is_visited = set_buffer0;
	is_visited.clear();

	is_visited.add(v);

	For(u, v)
	{
		is_visited.add(LBLB(u));
	}

	For(u, v)
	{
		For(w, LBLB(u))
		{
			if (!is_visited.get(LBLB(w)))
			{
				if (degree[LBLB(w)] <= MAX_SMALL_GRAPH_SIZE_H)
				{
					neighbors.push_back(LBLB(w));
				}
				is_visited.add(LBLB(w));
			}
		}
	}
}

weight_node GRAPH::double_heavy_set_reduction()
{
	if ((branching_stack.size() != 0) && (local_n < 100 || 9 * pre_n < 10 * local_n)) // && (branching_stack.size() % 3 != 0)
	{
		return 0;
	}

	weight_node reduction_offset = 0;

	vector<vertex> available_vertices;
	get_available_vertices(available_vertices);

	for (auto v : available_vertices)
	{
		if (!is_available[v])
			continue;

		if (weights[v] >= all_neighbors_weight[v])
		{
			reduction_offset += include_vertex(v);
			continue;
		}

		if (degree[v] < 2 || weights[v] == 0)
		{
			reduction_offset += degree01_process(v);
			continue;
		}

		if (degree[v] >= 2 && degree[v] <= MAX_SMALL_GRAPH_SIZE_H)
		{
			bool is_succes = false;
			vector<vertex> neighbors;
			SET &neighbors_set = set_buffer1;
			neighbors_set.clear();

			small_graph.heavy_v[0].reset();
			size_t idx_cnt = 0;
			small_graph.weight_v[0] = weights[v];
			For(i, v)
			{
				neighbors.push_back(LBLB(i));
				neighbors_set.add(LBLB(i));
				small_graph.heavy_v[0][idx_cnt++] = true;
			}

			vector<vertex> d2_neighbors;
			get_dis2_neighbors(v, d2_neighbors);

			for (auto w : d2_neighbors)
			{

				if (w <= v)
				{
					continue;
				}
				auto vertices = neighbors;
				For(k, w)
				{
					if (!neighbors_set.get(LBLB(k)))
					{
						vertices.push_back(LBLB(k));
						if (vertices.size() > MAX_SMALL_GRAPH_SIZE_H)
						{
							break;
						}
					}
				}
				if (vertices.size() > MAX_SMALL_GRAPH_SIZE_H)
				{
					continue;
				}

				SET &is_mapped = set_buffer4;
				is_mapped.clear();

				idx_cnt = 0;
				for (; idx_cnt < vertices.size(); idx_cnt++)
				{
					REMAP[vertices[idx_cnt]] = idx_cnt;
					is_mapped.add(vertices[idx_cnt]);
				}
				small_graph.build(idx_cnt, 0);
				small_graph.heavy_v[1].reset();
				For(l, w)
				{
					small_graph.heavy_v[1][REMAP[LBLB(l)]] = true;
				}
				small_graph.weight_v[1] = weights[w];
				for (vertex ver : vertices)
				{
					small_graph.assign_weight(REMAP[ver], weights[ver]);
					For(i, ver)
					{
						if (is_mapped.get(LBLB(i)))
						{
							small_graph.add_edge(REMAP[ver], REMAP[LBLB(i)]);
						}
					}
				}

				if (small_graph.is_heavy_set())
				{
					reduction_offset += include_vertex(v);
					reduction_offset += include_vertex(w);
					is_succes = true;
				}

				if (is_succes)
				{
					break;
				}
			}
		}
	}

	return reduction_offset;
}

weight_node GRAPH::MAX_IS()
{
	return is_weight;
}

void GRAPH::split_components(vector<vector<vertex>> &components)
{
	SET &is_visited = set_buffer0;
	is_visited.clear();
	size_t component_cnt = 0;
	for (size_t i = LBR(0); i != 0; i = LBR(i))
	{
		if (!is_visited.get(i))
		{
			components.push_back(vector<vertex>());
			vector<vertex> &now_component = components[component_cnt++];

			queue<vertex> Q;
			Q.push(i);
			is_visited.add(i);
			now_component.push_back(i);

			while (!Q.empty())
			{
				auto top = Q.front();
				Q.pop();

				For(j, top)
				{
					if (is_visited.get(LBLB(j)))
						continue;
					Q.push(LBLB(j));
					now_component.push_back(LBLB(j));
					is_visited.add(LBLB(j));
				}
			}
		}
	}
}

weight_node GRAPH::branch_in_component(GRAPH *componnent_alg, vector<vertex> &component, weight_node lower_weight, bool called_flag)
{

	weight_node max_is = 0;

	componnent_alg->build(component.size(), 10 * component.size());
	for (size_t i = 0; i < component.size(); i++)
	{
		REMAP[component[i]] = i + 1;
		componnent_alg->assign_weight(i + 1, weights[component[i]]);
	}
	size_t cnt_m = 0;
	for (size_t i = 0; i < component.size(); i++)
	{
		auto v = component[i];
		For(j, v)
		{
			componnent_alg->insert_edge(REMAP[v], REMAP[LBLB(j)]);
			cnt_m++;
		}
	}
	componnent_alg->is_called_from_split_direct = called_flag;

	componnent_alg->run_branch(lower_weight);
	max_is = componnent_alg->MAX_IS();

#ifdef DEBUG
	auto temp_ans = componnent_alg->export_best_is();
	for (size_t i = 0; i < component.size(); i++)
	{
		IS_STATUS[component[i]] = temp_ans[i + 1];
	}
#endif

	return max_is;
}

#ifdef DEBUG
dynamic_bitset<> GRAPH::export_best_is()
{
	local_solution = best_solution;
	IS_STATUS.reset();
	while (!local_solution.empty())
	{
		auto top = local_solution.top();
		local_solution.pop();

		if (top.vs == VS::VS_TOKEN)
		{
			continue;
		}

		if ((top.vs == VS::VS_INCLUDED) || (top.vs == VS::VS_IF && IS_STATUS[top.u]) || (top.vs == VS::VS_IFNOT && !IS_STATUS[top.u]))
		{
			IS_STATUS[top.v] = true;
		}
		else if ((top.vs == VS::VS_EXCLUDED) || (top.vs == VS::VS_NOTIF && IS_STATUS[top.u]) || (top.vs == VS::VS_NOTIFNOT && !IS_STATUS[top.u]))
		{
			IS_STATUS[top.v] = false;
		}
	}
	best_solution = local_solution;
	return IS_STATUS;
}
#endif

inline void GRAPH::update_solution(weight_node local_weight)
{
	if (local_weight > is_weight)
	{
		is_weight = local_weight;
#ifdef DEBUG
		best_solution = local_solution;
#endif
	}
}

inline void GRAPH::cout_best_is_now()
{
	std::chrono::duration<float> cost_time = std::chrono::system_clock::now() - start_time;
	cout << cost_time.count() << " " << is_weight << endl;
}

void GRAPH::local_search(weight_node local_weight) // Here we take 3 greedy strategies to find a lower bound fast.
{
	SET &is_included2ls = set_buffer0;
	is_included2ls.clear();
	vector<vertex> vertices;
	get_available_vertices(vertices);
#ifdef DEBUG
	auto temp_solution = local_solution;
#endif

	sort(vertices.begin(), vertices.end(), [&](const vertex lhs, const vertex rhs) {
		return (double)weights[lhs] / degree[lhs] > (double)weights[rhs] / degree[rhs];
	});
	weight_node ls_weight = 0;
	for (vertex v : vertices)
	{
		bool is_uncovered = false;

		For(i, v)
		{
			if (is_included2ls.get(LBLB(i)))
			{
				is_uncovered = true;
				break;
			}
		}
		if (is_uncovered)
		{
			continue;
		}

		ls_weight += weights[v];
		is_included2ls.add(v);
#ifdef DEBUG
		temp_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
	}

	ls_lower_weight = local_weight + ls_weight;
	if (ls_lower_weight > is_weight)
	{
		is_weight = ls_lower_weight;
#ifdef DEBUG
		best_solution = temp_solution;
#endif
	}

#ifdef DEBUG
	temp_solution = local_solution;
#endif
	is_included2ls.clear();

	sort(vertices.begin(), vertices.end(), [&](const vertex lhs, const vertex rhs) {
		return degree[lhs] < degree[rhs] || (degree[lhs] == degree[rhs] && weights[lhs] > weights[rhs]);
	});
	ls_weight = 0;
	for (vertex v : vertices)
	{
		bool is_uncovered = false;

		For(i, v)
		{
			if (is_included2ls.get(LBLB(i)))
			{
				is_uncovered = true;
				break;
			}
		}
		if (is_uncovered)
		{
			continue;
		}

		ls_weight += weights[v];
		is_included2ls.add(v);
#ifdef DEBUG
		temp_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
	}

	ls_lower_weight = local_weight + ls_weight;
	if (ls_lower_weight > is_weight)
	{
		is_weight = ls_lower_weight;
#ifdef DEBUG
		best_solution = temp_solution;
#endif
	}

#ifdef DEBUG
	temp_solution = local_solution;
#endif
	is_included2ls.clear();

	sort(vertices.begin(), vertices.end(), [&](const vertex lhs, const vertex rhs) {
		return weights[lhs] > weights[rhs];
	});
	ls_weight = 0;
	for (vertex v : vertices)
	{
		bool is_uncovered = false;

		For(i, v)
		{
			if (is_included2ls.get(LBLB(i)))
			{
				is_uncovered = true;
				break;
			}
		}
		if (is_uncovered)
		{
			continue;
		}

		ls_weight += weights[v];
		is_included2ls.add(v);
#ifdef DEBUG
		temp_solution.push(vertex_status(v, VS::VS_INCLUDED, 0));
#endif
	}

	ls_lower_weight = local_weight + ls_weight;
	if (ls_lower_weight > is_weight)
	{
		is_weight = ls_lower_weight;
#ifdef DEBUG
		best_solution = temp_solution;
#endif
	}
}

bool GRAPH::is_running_time_out()
{
	std::chrono::duration<float> cost_time = std::chrono::system_clock::now() - start_time;
	if (cost_time.count() >= OA.time_limit)
	{
		return true;
	}
	return false;
}

void GRAPH::reduce()
{
	weight_node local_weight = 0;
	pre_n = local_n;
	size_t re_cnt = 0;
	size_t old_local_weight;
	while (true)
	{
		old_local_weight = local_weight;

		for (size_t i = LBR(0); i != 0; i = LBR(i))
		{
			reduction_queue[0].push(i);
		}

		local_weight += apply_all_reductions();

		if (re_cnt == 0)
		{
			for (size_t i = 0; i < 7; i++)
			{
				first_reductions_cnt[i] = reductions_cnt[i];
				first_reductions_tims_cnt[i] = reductions_tims_cnt[i];
				first_reductions_weight_cnt[i] = reductions_weight_cnt[i];
			}
		}
		re_cnt++;
		if (old_local_weight == local_weight)
		{
			break;
		}
	}
	reduction_time = std::chrono::system_clock::now() - start_time;
	remaining_vertices = local_n;
	reduction_offset_all = local_weight;

	cout << reduction_time.count() << " " << local_n << " " << reduction_offset_all << endl;

	is_weight = local_weight;
}

void GRAPH::output_reduced_graph(string filepath)
{
	ofstream output(filepath.c_str());
	vector<vertex> vertices;
	get_available_vertices(vertices);
	size_t cnt_m = 0;
	long long total_weight = 0;
	for (size_t i = 0; i < vertices.size(); i++)
	{
		REMAP[vertices[i]] = i + 1;
		cnt_m += degree[vertices[i]];
		total_weight += weights[vertices[i]];
	}

	output << "p edge " << vertices.size() << " " << cnt_m / 2 << endl;
	for (size_t i = 0; i < vertices.size(); i++)
	{
		output << "v " << i + 1 << " " << weights[vertices[i]] << endl;
	}

	cnt_m = 0;
	for (size_t i = 0; i < vertices.size(); i++)
	{
		For(j, vertices[i])
		{
			if (REMAP[LBLB(j)] > REMAP[vertices[i]])
			{
				init_edges[cnt_m].first = REMAP[vertices[i]];
				init_edges[cnt_m].second = REMAP[LBLB(j)];
				cnt_m++;
			}
		}
	}
	sort(init_edges, init_edges + cnt_m);
	for (size_t i = 0; i < cnt_m; i++)
	{
		output << "e " << init_edges[i].first << " " << init_edges[i].second << endl;
	}
	output << "%" << total_weight << endl;
	cout << vertices.size() << " " << cnt_m << " " << total_weight << endl;
	output.close();
}

weight_node GRAPH::local_max_weight(vector<vertex> &vertices, weight_node lower_weight)
{
	SET &is_mapped = set_buffer4;
	is_mapped.clear();

	size_t idx_cnt = 0;
	sort(vertices.begin(), vertices.end(), [&](const vertex &lhs, const vertex &rhs) { return weights[lhs] < weights[rhs]; });
	for (; idx_cnt < vertices.size(); idx_cnt++)
	{
		REMAP[vertices[idx_cnt]] = idx_cnt;
		is_mapped.add(vertices[idx_cnt]);
	}
	weight_node min_weight = weights[vertices[0]];
	small_graph.build(idx_cnt, lower_weight);
	for (vertex u : vertices)
	{
		small_graph.assign_weight(REMAP[u], weights[u]);
		For(i, u)
		{
			if (is_mapped.get(LBLB(i)))
			{
				small_graph.add_edge(REMAP[u], REMAP[LBLB(i)]);
			}
		}
	}
	return small_graph.max_is();
}

void read_graph(size_t &n, size_t &m)
{
	ifstream in(OA.read_file);
	size_t flag;
	string line;
	in >> n >> m >> flag;

	init_edges = new pair<size_t, size_t>[2 * m];
	init_weights = new weight_node[2 * n];

	int cnt = 0;
	weight_node total_weight = 0;
	getline(in, line);
	for (size_t i = 1; i <= n; i++)
	{
		getline(in, line);
		stringstream ss(line);

		ss >> init_weights[i];
		total_weight += init_weights[i];
		while (ss >> init_edges[cnt].second)
		{
			init_edges[cnt].first = i;
			cnt++;
		}
	}
	in.close();

	cout << n << " " << m << " " << total_weight << endl;
}

void init_buffers(size_t n, size_t m)
{
	size_t graph_size = log(n);
	graph_size *= 5;
	GRAPH_BUFFER = new GRAPH[graph_size];
	for (size_t i = 0; i < graph_size; i++)
	{
		GRAPH_BUFFER[i].depth = i;
	}
	LARGE_BUFFER = new list_node[6 * (n + m)];
	WEIGHT_BUFFER = new weight_node[6 * n];
	NEIGHBORS_WEIGHT_BUFFER = new weight_node[6 * n];
	DEGREE_BUFFER = new size_t[6 * n];
	flow_graph = new ISAP(4 * n, 4 * (n + m));
	set_buffer0 = SET(2 * n);
	set_buffer1 = SET(2 * n);
	set_buffer2 = SET(2 * n);
	set_buffer3 = SET(2 * n);
	set_buffer4 = SET(2 * n);
	cap_weight = new weight_node[2 * n];
	REMAP = new size_t[2 * n];
}
void delete_buffers()
{
	delete[] LARGE_BUFFER;
	delete[] WEIGHT_BUFFER;
	delete[] NEIGHBORS_WEIGHT_BUFFER;
	delete[] DEGREE_BUFFER;
	delete[] GRAPH_BUFFER;
	delete flow_graph;
	delete[] cap_weight;
	delete[] REMAP;
}

int main(int argc, char *argv[])
{
	OA.get_arg(argc, argv);

	if (OA.help)
		return 0;

	size_t n, m;
	read_graph(n, m);
	init_buffers(n, m);

	auto ALG = GRAPH_BUFFER;
	ALG->reset(LARGE_BUFFER, WEIGHT_BUFFER, NEIGHBORS_WEIGHT_BUFFER, DEGREE_BUFFER);
	ALG->build(n, 2 * m, init_edges, init_weights);
	ALG->is_fisrt_runed = true;

	start_time = std::chrono::system_clock::now();
	if (OA.op)
		ALG->reduce();
	else
		ALG->run_branch(0);

	auto end_exact = std::chrono::system_clock::now();
	std::chrono::duration<float> branch_reduce_time = end_exact - start_time;
	cout << branch_reduce_time.count() << " " << GRAPH_BUFFER->MAX_IS() << endl;

	if (OA.op)
	{

		for (int i = 0; i < 7; i++)
		{
			cout << first_reductions_cnt[i] << " " << first_reductions_tims_cnt[i] << " " << first_reductions_weight_cnt[i] << endl;
		}

		for (int i = 0; i < 7; i++)
		{
			cout << reductions_cnt[i] << " " << reductions_tims_cnt[i] << " " << reductions_weight_cnt[i] << endl;
		}

		if (OA.output_flag)
			ALG->output_reduced_graph(OA.output_file);
	}
	else
	{
		for (int i = 0; i < 7; i++)
		{
			cout << reductions_tims_cnt[i] << " " << reductions_cnt[i] << " " << reductions_weight_cnt[i] << endl;
		}

#ifdef DEBUG
		auto ans = ALG->export_best_is(); // the maximum weighted independent set

		bool is_independent = true;
		for (auto i = 0; i < m; i++)
		{
			if (ans[init_edges[i].first] && ans[init_edges[i].second])
			{
				is_independent = false;
				cout << init_edges[i].first << " " << init_edges[i].second << endl;
				break;
			}
		}

		if (!is_independent)
		{
			cout << "WRONG! NOT INDEPENDENT!" << endl;
		}
		weight_node is_weight = 0;
		for (size_t i = 1; i <= n; i++)
		{
			if (ans[i])
				is_weight += init_weights[i];
		}

		cout << is_weight << endl;

		if (OA.output_flag)
		{
		}
#endif
	}

	delete_buffers();
	return 0;
}