#ifndef _OTHERS_H_
#define _OTHERS_H_

#include <bits/stdc++.h>
#include <boost/format.hpp>
#include <boost/dynamic_bitset.hpp>

using boost::dynamic_bitset;
using boost::format;
using std::bitset;
using std::cout;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::ofstream;
using std::pair;
using std::queue;
using std::stack;
using std::string;
using std::stringstream;
using std::swap;
using std::vector;

#define LBL(x) list_buffer[x].L
#define LBR(x) list_buffer[x].R
#define LBU(x) list_buffer[x].U
#define LBD(x) list_buffer[x].D
#define LBLB(x) list_buffer[x].label
#define size_t unsigned int
#define vertex unsigned int
#define weight_node int
#define INF 0x3f3f3f3f
#define SPLIT_LIMIT 100
#define LOCAL_SEARCH_LIMIT 10
#define MAX_SMALL_GRAPH_SIZE 8
#define MAX_SMALL_GRAPH_SIZE_H 8
#define rep(i, a, b) for (int i = a; i <= b; i++)
#define For(v, u) for (vertex v = list_buffer[u].D; v != u; v = list_buffer[v].D)
#define DEBUG

class SMALL_GRAPH
{
private:
	bitset<MAX_SMALL_GRAPH_SIZE> neighbors[MAX_SMALL_GRAPH_SIZE];
	weight_node weights[MAX_SMALL_GRAPH_SIZE];
	size_t n;
	weight_node lower_weight;

public:
	bitset<MAX_SMALL_GRAPH_SIZE> heavy_v[2];
	weight_node weight_v[2];

	void build(size_t &n_, weight_node lower_weight_);
	void assign_weight(vertex &v, weight_node &weight_v);
	void add_edge(vertex &u, vertex &v);
	weight_node max_is();
	bool is_clique_neighborhood();
	bool is_independent_neighborhood();
	bool is_heavy_set();
};

class SET
{
private:
	size_t *buffer;
	size_t valid_flag = 0;
	size_t valid_size = 0;
	size_t buffer_size = 0;

public:
	SET(size_t sz);
	SET();
	void clear();
	bool get(size_t ele);
	void add(size_t ele);
	void remove(size_t ele);
	size_t size();
};

class ISAP
{
private:
	const int inf = 1 << 28;
	int n, m, cnt, s, t;
	int *head, *pnt, *nxt, *flow, *pre, *d, *iter, *gap;
	bool *vis;

	void spfa();
	int augement();

public:
	ISAP(size_t max_n, size_t max_m);
	~ISAP();
	void init(int n_);
	void addedge(int u, int v, int f);
	int Max_flow();
	int ask_flow(size_t idx);
};
#endif