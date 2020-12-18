#ifndef DEFINITIONS_H
#define DEFINITIONS_H
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
#endif