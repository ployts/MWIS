#pragma once
#include<bits/stdc++.h>
#include<boost/format.hpp>
#include<boost/dynamic_bitset.hpp>

using std::cout;
using std::vector;
using std::endl;
using std::bitset;
using std::queue;
using std::min;
using std::max;
using std::pair;
using std::stack;
using std::string;
using std::swap;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using boost::dynamic_bitset;
using boost::format;


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
