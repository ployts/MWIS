#pragma once
#include "definitions.h"


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

	void build(size_t &n_, weight_node lower_weight_)
	{
		n = n_;
		lower_weight = lower_weight_;
		for (size_t i = 0; i < MAX_SMALL_GRAPH_SIZE; i++)
		{
			neighbors[i].reset();
		}
	}
	inline void assign_weight(vertex &v, weight_node &weight_v)
	{
		weights[v] = weight_v;
	}
	inline void add_edge(vertex &u, vertex &v)
	{
		neighbors[u][v] = true;
	}
	weight_node max_is()
	{
		weight_node max_weight = 0;
		weight_node is_weight;
		bool is_independent;
		for (size_t is_num = (1 << n) - 1; is_num > 0; is_num--)
		{
			bitset<MAX_SMALL_GRAPH_SIZE> is(is_num);
			is_independent = true;
			is_weight = 0;
			for (size_t i = 0; i < n; i++)
			{
				if (is[i])
				{
					is_weight += weights[i];
					if ((neighbors[i] & is).any())
					{
						is_independent = false;
					}
				}
				if (!is_independent)
				{
					break;
				}
			}
			if (is_independent && is_weight > max_weight)
			{
				max_weight = is_weight;
				if (max_weight > lower_weight)
				{
					return max_weight;
				}
			}
		}
		return max_weight;
	}

	bool is_clique_neighborhood()
	{
		for (size_t i = 0; i < n; i++)
		{
			if (neighbors[i].count() != n - 1)
			{
				return false;
			}
		}
		return true;
	}

	bool is_independent_neighborhood()
	{
		for (size_t i = 0; i < n; i++)
		{
			if (neighbors[i].count() != 0)
			{
				return false;
			}
		}
		return true;
	}

	bool is_heavy_set()
	{
		weight_node is_weight;
		bool is_independent;
		weight_node heavy_set_weight;
		for (size_t is_num = (1 << n) - 1; is_num > 0; is_num--)
		{
			bitset<MAX_SMALL_GRAPH_SIZE> is(is_num);
			is_independent = true;
			is_weight = 0;
			for (size_t i = 0; i < n; i++)
			{
				if (is[i])
				{
					is_weight += weights[i];
					if ((neighbors[i] & is).any())
					{
						is_independent = false;
					}
				}
				if (!is_independent)
				{
					break;
				}
			}
			if (is_independent)
			{
				heavy_set_weight = 0;
				for (size_t i = 0; i < 2; i++)
				{
					if ((heavy_v[i] & is).any())
					{
						heavy_set_weight += weight_v[i];
					}
				}
				if (heavy_set_weight < is_weight)
				{
					return false;
				}
			}
		}
		return true;
	}
};
