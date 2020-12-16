#pragma once
#include "definitions.h"

struct ISAP
{
	const int inf = 1 << 28;
	int n, m, cnt, s, t;
	int *head, *pnt, *nxt, *flow, *pre, *d, *iter, *gap;
	bool *vis;

	ISAP(size_t max_n, size_t max_m)
	{
		head = new int[max_n];
		pnt = new int[max_m];
		nxt = new int[max_m];
		flow = new int[max_m];
		pre = new int[max_n];
		d = new int[max_n];
		iter = new int[max_n];
		gap = new int[max_n];
		vis = new bool[max_n];
	}
	~ISAP()
	{
		delete[] head;
		delete[] pnt;
		delete[] nxt;
		delete[] flow;
		delete[] pre;
		delete[] d;
		delete[] iter;
		delete[] gap;
		delete[] vis;
	}

	void addedge(int u, int v, int f)
	{
		pnt[cnt] = v;
		nxt[cnt] = head[u];
		flow[cnt] = f;
		head[u] = cnt++;
		pnt[cnt] = u;
		nxt[cnt] = head[v];
		flow[cnt] = 0;
		head[v] = cnt++;
	}
	void init(int n_)
	{
		n = n_;
		cnt = 0;

		memset(head, -1, n << 2);
		memset(d, -1, n << 2);
		memset(vis, false, n);
		memset(gap, 0, n << 2);
		memset(pre, -1, n << 2);
		s = 0;
		t = n - 1;
	}
	void spfa()
	{
		queue<int> q;
		q.push(t);
		vis[t] = true;
		d[t] = 0;
		while (!q.empty())
		{
			int u = q.front();
			q.pop();
			for (int i = head[u]; i != -1; i = nxt[i])
				if (!vis[pnt[i]] && flow[i ^ 1] > 0)
				{
					vis[pnt[i]] = true;
					d[pnt[i]] = d[u] + 1;
					q.push(pnt[i]);
				}
		}
	}
	int augement()
	{
		int f = 0x3f3f3f3f;
		for (int i = pre[t]; i != -1; i = pre[pnt[i]])
			f = min(f, flow[i ^ 1]);
		for (int i = pre[t]; i != -1; i = pre[pnt[i]])
		{
			flow[i ^ 1] -= f;
			flow[i] += f;
		}
		return f;
	}
	int Max_flow()
	{
		int fl = 0;
		spfa();
		rep(i, 0, t) iter[i] = head[i];
		rep(i, 0, t) if (~d[i]) gap[d[i]]++;
		int u = 0;
		while (d[0] < (t + 1))
		{
			if (u == t)
			{
				fl += augement();
				u = 0;
			}
			bool adv = false;
			for (int &i = iter[u]; i != -1; i = nxt[i])
				if (flow[i] && d[u] == d[pnt[i]] + 1)
				{
					adv = true;
					pre[pnt[i]] = i ^ 1;
					u = pnt[i];
					break;
				}
			if (!adv)
			{
				int m = t + 1;
				for (int i = head[u]; i != -1; i = nxt[i])
					if (flow[i] && ~d[pnt[i]])
						m = min(m, d[pnt[i]]);
				if (--gap[d[u]] == 0)
					break;
				d[u] = m + 1;
				gap[d[u]]++;
				iter[u] = head[u];
				if (u != 0)
					u = pnt[pre[u]];
			}
		}
		return fl;
	}
};
