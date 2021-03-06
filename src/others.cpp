#include "others.h"

void SMALL_GRAPH::build(size_t &n_, weight_node lower_weight_)
{
    n = n_;
    lower_weight = lower_weight_;
    for (size_t i = 0; i < MAX_SMALL_GRAPH_SIZE; i++)
    {
        neighbors[i].reset();
    }
}
void SMALL_GRAPH::assign_weight(vertex &v, weight_node &weight_v)
{
    weights[v] = weight_v;
}
void SMALL_GRAPH::add_edge(vertex &u, vertex &v)
{
    neighbors[u][v] = true;
}
weight_node SMALL_GRAPH::max_is()
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

bool SMALL_GRAPH::is_clique_neighborhood()
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

bool SMALL_GRAPH::is_independent_neighborhood()
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

bool SMALL_GRAPH::is_heavy_set()
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

SET::SET(size_t sz)
{
    buffer = new size_t[sz];
    buffer_size = sz;
}
SET::SET()
{
}
void SET::clear()
{
    valid_flag++;
    valid_size = 0;
    if (valid_flag == 1 << 31)
    {
        memset(buffer, 0, sizeof(size_t) * buffer_size);
    }
}
const bool SET::operator[](const size_t &ele) const
{
    return buffer[ele] >= valid_flag;
}
void SET::operator+(const size_t &ele)
{
    if (buffer[ele] < valid_flag)
    {
        buffer[ele] = valid_flag;
        valid_size++;
    }
}
void SET::operator-(const size_t &ele)
{
    if (buffer[ele] >= valid_flag)
    {
        valid_size--;
        buffer[ele] = 0;
    }
}
size_t SET::size()
{
    return valid_size;
}

void ISAP::spfa()
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
int ISAP::augement()
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
ISAP::ISAP(size_t max_n, size_t max_m)
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
ISAP::~ISAP()
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
void ISAP::init(int n_)
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
void ISAP::addedge(int u, int v, int f)
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
int ISAP::Max_flow()
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
int ISAP::ask_flow(size_t idx)
{
    return flow[idx];
}

opt_arg::opt_arg()
{
    time_limit = 1000.0;
    output_flag = false;
    op = true;
    help = false;
    opt_string = "ho:t:r:w";
}

void opt_arg::get_arg(int argc, char *argv[])
{
    char opt;

    while ((opt = getopt(argc, argv, opt_string.c_str())) != -1)
    {
        switch (opt)
        {
        case 'h':
            help = true;
            cout << "-r the path of the input file" << endl;
            cout << "-o output_file, optional" << endl;
            cout << "-t the limit of running time, optional" << endl;
            cout << "-w run the whole algorithm" << endl;
            cout << "-h Helps" << endl;
            break;
        case 'r':
            read_file = string(optarg);
            break;

        case 'o':
            output_flag = true;
            output_file = string(optarg);
            break;
        case 'w':
            op = false;
            break;
        case 't':
            time_limit = atof(optarg);
            break;
        }
    }
}