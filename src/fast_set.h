#pragma once
#include "definitions.h"

class SET
{
private:
    size_t buffer[MAX_NUM_VERTICES];
    size_t valid_flag = 0;
    size_t valid_size = 0;

public:
    void clear()
    {
        valid_flag++;
        valid_size = 0;
        if (valid_flag == 1 << 31)
        {
            memset(buffer, 0, sizeof(buffer));
        }
    }
    bool get(size_t ele)
    {
        return buffer[ele] >= valid_flag;
    }
    void add(size_t ele)
    {
        if (!get(ele))
        {
            buffer[ele] = valid_flag;
            valid_size++;
        }
    }
    void remove(size_t ele)
    {
        if (get(ele))
        {
            valid_size--;
            buffer[ele] = 0;
        }
    }
    size_t size()
    {
        return valid_size;
    }
};