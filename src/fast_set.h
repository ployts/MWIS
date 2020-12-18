#pragma once
#include "definitions.h"

class SET
{
private:
    size_t *buffer;
    size_t valid_flag = 0;
    size_t valid_size = 0;
    size_t buffer_size = 0;

public:
    SET(size_t sz)
    {
        buffer = new size_t[sz];
        buffer_size = sz;
    }
    SET(){
        
    }
    void clear()
    {
        valid_flag++;
        valid_size = 0;
        if (valid_flag == 1 << 31)
        {
            memset(buffer, 0, sizeof(size_t)*buffer_size);
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