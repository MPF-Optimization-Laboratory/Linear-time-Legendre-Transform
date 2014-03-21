#ifndef BIVECTOR_H
#define BIVECTOR_H

#include <cmath>
#include <assert.h>
#include <vector>

using std::vector;

template<class T>
struct tuple
{
    T x;
    T y;

    tuple(T x, T y) : x(x), y(y) {}
};

// a wrapper to represent a univariate set
// of function parameters and values with all
// other variables held constant
template<class T>
class BiVector
{
public:
    BiVector(vector<T> &domain, vector<T> &range) : m_domain(&domain), m_range(&range) {}

    tuple<T> operator[](int i)
    {
        return tuple<T>((*m_domain)[i], (*m_range)[i]);
    }
    tuple<T> back()
    {
        return tuple<T>(m_domain->back(), m_range->back());
    }
    void push_back(T x, T y)
    {
        m_domain->push_back(x);
        m_range->push_back(y);
    }
    void push_back(tuple<T> p)
    {
        m_domain->push_back(p.x);
        m_range->push_back(p.y);
    }
    void pop_back()
    {
        m_domain->pop_back();
        m_range->pop_back();
    }
    float slope(int i)
    {
        if (i <= 0)
            return -3.40282347E+37F;
        if (i >= size())
            return 3.40282347E+37F;

        return ((*m_range)[i] - (*m_range)[i-1]) / ((*m_domain)[i] - (*m_domain)[i-1]);
    }
    float slope(tuple<T> &p)
    {
        return (p.y - m_range->back()) / (p.x - m_domain->back());
    }
    int size()
    {
        assert(m_domain->size() == m_range->size());
        return m_domain->size();
    }

    vector<T> *m_domain;
    vector<T> *m_range;
};

// a wrapper to represent a univariate set
// of function parameters and values with all
// other variables held constant
template<class T>
class BiArray
{
public:

    BiArray(int size, T *domain, T *range) : m_size(size), m_domain(domain), m_range(range) {}

    tuple<T> operator[](int i)
    {
        return tuple<T>(m_domain[i], m_range[i]);
    }
    tuple<T> back()
    {
        return tuple<T>(m_domain[m_size-1], m_range[m_size-1]);
    }
    float slope(int i)
    {
        if (i <= 0)
            return -3.40282347E+37F;
        if (i >= size())
            return 3.40282347E+37F;

        return (m_range[i] - m_range[i-1]) / (m_domain[i] - m_domain[i-1]);
    }
    float slope(tuple<T> &p)
    {
        return (p.y - m_range[m_size-1]) / (p.x - m_domain[m_size-1]);
    }
    int size()
    {
        return m_size;
    }

    int m_size;
    T  *m_domain;
    T  *m_range;
};

#endif // BIVECTOR_H