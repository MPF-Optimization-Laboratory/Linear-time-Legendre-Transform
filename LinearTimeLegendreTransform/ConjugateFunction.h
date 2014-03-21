#ifndef LUCET_LLT_H
#define LUCET_LLT_H

#define _USE_MATH_DEFINES

#include <vector>
#include <string.h>
#include "BiVector.h"

using std::vector;

struct iter_info
{
    /* constant state */
    int  N;        // the number of dimensions
    int  x_size;   // the size of domain (constant)
    int  s_size;   // the size of the slope array
    int  alloc;    // the amount of memory to allocate on each iteration for a partial solution

    /* iteration dependant state */
    int  subN;      // the subset of R^N currently being operated on
    int  iter_size;
    int  num_elem;  // the number of elements to iterate over this round
    int  s_elem;    // the number of s elements to skip in the solution
    int *index;     // a domain tuple (x1, x2, ..., xn), needs to be freed at the end
};

/* uses the beneath beyond algorithm to compute the lower convex
 * hull of a linear set of discrete points (this will be just the
 * end points for a concave function)
 *
 * inputs:
 * info: information regarding offsets, sizes and indices for pointer arithmatic
 * x: the set of x values in the domain
 * y: the set of values f(x)
 *
 * in/out:
 * hull_x: the set of x values in the hull, subset of x
 * hull_y: the set of y values in the hull, subset of y
 * hull_s: the set of slopes between each of the points in the hull */
void computeHull(iter_info &info, float *x, float *y, vector<float> &hull_x, vector<float> &hull_y, vector<float> &hull_s)
{
    BiArray <float> in(info.x_size, x, y);
    BiVector<float> out(hull_x, hull_y);

    out.push_back(in[0]);
    out.push_back(in[1]);
    hull_s.push_back(out.slope(1));

    float slope = 0;
    for (int i = 2; i < in.size(); i++)
    {
        slope = out.slope(in[i]); 
        while (hull_s.size() > 0 && slope < hull_s.back())
        {
            out.pop_back();
            hull_s.pop_back();
            slope = out.slope(in[i]);
        }
        if (hull_s.size() > 0 && slope == hull_s.back())
            out.pop_back();
        else
            hull_s.push_back(slope);
        out.push_back(in[i]);
        
    }
}

/* This function computes a pointer to the
 * array of function values for the given domain,
 * with all but the last value held fixed 
 * using the indices at which the domain values are stored
 *
 * input:
 * info:   information regarding offsets, sizes and indices for pointer arithmatic
 *         including the indices at which we are currently evaluating
 * values: the array of pointers to function values from which to retreive the
 *         answer
 *
 * output: 
 * the vector with all of the inputs held constant except x_n */
template <class C>
C *evaluate(iter_info &info, C *values)
{
    int increment = 0;
    for (int i = 0, j = info.subN; j >= 0; i++, j--)
    {
        increment += info.index[i] * pow((double)info.iter_size, j);
    }
    return values + increment;
}

/* inputs:
 * info:   information regarding offsets, sizes and indices for pointer arithmatic
 *            including the indices at which we are currently evaluating
 * domain: an array of N arrays of x values, each array representing one component
 *            of the input function f's domain
 * slopes: an array of N sets of values, the slopes, s_i, in each direction we are
 *            evaluating the conjugate at
 * prev_soln: a matrix of dimension N, with k dimensions of size x_size, and N-k
 *            dimensions of size s_size, as output from the previous iteration of
 *            function
 *
 * in/out:
 * soln: the location in memory to write the final solution (cannot be reused for
 *       each iteration due to packing issues) */
void recurse(iter_info info, float *domain, float *slopes, float *prev_soln, float *soln)
{
    // update step specific variables
    info.subN--;
    info.num_elem /= info.x_size;
    info.s_elem   /= info.s_size;

    // because the dimensions of s and x are different, the easiest solution
    // to ensure we don't overwrite each other is to allocate new memory
    float *part_soln = new float[info.alloc];

    // for each possible combination of x values with constant last element varying
    for (int i = 0; i < info.num_elem; i++)
    {
        vector<float> hull_x, hull_y, hull_s; // the domain and slope of each segment of the hull

        // acquire the function values varying only x_n
        float *f = evaluate(info, prev_soln);

        // compute the hull of {x_n, u(x1, x2, ..., xn) for fixed [x_1,x_{n-1}]
        computeHull(info, &domain[info.x_size*info.subN], f, hull_x, hull_y, hull_s);

        // compute partial solutions v_{s_n}(x_1,...,x_{n-1}) : x_n \in X_n
        for (int j = 0, k = 0; j < info.s_size; j++)
        {
            float s = slopes[info.s_size * info.subN + j];
            while (k < hull_s.size() && hull_s[k] < s)
                k++;

            // account for univariate and multi variate cases
            float value = hull_x[k] * s - hull_y[k];

            // in order to only evaluate the lower hull
            // we negate the values (meaning the next iteration will add them)
            // resulting in hull_x[k] * s + hull_y[k];
            // for the next iteration
            if (info.subN != 0)
                value *= -1;

            // if it is the final iteration, then write to the solution array
            // else index into our partial array in a way to make the next iteration's read sequential
            if (info.subN == 0)
                soln[i*info.s_size + j] = value;
            else 
                part_soln[j*info.num_elem + i] = value;
        }

        // get the next combination of x values f(x1, ..., x_n)
        //for (int j = 0; j < info.N && ++info.index[j] == info.x_size; j++)
        for (int j = info.subN-1; j >= 0 && ++info.index[j] == info.x_size; j--)
            info.index[j] = 0;
    }

    if (info.subN > 0)
        for (int i = 0; i < info.s_size; i++)
        {
            // evaluations in all iterations past the first are in terms of s_size
            info.iter_size = info.s_size;
            memset(info.index, 0, sizeof(int)*info.N);
            recurse(info, domain, slopes, &part_soln[i*info.num_elem], &soln[i*info.s_elem]);
        }

    delete [] part_soln;
}

/* inputs:
 * N: the dimension of the domain, x \in R^n
 * X: an array of n vectors of the same size, specifying the values in each
 *    dimension to combine to make the domain. (ie X[0][#] = x0, X[1] = x1 for all f(x0, x1, ...) )
 * Fx: a matrix \in R^n+1 containing the function evaluated at all the possible combinations of the
 *     array X of size x_size^n
 * x_size: the number of elements in the domain,
 *             where X is of size: N*x_size
 *             and Fx is of size x_size^N
 * s: the slope values at which to evaluate the conjugate function (of dimension N)
 * s_size: the size of the array s
 *             where s is of size: N*s_size
 *
 * outputs:
 * solution: the conjugate values at s. Must be of size s_size^n. */
void linearLegendre(int N, float *X, float *Fx, int x_size, float *s, int s_size, float *solution)
{
    int size_in  = pow((float)x_size, N);
    int size_out = pow((float)s_size, N);

    // this is info that will be needed for each stage
    iter_info info;
    info.N        = N;
    info.s_size   = s_size;
    info.x_size   = x_size;
    info.alloc    = (x_size > s_size) ? size_in : size_out;

    info.subN      = N;
    info.iter_size = x_size;
    info.num_elem  = size_in;
    info.s_elem    = size_out;
    info.index     = new int[N];                // the current index we are computing

    // set our index array to zeros
    memset(info.index, 0, sizeof(int)*N);

    // the pointer is nessecary because we don't know which of
    recurse(info, X,  s, Fx, solution);

    // cleanup resources
    delete [] info.index;
}

#endif // LUCET_LLT_H