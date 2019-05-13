#ifndef _FIBO_MIN_PROD_UTILS_H
#define _FIBO_MIN_PROD_UTILS_H

#include <iostream>
#include <list>
#include <vector>

/* type for vector of occurrences */
typedef std::vector<uint64_t> int_vec;

int cmp_lexico(int_vec&, int_vec&);
int cmp_prod(int_vec &u, int_vec &v);

uint64_t int_vec_sum(int_vec& u);
size_t word_length(uint64_t w);
void word_count(size_t& n0, size_t& n1, uint64_t w);

/* check whether the PAF (w, n0, n1) matches the word between begin and end */ 
int match(int_vec::iterator begin, int_vec::iterator end,
          uint64_t w, int_vec& n0, int_vec& n1);


inline bool int_vec_eq(int_vec& u, int_vec& v)
{
    if (u.size() != v.size()) return false;
    int_vec::iterator itu = u.begin();
    int_vec::iterator itv = v.begin();
    while (itu != u.end())
    {
        if (*itu != *itv) return false;
        itu++;
        itv++;
    }
    return true;
}

inline bool int_vec_is_zero(int_vec& u)
{
    for (int_vec::iterator i = u.begin(); i != u.end(); i++)
    {
        if (*i != 0) return false;
    }
    return true;

}

/*
class substitutive_word
{
    int_vec w;
    int_vec s0; // image of 0
    int_vec s1; // image of 1
    size_t i;   // position of the pointer for the next substitution

    substitutive_word(int_vec& w0, in_vec& w1);
    uint64_t operator[](size_t i);
};
*/

#endif

