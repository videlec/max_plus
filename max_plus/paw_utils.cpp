#include "paw_utils.h"

#include <stdexcept>

/* Words are encoded with an integer
 (the first 1 delimit the position)
 1 -> epsilon
 10 -> 0 
 11 -> 1 
 100 -> 00
 101 -> 01 
 110 -> 10
 111 -> 11
*/

/* compare for the product order */
/* (assuming they have the same length) */
/* -1: <  */
/* 0 : == */
/* 1 : >  */
/* 2 : incomparable */
int cmp_prod(int_vec &u, int_vec &v)
{
    int c = 0;
    int_vec::iterator itu = u.begin();
    int_vec::iterator itv = v.begin();
    while(itu != u.end() || itv != v.end())
    {
        if(*itu < *itv)
        {
            if(c == 0) c = -1;
            else if(c != -1) return 2;
        }
        else if (*itu > *itv)
        {
            if(c == 0) c = 1;
            else if(c != 1) return 2;
        }

        itu++;
        itv++;
    }

    return c;
}

int cmp_lexico(int_vec &u, int_vec& v)
{
    int_vec::iterator itu = u.begin();
    int_vec::iterator itv = v.begin();
    while (itu != u.end() || itv != v.end())
    {
        if (*itu < *itv) return -1;
        if (*itu > *itv) return 1;
        itu++;
        itv++;
    }
    if (itu != u.end()) return 1;
    if (itv != v.end()) return -1;

    return 0;
}

int cmp_paf(uint64_t u, int_vec& u_n0, int_vec &u_n1,
        uint64_t v, int_vec& v_n0, int_vec &v_n1)
{
    // compare lengths
    size_t lu = int_vec_sum(u_n0) + int_vec_sum(u_n1) + u_n0.size();
    size_t lv = int_vec_sum(v_n0) + int_vec_sum(v_n1) + v_n0.size();
    if (lu < lv) return -1;
    if (lu > lv) return 1;

    // compare sizes
    if (u_n0.size() < u_n1.size()) return -1;
    if (u_n0.size() > u_n1.size()) return 1;

    // don't care
    return -1;
}

uint64_t int_vec_sum(int_vec& u)
{
    uint64_t s = 0;
    for(int_vec::iterator itu = u.begin(); itu != u.end(); itu++)
        s += *itu;
    return s;
}

size_t word_length(uint64_t w)
{
    if (w == 0)
        throw std::domain_error("word can not be zero");
    size_t n = 0;
    while(w != 1)
    {
        w >>= 1;
        n += 1;
    }
    return n;
}

void word_count(size_t& n0, size_t& n1, uint64_t w)
{
    while(w != 1)
    {
        if (w & 1) n1 += 1;
        else n0 += 1;
        w >>= 1;
    }
}

int match(int_vec::iterator begin, int_vec::iterator end,
          uint64_t w, int_vec& n0, int_vec& n1)
{
    int_vec::iterator it0 = n0.begin();
    int_vec::iterator it1 = n1.begin();

    if (w <= 1)
        throw std::domain_error("empty word");
    if (begin == end)
       return 0;
    if ((w&1) != *begin)
        return 0;
    w >>= 1;
    begin++;
    
    while (it0 != n0.end() && it1 != n1.end() && w != 1 && begin != end)
    {
        // compare Abelian part
        size_t k0 = 0;
        size_t k1 = 0;
        for (size_t j = 0; j < *it0 + *it1; j++)
        {
            if (*begin == 0) k0++;
            else k1++;
            begin++;
            if (begin == end) return 0;
        }
        if (k0 != *it0 || k1 != *it1) return 0;
        it0++;
        it1++;

        // compare letter
        if ((w&1) != *begin) return 0;
        w >>= 1;
        begin++;
    }

    return it0 == n0.end() && it1 == n1.end() && w == 1 && begin == end;
}



