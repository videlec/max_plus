#include "paw.h"

#include <iostream>
#include <stdexcept>

PAW::PAW()
{
    this->n0 = new int_vec();
    this->n1 = new int_vec();
    this->w = 1;
}

PAW::PAW(uint64_t w)
{
    if (w <= 1)
        throw std::domain_error("too short word");
    this->w = w;
    size_t n = word_length(w);
    this->n0 = new int_vec(n, (uint64_t) 0);
    this->n1 = new int_vec(n, (uint64_t) 0);
}

PAW::PAW(uint64_t w, int_vec& n0, int_vec& n1)
{
    if (w <= 1)
        throw std::domain_error("too short word");
    this->n0 = new int_vec();
    this->n1 = new int_vec();
    this->set(w, n0, n1);
}


void PAW::set(uint64_t w, int_vec& n0, int_vec& n1)
{
    if (w <= 1)
        throw std::domain_error("too short word");
    size_t n = word_length(w);
    if (n0.size() != n-1 or n1.size() != n-1)
        throw std::length_error("wrong lengths");

    this->w = w;
    *(this->n0) = n0;
    *(this->n1) = n1;
}



PAW::~PAW()
{
    delete this->n0;
    delete this->n1;
}

void PAW::write()
{
    std::cout << this;
}

size_t PAW::size()
{
    return this->n0->size() + 1;
}

size_t PAW::length()
{
    return int_vec_sum(*(this->n0)) + int_vec_sum(*(this->n1)) + this->n0->size() + 1;
}


bool operator==(const PAW& u, const PAW& v)
{
    return (u.w == v.w &&
            int_vec_eq(*(u.n0), *(v.n0)) &&
            int_vec_eq(*(u.n1), *(v.n1)));
}

bool operator==(const uint64_t& u, const PAW& v)
{
    return (u == v.w &&
            int_vec_is_zero(*(v.n0)) &&
            int_vec_is_zero(*(v.n1)));
}

bool operator==(const PAW& u, const uint64_t& v)
{
    return (u.w == v &&
            int_vec_is_zero(*(u.n0)) &&
            int_vec_is_zero(*(u.n1)));
}

bool operator!=(const PAW& u, const PAW& v)
{
    return (u.w != v.w ||
            !int_vec_eq(*(u.n0), *(v.n0)) ||
            !int_vec_eq(*(u.n1), *(v.n1)));
}

bool operator!=(const uint64_t& u, const PAW& v)
{
    return (u != v.w ||
            !int_vec_is_zero(*(v.n0)) ||
            !int_vec_is_zero(*(v.n1)));
}

bool operator!=(const PAW& u, const uint64_t& v)
{
    return (u.w != v || 
            !int_vec_is_zero(*(u.n0)) ||
            !int_vec_is_zero(*(u.n1)));
}

std::ostream& operator<<(std::ostream& os, const PAW& dt)  
{
    uint64_t u = dt.w;
    os << (u&1);
    u >>= 1;
    int_vec::iterator it0 = dt.n0->begin();
    int_vec::iterator it1 = dt.n1->begin();
    while (u != 1)
    {
        if (*it0 || *it1)
            os << "(" << *(it0) << "," << *(it1) << ")";
        os << (u&1);
        it0++;
        it1++;
        u >>= 1;
    }
    return os;
}


/* we need to replace iterators
 *   product from itertools
 *   IntegerListsLex from Sage
 */
