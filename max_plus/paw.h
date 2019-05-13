/********************************************************************/
/* Partially Abelianized Words (PAW)                                */
/*                                                                  */
/* Only support binary alphabet and PAWS of size <= 64              */
/********************************************************************/

#ifndef _FIBO_MIN_PROD_PAW_H
#define _FIBO_MIN_PROD_PAW_H

#include "paw_utils.h"

#include <iostream>

class PAW
{
    public:
        uint64_t w;    // word
        int_vec * n0;  // abelianized zeros
        int_vec * n1;  // abelianized ones

        PAW();
        PAW(uint64_t w);
        PAW(uint64_t w, int_vec& n0, int_vec& n1);

        ~PAW();

        size_t size();
        size_t length();

        void set(uint64_t w, int_vec& n0, int_vec& n1);

        void write();

        friend std::ostream& operator<<(std::ostream& os, const PAW& u);
};


bool operator==(const PAW&, const PAW&);
bool operator==(const uint64_t&, const PAW&);
bool operator==(const PAW&, const uint64_t&);

bool operator!=(const PAW&, const PAW&);
bool operator!=(const uint64_t&, const PAW&);
bool operator!=(const PAW&, const uint64_t&);

#endif
