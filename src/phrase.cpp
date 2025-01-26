#include "phrase.h"

Phrase::Phrase()
    : exp(false), lrange(-1), rrange(-1), lneighbor(nullptr), rneighbor(nullptr) {}

Phrase::Phrase(int lrange, int rrange)
    : exp(false), lrange(lrange), rrange(rrange), lneighbor(nullptr), rneighbor(nullptr) {}

Phrase::Phrase(std::list<unsigned char> content)
    : exp(true), lrange(-1), rrange(-1), lneighbor(nullptr), rneighbor(nullptr), content(content) {}

Phrase::~Phrase() {}