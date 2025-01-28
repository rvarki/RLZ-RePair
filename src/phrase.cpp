#include "phrase.h"

Phrase::Phrase()
    : exp(false), lrange(-1), rrange(-1) {}

Phrase::Phrase(int lrange, int rrange)
    : exp(false), lrange(lrange), rrange(rrange) {}

Phrase::Phrase(std::list<unsigned char> content)
    : exp(true), lrange(-1), rrange(-1), content(content) {}

Phrase::~Phrase() {}