#include "phrase.h"

Phrase::Phrase()
    : non_explicit(true), lrange(-1), rrange(-1), lneighbor(nullptr), rneighbor(nullptr) {}

Phrase::Phrase(int lrange, int rrange)
    : non_explicit(true), lrange(lrange), rrange(rrange), lneighbor(nullptr), rneighbor(nullptr) {}

Phrase::Phrase(std::string content)
    : non_explicit(false), lrange(-1), rrange(-1), lneighbor(nullptr), rneighbor(nullptr), content(content) {}

Phrase::~Phrase() {}

int Phrase::get_lrange() const { return lrange; }
void Phrase::set_lrange(int lrange) { this->lrange = lrange; }

int Phrase::get_rrange() const { return rrange; }
void Phrase::set_rrange(int rrange) { this->rrange = rrange; }

const std::string& Phrase::get_content() const { return content; }
void Phrase::set_content(std::string content) { this->content = content; }

Phrase* Phrase::get_lneighbor() const { return lneighbor; }
void Phrase::set_lneighbor(Phrase* lneighbor) { this->lneighbor = lneighbor; }

Phrase* Phrase::get_rneighbor() const { return rneighbor; }
void Phrase::set_rneighbor(Phrase* rneighbor) { this->rneighbor = rneighbor; }
