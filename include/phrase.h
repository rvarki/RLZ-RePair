#ifndef PHRASE_H
#define PHRASE_H

#include <string>
#include <list>

class Phrase {

    public:
        bool exp;
        int lrange;
        int rrange;
        std::list<unsigned int> content;

        // Constructors
        Phrase();
        explicit Phrase(int lrange, int rrange);
        explicit Phrase(std::list<unsigned int> content);

        // Destructor
        ~Phrase();

};

#endif // PHRASE_H
