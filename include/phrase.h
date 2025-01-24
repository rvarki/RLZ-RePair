#ifndef PHRASE_H
#define PHRASE_H

#include <string>

class Phrase {

    private:
        bool non_explicit;
        int lrange;
        int rrange;
        Phrase* lneighbor;
        Phrase* rneighbor;
        std::string content;
    
    public:
        // Constructors
        Phrase();
        explicit Phrase(int lrange, int rrange);
        explicit Phrase(std::string content);

        // Destructor
        ~Phrase();

        // Getters and Setters
        int get_lrange() const;
        void set_lrange(int lrange);
        int get_rrange() const;
        void set_rrange(int rrange);
        const std::string& get_content() const;
        void set_content(std::string content);
        Phrase* get_lneighbor() const;
        void set_lneighbor(Phrase* lneighbor);
        Phrase* get_rneighbor() const;
        void set_rneighbor(Phrase* rneighbor);
};

#endif // PHRASE_H
