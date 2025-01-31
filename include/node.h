#ifndef NODE_H
#define NODE_H
#include <string>

class RefNode
{
    public:
        int val;
        bool deleted;
        RefNode* prev;
        RefNode* next;
        RefNode(int value) : val(value), deleted(false), prev(nullptr), next(nullptr) {}
};

class PhraseNode
{
    public: 
        std::string content;
        RefNode* lnode;
        RefNode* rnode;
        PhraseNode* prev;
        PhraseNode* next;

        PhraseNode(RefNode* lrange, RefNode* rrange) : content(""), lnode(lrange), rnode(rrange), prev(nullptr), next(nullptr) {}
        PhraseNode(std::string value) : content(value), lnode(nullptr), rnode(nullptr), prev(nullptr), next(nullptr) {}
};


#endif // NODE_H