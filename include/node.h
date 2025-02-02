#ifndef NODE_H
#define NODE_H
#include <list>

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
        std::list<unsigned int> content;
        bool exp;
        RefNode* lnode;
        RefNode* rnode;
        PhraseNode* prev;
        PhraseNode* next;
        int lrange;
        int rrange;

        PhraseNode(RefNode* lnode, RefNode* rnode, int lrange, int rrange) : content({}), exp(false), lnode(lnode), rnode(rnode), lrange(lrange), rrange(rrange), prev(nullptr), next(nullptr) {}
        PhraseNode(std::list<unsigned int> ilist) : content(ilist), exp(true), lnode(nullptr), rnode(nullptr), lrange(-1), rrange(-1), prev(nullptr), next(nullptr) {}
};


#endif // NODE_H