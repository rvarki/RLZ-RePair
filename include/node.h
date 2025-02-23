#ifndef NODE_H
#define NODE_H
#include <list>

class RefNode
{
    public:
        int val;
        int pos;
        bool deleted;
        RefNode* prev;
        RefNode* next;
        RefNode* forward;
        RefNode(int value) : val(value), pos(0), deleted(false), prev(nullptr), next(nullptr), forward(nullptr) {}
};

class PhraseNode
{
    public: 
        std::list<int> content;
        bool exp;
        RefNode* lnode;
        RefNode* rnode;
        PhraseNode* prev;
        PhraseNode* next;

        PhraseNode(RefNode* lnode, RefNode* rnode) : content({}), exp(false), lnode(lnode), rnode(rnode), prev(nullptr), next(nullptr) {}
        PhraseNode(std::list<int> ilist) : content(ilist), exp(true), lnode(nullptr), rnode(nullptr), prev(nullptr), next(nullptr) {}
};


#endif // NODE_H