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
        RefNode(int value) : val(value), pos(0), deleted(false), prev(nullptr), next(nullptr) {}
};

class PhraseNode
{
    public: 
        std::list<unsigned int> content;
        bool exp;
        bool leftReplaced;
        bool rightReplaced;
        RefNode* lnode;
        RefNode* rnode;
        PhraseNode* prev;
        PhraseNode* next;

        PhraseNode(RefNode* lnode, RefNode* rnode) : content({}), exp(false), leftReplaced(false), rightReplaced(false), lnode(lnode), rnode(rnode), prev(nullptr), next(nullptr) {}
        PhraseNode(std::list<unsigned int> ilist) : content(ilist), exp(true), leftReplaced(false), rightReplaced(false), lnode(nullptr), rnode(nullptr), prev(nullptr), next(nullptr) {}
};


#endif // NODE_H