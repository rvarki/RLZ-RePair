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
        std::list<int> content;
        bool exp;
        RefNode* lnode;
        RefNode* rnode;
        int ltmp; // Needed for case where consecutive non exp phrase are replacing their edges in the same step
        int rtmp; // Needed for case where consecutive non exp phrase are replacing their edges in the same step
        PhraseNode* prev;
        PhraseNode* next;

        PhraseNode(RefNode* lnode, RefNode* rnode) : content({}), exp(false), lnode(lnode), rnode(rnode), ltmp(-1), rtmp(-1), prev(nullptr), next(nullptr) {}
        PhraseNode(std::list<int> ilist) : content(ilist), exp(true), lnode(nullptr), rnode(nullptr), ltmp(-1), rtmp(-1), prev(nullptr), next(nullptr) {}
};


#endif // NODE_H