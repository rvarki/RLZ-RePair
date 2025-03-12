#ifndef NODE_H
#define NODE_H
#include <list>

class RefNode
{
    public:
        int val;
        int pos;
        int prev, next;
        bool deleted;
        RefNode() : val(-1), pos(-1), prev(-1), next(-1), deleted(false) {};
        RefNode(int value, int refPos) : val(value), pos(refPos), prev(-1), next(-1), deleted(false) {}
};

class PhraseNode
{
    public: 
        std::list<int> content;
        bool exp;
        int lnode;
        int rnode;
        int ltmp; // Needed for case where consecutive non exp phrase are replacing their edges in the same step
        int rtmp; // Needed for case where consecutive non exp phrase are replacing their edges in the same step
        PhraseNode* prev;
        PhraseNode* next;

        PhraseNode(int lnode, int rnode) : content({}), exp(false), lnode(lnode), rnode(rnode), ltmp(-1), rtmp(-1), prev(nullptr), next(nullptr) {}
        PhraseNode(std::list<int> ilist) : content(ilist), exp(true), lnode(-1), rnode(-1), ltmp(-1), rtmp(-1), prev(nullptr), next(nullptr) {}
};


#endif // NODE_H