#ifndef NODE_H
#define NODE_H
#include <list>

class RefNode
{
    public:
        int val;
        int prev, next;
        bool deleted;
        RefNode() : val(-1), prev(-1), next(-1), deleted(false) {};
        RefNode(int value) : val(value), prev(-1), next(-1), deleted(false) {}
};

class PhraseNode
{
    public: 
        std::list<int> content;
        bool exp;
        int lnode;
        int rnode;
        bool ltmp; // Needed for case where consecutive non exp phrase are replacing their edges in the same step
        bool rtmp; // Needed for case where consecutive non exp phrase are replacing their edges in the same step
        int prev;
        int next;

        PhraseNode(int lnode, int rnode) : content({}), exp(false), lnode(lnode), rnode(rnode), ltmp(false), rtmp(false), prev(-1), next(-1) {}
        PhraseNode(std::list<int> ilist) : content(ilist), exp(true), lnode(-1), rnode(-1), ltmp(false), rtmp(false), prev(-1), next(-1) {}
};


#endif // NODE_H