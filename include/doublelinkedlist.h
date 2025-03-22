#ifndef DOUBLE_LINKED_LIST_H
#define DOUBLE_LINKED_LIST_H

#include <iostream>
#include <vector>
#include "node.h"

class RefLinkedList 
{
    private:
        int head;
        int tail;
        int size;

    public:
        std::vector<RefNode> nodes;

        RefLinkedList(int capacity) : size(0), head(-1), tail(-1){
            nodes.reserve(capacity);
        }

        // Get the head and tail of linked list array
        int getHead(){ return head; }
        int getTail(){ return tail; }
        // Get size of reference
        int getSize(){ return size; }

        // Insert at the back
        int push_back(int value) {
            if (size >= nodes.capacity()){ // Overflow guard
                std::cerr << "ERROR: Tried inserting more elements than allowed into the reference linked-list array!" << std::endl;
                return -1; 
            } 
            nodes.emplace_back(value); // Updates size but node pos zero-index
            if (tail == -1) {
                head = tail = 0;
            } else {
                nodes[tail].next = size;
                nodes[size].prev = tail;
                tail = size;
            }
            return size++; // Return the pos in the vector of the most recent element pushed back
        }

        // Replace pair and cleanup
        void replacePair(int val, int leftIdx, int rightIdx)
        {
            nodes[leftIdx].val = val;
            nodes[rightIdx].deleted = true;
            nodes[rightIdx].prev = leftIdx;   
            if (nodes[rightIdx].next != -1){
                nodes[leftIdx].next = nodes[rightIdx].next;
                nodes[nodes[rightIdx].next].prev = leftIdx;
            }
        }

        // If the phrase endpoint on ref is deleted, 
        // then find the new ref endpoint
        int findNearestRef(int idx) 
        {
            while (idx != -1 && nodes[idx].deleted){
                idx = nodes[idx].prev;
            }
            return idx;
        }

        // If we need to find the forward adjacent RefNode
        int findForwardRef(int idx)
        {
            idx = nodes[idx].next;
            while (idx != -1 && nodes[idx].deleted){
                idx = nodes[idx].next;
            }
            return idx;
        }

        void printForward() const {
            for (int i = head; i != -1; i = nodes[i].next) {
                std::cout << nodes[i].val << " ";
            }
            std::cout << std::endl;
        }

        void printBackward() const {
            for (int i = tail; i != -1; i = nodes[i].prev) {
                std::cout << nodes[i].val << " ";
            }
            std::cout << std::endl;
        }
};

class PhraseLinkedList 
{
    private:
        int head;
        int tail;
        int size;

    public:
        std::vector<PhraseNode> phrases;

        PhraseLinkedList(int capacity) : head(-1), tail(-1), size(0) {
            phrases.reserve(2 * capacity);
        }

        int getHead(){ return head; }
        int getTail(){ return tail; }
        int getSize(){ return size; }

        // Insert before specified PhraseNode (using index)
        int insert(int nextNodeIndex, std::list<int> ilist)
        {
            if (nextNodeIndex == -1) return -1; // Invalid node index
            PhraseNode newNode(ilist);
            phrases.emplace_back(newNode);
            int newNodeIndex = phrases.size() - 1;

            // Setting new node's prev and next
            phrases[newNodeIndex].next = nextNodeIndex;
            phrases[newNodeIndex].prev = phrases[nextNodeIndex].prev;

            // Updating the previous node's next pointer
            if (phrases[nextNodeIndex].prev != -1) {
                phrases[phrases[nextNodeIndex].prev].next = newNodeIndex;
            } else {
                head = newNodeIndex;  // Update head if the new node is inserted at the front
            }

            // Updating the next node's prev pointer
            phrases[nextNodeIndex].prev = newNodeIndex;

            size++;
            return newNodeIndex;
        }

        // Lazy Delete at specified PhraseNode (using index)
        int remove(int currNodeIndex) 
        {
            if (currNodeIndex == -1) return -1; // Invalid node index

            int prevNodeIndex = phrases[currNodeIndex].prev;
            int nextNodeIndex = phrases[currNodeIndex].next;

            if (prevNodeIndex != -1) {
                phrases[prevNodeIndex].next = nextNodeIndex;
            } else {
                head = nextNodeIndex; // Update head if removing the first node
            }

            if (nextNodeIndex != -1) {
                phrases[nextNodeIndex].prev = prevNodeIndex;
            } else {
                tail = prevNodeIndex; // Update tail if removing the last node
            }

            //phrases.erase(phrases.begin() + currNodeIndex); // 
            //size--;

            // Return the index of the previous node or the next node if the previous is -1
            return (prevNodeIndex != -1) ? prevNodeIndex : nextNodeIndex;
        }

        // Insert at the back using list<int>
        int push_back(std::list<int> ilist) 
        {
            PhraseNode newNode(ilist);
            phrases.emplace_back(newNode);
            int newNodeIndex = phrases.size() - 1;
            
            if (tail == -1) {
                head = tail = newNodeIndex;
            } else {
                phrases[newNodeIndex].prev = tail;
                phrases[tail].next = newNodeIndex;
                tail = newNodeIndex;
            }

            size++;
            return newNodeIndex;
        }

        // Insert at the back using lnode, rnode
        int push_back(int lnode, int rnode) 
        {
            PhraseNode newNode(lnode, rnode);
            phrases.push_back(newNode);
            int newNodeIndex = phrases.size() - 1;

            if (tail == -1) {
                head = tail = newNodeIndex;
            } else {
                phrases[newNodeIndex].prev = tail;
                phrases[tail].next = newNodeIndex;
                tail = newNodeIndex;
            }

            size++;
            return newNodeIndex;
        }
};

#endif // DOUBLE_LINKED_LIST_H