/*
* RLZ-RePair - RePair compression using RLZ parse
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

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
        PhraseNode* head;
        PhraseNode* tail;
        int size;

    public:
        PhraseLinkedList() : size(0), head(nullptr), tail(nullptr) {}

        ~PhraseLinkedList() {
            while (head != nullptr) {
                PhraseNode* temp = head;
                head = head->next;
                delete temp;
            }
        }

        PhraseNode* getHead(){ return head; }
        PhraseNode* getTail(){ return tail; }
        int getSize(){ return size; }

        // Insert before specified PhraseNode
        PhraseNode* insert(PhraseNode* nextNode, std::list<int> ilist)
        {
            if (!nextNode) return nullptr; // Invalid node
            PhraseNode* newNode = new PhraseNode(ilist);
            newNode->next = nextNode;
            newNode->prev = nextNode->prev;
            if (nextNode->prev) {
                nextNode->prev->next = newNode;
            } else {
                head = newNode; 
            }
            nextNode->prev = newNode;
            size++;
            return newNode;
        }

        // Delete at specified PhraseNode
        PhraseNode* remove(PhraseNode* currNode) 
        {
            if (!currNode) return nullptr; // Invalid node

            if (currNode->prev) {
                currNode->prev->next = currNode->next;
            } else {
                head = currNode->next; // Update head if removing the first node
            }

            if (currNode->next) {
                currNode->next->prev = currNode->prev;
            } else {
                tail = currNode->prev; // Update tail if removing the last node
            }

            PhraseNode* prevNode = currNode->prev; // Store previous node before deleting
            PhraseNode* forwardNode = currNode->next; // Store previous node before deleting
            size--;
            delete currNode;
            if (prevNode != nullptr)
                return prevNode; // Return the prev node to allow easy iteration
            else
                return forwardNode; // If prev node is nullptr then give the next node
        }

        // Insert at the front
        PhraseNode* push_front(std::list<int> ilist) {
            PhraseNode* newNode = new PhraseNode(ilist);
            if (!head) {
                head = tail = newNode;
            } else {
                newNode->next = head;
                head->prev = newNode;
                head = newNode;
            }
            size++;
            return newNode;
        }

        // Insert at the front
        PhraseNode* push_front(int lnode, int rnode) {
            PhraseNode* newNode = new PhraseNode(lnode, rnode);
            if (!head) {
                head = tail = newNode;
            } else {
                newNode->next = head;
                head->prev = newNode;
                head = newNode;
            }
            size++;
            return newNode;
        }

        // Insert at the back
        PhraseNode* push_back(std::list<int> ilist) {
            PhraseNode* newNode = new PhraseNode(ilist);
            if (!tail) {
                head = tail = newNode;
            } else {
                tail->next = newNode;
                newNode->prev = tail;
                tail = newNode;
            }
            size++;
            return newNode;
        }

         // Insert at the back
        PhraseNode* push_back(int lnode, int rnode) {
            PhraseNode* newNode = new PhraseNode(lnode, rnode);
            if (!tail) {
                head = tail = newNode;
            } else {
                tail->next = newNode;
                newNode->prev = tail;
                tail = newNode;
            }
            size++;
            return newNode;
        }
};

#endif // DOUBLE_LINKED_LIST_H