#ifndef DOUBLE_LINKED_LIST_H
#define DOUBLE_LINKED_LIST_H

#include <iostream>
#include "node.h"

class RefLinkedList 
{
    private:
        RefNode* head;
        RefNode* tail;

    public:
        RefLinkedList() : head(nullptr), tail(nullptr) {}

        ~RefLinkedList() {
            while (head != nullptr) {
                RefNode* temp = head;
                head = head->next;
                delete temp;
            }
        }

        RefNode* getHead(){ return head; }
        RefNode* getTail(){ return tail; }

        // Insert at the front
        RefNode* push_front(int value) {
            RefNode* newNode = new RefNode(value);
            if (!head) {
                head = tail = newNode;
            } else {
                newNode->next = head;
                head->prev = newNode;
                head = newNode;
            }
            return newNode;
        }

        // Insert at the back
        RefNode* push_back(int value) {
            RefNode* newNode = new RefNode(value);
            if (!tail) {
                head = tail = newNode;
            } else {
                tail->next = newNode;
                newNode->prev = tail;
                tail = newNode;
            }
            return newNode;
        }

        // Replace pair and cleanup
        void replacePair(int val, RefNode* left, RefNode* right)
        {
            left->val = val;
            right->deleted = true;
            right->prev = left;
            if (!right->next){
                right->next->prev = left;
            }
        }

        // If the phrase endpoint on ref is deleted, 
        // then find the new ref endpoint
        RefNode* findNearestRef(RefNode* endpoint) 
        {
            if (!(endpoint->deleted))
                return endpoint;
            else{
                while(endpoint->deleted){
                    endpoint = endpoint->prev;
                }
                return endpoint;
            }
        }

        // If we need to find the forward adjacent RefNode
        RefNode* findForwardRef(RefNode* ref_elem)
        {
            ref_elem = ref_elem->next;
            while (ref_elem != nullptr && ref_elem->deleted){
                ref_elem = ref_elem->next;
            }
            return ref_elem;
        }

        void printForward() const {
            RefNode* temp = head;
            while (temp != nullptr) {
                std::cout << temp->val << " ";
                temp = temp->next;
            }
            std::cout << std::endl;
        }

        void printBackward() const {
            RefNode* temp = tail;
            while (temp != nullptr) {
                std::cout << temp->val << " ";
                temp = temp->prev;
            }
            std::cout << std::endl;
        }
};

class PhraseLinkedList 
{
    private:
        PhraseNode* head;
        PhraseNode* tail;

    public:
        PhraseLinkedList() : head(nullptr), tail(nullptr) {}

        ~PhraseLinkedList() {
            while (head != nullptr) {
                PhraseNode* temp = head;
                head = head->next;
                delete temp;
            }
        }

        PhraseNode* getHead(){ return head; }
        PhraseNode* getTail(){ return tail; }

        // Insert before specified PhraseNode
        PhraseNode* insert(PhraseNode* nextNode, std::list<unsigned int> ilist)
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
            delete currNode;
            return prevNode; // Return the next node to allow easy iteration
        }

        // Insert at the front
        PhraseNode* push_front(std::list<unsigned int> ilist) {
            PhraseNode* newNode = new PhraseNode(ilist);
            if (!head) {
                head = tail = newNode;
            } else {
                newNode->next = head;
                head->prev = newNode;
                head = newNode;
            }
            return newNode;
        }

        // Insert at the front
        PhraseNode* push_front(RefNode* lnode, RefNode* rnode, int lrange, int rrange) {
            PhraseNode* newNode = new PhraseNode(lnode, rnode, lrange, rrange);
            if (!head) {
                head = tail = newNode;
            } else {
                newNode->next = head;
                head->prev = newNode;
                head = newNode;
            }
            return newNode;
        }

        // Insert at the back
        PhraseNode* push_back(std::list<unsigned int> ilist) {
            PhraseNode* newNode = new PhraseNode(ilist);
            if (!tail) {
                head = tail = newNode;
            } else {
                tail->next = newNode;
                newNode->prev = tail;
                tail = newNode;
            }
            return newNode;
        }

         // Insert at the back
        PhraseNode* push_back(RefNode* lnode, RefNode* rnode, int lrange, int rrange) {
            PhraseNode* newNode = new PhraseNode(lnode, rnode, lrange, rrange);
            if (!tail) {
                head = tail = newNode;
            } else {
                tail->next = newNode;
                newNode->prev = tail;
                tail = newNode;
            }
            return newNode;
        }
};

#endif // DOUBLE_LINKED_LIST_H