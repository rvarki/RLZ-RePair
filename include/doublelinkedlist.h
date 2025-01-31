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

        // Insert at the front
        PhraseNode* push_front(std::string value) {
            PhraseNode* newNode = new PhraseNode(value);
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
        PhraseNode* push_front(RefNode* lrange, RefNode* rrange) {
            PhraseNode* newNode = new PhraseNode(lrange, rrange);
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
        PhraseNode* push_back(std::string value) {
            PhraseNode* newNode = new PhraseNode(value);
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
        PhraseNode* push_back(RefNode* lrange, RefNode* rrange) {
            PhraseNode* newNode = new PhraseNode(lrange, rrange);
            if (!tail) {
                head = tail = newNode;
            } else {
                tail->next = newNode;
                newNode->prev = tail;
                tail = newNode;
            }
            return newNode;
        }

        void printForward() const {
            PhraseNode* temp = head;
            while (temp != nullptr) {
                std::cout << temp->content << " ";
                temp = temp->next;
            }
            std::cout << std::endl;
        }

        void printBackward() const {
            PhraseNode* temp = tail;
            while (temp != nullptr) {
                std::cout << temp->content << " ";
                temp = temp->prev;
            }
            std::cout << std::endl;
        }
};

#endif // DOUBLE_LINKED_LIST_H