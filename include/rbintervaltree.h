#ifndef RB_INTERVAL_TREE_H
#define RB_INTERVAL_TREE_H

// C++ Program to Implement Interval Tree
// Modified from Geeks for Geeks RB-Tree using C++ (https://www.geeksforgeeks.org/red-black-tree-in-cpp/)

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

// Enumeration for colors of nodes in Interval Tree (Red-Black Tree)
enum Color{ RED, BLACK};

// Class template for IntervalRed-Black Tree
template <typename T>
class RBIntervalTree
{
    private:
        // Structure for a node in IntervalTree (Red-Black Tree)
        struct Node
        {
            int low, high, min, max;
            T data;
            Color color;
            Node *parent, *left, *right;
            Node(std::pair<int, int> interval, T value);
        };

        // Private members
        Node* root;
        std::vector<T> returnValues;

        // Private member function declarations
        void rotateLeft(Node*& node);
        void rotateRight(Node*& node);
        void fixInsert(Node*& node);
        void fixDelete(Node*& node);
        Node* minValueNode(Node*& node);
        void transplant(Node*& root, Node*& u, Node*& v);
        void printHelper(Node* root, std::string indent, bool last);
        void deleteTree(Node* node);
        void update(Node* node);
        void propagate(Node* node);
        void checkNode(std::pair<int, int>& spair, Node*& node);
    
    public:
        RBIntervalTree();
        ~RBIntervalTree();
        void insert(std::pair<int, int> interval, T data);
        void remove(std::pair<int, int> interval, T data);
        void clear();
        void printTree();
        std::vector<T> findContained(std::pair<int, int> spair);
        bool checkProperty(Node* node);
};

template <typename T>
RBIntervalTree<T>::Node::Node(std::pair<int, int> interval, T value)
        : low(interval.first), high(interval.second), min(interval.first), max(interval.second), data(value), color(RED), parent(nullptr), left(nullptr), right(nullptr){}

// Utility function: Left Rotation
template <typename T>
void RBIntervalTree<T>::rotateLeft(typename RBIntervalTree<T>::Node *&node)
{
    Node *child = node->right;
    node->right = child->left;
    if (node->right != nullptr)
        node->right->parent = node;
    child->parent = node->parent;
    if (node->parent == nullptr)
        root = child;
    else if (node == node->parent->left)
        node->parent->left = child;
    else
        node->parent->right = child;
    child->left = node;
    node->parent = child;
    // Fix max values
    update(node);   // First update the original node
    update(child);  // Then update the new subtree root
    propagate(child->parent);  // Ensure correctness for ancestors
}

// Utility function: Right Rotation
template <typename T>
void RBIntervalTree<T>::rotateRight(typename RBIntervalTree<T>::Node *&node)
{
    Node *child = node->left;
    node->left = child->right;
    if (node->left != nullptr)
        node->left->parent = node;
    child->parent = node->parent;
    if (node->parent == nullptr)
        root = child;
    else if (node == node->parent->left)
        node->parent->left = child;
    else
        node->parent->right = child;
    child->right = node;
    node->parent = child;
    update(node);   // First update the original node
    update(child);  // Then update the new subtree root
    propagate(child->parent);  // Ensure correctness for ancestors
}

// Utility function: Fixing Insertion Violation
template <typename T>
void RBIntervalTree<T>::fixInsert(typename RBIntervalTree<T>::Node *&node)
{
    Node *parent = nullptr;
    Node *grandparent = nullptr;
    while (node != root && node->color == RED && node->parent->color == RED)
    {
        parent = node->parent;
        grandparent = parent->parent;
        if (parent == grandparent->left)
        {
            Node *uncle = grandparent->right;
            // Case 1: Uncle is RED - Recoloring
            if (uncle != nullptr && uncle->color == RED)
            {
                grandparent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandparent;
            }
            else
            {
                // Case 2: Triangle (Left-Right) - Rotate Left at Parent
                if (node == parent->right)
                {
                    rotateLeft(parent);
                    node = parent;
                    parent = node->parent;
                }
                // Case 3: Line (Left-Left) - Rotate Right at Grandparent
                rotateRight(grandparent);
                std::swap(parent->color, grandparent->color);
                node = parent;
            }
        }
        else
        {
            Node *uncle = grandparent->left;
            // Case 1: Uncle is RED - Recoloring
            if (uncle != nullptr && uncle->color == RED)
            {
                grandparent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandparent;
            }
            else
            {
                // Case 2: Triangle (Right-Left) - Rotate Right at Parent
                if (node == parent->left)
                {
                    rotateRight(parent);
                    node = parent;
                    parent = node->parent;
                }
                // Case 3: Line (Right-Right) - Rotate Left at Grandparent
                rotateLeft(grandparent);
                std::swap(parent->color, grandparent->color);
                node = parent;
            }
        }
    }
    root->color = BLACK;
}

// Utility function: Fixing Deletion Violation
template <typename T>
void RBIntervalTree<T>::fixDelete(typename RBIntervalTree<T>::Node *&node)
{
    while (node != root && node != nullptr && node->color == BLACK)
    {
        if (node->parent != nullptr && node == node->parent->left)
        {
            Node *sibling = node->parent->right;
            
            if (sibling != nullptr && sibling->color == RED)
            {
                sibling->color = BLACK;
                node->parent->color = RED;
                rotateLeft(node->parent);
                sibling = node->parent->right;
            }

            if (sibling != nullptr &&
                (sibling->left == nullptr || sibling->left->color == BLACK) &&
                (sibling->right == nullptr || sibling->right->color == BLACK))
            {
                sibling->color = RED;
                node = node->parent;
            }
            else
            {
                if (sibling != nullptr && (sibling->right == nullptr || sibling->right->color == BLACK))
                {
                    if (sibling->left != nullptr)
                        sibling->left->color = BLACK;
                    sibling->color = RED;
                    rotateRight(sibling);
                    sibling = node->parent->right;
                }

                if (sibling != nullptr)
                {
                    sibling->color = node->parent->color;
                    node->parent->color = BLACK;
                    if (sibling->right != nullptr)
                        sibling->right->color = BLACK;
                    rotateLeft(node->parent);
                }
                node = root;
            }
        }
        else
        {
            Node *sibling = node->parent->left;
            
            if (sibling != nullptr && sibling->color == RED)
            {
                sibling->color = BLACK;
                node->parent->color = RED;
                rotateRight(node->parent);
                sibling = node->parent->left;
            }

            if (sibling != nullptr &&
                (sibling->left == nullptr || sibling->left->color == BLACK) &&
                (sibling->right == nullptr || sibling->right->color == BLACK))
            {
                sibling->color = RED;
                node = node->parent;
            }
            else
            {
                if (sibling != nullptr && (sibling->left == nullptr || sibling->left->color == BLACK))
                {
                    if (sibling->right != nullptr)
                        sibling->right->color = BLACK;
                    sibling->color = RED;
                    rotateLeft(sibling);
                    sibling = node->parent->left;
                }

                if (sibling != nullptr)
                {
                    sibling->color = node->parent->color;
                    node->parent->color = BLACK;
                    if (sibling->left != nullptr)
                        sibling->left->color = BLACK;
                    rotateRight(node->parent);
                }
                node = root;
            }
        }
    }

    if (node != nullptr)
        node->color = BLACK;  // Ensure the final node is black
}

// Utility function: Find Node with Minimum Value
template <typename T>
typename RBIntervalTree<T>::Node* RBIntervalTree<T>::minValueNode(typename RBIntervalTree<T>::Node *&node)
{
    Node *current = node;
    while (current->left != nullptr)
        current = current->left;
    return current;
}

// Utility function: Transplant nodes in IntervalTree (Red-Black Tree)
template <typename T>
void RBIntervalTree<T>::transplant(typename RBIntervalTree<T>::Node *&root, typename RBIntervalTree<T>::Node *&u, typename RBIntervalTree<T>::Node *&v)
{
    if (u->parent == nullptr)
        root = v;
    else if (u == u->parent->left)
        u->parent->left = v;
    else
        u->parent->right = v;
    if (v != nullptr)
        v->parent = u->parent;
}

// Utility function: Helper to print IntervalTree (Red-Black Tree)
template <typename T>
void RBIntervalTree<T>::printHelper(typename RBIntervalTree<T>::Node *root, std::string indent, bool last)
{
    if (root != nullptr)
    {
        std::cout << indent;
        if (last)
        {
            std::cout << "R----";
            indent += "   ";
        }
        else
        {
            std::cout << "L----";
            indent += "|  ";
        }
        std::string sColor = (root->color == RED) ? "RED" : "BLACK";
        std::cout << "[" << root->low << "," << root->high << "]" << " data: " << root->data << ", Max val: " << root->max << ", Min val: " << root->min << " (" << sColor << ")" << std::endl;
        printHelper(root->left, indent, false);
        printHelper(root->right, indent, true);
    }
}

// Utility function: Delete all nodes in the IntervalTree (Red-Black Tree)
template <typename T>
void RBIntervalTree<T>::deleteTree(typename RBIntervalTree<T>::Node *node)
{
    if (node != nullptr)
    {
        deleteTree(node->left);
        deleteTree(node->right);
        delete node;
    }
}

// Utility Function: Update the max and min value after insertion or deletion
template <typename T>
void RBIntervalTree<T>::update(typename RBIntervalTree<T>::Node *node)
{
    if (node)
    {
        node->max = node->high;
        node->min = node->low;
        if (node->left){
            node->max = std::max(node->max, node->left->max);
            node->min = std::min(node->min, node->left->min);
        }
        if (node->right){
            node->max = std::max(node->max, node->right->max);
            node->min = std::min(node->min, node->right->min);
        }
    }
}

// Utility Function: Propagate max and min value to parent
template <typename T>
void RBIntervalTree<T>::propagate(typename RBIntervalTree<T>::Node *node)
{
    while (node != nullptr)
    {
        update(node);
        node = node->parent;
    }
}

// Utility Function for findContained
template <typename T>
void RBIntervalTree<T>::checkNode(std::pair<int,int>& spair, typename RBIntervalTree<T>::Node*& node)
{
    if (spair.first >= node->low && spair.second <= node->high){
        returnValues.emplace_back(node->data);
    }
    if (node->left != nullptr && node->left->min <= spair.first && node->left->max >= spair.second){
        checkNode(spair, node->left);
    }
    if (node->right != nullptr && node->right->min <= spair.first && node->right->max >= spair.second){
        checkNode(spair, node->right);
    }
}

// Constructor: Initialize IntervalTree (Red-Black Tree)
template <typename T>
RBIntervalTree<T>::RBIntervalTree() : root(nullptr){}

// Destructor: Delete Interval Tree (Red-Black Tree)
template <typename T>
RBIntervalTree<T>::~RBIntervalTree() { deleteTree(root); }

// Public function: Clear the tree.
template <typename T>
void RBIntervalTree<T>::clear() {
    deleteTree(root);  // Call the helper function to recursively delete nodes
    root = nullptr;    // After deletion, set root to nullptr
}

// Public function: Insert a value into IntervalTree (Red-Black Tree)
template <typename T>
void RBIntervalTree<T>::insert(std::pair<int, int> interval, T data)
{
    Node *node = new Node(interval, data);
    Node* temp = node;
    Node *parent = nullptr;
    Node *current = root;
    while (current != nullptr)
    {
        parent = current;
        if (node->low < current->low)
            current = current->left;
        else
            current = current->right;
    }
    node->parent = parent;
    if (parent == nullptr)
        root = node;
    else if (node->low < parent->low)
        parent->left = node;
    else
        parent->right = node;
    
    fixInsert(node);
    propagate(temp);
}

// Public function: Remove a value from IntervalTree (Red-Black Tree)
template <typename T>
void RBIntervalTree<T>::remove(std::pair<int, int> interval, T data)
{
    Node *node = root;
    Node *z = nullptr;
    Node *x = nullptr;
    Node *y = nullptr;
    while (node != nullptr)
    {
        if (node->low == interval.first && node->high == interval.second && node->data == data)
        {
            z = node;
        }

        if (node->low <= interval.first)
        {
            node = node->right;
        }
        else
        {
            node = node->left;
        }
    }

    if (z == nullptr)
    {
        std::cout << "Key not found in the tree" << std::endl;
        return;
    }

    y = z;
    Color yOriginalColor = y->color;
    if (z->left == nullptr)
    {
        x = z->right;
        transplant(root, z, z->right);
    }
    else if (z->right == nullptr)
    {
        x = z->left;
        transplant(root, z, z->left);
    }
    else
    {
        y = minValueNode(z->right);
        yOriginalColor = y->color;
        x = y->right;
        if (y->parent == z)
        {
            if (x != nullptr)
                x->parent = y;
        }
        else
        {
            transplant(root, y, y->right);
            y->right = z->right;
            y->right->parent = y;
        }
        transplant(root, z, y);
        y->left = z->left;
        y->left->parent = y;
        y->color = z->color;
    }
    delete z;
    if (yOriginalColor == BLACK)
    {
        fixDelete(x);
    }
    Node *curr = (y != nullptr) ? y : x;
    propagate(curr);
}

// Public function: Print the Interval Tree (Red-Black Tree)
template <typename T>
void RBIntervalTree<T>::printTree()
{
    if (root == nullptr)
        std::cout << "Tree is empty." << std::endl;
    else
    {
        std::cout << "Interval (Red-Black) Tree:" << std::endl;
        printHelper(root, "", true);
    }
}

// Public function: Find all intervals containing pair
template <typename T>    
std::vector<T> RBIntervalTree<T>::findContained(std::pair<int,int> spair)
{
    returnValues.clear();
    Node* node = root;
    checkNode(spair, root);
    return returnValues;
}

// Public function: Checks if min and max property of tree are maintained
template <typename T> 
bool RBIntervalTree<T>::checkProperty(typename RBIntervalTree<T>::Node *node)
{
    if (node == nullptr)
        return true;

    // Compute expected max value
    int expectedMax = node->high;
    if (node->left)
        expectedMax = std::max(expectedMax, node->left->max);
    if (node->right)
        expectedMax = std::max(expectedMax, node->right->max);

    // Check if the max property holds for the current node
    if (node->max != expectedMax)
    {
        std::cout << "Max property violated at node [" << node->low << ", " << node->high << "]"
                << " Expected: " << expectedMax << ", Found: " << node->max << std::endl;
        return false;
    }

    // Compute expected min value
    int expectedMin = node->low;
    if (node->left)
        expectedMin = std::min(expectedMin, node->left->min);
    if (node->right)
        expectedMin = std::min(expectedMin, node->right->min);

    // Check if the min property holds for the current node
    if (node->min != expectedMin)
    {
        std::cout << "Min property violated at node [" << node->low << ", " << node->high << "]"
                << " Expected: " << expectedMin << ", Found: " << node->min << std::endl;
        return false;
    }

    // Recursively check left and right subtrees
    return checkProperty(node->left) && checkProperty(node->right);
}

#endif // RB_INTERVAL_TREE_H