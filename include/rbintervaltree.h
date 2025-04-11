/*
* RLZ-RePair - RePair compression using RLZ parse
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef RB_INTERVAL_TREE_H
#define RB_INTERVAL_TREE_H

// C++ Program to Implement Interval Tree
// Modified from Geeks for Geeks RB-Tree using C++ (https://www.geeksforgeeks.org/red-black-tree-in-cpp/)

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>

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
            int low, min, max;
            std::vector<int> highVec;
            std::vector<T> dataVec;
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
        bool isBST(Node* root, std::string message, Node* minNode, Node* maxNode);
        bool validateRB(Node* root, int blackHeight);
        bool checkProperty(Node* node);
  
    public:
        RBIntervalTree();
        ~RBIntervalTree();
        void insert(std::pair<int, int> interval, T data);
        void remove(std::pair<int, int> interval, T data);
        void clear();
        void printTree();
        std::vector<T> findContained(std::pair<int, int> spair);
        bool isValidIT();
        bool isValidRB();
};

template <typename T>
RBIntervalTree<T>::Node::Node(std::pair<int, int> interval, T value)
    : low(interval.first), highVec{interval.second}, min(interval.first), max(interval.second), dataVec{value}, color(RED), parent(nullptr), left(nullptr), right(nullptr){}

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
    // Fix min and max values
    update(node->right);   // First update the right child of the original node since it is swapped.
    propagate(node);  // Propagate node since left child same and right child updated
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
    // Fix min and max values
    update(node->left);   // First update the left child of the original node since it is swapped.
    propagate(node);  // Propagate node since right child same and left child updated
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

        std::string highValues = "";
        for (int i = 0; i < root->highVec.size(); i++){
            highValues += std::to_string(root->highVec[i]);
            if (i != root->highVec.size() - 1){
                highValues += ",";
            }
        }

        std::cout << "[" << root->low << "," << "(" << highValues << ")" << "]" << " Min val: " << root->min << ", Max val:" << root->max << " (" << sColor << ")" << std::endl;
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
        node->max = *std::max_element(node->highVec.begin(), node->highVec.end());
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
    if (spair.first >= node->low){
        // 'it' is an iterator to high vals in highVec  
        for (auto it = node->highVec.begin(); it != node->highVec.end(); it++)
        {
            if (spair.second <= *it){
                int dataPos = std::distance(node->highVec.begin(), it);
                returnValues.emplace_back(node->dataVec[dataPos]);
            }
        }
    }
    if (node->left != nullptr && node->left->min <= spair.first && node->left->max >= spair.second){
        checkNode(spair, node->left);
    }
    if (node->right != nullptr && node->right->min <= spair.first && node->right->max >= spair.second){
        checkNode(spair, node->right);
    }
}

// Utility Function: Checks if the Interval Tree (Red-Black Tree) is a valid BST 
template <typename T>
bool RBIntervalTree<T>::isBST(typename RBIntervalTree<T>::Node* node, std::string message, typename RBIntervalTree<T>::Node* minNode, typename RBIntervalTree<T>::Node* maxNode)
{
    if (!node) return true;

    if ((minNode != nullptr && node->low < minNode->low) || (maxNode != nullptr && node->low >= maxNode->low)){
        std::cout << "BST Violated" << std::endl;
        return false;
    } 

    return isBST(node->left, "left", minNode, node) && isBST(node->right, "right", node, maxNode);
}

// Utility Function: Checks Red-Black properties of the Interval Tree (Red-Black Tree)
template <typename T>
bool RBIntervalTree<T>::validateRB(typename RBIntervalTree<T>::Node* node, int blackHeight)
{
    // If empty tree then automatically valid
    if (node == nullptr){
        blackHeight = 1;
        return true;
    }
    // Checks if root is black
    if (node->parent == nullptr && node->color != BLACK){
        std::cout << "ROOT IS NOT BLACK!" << std::endl;
    }
    // Checks if there are two consecutive red nodes
    if (node->color == RED){
        if (node->left != nullptr && node->left->color == RED || node->right != nullptr && node->right->color == RED){
            std::cout << "RED NODE HAS A RED CHILD" << std::endl;
            return false;
        }
    }
    // Keep track of the number of black nodes in the left and right subtrees
    int leftBlackHeight = 0;
    int rightBlackHeight = 0;

    if (!validateRB(node->left, leftBlackHeight) || !validateRB(node->right, rightBlackHeight)){
        return false;
    }

    // Ensure all paths have the same black height
    if (leftBlackHeight != rightBlackHeight) {
        std::cout << "BLACK HEIGHT MISMATCH!" << std::endl;
        return false;
    }
    // Compute the black height for the current node
    blackHeight = leftBlackHeight + (node->color == BLACK ? 1 : 0);

    return true;
}

// Public function: Checks if min and max property of tree are maintained
template <typename T> 
bool RBIntervalTree<T>::checkProperty(typename RBIntervalTree<T>::Node *node)
{
    if (node == nullptr)
        return true;

    // Compute expected max value
    int expectedMax = *std::max_element(node->highVec.begin(), node->highVec.end());
    if (node->left)
        expectedMax = std::max(expectedMax, node->left->max);
    if (node->right)
        expectedMax = std::max(expectedMax, node->right->max);

    // Check if the max property holds for the current node
    if (node->max != expectedMax)
    {
        std::string highValues = "";
        for (int i = 0; i < node->highVec.size(); i++){
            highValues += std::to_string(node->highVec[i]);
            if (i != node->highVec.size() - 1){
                highValues += ",";
            }
        }
        std::cout << "Max property violated at node [" << node->low << ", (" << highValues << ")]"
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
        std::string highValues = "";
        for (int i = 0; i < node->highVec.size(); i++){
            highValues += std::to_string(node->highVec[i]);
            if (i != node->highVec.size() - 1){
                highValues += ",";
            }
        }
        std::cout << "Min property violated at node [" << node->low << ",(" << highValues << ")]"
                << " Expected: " << expectedMin << ", Found: " << node->min << std::endl;
        return false;
    }

    // Recursively check left and right subtrees
    return checkProperty(node->left) && checkProperty(node->right);
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
    bool exist = false;
    Node *parent = nullptr;
    Node *current = root;
    while (current != nullptr)
    {
        parent = current;
        // The low value of the interval could exist in the red black tree already
        if (interval.first == current->low){
            exist = true;
            current->highVec.emplace_back(interval.second);
            current->dataVec.emplace_back(data);
            break;
        }
        else if (interval.first < current->low){
            current = current->left;
        }
        else{
            current = current->right;
        }
    }
    // If the low value of the interval did not exist in the red black tree already 
    if (!exist)
    {
        Node *node = new Node(interval, data);
        Node* temp = node;
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
    else{
        propagate(current);
    }
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
        // If node low value is present 
        if (node->low == interval.first)
        {
            // Get the iterator of the data value if it exists
            auto dataIt = std::find(node->dataVec.begin(), node->dataVec.end(), data);
            if (dataIt != node->dataVec.end())
            {
                // Calculate the position in the vector
                int dataPos = std::distance(node->dataVec.begin(), dataIt);
                // The high value should exist in this position if the other two conditions are true
                if (interval.second == node->highVec[dataPos]){
                    // Erase the value from both highVec and dataVec
                    node->highVec.erase(node->highVec.begin() + dataPos);
                    node->dataVec.erase(dataIt);
                    if (node->highVec.size() != node->dataVec.size()){
                        std::cout << "Node content is out of sync" << std::endl;
                    }
                    // If size is zero then delete from the tree
                    if (node->highVec.size() == 0 && node->dataVec.size() == 0){
                        z = node;
                        break;
                    }
                    // Update the tree if needed and exit
                    else{
                        propagate(node);
                        return;
                    }
                }
                else{
                    std::cout << "Something is wrong. Data is present but not the corresponding high value" << std::endl;
                }
            }
            else{
                std::cout << "Data does not exist in tree" << std::endl;
            }
        }

        // If node value is less than or equal then go right.
        else if (node->low < interval.first)
        {
            node = node->right;
        }
        // Otherwise go left.
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
        if (z->right == nullptr) // If z->right does not exist then propagate from z->parent itself
            propagate(z->parent);
        else
            propagate(z->right); // Have only reattached z->right to parent of z.
    }
    else if (z->right == nullptr)
    {
        x = z->left;
        transplant(root, z, z->left);
        if (z->left == nullptr)
            propagate(z->parent); // If z->left does not exist then propagate from z->parent itself
        else
            propagate(z->left); // Have only reattached z->left to parent of z.
    }
    else
    {
        y = minValueNode(z->right);
        yOriginalColor = y->color;
        x = y->right;
        bool yParentz = false;
        if (y->parent == z) // z->right does not have left subtree
        {
            yParentz = true;
            if (x != nullptr)
                x->parent = y;
        }
        else
        {
            transplant(root, y, y->right);
            y->right = z->right;
            y->right->parent = y;
        }
        Node* y_orig_parent = y->parent;
        transplant(root, z, y);
        y->left = z->left;
        y->left->parent = y;
        y->color = z->color;
        if (yParentz){
            propagate(y);
        }
        else{ 
            if (x == nullptr) // If y moves far and x is nullptr then need to call propagate on y's original parent.
                propagate(y_orig_parent);
            else
                propagate(x);
        }
    }
    delete z;
    if (yOriginalColor == BLACK)
    {
        fixDelete(x);
    }
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
    if (node == nullptr){
        return returnValues;
    }
    checkNode(spair, root);
    return returnValues;
}

// Public function: Check if Red Black tree is valid
template <typename T>
bool RBIntervalTree<T>::isValidRB()
{
    int blackHeight = 0;
    return isBST(root, "Root", nullptr, nullptr) && validateRB(root, blackHeight);
}

// Public function: Checks if the Red Black tree is valid Interval Tree 
template <typename T>
bool RBIntervalTree<T>::isValidIT()
{
    bool valid;
    valid = checkProperty(root);
    if (!valid){
        std::cout << "Min and Max Ranges in Subtree Not Correct" << std::endl;
    }
    return valid;
}


#endif // RB_INTERVAL_TREE_H