// C++ Program to Implement Interval Tree
// Modified from Geeks for Geeks RB-Tree using C++ (https://www.geeksforgeeks.org/red-black-tree-in-cpp/)

#include <iostream>

// Enumeration for colors of nodes in Interval Tree (Red-Black Tree)
enum Color
{
    RED,
    BLACK
};

// Class template for IntervalRed-Black Tree
template <typename T>
class RBIntervalTree
{
private:
    // Structure for a node in IntervalTree (Red-Black Tree)
    struct Node
    {
        int low;
        int high;
        int max;
        T data;
        Color color;
        Node *parent;
        Node *left;
        Node *right;

        // Constructor to initialize node with low and color
        Node(std::pair<int, int> interval, T value)
            : low(interval.first), high(interval.second), max(interval.second), data(value), color(RED), parent(nullptr), left(nullptr), right(nullptr)
        {
        }
    };

    //Node *root; // Root of the Interval Tree (Red-Black Tree)

    // Utility function: Left Rotation
    void rotateLeft(Node *&node)
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
        updateMax(node);   // First update the original node
        updateMax(child);  // Then update the new subtree root
        propagateMax(child->parent);  // Ensure correctness for ancestors
    }

    // Utility function: Right Rotation
    void rotateRight(Node *&node)
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
        updateMax(node);   // First update the original node
        updateMax(child);  // Then update the new subtree root
        propagateMax(child->parent);  // Ensure correctness for ancestors
    }

    // Utility function: Fixing Insertion Violation
    void fixInsert(Node *&node)
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
    void fixDelete(Node *&node)
    {
        while (node != root && node->color == BLACK)
        {
            if (node == node->parent->left)
            {
                Node *sibling = node->parent->right;
                if (sibling->color == RED)
                {
                    sibling->color = BLACK;
                    node->parent->color = RED;
                    rotateLeft(node->parent);
                    sibling = node->parent->right;
                }
                if ((sibling->left == nullptr || sibling->left->color == BLACK) && (sibling->right == nullptr || sibling->right->color == BLACK))
                {
                    sibling->color = RED;
                    node = node->parent;
                }
                else
                {
                    if (sibling->right == nullptr || sibling->right->color == BLACK)
                    {
                        if (sibling->left != nullptr)
                            sibling->left->color = BLACK;
                        sibling->color = RED;
                        rotateRight(sibling);
                        sibling = node->parent->right;
                    }
                    sibling->color = node->parent->color;
                    node->parent->color = BLACK;
                    if (sibling->right != nullptr)
                        sibling->right->color = BLACK;
                    rotateLeft(node->parent);
                    node = root;
                }
            }
            else
            {
                Node *sibling = node->parent->left;
                if (sibling->color == RED)
                {
                    sibling->color = BLACK;
                    node->parent->color = RED;
                    rotateRight(node->parent);
                    sibling = node->parent->left;
                }
                if ((sibling->left == nullptr || sibling->left->color == BLACK) && (sibling->right == nullptr || sibling->right->color == BLACK))
                {
                    sibling->color = RED;
                    node = node->parent;
                }
                else
                {
                    if (sibling->left == nullptr || sibling->left->color == BLACK)
                    {
                        if (sibling->right != nullptr)
                            sibling->right->color = BLACK;
                        sibling->color = RED;
                        rotateLeft(sibling);
                        sibling = node->parent->left;
                    }
                    sibling->color = node->parent->color;
                    node->parent->color = BLACK;
                    if (sibling->left != nullptr)
                        sibling->left->color = BLACK;
                    rotateRight(node->parent);
                    node = root;
                }
            }
        }
        node->color = BLACK;
    }

    // Utility function: Find Node with Minimum Value
    Node *minValueNode(Node *&node)
    {
        Node *current = node;
        while (current->left != nullptr)
            current = current->left;
        return current;
    }

    // Utility function: Transplant nodes in IntervalTree (Red-Black Tree)
    void transplant(Node *&root, Node *&u, Node *&v)
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
    void printHelper(Node *root, std::string indent, bool last)
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
            std::cout << "[" << root->low << "," << root->high << "]" << " data: " << root->data << ", Max val: " << root->max << " (" << sColor << ")" << std::endl;
            printHelper(root->left, indent, false);
            printHelper(root->right, indent, true);
        }
    }

    // Utility function: Delete all nodes in the IntervalTree (Red-Black Tree)
    void deleteTree(Node *node)
    {
        if (node != nullptr)
        {
            deleteTree(node->left);
            deleteTree(node->right);
            delete node;
        }
    }

    // Utility Function: Update the max value after insertion or deletion
    void updateMax(Node *node)
    {
        if (node)
        {
            node->max = node->high;
            if (node->left)
                node->max = std::max(node->max, node->left->max);
            if (node->right)
                node->max = std::max(node->max, node->right->max);
        }
    }

    // Utility Function: Propagate max value to parent
    void propagateMax(Node *node)
    {
        while (node != nullptr)
        {
            updateMax(node);
            node = node->parent;
        }
    }

public:
    Node *root; // Root of the Interval Tree (Red-Black Tree)

    // Constructor: Initialize IntervalTree (Red-Black Tree)
    RBIntervalTree()
        : root(nullptr)
    {
    }

    // Destructor: Delete Interval Tree (Red-Black Tree)
    ~RBIntervalTree() { deleteTree(root); }

    // Public function: Insert a value into IntervalTree (Red-Black Tree)
    void insert(std::pair<int, int> interval, T data)
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
        propagateMax(temp);
    }

    // Public function: Remove a value from IntervalTree (Red-Black Tree)
    void remove(std::pair<int, int> interval)
    {
        Node *node = root;
        Node *z = nullptr;
        Node *x = nullptr;
        Node *y = nullptr;
        while (node != nullptr)
        {
            if (node->low == interval.first)
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
    }

    // Public function: Print the Interval Tree (Red-Black Tree)
    void printTree()
    {
        if (root == nullptr)
            std::cout << "Tree is empty." << std::endl;
        else
        {
            std::cout << "Interval (Red-Black) Tree:" << std::endl;
            printHelper(root, "", true);
        }
    }

    bool checkMaxProperty(Node *node)
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

        // Recursively check left and right subtrees
        return checkMaxProperty(node->left) && checkMaxProperty(node->right);
    }
};

// Driver program to test Red-Black Tree
int main()
{
    RBIntervalTree<int> rbtree;

    // Inserting values into Red-Black Tree
    for (int i = 1; i < 20; i++){
        rbtree.insert({i,i+1},i+2);
        rbtree.checkMaxProperty(rbtree.root);
    }

    rbtree.checkMaxProperty(rbtree.root);
    rbtree.printTree();
    
    // Deleting nodes from Red-Black Tree
    // std::cout << "After deleting 18:" << std::endl;
    //rbtree.remove(18);
    //rbtree.printTree();

    // std::cout << "After deleting 11:" << std::endl;
    // rbtree.remove(11);
    // rbtree.printTree();

    // std::cout << "After deleting 3:" << std::endl;
    // rbtree.remove(3);
    // rbtree.printTree();

    return 0;
}