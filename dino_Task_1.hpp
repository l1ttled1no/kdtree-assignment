#include "main.hpp"
// #include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct Node
{
    vector<int> data;
    Node *left;
    Node *right;
    Node(vector<int> data, Node *left = nullptr, Node *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }

    void print() const
    {
        OUTPUT << "(";
        for(int i = 0; i < data.size(); i++)
        {
            OUTPUT << data[i];
            if (i == data.size() - 1) OUTPUT << ")";
            else OUTPUT << ",";
        }
        OUTPUT << " ";
    }
};

class kDTree
{
private:
    int k;
    Node *root;
    int count;
private:
    // hàm phụ nếu cần
public:
    kDTree(int k = 2) 
    {
        this->k = k; 
        this->root = nullptr;
        this->count = 0;
    }
    //recursive destructor
    void recDestructor(Node *root){
        if (root == nullptr) return;
        recDestructor(root->left);
        recDestructor(root->right);
        delete root;
    }
    ~kDTree()
    {
        recDestructor(root);
    }
    Node *copy(const Node *root)
    {
        if (root == nullptr) return nullptr;
        Node *node = new Node(root->data);
        node->left = copy(root->left);
        node->right = copy(root->right);
        return node;
    }

    const kDTree &operator=(const kDTree &other){
        // Deep copy
        if (this == &other) return *this;
        this->k = other.k;
        this->root = this->copy(other.root);
        this->count = other.count;
        return *this;
    }
    kDTree(const kDTree &other) {
        //deep copy
        this->k = other.k;
        this->root = this->copy(other.root);
        this->count = other.count;
    }
    int nodeCount() const{
        return count;
    }
    //recursive height
    int recHeight(Node *root) const{
        if (root == nullptr) return 0;
        return 1 + max(recHeight(root->left), recHeight(root->right));
    }
    int height() const{
        return recHeight(root);
    }
    //recursive leafCount
    int recLeafCount(Node *root) const{
        if (root == nullptr) return 0;
        if (root->left == nullptr && root->right == nullptr) return 1;
        return recLeafCount(root->left) + recLeafCount(root->right);
    }
    int leafCount() const{
        return recLeafCount(root);
    }

    void recInorder(Node *root) const{
        if (root == nullptr) return;
        recInorder(root->left);
        root->print();
        recInorder(root->right);
    }
    void inorderTraversal() const{
        this->recInorder(root);
    }
    void recPreorder(Node *root) const{
        if (root == nullptr) return;
        root->print();
        recPreorder(root->left);
        recPreorder(root->right);
    }
    void preorderTraversal() const{
        this->recPreorder(root);
    }
    void recPostorder(Node *root) const{
        if (root == nullptr) return;
        recPostorder(root->left);
        recPostorder(root->right);
        root->print();
    }
    void postorderTraversal() const{
        this->recPostorder(root);
    }
    //recursive insert
    Node *recInsert(Node *root, const vector<int> &point, int depth){
        if (root == nullptr){
            count++;
            return new Node(point);
        }
        int axis = depth % k;
        if (point[axis] < root->data[axis]){
            root->left = recInsert(root->left, point, depth + 1);
        } else {
            root->right = recInsert(root->right, point, depth + 1);
        }
        return root;
    }

    void insert(const vector<int> &point)
    {
        root = recInsert(root, point, 0);
    }
    //remove helper
    
    void copyData(vector<int> &a, const vector<int> &b)
    {
        for (int i = 0; i < a.size(); i++){
            a[i] = b[i];
        }
    }

    bool isEqual(const vector<int> &a, const vector<int> &b)
    {
        for (int i = 0; i < a.size(); i++){
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    //recursive remove
    Node *recRemove(Node *root, const vector<int> &point, int depth){
        if (root == nullptr) return nullptr;
        int axis = depth % k;
        if (isEqual(root->data, point)){
            if (root->right == nullptr){
                Node *temp = root->left;
                delete root;
                count--;
                return temp;
            }
            if (root->left == nullptr){
                Node *temp = root->right;
                delete root;
                count--;
                return temp;
            }
            Node *temp = root->right;
            while (temp->left != nullptr){
                temp = temp->left;
            }
            copyData(root->data, temp->data);
            root->right = recRemove(root->right, temp->data, depth + 1);
        } else if (point[axis] < root->data[axis]){
            root->left = recRemove(root->left, point, depth + 1);
        } else {
            root->right = recRemove(root->right, point, depth + 1);
        }
        return root;
    }



    void remove(const vector<int> &point){
        root = recRemove(root, point, 0);
    }

    bool recSrc(Node *node, const vector<int> &point, int depth){
        if (node == nullptr) return false;
        if (isEqual(node->data, point)) return true;
        int axis = depth % k;
        if (point[axis] < node->data[axis]){
            return recSrc(node->left, point, depth + 1);
        } else {
            return recSrc(node->right, point, depth + 1);
        }
    }
    bool search(const vector<int> &point){
        return recSrc(root, point, 0);
    }
    void buildTree(const vector<vector<int>> &pointList){
        for (int i = 0; i < pointList.size(); i++){
            insert(pointList[i]);
        }
    }
    // void nearestNeighbour(const vector<int> &target, Node *best);
    // void kNearestNeighbour(const vector<int> &target, int k, vector<Node *> &bestList);

};



// Please add more or modify as needed
