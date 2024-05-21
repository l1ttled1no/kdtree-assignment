#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */

const int max_int = 2147483647;

struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    int label;
    bool isChosen = false; //for k nearest neighbour
    kDTreeNode(vector<int> data,
               int label = -1, 
               kDTreeNode *left  = nullptr, 
               kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
        this->label = label;
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
    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for(int i = 0; i < node.data.size(); i++)
        {
            os << node.data[i];
            if (i != node.data.size() - 1) os << ", ";
        }
        os << ")";
        return os;
    }
    void reset(){ isChosen = false; }
    void set(){ isChosen = true; }
    
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;
    int count;
private:
    void kValue(int k){
        this->k = k;
    }
public:
    kDTree(int k = 2) 
    {
        this->k = k; 
        this->root = nullptr;
        this->count = 0;
    }
    //recursive destructor
    void recDestructor(kDTreeNode *root){
        if (root == nullptr) return;
        recDestructor(root->left);
        recDestructor(root->right);
        delete root;
    }
    ~kDTree()
    {
        recDestructor(root);
    }
    
    kDTreeNode *copy(const kDTreeNode *root)
    {
        if (root == nullptr) return nullptr;
        kDTreeNode *node = new kDTreeNode(root->data);
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
    //recursive count
    void recCount(kDTreeNode *root){
        if (root == nullptr) return;
        count++;
        recCount(root->left);
        recCount(root->right);
    }

    void updateCount(){
        count = 0;
        recCount(root);
    }
    int nodeCount() const{
        return count;
    }
    //recursive height
    int recHeight(kDTreeNode *root) const{
        if (root == nullptr) return 0;
        return 1 + max(recHeight(root->left), recHeight(root->right));
    }
    int height() const{
        return recHeight(root);
    }
    //recursive leafCount
    int recLeafCount(kDTreeNode *root) const{
        if (root == nullptr) return 0;
        if (root->left == nullptr && root->right == nullptr) return 1;
        return recLeafCount(root->left) + recLeafCount(root->right);
    }
    int leafCount() const{
        return recLeafCount(root);
    }

    void recInorder(kDTreeNode *root) const{
        if (root == nullptr) return;
        recInorder(root->left);
        OUTPUT << *root << " ";
        recInorder(root->right);
    }
    void inorderTraversal() const{
        this->recInorder(root);
    }
    void recPreorder(kDTreeNode *root) const{
        if (root == nullptr) return;
        OUTPUT << *root << " ";
        recPreorder(root->left);
        recPreorder(root->right);
    }
    void preorderTraversal() const{
        this->recPreorder(root);
    }
    void recPostorder(kDTreeNode *root) const{
        if (root == nullptr) return;
        recPostorder(root->left);
        recPostorder(root->right);
        OUTPUT << *root << " ";
    }
    void postorderTraversal() const{
        this->recPostorder(root);
    }
    //recursive insert
    kDTreeNode *recInsert(kDTreeNode *root, const vector<int> &point, int depth){
        if (root == nullptr){
            count++;
            return new kDTreeNode(point);
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

    void insertLabel(const vector<int> &point, int &label){
        if (search(point)){
            kDTreeNode *node = nodeSearch(point);
            node->label = label;
        }
        else OUTPUT << "Insert failed: Point not found\n";
    }
    //remove helper
    void copyData(kDTreeNode *a, kDTreeNode *b){
        for (int i = 0; i < a->data.size(); i++){
            a->data[i] = b->data[i];
        }
    }
    bool isEqual(const vector<int> &a, const vector<int> &b)
    {
        for (int i = 0; i < a.size(); i++){
            if (a[i] != b[i]) return false;
        }
        return true;
    }
    kDTreeNode* minNode(kDTreeNode* x, kDTreeNode* y, kDTreeNode* z, int d) {
        kDTreeNode* res = x;

        if (y != nullptr && y->data[d] < res->data[d])
            res = y;
        if (z != nullptr && z->data[d] < res->data[d])
            res = z;

        return res;
    }
    kDTreeNode* findMinRec(kDTreeNode* node, int d, int depth) {
        if (node == nullptr)
            return nullptr;
    
        int cd = depth % k;
    
        if (cd == d) {
            if (node->left == nullptr)
                return node;
            return findMinRec(node->left, d, depth + 1);
        }
    
        return minNode(node,
                       findMinRec(node->left, d, depth + 1),
                       findMinRec(node->right, d, depth + 1), 
                       d);
    }
    kDTreeNode* recRemove(kDTreeNode* root, const vector<int>& point, int depth) {
        if (root == nullptr) return nullptr;
        if (!search(point)) return root;
        int cd = depth % k;

        if (isEqual(root->data, point)) {
            if (root->right != nullptr) {
                kDTreeNode* min = findMinRec(root->right, cd, depth + 1);
                copyData(root, min);
                root->right = recRemove(root->right, min->data, depth + 1);
            } else if (root->left != nullptr) {
                kDTreeNode* min = findMinRec(root->left, cd, depth + 1);
                copyData(root, min);
                root->right = root->left;
                root->left = nullptr;
                root->right = recRemove(root->right, min->data, depth + 1);
            } else {
                delete root;
                return nullptr;
            }
        }

        else if (point[cd] < root->data[cd]) root->left = recRemove(root->left, point, depth + 1);
        else root->right = recRemove(root->right, point, depth + 1);

        return root;
    }
    void remove(const vector<int> &point){
        root = recRemove(root, point, 0);
        updateCount();
    }
    //search helper
    bool recSrc(kDTreeNode *node, const vector<int> &point, int depth){
        if (node == nullptr) return false;
        if (isEqual(node->data, point)) return true;
        int axis = depth % k;
        if (point[axis] < node->data[axis]){
            return recSrc(node->left, point, depth + 1);
        } else {
            return recSrc(node->right, point, depth + 1);
        }
    }
    kDTreeNode *nodeSrcRec(kDTreeNode *node, const vector<int> &point, int depth){
        if (node == nullptr) return nullptr;
        if (isEqual(node->data, point)) return node;
        int axis = depth % k;
        if (point[axis] < node->data[axis]){
            return nodeSrcRec(node->left, point, depth + 1);
        } else {
            return nodeSrcRec(node->right, point, depth + 1);
        }
    }

    kDTreeNode *nodeSearch(const vector<int> &point){
        return nodeSrcRec(root, point, 0);  
    }

    bool search(const vector<int> &point){
        return recSrc(root, point, 0);
    }
    //build tree helper
    /// use merge sort to sort the pointList
    void merge(vector<vector<int>>& pointList, int left, int mid, int right, int cd) {
        int i, j, k;
        int n1 = mid - left + 1;
        int n2 = right - mid;
    
        vector<vector<int>> L(n1), R(n2);
    
        for (i = 0; i < n1; i++)
            L[i] = pointList[left + i];
        for (j = 0; j < n2; j++)
            R[j] = pointList[mid + 1 + j];
    
        i = 0;
        j = 0;
        k = left;
        while (i < n1 && j < n2) {
            if (L[i][cd] <= R[j][cd]) {
                pointList[k] = L[i];
                i++;
            } else {
                pointList[k] = R[j];
                j++;
            }
            k++;
        }
    
        while (i < n1) {
            pointList[k] = L[i];
            i++;
            k++;
        }
    
        while (j < n2) {
            pointList[k] = R[j];
            j++;
            k++;
        }
    }
    
    void mergeSort(vector<vector<int>>& pointList, int left, int right, int cd) {
        if (left < right) {
            int mid = left + (right - left) / 2;
    
            mergeSort(pointList, left, mid, cd);
            mergeSort(pointList, mid + 1, right, cd);
            if (pointList[mid][cd] > pointList[mid + 1][cd]) 
                merge(pointList, left, mid, right, cd);
        }
        else return;
    }

    int medianPointList(vector<vector<int>>& pointList, int cd) {
        int median = (pointList.size() - 1) / 2;
        while (median > 0 && pointList[median][cd] == pointList[median - 1][cd]) {
            median--;
        }
        return median;
    }
    kDTreeNode* buildTreeRec(vector<vector<int>>& pointList, int depth) {
        if (pointList.empty())
            return nullptr;

        int cd = depth % this->k; // Current dimension

        // Use mergeSort instead of sort
        mergeSort(pointList, 0, pointList.size() - 1, cd); // Sort pointList by current dimension

        // Find median
        int med = medianPointList(pointList, cd);
        // int median = pointList.size() / 2;  
        // Create new node and construct subtrees
        kDTreeNode* node = new kDTreeNode(pointList[med]);
        vector<vector<int>> leftPoints(pointList.begin(), pointList.begin() + med);
        vector<vector<int>> rightPoints(pointList.begin() + med + 1, pointList.end());
        node->left = buildTreeRec(leftPoints, depth + 1);
        node->right = buildTreeRec(rightPoints, depth + 1);
        return node;
    }

    void buildTree(vector<vector<int>>& pointList) {
        this->root = buildTreeRec(pointList, 0);
        this->count = pointList.size();
    }
    //supreme buildtree: build tree with label, 
    //nearest neighbour helper

    double distance(const vector<int> &a, const vector<int> &b){
        double sum = 0;
        for (int i = 0; i < a.size(); i++){
            sum += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return sqrt(sum);
    }
    void nearestNeigbourRec (kDTreeNode *&root, const vector<int> &target, kDTreeNode *&best, int depth, double &bestDist){
        if (root == nullptr) return;
        //Get to leaf node
        if (root->left == nullptr && root->right == nullptr){
            double dist = distance(root->data, target);
            if (dist < bestDist && !root->isChosen){
                bestDist = dist;
                best = root;
            }
            return;
        }
        int axis = depth % k;
        if (target[axis] < root->data[axis]){ // if target < root
            nearestNeigbourRec(root->left, target, best, depth + 1, bestDist); //go left
            double dist = distance(root->data, target); //check distance
            if (dist < bestDist && !root->isChosen){ //update best
                bestDist = dist;
                best = root;
            }
            double distAxis = abs(target[axis] - root->data[axis]); //distance to axis
            if (bestDist > distAxis){ //check if need to go right
                nearestNeigbourRec(root->right, target, best, depth + 1, bestDist);
            }

        } else {
            nearestNeigbourRec(root->right, target, best, depth + 1, bestDist); //go right
            double dist = distance(root->data, target); //check distance
            if (dist < bestDist && !root->isChosen){ //update best
                bestDist = dist;
                best = root;
            }
            double distAxis = abs(target[axis] - root->data[axis]); //distance to axis
            if (bestDist > distAxis){ //check if need to go left
                nearestNeigbourRec(root->left, target, best, depth + 1, bestDist);
            }
        }        
    }
    
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best){
        best = nullptr; //reset best
        // nearestNeighbourRec_(root, target, best, 0)->print();
        double maxDist = max_int;
        nearestNeigbourRec(root, target, best, 0, maxDist);
        // if (best != nullptr) OUTPUT << *best << " ";

    }

    //k nearest neighbour helper
    
    void resetChosen(kDTreeNode *root){
        if (root == nullptr) return;
        root->reset();
        resetChosen(root->left);
        resetChosen(root->right);
    }
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList){
        resetChosen(root);
        for (int i = 0; i < k; i++){
            kDTreeNode *best = nullptr;
            double maxDist = max_int;
            nearestNeigbourRec(root, target, best, 0, maxDist);
            if (best != nullptr){
                best->set();
                bestList.push_back(best);
            }
            else break;
        }
        resetChosen(root);
    }
    friend class kNN;
};

class kNN
{
private:
    int k;
    Dataset *X_train;
    Dataset *y_train;
    kDTree tree; 

public:
    kNN(int k = 5){ this->k = k; }
    vector<vector<int>> convertList(const list<list<int>> &l){
        vector<vector<int>> res;
        for (auto i : l){
            vector<int> temp;
            for (auto j : i){
                temp.push_back(j);
            }
            res.push_back(temp);
        }
        return res;
    }
    //destructor

    list<int> convertVector(const vector<int> &v){
        list<int> res;
        for (auto i : v){
            res.push_back(i);
        }
        return res;
    }
    void fit(Dataset &X_train, Dataset &y_train){
        this->X_train = &X_train;
        this->y_train = &y_train;
        tree.k = X_train.columnName.size();
        vector<vector<int>> pointList = convertList(this->X_train->data);
        tree.buildTree(pointList);
        
        vector<vector<int>> labelList = convertList(this->y_train->data);
        for (int i = 0; i < labelList.size(); i++){
            tree.insertLabel(pointList[i], labelList[i][0]);
        }
    }

    Dataset predict(Dataset &X_test)
    {
        Dataset res; 
        res.columnName.push_back("label");
        vector<vector<int>> pointList = convertList(X_test.data);

        for (int i = 0; i < pointList.size(); i++){
            vector <kDTreeNode *> bestList;
            tree.kNearestNeighbour(pointList[i], k, bestList);
            // vector 10 label count at 0
            int labelCount[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (int j = 0; j < bestList.size(); ++j){
                labelCount[bestList[j]->label]++;
            }
            //find max label
            int max = -712, label = 0;
            for (int j = 0; j < 10; j++){
                if (labelCount[j] > max){
                    max = labelCount[j];
                    label = j;
                }
            }

            list<int> temp; 
            temp.push_back(label);
            res.data.push_back(temp);
        }
        return res;
    }
    
    double score(const Dataset &y_test, const Dataset &y_pred){
        vector<vector<int>> y_test_list = convertList(y_test.data);
        vector<vector<int>> y_pred_list = convertList(y_pred.data); 
        int correct = 0;
        for (int i = 0; i < y_test_list.size(); i++){
            if (y_test_list[i][0] == y_pred_list[i][0]){
                correct++;
            }
        }
        return (double)correct / y_test.data.size();
    }
    void print_Y(const Dataset &y){
        vector<vector<int>> y_list = convertList(y.data);
        OUTPUT << "label: "; 
        for (int i = 0; i < y_list.size(); i++){
            OUTPUT << y_list[i][0] << " ";
        }
        // OUTPUT << endl;
    }
};


// Please add more or modify as needed
