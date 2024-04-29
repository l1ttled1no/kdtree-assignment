#include "main.hpp"

class Dataset
{
private:
    list<list<int>> data;
    vector<string> columnName;

public:
    Dataset();
    ~Dataset();
    bool loadFromCSV(const char *fileName);
    void printHead(int nRows = 5, int nCols = 5) const;
    void printTail(int nRows = 5, int nCols = 5) const;
    void getShape(int &nRows, int &nCols) const;
    void columns() const;
    bool drop(int axis = 0, int index = 0, std::string columns = "");
    Dataset extract(int startRow = 0, int endRow = -1, int startCol = 0, int endCol = -1) const;
    Dataset(const Dataset &other);
    Dataset &operator=(const Dataset &other);
    friend void train_test_split(Dataset &X, Dataset &y, double test_size,
                                 Dataset &X_train, Dataset &X_test, Dataset &y_train, Dataset &y_test);
    friend class kNN;
};

void train_test_split(Dataset &X, Dataset &y, double test_size,
                      Dataset &X_train, Dataset &X_test, Dataset &y_train, Dataset &y_test);
