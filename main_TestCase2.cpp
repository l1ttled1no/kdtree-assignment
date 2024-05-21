#include "kDTree.hpp"



/*
 * Author: Võ Tiến
 * Date: 30.04.2023
 * FB: https://www.facebook.com/profile.php?id=100056605580171
 * FB nhóm: https://www.facebook.com/groups/211867931379013
 */
/*
* cách dùng hàm main.cpp
* b2: mở Terminal trên thanh công cụ mở new Terminal
* b3: chạy lệnh g++ -o main main_TestCase2.cpp ID_BTL2_DSA.cpp nếu không gì hiện ra thì không bị lỗi code
* b4: chạy ./main KEYWORDS (./main task1, ./main task2)
* b5: nếu hiện Success: test 1 to 200 thì đã đúng hết

! KEYWORDS có thể là
* ./main 2 này chạy phần input2
* ./main 1 5 này sẽ chạy phần input1 -> input5
* output là file của bạn, expect là đáp án
*/


#include <algorithm>

#define folder_input "TestCase_2/input/input"
#define folder_output "TestCase_2/output/output"
#define folder_expect "TestCase_2/expect/expect"

bool isNumber(string str);
void comparefile(int start, int end);


void run_task_2(string fileName)
{
    ifstream file(folder_input + fileName + ".txt");
    if(!file.is_open())
    {
        cout << "Error: " << fileName << endl;
        return;
    }

    int toRow = 0, fromRow = -1;
    int toColFeature = 1, fromColFeature = -1;
    int K_knn = 5;
    double test_size = 0.2;

    string s;
    file >> s; OUTPUT << "input + " << s << endl;
    file >> s; OUTPUT << s << endl; toRow = stoi(s.substr(6));
    file >> s; OUTPUT << s << endl; fromRow = stoi(s.substr(8));
    file >> s; OUTPUT << s << endl; toColFeature = stoi(s.substr(13));
    file >> s; OUTPUT << s << endl; fromColFeature = stoi(s.substr(15));
    file >> s; OUTPUT << s << endl; K_knn = stoi(s.substr(6));
    file >> s; OUTPUT << s << endl; test_size = stod(s.substr(10));
    OUTPUT << "------------------------------------" << endl;

    Dataset dataset;
    dataset.loadFromCSV("TestCase_2/mnist.csv");
    int nRows, nCols;

    kNN knn(K_knn);
    Dataset X_train, X_test, y_train, y_test;
    Dataset feature = dataset.extract(toRow, fromRow, toColFeature , fromColFeature);
    Dataset label = dataset.extract(toRow, fromRow, 0, 0);

    train_test_split(feature, label, test_size, X_train, X_test, y_train, y_test);
    knn.fit(X_train, y_train);
    Dataset y_pred = knn.predict(X_test);
    double accuracy = knn.score(y_test, y_pred);
    OUTPUT << "y_pred" << endl;
    knn.print_Y(y_pred);
    OUTPUT << "y_test" << endl;
    knn.print_Y(y_test);
    OUTPUT << "Accuracy: " << accuracy << endl;   
}



int main(int argc, char *argv[])
{
    cout << "Start program assignments 1 DSA" << endl;

    int START = 0, END = 0;
    int fullTask = 2000;
    //* check all test cases

    if (argc == 1)
    {
        cout << "checking test : ";
        START = 1;
        END = 2000;
        for (int i = START; i <= END; i++)
        {
            cout << i << " ";
            OUTPUT.open(folder_output + to_string(i) + ".txt");
            run_task_2(to_string(i));
          
            OUTPUT.close();
        }
    }

    else if (argc == 2)
    {
    
        string stringTask = argv[1];
        if (isNumber(stringTask))
        {
            cout << "checking test : " << argv[1];
            START = stoi(stringTask);
            END = stoi(stringTask);
            OUTPUT.open(folder_output + stringTask + ".txt");
            run_task_2(stringTask);
            OUTPUT.close();
        }
        else if (stringTask == "me")
        {
            cout << "checking test : ";
            for (int i = 1; i <= 2000; i++)
            {
                cout << i << " ";
                OUTPUT.open(folder_expect + to_string(i) + ".txt");
                run_task_2(to_string(i));
                    OUTPUT.close();
            }
        }
        else if (stringTask.substr(0, 4) != "task" || stringTask[4] <= '0' || stringTask[4] > '3')
        {
            cout << "Please re-enter the correct command" << endl;
            cout << "./main task1\n./main task2\n./main task3\n./main number" << endl;
            return 0;
        }
        else
        {
            cout << "checking test : ";
            for (int i = 1; i <= 2000; i++)
            {
                cout << i << " ";
                OUTPUT.open(folder_output + to_string(i) + ".txt");
                run_task_2(to_string(i));
                OUTPUT.close();
            }
        }
    }
    else if (argc == 3 && isNumber(argv[1]) && isNumber(argv[2]))
    {
        cout << "checking test : ";
        START = stoi(argv[1]);
        END = stoi(argv[2]);
        for (int i = stoi(argv[1]); i <= stoi(argv[2]); i++)
        {
            cout << i << " ";
            OUTPUT.open(folder_output + to_string(i) + ".txt");
            run_task_2(to_string(i));
            OUTPUT.close();
        }
    }
    else
    {
        cout << "Please re-enter the correct command" << endl;
        cout << "./main task1\n./main task2\n./main task3OUTPUT\n./main number" << endl;
        return 0;
    }

    cout << "\nOK: without errors\n"
         << endl;
    comparefile(START, END);
    return 1;
}

bool isNumber(string str)
{
    for (int i = 0; i < str.length(); i++)
    {
        if (!isdigit(str[i]))
            return false;
    }
    return true;
}

void comparefile(int start, int end)
{
    vector<int> result;
    for (int i = start; i <= end; i++)
    {
        ifstream read_output(folder_output + to_string(i) + ".txt");
        ifstream read_expect(folder_expect + to_string(i) + ".txt");
        if (read_output.fail() || read_expect.fail())
        {
            cout << "Error reading file" << end;
            return;
        }
        string s1, s2;
        while (read_output >> s1 && read_expect >> s2)
        {
            if (s1 != s2)
            {
                result.push_back(i);
                break;
            }
        }
        if (read_output >> s1 || read_expect >> s2)
        {
            if (result.size() == 0 || result.size() > 0 && result[result.size() - 1] != i)
                result.push_back(i);
        }
    }

    if (result.size() == 0)
    {
        cout << "Success: test " << start << " to " << end << endl;
    }
    else
    {
        cout << "percent success : " << (1 - result.size() * 1.0 / (end - start + 1)) * 100 << "%" << endl;
        cout << "Fail : test [";
        for (int i = 0; i < result.size() - 1; i++)
        {
            cout << result[i] << ", ";
        }
        cout << result[result.size() - 1] << "]\n";
        cout << "link check comparefile: https://www.diffchecker.com/text-compare/" << endl;
    }
}