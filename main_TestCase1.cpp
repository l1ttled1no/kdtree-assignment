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
* b3: chạy lệnh g++ -o main main_TestCase1.cpp ID_BTL2_DSA.cpp Dataset.o nếu không gì hiện ra thì không bị lỗi code
* b4: chạy ./main KEYWORDS (./main task1, ./main task2)
* b5: nếu hiện Success: test 1 to 200 thì đã đúng hết

! KEYWORDS có thể là
* ./main 2 này chạy phần input2
* ./main 1 5 này sẽ chạy phần input1 -> input5
* output là file của bạn, expect là đáp án
*/

#include <algorithm>

#define folder_input "TestCase_1/input/input"
#define folder_output "TestCase_1/output/output"
#define folder_expect "TestCase_1/expect/expect"

bool isNumber(string str);
void comparefile(int start, int end);

// Hàm so sánh để sắp xếp các Node dựa trên data
bool compareNodes(const kDTreeNode *a, const kDTreeNode *b)
{
    std::vector<int> vec1 = a->data;
    std::vector<int> vec2 = b->data;
    return std::lexicographical_compare(vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
}
void run_task_1(string fileName)
{
    ifstream file(folder_input + fileName + ".txt");
    if (!file.is_open())
    {
        cout << "Error: " << fileName << endl;
        return;
    }
    string s;
    file >> s;
    OUTPUT << "line " << 1 << ": " << s << endl;
    int k = stoi(s.substr(2));
    kDTree *VOTIEN = new kDTree(k);

    int value, index, i = 2, number;
    while (file >> s)
    {
        if (s == "nodeCount_height_leafCount")
        {
            OUTPUT << "line " << i << ": " << s << endl;
            OUTPUT << "nodeCount=" << VOTIEN->nodeCount() << ";" << "height=" << VOTIEN->height() << ";" << "leafCount=" << VOTIEN->leafCount() << "\n";
        }
        else if (s == "inorderTraversal")
        {
            OUTPUT << "line " << i << ": " << s << endl;
            VOTIEN->inorderTraversal();
            OUTPUT << endl;
        }
        else if (s == "preorderTraversal")
        {
            OUTPUT << "line " << i << ": " << s << endl;
            VOTIEN->preorderTraversal();
            OUTPUT << endl;
        }
        else if (s == "postorderTraversal")
        {
            OUTPUT << "line " << i << ": " << s << endl;
            VOTIEN->postorderTraversal();
            OUTPUT << endl;
        }
        else if (s == "insert")
        {
            vector<int> data;
            OUTPUT << "line " << i << ": " << s << " ";
            file >> s;
            OUTPUT << s << endl;
            stringstream ss(s);

            ss.ignore(); // Bỏ qua ký tự '(' đầu tiên
            int number;
            while (ss >> number)
            {
                data.push_back(number);
                if (ss.peek() == ',')
                    ss.ignore(); // Bỏ qua ký tự ',' sau mỗi số
            }

            VOTIEN->insert(data);
        }
        else if (s == "remove")
        {
            vector<int> data;
            OUTPUT << "line " << i << ": " << s << " ";
            file >> s;
            OUTPUT << s << endl;
            stringstream ss(s);

            ss.ignore(); // Bỏ qua ký tự '(' đầu tiên
            int number;
            while (ss >> number)
            {
                data.push_back(number);
                if (ss.peek() == ',')
                    ss.ignore(); // Bỏ qua ký tự ',' sau mỗi số
            }
            VOTIEN->remove(data);
        }
        else if (s == "search")
        {
            vector<int> data;
            OUTPUT << "line " << i << ": " << s << " ";
            file >> s;
            OUTPUT << s << endl;
            stringstream ss(s);

            ss.ignore(); // Bỏ qua ký tự '(' đầu tiên
            int number;
            while (ss >> number)
            {
                data.push_back(number);
                if (ss.peek() == ',')
                    ss.ignore(); // Bỏ qua ký tự ',' sau mỗi số
            }

            if (VOTIEN->search(data))
            {
                OUTPUT << "TRUE" << endl;
            }
            else
                OUTPUT << "FALSE" << endl;
        }
        else if (s == "buildTree")
        {
            vector<vector<int>> pointList;
            int number = 0;
            file >> number;
            OUTPUT << "line " << i << ": " << s << " " << number << endl;

            for (int i = 0; i < number; i++)
            {
                file >> s;
                OUTPUT << s << " ";
                stringstream ss(s);
                vector<int> data;
                ss.ignore(); // Bỏ qua ký tự '(' đầu tiên
                int number;
                while (ss >> number)
                {
                    data.push_back(number);
                    if (ss.peek() == ',')
                        ss.ignore(); // Bỏ qua ký tự ',' sau mỗi số
                }
                pointList.push_back(data);
            }
            OUTPUT << endl;
            VOTIEN->buildTree(pointList);
        }
        else if (s == "clear")
        {
            OUTPUT << "line " << i << ": " << s << endl;
            delete VOTIEN;
            VOTIEN = new kDTree(k);
        }
        else if (s == "nearestNeighbour")
        {
            vector<int> data;
            OUTPUT << "line " << i << ": " << s << " ";
            file >> s;
            OUTPUT << s << endl;
            stringstream ss(s);

            ss.ignore(); // Bỏ qua ký tự '(' đầu tiên
            int number;
            while (ss >> number)
            {
                data.push_back(number);
                if (ss.peek() == ',')
                    ss.ignore(); // Bỏ qua ký tự ',' sau mỗi số
            }
            kDTreeNode *best = nullptr;
            VOTIEN->nearestNeighbour(data, best);
            if (best)
                OUTPUT << *best;
            OUTPUT << endl;
        }
        else if (s == "kNearestNeighbour")
        {
            vector<int> data;
            int kx;
            OUTPUT << "line " << i << ": " << s << " ";
            file >> kx >> s;
            OUTPUT << kx << " " << s << endl;
            stringstream ss(s);

            ss.ignore(); // Bỏ qua ký tự '(' đầu tiên
            int number;
            while (ss >> number)
            {
                data.push_back(number);
                if (ss.peek() == ',')
                    ss.ignore(); // Bỏ qua ký tự ',' sau mỗi số
            }
            vector<kDTreeNode *> best;
            VOTIEN->kNearestNeighbour(data, kx, best);

            // for (int i = 0; i < best.size(); i++)
            // {
            //     for(int j = 0; j < best.size() - 1; j++)
            //     {
            //         if(compareNodes(best[j+1], best[j])) swap(best[j], best[j+1]);
            //     }
            // }
            for (kDTreeNode *i : best)
            {
                OUTPUT << *i;
            }
            OUTPUT << endl;
        }

        i++;
    }
    delete VOTIEN;
}

void run_task_2(string fileName) {}
void run_task_3(string fileName) {}

//* ./main task1; ./main task2; ./main task3;
//* ./main i j (với i với j là số từ test i -> j)
//* ./main i (i là test cần test)
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
            run_task_1(to_string(i));

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
            run_task_1(stringTask);
            OUTPUT.close();
        }
        else if (stringTask == "me")
        {
            cout << "checking test : ";
            for (int i = 1; i <= 2000; i++)
            {
                cout << i << " ";
                OUTPUT.open(folder_expect + to_string(i) + ".txt");
                run_task_1(to_string(i));
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
                run_task_1(to_string(i));
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
            run_task_1(to_string(i));
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