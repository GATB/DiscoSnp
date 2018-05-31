#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <map>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <stack>          // std::stack
//#include <mutex>
//#include <functional>
//#include <utility>
//typedef	unsigned long long	u_int64_t;

using namespace std;
unordered_set <string> visited; // AVOIDS TO PUT IT IN THE RECURSION STACK.
unordered_map <string, unordered_set<string> > nodeToNeighbors; // AVOIDS TO PUT IT IN THE RECURSION STACK.


// DEPRECATED. Recursion depth is to high and causes the computation fails. 
void DFS(const string n ){
    visited.insert(n);
    cout<<" "<<n;
    for (auto&& neigh : nodeToNeighbors[n]){
        if (not visited.count(neigh)) DFS(neigh);
    }
}


void nonrecursive_DFS(const string n){
    stack<string> mystack;
    mystack.push(n);
    visited.insert(n);
    while (not mystack.empty()){
        auto m = mystack.top();
        mystack.pop();
        cout<<" "<<m;
        for (auto&& neigh : nodeToNeighbors[m]){
            if (visited.count(neigh)==0)
                mystack.push(neigh);
                visited.insert(neigh);
        }
    }
}

// awaits  input lines like this 12 85834
void parsingPairsOfNodes(ifstream & refFile){
    string listNodes;
    while (not refFile.eof()){
        getline(refFile, listNodes);
        if (listNodes.size() == 0) break;
        if (listNodes[0]=='#') continue; // HEADER.
        
        string node1, node2;
        bool firstnode=true;
        for (char c : listNodes){
            if (c=='\t' || c==' ') {firstnode=false;continue;}
            if (firstnode)  node1+=c;
            else            node2+=c;
        }
        if (nodeToNeighbors.count(node1))  nodeToNeighbors[node1].insert(node2);
        else nodeToNeighbors.insert({node1, {node2}});
        if (nodeToNeighbors.count(node2))  nodeToNeighbors[node2].insert(node1);
        else nodeToNeighbors.insert({node2, {node1}});
    }
}




int main(int argc, char** argv){

    if (argc > 1){
		
        string fileName(argv[1]);
        ifstream refFile(fileName);
        // parse SRC's output, for each line we get a node and its neighbors and fill the map nodeToNeighbors
        cerr << "Parsing..." << endl;
        parsingPairsOfNodes(refFile);
        cerr << "Compute CCs..." << endl;
        for (auto node(nodeToNeighbors.begin()); node != nodeToNeighbors.end(); ++node){
            if (visited.count(node->first)==0){
//                DFS(node->first);
                nonrecursive_DFS(node->first);
                cout << endl;
            }
        }
    }
    return 0;
}
