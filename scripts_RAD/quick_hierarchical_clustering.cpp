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
#include <mutex>
#include <functional>
#include <utility>

uint maxdepth(0);
using namespace std;
unordered_set<uint> visited; // AVOIDS TO PUT IT IN THE RECURSION STACK.
unordered_map <uint, unordered_set<uint> > nodeToNeighbors; // AVOIDS TO PUT IT IN THE RECURSION STACK.
void DFS(uint n, unordered_set<uint>& nodesInConnexComp){
    if (not visited.count(n)){
        visited.insert(n);
        nodesInConnexComp.insert(n);
            for (auto&& neigh : nodeToNeighbors[n]){
                DFS(neigh, nodesInConnexComp);
            }
        
    }
}

vector<string> split(const string &s, char delim){
    stringstream ss(s);
    string item;
    vector<string> elems;
    while (getline(ss, item, delim)) {
        elems.push_back(move(item)); 
    }
    return elems;
}

// awaits  input lines like this 14:12 85834 13 85835 14 12319 15 19508 154886 14536 19509 154887 80595 80603 
void parsingSRC(ifstream & refFile){
    string listNodes;
    // header
    vector<string> splitted1, splitted2, splitted3;
    uint read, source;
    while (not refFile.eof()){
        getline(refFile, listNodes);
        if (listNodes[0]=='#') continue; // HEADER.
        splitted1 = split(listNodes, ':');
        if (splitted1.size() > 1){
            splitted2 = split(splitted1[1], ' ');
            unordered_set<uint> reads;
            source = stoi(splitted1[0]);  // source read
            if (not splitted2.empty()){
                for (uint i(0); i < splitted2.size(); ++i){
                    read = stoi(splitted2[i]);  // recruited read
                    if (read != source){
                        reads.insert(read);
                        if (nodeToNeighbors.count(read)){
                            nodeToNeighbors[read].insert(source);
                        } else {
                            nodeToNeighbors.insert({read, {source}});
                        }
                    }
                }
            }
            if (nodeToNeighbors.count(source)){
                for (auto&& n : reads){
                    nodeToNeighbors[source].insert(n);
                }
            } else {
                nodeToNeighbors.insert({source, reads});
            }
        }   
    }
}



int main(int argc, char** argv){

    if (argc > 1){
        bool approx(false);
        //~ string outFileName("final_g_clusters.txt");
		
        string fileName(argv[1]);
        ifstream refFile(fileName);
        // parse SRC's output, for each line we get a node and its neighbors and fill the map nodeToNeighbors
        cerr << "Parsing..." << endl;
        parsingSRC(refFile);
        cerr << "Compute CCs..." << endl;
        uint nbConnexComp(0);
        vector<unordered_set<uint> > nodesInConnexComp;
        for (auto node(nodeToNeighbors.begin()); node != nodeToNeighbors.end(); ++node){
            if (not (visited.count(node->first))){
                unordered_set<uint> s;
                nodesInConnexComp.push_back(s);
                DFS(node->first, nodesInConnexComp.back());
                ++ nbConnexComp;
            }
        }
        cerr<<"Print CCs"<<endl;
        for (uint i(0); i < nodesInConnexComp.size(); ++i){
            for (auto&& n : nodesInConnexComp[i]){
                cout << n << " ";
            }
            cout << endl;
        }
        cerr << "Connected components: " << nbConnexComp << endl;
    }
}
