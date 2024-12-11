#include <string.h>
#include <string>
#include <fstream>
#include <vector>
#include "read.h"

bool getlmsg(std::ifstream &file_stream, std::string &str){
    std::getline(file_stream, str);
    while (true){
        if (!(str.length() == 0 || (str.at(0) == '#'))){
            break;
        }
        else{
            if(file_stream.eof()){
                return true;
            }
            std::getline(file_stream, str);
        }
    }
    std::string str_f, str_b;
    split(str, str_f, str_b, '#');
    str = str_f;
    return false;
};

void first_wd(std::string &str, std::string &word){
    bool stop = false;
    word.clear();
    std::vector<std::string::iterator> del_list;
    for(std::string::iterator i = str.begin(); i != str.end(); i++){
        if(*i == ','){
            del_list.push_back(i);
            break;
        }
        else if(*i == ' '){
            del_list.push_back(i);
            // stop = true;
        }
        else{
            // if(stop) break;
            word.push_back(*i);
            del_list.push_back(i);
        }
    }
    for (std::vector<std::string::iterator>::reverse_iterator i = del_list.rbegin(); i != del_list.rend(); i++){
        str.erase(*i);
    }
    str.shrink_to_fit();
};

void int2str(const int* i_arr, const int len, std::string& str){
    str.clear();
    str.shrink_to_fit();
    char cdum;
    for(int i = 0; i < len; i++){
        if(i_arr[i] > 127 || i_arr[i] < 0){
            printf("\033[33mWARNING\033[0m: the value %i in Int array not in [0,127], ignored!\n", i_arr[i]);
        }
        else{
            cdum = i_arr[i];
        }
        str.append({cdum});
    }
};

void str2int(const std::string& str, const int len, int* i_arr){
    int i = 0;
    for(std::string::const_iterator ichar = str.begin(); ichar != str.end(); ichar++){
        i_arr[i] = *ichar;
        i++;
    }
    while(i < len){
        i_arr[i] = 32;
        i++;
    }
};

void trim(std::string& str){
    int start = str.find_first_not_of(' ');
    int end = str.find_last_not_of(' ');
    if (start == std::string::npos) {
        str = "";
    }
    else {
        str = str.substr(start, end - start + 1);
    }
};

bool split(const std::string& str, std::string& front, std::string& back, char charater){
    int posi = str.find_first_of(charater);
    int len = str.size();
    if(posi == std::string::npos){
        front = str;
        back = "";
        return false;
    }
    else{
        front = str.substr(0, posi);
        back = str.substr(posi + 1, len - posi);
        return true;
    }
};