#include "elements.h"

void elements::resize(const std::vector<int> &sot){
    int size = sot.size();

    if(size != 4){
        printf("error size %i of input vector!(required 4)", size);
    }
    else{
        size_of_eletype = sot;
    }
    alloc();
};

void elements::set_sont(const std::vector<std::pair<int,int>> &sont){
    size_of_nodetype = sont;
};

void elements::alloc(){
    B2TS_eles.resize(size_of_eletype[0]);
    B2TS_eles.shrink_to_fit();

    B2TSR_eles.resize(size_of_eletype[1]);
    B2TSR_eles.shrink_to_fit();

    B3TS_eles.resize(size_of_eletype[2]);
    B3TS_eles.shrink_to_fit();

    B3TSR_eles.resize(size_of_eletype[3]);
    B3TSR_eles.shrink_to_fit();
};

int elements::match_type(std::string &name){
    if(!_strcmpi(name.data(), "B2TS")){
        return 0;
    }
    else if(!_strcmpi(name.data(), "B2TSR")){
        return 1;
    }
    else if(!_strcmpi(name.data(), "B3TS")){
        return 2;
    }
    else if(!_strcmpi(name.data(), "B3TSR")){
        return 3;
    }
    return -1;
};

std::vector<B2TS>* elements::eleset_ptr(B2TS &std_ele){
    return &B2TS_eles;
};

std::vector<B2TSR>* elements::eleset_ptr(B2TSR &std_ele){
    return &B2TSR_eles;
};

std::vector<B3TS>* elements::eleset_ptr(B3TS &std_ele){
    return &B3TS_eles;
};

std::vector<B3TSR>* elements::eleset_ptr(B3TSR &std_ele){
    return &B3TSR_eles;
};

const std::vector<int>& elements::sot(){
    return size_of_eletype;
}

const std::vector<std::pair<int, int>>& elements::sont(){
    return size_of_nodetype;
};