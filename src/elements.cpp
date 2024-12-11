#include "elements.h"

void elements::resize(const std::vector<int> &sot){
    int size = sot.size(), i;
    bool changed = false;

    if(size > 5){
        printf("Warning: input vetor size %i is greater than expected!(only req-\nuired 5)\n", size);
    }
    
    for(i = 0; i < size; i++){
        if(size_of_eletype.at(i) != sot.at(i)){
            size_of_eletype.at(i) = sot.at(i);
            changed = true;
        }
    }
    if(changed) alloc();
};

void elements::set_sont(const std::vector<std::pair<int,int>> &sont){
    size_of_nodetype = sont;
};

void elements::alloc(){
    B2TS_eles.resize(size_of_eletype[0]);
    B2TS_eles.shrink_to_fit();

    B3TS_eles.resize(size_of_eletype[1]);
    B3TS_eles.shrink_to_fit();

    PQL__eles.resize(size_of_eletype[2]);
    PQL__eles.shrink_to_fit();
    
    PQS__eles.resize(size_of_eletype[3]);
    PQL__eles.shrink_to_fit();

    PSL__eles.resize(size_of_eletype[4]);
    PSL__eles.shrink_to_fit();
};

int elements::match_type(std::string &name){
    if(!_strcmpi(name.data(), "B2TS")){
        return 0;
    }
    else if(!_strcmpi(name.data(), "B3TS")){
        return 1;
    }
    else if(!_strcmpi(name.data(), "PQL_")){
        return 2;
    }
    else if(!_strcmpi(name.data(), "PQS_")){
        return 3;
    }
    else if(!_strcmpi(name.data(), "PSL_")){
        return 4;
    }
    return -1;
};

std::vector<B2TS>* elements::eleset_ptr(B2TS &std_ele){
    return &B2TS_eles;
};

std::vector<B3TS>* elements::eleset_ptr(B3TS &std_ele){
    return &B3TS_eles;
};

std::vector<PQL_>* elements::eleset_ptr(PQL_ &std_ele){
    return &PQL__eles;
};

std::vector<PQS_>* elements::eleset_ptr(PQS_ &std_ele){
    return &PQS__eles;
};

std::vector<PSL_>* elements::eleset_ptr(PSL_ &std_ele){
    return &PSL__eles;
};

const std::vector<int>& elements::sot(){
    return size_of_eletype;
}

const std::vector<std::pair<int, int>>& elements::sont(){
    return size_of_nodetype;
};