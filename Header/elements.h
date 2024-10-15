#ifndef ELEMENTS_H
#define ELEMENTS_H
#include <string>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "beam_TMSK2.h"
#include "beam_TMSK3.h"
#include "PQL_.h"
#include "PQS_.h"

class elements{
private:
    std::vector<int> size_of_eletype;
    std::vector<std::pair<int, int>> size_of_nodetype;
    std::vector<B2TS> B2TS_eles;
    std::vector<B3TS> B3TS_eles;
    std::vector<PQL_> PQL__eles;
    std::vector<PQS_> PQS__eles;

public:
    elements(std::vector<std::pair<std::string, int>> eles_num){
        size_of_eletype.resize(4, 0);
        for(std::vector<std::pair<std::string, int>>::iterator iten = eles_num.begin(); iten != eles_num.end(); iten++){
            size_of_eletype.at(match_type(iten->first)) = iten->second;
        }
        alloc();
    };

    elements(const std::vector<int> &sot = {0,0,0,0}){
        int size = sot.size();

        if(size != 4){
            printf("error size %i of input vector!(required 4)", size);
        }
        else{
            size_of_eletype = sot;
        }
        alloc();
    }

    void resize(const std::vector<int> &sot);

    void set_sont(const std::vector<std::pair<int,int>> &sont);

    /// @brief allocated memory according to size_of_eletype vector.
    void alloc();

    /// @brief find the name {element type} 's position
    /// @param name means element type
    /// @return position
    int match_type(std::string &name);

    /// @brief return a pointer of type's element set
    /// @param std_ele an element type's empty object
    /// @return element set reference
    std::vector<B2TS>* eleset_ptr(B2TS &std_ele);
    std::vector<B3TS>* eleset_ptr(B3TS &std_ele);
    std::vector<PQL_>* eleset_ptr(PQL_ &std_ele);
    std::vector<PQS_>* eleset_ptr(PQS_ &std_ele);

    /// @brief return an access to size_of_eletype vector(read only)
    /// @return const reference of size_of_eletype vector
    const std::vector<int>& sot();

    /// @brief return an access to size_of_nodetype vector(read only)
    /// @return const reference of size_of_nodetype vector
    const std::vector<std::pair<int, int>>& sont();

    ~elements(){};
};

#endif