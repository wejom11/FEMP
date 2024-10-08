#ifndef ELEMENTS_H
#define ELEMENTS_H
#include <string>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "beam_TMSK2.h"
#include "beam_TMSK2R.h"
#include "beam_TMSK3.h"
#include "beam_TMSK3R.h"

class elements{
private:
    std::vector<int> size_of_eletype;
    std::vector<std::pair<int, int>> size_of_nodetype;
    std::vector<B2TS> B2TS_eles;
    std::vector<B2TSR> B2TSR_eles;
    std::vector<B3TS> B3TS_eles;
    std::vector<B3TSR> B3TSR_eles;

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

    /// @brief return a reference of type's element set
    /// @param type element type
    /// @return element set reference
    std::vector<B2TS>* eleset_ptr(B2TS &std_ele);
    std::vector<B2TSR>* eleset_ptr(B2TSR &std_ele);
    std::vector<B3TS>* eleset_ptr(B3TS &std_ele);
    std::vector<B3TSR>* eleset_ptr(B3TSR &std_ele);

    /// @brief return an access to size_of_eletype vector(read only)
    /// @return const reference of size_of_eletype vector
    const std::vector<int>& sot();

    /// @brief return an access to size_of_nodetype vector(read only)
    /// @return const reference of size_of_nodetype vector
    const std::vector<std::pair<int, int>>& sont();

    ~elements(){};
};

#endif