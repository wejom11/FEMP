#include "assemble.h"

bool set_seq(int** paras, int ptr, const std::string& val){
    std::string wd_front, wd_back;

    if((ptr == 0) || (ptr >= 5 && ptr <= 7) || (ptr >= 9 && ptr <= 11) || (ptr == 14) || (ptr >= 17 && ptr <= 18)){
        bool is_num = true;
        for(std::string::const_iterator ich = val.begin(); ich != val.end(); ich++){
            if((*ich < 48 || *ich > 57) && (*ich != 45) && (*ich != 43)){
                is_num = false;
                break;
            }
        }
        if(!is_num) printf("\033[33mWARNING\033[0m: [\033[1;34m%s\033[0m] isn't a Integer, recommanded reinput!\n", val.data());
        
        try{
            paras[ptr] = new int{std::stoi(val)};
        }
        catch(std::invalid_argument& e){
            std::cerr << "Please input Integer!" << std::endl;
            return false;
        }
    }
    else if(ptr == 8 || (ptr >= 1 && ptr <= 4) || (ptr >= 19 && ptr <= 21)){
        bool is_num = true;
        for(std::string::const_iterator ich = val.begin(); ich != val.end(); ich++){
            if((*ich < 48 || *ich > 57) && (*ich != 45) && (*ich != 43) && (*ich != 'E') && (*ich != 'e')
                && (*ich != '.')){
                is_num = false;
                break;
            }
        }
        if(!is_num) printf("\033[33mWARNING\033[0m: [\033[1;34m%s\033[0m] isn't a number, recommanded reinput!\n", val.data());
        if(!*paras){
            printf("\033[33mWARNING\033[0m: exponent digits(Exp_Digit) should configure first\nUse default value -2");
            *paras = new int{-2};
        }
        
        bool is_ed = split(val, wd_front, wd_back, 'E');
        bool is_d = false;
        if(!is_ed) is_ed = split(val, wd_front, wd_back, 'e');
        if(!is_ed) is_d = split(val, wd_front, wd_back, '.');
        
        if(is_d){
            int size = wd_back.size();
            int sign = *paras[0] > 0 ? -1 : 1;
            std::string exper = std::to_string(size);
            if(exper.size() > abs(*paras[0])){
                printf("\033[1;31mERROR\033[0m: the exponent digits not enough for %s\n\
Please ste [Exp_Digit] larger\n", val.data());
            }
            else{
                for(int i = exper.size(); i < abs(*paras[0]); i++){
                    exper.insert(exper.begin(), '0');
                }
            }

            try{
                paras[ptr] = new int{sign * std::stoi(wd_front + wd_back + exper)};
            }
            catch(std::invalid_argument& e){
                std::cerr << "Please input value!" << std::endl;
                return false;
            }
        }
        else{
            std::string subwd_f, subwd_b;
            int sign1 = *paras[0] > 0 ? -1 : 1;
            if(is_ed){
                split(wd_front, subwd_f, subwd_b, '.');
                int size = subwd_b.size();
                int exp = 0;
                if(wd_back.size() > 0){
                    try{
                        exp = size - std::stoi(wd_back);
                    }
                    catch(std::invalid_argument& e){
                        std::cerr << "Please input value!" << std::endl;
                        return false;
                    }
                }
                else{
                    exp = size;
                }
                int sign2 = exp > 0 ? 1 : -1;

                std::string exper = std::to_string(abs(exp));
                if(exper.size() > abs(*paras[0])){
                    printf("\033[1;31mERROR\033[0m: the exponent digits not enough for %s\n\
Please ste [Exp_Digit] larger\n", val.data());
                }
                else{
                    for(int i = exper.size(); i < abs(*paras[0]); i++){
                        exper.insert(exper.begin(), '0');
                    }
                }

                try{
                    paras[ptr] = new int{sign1 * sign2 * std::stoi(subwd_f + subwd_b + exper)};
                }
                catch(std::invalid_argument& e){
                    std::cerr << "Please input value!" << std::endl;
                    return false;
                }
            }
            else{
                int sign = *paras[0] > 0 ? -1 : 1;
                std::string exper = "0";
                if(exper.size() > abs(*paras[0])){
                    printf("\033[1;31mERROR\033[0m: the exponent digits not enough for %s\n\
Please ste [Exp_Digit] larger\n", val.data());
                }
                else{
                    for(int i = exper.size(); i < abs(*paras[0]); i++){
                        exper.insert(exper.begin(), '0');
                    }
                }
                try{
                    paras[ptr] = new int{sign * std::stoi(val + exper)};
                }
                catch(std::invalid_argument& e){
                    std::cerr << "Please input value!" << std::endl;
                    return false;
                }
            }
        }
    }
    else if((ptr >= 12 && ptr <= 13) || (ptr >= 15 && ptr <= 16)){
        if(!paras[10]){
            printf("\033[33mWARNING\033[0m: string length(Str_Len) should configure first\nUse default value 64\n");
            paras[10] = new int{64};            
        }
        int* istr = new int[*paras[10]];
        for(int i = 0; i < *paras[10]; i++){
            istr[i] = 32;
        }

        str2int(val, *paras[10], istr);
        paras[ptr] = istr;
    }
    else{
        printf("\033[33mWARNING\033[0m: parameters[%i] don't need configure", ptr);
        return false;
    }
    return true;
};

void show_seq(int** paras, int ptr, std::string& val){
    std::string wd_front, wd_back;
    if(!paras[ptr]){
        val = "\033[0m\033[31mUnset   ";
        return;
    }

    if((ptr == 0) || (ptr >= 5 && ptr <= 7) || (ptr >= 9 && ptr <= 11) || (ptr == 14) || (ptr >= 17 && ptr <= 18)){
        val = std::to_string(*paras[ptr]);
    }
    else if(ptr == 8 || (ptr >= 1 && ptr <= 4) || (ptr >= 19 && ptr <= 21)){
        int exp_digit = pow(10, abs(**paras));
        int sign = **paras > 0 ? 1 : -1;
        int a, b;
        a = *paras[ptr] / exp_digit;
        b = *paras[ptr] % exp_digit * sign;
        val = std::to_string(abs(a)) + "E" + std::to_string(b);
    }
    else if((ptr >= 12 && ptr <= 13) || (ptr >= 15 && ptr <= 16)){
        int2str(paras[ptr], *paras[10], val);
        trim(val);
    }
    else{
        printf("\033[33mWARNING\033[0m: parameters[%i] don't need configure", ptr);
    }
};

void set_appoint(int** paras, const std::string& para_name, const std::string& val){
    std::vector<std::string> para_names;
    para_names = {"Exp_Digit", "Alpha", "Initial_ArcLen", "Max_ArcLen_Inc_Times", "Min_ArcLen", "Optimal_Iter_Num",
                  "Max_Iter_Num_In_NR", "Max_Iter_Num", "Tolerance", "Mode", "Str_Len", "Write_Equi_Path_File",
                  "Path_to_Equi_Path_File", "Name_Of_Equi_Path_File", "Write_Disp_File",
                  "Path_to_Disp_File", "Name_Of_Disp_File", "Element_Type", "Disturb_Order",
                  "Disturb_Magnitude", "Width", "Length"};
    int ptr = -1;
    for(int i = 0; i < para_names.size(); i++){
        if(!_strcmpi(para_name.data(), para_names.at(i).data())){
            ptr = i;
            break;
        }
    }

    if(ptr >= 0){
        set_seq(paras, ptr, val);
    }
    else{
        printf("\033[1;31mERROR\033[0m: no such parameter [%s]\n", para_name.data());
    }
};

void show(std::vector<int>& vec, std::vector<int>& opt_vec){
    std::vector<std::string> para_names;
    para_names = {"Exp_Digit             ", "Alpha                 ", "Initial_ArcLen        ",
                  "Max_ArcLen_Inc_Times  ", "Min_ArcLen            ", "Optimal_Iter_Num      ",
                  "Max_Iter_Num_In_NR    ", "Max_Iter_Num          ", "Tolerance             ",
                  "Mode                  ", "Str_Len               ", "Write_Equi_Path_File  ",
                  "Path_to_Equi_Path_File", "Name_Of_Equi_Path_File", "Write_Disp_File       ",
                  "Path_to_Disp_File     ", "Name_Of_Disp_File     ", "Element_Type          ",
                  "Disturb_Order         ", "Disturb_Magnitude     ", "Width                 ", 
                  "Length                "};
    
    printf("Remain parameters to configure:\n");
    int j = 1;
    for(std::vector<int>::iterator ipara = vec.begin(); ipara != vec.end(); ipara++){
        printf("\033[1;36m%s\033[0m, ", para_names.at(*ipara).data());
        if(j == 4){
            printf("\n");
            j = 0;
        }
        j++;
    }

    printf("\nOptional parameters to configure:\n");
    j = 1;
    for(std::vector<int>::iterator ipara = opt_vec.begin(); ipara != opt_vec.end(); ipara++){
        printf("\033[1;36m%s\033[0m, ", para_names.at(*ipara).data());
        if(j == 4){
            printf("\n");
            j = 0;
        }
        j++;
    }
};

bool read_cfg(int** paras){
    std::ifstream cfg_fio("../../nl_demo.cfg");
    std::string wkstr;
    std::string wkwd;
    std::string wd_front, wd_back;
    bool is_done = false;
    bool is_end = false;
    bool mode_set = false;
    int mode = 0; //0 - appoint; 1 - sequence
    int where;
    int para_ptr = 0;

    is_end = getlmsg(cfg_fio, wkstr);
    while(!is_done && !is_end){
        trim(wkstr);

        while(wkstr.size() > 0){
            first_wd(wkstr, wkwd);

            if(wkwd.size() > 0){
                bool is_ap = split(wkwd, wd_front, wd_back, '=');

                if(wd_back.size() == 0 && !is_ap){
                    if(!mode_set){
                        mode = 1;
                        mode_set = true;
                    }
                    else{
                        if(mode == 0){
                            printf("\033[1;31mERROR\033[0m: unexpected configure command %s\n", (wd_front + "=").data());
                        }
                    }
                    set_seq(paras, para_ptr, wd_front);
                    para_ptr++;
                }
                else if(wd_back.size() != 0 && is_ap){
                    if(!mode_set){
                        mode = 0;
                        mode_set = true;
                    }
                    else{
                        if(mode == 1){
                            printf("\033[1;31mERROR\033[0m: expected configure command ${Para's name}=%s\n", wd_front.data());
                        }
                    }
                    set_appoint(paras, wd_front, wd_back);
                    para_ptr++;
                }
                else{
                    printf("\033[1;31mERROR\033[0m: need an input after %s=\n", wd_front.data());
                }
            }
            else{
                continue;
            }

            if(para_ptr > 21){
                break;
            }
        }

        if(para_ptr > 21){
            is_done = true;
        }
        is_end = getlmsg(cfg_fio, wkstr);
    }
    cfg_fio.close();
    return is_done;
};

int main(){
    char buffer[64]{0};
    std::vector<std::string> para_names, para_names_ws;
    para_names = {"Exp_Digit", "Alpha", "Initial_ArcLen", "Max_ArcLen_Inc_Times", "Min_ArcLen", "Optimal_Iter_Num",
                  "Max_Iter_Num_In_NR", "Max_Iter_Num", "Tolerance", "Mode", "Str_Len", "Write_Equi_Path_File",
                  "Path_to_Equi_Path_File", "Name_Of_Equi_Path_File", "Write_Disp_File",
                  "Path_to_Disp_File", "Name_Of_Disp_File", "Element_Type", "Disturb_Order",
                  "Disturb_Magnitude", "Width", "Length"};

    para_names_ws = {"Exp_Digit             ", "Alpha                 ", "Initial_ArcLen        ",
                     "Max_ArcLen_Inc_Times  ", "Min_ArcLen            ", "Optimal_Iter_Num      ",
                     "Max_Iter_Num_In_NR    ", "Max_Iter_Num          ", "Tolerance             ",
                     "Mode                  ", "Str_Len               ", "Write_Equi_Path_File  ",
                     "Path_to_Equi_Path_File", "Name_Of_Equi_Path_File", "Write_Disp_File       ",
                     "Path_to_Disp_File     ", "Name_Of_Disp_File     ", "Element_Type          ",
                     "Disturb_Order         ", "Disturb_Magnitude     ", "Width                 ", 
                     "Length                "};

    int method;
    printf("The method to use:\n(\033[35m0-GALM\033[0m: General Arc Length Method;\n \033[35m1-MLSM\033[0m: Multiple Load Step Method)\n");
    std::cin >> method;
    std::string m_name = method ? "MLSM" : "GALM";
    asb_manager asb("SN", m_name, 2);
    //asb_manager asb("SN", "BAEA", 2);
    int wide_elenum = 2;
    printf("input wide element number:\n");
    std::cin >> wide_elenum;
    int bnd_type = 22;
    printf("input boundary condition:\n(\033[35m02\033[0m - one side anchorage and one side free;\n \
\033[35m11\033[0m - Simply supported at both ends;\n \033[35m12\033[0m - one side anchorage and one side simply supported;\n \
\033[35m22\033[0m - anchoraged in two sides)\n");
    std::cin >> bnd_type;

    // -- default parameters
    std::string exp_digits = "-2", alpha = "0.E3", dl0 = "5.0E1", max_inc_t = "1.2", dl_min = "0.5E-2",
        r0 = "14", max_itr = "16", max_incr = "40", tol = "1.E-8", mode = "1", str_len = "64",
        plot = "1", write = "1";
    std::string ipath = "", iname_cur = "", iname_dis = "", ipath2 = "";
    // std::string path = ".";
    // ipath = new int[str_len]{0};
    // for(int i = 0; i < str_len; i++){
    //     ipath[i] = 32;
    // };
    // str2int(path, str_len, ipath);

    std::string ele_type = "1", order = "1";
    std::string disturb_magni = "1E-4";        // 1E-4
    std::string side_len = "1.";             // 1
    std::string len = "10.";                 // 10
    std::vector<std::string> paras_dft = {exp_digits, alpha, dl0, max_inc_t, dl_min, r0, max_itr, max_incr, tol, mode,
                          str_len, plot, ipath, iname_cur, write, ipath2, iname_dis, ele_type, order, disturb_magni,
                          side_len, len};
    // -- default parameters

    int *solver_paras[64]{nullptr};

    bool ask_for_cfg = read_cfg(solver_paras);

    std::vector<int> non_cfg_list, opt_cfg_list;
    for(int i = 0; i < 22; i++){
        if(!solver_paras[i]){
            if(i == 12 || i == 13 || i == 15 || i == 16){
                opt_cfg_list.push_back(i);
            }
            else{
                non_cfg_list.push_back(i);
            }
        }
    }

    char set_in_cmd = 'n', confirm = 'n';
    if(non_cfg_list.size() != 0){
        printf("There are some parametres do not configure in 'nl_demo.cfg' file\n\
Do you want set these parameters in cmdline?(\033[32mY\033[0m/\033[31mn\033[0m)\n");
        std::cin >> set_in_cmd;
    }

    std::string input_msg, wd_f, wd_b, wd;
    if(set_in_cmd == 'Y' || set_in_cmd == 'y'){
        printf("Please configure [\033[36m%s\033[0m]:\n\
(or you can use {para_name}={val} to configure specified parameter;\n\
 input \033[1;33muncfg\033[0m to view remain parameters need to be configured;\n\
 input \033[1;33mexit\033[0m to exit parameter configuration;\n\
 input \033[1;33mshow\033[0m to show current parameters' value;\n\
 input \033[1;33mhelp\033[0m to show this message again)\n", para_names.at(non_cfg_list.front()).data());

        std::cin.getline(buffer, 64, '\n');
        while(non_cfg_list.size() > 0){
            // std::cin >> std::noskipws >> input_msg;
            printf("--------------------------------------------------------------------------------v[Don't go past here]\n");
            std::cin.getline(buffer, 64, '\n');
            input_msg.clear();
            input_msg.resize(64, ' ');
            for(int i = 0; i < 64; i++){
                if(buffer[i] != '\n' && buffer[i] != 0) {
                    input_msg.at(i) = buffer[i];
                }
                else{
                    break;
                }
            }
            //for (int i = 0; i < 64; i++) {
            //    buffer[i] = 0;
            //}

            trim(input_msg);
            if(!_strcmpi(input_msg.data(), "uncfg")){
                show(non_cfg_list, opt_cfg_list);
                printf("\n");
            }
            else if(!_strcmpi(input_msg.data(), "exit")){
                break;
            }
            else if(!_strcmpi(input_msg.data(), "help")){
                printf("(input \033[1;33muncfg\033[0m to view remain parameters need to be configured;\n\
 input \033[1;33mexit\033[0m to exit parameter configuration;\n\
 input \033[1;33mshow\033[0m to show current parameters' value;\n\
 input \033[1;33mhelp\033[0m to show this message again)\n");
            }
            else if(!_strcmpi(input_msg.data(), "show")){
                int j = 1;
                std::string val;
                for(int i = 0; i < 22; i++){
                    show_seq(solver_paras, i, val);
                    for(int k = val.size(); k < 8; k++){
                        val.append(" ");
                    }
                    printf("\033[36m%s\033[0m: \033[32m%s\033[0m,  ", para_names_ws.at(i).data(), val.data());
                    if(j == 3){
                        printf("\n");
                        j = 0;
                    }
                    j++;
                }
                printf("\n");
            }
            else{
                while(input_msg.size() > 0){
                    first_wd(input_msg, wd);
                    
                    bool is_speci = split(wd, wd_f, wd_b, '=');
                    if(wd_b.size() == 0 && !is_speci){
                        if(set_seq(solver_paras, non_cfg_list.front(), wd_f)){
                            non_cfg_list.erase(non_cfg_list.begin());
                        }
                    }
                    else if(wd_b.size() != 0 && is_speci){
                        int ptr = -1;
                        int where = -1;
                        int where_opt = -1;
                        for(int i = 0; i < para_names.size(); i++){
                            if(!_strcmpi(wd_f.data(), para_names.at(i).data())){
                                ptr = i;
                                break;
                            }
                        }
                        for(int i = 0; i < non_cfg_list.size(); i++){
                            if(ptr == non_cfg_list.at(i)){
                                where = i;
                                break;
                            }
                        }
                        if(where < 0){
                            for(int i = 0; i < opt_cfg_list.size(); i++){
                                if(ptr == opt_cfg_list.at(i)){
                                    where_opt = i;
                                    break;
                                }
                            }
                        }

                        if(where >= 0 && ptr >= 0){
                            if(set_seq(solver_paras, ptr, wd_b)){
                                non_cfg_list.erase(non_cfg_list.begin() + where);
                            }
                        }
                        else if(where_opt >= 0){
                            if(set_seq(solver_paras, ptr, wd_b)){
                                opt_cfg_list.erase(opt_cfg_list.begin() + where_opt);
                            }
                        }
                        else if(ptr < 0){
                            printf("\033[1;31mERROR\033[0m: no such parameter [\033[36m%s\033[0m]\n", wd_f.data());
                        }
                        else{
                            printf("parameter [%s] is already configured, want to reconfigure?(\033[32mY\033[0m/\033[31mn\033[0m)\n", wd_f.data());
                            std::cin >> confirm;
                            std::cin.getline(buffer, 64, '\n');
                            if(confirm == 'Y' || confirm == 'y'){
                                set_seq(solver_paras, ptr, wd_b);
                            }
                        }
                    }
                    else{
                        printf("\033[1;31mERROR\033[0m: need an input after '\033[35m%s=\033[0m'\n", wd_f.data());
                    }
                }
            }

            if(non_cfg_list.size() > 0) printf("Please configure [\033[36m%s\033[0m]:\n", para_names.at(non_cfg_list.front()).data());
        }
    }

    if(non_cfg_list.size() > 0){
        printf("There are still some parametres do not configure\n\
program will used default values\n");

        for(std::vector<int>::iterator itncl = non_cfg_list.begin(); itncl != non_cfg_list.end(); itncl++){
            delete[] solver_paras[*itncl]; solver_paras[*itncl] = nullptr;
            set_seq(solver_paras, *itncl, paras_dft[*itncl]);
        }
    }

    printf("Current parameters:\n");
    int j = 1;
    std::string val;
    for(int i = 0; i < 22; i++){
        show_seq(solver_paras, i, val);
        for(int k = val.size(); k < 8; k++){
            val.append(" ");
        }
        printf("\033[36m%s\033[0m: \033[32m%s\033[0m,  ", para_names_ws.at(i).data(), val.data());
        if(j == 3){
            printf("\n");
            j = 0;
        }
        j++;
    }

    printf("\nStart to solve?(\033[32mY\033[0m/\033[31mn\033[0m)?\n");
    std::cin >> confirm;

    if(confirm == 'Y' || confirm == 'y'){
        int exp_digit = pow(10, abs(**solver_paras));
        int sign = **solver_paras > 0 ? 1 : -1;
        double a, b;
        a = *solver_paras[19] / exp_digit;
        b = *solver_paras[19] % exp_digit;
        double db = abs(a) * pow(10, sign * b);

        a = *solver_paras[20] / exp_digit;
        b = *solver_paras[20] % exp_digit;
        double seca = abs(a) * pow(10, sign * b);

        a = *solver_paras[21] / exp_digit;
        b = *solver_paras[21] % exp_digit;
        double l = abs(a) * pow(10, sign * b);
        if (method || *solver_paras[9] == 0) db = 0;

        asb.init_mesh(wide_elenum, *solver_paras[17], bnd_type, db, seca, l, *solver_paras[18]);
        asb.solve(solver_paras);
    }
    else{
        printf("exit\n");
    }
    //for (int i = 0; i < asb.KF.K_sparse.rows; i++) {
    //    printf("% .10f ", asb.sln.Var[i]);
    //    if (j == (ele_type + 1) * wide_elenum + 1) {
    //        printf("\n");
    //        j = 0;
    //    }
    //    j++;
    //}

    
}