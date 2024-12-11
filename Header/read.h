#include <string.h>
#include <fstream>
#include <vector>

#ifndef READ_H
#define READ_H

/// @brief read a line of the specified file, skip the comment.
/// @param file_stream specified file stream
/// @param str line string
/// @return is reach end of file
bool getlmsg(std::ifstream &file_stream, std::string &str);

/// @brief read the first word of a string
/// @param str string without first word (output)
/// @param word first word (output)
void first_wd(std::string &str, std::string &word);

/// @brief convet int array to string
/// @param i_arr int array
/// @param len length of i_arr
/// @param str overwrite in str
void int2str(const int* i_arr, const int len, std::string& str);

/// @brief convet string to int array
/// @param i_arr int array
/// @param str overwrite in str
void str2int(const std::string& str, const int len, int* i_arr);

/// @brief delete the front and back space
/// @param str string which to trim
void trim(std::string& str);

/// @brief split a string into two parts from the first position of charater in string
bool split(const std::string& str, std::string& front, std::string& back, char charater = ' ');

#endif