#ifndef UTIL_SYSTEM_SYSTEM_H_
#define UTIL_SYSTEM_SYSTEM_H_

#include <stdio.h>
#include <string>

enum class Color { RED, GREEN };

void set_color(Color color);
void reset_color();
std::string executable_path();
bool exists(const std::string &file_name);
void auto_append_extension(std::string &str, const char *ext);
void auto_append_extension_if_exists(std::string &str, const char *ext);
size_t getCurrentRSS();
size_t getPeakRSS();
void log_rss();
size_t file_size(const char* name);

#ifdef _MSC_VER
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

#endif