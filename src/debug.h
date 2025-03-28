#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>

#ifndef LOG_LEVEL
#define LOG_LEVEL 4  
#endif

#define LOG_LEVEL_NONE   0
#define LOG_LEVEL_ERROR  1
#define LOG_LEVEL_WARN   2
#define LOG_LEVEL_INFO   3
#define LOG_LEVEL_DEBUG  4

#if LOG_LEVEL >= LOG_LEVEL_DEBUG
    #define LOG_DEBUG(...)    printf("[DEBUG] " __VA_ARGS__)
#else
    #define LOG_DEBUG(...)    ((void)0)
#endif

#if LOG_LEVEL >= LOG_LEVEL_INFO
    #define LOG_INFO(...)     printf("[INFO ] " __VA_ARGS__)
#else
    #define LOG_INFO(...)     ((void)0)
#endif

#if LOG_LEVEL >= LOG_LEVEL_WARN
    #define LOG_WARN(...)     printf("[WARN ] " __VA_ARGS__)
#else
    #define LOG_WARN(...)     ((void)0)
#endif

#if LOG_LEVEL >= LOG_LEVEL_ERROR
    #define LOG_ERROR(...)    printf("[ERROR] " __VA_ARGS__)
#else
    #define LOG_ERROR(...)    ((void)0)
#endif

#define debug_printf(...) LOG_DEBUG(__VA_ARGS__)

#endif
