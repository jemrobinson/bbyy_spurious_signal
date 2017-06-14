#pragma once
#include <iostream>
#include <iomanip>

#define MSG_LVL_NONE      0x00
#define MSG_LVL_ERROR     0x01
#define MSG_LVL_WARNING   0x02
#define MSG_LVL_INFO      0x03
#define MSG_LVL_DEBUG     0x04
#define MSG_LVL_VERBOSE   0x05

#ifndef MSG_LEVEL
#define MSG_LEVEL   MSG_LVL_DEBUG
#endif

#define ERROR_TAG   "  \033[1;31mERROR\033[1;0m"
#define WARNING_TAG "\033[1;35mWARNING\033[1;0m"
#define INFO_TAG    "   \033[1;32mINFO\033[1;0m"
#define DEBUG_TAG   "  \033[1;34mDEBUG\033[1;0m"
#define VERBOSE_TAG "\033[1;34mVERBOSE\033[1;0m"

#define LOGFUNCTION(lvl, message)      std::cout << lvl << ": " << message << std::endl;

#if MSG_LEVEL >= MSG_LVL_VERBOSE
#define MSG_VERBOSE(message)   LOGFUNCTION(VERBOSE_TAG, message)
#else
#define MSG_VERBOSE(message)
#endif

#if MSG_LEVEL >= MSG_LVL_DEBUG
#define MSG_DEBUG(message)     LOGFUNCTION(DEBUG_TAG, message)
#else
#define MSG_DEBUG(message)
#endif

#if MSG_LEVEL >= MSG_LVL_INFO
#define MSG_INFO(message)      LOGFUNCTION(INFO_TAG, message)
#else
#define MSG_INFO(message)
#endif

#if MSG_LEVEL >= MSG_LVL_WARNING
#define MSG_WARNING(message)   LOGFUNCTION(WARNING_TAG, message)
#else
#define MSG_WARNING(message)
#endif

#if MSG_LEVEL >= MSG_LVL_ERROR
#define MSG_ERROR(message)     LOGFUNCTION(ERROR_TAG, message)
#else
#define MSG_ERROR(message)
#endif