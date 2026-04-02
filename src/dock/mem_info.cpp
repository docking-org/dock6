#ifndef __APPLE__

#include <stdio.h>
#include <stdint.h>
#include <sys/sysinfo.h>
#include "string.h"
#include "stdlib.h"
#include <iostream>

// MEMORY FUNCTIONS
int parseMemLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++; 
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

// this is VIRT RAM MEMORY
int getVirtValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseMemLine(line);
            break;
        }
    }    
    fclose(file);
    return result;
}

// this is PHYSICAL RAM memory
int getPhysValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseMemLine(line);
            break;
        }
    }    
    fclose(file);
    return result;
}
#endif
