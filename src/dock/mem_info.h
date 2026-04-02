#ifndef _MEM_INFO_H_
#define _MEM_INFO_H_

#include <stdio.h>
#include <stdint.h>
#include "string.h"
#include "sys/types.h"


#ifndef __APPLE__
#include <sys/sysinfo.h>
#endif


int parseMemLine(char*);
int getVirtValue();
int getPhysValue();


#endif  
