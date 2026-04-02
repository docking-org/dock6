/* call to the C function realloc is necessary as memory is adjusted at
   run time 			S. Sridharan. Sep 95 */
#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 

#ifdef IRIX
#define memalloc memalloc_
#define memcopyit memcopyit_
#elif defined (LINUX)
#define memalloc memalloc_
#define memcopyit memcopyit_
#elif defined (CRAY)
#define memalloc MEMALLOC
#define memcopyit MEMCOPYIT
#elif defined (PC)
#define memalloc (_stdcall MEMALLOC)
#define memcopyit (_stdcall MEMCOPYIT)
#endif
void *memalloc( void **ptr, int *new_size )
{	void *newptr;   

        //printf("new_size = %d\n", new_size);
        //printf("new_size = %d\n", *new_size);

	if(!*new_size) 
        {
	  if(*ptr) free(*ptr);
	  return NULL;
	}
	if(!*ptr) 
        {
          if((*new_size)<0)
          {
            *new_size=-(*new_size);
            newptr=calloc(*new_size,1);
          }else
          {
            //printf("I AM HERE 1\n");
            //printf("size = %d\n",*new_size);
            //newptr=calloc((*new_size)/4,4);
            //newptr=calloc((*new_size)/8,4);
            //newptr=calloc((*new_size),4);
            //newptr=calloc((*new_size),8);
            //printf("%d\n",newptr);
	    /*newptr=malloc(*new_size);*/
	    //newptr=malloc(*new_size*4);
	    //newptr=malloc(*new_size*8);
            //newptr=calloc((*new_size),4);
            newptr=calloc((*new_size),8);
            //printf("%d\n",newptr);
            //printf("I AM HERE 1 e\n");
          }
	} else 
        {
          //printf("I AM HERE 2\n");
          if(*new_size<0) { *new_size=-(*new_size);}
          //printf("size = %d\n",*new_size);
          //*new_size=*new_size*4;
	  //newptr=realloc(*ptr, *new_size);
	  //newptr=realloc(*ptr, *new_size*4);
	  newptr=realloc(*ptr, *new_size*8);
          //printf("I AM HERE 2 e\n");
	}
	if (newptr == 0 && *new_size != 0) 
        {
#ifdef PC
          perror("memalloc");
#else
          //perror("memalloc", 8L);
          perror("memalloc");
#endif
	  exit(EXIT_FAILURE);
	}

	return newptr;
}

void *memcopyit(void **dest, void **src, int *num) {
    void *newptr;
    newptr = memcpy(*dest, *src, *num);
    return newptr;
}
