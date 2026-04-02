# define COLORTEXT "YES" 
# define REDUCEFACTOR 1
# define PI 3.1415926
# define MAXCHAR 2048
# define MAXATOM 128
# define MAXBOND 512
# define MAXRING 500 
# define MAXGAS  500 /*maximum gasiteger parameters*/
/*
For MAXRING, no dynamic memory applied since the actuall number is determined 
using a recursive function. However, for small and middle-sized molecules, 
it is unlikely that the ring num is larger than 1000
*/
# define MAXCYCLE 100
# define OUTPUTSTEP 10
# define MAXTWIST 10
# define ECSLONG 2
# define COSCUT 120
# define DEGRAD 3.1415926/180
# define VDWIDIST 10
# define ESIDIST 14
# define THETACUT 15
# define CUBE 2.0
# define MAXWILDATOM 20
# define MAXSCHAIN 100 
# define MAXCES 20
# define MAXBEED 20
# define MAXATOMTYPE 250
# define MAXVASTATE 2000
# define PSCUTOFF  7 
# define MAX_CES_BOND 100
