#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>



double user_time(){
 	struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec/1000000.0;
    }

double user_and_system_time(){
 	struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec/1000000.0 + (double)ru.ru_stime.tv_sec + (double)ru.ru_stime.tv_usec/1000000.0;
    }    






