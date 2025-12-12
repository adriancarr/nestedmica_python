#if defined (NMICA_TEST_ENDIAN)
#if defined (__BIG_ENDIAN__)
#define NMICA_BIG_ENDIAN
#else
#define NMICA_LITTLE_ENDIAN
#endif
#endif

#if defined (__SVR4) && defined (__sun)

#include <values.h>
#define NMICA_NEGATIVE_INFINITY -MAXDOUBLE

#else

#define NMICA_NEGATIVE_INFINITY -INFINITY

#endif

