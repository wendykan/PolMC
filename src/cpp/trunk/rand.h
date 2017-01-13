
#define InitRandomGen    (double) RandomGen(0, 1, NULL)
     /* Initializes the seed for the random number generator. */     
#define RandomNum        (double) RandomGen(1, 0, NULL)
     /* Calls for a random number from the randum number generator. */

/* DECLARE FUNCTION */
double RandomGen(char Type, long Seed, long *Status);  
     /* Random number generator */

