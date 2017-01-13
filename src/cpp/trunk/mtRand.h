/* A C++-program for MT19937: Real number version  (1998/4/6)    */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

/* Period parameters */  
#ifndef mtRand_N

#define mtRand_N 624

#ifndef mtRand_CC
extern unsigned short mtRand_xsubi[3];
#endif

class mtRand
{


/* initializing the array with a NONZERO seed */
public:
    void sgenrand(unsigned long seed);  

    mtRand(unsigned long seed)
    {
        if (seed)
        {
        sgenrand(seed);
        } else
        {
        //sgenrand((unsigned long)nrand48(mtRand_xsubi));
        sgenrand((unsigned long)nrand48(mtRand_xsubi));
        }
    }

    mtRand(void)
    {
        sgenrand((unsigned long)nrand48(mtRand_xsubi));
    }


    double gendrand();

    unsigned long genlrand();
protected:
     unsigned long mt[mtRand_N]; /* the array for the state vector  */
     int mti; /* mti==N+1 means mt[N] is not initialized */
};
#endif
