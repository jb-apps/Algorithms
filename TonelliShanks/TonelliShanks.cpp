//  Created by Jonathan Benavides Vallejo on 10/04/2014.

#include "TonelliShanks.hpp"

/*
 *  Find a quadratic residue (mod p) of 'a'. p must be an odd prime.
 *
 *  Solve the congruence of the form: x^2 = a (mod p) and returns x. 
 *  Note that p - x is also a root.
 *
 *  0 is returned if no square root exists for these a and p.
 *
 *  The Tonelli-Shanks algorithm is used (except for some simple cases in which the solution
 *  is known from an identity). This algorithm runs in polynomial time
 *  (unless the generalized Riemann hypothesis is false).
 */
void modular_sqrt(mpz_t result, mpz_t a, mpz_t p)
{
    mpz_t aux;
    mpz_init(aux);
    mpz_t aux_2;
    mpz_init(aux_2);
    mpz_powm_ui(aux, a, 4, p);
    
    // Simple Cases
    if (mpz_legendre(a, p) != 1)
    {
        mpz_init_set_ui(result, 0);
        return;
    }
    else if ( mpz_cmp_ui(a, 0) == 0 )
    {
        mpz_init_set_ui(result, 0);
        return;
    }
    else if ( mpz_cmp_ui(p, 2) == 0 )
    {
        mpz_init_set(result, a);
        return;
    }
    else if ( mpz_cmp_ui(aux, 3) == 0 )
    {
        mpz_init_set(aux_2, p);
        mpz_add_ui(aux_2,aux_2,1);
        
        mpz_init_set(aux, aux_2);
        
        mpz_fdiv_q_ui(aux, aux, 4);
        mpz_powm(result, a, aux, p);
        return;
    }
    
    // Partition p-1 to s * 2^e for an odd s
    mpz_init_set(aux_2, p);
    mpz_sub_ui(aux_2, aux_2, 1);
    
    mpz_t s;
    mpz_init_set(s, aux_2);
    
    mpz_t e;
    mpz_init_set_ui(e,0);
    
    mpz_t mod;
    mpz_init(mod);
    mpz_mod_ui(mod, s, 2);
    
    while ( mpz_cmp_ui(mod, 0) == 0)
    {
        mpz_fdiv_q_ui(s, s, 2);
        mpz_add_ui(e, e, 1);
        mpz_mod_ui(mod, s, 2);
    }
    
    // Find some 'n' with a legendre symbol n|p = -1.
    mpz_t n;
    mpz_init(n);
    mpz_init_set_ui(n,2);
    
    while (mpz_legendre(n, p) != -1)
    {
        mpz_add_ui(n,n,1);
    }
    
    // Read the paper "Square roots from 1; 24, 51,
    // 10 to Dan Shanks" by Ezra Brown for more
    // information
    // https://www.maa.org/sites/default/files/pdf/upload_library/22/Polya/07468342.di020786.02p0470a.pdf
    
    mpz_add_ui(s,s,1);
    mpz_set(aux, s);
    mpz_cdiv_q_ui(aux, aux, 2);
    mpz_powm(result, a, aux, p);
    
    mpz_t b;
    mpz_init(b);
    mpz_powm(b, a, s, p);
    
    mpz_t g;
    mpz_init(g);
    mpz_powm(g, a, s, p);
    ///////////////////////////////
    
    mpz_t r;
    mpz_init(r);
    mpz_init_set(r,e);
    
    while (true)
    {
        mpz_t t;
        mpz_init(t);
        mpz_init_set(t,b);
        mpz_t m;
        mpz_init(m);
        
        for (mpz_init_set_ui(m,0); mpz_cmp(m,r) < -1; mpz_add_ui(m,m,1))
        {
            if ( mpz_cmp_ui(t, 1) == 0 ) break;
            mpz_powm_ui(t, t, 2, p);
        }
        
        if (mpz_cmp_ui(m,0) == 0)
        {
            // It returns result with the current value
            return;
        }
        
        mpz_t gs;
        mpz_init(gs);
        mpz_sub(aux_2, r, m);
        mpz_sub_ui(aux_2, aux_2, 1);
        
        mpz_powm_ui(aux, aux, 2, p);
        mpz_powm(gs, g, aux, p);
        
        mpz_t g;
        mpz_init(g);
        mpz_mul(g, gs, gs);
        mpz_mod(g, gs, p);
        
        mpz_mul(result, result, gs);
        mpz_mod(result, result, p);
        
        
        mpz_mul(b, b, g);
        mpz_mod(b, b, p);
        mpz_init_set(r,m);
    }
}