#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Change here the number of steps, the cell geometry, etc
#define NITER 5000
#define STEPITER 1000
#define delx 0.5
#define dely 0.25
// end change here.
/***/
       int main( int argc, char **argv )
{
       int i, j, n, mx1, my1, mx, my;
       FILE *fp;
       double rdx2, rdy2, beta;
/***/
       if ( argc != 3 )
        { 
         fprintf( stderr, "Usage: %s <nrows> <ncolumns>\n", argv[0] );
         return (-1);
          }
       mx = (int)atol( argv[1] );
       my = (int)atol( argv[2] );
       mx += 2;
       my += 2;
/* Here we know the array sizes, so make the arrays themselves: */
       {
        double (*f)[my];
        double (*newf)[my];
        double (*r)[my];
        void *ptmp; 
/***/
        f = (typeof(f))malloc( mx*sizeof(*f) );
          newf = (typeof(newf))malloc( mx*sizeof(*newf) );
        r = (typeof(r))malloc( mx*sizeof(*r) );
        if ( (!f) || (!newf) || (!r) )
         {
          fprintf( stderr, "Cannot allocate, exiting\n" );
          return( -1 );
         }
        rdx2 = 1./delx/delx;
        rdy2 = 1./dely/dely;
        beta = 1.0/(2.0*(rdx2+rdy2)); 
        printf( "Solving task on %d by %d grid\n", mx-2, my-2 );     
        fflush( stdout );
        for ( i = 0; i < mx; i++ )
         {
          for ( j = 0; j < my; j++ )
           {
            if ( (i==0) || (j==0) || (i==(mx-1)) || (j==(my-1)) )
             newf[i][j] = f[i][j] = 1.0;
            else
             newf[i][j] = f[i][j] = 0.0;
            r[i][j] = 0.0;        
           }
         }    
        mx1 = mx - 1;
        my1 = my - 1; 
/* Iteration loop: */
        for ( n = 0; n < NITER; n++ )
         {
           if ( !(n%STEPITER) ) printf( "Iteration %d\n", n );
/* Step of calculation starts here: */ 
          for ( i = 1; i < mx1; i++ )
           {
            for ( j = 1; j < my1; j++ )
             {
              newf[i][j] = ((f[i-1][j]+f[i+1][j])*rdx2
                       +(f[i][j-1]+f[i][j+1])*rdy2-r[i][j])*beta;
             }
            }
/* swap f & newf: */      
             ptmp = f;
              f = newf;
              newf = (typeof(newf))ptmp;
           }    
          fp = fopen("foutput.txt", "w" );
          for ( i = 1; i < (mx-1); i++ )
           {
            //fwrite( &(f[i][1]), my-2, sizeof(f[0][0]), fp );
            for (j = 1; j < (my-1); j++ )
            {
              fwrite(&(f[i][j]), 1, sizeof(f[0][0]), fp);
            }
           } 
          fclose( fp );    
          }   
       return 0;
}
