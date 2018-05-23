/* Program to compute coordinate values at data points
 *
 * Written by cthree-40
 * The Johns Hopkins University
 * NOTE (05-04-2015): much of the this program is hard-coded 
 * for the hydroxymethyl problem. Please be careful.
 *
 * Log:
 *  05-04-2015 Added OOP coordinate computation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MATH_PI 3.14159265358979323846

struct point {
  int          index;
  double      *coord;
  double        rmse;
  double        rmsg;
  struct point *next;
};

struct xyz {
  double     x;
  double     y;
  double     z;
};

struct geometry {
  int                 index;
  struct xyz          *atom;
  struct geometry     *next;
};

/* subfunctions */
void   readinput   ( int natoms, int nrgeoms, struct geometry *refgeoms );
void   computedist ( int  *dist, int nrgeoms, struct geometry *refgeoms,
		     struct point *data );
void   computerms  ( int nrgeoms, struct point *data );
double dist3d      ( struct xyz *vec1, struct xyz *vec2 );
void   printoutput ( int nrgeoms,  struct point *data );
double oop3d       ( struct xyz *vec1, struct xyz *vec2, struct xyz *vec3,
		     struct xyz *vec4 );
struct xyz subtractxyz ( struct xyz *vec1, struct xyz *vec2 );
double dotprod ( struct xyz vec1, struct xyz vec2 );
double angl2vecxyz ( struct xyz vec1, struct xyz vec2 );
struct xyz crossprod ( struct xyz vec1, struct xyz vec2 );

/* main */
main( int argc, char *argv[] )
{

  int                natoms, nrgeoms;
  int                          *dist;
  struct geometry           *refgeom, *ptr2, *tmp2;
  struct point               *ptdata, *ptr1, *tmp1;

  /* process command line arguments */
  if ( argc == 3 ) { /* argc should be 3 */
    natoms  = strtol( argv[1], NULL, 10 ); /* number of atoms */
    nrgeoms = strtol( argv[2], NULL, 10 ); /* number of geoms */
  } else {
    printf( "\n usage: %s [num. atoms] [num. geoms]\n", argv[0] );
    exit(1);
  }
  /* allocate first node of refgeom */
  refgeom = malloc( sizeof( struct geometry ) );
  refgeom->atom = malloc( natoms * sizeof( struct xyz ) );
  refgeom->next = 0;
  /* read in refgeom file, building refgeom array */
  readinput( natoms, nrgeoms, refgeom );
  
  /* allocate distances array */
  dist = malloc( (5 * 2) * sizeof( int ) );
  dist[0] = 1;
  dist[1] = 2; /* O-H1  bond distance */
  dist[2] = 1;
  dist[3] = 3; /* O-H2  bond distance */
  dist[4] = 1;
  dist[5] = 4; /* O-H3 bond distance */
  dist[6] = 2;
  dist[7] = 3; /* H1-H2 bond distance */
  dist[8] = 2; 
  dist[9] = 4; /* H1-H3 bond distance */
  /* allocate first node of data aray */
  ptdata = malloc( sizeof( struct point ) );
  /* 5 distances + 2 OOP angles */
  ptdata->coord = malloc( (5 + 2) * sizeof( double ) );
  ptdata->next = 0;
  /* compute distances */
  computedist( dist, nrgeoms, refgeom, ptdata );
  /* compute rms energy and gradient errors */
  computerms( nrgeoms, ptdata );
  /* print distances */
  printoutput( nrgeoms, ptdata );
  
  /* dellocate memory */
  ptr1 = ptdata;
  while ( ptr1 != NULL ) {
    tmp1 = ptr1;
    ptr1 = ptr1->next;
    free( tmp1 );
  }
  ptdata = NULL;
  ptr2 = refgeom;
  while (ptr2 != NULL ) {
    tmp2 = ptr2;
    ptr2 = ptr2->next;
    free( tmp2 );
  }
  refgeom = NULL;
}

void readinput( int natoms, int nrgeoms, struct geometry *rgeom )
{
  /* read input file, generating rgeom array */
  FILE                       *rgeomfl;
  char                    scrstr[100];
  struct geometry             *crgeom;
  int                       i, gindex;

  /* set crgeom pointer */
  crgeom = rgeom;
  
  /* open refgeom file */
  rgeomfl = fopen( "refgeom", "r" );
  if ( rgeomfl == NULL ) { /* file could not be opened */
    printf( " Could not open refgeom. Exiting...\n " );
    exit(1);
  }

  /* read to end of geom file */
  /* geometries are separated by empty line */
  gindex = 1;
  while ( fscanf( rgeomfl, "%lf %lf %lf\n", &crgeom->atom[0].x,
	       &crgeom->atom[0].y, &crgeom->atom[0].z ) != EOF ) {
    for ( i=1; i < natoms; i++ ) {
      fscanf( rgeomfl, "%lf %lf %lf\n", &crgeom->atom[i].x,
	      &crgeom->atom[i].y, &crgeom->atom[i].z );
    }

#ifdef debugging
    /* print geometry */
    for ( i = 0; i < natoms; i++ ) {
      printf( " %lf %lf %lf\n", crgeom->atom[i].x, crgeom->atom[i].y,
	      crgeom->atom[i].z );
    }
#endif
       
    /* iterate index */
    gindex = gindex + 1;
    /* set next pointer */
    crgeom->next = malloc( sizeof( struct geometry ) );
    crgeom = crgeom->next;
    crgeom->atom = malloc( natoms * sizeof( struct xyz ) );
    crgeom->next = 0;
  } /* while */
  
  /* check if gindex == nrgeoms */
  gindex = gindex - 1;
  if ( gindex != nrgeoms ) {
    printf( " nrgeoms != gindex : %d != %d \n", nrgeoms, gindex );
  }
  
  /* close file */
  fclose( rgeomfl );
}
void computedist( int *dist, int nrgeoms, struct geometry *rgeoms,
		  struct point *ptdata )
{
  struct point           *currptd;
  struct geometry        *currrgm;
  int                   indexp, i;
  
  /* set pointers */
  currrgm = rgeoms;
  currptd = ptdata;
  indexp = 1;
  /* loop over rgeoms */
  while ( currrgm-> next != 0 ) {
    /* set index */
    currptd->index = indexp;
    /* compute each two atom distance */
    for ( i = 0; i < 5; i++ ) {
      currptd->coord[i] = dist3d( &(currrgm->atom[dist[(i * 2)]-1]),
				  &(currrgm->atom[dist[(i * 2) + 1 ]-1]));
    } /* i */
    /* compute OOP angles */
    /* ( 5, 1, 2, 3 ) */
    currptd->coord[i] = oop3d( &(currrgm->atom[3]), &(currrgm->atom[0]),
			       &(currrgm->atom[1]), &(currrgm->atom[2]));
    i++;
    /* ( 3, 1, 2, 4 ) */
    currptd->coord[i] = oop3d( &(currrgm->atom[3]), &(currrgm->atom[0]),
			       &(currrgm->atom[1]), &(currrgm->atom[2]));
    
#ifdef debugging
    /* print CO bond distances */
    printf( " %d C-O distance: %lf \n", currptd->index, currptd->coord[0] );
#endif

    /* allocate new ptdata node */
    currptd->next  = malloc( sizeof( struct point ) );
    currptd = currptd->next;
    /* 5 distances + 2 OOP angles */
    currptd->coord = malloc( (5 + 2) * sizeof( double ) );
    currptd->next = 0;
    /* move to new geometry */
    currrgm = currrgm->next;
    indexp = indexp + 1;
  } /* while */

}
void computerms( int nrgeoms, struct point *data )
{
  /* compute rms errors for energies and gradiets from fitinfo.csv file. */
  
  FILE            *fitinfo;
  char            *fistream;
  int              nstates, ngrads;
  double           ecutoff, weight;
  double          *enerror;
  double          *gderror;
  int                  i,j;
  struct point      *cdata;

  /* open fitinfo.csv */
  fitinfo = fopen( "fitinfo.csv", "r" );
  if ( fitinfo == NULL ) { /* if file cannot be opened... */
    printf( " fitinfo.csv file cannot be opened. Exiting...\n" );
    exit(1);
  }
  
  /* set cdata pointer to first node of data */
  cdata = data;
  
  /* read first line of fitinfo, to obtain state information
   *  and energy cutoff information. */
  fscanf( fitinfo, " %d, %lf ", &nstates, &ecutoff );
#ifdef debugging
  printf( "%d %lf \n", nstates, ecutoff );
#endif
  /* allocate energy error array */
  ngrads = ( nstates + 1) * nstates / 2 ;
  enerror = malloc( ngrads * sizeof( double ) );

  /* read in rms data for each point */
  fistream = malloc( ngrads * 20 * sizeof( char ) );
  for ( i = 0; i < nrgeoms; i++ ) {
    /* get line */
    fgets( fistream, (ngrads * 20), fitinfo );
    /* read weight */
    sscanf( fistream, "%*d, %lf", &weight );
    if ( weight > 0.0 ) { /* if included in fit */
      // add rms terms
      
    }
      
  }
      
}
double dist3d( struct xyz *vec1, struct xyz *vec2 )
{
  /* compute distance between 2 points in 3D */
  
  double      xdiff, ydiff, zdiff; 
  double            xsq, ysq, zsq;
  double                   result;
  
  /* compute distances */
  xdiff = vec1->x - vec2->x;
  ydiff = vec1->y - vec2->y;
  zdiff = vec1->z - vec2->z;
  
  /* square distances and take square root of sum */
  result = sqrt( pow( xdiff, 2 ) + pow( ydiff, 2 ) +
		 pow( zdiff, 2 ) );
  return result;

}
double oop3d ( struct xyz *vec1, struct xyz *vec2, struct xyz *vec3,
	       struct xyz *vec4 )
{
  /* compute oop angle between plane formed by 4 atoms */
  struct xyz plane1vec1, plane1vec2;
  struct xyz plane2vec1, plane2vec2;
  struct xyz normvec1, normvec2;
  double angle;
  
  /* plane1 is atoms 1, 2, 3
   * plane2 is atoms 2, 3, 4 */
  plane1vec1 = subtractxyz( vec1, vec2 );
  plane1vec2 = subtractxyz( vec3, vec2 );
  plane2vec1 = subtractxyz( vec4, vec3 );

  /* compute normal vectors */
  plane2vec2.x = plane1vec2.x * (-1);
  plane2vec2.y = plane1vec2.y * (-1);
  plane2vec2.z = plane1vec2.z * (-1);
  normvec1 = crossprod( plane1vec1, plane1vec2 );
  normvec2 = crossprod( plane2vec1, plane2vec2 );
  
  /* compute angle between normal vectors */
  angle = angl2vecxyz( normvec1, normvec2 ) * 180.0 / MATH_PI;

  /* return angle */
  return angle;
}
double dotprod( struct xyz vec1 ,struct xyz vec2 )
{
  double result;
  
  result = 0.0;
  
  result = result + vec1.x * vec2.x + vec1.y * vec2.y 
    + vec1.z * vec2.z;
  
  return result;
}
/*>angl2vecxyz
 * Compute angle between two vectors by
 * cos q = v.w / (||v|| ||w||)
 */
double angl2vecxyz( struct xyz vec1, struct xyz vec2 )
{
  double angle, acosin;
  double dotvw, normv, normw;
  
  dotvw = dotprod( vec1, vec2 );
  normv = sqrt( dotprod( vec1, vec1 ) );
  normw = sqrt( dotprod( vec2, vec2 ) );
  
  acosin = dotvw / (normv * normw);
  angle = acos( acosin );

  return angle;
}
/*>crossprod
 * Compute the cross product 
 */
struct xyz crossprod( struct xyz vec1, struct xyz vec2 )
{
  struct xyz result;
  
  result.x = ( vec1.y * vec2.z - vec1.z * vec2.y );
  result.y = ( vec1.z * vec2.x - vec1.x * vec2.z );
  result.z = ( vec1.x * vec2.y - vec1.y * vec2.x );
  
  return result;
}
void printoutput( int nrgeoms, struct point *data )
{
  /* print bond distance information in .csv format */

  FILE         *outputfl;
  struct point    *currd;
  int                  i;
  
  /* open output file */
  outputfl = fopen( "bnddist.csv", "w" );
  if ( outputfl == NULL ) { /* file could not be opened */
    printf( " Could not open output file. Exiting..\n" );
    exit(1);
  }
  
  /* start at first data point */
  currd = data;
  while ( currd->next != 0 ) {
    fprintf( outputfl, "%lf, %lf, %lf, %lf, %lf, %lf, %lf \n",
	     currd->coord[0], currd->coord[1], currd->coord[2],
	     currd->coord[3], currd->coord[4], currd->coord[5],
	     currd->coord[6] );
    currd = currd->next;
  }

  /* close output file */
  fclose(outputfl);

}
/*>subtractxyz
 * Subtract to xyz points 
 */
struct xyz subtractxyz ( struct xyz *vec1, struct xyz *vec2 )
{
  struct xyz result;
  
  result.x = ( vec1->x - vec2->x );
  result.y = ( vec1->y - vec2->y );
  result.z = ( vec1->z - vec2->z );
  
  return result;
}
