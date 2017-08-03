/*
* T-Dtransverse4velocity.cpp
*	Created on  2017.8.3
*		Author  weiyw
*/
#include <rsf.h>
int main(int argc, char* argv[])
{
	/* parameters defination */
	bool verb;						/* verbose flag */
	sf_file Fin=NULL,Fot=NULL;		/* I/O files */
	sf_axis ax,az;					/* axis */
	int iz,ix;						/* index variables */
	int nz,nx;
	float dt,dz,dx;

	float **vv, **tt;


	sf_init(argc,argv);
	if(! sf_getbool("verb", &verb)) verb=false;

	/* setuo I/O files */
	Fin = sf_input("in");
	Fot = sf_output("out");

	/* prapare variables */
	ax = sf_iaxa(Fin,2); nx = sf_n(ax); dx = sf_d(ax);
	az = sf_iaxa(Fin,1); nz = sf_n(az); dz = sf_d(az);

	sf_oaxa(Fot,az,1);
	sf_oaxa(Fot,ax,2);

	vv = sf_floatalloc2(nz,nx); sf_floatread(vv[0],nz*nx,Fin);
	tt = sf_floatalloc2(nz,nx);

	for(ix=0; ix < nx; ix++) {
		for(iz=0; iz < nz; iz++) {
			tt[ix][iz] = 0.0;
		}
	}

	/* get time */
	for(ix = 0; ix < nx; ix++) {
		for(iz = 1; iz < nz; iz++) {
			tt[ix][iz] = tt[ix][iz-1] + dz/vv[ix][iz];
		}
	}
	
	/* write out result */
	sf_floatwrite(tt[0], nz*nx, Fot);


	/* end */
	if(verb) fprintf(stderr,"/n");
	sf_close();
	exit(0);
	
}
