/*__________________________________________________________________

	Title: Gray-Scott Reaction Diffusion Model
	''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
	Author: Jolyon Troscianko
	Date: 8/1/2018
............................................................................................................

This code generates Reaction Diffusion patterns.

The algorithm joins opposite edges, meaning the patterns can be tiled, and the growth
is effectively on a torus-type surface.

The starting condition is: all "A" concentrations to 1, and all "B" concentrations to 0,
except for a randomly located vertical line 1 px thick, where "B" concentrations are 1.

The concentrations of A and B are limited between 0 and 1.

Random noise is added to A and B, which helps stabilise pattern growth with low values,
and at high values breaks up the patterns

The patterns grow from their random start column to their surroundings. Once the
patterns reach the opposite end of the image the growth continues for a length
of time dependant on the width of the image. Small images require additional extra
growth time for the patterns to start to stabilise. An exponential decay is used to 
calculate the number of additional steps y = 5000*0.996^w where w is image width px.

Useful descriptions:
http://www.algosome.com/articles/reaction-diffusion-gray-scott.html
http://www.mrob.com/pub/comp/xmorphia/

Useful example code:
https://codepen.io/eskimoblood/pen/sFKkc


____________________________________________________________________
*/

import ij.*;
import ij.plugin.PlugIn;
import ij.process.*;
import ij.gui.*;


public class Reaction_Diffusion_Pattern_Generator_gradient_smlap_adaptivek implements PlugIn {
public void run(String arg) {

IJ.showStatus("Reaction Diffusion");

int w = 1000;
int h = 1000;
int nSteps = 3000; // maximum number of steps
double randomFeed = 0.01;

double fy =  0.0;
double fmin =  0.01;
double fmax = 0.101;

double kx = 0.0;
double kmin = 0.045;
double kmax = 0.065;

double Da = 0.2;
double Db = Da/2;

int recSlice = 100;

/*

Keeping k around 0.77 works well.
f = 0.020 = fail
f = 0.0215 = spots (homogeneous, fails with small images)
f = 0.022 = spots (wobbly)
f = 0.023 = stripes (less connected)
f = 0.024 = stripes (quite connected)
f = 0.025 = inverse spots (only form when image is large enough, around 200)
f = 0.026 = chaotic markings
f = 0.027 = mottle (due to random values)
above with noise around 0.01, values below this don't seem to make a difference
highest noise should  be around 0.2, 0.05 is intermediate/high, 0.01 is intermediate, 0.001 is low

*/


//-------f-k modelling
//Formula: y = a+bx+cx^2+dx^3+ex^4

double akmin = 0.01144;
double bkmin = 3.57996;
double ckmin = -53.74636;
double dkmin = 553.13123;
double ekmin = -2454.95711;

double akmax = 0.035563;
double bkmax = 2.95502;
double ckmax = -49.99330;
double dkmax = 625.35638;
double ekmax = -3263.73676;




GenericDialog gd = new GenericDialog("Reaction Diffusion Pattern Generation");
		gd.addMessage("Pattern Generation Values");
		gd.addNumericField("fmin Feed rate", fmin, 5);
		gd.addNumericField("fmax Feed rate", fmax, 5);
		//gd.addNumericField("kmin Kill rate", kmin, 5);
		//gd.addNumericField("kmax Kill rate", kmax, 5);
		gd.addNumericField("Da Diffusion rate", Da, 5);
		gd.addNumericField("Db Diffusion rate", Db, 5);
		gd.addNumericField("noise", randomFeed, 5);
		gd.addNumericField("Max steps", nSteps, 0);
		gd.addMessage("Image Dimensions");
		gd.addNumericField("w Width", w, 0);
		gd.addNumericField("h Height", h, 0);
		gd.addNumericField("Output slice frequency", recSlice, 0);
gd.showDialog();
	if (gd.wasCanceled())
		return;

	
fmin = gd.getNextNumber();
fmax = gd.getNextNumber();
//kmin = gd.getNextNumber();
//kmax = gd.getNextNumber();

Da = gd.getNextNumber();
Db = gd.getNextNumber();
randomFeed = gd.getNextNumber();
nSteps = (int) Math.round( gd.getNextNumber());
w = (int) Math.round( gd.getNextNumber());
h = (int) Math.round( gd.getNextNumber());
recSlice = (int) Math.round( gd.getNextNumber());


//---------------------Set up environment------------------

int dimension = w*h;

double[] outA = new double[dimension];
double[] outB = new double[dimension];
double[] A = new double[dimension];
double[] B = new double[dimension];

// initialise image
for(int i=0; i<dimension; i++){
	A[i] = 1;
	B[i] = 0;
}

/*
// draw a vertical line in a random column of the image as a starting point
//int midRow = (int) (Math.round(Math.random()*(w-1)));
int midRow = Math.round(w/2);

for(int y=0; y<h; y++)
	B[(y*w)+midRow] = 1;

if(midRow >= Math.round(w/2)) // prepare midRow for finding the opposite side of the image (i.e last pixels to be filled)
	midRow = midRow - Math.round(w/2);
else midRow = midRow + Math.round(w/2);


//for(int i=0; i<dimension; i++)
//	B[i] = Math.round(Math.random()*Math.random());

*/


int randLoc = 0;

// random start points
//for(int i=0; i<(w+h)/2; i++){
for(int i=0; i<(w+h); i++){

	randLoc = (int) Math.round(Math.random()*(dimension-(2*h)-(2*w))+1+w);

	B[randLoc] = 1.0;

	B[randLoc - w] = 1.0;
	B[randLoc - 1] = 1.0;
	B[randLoc + w] = 1.0;
	B[randLoc + 1] = 1.0;

	B[randLoc - w -1] = 0.2;
	B[randLoc - w +1] = 0.2;
	B[randLoc + w -1] = 0.2;
	B[randLoc + w +1] = 0.2;

/*
	B[randLoc - w] = 0.2;
	B[randLoc - 1] = 0.2;
	B[randLoc + w] = 0.2;
	B[randLoc + 1] = 0.2;

	B[randLoc - w -1] = 0.05;
	B[randLoc - w +1] = 0.05;
	B[randLoc + w -1] = 0.05;
	B[randLoc + w +1] = 0.05;
*/
}



double Aval = 0.0;
double Bval = 0.0;

double convCount = 0.0; // used to work with image edges
double cA = 0.0; // convolution sum A
double cB = 0.0; // convolution sum A

double maxB = 0.0;
double minB = 10E10;
double xBSum = 0.0;
int rowsFilled = 0;
int finalFlag = 0;

int pxid = 0;
double d2 = 0.0;
double randVal = 0.0;

double fyMOD1=350;
double fyMOD2=1.315;
double fyMOD3=4.2;

double fyOLD=0;

double fyDIF=0;

int recSliceCounter = 0;

ImageStack outStack = new ImageStack(w, h);

for(int j=0; j<nSteps; j++){

recSliceCounter ++;

for (int y=1; y<h-1; y++){



	fy = (double) (fmin + (fmax-fmin)*((y*1.0)/(h*1.0)));

  fyDIF = (double)  (fmax-fmin)- (fy-fmin);

  fy = (double) (fmin + (fmax-fmin)*((Math.sqrt(y*fyMOD1))/(h*1.0)));

  fy =(double) fy*(fyMOD2-(fyMOD3*fyDIF));



fyOLD=(double)fy;


kmin = akmin+bkmin * fy+ckmin * Math.pow(fy,2)+dkmin * Math.pow(fy,3) + ekmin * Math.pow(fy,4);

kmax = akmax+bkmax * fy+ckmax * Math.pow(fy,2)+dkmax * Math.pow(fy,3) + ekmax * Math.pow(fy,4);


//kmin = -16.6329*fy*fy+2.77755*fy+0.013256;
//kmax = -11.1913421*fy*fy+2.181609*fy+0.04310596;


for (int x=1; x<w-1; x++){

	//Laplacian Convolution

	cA = 0.0;
	cB = 0.0;
	pxid = (y*w)+x;
	convCount = 0.0;
/*
	//TOP
		cA += A[pxid - w];
		cB += B[pxid - w];

	//LEFT
		cA += A[pxid - 1];
		cB += B[pxid - 1];

	//RIGHT
		cA += A[pxid + 1];
		cB += B[pxid + 1];

	//BOTTOM
		cA += A[pxid + w];
		cB += B[pxid + w];

	// CENTRE
	Aval = A[pxid];
	Bval = B[pxid];

*/
	//TOP-LEFT
		cA += 0.25*A[pxid - w - 1];
		cB += 0.25*B[pxid - w - 1];

	//TOP
		cA += A[pxid - w];
		cB += B[pxid - w];

	//TOP-RIGHT
		cA += 0.25*A[pxid - w + 1];
		cB += 0.25*B[pxid - w + 1];

	//LEFT
		cA += A[pxid - 1];
		cB += B[pxid - 1];

	//RIGHT
		cA += A[pxid + 1];
		cB += B[pxid + 1];

	//BOTTOM-LEFT
		cA += 0.25*A[pxid + w - 1];
		cB += 0.25*B[pxid + w - 1];

	//BOTTOM
		cA += A[pxid + w];
		cB += B[pxid + w];

	//BOTTOM-RIGHT
		cA += 0.25*A[pxid + w + 1];
		cB += 0.25*B[pxid + w + 1];


	// CENTRE


	Aval = A[pxid];
	Bval = B[pxid];
	
	d2 = Aval*Bval*Bval;

	kx = kmin + (kmax-kmin)*((x*1.0)/(w*1.0));

	outA[pxid] =    Aval + ((Da* ( cA - 5 * Aval) - d2) + fy * (1 - Aval)) ;
	outB[pxid] =    Bval + ((Db* ( cB - 5 * Bval) + d2) - kx *Bval) ;
	
	//outA[pxid] =    Aval + ((Da* ( cA - 4 * Aval) - d2) + fy * (1 - Aval)) ;
	//outB[pxid] =    Bval + ((Db* ( cB - 4 * Bval) + d2) - kx *Bval) ;

/*
	A[pxid] =    Aval + ((Da* ( cA - 5 * Aval) - d2) + fy * (1 - Aval)) ;
	B[pxid] =    Bval + ((Db* ( cB - 5 * Bval) + d2) - kx *Bval) ;
	
	//A[pxid] =    Aval + ((Da*  cA - d2) + fy * (1 - Aval)) ;
	//B[pxid] =    Bval + ((Db*  cB + d2) - kx *Bval) ;
	
	if(A[pxid] > 1.0)
		A[pxid] = 1.0;
	else if(A[pxid] < 0.0)
		A[pxid] = 0.0;

	if(B[pxid] > 1.0)
		B[pxid] = 1.0;
	else if(B[pxid] < 0.0)
		B[pxid] = 0.0;	



	A[pxid] +=  randomFeed*(Math.random()-0.5);
	if(j>nSteps-1) // add noise to B channel except fot the last go
		B[pxid] += randomFeed*(Math.random()-0.5);
*/


}//x
}//y


//xBSum = 0.0;

for (int y=0; y<h; y++){
for (int x=0; x<w; x++){

	pxid = (y*w)+x;

	if(outA[pxid] > 1.0)
		A[pxid] = 1.0;
	else if(outA[pxid] < 0.0)
		A[pxid] = 0.0;	
	else A[pxid] = outA[pxid];

	if(outB[pxid] > 1.0)
		B[pxid] = 1.0;
	else if(outB[pxid] < 0.0)
		B[pxid] = 0.0;	
	else B[pxid] = outB[pxid];

	A[pxid] +=  randomFeed*(Math.random()-0.5);
	if(j>nSteps-1) // add noise to B channel except fot the last go
		B[pxid] += randomFeed*(Math.random()-0.5);

	xBSum += B[pxid];

}//x
}//y



//if(xBSum/dimension < 0.01/dimension) // stop when the pattern growth has failed
//	j = nSteps;

//--------------is image filled yet?---------------
/*
if(finalFlag == 0){
	xBSum = 0.0;
	for (int y=0; y<h; y++)
		xBSum += B[(y*w)+midRow];

	xBSum = xBSum/h;

	if(xBSum > 0.1){ // final row filled
		// use equation for exponential decay to increase the numebr of additional passes
		// with smaller images requiring additional extra passes.
		nSteps = (int) (j + Math.round(    Math.pow(0.996,w)*5000));

		finalFlag = 1;	
	}

}// final flag
*/

IJ.showProgress((float) j/nSteps);


if(recSliceCounter >= recSlice){
	float[] Bfloat = new float[dimension];
	for(int i=0; i<dimension; i++)
		Bfloat[i] = (float) (B[i]);

	outStack.addSlice("Steps: " + (j+1), Bfloat);
	recSliceCounter = 0;
}


}//j steps

IJ.showProgress((float) 1);

//new ImagePlus("fmin"+ fmin + " fmax"+ fmax + " kmin"+ kmin + " kmax"+ kmax + " rand" + randomFeed, outStack).show();
new ImagePlus("fmin"+ fmin + " fmax"+ fmax + " rand" + randomFeed, outStack).show();


//-------create FK map--------------

float[] kmap = new float[dimension];
float[] fmap = new float[dimension];

ImageStack fkStack = new ImageStack(w, h);

for (int y=1; y<h-1; y++){

	fy = (double) (fmin + (fmax-fmin)*((y*1.0)/(h*1.0)));

	fyDIF = (double)  (fmax-fmin)- (fy-fmin);

	fy = (double) (fmin + (fmax-fmin)*((Math.sqrt(y*fyMOD1))/(h*1.0)));

	fy =(double) fy*(fyMOD2-(fyMOD3*fyDIF));

	kmin = -16.6329*fy*fy+2.77755*fy+0.013256;
	kmax = -11.1913421*fy*fy+2.181609*fy+0.04310596;


	for (int x=1; x<w-1; x++){

		kx = kmin + (kmax-kmin)*((x*1.0)/(w*1.0));

		kmap[x+w*y] = (float) kx;
		fmap[x+w*y] = (float) fy;

	}//x

}//y

fkStack.addSlice("Feedrate", fmap);
fkStack.addSlice("Killrate", kmap);

new ImagePlus("Feed Kill maps", fkStack).show();


} // void
} //class
