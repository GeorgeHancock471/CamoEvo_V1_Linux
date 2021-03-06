// Code automatically generated by 'Generate Cone Mapping Model' script by Jolyon Troscianko

//Model fits:
//lw 0.9994274485823442
//mw 0.9999840285946737
//uv 0.9975966112912081


// Generated: 2019/2/18   21:7:29


import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

public class Samsung_NX1000_Nikkor_EL_80mm_D65_to_Jumping_Spider_H_pyrrithrix_D65 implements PlugInFilter {

ImageStack stack;
	public int setup(String arg, ImagePlus imp) { 
	stack = imp.getStack(); 
	return DOES_32 + STACK_REQUIRED; 
	}
public void run(ImageProcessor ip) {

IJ.showStatus("Cone Mapping");
float[] vR;
float[] vG;
float[] vB;
float[] uB;
float[] uR;
int w = stack.getWidth();
int h = stack.getHeight();
int dimension = w*h;

float[] lw = new float[dimension];
float[] mw = new float[dimension];
float[] uv = new float[dimension];

vR = (float[]) stack.getPixels(1);
vG = (float[]) stack.getPixels(2);
vB = (float[]) stack.getPixels(3);
uB = (float[]) stack.getPixels(4);
uR = (float[]) stack.getPixels(5);

for (int i=0;i<dimension;i++) {
lw[i] = (float) (-0.002297610431565501 +(vR[i]*0.008999485172892757)+(vG[i]*-0.0012853952689448902)+(vB[i]*0.0011735971779246726)+(uB[i]*2.0847099486195862E-4)+(uR[i]*9.317843221450548E-4)+(vR[i]*vG[i]*-1.847865788453742E-5)+(vR[i]*vB[i]*2.6002687339557943E-5)+(vR[i]*uB[i]*5.090696751004094E-6)+(vR[i]*uR[i]*-1.217973236191995E-5)+(vG[i]*vB[i]*-7.579625443810254E-6)+(vG[i]*uB[i]*-1.8443416527454413E-5)+(vG[i]*uR[i]*3.099937833567382E-5)+(vB[i]*uB[i]*2.659921662818464E-5)+(vB[i]*uR[i]*-3.126872240526826E-5)+(uB[i]*uR[i]*-9.640664830979268E-7));
mw[i] = (float) (-3.8895528941290864E-4 +(vR[i]*-0.001001243632415447)+(vG[i]*0.00846420305845785)+(vB[i]*0.0020565273768294633)+(uB[i]*-5.051104380676874E-5)+(uR[i]*5.432926297591702E-4)+(vR[i]*vG[i]*-1.3674170095166513E-6)+(vR[i]*vB[i]*1.3549322627180115E-6)+(vR[i]*uB[i]*1.033117166931676E-5)+(vR[i]*uR[i]*-1.2077252198461831E-5)+(vG[i]*vB[i]*-2.1233480357656027E-7)+(vG[i]*uB[i]*-1.5495854569162447E-5)+(vG[i]*uR[i]*1.9388729053977932E-5)+(vB[i]*uB[i]*7.249389191160168E-6)+(vB[i]*uR[i]*-8.788274063288446E-6)+(uB[i]*uR[i]*-5.442491421535315E-7));
uv[i] = (float) (0.0010779158691523844 +(vR[i]*3.630532579682056E-4)+(vG[i]*-8.755033006055652E-4)+(vB[i]*0.002137149421107564)+(uB[i]*0.003787207645764185)+(uR[i]*0.0043409189616341265)+(vR[i]*vG[i]*-3.4123731880978737E-6)+(vR[i]*vB[i]*7.1144789390601836E-6)+(vR[i]*uB[i]*-1.1562264537520436E-5)+(vR[i]*uR[i]*-6.895998742432142E-6)+(vG[i]*vB[i]*-6.822250660117221E-6)+(vG[i]*uB[i]*-7.68734731503483E-5)+(vG[i]*uR[i]*1.0393536458708838E-4)+(vB[i]*uB[i]*1.386752467614376E-4)+(vB[i]*uR[i]*-1.3990346253711557E-4)+(uB[i]*uR[i]*-2.2387694620362043E-6));
IJ.showProgress((float) i/dimension);
}

ImageStack outStack = new ImageStack(w, h);
outStack.addSlice("lw", lw);
outStack.addSlice("mw", mw);
outStack.addSlice("uv", uv);
new ImagePlus("Output", outStack).show();

}
}
