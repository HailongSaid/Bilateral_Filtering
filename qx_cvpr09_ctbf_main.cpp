#include "qx_ppm.h"
#include "qx_cvpr09_ctbf_basic.h"
#include "qx_constant_time_bilateral_filter.h"
#define QX_DEF_ARGC1								4
#define QX_DEF_ARGC2								6


int main(int argc,char *argv[])
{
	//if(argc!=QX_DEF_ARGC1&&argc!=QX_DEF_ARGC2)
	//{		
	//	printf("Usage 1: \n");
	//	printf("      Gaussian bilateral filter: qx_cvpr09_ctbf.exe 1.pgm girl.pgm  0\n");
	//	printf("      Box bilateral filter: qx_cvpr09_ctbf.exe fileout filein 1 \n\n");
	//	printf("Usage 2: \n");
	//	printf("      Gaussian bilateral filter: qx_cvpr09_ctbf.exe fileout filein 0 sigma_spatial(0 to 1) sigma_range(0 to 1)\n");
	//	printf("      Box bilateral filter: qx_cvpr09_ctbf.exe fileout filein 1 sigma_spatial(0 to 1) sigma_range(0 to 1)\n\n");
	//	getchar();
	//	return(-1);
	//}
	//char *fileout = argv[1];
	//char *filein = argv[2];
	
	//spatial_filter = QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER;
	//argc = QX_DEF_ARGC2;
	//if(argc==QX_DEF_ARGC2)
	//{
	//	if(spatial_filter==QX_DEF_CTBF_BOX_BILATERAL_FILTER) sigma_spatial=atof(argv[4]);
	//	else if(spatial_filter==QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER) sigma_spatial=atof(argv[4]);
	//	else 
	//	{
	//		printf("Note: ONLY support box and Gaussian spatial filter!");
	//		printf("Switching to Gaussian spatial filter automatically.....");
	//		spatial_filter=QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER;
	//		sigma_spatial=atof(argv[4]);
	//	}
	//	sigma_range=sqrt(atof(argv[5]));
	//}

	printf("%lf ", 3 / (double)2.0);
	char *filein = "girl.pgm";
	char *fileout = "1.jpg";
	int spatial_filter = QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER; double sigma_spatial = 0.1, sigma_range = 0.1;

	qx_constant_time_bilateral_filter m_bf; qx_timer timer;
	unsigned char **image,**image_filtered; int h,w;
	image=loadimage_pgm_u(filein,h,w);
	image_filtered=qx_allocu(h,w);
	//if(argc==QX_DEF_ARGC2)
	//{
	//	m_bf.init(h,w,spatial_filter,sigma_spatial,sigma_range); /*initialization*/
	//}
	//else
	//{
	//	m_bf.init(h,w,spatial_filter); /*initialization*/
	//}

	m_bf.init(h, w, spatial_filter, sigma_spatial, sigma_range); /*initialization*/
	timer.start();
	
	m_bf.bilateral_filter(image_filtered,image); /*bilateral filtering*/
	timer.time_display("CTBF");

	saveimage_pgm(fileout,image_filtered,h,w); /*save the filtered image as an unsigned char pgm image*/
	qx_freeu(image); image=NULL;
	qx_freeu(image_filtered); image_filtered=NULL;
	return(0);
} 