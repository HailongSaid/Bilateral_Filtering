/*****************************************************************************************
\Author:	Qingxiong Yang (liiton)
\Function:	Constant time (joint/cross) bilateral filtering
\paper:     Q. Yang, K.-H. Tan, N. Ahuja, Real-time O(1) Bilateral Filtering, CVPR, 2009. 
*****************************************************************************************/
#ifndef QX_CONSTANT_TIME_BILATERAL_FILTER_H
#define QX_CONSTANT_TIME_BILATERAL_FILTER_H
#define QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER		0
#define QX_DEF_CTBF_BOX_BILATERAL_FILTER			1
#define QX_DEF_CTBF_INTENSITY_RANGE					256
#define QX_DEF_CTBF_SIGMA_SPATIAL_DEFAULT			0.03//0.03~16/512
#define QX_DEF_CTBF_SIGMA_RANGE_DEFAULT				0.08//0.08~20/255
class qx_constant_time_bilateral_filter
{
public:
    qx_constant_time_bilateral_filter();
    ~qx_constant_time_bilateral_filter();
    void clean();
	int init(int h,int w,
		int spatial_filter=QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER,//default: Gaussian BF
		double sigma_spatial=QX_DEF_CTBF_SIGMA_SPATIAL_DEFAULT,//0.03~16/512
		double sigma_range=QX_DEF_CTBF_SIGMA_RANGE_DEFAULT);//0.08~20/255
	int bilateral_filter(unsigned char **image_filtered,unsigned char **image,
		double sigma_range=0,//0.08~20/255
		unsigned char **texture=NULL);
	int bilateral_filter(float **image_filtered,float **image,
		double sigma_range=0,//0.08~20/255
		unsigned char **texture=NULL);
private:
	char m_str[300];
	int m_h,m_w,m_radius; double m_sigma_range,m_sigma_spatial; int m_nr_scale; int m_spatial_filter;
	unsigned char**m_image,**m_image_filtered; double ***m_jk,**m_wk; double **m_box; 
	double *m_grayscale;
	double *m_table;
};
#endif