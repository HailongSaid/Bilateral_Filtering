#include "qx_cvpr09_ctbf_basic.h"
#include "qx_constant_time_bilateral_filter.h"
qx_constant_time_bilateral_filter::qx_constant_time_bilateral_filter()
{
	m_box=NULL;
	m_jk=NULL;
	m_wk=NULL;
	m_grayscale=NULL;
}
qx_constant_time_bilateral_filter::~qx_constant_time_bilateral_filter()
{
	clean();
}
void qx_constant_time_bilateral_filter::clean()
{
	qx_freed(m_box); m_box=NULL;
	qx_freed_3(m_jk); m_jk=NULL;
	qx_freed(m_wk); m_wk=NULL;
	if(m_grayscale!=NULL) delete [] m_grayscale; m_grayscale=NULL;
}
int qx_constant_time_bilateral_filter::init(int h,int w,int spatial_filter,double sigma_spatial,double sigma_range)
{
	m_h=h; m_w=w; m_spatial_filter=spatial_filter; m_sigma_spatial=sigma_spatial; m_sigma_range=sigma_range; m_nr_scale=QX_DEF_CTBF_INTENSITY_RANGE;
	if(m_spatial_filter==QX_DEF_CTBF_BOX_BILATERAL_FILTER) m_radius=int(m_sigma_spatial*qx_min(m_h,m_w)+0.5);
	else if(m_spatial_filter!=QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER) 
	{
		printf("Note: ONLY support box and Gaussian spatial filter!");
		printf("Switching to Gaussian spatial filter automatically.....");
		m_spatial_filter=QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER;
	}
	/*memory allocation*/
	m_box=qx_allocd(h,w);
	m_jk=qx_allocd_3(2,h,w);
	m_wk=qx_allocd(h,w);
	m_table=get_color_weighted_table(m_sigma_range*QX_DEF_CTBF_INTENSITY_RANGE,QX_DEF_CTBF_INTENSITY_RANGE);
	m_grayscale=new double [m_nr_scale];
	return(0);
}
int qx_constant_time_bilateral_filter::bilateral_filter(unsigned char **image_filtered,unsigned char **image,double sigma_range,unsigned char **texture)
{
	unsigned char image_min,image_max; 
	int y,x,jk_0,jk_1;
	if(sigma_range>QX_DEF_THRESHOLD_ZERO) 
	{
		m_sigma_range=sigma_range;
		color_weighted_table_update(m_table,m_sigma_range*QX_DEF_CTBF_INTENSITY_RANGE,QX_DEF_CTBF_INTENSITY_RANGE);
	}
	qx_timer timer;
	timer.start();
	if(texture==NULL)
	{
		vec_min_val(image_min,image[0],m_h*m_w);
		vec_max_val(image_max,image[0],m_h*m_w);
	}
	else
	{
		vec_min_val(image_min,texture[0],m_h*m_w);
		vec_max_val(image_max,texture[0],m_h*m_w);
	}
	m_nr_scale=qx_max(1,int(double(image_max-image_min)/(255*m_sigma_range)+0.5));
	//printf("[qx_max,qx_min]:[%5.5f,%5.5f]\n",(float)image_max,(float)image_min);
	//printf("[sigma_range: %1.3f]\n",m_sigma_range);
	//printf("[nr_scale: %d]\n",m_nr_scale);
	m_grayscale[0]=(double)image_min;
	m_grayscale[m_nr_scale-1]=(double)image_max;
	double delta_scale=double(image_max-image_min)/(m_nr_scale-1);
	for(int i=1;i<m_nr_scale-1;i++) m_grayscale[i]=(double)image_min+delta_scale*i;
	for(int i=0;i<m_nr_scale;i++)
	{
		double **jk;
		if(i==0)
		{
			jk_0=0;
			jk_1=1;
			jk=m_jk[jk_0];
		}
		else 
			jk=m_jk[jk_1];
		for(y=0;y<m_h;y++)
		{
			for(x=0;x<m_w;x++)
			{
				int index;
				if(texture==NULL) index=int(abs(m_grayscale[i]-image[y][x])+0.5f);
				else index=int(abs(m_grayscale[i]-texture[y][x])+0.5f); /*cross/joint bilateral filtering*/
				jk[y][x]=m_table[index]*image[y][x];
				m_wk[y][x]=m_table[index];
			}
		}
		if(m_spatial_filter==QX_DEF_CTBF_BOX_BILATERAL_FILTER)
		{
			boxcar_sliding_window(jk,jk,m_box,m_h,m_w,m_radius);
			boxcar_sliding_window(m_wk,m_wk,m_box,m_h,m_w,m_radius);
		}
		else if(m_spatial_filter==QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER)
		{
			gaussian_recursive(jk,m_box,m_sigma_spatial*qx_min(m_h,m_w),0,m_h,m_w);
			gaussian_recursive(m_wk,m_box,m_sigma_spatial*qx_min(m_h,m_w),0,m_h,m_w);
		}
		for(y=0;y<m_h;y++)
		{
			for(x=0;x<m_w;x++)
			{
				jk[y][x]/=m_wk[y][x];
			}
		}
		//image_display(jk,m_h,m_w);
		if(i>0)
		{
			for(y=0;y<m_h;y++)
			{
				for(x=0;x<m_w;x++)
				{
					double kf;
					if(texture==NULL) kf=double(image[y][x]-image_min)/delta_scale;
					else kf=double(texture[y][x]-image_min)/delta_scale; /*cross/joint bilateral filtering*/
					int k=int(kf); 
					if(k==(i-1))
					{
						double alpha=(k+1)-kf;
						image_filtered[y][x]=(unsigned char)qx_min(qx_max(alpha*m_jk[jk_0][y][x]+(1.f-alpha)*m_jk[jk_1][y][x],0.f)+0.5f,255.f);
					}
					else if(k==i&&i==(m_nr_scale-1)) image_filtered[y][x]=(unsigned char)(m_jk[jk_1][y][x]+0.5f);
				}
			}
			jk_1=jk_0;
			jk_0=(jk_0+1)%2;
		}
	}
	//timer.time_display("bilateral filter");
	return(0);
}


int qx_constant_time_bilateral_filter::bilateral_filter(float **image_filtered,float**image,double sigma_range,unsigned char **texture)
{
	unsigned char image_min_u,image_max_u; 
	float image_min,image_max; 
	int y,x,jk_0,jk_1;
	if(sigma_range>QX_DEF_THRESHOLD_ZERO) 
	{
		m_sigma_range=sigma_range;
		color_weighted_table_update(m_table,m_sigma_range*QX_DEF_CTBF_INTENSITY_RANGE,QX_DEF_CTBF_INTENSITY_RANGE);
	}
	qx_timer timer;
	timer.start();
	if(texture==NULL)
	{
		vec_min_val(image_min,image[0],m_h*m_w);
		vec_max_val(image_max,image[0],m_h*m_w);
	}
	else
	{
		vec_min_val(image_min_u,texture[0],m_h*m_w); image_min=(float)image_min_u;
		vec_max_val(image_max_u,texture[0],m_h*m_w); image_max=(float)image_max_u;
	}
	m_nr_scale=qx_max(1,int(double(image_max-image_min)/(255*m_sigma_range)+0.5));
	//printf("[qx_max,qx_min]:[%5.5f,%5.5f]\n",image_max,image_min);
	//printf("[sigma_range: %1.3f]\n",m_sigma_range);
	//printf("[nr_scale: %d]\n",m_nr_scale);
	m_grayscale[0]=(double)image_min;
	m_grayscale[m_nr_scale-1]=(double)image_max;
	double delta_scale=double(image_max-image_min)/(m_nr_scale-1);
	for(int i=1;i<m_nr_scale-1;i++) m_grayscale[i]=(double)image_min+delta_scale*i;
	for(int i=0;i<m_nr_scale;i++)
	{
		double **jk;
		if(i==0)
		{
			jk_0=0;
			jk_1=1;
			jk=m_jk[jk_0];
		}
		else 
			jk=m_jk[jk_1];
		for(y=0;y<m_h;y++)
		{
			for(x=0;x<m_w;x++)
			{
				int index;
				if(texture==NULL) 
					index=int(abs(m_grayscale[i]-image[y][x])+0.5f);
				else 
					index=int(abs(m_grayscale[i]-texture[y][x])+0.5f); /*cross/joint bilateral filtering*/
				jk[y][x]=m_table[index]*image[y][x];
				m_wk[y][x]=m_table[index];
			}
		}
		if(m_spatial_filter==QX_DEF_CTBF_BOX_BILATERAL_FILTER)
		{
			boxcar_sliding_window(jk,jk,m_box,m_h,m_w,m_radius);
			boxcar_sliding_window(m_wk,m_wk,m_box,m_h,m_w,m_radius);
		}
		else if(m_spatial_filter==QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER)
		{
			gaussian_recursive(jk,m_box,m_sigma_spatial*qx_min(m_h,m_w),0,m_h,m_w);
			gaussian_recursive(m_wk,m_box,m_sigma_spatial*qx_min(m_h,m_w),0,m_h,m_w);
			//qx_filtering_guassian(jk,jk,m_box,m_sigma_spatial*qx_min(m_h,m_w),m_h,m_w); /*Exact gaussian spatial filtering*/
			//qx_filtering_guassian(m_wk,m_wk,m_box,m_sigma_spatial*qx_min(m_h,m_w),m_h,m_w); /*Exact gaussian spatial filtering*/
			//printf("[%d] ",i);
		}
		for(y=0;y<m_h;y++)
		{
			for(x=0;x<m_w;x++)
			{
				jk[y][x]/=m_wk[y][x];
			}
		}
		//image_display(jk,m_h,m_w);
		if(i>0)
		{
			for(y=0;y<m_h;y++)
			{
				for(x=0;x<m_w;x++)
				{
					double kf;
					if(texture==NULL)
						kf=double(image[y][x]-image_min)/delta_scale;
					else 
						kf=double(texture[y][x]-image_min)/delta_scale; /*cross/joint bilateral filtering*/
					int k=int(kf); 
					if(k==(i-1))
					{
						double alpha=(k+1)-kf; 
						image_filtered[y][x]=(float)(alpha*m_jk[jk_0][y][x]+(1.f-alpha)*m_jk[jk_1][y][x]);
					}
					else if(k==i&&i==(m_nr_scale-1)) 
						image_filtered[y][x]=(float)m_jk[jk_1][y][x];
				}
			}
			//image_display(image_filtered,m_h,m_w);
			jk_1=jk_0;
			jk_0=(jk_0+1)%2;
		}
	}
	//timer.time_display("bilateral filter");
	return(0);
}

 