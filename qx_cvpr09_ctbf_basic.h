/*******************************************************
\Author:	Qingxiong Yang (liiton)
\Function:	Basic functions.
********************************************************/
#ifndef QX_CVPR09_CTBF_BASIC_H
#define QX_CVPR09_CTBF_BASIC_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <vector>
#include <process.h>
#include <direct.h>
#include <io.h>
#include <time.h>
#include <string>
#include <memory.h>
#include <algorithm>
#include <functional>      // For greater<int>()
#include <iostream>
using namespace std;
#define QX_DEF_PADDING					10
#define QX_DEF_THRESHOLD_ZERO			1e-6
class   qx_timer		{public: void start();	float stop(); void time_display(char *disp=""); void fps_display(char *disp=""); private: clock_t m_begin; clock_t m_end;};

/*Box filter*/
void boxcar_sliding_window(double **out,double **in,double **temp,int h,int w,int radius);
void boxcar_sliding_window_x(double *out,double *in,int h,int w,int radius);
void boxcar_sliding_window_y(double *out,double *in,int h,int w,int radius);
/*Gaussian filter*/
int gaussian_recursive(double **image,double **temp,double sigma,int order,int h,int w);
void gaussian_recursive_x(double **od,double **id,int w,int h,double a0,double a1,double a2,double a3,double b1,double b2,double coefp,double coefn);
void gaussian_recursive_y(double **od,double **id,int w,int h,double a0,double a1,double a2,double a3,double b1,double b2,double coefp,double coefn);
/*basic functions*/

inline float qx_min(float a,float b){if(a<b) return(a); else return(b);}
inline float qx_max(float a,float b){if(a>b) return(a); else return(b);}
inline int qx_sum_u3(unsigned char *a) {return(a[0]+a[1]+a[2]);}
inline double qx_sum_d3(double*a){return(a[0]+a[1]+a[2]);}
inline unsigned char qx_min_u3(unsigned char *a){return(qx_min(qx_min(a[0],a[1]),a[2]));}
inline unsigned char qx_max_u3(unsigned char *a){return(qx_max(qx_max(a[0],a[1]),a[2]));}
inline void image_zero(float **in,int h,int w,float zero=0){memset(in[0],zero,sizeof(float)*h*w);}
inline void image_zero(double **in,int h,int w,double zero=0){memset(in[0],zero,sizeof(double)*h*w);}
inline void image_zero(unsigned char**in,int h,int w,unsigned char zero=0){memset(in[0],zero,sizeof(unsigned char)*h*w);}
inline void image_zero(double ***in,int h,int w,int d,double zero=0){memset(in[0][0],zero,sizeof(double)*h*w*d);}
inline unsigned char rgb_2_gray(unsigned char*in){return(unsigned char(0.299*in[0]+0.587*in[1]+0.114*in[2]+0.5));}
inline int qx_square_difference_u3(unsigned char *a,unsigned char *b){int d1,d2,d3; d1=(*a++)-(*b++); d2=(*a++)-(*b++);	d3=(*a++)-(*b++); return(int(d1*d1+d2*d2+d3*d3));}
void qx_specular_free_image(unsigned char ***image_specular_free,unsigned char ***image_normalized,float **diffuse_chromaticity_max,int h,int w);
inline double *get_color_weighted_table(double sigma_range,int len)
{
	double *table_color,*color_table_x; int y;
	table_color=new double [len];
	color_table_x=&table_color[0];
	for(y=0;y<len;y++) (*color_table_x++)=exp(-double(y*y)/(2*sigma_range*sigma_range));
	return(table_color);
}
inline void color_weighted_table_update(double *table_color,double dist_color,int len)
{
	double *color_table_x; int y;
	color_table_x=&table_color[0];
	for(y=0;y<len;y++) (*color_table_x++)=exp(-double(y*y)/(2*dist_color*dist_color));
}

inline void vec_min_val(float &min_val,float *in,int len)
{
	min_val=in[0];
	for(int i=1;i<len;i++) if(in[i]<min_val) min_val=in[i];	
}
inline void vec_min_val(unsigned char &min_val,unsigned char *in,int len)
{
	min_val=in[0];
	for(int i=1;i<len;i++) if(in[i]<min_val) min_val=in[i];	
}
inline void vec_max_val(float &max_val,float *in,int len)
{
	max_val=in[0];
	for(int i=1;i<len;i++) if(in[i]>max_val) max_val=in[i];	
}
inline void vec_max_val(unsigned char &max_val,unsigned char *in,int len)
{
	max_val=in[0];
	for(int i=1;i<len;i++) if(in[i]>max_val) max_val=in[i];	
}
/*memory*/
inline double *** qx_allocd_3(int n,int r,int c,int padding=QX_DEF_PADDING)
{
	double *a,**p,***pp;
    int rc=r*c;
    int i,j;
	a=(double*) malloc(sizeof(double)*(n*rc+padding));
	if(a==NULL) {printf("qx_allocd_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
    p=(double**) malloc(sizeof(double*)*n*r);
    pp=(double***) malloc(sizeof(double**)*n);
    for(i=0;i<n;i++) 
        for(j=0;j<r;j++) 
            p[i*r+j]=&a[i*rc+j*c];
    for(i=0;i<n;i++) 
        pp[i]=&p[i*r];
    return(pp);
}
inline void qx_freed_3(double ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline unsigned char** qx_allocu(int r,int c,int padding=QX_DEF_PADDING)
{
	unsigned char *a,**p;
	a=(unsigned char*) malloc(sizeof(unsigned char)*(r*c+padding));
	if(a==NULL) {printf("qx_allocu() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(unsigned char**) malloc(sizeof(unsigned char*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freeu(unsigned char **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline unsigned char *** qx_allocu_3(int n,int r,int c,int padding=QX_DEF_PADDING)
{
	unsigned char *a,**p,***pp;
    int rc=r*c;
    int i,j;
	a=(unsigned char*) malloc(sizeof(unsigned char )*(n*rc+padding));
	if(a==NULL) {printf("qx_allocu_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
    p=(unsigned char**) malloc(sizeof(unsigned char*)*n*r);
    pp=(unsigned char***) malloc(sizeof(unsigned char**)*n);
    for(i=0;i<n;i++) 
        for(j=0;j<r;j++) 
            p[i*r+j]=&a[i*rc+j*c];
    for(i=0;i<n;i++) 
        pp[i]=&p[i*r];
    return(pp);
}
inline void qx_freeu_3(unsigned char ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline float** qx_allocf(int r,int c,int padding=QX_DEF_PADDING)
{
	float *a,**p;
	a=(float*) malloc(sizeof(float)*(r*c+padding));
	if(a==NULL) {printf("qx_allocf() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(float**) malloc(sizeof(float*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freef(float **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline float *** qx_allocf_3(int n,int r,int c,int padding=QX_DEF_PADDING)
{
	float *a,**p,***pp;
    int rc=r*c;
    int i,j;
	a=(float*) malloc(sizeof(float)*(n*rc+padding));
	if(a==NULL) {printf("qx_allocf_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
    p=(float**) malloc(sizeof(float*)*n*r);
    pp=(float***) malloc(sizeof(float**)*n);
    for(i=0;i<n;i++) 
        for(j=0;j<r;j++) 
            p[i*r+j]=&a[i*rc+j*c];
    for(i=0;i<n;i++) 
        pp[i]=&p[i*r];
    return(pp);
}
inline void qx_freef_3(float ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline int** qx_alloci(int r,int c,int padding=QX_DEF_PADDING)
{
	int *a,**p;
	a=(int*) malloc(sizeof(int)*(r*c+padding));
	if(a==NULL) {printf("qx_alloci() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(int**) malloc(sizeof(int*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freei(int **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline double** qx_allocd(int r,int c,int padding=QX_DEF_PADDING)
{
	double *a,**p;
	a=(double*) malloc(sizeof(double)*(r*c+padding));
	if(a==NULL) {printf("qx_allocd() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(double**) malloc(sizeof(double*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freed(double **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
#endif