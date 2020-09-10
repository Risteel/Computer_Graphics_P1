using namespace std;
#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <algorithm>

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene(void)
{

	ui_instance = Qt_Opengl_Framework::getInstance();

}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage(QString filePath)
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);
	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));
	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath)
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB(void)
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (!img_data)
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0; j < img_width; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j * 4), rgb + (out_offset + j * 3));
		}
	}

	return rgb;
}


void Application::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0; i < 3; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = gray;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			img_data[offset_rgba + rr] = rgb[offset_rgb + rr] & 224;//取最高三位-->二進位11100000
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg] & 224;//取最高三位-->二進位11100000
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb] & 192;//取最高兩位-->二進位11000000
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();
	int *color = new int[32768](), *count = new int[32768]();
	int *colorTable = new int[256]();
	int index = 0;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int index = ((rgb[offset_rgb + rr] >> 3) << 10) + ((rgb[offset_rgb + gg] >> 3) << 5) + (rgb[offset_rgb + bb] >> 3);
			color[index]++;
			count[index]++;
		}
	}
	index = 0;
	sort(count, count + 32768);
	for (int i = 0; i < 32768 && index < 256; i++)
	{
		if (color[i] >= count[32512]) {
			colorTable[index++] = i;
		}
	}
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			int min = 196609, newBlue, newRed, newGreen;
			for (int k = 0; k < 256; k++) {
				int blue = ((colorTable[k] & 31) << 3), green = (((colorTable[k] >> 5) & 31) << 3), red = (((colorTable[k] >> 10) & 31) << 3);
				int distance = (rgb[offset_rgb + rr] - red)*(rgb[offset_rgb + rr] - red) + (rgb[offset_rgb + gg] - green)*(rgb[offset_rgb + gg] - green) + (rgb[offset_rgb + bb] - blue)*(rgb[offset_rgb + bb] - blue);
				if (distance < min) {
					min = distance;
					newBlue = blue;
					newRed = red;
					newGreen = green;
				}
			}
			img_data[offset_rgba + rr] = newRed;
			img_data[offset_rgba + gg] = newGreen;
			img_data[offset_rgba + bb] = newBlue;
		}
	}
	delete[] rgb;
	delete[] color;
	delete[] count;
	delete[] colorTable;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			double gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			int color = (gray > 127) * 255;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = color;
			}
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			double gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			int sign = rand() % 2, color;
			float  value = (rand() % 21) / 100.0;
			if (sign) {
				value += 0.5;
			}
			else {
				value = 0.5 - value;
			}
			color = (gray > (value * 255)) * 255;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = color;
			}
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();
	float *nrgb = new float[img_height*img_width * 3];
	for (int i = 0; i < img_height*img_width * 3; i++)
	{
		nrgb[i] = rgb[i];
	}
	for (int i = 0; i < img_height; i++)
	{
		for (int j = (i % 2 ? img_width - 1 : 0); j >= 0 && j < img_width; j += (i % 2 ? -1 : 1))
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			double oldpixel = 0.299 * nrgb[offset_rgb + rr] + 0.587 * nrgb[offset_rgb + gg] + 0.114 * nrgb[offset_rgb + bb], newpixel;
			double error;
			newpixel = ((oldpixel / 255.0) > 0.5) * 255;
			error = oldpixel - newpixel;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = newpixel;
				if (i % 2) {
					if (i + 1 < img_height) {
						auto down = (i + 1)*img_width * 3 + j * 3;
						nrgb[down + k] += (error*(5 / 16.0));
						if (j > 0) {
							nrgb[down - 3 + k] += (error*(1 / 16.0));
						}
						if (j + 1 < img_width) {
							nrgb[down + 3 + k] += (error*(3 / 16.0));
						}
					}
					if (j > 0) {
						nrgb[offset_rgb - 3 + k] += (error*(7 / 16.0));
					}
				}
				else {
					if (i + 1 < img_height) {
						auto down = (i + 1)*img_width * 3 + j * 3;
						nrgb[down + k] += (error*(5 / 16.0));
						if (j > 0) {
							nrgb[down - 3 + k] += (error*(3 / 16.0));
						}
						if (j + 1 < img_width) {
							nrgb[down + 3 + k] += (error*(1 / 16.0));
						}
					}
					if (j + 1 < img_width) {
						nrgb[offset_rgb + 3 + k] += (error*(7 / 16.0));
					}
				}
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	delete[] nrgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
	double *newRGB = new double[img_height*img_width](), *rgba = new double[img_height*img_width * 4];
	int counts = 0;
	double graySum;
	double grayAverage;
	graySum = 0;
	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			double gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			newRGB[i*img_width + j] = gray;
			rgba[offset_rgba] = gray;
			graySum += gray;
		}
	}
	sort(newRGB, newRGB + img_width*img_height);
	counts = newRGB[img_width*img_height - 1];
	grayAverage = graySum / ((float)(img_height*img_width));
	grayAverage /= 255.0;
	counts = (1 - grayAverage) *(img_height*img_width);
	newRGB[counts];
	graySum = 0;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4, color = 0;
			if (rgba[offset_rgba] >= newRGB[counts]) {
				color = 255;
			}
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = color;
			}
			graySum += color;
		}
	}
	graySum /= ((float)(img_height*img_width));
	delete[] rgb;
	delete[] newRGB;
	delete[] rgba;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();
	double mask[4][4] = { {0.7509,0.3529,0.5882,0.2353},{0.0588,0.9412,0.8235,0.4118},{0.4706,0.7647,0.8824,0.1176},{0.1765,0.5294,0.2941,0.6471} };
	Gray();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			for (int k = 0; k < 3; k++) {
				if (img_data[offset_rgba + k] > 255 * mask[i % 4][j % 4]) {
					img_data[offset_rgba + k] = 255;
					continue;
				}
				img_data[offset_rgba + k] = 0;
			}
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();
	double *nrgb = new double[img_height*img_width * 3];
	for (int i = 0; i < img_height*img_width * 3; i++)
	{
		nrgb[i] = rgb[i];
	}
	int rg[8] = { 0,36,73,109,146,182,219,255 }, b[4] = { 0,85,170,255 }, rgbCounts[3] = { 4,8,8 };
	for (int i = 0; i < img_height; i++)
	{
		for (int j = (i % 2 ? img_width - 1 : 0); j >= 0 && j < img_width; j += (i % 2 ? -1 : 1))
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			for (int k = 0; k < 3; k++) {
				double min = 256, error, pixel;
				for (int l = 0; l < rgbCounts[k]; l++) {
					if (k > 0) {
						if (abs(nrgb[offset_rgb + k] - rg[l]) < min) {
							pixel = rg[l];
							min = abs(nrgb[offset_rgb + k] - rg[l]);
						}
					}
					else {
						if (abs(nrgb[offset_rgb + k] - b[l]) < min) {
							pixel = b[l];
							min = abs(nrgb[offset_rgb + k] - rg[l]);
						}
					}
				}
				error = nrgb[offset_rgb + k] - pixel;
				img_data[offset_rgba + k] = pixel;
				if (i % 2) {
					if (i + 1 < img_height) {
						auto down = (i + 1)*img_width * 3 + j * 3;
						nrgb[down + k] += (error*(5 / 16.0));
						if (j > 0) {
							nrgb[down - 3 + k] += (error*(1 / 16.0));
						}
						if (j + 1 < img_width) {
							nrgb[down + 3 + k] += (error*(3 / 16.0));
						}
					}
					if (j > 0) {
						nrgb[offset_rgb - 3 + k] += (error*(7 / 16.0));
					}
				}
				else {
					if (i + 1 < img_height) {
						auto down = (i + 1)*img_width * 3 + j * 3;
						nrgb[down + k] += (error*(5 / 16.0));
						if (j > 0) {
							nrgb[down - 3 + k] += (error*(3 / 16.0));
						}
						if (j + 1 < img_width) {
							nrgb[down + 3 + k] += (error*(1 / 16.0));
						}
					}
					if (j + 1 < img_width) {
						nrgb[offset_rgb + 3 + k] += (error*(7 / 16.0));
					}
				}
			}
		}
	}
	delete[] rgb;
	delete[] nrgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering(double filter[][5])
{
	unsigned char *rgb = this->To_RGB();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4, x, y;
			for (int k = 0; k < 3; k++) {
				double pixel = 0;
				for (int m = -2; m < 3; m++) {
					for (int n = -2; n < 3; n++) {
						y = i + m; x = j + n;
						if (x < 0 || x >= img_width) {
							x = j - n;
						}
						if (y < 0 || y >= img_height) {
							y = i - m;
						}
						pixel += (double)rgb[y*img_width * 3 + x * 3 + k] * filter[m + 2][n + 2];
					}
				}
				if (pixel < 0) {
					pixel = 0;
				}
				else if (pixel > 255) {
					pixel = 255;
				}
				img_data[offset_rgba + k] = pixel;
			}
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::filtering(double *filter, int N)
{
	unsigned char *rgb = this->To_RGB();
	unsigned int sum = 0;
	sum = (1 << (N - 1 << 1));
	int start = -(N >> 1), end = N - (N >> 1);
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4, x, y;
			for (int k = 0; k < 3; k++) {
				double pixel = 0;
				for (int m = start; m < end; m++) {
					for (int n = start; n < end; n++) {
						y = i + m; x = j + n;
						if (x < 0 || x >= img_width) {
							x = j - n;
						}
						if (y < 0 || y >= img_height) {
							y = i - m;
						}
						pixel += rgb[y*img_width * 3 + x * 3 + k] * filter[n - start] * filter[m - start];
					}
				}
				img_data[offset_rgba + k] = (pixel / sum);
				if (pixel / sum > 0) {
					sum = sum;
				}
			}
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	double filter[5][5] = { {1,1,1,1,1},{ 1,1,1,1,1 },{ 1,1,1,1,1 },{ 1,1,1,1,1 },{ 1,1,1,1,1 } };
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			filter[i][j] /= 25.0;
		}
	}
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	double filter[5][5] = { {1,2,3,2,1},{ 2,4,6,4,2 },{ 3,6,9,6,3 },{ 2,4,6,4,2 },{ 1,2,3,2,1 } };
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			filter[i][j] /= 81.0;
		}
	}
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	double filter[5][5] = { { 1,4,6,4,1},{ 4,16,24,16,4 },{ 6,24,36,24,6 },{ 4,16,24,16,4 },{ 1,4,6,4,1 } };
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			filter[i][j] /= 256.0;
		}
	}
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N(unsigned int N)
{
	if (N <= 1) return;
	double **bufFilter = new double*[N - 1], *filter = new double[N]();
	bufFilter[0] = new double[N]();
	bufFilter[0][0] = bufFilter[0][1] = 1;
	for (int i = 1; i < N - 1; i++) {
		bufFilter[i] = new double[N]();
		bufFilter[i][0] = bufFilter[i][i + 1] = 1;
		for (int j = 1; j <= i; j++) {
			bufFilter[i][j] = bufFilter[i - 1][j - 1] + bufFilter[i - 1][j];
		}
	}
	for (int i = 0; i < N; i++) {
		filter[i] = bufFilter[N - 2][i];
	}
	filtering(filter, N);
	for (int i = 0; i < N - 1; i++) {
		delete[] bufFilter[i];
	}
	delete[] bufFilter;
	delete[] filter;
}
//------------
void Application::Filter_Gaussian_N(unsigned int N, unsigned char *img)
{
	if (N <= 1) return;
	double **bufFilter = new double*[N - 1], *filter = new double[N]();
	unsigned char *rgb = this->To_RGB();
	unsigned int sum = 0;
	bufFilter[0] = new double[N]();
	bufFilter[0][0] = bufFilter[0][1] = 1;
	for (int i = 1; i < N - 1; i++) {
		bufFilter[i] = new double[N]();
		bufFilter[i][0] = bufFilter[i][i + 1] = 1;
		for (int j = 1; j <= i; j++) {
			bufFilter[i][j] = bufFilter[i - 1][j - 1] + bufFilter[i - 1][j];
		}
	}
	for (int i = 0; i < N; i++) {
		filter[i] = bufFilter[N - 2][i];
	}
	sum = (1 << (N - 1 << 1));
	int start = -(N >> 1), end = N - (N >> 1);
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4, x, y;
			for (int k = 0; k < 3; k++) {
				double pixel = 0;
				for (int m = start; m < end; m++) {
					for (int n = start; n < end; n++) {
						y = i + m; x = j + n;
						if (x < 0 || x >= img_width) {
							x = j - n;
						}
						if (y < 0 || y >= img_height) {
							y = i - m;
						}
						pixel += rgb[y*img_width * 3 + x * 3 + k] * filter[n - start] * filter[m - start];
					}
				}
				img[offset_rgba + k] = (pixel / sum);
			}
		}
	}
	for (int i = 0; i < N - 1; i++) {
		delete[] bufFilter[i];
	}
	delete[] rgb;
	delete[] bufFilter;
	delete[] filter;
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	double filter[5][5] = { { 1,4,6,4,1 },{ 4,16,24,16,4 },{ 6,24,-220,24,6 },{ 4,16,24,16,4 },{ 1,4,6,4,1 } };
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			filter[i][j] /= -256.0;
		}
	}
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();
	Filter_Edge();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_rgb = i*img_width * 3 + j * 3;
			for (int k = 0; k < 3; k++) {
				int pixel = img_data[offset_rgba + k] + rgb[offset_rgb + k];
				if (pixel > 255) {
					pixel = 255;
				}
				img_data[offset_rgba + k] = pixel;
			}
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	Resize(0.5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	Resize(2);
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize(float scale)
{
	int reImg_width = img_width * scale, reImg_height = img_height * scale;
	if (scale > 1 || scale < 1) {
		unsigned char *reImg_data = new unsigned char[reImg_width* reImg_height * 4];
		if (scale > 1) {
			double filter[4][4][4] = { { { 1 / 16.0,1 / 8.0,1 / 16.0 ,0 },{ 1 / 8.0,1 / 4.0,1 / 8.0,0 },{ 1 / 16.0,1 / 8.0,1 / 16.0,0 },{ 0,0,0,0 } },
			{ { 1.0 / 64,3.0 / 64,3.0 / 64,1.0 / 64 },{ 3.0 / 64,9.0 / 64,9.0 / 64,3.0 / 64 },{ 3.0 / 64,9.0 / 64,9.0 / 64,3.0 / 64 },{ 1.0 / 64,3.0 / 64,3.0 / 64,1.0 / 64 } } ,
			{ { 1.0 / 32,2.0 / 32,1.0 / 32,0 },{ 3.0 / 32,6.0 / 32,3.0 / 32,0 },{ 3.0 / 32,6.0 / 32,3.0 / 32,0 },{ 1.0 / 32,2.0 / 32,1.0 / 32,0 } },
			{ { 1.0 / 32,3.0 / 32,3.0 / 32,1.0 / 32 },{ 2.0 / 32,6.0 / 32,6.0 / 32,2.0 / 32 },{ 1.0 / 32,3.0 / 32,3.0 / 32,1.0 / 32 },{ 0,0,0,0 } } };
			for (int i = 0; i < reImg_height; i++) {
				for (int j = 0; j < reImg_width; j++) {
					int y = i / scale, x = j / scale, oddI = i % 2, oddJ = j % 2, upX, upY, index;
					double pixel = 0;
					if (oddI == 0 && oddJ == 0) {
						upX = 2; upY = 2;
						index = 0;
					}
					if (oddI == 1 && oddJ == 1) {
						upX = 3; upY = 3;
						index = 1;
					}
					if (oddI == 0 && oddJ == 1) {
						upX = 3; upY = 2;
						index = 3;
					}
					if (oddI == 1 && oddJ == 0) {
						upX = 2; upY = 3;
						index = 2;
					}
					for (int k = 0; k < 4; k++) {
						pixel = 0;
						for (int m = -1; m < upY; m++) {
							for (int n = -1; n < upX; n++) {
								int nowX = x + n, nowY = y + m;
								if (nowX >= img_width || nowX < 0) {
									nowX = x - n;
								}
								if (nowY >= img_height || nowY < 0) {
									nowY = y - m;
								}
								pixel += img_data[nowY*img_width * 4 + nowX * 4 + k] * filter[index][m + 1][n + 1];
							}
						}
						reImg_data[i* reImg_width * 4 + j * 4 + k] = pixel;
					}
				}
			}
		}
		else {
			float filter[3][3] = { { 1 / 16.0,1 / 8.0,1 / 16.0 },{ 1 / 8.0,1 / 4.0,1 / 8.0 },{ 1 / 16.0,1 / 8.0,1 / 16.0 } };
			for (int i = 0; i < reImg_height; i++) {
				for (int j = 0; j < reImg_width; j++) {
					int y = i / scale, x = j / scale;
					float pixel;
					for (int k = 0; k < 4; k++) {
						pixel = 0;
						for (int m = -1; m < 2; m++) {
							for (int n = -1; n < 2; n++) {
								int nowX = x + n, nowY = y + m;
								if (nowX >= img_width || nowX < 0) {
									nowX = x - n;
								}
								if (nowY >= img_height || nowY < 0) {
									nowY = y - m;
								}
								pixel += img_data[nowY*img_width * 4 + nowX * 4 + k] * filter[m + 1][n + 1];
							}
						}
						reImg_data[i* reImg_width * 4 + j * 4 + k] = pixel;
					}
				}
			}
		}
		img_data = new unsigned char[reImg_width*reImg_height * 4];
		for (int i = 0; i < reImg_width*reImg_height * 4; i++) {
			img_data[i] = reImg_data[i];
		}
		delete[] reImg_data;
		mImageDst = QImage(img_data, img_width * scale, img_height * scale, QImage::Format_ARGB32);
		img_width *= scale;
		img_height *= scale;
		renew();
	}

}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate(float angleDegrees)
{
	unsigned char *rgb = this->To_RGB();
	double *reImg_data = new double[img_height*img_width * 4];
	double PI = acos(-1), rad = -(angleDegrees*PI / 180.0);
	double matrix[2][2] = { {cos(rad),-sin(rad)},{sin(rad),cos(rad)} };
	int y = img_height / 2.0, x = img_width / 2.0;
	float offset_x = x - ((float)x*cos(rad) - (float)y*sin(rad)), offset_y = y - ((float)y*cos(rad) + (float)x*sin(rad));
	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			float x = (float)j*cos(rad) - (float)i*sin(rad), y = (float)i*cos(rad) + (float)j*sin(rad);
			x += offset_x;
			y += offset_y;
			if (y >= 0 && y < img_height && x >= 0 && x < img_width) {
				reImg_data[i*img_width * 4 + j * 4 + aa] = WHITE;
				int left = (int)x, right = (int)(x + 1), down = (int)(y + 1), up = (int)y;
				reImg_data[i* img_width * 4 + j * 4 + aa] = WHITE;
				if ((x == (int)x && y == (int)y) || down >= img_height || right >= img_width || left < 0 || up < 0) {
					for (int k = 0; k < 3; k++) {
						reImg_data[i* img_width * 4 + j * 4 + k] = img_data[(int)y*img_width * 4 + (int)x * 4 + k];
					}
				}
				else {
					double pixel = 0, r1, r2;
					for (int k = 0; k < 3; k++) {
						pixel = 0;
						r1 = (float)(right - x)*img_data[down*img_width * 4 + left * 4 + k] + (float)(x - left)*img_data[down*img_width * 4 + right * 4 + k];
						r2 = (float)(right - x)*img_data[up*img_width * 4 + left * 4 + k] + (float)(x - left) *img_data[up*img_width * 4 + right * 4 + k];
						pixel = (float)(down - y)*r2 + (float)(y - up)*r2;
						reImg_data[i* img_width * 4 + j * 4 + k] = pixel;
					}
					reImg_data[i* img_width * 4 + j * 4 + aa] = WHITE;
				}
			}
			else {
				reImg_data[i*img_width * 4 + j * 4 + rr] = reImg_data[i*img_width * 4 + j * 4 + gg] = reImg_data[i*img_width * 4 + j * 4 + bb] = 0;
				reImg_data[i*img_width * 4 + j * 4 + aa] = 255;
			}
		}
	}
	for (int i = 0; i < img_height* img_width * 4; i++) {
		img_data[i] = reImg_data[i];
	}
	delete[] reImg_data;
	delete[] rgb;
	this->mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge(QString filePath)
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image(int tMethod)
{
	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgba = i*img_width * 4 + j * 4;
			double AF, AB;
			switch (tMethod) {
			case over:AF = img_data[offset_rgba + aa] / 255.0, AB = (1 - AF)*(img_data2[offset_rgba + aa] / 255.0); break;
			case in:AF = (img_data[offset_rgba + aa] / 255.0)*(img_data2[offset_rgba + aa] / 255.0), AB = 0; break;
			case out:AF = (img_data[offset_rgba + aa] / 255.0)*(1 - (img_data2[offset_rgba + aa] / 255.0)), AB = 0; break;
			case atop:AF = (img_data[offset_rgba + aa] / 255.0)*(img_data2[offset_rgba + aa] / 255.0), AB = (1 - (img_data[offset_rgba + aa] / 255.0))*(img_data2[offset_rgba + aa] / 255.0); break;
			case xor:AF = (img_data[offset_rgba + aa] / 255.0)*(1 - (img_data2[offset_rgba + aa] / 255.0)), AB = (1 - (img_data[offset_rgba + aa] / 255.0))*(img_data2[offset_rgba + aa] / 255.0); break;
			}
			for (int k = 0; k < 3; k++) {
				double data = (double)(AF*img_data[offset_rgba + k]) + (double)(AB*img_data2[offset_rgba + k]);
				if (data > 255) {
					data = 255;
				}
				img_data[offset_rgba + k] = data;
			}
			img_data[offset_rgba + aa] = (AF + AB) * 255;
		}
	}
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{
		Comp_image(over);
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{
		Comp_image(in);
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{
		Comp_image(out);
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{
		Comp_image(atop);
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{
		Comp_image(xor);
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	int brushSize[3] = { 7,3,1 };
	unsigned char *sourceImage = new unsigned char[img_width*img_height * 4], *referenceImage = new unsigned char[img_width*img_height * 4];
	for (int i = 0; i < img_width*img_height * 4; i++) {
		sourceImage[i] = img_data[i];
		referenceImage[i] = img_data[i];
	}
	for (int i = 0; i < 3; i++) {
		Filter_Gaussian_N(brushSize[i]+1 , referenceImage);
		NPR_Paint_Layer(img_data, referenceImage, brushSize[i]);
		for (int j = 0; j < img_width*img_height * 4; j++) {
			referenceImage[j] = sourceImage[j];
		}
	}
	delete[] sourceImage;
	delete[] referenceImage;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::NPR_Paint_Layer(unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize)
{
	double *difference = new double[img_width*img_height]();
	vector<Stroke> S;
	for (int i = 0; i < img_height; i++) {
		for (int j = 0; j < img_width; j++) {
			int offset_rgba = i*img_width * 4 + j * 4;
			double dis = 0;
			for (int k = 0; k < 3; k++) {
				dis += (tCanvas[offset_rgba + k] - tReferenceImage[offset_rgba + k])*(tCanvas[offset_rgba + k] - tReferenceImage[offset_rgba + k]);
			}
			difference[i*img_width + j] = sqrt(dis);
		}
	}
	for (int i = 0; i < img_height; i += tBrushSize)
	{
		for (int j = 0; j < img_width; j += tBrushSize)
		{
			int n = tBrushSize / 2;
			int startX = j - n / 2, startY = i - n / 2, endX = j + tBrushSize - n / 2, endY = i + tBrushSize - n / 2;
			double sum = 0;
			if (startX >= 0 && endX < img_width && startY >= 0 && endY < img_height) {
				int setX, setY, max = -1, offset_rgba;
				for (int y = startY; y <= endY; y++) {
					for (int x = startX; x <= endX; x++) {
						sum += difference[y*img_width + x];
						if (difference[y*img_width + x] > max) {
							max = difference[y*img_width + x];
							setX = x;
							setY = y;
						}
					}
				}
				sum /= (float)(tBrushSize*tBrushSize);
				double T = 1.1 - (tBrushSize / 10.0);
				if (sum > T*difference[setY*img_width + setX]) {
					offset_rgba = setY*img_width * 4 + setX * 4;
					S.push_back(Stroke(tBrushSize, setX + rand() % 3 - 1, setY + rand() % 3 - 1, tReferenceImage[offset_rgba + rr], tReferenceImage[offset_rgba + gg], tReferenceImage[offset_rgba + bb], tReferenceImage[offset_rgba + aa]));
				}
			}
		}
	}
	for (int i = 0; i < S.size(); i++) {
		int s1 = rand() % S.size(), s2 = rand() % S.size();
		Stroke temp = S[s1];
		S[s1] = S[s2];
		S[s2] = temp;
	}
	for (int i = 0; i < S.size(); i++) {
		Paint_Stroke(S[i]);
	}
	S.clear();
	delete[] difference;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke(const Stroke& s)
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++)
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++)
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height))
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared)
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				}
				else if (dist_squared == radius_squared + 1)
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
	radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}



