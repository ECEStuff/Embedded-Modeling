//////////////////////////////////////////////////////////////////////
// canny.cpp: prepared for shortened A7 (02/10/20, RD)
//////////////////////////////////////////////////////////////////////

#include <specc.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned char img[4110080];
typedef short int simg[4110080];

class i_img_sender;
class i_img_receiver;
class i_img_tranceiver;
class c_img_queue;
class Stimulus;
class Monitor;
class DataIn;
class DataOut;
class Receive_Image;
class Gaussian_Kernel;
class BlurX;
class BlurY;
class Gaussian_Smooth;
class Derivative_X_Y;
class Magnitude_X_Y;
class Non_Max_Supp;
class Apply_Hysteresis;
class DUT;
class Platform;
class Main;

class i_img_sender
{
private:
public:
    virtual ~i_img_sender(void){};
    virtual void send(unsigned char [4110080]) = 0;
private:
};

class i_img_receiver
{
private:
public:
    virtual ~i_img_receiver(void){};
    virtual void receive(unsigned char (*)[4110080]) = 0;
private:
};

class i_img_tranceiver
{
private:
public:
    virtual ~i_img_tranceiver(void){};
    virtual void receive(unsigned char (*)[4110080]) = 0;
    virtual void send(unsigned char [4110080]) = 0;
private:
};

class c_img_queue : public _specc::channel, public i_img_sender, public i_img_receiver, public i_img_tranceiver
{
private:
    const unsigned long int (&size);
public:
    c_img_queue(const unsigned long int (&size));
    virtual ~c_img_queue(void);
    void cleanup(void);
    void receive(unsigned char (*)[4110080]);
    void send(unsigned char [4110080]);
    void setup(void);
private:
    unsigned char (*buffer)[4110080];
    unsigned long int n;
    unsigned long int p;
    _specc::event r;
    _specc::event s;
    unsigned long int wr;
    unsigned long int ws;
};

class Stimulus : public _specc::behavior
{
private:
    i_img_sender (&ImgOut);
public:
    Stimulus(unsigned int _idcnt, i_img_sender (&ImgOut));
    virtual ~Stimulus(void);
    void main(void);
    int read_pgm_image(const char *, unsigned char *, int, int);
private:
    unsigned char Image[4110080];
};

class Monitor : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
public:
    Monitor(unsigned int _idcnt, i_img_receiver (&ImgIn));
    virtual ~Monitor(void);
    void main(void);
    int write_pgm_image(const char *, unsigned char *, int, int, const char *, int);
private:
    unsigned char EdgeImage[4110080];
};

class DataIn : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
    i_img_sender (&ImgOut);
public:
    DataIn(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut));
    virtual ~DataIn(void);
    void main();
private:
    unsigned char Image[4110080];
};

class DataOut : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
    i_img_sender (&ImgOut);
public:
    DataOut(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut));
    virtual ~DataOut(void);
    void main();
private:
    unsigned char Image[4110080];
};

class Receive_Image : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
    unsigned char (&image)[4110080];
public:
    Receive_Image(unsigned int _idcnt, i_img_receiver (&ImgIn), unsigned char (&image)[4110080]);
    virtual ~Receive_Image(void);
    void main(void);
private:
};

class Gaussian_Kernel : public _specc::behavior
{
private:
    float (&gaussian_kernel)[21];
    int (&kernel_center);
public:
    Gaussian_Kernel(unsigned int _idcnt, float (&gaussian_kernel)[21], int (&kernel_center));
    virtual ~Gaussian_Kernel(void);
    void main(void);
    void make_gaussian_kernel(float, float *, int *);
private:
};

class BlurX : public _specc::behavior
{
private:
    unsigned char (&image)[4110080];
    float (&kernel)[21];
    int (&center);
    float (&tempim)[4110080];
public:
    BlurX(unsigned int _idcnt, unsigned char (&image)[4110080], float (&kernel)[21], int (&center), float (&tempim)[4110080]);
    virtual ~BlurX(void);
    void blur_x(int, int);
    void main(void);
private:
};

class BlurY : public _specc::behavior
{
private:
    float (&tempim)[4110080];
    float (&kernel)[21];
    int (&center);
    short int (&smoothedim)[4110080];
public:
    BlurY(unsigned int _idcnt, float (&tempim)[4110080], float (&kernel)[21], int (&center), short int (&smoothedim)[4110080]);
    virtual ~BlurY(void);
    void blur_y(int, int);
    void main(void);
private:
};

class Gaussian_Smooth : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
    short int (&smoothedim)[4110080];
public:
    Gaussian_Smooth(unsigned int _idcnt, i_img_receiver (&ImgIn), short int (&smoothedim)[4110080]);
    virtual ~Gaussian_Smooth(void);
    void main(void);
private:
    int center;
    unsigned char image[4110080];
    float kernel[21];
    float tempim[4110080];
    BlurX blurX;
    BlurY blurY;
    Gaussian_Kernel gauss;
    Receive_Image receive;
};

class Derivative_X_Y : public _specc::behavior
{
private:
    short int (&smoothedim)[4110080];
    short int (&delta_x)[4110080];
    short int (&delta_y)[4110080];
public:
    Derivative_X_Y(unsigned int _idcnt, short int (&smoothedim)[4110080], short int (&delta_x)[4110080], short int (&delta_y)[4110080]);
    virtual ~Derivative_X_Y(void);
    void derivative_x_y(int, int);
    void main(void);
private:
};

class Magnitude_X_Y : public _specc::behavior
{
private:
    short int (&delta_x)[4110080];
    short int (&delta_y)[4110080];
    short int (&magnitude)[4110080];
public:
    Magnitude_X_Y(unsigned int _idcnt, short int (&delta_x)[4110080], short int (&delta_y)[4110080], short int (&magnitude)[4110080]);
    virtual ~Magnitude_X_Y(void);
    void magnitude_x_y(int, int);
    void main(void);
private:
};

class Non_Max_Supp : public _specc::behavior
{
private:
    short int (&gradx)[4110080];
    short int (&grady)[4110080];
    short int (&mag)[4110080];
    unsigned char (&nms)[4110080];
public:
    Non_Max_Supp(unsigned int _idcnt, short int (&gradx)[4110080], short int (&grady)[4110080], short int (&mag)[4110080], unsigned char (&nms)[4110080]);
    virtual ~Non_Max_Supp(void);
    void main(void);
    void non_max_supp(int, int, unsigned char *);
private:
};

class Apply_Hysteresis : public _specc::behavior
{
private:
    short int (&mag)[4110080];
    unsigned char (&nms)[4110080];
    i_img_sender (&ImgOut);
public:
    Apply_Hysteresis(unsigned int _idcnt, short int (&mag)[4110080], unsigned char (&nms)[4110080], i_img_sender (&ImgOut));
    virtual ~Apply_Hysteresis(void);
    void apply_hysteresis(int, int, float, float, unsigned char *);
    void follow_edges(unsigned char *, short int *, short int, int);
    void main(void);
private:
    unsigned char EdgeImage[4110080];
};

class DUT : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
    i_img_sender (&ImgOut);
public:
    DUT(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut));
    virtual ~DUT(void);
    void main(void);
private:
    short int delta_x[4110080];
    short int delta_y[4110080];
    short int magnitude[4110080];
    unsigned char nms[4110080];
    short int smoothedim[4110080];
    Apply_Hysteresis apply_hysteresis;
    Derivative_X_Y derivative_x_y;
    Gaussian_Smooth gaussian_smooth;
    Magnitude_X_Y magnitude_x_y;
    Non_Max_Supp non_max_supp;
};

class Platform : public _specc::behavior
{
private:
    i_img_receiver (&ImgIn);
    i_img_sender (&ImgOut);
public:
    Platform(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut));
    virtual ~Platform(void);
    void main();
private:
    unsigned long int _scc_const_port_0;
    unsigned long int _scc_const_port_1;
    DUT canny;
    DataIn din;
    DataOut dout;
    c_img_queue q1;
    c_img_queue q2;
};

class Main : public _specc::class_type
{
private:
public:
    Main(unsigned int _idcnt);
    virtual ~Main(void);
    int main(void);
private:
    unsigned long int _scc_const_port_0;
    unsigned long int _scc_const_port_1;
    Monitor monitor;
    Platform platform;
    c_img_queue q1;
    c_img_queue q2;
    Stimulus stimulus;
};
void _scc_bit4_err_handle(const _bit4&);

//////////////////////////////////////////////////////////////////////

unsigned int _IDcnt = 0;

c_img_queue::c_img_queue(const unsigned long int (&size))
    : _specc::channel(), size(size),
    buffer(0),
    n(0ul),
    p(0ul),
    wr(0ul),
    ws(0ul)
{   
}

c_img_queue::~c_img_queue(void)
{   
}

void c_img_queue::cleanup(void)
{   
    if ( !n)
    {   
	free(buffer);
	buffer = 0;
    }
}

void c_img_queue::receive(unsigned char (*d)[4110080])
{   
    while( !n)
    {   
	wr++ ;
	_specc::wait(event(&r), ((void*)0));
	wr-- ;
    }
    if (n <= p)
    {   
	{ unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) ( *d)[_scc_index_0] = (buffer[p - 
	    n])[_scc_index_0]; }
    }
    else 
    {   
	{ unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) ( *d)[_scc_index_0] = (buffer[p + 
	    size - n])[_scc_index_0]; }
    }
    n-- ;
    if (ws)
    {   
	_specc::notify(event(&s), ((void*)0));
    }
    cleanup();
}

void c_img_queue::send(unsigned char d[4110080])
{   
    while(n >= size)
    {   
	ws++ ;
	_specc::wait(event(&s), ((void*)0));
	ws-- ;
    }
    setup();
    { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) (buffer[p])[_scc_index_0] = (d)[_scc_index_0]; }
    p++ ;
    if (p >= size)
    {   
	p = 0;
    }
    n++ ;
    if (wr)
    {   
	_specc::notify(event(&r), ((void*)0));
    }
}

void c_img_queue::setup(void)
{   
    if ( !buffer)
    {   
	unsigned char dummy[4110080];
	unsigned long int i;

	if ( !(buffer = (unsigned char (*)[4110080])malloc(sizeof(unsigned char [4110080]) * 
		    size)))
	{   
	    perror("c_typed_queue");
	    abort();
	}
	for(i = 0; i < size; i++ )
	{   
	    memcpy( &buffer[i],  &dummy, sizeof(unsigned char [4110080]));
	}
    }
}

Stimulus::Stimulus(unsigned int _idcnt, i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgOut(ImgOut)
{   
}

Stimulus::~Stimulus(void)
{   
}

void Stimulus::main(void)
{   
    int i = 0;
    int n = 0;
    char infilename[40];

    for(i = 0; i < 20; i++ )
    {   
	n = i % 20;
	sprintf(infilename, "video/EngPlaza%03d.pgm", n + 1);
	if (0)
	    printf("Reading the image %s.\n", infilename);
	if (read_pgm_image(infilename, Image, 1520, 2704) == 0)
	{   
	    fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	    exit(1);
	}
	ImgOut.send(Image);
	printf("Stimulus sent frame%3d.\n", n + 1);
    }
}

int Stimulus::read_pgm_image(const char *infilename, unsigned char *image, 
    int rows, int cols)
{   
    struct _IO_FILE *fp;
    char buf[71];
    int c;
    int r;

    if (infilename == ((void *)0))
	fp = stdin;
    else 
    {   
	if ((fp = fopen(infilename, "r")) == ((void *)0))
	{   
	    fprintf(stderr, "Error reading the file %s in read_pgm_image().\n", 
		infilename);
	    return (0);
	}
    }
    fgets(buf, 70, fp);
    if (strncmp(buf, "P5", 2) != 0)
    {   
	fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	fprintf(stderr, "read_pgm_image().\n");
	if (fp != stdin)
	    fclose(fp);
	return (0);
    }
    do 
    {   
	fgets(buf, 70, fp);
    }
    while(buf[0] == '#');
    __isoc99_sscanf(buf, "%d %d",  &c,  &r);
    if (c != cols || r != rows)
    {   
	fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, 
	    cols, rows);
	fprintf(stderr, "read_pgm_image().\n");
	if (fp != stdin)
	    fclose(fp);
	return (0);
    }
    do 
    {   
	fgets(buf, 70, fp);
    }
    while(buf[0] == '#');
    if ((unsigned int)rows != fread(image, cols, rows, fp))
    {   
	fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	if (fp != stdin)
	    fclose(fp);
	return (0);
    }
    if (fp != stdin)
	fclose(fp);
    return (1);
}

Monitor::Monitor(unsigned int _idcnt, i_img_receiver (&ImgIn))
    : _specc::behavior(_idcnt), ImgIn(ImgIn)
{   
}

Monitor::~Monitor(void)
{   
}

void Monitor::main(void)
{   
    char outfilename[128];
    int i;
    int n;

    for(i = 0; i < 20; i++ )
    {   
	ImgIn.receive( &EdgeImage);
	n = i % 20;
	printf("Monitor received frame%3d.\n", n + 1);
	sprintf(outfilename, "EngPlaza%03d_edges.pgm", n + 1);
	if (0)
	    printf("Writing the edge image in the file %s.\n", outfilename);
	if (write_pgm_image(outfilename, EdgeImage, 1520, 2704, "", 255) == 
	    0)
	{   
	    fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	    exit(1);
	}
    }
    if (0)
	printf("Monitor exits simulation.\n");
    sim_exit(0);
}

int Monitor::write_pgm_image(const char *outfilename, unsigned char *image, 
    int rows, int cols, const char *comment, int maxval)
{   
    struct _IO_FILE *fp;

    if (outfilename == ((void *)0))
	fp = stdout;
    else 
    {   
	if ((fp = fopen(outfilename, "w")) == ((void *)0))
	{   
	    fprintf(stderr, "Error writing the file %s in write_pgm_image().\n", 
		outfilename);
	    return (0);
	}
    }
    fprintf(fp, "P5\n%d %d\n", cols, rows);
    if (comment != ((void *)0))
	if (strlen(comment) <= 70)
	    fprintf(fp, "# %s\n", comment);
    fprintf(fp, "%d\n", maxval);
    if ((unsigned int)rows != fwrite(image, cols, rows, fp))
    {   
	fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	if (fp != stdout)
	    fclose(fp);
	return (0);
    }
    if (fp != stdout)
	fclose(fp);
    return (1);
}

DataIn::DataIn(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut)
{   
}

DataIn::~DataIn(void)
{   
}

void DataIn::main()
{   
    while(1)
    {   
	ImgIn.receive( &Image);
	ImgOut.send(Image);
    }
}

DataOut::DataOut(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut)
{   
}

DataOut::~DataOut(void)
{   
}

void DataOut::main()
{   
    while(1)
    {   
	ImgIn.receive( &Image);
	ImgOut.send(Image);
    }
}

Receive_Image::Receive_Image(unsigned int _idcnt, i_img_receiver (&ImgIn), unsigned char (&image)[4110080])
    : _specc::behavior(_idcnt), ImgIn(ImgIn), image(image)
{   
}

Receive_Image::~Receive_Image(void)
{   
}

void Receive_Image::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    unsigned char Image[4110080];
    Tstart = clock();

    ImgIn.receive( &Image);
    { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) (image)[_scc_index_0] = (Image)[_scc_index_0]; }
    
    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Receive_Image time: %f\n", T1);
}

Gaussian_Kernel::Gaussian_Kernel(unsigned int _idcnt, float (&gaussian_kernel)[21], int (&kernel_center))
    : _specc::behavior(_idcnt), gaussian_kernel(gaussian_kernel), kernel_center(kernel_center)
{   
}

Gaussian_Kernel::~Gaussian_Kernel(void)
{   
}

void Gaussian_Kernel::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    int windowsize;
    int center;
    float kernel[21] = 
    { 0.000000e+00f };
    Tstart = clock();

    if (0)
	printf("   Computing the gaussian smoothing kernel.\n");
    make_gaussian_kernel(6.000000000000000e-01, kernel,  &windowsize);
    center = windowsize / 2;
    { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<21;_scc_index_0++) (gaussian_kernel)[_scc_index_0] = (kernel)[_scc_index_0]; }
    kernel_center = center;
    
    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Gaussian_Kernel time: %f\n", T1);
}

void Gaussian_Kernel::make_gaussian_kernel(float sigma, float *kernel, 
    int *windowsize)
{   
    int center;
    int i;
    float fx;
    float sum = 0.000000e+00f;
    float x;

     *windowsize = 1 + 2 * ceil(2.500000000000000e+00 * sigma);
    center = ( *windowsize) / 2;
    if (0)
	printf("      The kernel has %d elements.\n",  *windowsize);
    for(i = 0; i < ( *windowsize); i++ )
    {   
	x = (float)(i - center);
	fx = pow(2.718280000000000e+00,  -5.000000000000000e-01 * x * 
	    x / (sigma * sigma)) / (sigma * sqrt(6.283185300000000e+00));
	kernel[i] = fx;
	sum += fx;
    }
    for(i = 0; i < ( *windowsize); i++ )
	kernel[i] /= sum;
    if (0)
    {   
	printf("The filter coefficients are:\n");
	for(i = 0; i < ( *windowsize); i++ )
	    printf("kernel[%d] = %f\n", i, kernel[i]);
    }
}

BlurX::BlurX(unsigned int _idcnt, unsigned char (&image)[4110080], float (&kernel)[21], int (&center), float (&tempim)[4110080])
    : _specc::behavior(_idcnt), image(image), kernel(kernel), center(center), tempim(tempim)
{   
}

BlurX::~BlurX(void)
{   
}

void BlurX::blur_x(int rows, int cols)
{   
    int c;
    int cc;
    int r;
    float dot;
    float sum;

    if (0)
	printf("   Blurring the image in the X-direction.\n");
    for(r = 0; r < rows; r++ )
    {   
	for(c = 0; c < cols; c++ )
	{   
	    dot = 0.000000000000000e+00;
	    sum = 0.000000000000000e+00;
	    for(cc = ( -center); cc <= center; cc++ )
	    {   
		if (((c + cc) >= 0) && ((c + cc) < cols))
		{   
		    dot += (float)image[r * cols + (c + cc)] * kernel[center + 
		    cc];
		    sum += kernel[center + cc];
		}
	    }
	    tempim[r * cols + c] = dot / sum;
	}
    }
}

void BlurX::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    Tstart = clock();    

    blur_x(1520, 2704);

    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("BlurX time: %f\n", T1);
}

BlurY::BlurY(unsigned int _idcnt, float (&tempim)[4110080], float (&kernel)[21], int (&center), short int (&smoothedim)[4110080])
    : _specc::behavior(_idcnt), tempim(tempim), kernel(kernel), center(center), smoothedim(smoothedim)
{   
}

BlurY::~BlurY(void)
{   
}

void BlurY::blur_y(int rows, int cols)
{   
    int c;
    int r;
    int rr;
    float dot;
    float sum;

    if (0)
	printf("   Blurring the image in the Y-direction.\n");
    for(c = 0; c < cols; c++ )
    {   
	for(r = 0; r < rows; r++ )
	{   
	    sum = 0.000000000000000e+00;
	    dot = 0.000000000000000e+00;
	    for(rr = ( -center); rr <= center; rr++ )
	    {   
		if (((r + rr) >= 0) && ((r + rr) < rows))
		{   
		    dot += tempim[(r + rr) * cols + c] * kernel[center + 
		    rr];
		    sum += kernel[center + rr];
		}
	    }
	    smoothedim[r * cols + c] = (short int)(dot * 9.000000000000000e+01 / 
		sum + 5.000000000000000e-01);
	}
    }
}

void BlurY::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    Tstart = clock();

    blur_y(1520, 2704);

    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("BlurY time: %f\n", T1);
}

Gaussian_Smooth::Gaussian_Smooth(unsigned int _idcnt, i_img_receiver (&ImgIn), short int (&smoothedim)[4110080])
    : _specc::behavior(_idcnt), ImgIn(ImgIn), smoothedim(smoothedim),
    blurX(_IDcnt, image, kernel, center, tempim),
    blurY(_IDcnt, tempim, kernel, center, smoothedim),
    gauss(_IDcnt, kernel, center),
    receive(_IDcnt, ImgIn, image)
{   
}

Gaussian_Smooth::~Gaussian_Smooth(void)
{   
}

void Gaussian_Smooth::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    Tstart = clock();    

    receive.main();
    gauss.main();
    blurX.main();
    blurY.main();

    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Total Gaussian_Smooth time: %f\n", T1);
}

Derivative_X_Y::Derivative_X_Y(unsigned int _idcnt, short int (&smoothedim)[4110080], short int (&delta_x)[4110080], short int (&delta_y)[4110080])
    : _specc::behavior(_idcnt), smoothedim(smoothedim), delta_x(delta_x), delta_y(delta_y)
{   
}

Derivative_X_Y::~Derivative_X_Y(void)
{   
}

void Derivative_X_Y::derivative_x_y(int rows, int cols)
{   
    int c;
    int pos;
    int r;

    if (0)
	printf("   Computing the X-direction derivative.\n");
    for(r = 0; r < rows; r++ )
    {   
	pos = r * cols;
	delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos];
	pos++ ;
	for(c = 1; c < (cols - 1); c++  , pos++ )
	{   
	    delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
	}
	delta_x[pos] = smoothedim[pos] - smoothedim[pos - 1];
    }
    if (0)
	printf("   Computing the Y-direction derivative.\n");
    for(c = 0; c < cols; c++ )
    {   
	pos = c;
	delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos];
	pos += cols;
	for(r = 1; r < (rows - 1); r++  , pos += cols)
	{   
	    delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
	}
	delta_y[pos] = smoothedim[pos] - smoothedim[pos - cols];
    }
}

void Derivative_X_Y::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    Tstart = clock();

    derivative_x_y(1520, 2704);

    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Derivative_X_Y time: %f\n", T1);
}

Magnitude_X_Y::Magnitude_X_Y(unsigned int _idcnt, short int (&delta_x)[4110080], short int (&delta_y)[4110080], short int (&magnitude)[4110080])
    : _specc::behavior(_idcnt), delta_x(delta_x), delta_y(delta_y), magnitude(magnitude)
{   
}

Magnitude_X_Y::~Magnitude_X_Y(void)
{   
}

void Magnitude_X_Y::magnitude_x_y(int rows, int cols)
{   
    int c;
    int pos;
    int r;
    int sq1;
    int sq2;

    for(r = 0 , pos = 0; r < rows; r++ )
    {   
	for(c = 0; c < cols; c++  , pos++ )
	{   
	    sq1 = (int)delta_x[pos] * (int)delta_x[pos];
	    sq2 = (int)delta_y[pos] * (int)delta_y[pos];
	    magnitude[pos] = (short int)(5.000000000000000e-01 + sqrt((float)sq1 + 
		    (float)sq2));
	}
    }
}

void Magnitude_X_Y::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    Tstart = clock();

    magnitude_x_y(1520, 2704);

    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Magnitude_X_Y time: %f\n", T1);
}

Non_Max_Supp::Non_Max_Supp(unsigned int _idcnt, short int (&gradx)[4110080], short int (&grady)[4110080], short int (&mag)[4110080], unsigned char (&nms)[4110080])
    : _specc::behavior(_idcnt), gradx(gradx), grady(grady), mag(mag), nms(nms)
{   
}

Non_Max_Supp::~Non_Max_Supp(void)
{   
}

void Non_Max_Supp::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    unsigned char result[4110080];
    
    Tstart = clock();

    non_max_supp(1520, 2704, result);
    { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<4110080;_scc_index_0++) (nms)[_scc_index_0] = (result)[_scc_index_0]; }
    
    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Non_Max_Supp time: %f\n", T1);
}

void Non_Max_Supp::non_max_supp(int nrows, int ncols, unsigned char *result)
{   
    int colcount;
    int count;
    int rowcount;
    short int *magptr;
    short int *magrowptr;
    short int *gxptr;
    short int *gxrowptr;
    short int *gyptr;
    short int *gyrowptr;
    short int z1;
    short int z2;
    short int gx;
    short int gy;
    short int m00;
    float mag1;
    float mag2;
    float xperp;
    float yperp;
    unsigned char *resultptr;
    unsigned char *resultrowptr;

    for(count = 0 , resultrowptr = result , resultptr = result + ncols * 
	(nrows - 1); count < ncols; resultptr++  , resultrowptr++  , count++ )
    {   
	 *resultrowptr =  *resultptr = (unsigned char)0;
    }
    for(count = 0 , resultptr = result , resultrowptr = result + ncols - 
	1; count < nrows; count++  , resultptr += ncols , resultrowptr += 
	ncols)
    {   
	 *resultptr =  *resultrowptr = (unsigned char)0;
    }
    for(rowcount = 1 , magrowptr = mag + ncols + 1 , gxrowptr = gradx + 
	ncols + 1 , gyrowptr = grady + ncols + 1 , resultrowptr = result + 
	ncols + 1; rowcount <= nrows - 2; rowcount++  , magrowptr += ncols , 
	gyrowptr += ncols , gxrowptr += ncols , resultrowptr += ncols)
    {   
	for(colcount = 1 , magptr = magrowptr , gxptr = gxrowptr , gyptr = 
	    gyrowptr , resultptr = resultrowptr; colcount <= ncols - 2; colcount++  , 
	    magptr++  , gxptr++  , gyptr++  , resultptr++ )
	{   
	    m00 =  *magptr;
	    if (m00 == 0)
	    {   
		 *resultptr = (unsigned char)255;
	    }
	    else 
	    {   
		xperp =  -(gx =  *gxptr) / ((float)m00);
		yperp = (gy =  *gyptr) / ((float)m00);
	    }
	    if (gx >= 0)
	    {   
		if (gy >= 0)
		{   
		    if (gx >= gy)
		    {   
			z1 =  *(magptr - 1);
			z2 =  *(magptr - ncols - 1);
			mag1 = (m00 - z1) * xperp + (z2 - z1) * yperp;
			z1 =  *(magptr + 1);
			z2 =  *(magptr + ncols + 1);
			mag2 = (m00 - z1) * xperp + (z2 - z1) * yperp;
		    }
		    else 
		    {   
			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols - 1);
			mag1 = (z1 - z2) * xperp + (z1 - m00) * yperp;
			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols + 1);
			mag2 = (z1 - z2) * xperp + (z1 - m00) * yperp;
		    }
		}
		else 
		{   
		    if (gx >=  -gy)
		    {   
			z1 =  *(magptr - 1);
			z2 =  *(magptr + ncols - 1);
			mag1 = (m00 - z1) * xperp + (z1 - z2) * yperp;
			z1 =  *(magptr + 1);
			z2 =  *(magptr - ncols + 1);
			mag2 = (m00 - z1) * xperp + (z1 - z2) * yperp;
		    }
		    else 
		    {   
			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols - 1);
			mag1 = (z1 - z2) * xperp + (m00 - z1) * yperp;
			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols + 1);
			mag2 = (z1 - z2) * xperp + (m00 - z1) * yperp;
		    }
		}
	    }
	    else 
	    {   
		if ((gy =  *gyptr) >= 0)
		{   
		    if ( -gx >= gy)
		    {   
			z1 =  *(magptr + 1);
			z2 =  *(magptr - ncols + 1);
			mag1 = (z1 - m00) * xperp + (z2 - z1) * yperp;
			z1 =  *(magptr - 1);
			z2 =  *(magptr + ncols - 1);
			mag2 = (z1 - m00) * xperp + (z2 - z1) * yperp;
		    }
		    else 
		    {   
			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols + 1);
			mag1 = (z2 - z1) * xperp + (z1 - m00) * yperp;
			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols - 1);
			mag2 = (z2 - z1) * xperp + (z1 - m00) * yperp;
		    }
		}
		else 
		{   
		    if ( -gx >  -gy)
		    {   
			z1 =  *(magptr + 1);
			z2 =  *(magptr + ncols + 1);
			mag1 = (z1 - m00) * xperp + (z1 - z2) * yperp;
			z1 =  *(magptr - 1);
			z2 =  *(magptr - ncols - 1);
			mag2 = (z1 - m00) * xperp + (z1 - z2) * yperp;
		    }
		    else 
		    {   
			z1 =  *(magptr + ncols);
			z2 =  *(magptr + ncols + 1);
			mag1 = (z2 - z1) * xperp + (m00 - z1) * yperp;
			z1 =  *(magptr - ncols);
			z2 =  *(magptr - ncols - 1);
			mag2 = (z2 - z1) * xperp + (m00 - z1) * yperp;
		    }
		}
	    }
	    if ((mag1 > 0.000000000000000e+00) || (mag2 > 0.000000000000000e+00))
	    {   
		 *resultptr = (unsigned char)255;
	    }
	    else 
	    {   
		if (mag2 == 0.000000000000000e+00)
		     *resultptr = (unsigned char)255;
		else 
		     *resultptr = (unsigned char)128;
	    }
	}
    }
}

Apply_Hysteresis::Apply_Hysteresis(unsigned int _idcnt, short int (&mag)[4110080], unsigned char (&nms)[4110080], i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), mag(mag), nms(nms), ImgOut(ImgOut)
{   
}

Apply_Hysteresis::~Apply_Hysteresis(void)
{   
}

void Apply_Hysteresis::apply_hysteresis(int rows, int cols, float tlow, 
    float thigh, unsigned char *edge)
{   
    int c;
    int highcount;
    int highthreshold;
    int hist[32768];
    int lowthreshold;
    int numedges;
    int pos;
    int r;
    short int maximum_mag;

    for(r = 0 , pos = 0; r < rows; r++ )
    {   
	for(c = 0; c < cols; c++  , pos++ )
	{   
	    if (nms[pos] == 128)
		edge[pos] = 128;
	    else 
		edge[pos] = 255;
	}
    }
    for(r = 0 , pos = 0; r < rows; r++  , pos += cols)
    {   
	edge[pos] = 255;
	edge[pos + cols - 1] = 255;
    }
    pos = (rows - 1) * cols;
    for(c = 0; c < cols; c++  , pos++ )
    {   
	edge[c] = 255;
	edge[pos] = 255;
    }
    for(r = 0; r < 32768; r++ )
	hist[r] = 0;
    for(r = 0 , pos = 0; r < rows; r++ )
    {   
	for(c = 0; c < cols; c++  , pos++ )
	{   
	    if (edge[pos] == 128)
		hist[mag[pos]]++ ;
	}
    }
    for(r = 1 , numedges = 0; r < 32768; r++ )
    {   
	if (hist[r] != 0)
	    maximum_mag = r;
	numedges += hist[r];
    }
    highcount = (int)(numedges * thigh + 5.000000000000000e-01);
    r = 1;
    numedges = hist[1];
    while((r < (maximum_mag - 1)) && (numedges < highcount))
    {   
	r++ ;
	numedges += hist[r];
    }
    highthreshold = r;
    lowthreshold = (int)(highthreshold * tlow + 5.000000000000000e-01);
    if (0)
    {   
	printf("The input low and high fractions of %f and %f computed to\n", 
	    tlow, thigh);
	printf("magnitude of the gradient threshold values of: %d %d\n", 
	    lowthreshold, highthreshold);
    }
    for(r = 0 , pos = 0; r < rows; r++ )
    {   
	for(c = 0; c < cols; c++  , pos++ )
	{   
	    if ((edge[pos] == 128) && (mag[pos] >= highthreshold))
	    {   
		edge[pos] = 0;
		follow_edges((edge + pos), (mag + pos), lowthreshold, 
		    cols);
	    }
	}
    }
    for(r = 0 , pos = 0; r < rows; r++ )
    {   
	for(c = 0; c < cols; c++  , pos++ )
	    if (edge[pos] != 0)
		edge[pos] = 255;
    }
}

void Apply_Hysteresis::follow_edges(unsigned char *edgemapptr, short int *edgemagptr, 
    short int lowval, int cols)
{   
    short int *tempmagptr;
    unsigned char *tempmapptr;
    int i;
    int x[8] = 
    { 1,1,0,-1,-1,-1,0,1 };
    int y[8] = 
    { 0,1,1,1,0,-1,-1,-1 };

    for(i = 0; i < 8; i++ )
    {   
	tempmapptr = edgemapptr - y[i] * cols + x[i];
	tempmagptr = edgemagptr - y[i] * cols + x[i];
	if (( *tempmapptr == 128) && ( *tempmagptr > lowval))
	{   
	     *tempmapptr = (unsigned char)0;
	    follow_edges(tempmapptr, tempmagptr, lowval, cols);
	}
    }
}

void Apply_Hysteresis::main(void)
{   
    clock_t Tstart, Tstop;
    double T1 = 0.0;
    Tstart = clock();
    
    apply_hysteresis(1520, 2704, 3.000000000000000e-01, 8.000000000000000e-01, 
	EdgeImage);
    ImgOut.send(EdgeImage);
 
    Tstop = clock();
    T1 = (double) (Tstop-Tstart)/CLOCKS_PER_SEC;
    printf("Apply_Hysteresis time: %f\n", T1);
}

DUT::DUT(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut),
    apply_hysteresis(_IDcnt, magnitude, nms, ImgOut),
    derivative_x_y(_IDcnt, smoothedim, delta_x, delta_y),
    gaussian_smooth(_IDcnt, ImgIn, smoothedim),
    magnitude_x_y(_IDcnt, delta_x, delta_y, magnitude),
    non_max_supp(_IDcnt, delta_x, delta_y, magnitude, nms)
{   
}

DUT::~DUT(void)
{   
}

void DUT::main(void)
{   
    {
	enum { _scc_state_0, _scc_state_gaussian_smooth, _scc_state_derivative_x_y, _scc_state_magnitude_x_y, _scc_state_non_max_supp, _scc_state_apply_hysteresis } _scc_next_state = _scc_state_gaussian_smooth;
	do switch(_scc_next_state)
	{
	    case _scc_state_gaussian_smooth: { gaussian_smooth.main();
		{ _scc_next_state = _scc_state_derivative_x_y; break; }
		}
	    case _scc_state_derivative_x_y: { derivative_x_y.main();
		{ _scc_next_state = _scc_state_magnitude_x_y; break; }
		}
	    case _scc_state_magnitude_x_y: { magnitude_x_y.main();
		{ _scc_next_state = _scc_state_non_max_supp; break; }
		}
	    case _scc_state_non_max_supp: { non_max_supp.main();
		{ _scc_next_state = _scc_state_apply_hysteresis; break; }
		}
	    case _scc_state_apply_hysteresis: { apply_hysteresis.main();
		{ _scc_next_state = _scc_state_gaussian_smooth; break; }
		}
	    case _scc_state_0: { _scc_next_state = _scc_state_0; break; }
	}
	while(_scc_next_state != _scc_state_0);
    }
}

Platform::Platform(unsigned int _idcnt, i_img_receiver (&ImgIn), i_img_sender (&ImgOut))
    : _specc::behavior(_idcnt), ImgIn(ImgIn), ImgOut(ImgOut),
    _scc_const_port_0(1ul),
    _scc_const_port_1(1ul),
    canny(++_IDcnt, q1, q2),
    din(++_IDcnt, ImgIn, q1),
    dout(++_IDcnt, q2, ImgOut),
    q1(_scc_const_port_0),
    q2(_scc_const_port_1)
{   
}

Platform::~Platform(void)
{   
}

void Platform::main()
{   
    {
	_specc::fork _scc_fork_0(&din), _scc_fork_1(&canny), _scc_fork_2(&dout);
	_specc::par(&_scc_fork_0,
	    &_scc_fork_1,
	    &_scc_fork_2,
	    ((_specc::fork*)0));
    }
}

Main::Main(unsigned int _idcnt)
    : _specc::class_type(_idcnt),
    _scc_const_port_0(1ul),
    _scc_const_port_1(1ul),
    monitor(_IDcnt, q2),
    platform(_IDcnt, q1, q2),
    q1(_scc_const_port_0),
    q2(_scc_const_port_1),
    stimulus(_IDcnt, q1)
{   
}

Main::~Main(void)
{   
}

int Main::main(void)
{   
    {
	_specc::fork _scc_fork_0(&stimulus), _scc_fork_1(&platform), _scc_fork_2(&monitor);
	_specc::par(&_scc_fork_0,
	    &_scc_fork_1,
	    &_scc_fork_2,
	    ((_specc::fork*)0));
    }
    return 0;
}

Main _scc_main(_IDcnt);

int main(void)
{   
    int _scc_main_return;
    
    _specc::start();
    _scc_main_return = _scc_main.main();
    _specc::end();
    return(_scc_main_return);
}

void _scc_bit4_err_handle(
    const _bit4& bit4vec)
{   
    char temp_bits[1024], *p;
    p=bit2str(2,&temp_bits[1023], bit4vec);
    _specc::abort(
	"ERROR:\t Casting a bit4 vector failed \n"
	"Bit4 vector contains X/Z values %s\n"
	"Simulation aborted.\n", p);
	
}

// EOF canny.cpp
