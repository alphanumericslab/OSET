//#define SWAP(a,b) temp = (a);(a) = (b);(b) = temp; 

class item
{
public:
	double value;
	int from;
	int to;
};

class list
{
private:
	item* itmlist;
	int ptr;
	int len;
	int start,stop;

	int resort(void);

public:
	list(int N,double init=0);
	~list(void);
	
	int initialize(double x);
	int insert(double data);
	
	double mean(int alpha, int beta);
	double wmedian(const double *weight, int m, int alpha, int beta);
	double medianOdd(void);
	double medianEven(void);
};