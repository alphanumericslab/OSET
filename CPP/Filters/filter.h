#ifndef FILTER_H_
#define FILTER_H_

//#define SWAP(a,b) temp = (a);(a) = (b);(b) = temp;
template <class T>
void SWAP(T &a,T &b)
{
    T temp;
    
    temp = a;
    a = b;
    b = temp;
}
//abstract filter class
class Filter
{
public:
	virtual int process(const double *input_vector, double *output_vector,const int length) = 0;
	virtual ~Filter(){}
	
};

#endif /*FILTER_H_*/
