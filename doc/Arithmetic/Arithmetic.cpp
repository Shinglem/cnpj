// Arithmetic.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"
#include "Arithmetic.h"
#include "malloc.h"
#include "math.h"
const double PI=3.1415926535897932384626;
const double MinDoubleNum = 10e-10;

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}

// 导数计算公式程序
int __stdcall Derivate(double *pfX, double *pfY, int iNum, int bLn, int iN, double *pfdY)
              //bLn = 0    直接计算导数dY/dX    ----X = Ln(T)
              //bLn = 1    计算LnX导数dY/dLnX   ----X = T
              //bLn = 2    计算导数dY/dX*2.303  ----X = Log10(T)
{
	int    i, j, k, m, iBasej;
	double fA, fB, fC, fD, fX;

	if( iNum < 3 ) return 1;
	if( iNum < iN ) iN = iNum;
	for( k=0; k<iNum; k++ )
	{
		if( k < (iN-1)/2 )
		{
			//计算导数
			iBasej = 0;
			fX = pfX[k];
			fD = 0;
			for( j=0; j<iN; j++ )//和式
			{
				//计算插值基函数的导数
				fB = 0;
				for( m=0; m<iN; m++ )//和式
				{//m!=j---所有积式
					//计算积式
					if(m==j) continue;
					fC = 1;
					for(i=0; i<iN;i++)
					{//i!=j, i!=m    x-xi
						if(i==j) continue;
						if(i==m) continue;
						fC = fC*(fX-pfX[iBasej+i]);
					}
					fB = fB+fC;
				}
				//计算积式   xj-xi
				fA = 1;
				for( i=0; i<iN; i++ )
				{//i!=j
					if(i==j) continue;
					fA=fA*(pfX[iBasej+j]-pfX[iBasej+i]);
				}
				fD = fD+fB/fA*pfY[iBasej+j];
			}
			pfdY[k] = (float)fD;
			continue;
		}
		if( iNum-1-k < (iN-1)/2 )
		{
			//计算导数
			iBasej = iNum-iN;
			fX = pfX[k];
			fD = 0;
			for( j=0; j<iN; j++ )//和式
			{
				//计算插值基函数的导数
				fB = 0;
				for( m=0; m<iN; m++ )//和式
				{//m!=j---所有积式
					//计算积式
					if(m==j) continue;
					fC = 1;
					for(i=0;i<iN;i++)
					{//i!=j, i!=m    x-xi
						if(i==j) continue;
						if(i==m) continue;
						fC = fC*(fX-pfX[iBasej+i]);
					}
					fB = fB+fC;
				}
				//计算积式   xj-xi
				fA = 1;
				for( i=0; i<iN; i++ )
				{//i!=j
					if(i==j) continue;
					fA=fA*(pfX[iBasej+j]-pfX[iBasej+i]);
				}
				fD = fD+fB/fA*pfY[iBasej+j];
			}
			pfdY[k] = (float)fD;
			continue;
		}

		//计算导数
		iBasej = k-(iN-1)/2;
		fX = pfX[k];
		fD = 0;
		for( j=0; j<iN; j++ )//和式
		{
			//计算插值基函数的导数
			fB = 0;
			for( m=0; m<iN; m++ )//和式
			{//m!=j---所有积式
				//计算积式
				if(m==j) continue;
				fC = 1;
				for( i=0; i<iN; i++ )
				{//i!=j, i!=m    x-xi
					if(i==j) continue;
					if(i==m) continue;
					fC = fC*(fX-pfX[iBasej+i]);
				}
				fB = fB+fC;
			}
			//计算积式   xj-xi
			fA = 1;
			for( i=0; i<iN; i++ )
			{//i!=j
				if(i==j) continue;
				fA=fA*(pfX[iBasej+j]-pfX[iBasej+i]);
			}
			fD = fD+fB/fA*pfY[iBasej+j];
		}
		pfdY[k] = (float)fD;
	}

	if( bLn == 1)
		for( k=0; k<iNum; k++ )
			pfdY[k] = (float)fabs((float)pfdY[k]*pfX[k]);
	else if( bLn == 2)
		for( k=0; k<iNum; k++ )
			pfdY[k] = (float)fabs((float)pfdY[k]/log(10.0));

	for( k=0; k<iNum; k++ )
	{
		pfdY[k] = (float)fabs(pfdY[k]);
		if( pfdY[k] < 1e-4 )
			pfdY[k] = float(1e-4);
		if( pfdY[k] > 1e6 )
			pfdY[k] = float(1e6);
	}
	return 0;
}
//I0计算公式程序 
double __stdcall I0(double x)
{
	double p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9;
	double y,aaa,bbb,ax,temp;
    p1 = 1.0;             p2 = 3.5156229;
    p3 = 3.0899424;       p4 = 1.2067492;
    p5 = 0.2659732;       p6 = 0.0360768;
    p7 = 0.0045813;
    q1 = 0.39894228;      q2 = 0.01328592;
    q3 = 0.00225319;      q4 = -0.00157565;
    q5 = 0.00916281;      q6 = -0.02057706;
    q7 = 0.02635537;      q8 = -0.01647633;
    q9 = 0.00392377;
    if( fabs(x) < 3.75)
	{
        y = (x / 3.75)*(x / 3.75);
        aaa = y * (p5 + y * (p6 + y * p7));
        temp = p1 + y * (p2 + y * (p3 + y * (p4 + aaa)));
	}
    else
	{
        ax = fabs(x);
        y = 3.75 / ax;
		if(ax<709)
		{
           aaa = exp(ax) / sqrt(ax);
		}
		if(ax>=709)
		{
           aaa = exp(709) / sqrt(ax);
		}
        bbb = q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))));
        temp = aaa * (q1 + y * (q2 + y * (q3 + y * bbb)));
	}
	return temp;
}
//K0计算公式程序 
double __stdcall K0(double x)
{
	double p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7;
	double y,aaa,bbb,temp;
    p1 = -0.57721566;      p2 = 0.4227842;
    p3 = 0.23069756;       p4 = 0.0348859;
    p5 = 0.00262698;       p6 = 0.0001075;
    p7 = 0.0000074;
    q1 = 1.25331414;       q2 = -0.07832358;
    q3 = 0.02189568;       q4 = -0.01062446;
    q5 = 0.00587872;       q6 = -0.0025154;
    q7 = 0.00053208;
    if( x <= 2.0)
	{
        y = x * x / 4.0;
        bbb = y * (p5 + y * (p6 + y * p7));
        aaa= p1 + y * (p2 + y * (p3 + y * (p4 + bbb)));
        temp = (-log(x / 2.0) * I0(x)) + aaa;
	}
    else
	{
        y = (2.0 / x);
        bbb = y * (q5 + y * (q6 + y * q7));
        aaa = q1 + y * (q2 + y * (q3 + y * (q4 + bbb)));
        temp = (exp(-x) / sqrt(x)) * aaa;
	}
	return temp;
}
//I1计算公式程序 
double __stdcall I1(double x)
{
	double p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9;
	double y,aaa,bbb,ax,temp;
    p1 = 0.5;              p2 = 0.87890594;
    p3 = 0.51498869;       p4 = 0.15084934;
    p5 = 0.02658733;       p6 = 0.00301532;
    p7 = 0.00032411;
    q1 = 0.39894228;       q2 = -0.03988024;
    q3 = -0.00362018;      q4 = 0.00163801;
    q5 = -0.01031555;      q6 = 0.02282967;
    q7 = -0.02895312;      q8 = 0.01787654;
    q9 = -0.00420059;
    if (fabs(x) < 3.75 )
	{
        y = (x / 3.75) *(x / 3.75) ;
        aaa = y * (p4 + y * (p5 + y * (p6 + y * p7)));
        temp = x * (p1 + y * (p2 + y * (p3 + aaa)));
	}
    else
	{
        ax = fabs(x);
        y = 3.75 / ax;
        if(ax<709)
		{
           aaa = exp(ax) / sqrt(ax);
		}
		if(ax>=709)
		{
           aaa = exp(709) / sqrt(ax);
		}
        bbb = y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))));
        temp = aaa * (q1 + y * (q2 + y * (q3 + y * (q4 + bbb))));
	}
    return temp;
}
//K1计算公式程序 
double __stdcall K1(double x)
{
	double p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7;
	double y,aaa,bbb,ccc,temp; 
    p1 = 1.0;               p2 = 0.15443144;
    p3 = -0.67278579;       p4 = -0.18156897;
    p5 = -0.01919402;       p6 = -0.00110404;
    p7 = -0.00004686;
    q1 = 1.25331414;        q2 = 0.23498619;
    q3 = -0.0365562;        q4 = 0.01504268;
    q5 = -0.00780353;       q6 = 0.00325614;
    q7 = -0.00068245;
    if (x <= 2.0)
	{
        y = x * x / 4;
        aaa = log(x / 2.0) * I1(x);
        ccc = y * (p5 + y * (p6 + y * p7));
        bbb = p1 + y * (p2 + y * (p3 + y * (p4 + ccc)));
        temp = aaa + (1.0 / x) * bbb;
	}
    else
	{
        y = (2.0 / x);
        bbb = y * (q5 + y * (q6 + y * q7));
        aaa = q1 + y * (q2 + y * (q3 + y * (q4 + bbb)));
        temp = (exp(-x) / sqrt(x)) * aaa;
	}    
	return temp;
}
//计算K0从x到无穷积分
double __stdcall IntegelK0x(double fx)
{
	if(fx<=2)
	{
		double fck[7];
		double fdk[7];
		double fa=0.0f;
		double fb=0.0f;
		fck[0]=2;fck[1]=0.666666667;fck[2]=0.100000003;fck[3]=0.007936494;
		fck[4]=0.000385833;fck[5]=0.000012590;fck[6]=0.000000319;
		fdk[0]=0.84556868;fdk[1]=0.50407836;fdk[2]=0.11227902;fdk[3]=0.01110118;
		fdk[4]=0.00062664;fdk[5]=0.00002069;fdk[6]=0.00000116;
        for(int i=0;i<=6;i++)
		{
			fa+=fck[i]*pow(fx/2.0,2*i+1);
		}
		for(i=0;i<=6;i++)
		{
			fb+=fdk[i]*pow(fx/2.0,2*i+1);
		}
		return PI/2.0+log(fx/2.0)*fa-fb;
	}
	else if(fx>2 && fx<=4)
	{
		double fak[5];
		double fa=0.0f;
		fak[0]=1.2494934;fak[1]=0.3584641;fak[2]=0.1859840;
		fak[3]=0.0781715;fak[4]=0.0160395;
		for(int i=0;i<=4;i++)
		{
			fa+=fak[i]*pow(fx/2.0,-i)*pow(-1,i);
		}
		return pow(fx,-0.5)*exp(-fx)*fa;
	}
	else if(fx>4 && fx<=7)
	{
		double fek[7];
		double fe=0.0f;
		fek[0]=1.25331414;fek[1]=0.11190289;fek[2]=0.02576646;fek[3]=0.00933994;
		fek[4]=0.00417454;fek[5]=0.00163271;fek[6]=0.00033934;
		for(int i=0;i<=6;i++)
		{
			fe+=pow(-1,i)*fek[i]*pow(fx/7.0,-i);
		}
		return pow(fx,-0.5)*exp(-fx)*fe;
	}
	else
	{
		double fck[7];
		double fc=0.0f;
		fck[0]=1.25331414;fck[1]=0.11190289;fck[2]=0.02576646;fck[3]=0.00933994;
		fck[4]=0.00417454;fck[5]=0.00163271;fck[6]=0.00033934;
		for(int i=0;i<=6;i++)
		{
			fc+=pow(-1,i)*fck[i]*pow(fx/7.0,-i);
		}
		return pow(fx,-0.5)*exp(-fx)*fc;
	}
}
//累乘计算
double __stdcall multp(int nx)
{
     int num,i;
	 num=1;
	 for(i=1;i<=nx;i++)
	 {
	     num*=i;   
	 }
	 return num;
}
//计算K0从0到x积分
double __stdcall IntegelK0(double fx)
{
	if(fx>20)
		return PI/2;
    else
	{
		double fKresult=0;
		double fMresult1=0;
		double fMresult2=0;
		double fMresult3=0;
		double ferror=0.000001;
		double fa=0;
		double fb=0;
		double fc=0;
		int nN=0;
		int nStop=0;
		int i=0;
		double fMR;
		nStop=int(fx/2)+1;
		fa=1;
		fMresult1=fa;
		i=0;
		do
		{
			i++;
			fa=0;
			fa=pow(fx/2,2*i)/(pow(multp(i),2)*(2*i+1));
			fMresult1+=fa;
		}while(fabs((log(fx/2)+0.5772)*fx*fa)>ferror && i<nStop);
		fb=1;
		fMresult2=fb;
		i=0;
		do
		{
			i++;
			fb=0;
			fb=pow(fx/2,2*i)/(pow(multp(i),2)*pow(2*i+1,2));
			fMresult2+=fb;		
		}while(fabs(fx*fb)>ferror  && i<nStop);
		i=0;
		do
		{
			i++;
			fMR=0;
			for(nN=1;nN<=i;nN++)
			{
				fMR=fMR+1/((double)(nN));
			}
			fc=0;
			fc=fMR*pow(fx/2,2*i)/(pow(multp(i),2)*(2*i+1));
			fMresult3+=fc;
		}while(fabs(fx*fc)>ferror  && i<nStop);
		fKresult=fx*(fMresult2+fMresult3)-(log((double)(fx/2))+0.5772)*fx*fMresult1;
		return fKresult;
	}
}

//Ki2积分
double __stdcall Ki2(double fu)
{
	double fa=0.0f;
	fa=-fu*IntegelK0x(fu)+fu*K1(fu);
	return fa;
}

//计算I0从0到x积分
double __stdcall IntegelI0(double fx)
{
	double fak[5];
	double fbk[7];
	double fck[7];
	double fa=0.0f;
	if(fx>0 && fx<=2)
	{
		fck[0]=2;fck[1]=0.666666667;fck[2]=0.100000003;fck[3]=0.007936494;
    	fck[4]=0.000385833;fck[5]=0.000012590;fck[6]=0.000000319;
		for(int i=0;i<=6;i++)
		{
			fa+=fck[i]*pow(fx/2.0,2*i+1);
		}
		return fa;
	}
	else if(fx>2 && fx<=5)
	{
		fck[0]=5.000000;fck[1]=10.416666367;fck[2]=9.765629849;fck[3]=4.844024624;
		fck[4]=1.471860153;fck[5]=0.300704878;fck[6]=0.044686921;fck[7]=0.004500642;
		fck[8]=0.000594340;
		for(int i=0;i<=6;i++)
		{
			fa+=fck[i]*pow(fx/5.0,2*i+1);
		}
		return fa;
	}
	else if(fx>5 && fx<8)
	{
		fak[0]=0.41612;fak[1]=-0.0302912;fak[2]=0.1294122;
		fak[3]=-0.0202292;fak[4]=-0.0151660;
		for(int i=0;i<=4;i++)
		{
			fa+=fak[i]*pow(fx/5.0,-i);
		}
		return pow(fx,-0.5)*exp(fx)*fa;
	}
	else
	{
		fbk[0]=0.3989423;fbk[1]=0.0311734;fbk[2]=0.0059191;fbk[3]=0.0055956;
		fbk[4]=-0.0114858;fbk[5]=0.0177440;fbk[6]=-0.0073995;
		for(int i=0;i<=6;i++)
		{
			fa+=fbk[i]*pow(fx/8.0,-i);
		}
		return pow(fx,-0.5)*exp(fx)*fa;
	}
}

//自变量范围为1——2的Gama函数(1<=x<=2)
double __stdcall Gam1(double x )       
{
	 double tmp, gama, ser;

	 tmp  = x + 4.5;
	 tmp  = ( x - 0.5 ) * logl( tmp ) - tmp;
	 ser  = 1.0+76.18009173/x-86.50532033/(x+1.0)+24.01409822/(x+2.0)-1.231739516/(x+3.0)+0.00120858003/(x+4.0)-5.36382e-6/(x+5.0);
	 gama = expl( tmp + logl( 2.50662827465*ser ) );

	 return(gama);
}

//Gama函数
double __stdcall Gamm(double x )
{
	 double gama, Sum = 1.0;

	 double xb = x;
	 if( x < 0 )
		 x = -x;

	 if( x >= 2.0 )
	 {
		  do
		  {
				x   = x - 1;
				Sum = Sum * x;
		  }while( x >= 2.0 );
		  gama = Sum * Gam1( x );
	 }
	 else
	     if( x < 1.0 )
            gama = Gam1( 1+x ) / x;
	     else
            gama = Gam1( x );

	if( xb < 0 )
		gama = -PI/(x*sin(PI*x)*gama);

	return( gama );
}

//自变量范围为2——3的Gama函数(2<=x<=3)
double __stdcall Gam2To3(double x )       
{
	 double fSum = 0;
	 double pfAi[11]={0.0000677106, -0.0003442342, 0.0015397681, -0.0024467480, 0.0109736958, -0.0002109075,
	                  0.0742379071,  0.0815782188, 0.4118402518,  0.4227843370, 1.0 };
	 for(int i=0; i<11; i++)
		 fSum += (pfAi[i]*pow(x-2, 10-i)); 
	 return fSum;
}

//分数阶Bessel函数Iv
double __stdcall Iv(double v, double z)
{
	 int    k = 0;
	 long double Iv, tmp, fGamm, sum = 0;

	 do
	 {
		  fGamm = Gamm( v+k+1 );
		  tmp = powl( z/2, 2*k ) / multp(k) / fGamm;
		  sum = sum + tmp;
		  k   = k + 1;
	 }while( tmp/sum>1e-18 );
	 Iv = powl( z/2.0, v ) * sum;

	 return( Iv );
}

//分数阶Bessel函数Kv
double __stdcall Kv(double v, double z)
{
	double fV2 = Iv( v, z );
	double fV1 = Iv( (-1.0)*v, z );
	double fV3 = 2*sin( v*PI );
	return PI * (fV1-fV2) / fV3;
}

//stefest计算
double __stdcall V(int nx)
{
	int uplimit,downlimit;
	uplimit = 0;
	downlimit = 0;
	downlimit = floor((nx+1)/2);
	if (nx>=4)
	{
	    uplimit=4;
	}
	else
	{
	    uplimit=nx;
	}
	int num;
	num=0;
	double numv;
	numv=0.0f;
	for(num=downlimit;num<=uplimit;num++)
	{
	     numv+=pow(num,5)*multp(2*num)/(multp(4-num)*pow(multp(num),2)*multp(nx-num)*multp(2*num-nx));  
	}
	return pow(-1,4+nx)*numv;

}

//数据磨光处理求导
double __stdcall Spole(double x, double x1, double x2, double y1, double y2)
{
	double	y;

	y=(y2-y1)*(x-x1)/(x2-x1)+y1;

	return(y);
}

double __stdcall Spole1(int jl1, int jr1, int jl2, int jr2, double *x, double *y)
{
	double	dx1,dx2,dy1,dy2,m1,m2;
	double aaa;

	dx1=x[jr1]-x[jl1];
	dx2=x[jr2]-x[jl2];
	dy1=y[jr1]-y[jl1];
	dy2=y[jr2]-y[jl2];
	m1=dy1*dx2/dx1;
	m2=dy2*dx1/dx2;

	aaa=(m1+m2)/(dx1+dx2);
	return(aaa);
}

//数据磨光处理求导
int __stdcall DataDerivateSmooth(	double *td,      //输入时间数据序列点
									double *pd,      //输入压力数据序列点
									int    npoint,   //数据序列点数
									int    Flag,     //0为对数求导,1为自然求导
									double l,
									double *pdd)     //输出压力导数数据序列点        //磨光系数
{
	int	i,jl1,jr1,jl2,jr2,k,k1;
	double	*h;
	double x1,x2;

	h=new double[npoint+1];

	jl1 = 0;
	jr1 = 1;
	jl2 = 1;
	jr2 = 2;
	k1  = 0;
	 	 
	for(i=1; i<npoint-1; i++)
	{
		if(l<=0.0 || fabs(log(td[i])-log(td[0]))<=l)
		{
			jl1=i-1;
			jr1=i;
			jl2=i;
			jr2=i+1;
			h[i]=Spole1(jl1,jr1,jl2,jr2,td,pd);
			continue;
		}
		else
		{
			if(fabs(log(td[npoint-1])-log(td[i]))>l)
			{
		 	    x1=log(td[i])-l;
				x2=log(td[i])+l;
				for(k=1; k<npoint-1; k++)
				{
				    if(x1>=log(td[k]) && x1<=log(td[k+1]))
					{
						jl1=k;
						jr1=i;
					}
					if(x2>=log(td[k]) && x2<=log(td[k+1]))
					{
						jl2=i;
						jr2=k+1;
					}
			    }
				h[i]=Spole1(jl1,jr1,jl2,jr2,td,pd);
				continue;
			}
			else
			{
				for(k=1; k<npoint-1; k++)
				{
				    if(x1>=log(td[k]) && x1<=log(td[k+1]))
					{
						jl1=k;
						jr1=i;
					}
				}
				k1=k1+1;
				jl2=i-k1;
				jr2=npoint-1;
				h[i]=Spole1(jl1,jr1,jl2,jr2,td,pd);
				continue;
			}
		}
	}
	pdd[0]=(float)Spole(td[0],td[1],td[2],h[1]*td[1],h[2]*td[2]);
	pdd[npoint-1]=(float)Spole(td[npoint-1],td[npoint-2],td[npoint-3],h[npoint-2]*td[npoint-2],h[npoint-3]*td[npoint-3]);
	
	for(i=1; i<npoint-1; i++)
	{
		pdd[i]=float(h[i]*td[i]);
	}

	if(Flag == 1)
	{
        for(i=0; i<npoint; i++)
	    {
		    pdd[i]=float(pdd[i]/td[i]);
	    }
	}
	 
    for(i=0; i<npoint; i++)
	{
		//pdd[i] = (float)fabs(pdd[i]);
		if(pdd[i]<1e-3)
			pdd[i] = float(1e-3);
	}

	delete []h;
	return 0;
}
//转置矩阵
int __stdcall funJZZZ(double** fA, double** fB, int NumX, int NumY)
{
	for(int i=0;i<NumX;i++)
	{
		for(int j=0;j<NumY;j++)
		{
			fB[i][j]=fA[j][i];
		}
	}
	return 0;
}

//矩阵的乘积
int __stdcall funJZCJ(double** fA, double** fB, double** fC, int NumAX, int NumAY, int NumBX, int NumBY)
{
	for(int i=0;i<NumAY;i++)
	{
		for(int j=0;j<NumBX;j++)
		{
			double fa;
			fa=0;
			for(int z=0;z<NumAX;z++)
			{
				fa=fa+fA[z][i]*fB[j][z];
			}
			fC[j][i]=fa;
		}
	}
	return 0;
}

//矩阵的逆
int __stdcall funJZN(double** fA, double** fB, int NumX)
{
	double fa=0.0f;
    for(int i=0;i<NumX;i++)
	{
		for(int j=0;j<NumX;j++)
		{
			if(i==j)
			{
				fB[i][j]=1;
			}
			else
			{
				fB[i][j]=0;
			}
		}
	}
 	for(i=0;i<NumX;i++)//向下
	{
		if(fA[i][i]==0)
		{
			for(int j=0;j<NumX;j++)
			{
				if(fA[i][j]!=0)
				{
					for(int z=0;z<NumX;z++)
					{
						fB[z][i]=fB[z][i]+fB[z][j];
						fA[z][i]=fA[z][i]+fA[z][j];
					}
					j=NumX;
				}
			}
		}
		fa=fA[i][i];
		for(int j=0;j<NumX;j++)//向下
		{
			fB[j][i]=fB[j][i]/fa;
			fA[j][i]=fA[j][i]/fa;//将对角系数化为1
		}
		for(j=i+1;j<NumX;j++)//向下
		{
			fa=fA[i][j];
			for(int z=0;z<NumX;z++)
			{
				fB[z][j]=fB[z][j]-fB[z][i]*fa;
				fA[z][j]=fA[z][j]-fA[z][i]*fa;//下三角化为0
			}
		}
	}
	for(i=NumX-1;i>0;i--)//上三角化为0
	{
		for(int j=i-1;j>-1;j--)//向上
		{
			fa=fA[i][j];
			for(int z=0;z<NumX;z++)
			{
				fB[z][j]=fB[z][j]-fB[z][i]*fa;
				fA[z][j]=fA[z][j]-fA[z][i]*fa;//上三角化为0
			}
		}
	}
	return 0;	
}

// 计算拟稳态串流
double __stdcall ffs(double fu, GJJS_Parameter* GJJSParameter, double fCDe2S, double fS)
{
	double fw;
	double fnamda;
	fw = GJJSParameter->fResWomiga;
	fnamda = GJJSParameter->fResNamuda;
	double fa;
	double fb;
	double fc;
	fa=exp(-2.0*fS)*fnamda+fw*(1-fw)*fu/fCDe2S;
	fb=exp(-2.0*fS)*fnamda+(1-fw)*fu/fCDe2S;
	fc=fa/fb;
	return fc;
}

double __stdcall ffs3(double fu, GJJS_Parameter* GJJSParameter, double fCDe2S, double fS)
{
	double fw1;
	double fw2;
	double fnamda1;
	double fnamda2;
	double fkesai1;
	double fkesai2;
	double fdelta1;
	double fdelta2;
	double fa;
	double fb;
	double fc;
	double fd;
	fw1=GJJSParameter->fResWomiga1;
	fw2=GJJSParameter->fResWomiga2;
	fnamda1=GJJSParameter->fResNamuda1;
	fnamda2=GJJSParameter->fResNamuda2;
	fkesai1=0.0f;
	fkesai2=0.0f;
	fdelta1=0.0f;
	fdelta2=0.0f;
	fa=0.0f;
	fb=0.0f;
	fc=0.0f;
	fd=0.0f;
	fa=1-fw1-fw2;
	fb=(exp(-2.0*fS)*fnamda2*(1-fw1)/fw2)+(exp(-2.0*fS)*fnamda1*(1-fw2)/fw1);
	fc=exp(-2.0*fS)*fnamda1*exp(-2.0*fS)*fnamda2/(fw1*fw2);
	fkesai1=(fb-sqrt(fb*fb-4*fa*fc))/(2*fa);
	fkesai2=(fb+sqrt(fb*fb-4*fa*fc))/(2*fa);
	fdelta1=exp(-2.0*fS)*fnamda1/fw1;
	fdelta2=exp(-2.0*fS)*fnamda2/fw2;
	fd=fa*(fu/fCDe2S+fkesai1)*(fu/fCDe2S+fkesai2)/((fu/fCDe2S+fdelta1)*(fu/fCDe2S+fdelta2));
	return fd;
}

// 在拉氏空间中压力计算公式程序
double __stdcall laplas(double fu, double fCD, double fS, GJJS_Parameter* GJJSParameter, double frd)
{
	double fa=0.0f;
	double fb=0.0f;
	double fc=0.0f;
	double fd=0.0f;
	double fe=0.0f;
	double ff=0.0f;
	double fg=0.0f;
	double fcde2s=0.0f;
    double fvalue=0.0f;
	fcde2s=fCD*exp(2.0*fS);
	if(GJJSParameter->modeltype==1)//均质
	{
		fa=K0(frd*sqrt(ffsy(fu,GJJSParameter)/fcde2s));
    	fb=fu*(fu*K0(sqrt(ffsy(fu,GJJSParameter)/fcde2s))+sqrt(ffsy(fu,GJJSParameter)/fcde2s)*K1(sqrt(ffsy(fu,GJJSParameter)/fcde2s)));
    	fvalue=fa/fb;
	}
	else if(GJJSParameter->modeltype==2)//双重
	{
		fa=K0(frd*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s));
    	fb=fu*(fu*K0(sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s))+
			sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s)*K1(sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s)));
    	fvalue=fa/fb;
	}
	else if(GJJSParameter->modeltype==3)//无限
	{
		fa=-(2*frd*sqrt(ffsy(fu,GJJSParameter)))*IntegelK0x(2*frd*sqrt(ffsy(fu,GJJSParameter)))+2*frd*sqrt(ffsy(fu,GJJSParameter))*K1(2*frd*sqrt(ffsy(fu,GJJSParameter)));
        fb=(PI-(1-fa)/(frd*sqrt(ffsy(fu,GJJSParameter))))/(2*frd*fu*sqrt(ffsy(fu,GJJSParameter)));

		double fpa=0.0f;//stephest算法
    	double fpb=0.0f;
    	double fpc=0.0f;
    	double fpd=0.0f;
    	fpa=fb;
    	fpb=fu*fpa+fS;
    	fpc=fu*(1.0+fu*fCD*(fu*fpa+fS));
    	fvalue=fpb/fpc;
	}
	else if(GJJSParameter->modeltype==9)//无限-双重
	{
		fa=-(2*frd*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu))*IntegelK0x(2*frd*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu))+2*frd*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu)*K1(2*frd*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu));
        fb=(PI-(1-fa)/(frd*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu)))/(2*frd*fu*pow(ffs(fu,GJJSParameter,fcde2s,fS)*fu,0.5));

		double fpa=0.0f;//stephest算法
    	double fpb=0.0f;
    	double fpc=0.0f;
    	double fpd=0.0f;
    	fpa=fb;
    	fpb=fu*fpa+fS;
    	fpc=fu*(1.0+fu*fCD*(fu*fpa+fS));
    	fvalue=fpb/fpc;
	}
	else if(GJJSParameter->modeltype==4)//有限
	{
		fa=1.0/(fu*sqrt(GJJSParameter->fResFcd)*sqrt(1.0/(PI*fS+PI*PI/(2.0*sqrt(ffsy(fu,GJJSParameter)))))+fCD*fu*fu);
		fvalue=fa;
	}
    else if(GJJSParameter->modeltype==10)//有限-双重
	{
		fa=1.0/(fu*sqrt(GJJSParameter->fResFcd)*sqrt(1.0/(PI*fS+PI*PI/(2.0*sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu))))+fCD*fu*fu);
		fvalue=fa;
	}
	else if(GJJSParameter->modeltype==5)//三重
	{
		fa=K0(frd*sqrt(ffs3(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s));
    	fb=fu*(fu*K0(sqrt(ffs3(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s))+
			sqrt(ffs3(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s)*K1(sqrt(ffs3(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s)));
    	fvalue=fa/fb;
	}
	else if(GJJSParameter->modeltype==6)//部分射孔底水
	{
		double fe=0.0f;
       	double fg=0.0f;
    	double ftop=0.0f;
    	double fdown=0.0f;
    	double fnetpay=0.0f;
		double fLD=GJJSParameter->BasicProH/GJJSParameter->BasicProRw;
    	ftop=GJJSParameter->fResHpt;
    	fdown=GJJSParameter->fResHpd;
    	fnetpay=ftop-fdown;
     	int i=0;
    	do
		{
	    	i++;
    		fe = pow(cos((2*i-1)*PI*ftop/2)-cos((2*i-1)*PI*fdown/2),2)*
				 K0(sqrt(ffsy(fu,GJJSParameter)+(2*i-1) * (2*i-1)*PI*PI/(4*fLD*fLD)))/((2*i-1)*(2*i-1));
    		fg=fg+fe;
		}while((fabs(fe) > 10e-10) && i < 300);
        fg=8*fg/(PI*PI*fnetpay*fnetpay*fu);

		double fpa=0.0f;//stephest算法
    	double fpb=0.0f;
    	double fpc=0.0f;
    	double fpd=0.0f;
    	fpa=fg;
    	fpb=fu*fpa+fS;
    	fpc=fu*(1.0+fu*fCD*(fu*fpa+fS));
    	fvalue=fpb/fpc;
	}
	else if(GJJSParameter->modeltype==11)//部分射孔底水--双重
	{
		double fe=0.0f;
       	double fg=0.0f;
    	double ftop=0.0f;
    	double fdown=0.0f;
    	double fnetpay=0.0f;
		double fLD=GJJSParameter->BasicProH/GJJSParameter->BasicProRw;
    	ftop=GJJSParameter->fResHpt;
    	fdown=GJJSParameter->fResHpd;
    	fnetpay=ftop-fdown;
     	int i=0;
    	do
		{
	    	i++;
    		fe = pow(cos((2*i-1)*PI*ftop/2)-cos((2*i-1)*PI*fdown/2),2)*
				 K0(sqrt(ffs(fu,GJJSParameter,fcde2s,fS)*fu+(2*i-1) * (2*i-1)*PI*PI/(4*fLD*fLD)))/((2*i-1)*(2*i-1));
    		fg=fg+fe;
		}while((fabs(fe) > 10e-10) && i < 300);
        fg=8*fg/(PI*PI*fnetpay*fnetpay*fu);

		double fpa=0.0f;//stephest算法
    	double fpb=0.0f;
    	double fpc=0.0f;
    	double fpd=0.0f;
    	fpa=fg;
    	fpb=fu*fpa+fS;
    	fpc=fu*(1.0+fu*fCD*(fu*fpa+fS));
    	fvalue=fpb/fpc;
	}
	else if(GJJSParameter->modeltype==7)//三重边水
	{		
		double fa=0.0f;
		double fb=0.0f;
		double fc=0.0f;
		double fd=0.0f;
		double fred=0.0f;
		fd=sqrt(ffs3(fu,GJJSParameter,fcde2s,fS)*fu/fcde2s);
		fred=GJJSParameter->fResBD1*exp(fS)/GJJSParameter->BasicProRw;
		fa=fd*(I1(fd)*K1(fred*fd)+K1(fd)*I0(fred*fd));
		fb=K0(fd)*I0(fred*fd)-I0(fd)*K0(fred*fd);
		fc=1.0/(fu*(fa/fb+fu));
		fvalue=fc;
	}
	else if(GJJSParameter->modeltype==8)//两区复合油藏
	{		
	    double cd=0.0f,fs=0.0f,fm=3.0,nam=0.0f;
    	double rd1=0.0f,rd2=0.0f;
    	
		cd=fCD;fs=fS;fm=GJJSParameter->fResM12;nam=GJJSParameter->fResEta12;
		rd1=GJJSParameter->fResRdi/GJJSParameter->BasicProRw;
		rd2=GJJSParameter->fResRdo/GJJSParameter->BasicProRw;

    	double a11=0.0f,a12=0.0f,a13=0.0f,a14=0.0f;
    	double a21=0.0f,a22=0.0f,a23=0.0f,a24=0.0f;
    	double a31=0.0f,a32=0.0f,a33=0.0f,a34=0.0f;
    	double a43=0.0f,a44=0.0f;

    	a11=cd*fu*fu*I0(sqrt(ffsy(fu,GJJSParameter)))-fu*sqrt(ffsy(fu,GJJSParameter))*(cd*fu*fs+1.0)*I1(sqrt(ffsy(fu,GJJSParameter)));
    	a12=cd*fu*fu*K0(sqrt(ffsy(fu,GJJSParameter)))+fu*sqrt(ffsy(fu,GJJSParameter))*(cd*fu*fs+1.0)*K1(sqrt(ffsy(fu,GJJSParameter)));
    	a22=K0(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a24=-K0(rd1*sqrt(nam*ffsy(fu,GJJSParameter)));
		a21=I0(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a31=fm*sqrt(ffsy(fu,GJJSParameter))*I1(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a32=-fm*sqrt(ffsy(fu,GJJSParameter))*K1(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a34=sqrt(nam*ffsy(fu,GJJSParameter))*K1(rd1*sqrt(nam*ffsy(fu,GJJSParameter)));
        //无限大
    	a23=0.0f;
    	a33=0.0f;
    	a43=0.0f;
    	a44=0.0f;

    	double b10=0.0f,b11=0.0f,b12=0.0f;
    	double b20=0.0f,b21=0.0f,b22=0.0f;

    	b10=-a21*1.0/a11;
    	b11=a22-a21*a12/a11;
       	b12=a24;
    	b20=-a31*1.0/a11;
    	b21=a32-a31*a12/a11;
      	b22=a34;

    	double a1=0.0f,a2=0.0f;
    	a1=(1.0/a11)-(a12*(b20*b12-b22*b10))/(a11*(b12*b21-b22*b11));
    	a2=(b20*b12-b22*b10)/(b12*b21-b22*b11);

    	double pwd=0.0f;
    	pwd=a1*(I0(sqrt(ffsy(fu,GJJSParameter)))-fs*sqrt(ffsy(fu,GJJSParameter))*I1(sqrt(ffsy(fu,GJJSParameter))))+
    		a2*(K0(sqrt(ffsy(fu,GJJSParameter)))+fs*sqrt(ffsy(fu,GJJSParameter))*K1(sqrt(ffsy(fu,GJJSParameter))));

    	fvalue=pwd;
	}
	else if(GJJSParameter->modeltype==9)//两区复合油藏封闭边界
	{		
	    double cd=0.0f,fs=0.0f,fm=3.0,nam=0.0f;
    	double rd1=0.0f,rd2=0.0f;
    	
		cd=fCD;fs=fS;fm=GJJSParameter->fResM12;nam=GJJSParameter->fResEta12;
		rd1=GJJSParameter->fResRdi/GJJSParameter->BasicProRw;
		rd2=GJJSParameter->fResRdo/GJJSParameter->BasicProRw;

    	double a11=0.0f,a12=0.0f,a13=0.0f,a14=0.0f;
    	double a21=0.0f,a22=0.0f,a23=0.0f,a24=0.0f;
    	double a31=0.0f,a32=0.0f,a33=0.0f,a34=0.0f;
    	double a43=0.0f,a44=0.0f;

    	a11=cd*fu*fu*I0(sqrt(ffsy(fu,GJJSParameter)))-fu*sqrt(ffsy(fu,GJJSParameter))*(cd*fu*fs+1.0)*I1(sqrt(ffsy(fu,GJJSParameter)));
    	a12=cd*fu*fu*K0(sqrt(ffsy(fu,GJJSParameter)))+fu*sqrt(ffsy(fu,GJJSParameter))*(cd*fu*fs+1.0)*K1(sqrt(ffsy(fu,GJJSParameter)));
    	a21=I0(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a22=K0(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a24=-K0(rd1*sqrt(nam*ffsy(fu,GJJSParameter)));
    	a31=fm*sqrt(ffsy(fu,GJJSParameter))*I1(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a32=-fm*sqrt(ffsy(fu,GJJSParameter))*K1(rd1*sqrt(ffsy(fu,GJJSParameter)));
    	a34=sqrt(nam*ffsy(fu,GJJSParameter))*K1(rd1*sqrt(nam*ffsy(fu,GJJSParameter)));
        //封闭
       	a23=-I0(rd1*sqrt(nam*ffsy(fu,GJJSParameter)));
       	a33=-sqrt(nam*ffsy(fu,GJJSParameter))*I1(rd1*sqrt(nam*ffsy(fu,GJJSParameter)));
       	a43=sqrt(nam*ffsy(fu,GJJSParameter))*I1(rd2*sqrt(nam*ffsy(fu,GJJSParameter)));
       	a44=-sqrt(nam*ffsy(fu,GJJSParameter))*K1(rd2*sqrt(nam*ffsy(fu,GJJSParameter)));

      

    	double b10=0.0f,b11=0.0f,b12=0.0f;
    	double b20=0.0f,b21=0.0f,b22=0.0f;

    	b10=-a21*1.0/a11;
    	b11=a22-a21*a12/a11;
       	b12=a24-a23*a44/a43;
       	b20=-a31*1.0/a11;
    	b21=a32-a31*a12/a11;
       	b22=a34-a33*a44/a43;
     	
    	double a1=0.0f,a2=0.0f;
    	a1=(1.0/a11)-(a12*(b20*b12-b22*b10))/(a11*(b12*b21-b22*b11));
    	a2=(b20*b12-b22*b10)/(b12*b21-b22*b11);

    	double pwd=0.0f;
    	pwd=a1*(I0(sqrt(ffsy(fu,GJJSParameter)))-fs*sqrt(ffsy(fu,GJJSParameter))*I1(sqrt(ffsy(fu,GJJSParameter))))+
    		a2*(K0(sqrt(ffsy(fu,GJJSParameter)))+fs*sqrt(ffsy(fu,GJJSParameter))*K1(sqrt(ffsy(fu,GJJSParameter))));

    	fvalue=pwd;
	}
	else if(GJJSParameter->modeltype==10)//顶底封闭部分射孔
	{
		double fa=0.0f;
    	double fb=0.0f;
    	double fc=0.0f;
    	double fd=0.0f;
    	double fe=0.0f;
    	double ff=0.0f;
    	double fg=0.0f;
    	fa=K0(sqrt(fu))/(fu);
    	double ftop=0.0f;
    	double fdown=0.0f;
    	double fnetpay=0.0f;
		double fLD=GJJSParameter->BasicProH/GJJSParameter->BasicProRw;
    	ftop=GJJSParameter->fResHpt;
    	fdown=GJJSParameter->fResHpd;
       	fnetpay=ftop-fdown;
      	int i=0;
	    do
		{
    		i++;
    		fe=pow(sin(i*PI*ftop)-sin(i*PI*fdown),2)*K0(sqrt(ffsy(fu,GJJSParameter)+i*i*PI*PI/(fLD*fLD)))/(i*i);
    		fg=fg+fe;
		}while(fabs(fe) > MinDoubleNum && i<3000);
        fg=2*fg/(PI*PI*fnetpay*fnetpay*fu);

    	double fpa=0.0f;//stephest算法
    	double fpb=0.0f;
    	double fpc=0.0f;
    	double fpd=0.0f;
    	fpa=fa+fg;
    	fpb=fu*fpa+fS;
    	fpc=fu*(1.0+fu*fCD*(fu*fpa+fS));
    	fpd=fpb/fpc;
		fvalue=fpd;
	}
	if(GJJSParameter->modeltype==12)//均质启动压力梯度
	{
		double fa=0.0f;
    	double fb=0.0f;
    	double fc=0.0f;
    	double fd=0.0f;
        double fNambd=0.0f;

		fNambd=((GJJSParameter->fResK)*(GJJSParameter->BasicProH)*(GJJSParameter->BasicProRw)*(GJJSParameter->fQidongyali))/
			(1.842*(GJJSParameter->BasicProQ)*(GJJSParameter->BasicProBo)*(GJJSParameter->BasicProMiuo));
	
		fNambd=fNambd*exp(-fS);

		fa=fNambd+PI*fNambd*I1(sqrt(ffsy(fu,GJJSParameter)/fcde2s))/2.0-PI*fNambd*fu*I0(sqrt(ffsy(fu,GJJSParameter)/fcde2s))/(2.0*sqrt(ffsy(fu,GJJSParameter)/fcde2s));

		fb=(1+fa)*K0(sqrt(ffsy(fu,GJJSParameter)/fcde2s))/(fu*(sqrt(ffsy(fu,GJJSParameter)/fcde2s)*K1(sqrt(ffsy(fu,GJJSParameter)/fcde2s))+fu*K0(sqrt(ffsy(fu,GJJSParameter)/fcde2s))));

		fc=PI*fNambd*I0(sqrt(ffsy(fu,GJJSParameter)/fcde2s))/(2.0*fu*sqrt(ffsy(fu,GJJSParameter)/fcde2s));

		fd=fb+fc;

    	fvalue=fd;
	}
	return fvalue;	      
}

// 边界叠加处理
double __stdcall BoundaryLap(double fu,
							 double fCD,
							 double fS,
							 GJJS_Parameter* GJJSParameter,
							 double fd1,
							 double fd2, 
							 double fd3,
							 double fd4,
							 double fvnum)
{
	double numpwd=0.0f;
	double flinelong;
	int i;
	if(GJJSParameter->mboundary==1)//直线断层
	{
	    flinelong=2.0*fd1;
	    numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
	}
	if(GJJSParameter->mboundary==2)//60断层
	{
		flinelong=2.0*fd1;
	    numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
		flinelong=2.0*fd2;
	    numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
		flinelong=2.0*sqrt(fd1*fd1+fd2*fd2+fd1*fd2);
	    numpwd+=2.0*fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
		flinelong=2.0*2.0*sqrt(fd1*fd1+fd2*fd2+fd1*fd2)/sqrt(3.0);
	    numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
	}
	if(GJJSParameter->mboundary==3)//90断层
	{
		flinelong=2.0*fd1;
	    numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
		flinelong=2.0*fd2;
	   	numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
		flinelong=2.0*sqrt(fd1*fd1+fd2*fd2);
	   	numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
	 }
	 if(GJJSParameter->mboundary==4)//120断层
	 {
		flinelong=2.0*fd1;
	   	numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
		flinelong=2.0*fd2;
	    numpwd+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
	 }
	 if(GJJSParameter->mboundary==5)//180断层
	 {
		i=1;
		double fa=0.0f;
		double fb=0.0f;
		double fc=0.0f;
	    do
		{
			fa=0.0f;
			flinelong=2.0*i*(fd1+fd2);
			fa=2.0*fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=2.0*i*fd1+2.0*(i-1)*fd2;
			fa+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=2.0*(i-1)*fd1+2.0*i*fd2;
			fa+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			numpwd+=fa;
			i++;
		}while(fabs(fa/numpwd)>0.0000001 && i>20);
	 }
	 if(GJJSParameter->mboundary==6)//厂型断层
	 {
		i=1;
		double fa=0.0f;
		double fb=0.0f;
		double fc=0.0f;
	    do
		{
			fa=0.0f;
			flinelong=2.0*i*(fd1+fd2);
			fa=2.0*fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=sqrt(fd3*fd3+4.0*i*(fd1+fd2)*i*(fd1+fd2));
			fa=2.0*fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=2.0*i*fd1+2.0*(i-1)*fd2;
			fa+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=sqrt(fd3*fd3+(2.0*i*fd1+2.0*(i-1)*fd2)*(2.0*i*fd1+2.0*(i-1)*fd2));
			fa+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=2.0*(i-1)*fd1+2.0*i*fd2;
			fa+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			flinelong=sqrt(fd3*fd3+(2.0*(i-1)*fd1+2.0*i*fd2)*(2.0*(i-1)*fd1+2.0*i*fd2));
			fa+=fvnum*laplas(fu,fCD,fS,GJJSParameter,flinelong);
			numpwd+=fa;
			i++;
		}while(fabs(fa/numpwd)>0.0000001 && i>20);
	 }
  	return numpwd;

}
// 实空间压力计算公式程序----------
int __stdcall Calfinit(double* pfPwd, double* pfTD, int& number, GJJS_Parameter* GJJSParameter)
{
	//无因次化
	double fa=GJJSParameter->InputDiaPP[0];
	double fb=GJJSParameter->InputDiaTP[0];
	double InputDiaNum=GJJSParameter->InputDiaNum;
	double ftp=GJJSParameter->ftp;
	double fS=GJJSParameter->fResS;
	double mintime=0.0f;
	double maxtime=0.0f;
	int    TheoryNum=0;
	double fCD=0.0f;
	int i=0;
    for(i=0;i<=InputDiaNum-2;i++)
	{
		if(GJJSParameter->innerboundary==0)
		{
        	pfTD[i]=GJJSParameter->InputDiaTP[i+1]-fb;
        	pfTD[i]=0.0036*(GJJSParameter->fResK)*pfTD[i]/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
        		(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
    	}
    	if(GJJSParameter->innerboundary==4 ||GJJSParameter->innerboundary==5)//表皮叠加公式
		{
         	pfTD[i]=GJJSParameter->InputDiaTP[i+1]-fb;
        	pfTD[i]=0.0036*(GJJSParameter->fResK)*pfTD[i]/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
    		(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
    	}
    	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
		{
	    	if(fS<0.0)
			{
		    	GJJSParameter->fResS=0.0f;
	    		fS=GJJSParameter->fResS;
			}
	    	pfTD[i]=GJJSParameter->InputDiaTP[i+1]-fb;
        	pfTD[i]=0.0036*(GJJSParameter->fResK)*pfTD[i]/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
        		(GJJSParameter->BasicProCt)*(GJJSParameter->fResXf)*(GJJSParameter->fResXf));
        }

	}


	if(GJJSParameter->innerboundary==0)
	{
    	fCD=(GJJSParameter->fResC)/
	    	(2.0*3.1415926*(GJJSParameter->BasicProFai)*(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProH)*
    		(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
	    ftp=(1.0/fCD)*0.0036*(GJJSParameter->fResK)*ftp/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
	    	(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
	}
	if(GJJSParameter->innerboundary==4 ||GJJSParameter->innerboundary==5)//表皮叠加公式
	{
    	fCD=(GJJSParameter->fResC)/
	    	(2.0*3.1415926*(GJJSParameter->BasicProFai)*(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProH)*
    		(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
 	    ftp=0.0036*(GJJSParameter->fResK)*ftp/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
	    	(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
	}
	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
	{
		if(fS<0.0)
		{
			GJJSParameter->fResS=0.0f;
			fS=GJJSParameter->fResS;
		}
		fCD=(GJJSParameter->fResC)/
	    	(2.0*3.1415926*(GJJSParameter->BasicProFai)*(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProH)*
    		(GJJSParameter->fResXf)*(GJJSParameter->fResXf));
        ftp=0.0036*(GJJSParameter->fResK)*ftp/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
	    	(GJJSParameter->BasicProCt)*(GJJSParameter->fResXf)*(GJJSParameter->fResXf));
	}
	////////////////////////////////////////////////////////////
	//计算无因次理论压力（压恢）
    int j;
	int pnumber=0;
    double ff=0.07f;
	int nmin;
	int nmax;
	double timestep;
	double fd1=0.0f,fd2=0.0f;
	double fd3=0.0f,fd4=0.0f;
	double fretype=0.0f;
	double flinelong=0.0f;
	if(GJJSParameter->innerboundary==0 )//无因此化
	{
		fretype=exp(fS)/GJJSParameter->BasicProRw;
	}
	if(GJJSParameter->innerboundary==4 || GJJSParameter->innerboundary==5)//无因此化
	{
		fretype=1.0/GJJSParameter->BasicProRw;
	}
	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
	{
		fretype=1.0/GJJSParameter->fResXf;
	}

	if(GJJSParameter->mboundary==1 || GJJSParameter->mboundary==8)//一条边界
	{
	    fd1=GJJSParameter->fResBD1*fretype;
	}
	if(GJJSParameter->mboundary==2 || GJJSParameter->mboundary==3
			|| GJJSParameter->mboundary==4 || GJJSParameter->mboundary==5)//2条边界
	{
		fd1=GJJSParameter->fResBD1*fretype;
		fd2=GJJSParameter->fResBD2*fretype;
	}
	if(GJJSParameter->mboundary==6)//3条边界
	{
		fd1=GJJSParameter->fResBD1*fretype;
		fd2=GJJSParameter->fResBD2*fretype;
		fd3=GJJSParameter->fResBD3*fretype;
	}
	if(GJJSParameter->mboundary==7)//4条边界
	{
		fd1=GJJSParameter->fResBD1*fretype;
		fd2=GJJSParameter->fResBD2*fretype;
		fd3=GJJSParameter->fResBD3*fretype;
		fd4=GJJSParameter->fResBD4*fretype;
	}

	if(GJJSParameter->innerboundary==0 )//有效井井
	{
    	for(i=0;i<=InputDiaNum-2;i++)
		{
	 	          int num;
		          num=0;
		          double numpwd;
		          double s[8];
		          double nump;
		          nump=0;
		          numpwd=0.0f;
	              for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/(pfTD[i]*(1.0/fCD));
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
	    	      pfPwd[i]=log(2)*numpwd/(pfTD[i]*(1.0/fCD));
 
		    	  numpwd=0.0f;
		    	  for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/((pfTD[i]*(1.0/fCD))+ftp);
		    	       numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
		          pfPwd[i]=pfPwd[i]-log(2)*numpwd/((pfTD[i]*(1.0/fCD))+ftp);

		    	  numpwd=0.0f;
		    	  for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/ftp;
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
		          pfPwd[i]=pfPwd[i]+log(2)*numpwd/ftp;	

		  }
		
	}
	if(GJJSParameter->innerboundary==4 || GJJSParameter->innerboundary==5)//井
	{
    	for(i=0;i<=InputDiaNum-2;i++)
		{
	 	          int num;
		          num=0;
		          double numpwd;
		          double s[8];
		          double nump;
		          nump=0;
		          numpwd=0.0f;
	              for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/(pfTD[i]);
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
	    	      pfPwd[i]=log(2)*numpwd/(pfTD[i]);
 
		    	  numpwd=0.0f;
		    	  for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/((pfTD[i])+ftp);
		    	       numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
		          pfPwd[i]=pfPwd[i]-log(2)*numpwd/((pfTD[i])+ftp);

		    	  numpwd=0.0f;
		    	  for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/ftp;
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
		          pfPwd[i]=pfPwd[i]+log(2)*numpwd/ftp;	

		    	  
		}
	}
	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)//裂缝
	{
		for(i=0;i<=InputDiaNum-2;i++)
		{
	     	      int num;
	    	      num=0;
		          double numpwd;
		          double s[8];
	    	      double nump;
	    	      nump=0;
		          numpwd=0.0f;
	              for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/(pfTD[i]);
	    		       numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
					   
				  }
		          pfPwd[i]=log(2)*numpwd/(pfTD[i]);

		    	  numpwd=0.0f;
			      for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/((pfTD[i])+ftp);
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
		          pfPwd[i]=pfPwd[i]-log(2)*numpwd/((pfTD[i])+ftp);

		    	  numpwd=0.0f;
		    	  for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/ftp;
		    	       numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
                       if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
		          pfPwd[i]=pfPwd[i]+log(2)*numpwd/ftp;	
			 
		}
	}
	
	//有因次化
	for(i=0;i<InputDiaNum-2;i++)
	{
		pfPwd[i]=pfPwd[i]*1.842*(GJJSParameter->BasicProQ)*(GJJSParameter->BasicProBo)*(GJJSParameter->BasicProMiuo)/
		    	((GJJSParameter->fResK)*(GJJSParameter->BasicProH));
		if(GJJSParameter->innerboundary==0 || GJJSParameter->innerboundary==4 || GJJSParameter->innerboundary==5)
		{
	    	pfTD[i]=pfTD[i]*((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*(GJJSParameter->BasicProCt)*
    		    	(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw))/(0.0036*(GJJSParameter->fResK));
		}
		if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
		{
    		pfTD[i]=pfTD[i]*((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*(GJJSParameter->BasicProCt)*
       		    	(GJJSParameter->fResXf)*(GJJSParameter->fResXf))/(0.0036*(GJJSParameter->fResK));
		}
	}
	number=InputDiaNum-2;
	return 0;
}


int __stdcall CalAutofix(GJJS_Parameter* GJJSParameter)
{
	double fK=0;
	double fC=0;
	double fS=0;
	double fNumda=0;
	double fWomiga=0;
	double ftemp=0;
	int  Paranum=0;
	int i =0;
	double fKstep[1000];
	double fCstep[1000];
	double fSstep[1000];
	double fPara[100][100];
	fK=GJJSParameter->fResK;
    ftemp=fK/10;
	for(i=0;i<20;i++)
	{
		fPara[0][i]=fK/10.0+i*ftemp;
	}
	fC=GJJSParameter->fResC;
    ftemp=fC/10;
	for(i=0;i<20;i++)
	{
		fPara[1][i]=fC/10.0+i*ftemp;
	}
    fS=GJJSParameter->fResS;
    ftemp=fS/10;
	for(i=0;i<20;i++)
	{
		fPara[2][i]=fS/10.0+i*ftemp;
	}
    Paranum=3;
	//双重介质
	if(GJJSParameter->modeltype==2)
	{
		Paranum=Paranum+2;
		fNumda=GJJSParameter->fResNamuda;
		ftemp=fNumda/10;
     	for(i=0;i<20;i++)
		{
    		fPara[Paranum-2][i]=fS/10.0+i*ftemp;
		}

		fNumda=GJJSParameter->fResWomiga;
		ftemp=fWomiga/10;
     	for(i=0;i<20;i++)
		{
    		fPara[Paranum-1][i]=fS/10.0+i*ftemp;
		}


	}

	if(GJJSParameter->modeltype==3)
	{
		Paranum=Paranum+2;
		fNumda=GJJSParameter->fResNamuda;
		ftemp=fNumda/10;
     	for(i=0;i<20;i++)
		{
    		fPara[Paranum-2][i]=fS/10.0+i*ftemp;
		}

		fNumda=GJJSParameter->fResWomiga;
		ftemp=fWomiga/10;
     	for(i=0;i<20;i++)
		{
    		fPara[Paranum-1][i]=fS/10.0+i*ftemp;
		}


	}


    return 0;
}

// 井筒压降计算程序
double  H_Hagedorn_Brown(double Panbie, double Qg_sc, double Qw_sc, double Qg_ing,double  D_youg_nj,double  Tjingkou_sc,double  Tjingdi_sc,double  Hyouc_sc, 
						 double bz_scg,double  Bz_ing,double  Ppr_sc,double  Ttr_sc,double s, double Md_w_sc,double  Czd_xd, double P_qid,double  Bmyz_w,double Md_o_sc,double Qo_sc)
{
	double P_pj=0.0f; double T_pj=0.0f; double Bzg=0.0f; double PPr=0.0f; double Tpr=0.0f; double Khd=0.0f; double P1=0.0f;
	double T1=0.0f; double  P2=0.0f; double T2=0.0f; double Nd_o=0.0f; double Md_o=0.0f; double Pb=0.0f; double  Rso=0.0f;
	double Zg=0.0f; double  Rso_1=0.0f; double Rso_2=0.0f; double Zg_1=0.0f; double Zg_2=0.0f; double e = 2.718282f; double  Pi = 3.14159265359f;
	double g = 9.81f;
	int i=0 ; int K_zus=0;
   	Md_o = Md_o_sc;
	double fw=0.0f;
	fw = Qw_sc / (Qw_sc + Qo_sc);
	double  DeltaH0=0.0f;
	DeltaH0 = 50; //M迭代高度;
	int i1=0;
    double i2=0.0f;
	int Ni=0;
	i1 = int(Hyouc_sc / DeltaH0);
	i2 = Hyouc_sc / DeltaH0;
	Ni = i1 + 1;
	if(i1 = i2 )
	{
		 Ni = i1 + 1;
	}
	K_zus = int(Hyouc_sc / DeltaH0) + 1;

	double H[1000];
	double P[1000];
	PPr = Ppr_sc;
	Tpr = Ttr_sc;
    Qg_sc = Qg_sc;
	Qg_ing = Qg_ing;// '二者相加！！！
	int I=0;
	int J=0;
    double DeltaH_1=0.0f;
	double H_sd=0.0f;
	double DeltaH_i=0.0f;
	double Qg=0.0f;
	double DeltaT=0.0f;
	double DeltaP=0.0f;
	double B_o=0.0f;
	double Bmyz_o=0.0f;//如何赋值


	if( Panbie = 1)
	{
		H_sd = 0; I = 0; T1 = Tjingkou_sc;
        DeltaH_1 = DeltaH0; DeltaH_i = DeltaH0;
        H[0]= 0; 
		P[0]= P_qid; 
		J = 0;
        Qg = Qg_sc + Qg_ing;
        Bzg = (Qg_sc * bz_scg + Qg_ing * Bz_ing) / Qg;
	}
	else
	{
		H_sd = Hyouc_sc;  T1 = Tjingdi_sc;
       if(i1 = i2)
	   {
		    DeltaH_i = 100;
			DeltaH_1 = 100;
	   }
	   else
	   {
		   DeltaH_i = 100; 
		   DeltaH_1 = Hyouc_sc - i1 * DeltaH_i;
	   }
		H[Ni] = H_sd; P[Ni] = P_qid;J = Ni;
        Qg = Qg_sc; Bzg = bz_scg;
	}
	double A_youg=0.0f;
	double DeltaH=0.0f;
	double B_w=0.0f;
	double Bg=0.0f;
	double Rsw=0.0f;
	double Nd_g=0.0f;
	double Ndw=0.0f;
	double Md_g=0.0f;
	double Md_w=0.0f;
	double Mdl=0.0f;
	double Bmyzl=0.0f;
	double Qsl=0.0f;
	double Qsg=0.0f;
	double Vsg_1=0.0f;
	double Vsl=0.0f;
	double Vsg=0.0f;
	double Vm=0.0f;
	double Qm=0.0f;
	double Lamda=0.0f;
	double NL=0.0f;
	double Nvl=0.0f;
	double nvg=0.0f;
	double Ndl=0.0f;
	double L1=0.0f;
	double L2=0.0f;
	double Vs=0.0f;
	double Hg=0.0f;
	double Md_m=0.0f;
	double Nrem=0.0f;
	double F_m=0.0f;
	double Nd=0.0f;
	double Fai=0.0f;
	double CNl=0.0f;
	double Fai_1=0.0f;
	double  Fai_xz=0.0f;
	double Xie=0.0f;
	double H_l=0.0f;
	double Hl_Fai=0.0f;
	double Nd_m=0.0f;
	double Mt=0.0f;
	double F_m_1=0.0f;
	double DeltaP2_1=0.0f;
	double DeltaP2=0.0f;
	double DP=0.0f;
	int panbiekk=0;
	A_youg =( Pi / 4 )*(D_youg_nj / 1000)*(D_youg_nj / 1000);
	P1 = P_qid;
    DeltaT = (-Tjingkou_sc + Tjingdi_sc) / Hyouc_sc;
    DeltaH = DeltaH_1;
	do
	{
		DeltaP = 1;
        H_sd = H_sd + DeltaH * Panbie;
        if( H_sd > Hyouc_sc )
		{
			DeltaH = Hyouc_sc - H_sd + DeltaH;
            H_sd = Hyouc_sc;
		}
		if(H_sd <= 0)
		{
			DeltaH = H_sd + DeltaH;
            H_sd = 0;
		}
		do
		{
			T2 = T1 + DeltaT * DeltaH * Panbie; P2 = P_qid + DeltaP * Panbie;
            T_pj = T1 + DeltaT * DeltaH * Panbie / 2;
            P_pj = P1 + DeltaP * Panbie / 2;
            Zg = Z_yasyz(P_pj, T_pj, Bzg);
            B_w = B_water(T_pj, P_pj) ;//'T_pj----C
            
            Bg = 3.458 * pow(10,-4) * Zg * (T_pj + 273) / P_pj;
            Khd = s;
            Rsw = Rsg_w(T_pj, P_pj, Khd);// 'S--矿化度
			if(fw < 1)
			{
				Rso = Rsg_o(T_pj, P_pj, Md_o, Bzg);
                B_o = B_oil(T_pj, P_pj, Bzg, Md_o, Rso);
                Nd_o = Nd_oil(T_pj, P_pj, Pb, Md_o, Bzg, Rso);
			}
			Nd_g = Nd_gas(Bzg, P_pj, T_pj, Zg);
            Ndw = Nd_water(T_pj);
			Nd_g = 0.02;
            Ndw = 0.789;
            
            Md_g = Md_gas(Bzg, P_pj, T_pj, Zg);  //'Kg/m^3
            Md_w = Md_w_sc;
			if(fw < 1)
			{
				Md_o = 1;
			}
            Mdl = Md_w * fw + Md_o * (1 - fw);
            Md_g = 1.205 * Bzg * (P_pj / 0.101325) * (293 / (T_pj + 273)) / Zg;
            Ndl = Ndw * fw + Nd_o * (1 - fw);
            
            Bmyzl = Bmyz_w * fw + Bmyz_o * (1 - fw);
                
            Qsl = Qw_sc * B_w + Qo_sc * B_o;
            Qsg = 4.084 * pow(10,-9) * (Qg) * (T_pj + 273) * Zg / (P_pj);
            Vsg_1 = Qsg / A_youg;
            Vsl = Qsl / (86400 * A_youg);
            Vsg = Vsl * (Qg / (Qw_sc + Qo_sc)) * (0.101325 / P_pj) * ((T_pj + 273) / 293) * Zg;
            Vm = Vsl + Vsg;
            Qm = Qsg + Qsl / 86400;
            Lamda = Qsg / Qm;
                        
            NL = 0.3147 * Ndl * pow((Mdl * pow(Bmyzl , 3)),-0.25) ;// '液相粘度准数
            Nvl = 3.178 * Vsl *pow( (Mdl / Bmyzl), 0.25);//  '液相速度准数
            nvg = 3.178 * Vsg *pow( (Mdl / Bmyzl), 0.25);// '气流速度准数

			L1 = 1.071 - 0.7277 * pow(Vm,2) / (D_youg_nj / 1000);// 'D_youg_nj__单位？？？？
            if( L1 < 0.13)
			{
				 L1 = 0.13;
			}
			L2 = Vsg / (Vsg + Vsl);
			if( L1 > L2)
			{
				Vs = 0.3048 * 0.8;
                Hg = 0.5 * (1 + Qm / (Vs * A_youg) -pow( (pow((1 + Qm / (Vs * A_youg)) , 2) - 4 * Qsg / (Vs * A_youg)) ,0.5));
                Md_m = (1 - Hg) * Mdl + Hg * Md_g;
                Nrem = 1000 * Mdl * D_youg_nj / 1000 * Vsl / (Ndl * (1 - Hg));
                F_m = fm(Czd_xd, Nrem, Qg);
			}
			else
			{
				Nd = 99.405 * (D_youg_nj / 1000) *pow( (Mdl / Bmyzl) , 0.5);// '管径准数
                Fai = Nvl * pow((P_pj / 0.101325) , 0.1) * CNl * pow(nvg, (-0.575) )/ Nd;//shifou youwent //cnl
                Fai_1 = Fai * 1000;
//'                Nname$ = App.Path & "\H&B_HlFai.dat"
//'                Hl_Fai = CHazhi(Nname$, Fai_1)  '？？？？
                
               Fai_xz = nvg *pow( NL, 0.38 )/pow( Nd , 2.14);
               if(Fai_xz <= 0.01 )
			   {
				   Xie=1.0f;
			   }
			   else
			   {
				   Xie=1.0f;///youwenti
				//   Nname$ = App.Path & "\H&B_Fai.dat"
                 //   Xie = CHazhi(Nname$, Fai)
			   }

			    H_l = Hl_Fai * Xie;//   'Xie=1///Hl_Fai
                Md_m = Mdl * H_l + Md_g * (1 - H_l);
                Nd_m =pow( Ndl, H_l) *pow( Nd_g ,(1 - H_l));
                Mt = (Md_w_sc * fw + Md_o_sc * (1 - fw)) + 1.025 * Bzg * Qg / (Qw_sc + Qo_sc);
               
                Nrem = 1.474 *pow( 10 , -2 )* (Qw_sc * B_w * fw + Qo_sc * B_o * (1 - fw)) * Mt / (D_youg_nj / 1000 * Nd_m);
             //   '杨川东  < Jain 方法
                F_m_1 = fm(Czd_xd, Nrem, Qg);
         //       '刘建议——Jain 方法
                F_m =pow( (1.14 - 2 * log10(Czd_xd + 21.25 / pow(Nrem,0.9))) , -2 );//'刘建议????


			}
			////////

			DeltaP2_1 =pow( 10 ,-6) * (Md_m * g + F_m_1 * pow(Qsl, 2) * pow(Mt , 2) / (9.21 * pow(10 , 9) * Md_m * pow((D_youg_nj / 1000) , 5))) * DeltaH;
            DeltaP2 = pow(10,-6) * (Md_m * g + F_m * pow(Qsl ,2) *pow( Mt , 2) / (9.21 * pow(10 , 9) * Md_m *pow( (D_youg_nj / 1000) ,5))) * DeltaH;
            
            DP = fabs(DeltaP - DeltaP2) / DeltaP2;
            if(DP <= 0.01 )
			{
				T1 = T2; 
				P1 = P2;
			}
			else
			{
				DeltaP = DeltaP2;
			}
			
			
		}while(DP > 0.01  );
		DeltaH = DeltaH_i;
        P2 = P_qid + DeltaP * Panbie;
		P_qid = P2;
       if(Panbie = 1 )
	   {
		   J = J + 1;
	   }
	   else
	   {
		   J = J - 1;
	   }
		
		H[J] = H_sd; 
		P[J] = P2;
       
	}while(H_sd < Hyouc_sc);
	return P2;
}


double Z_yasyz(double P_pj ,double T_pj,double Bz_g) //偏差系数计算//ok
{
	//以下计算P_平均，T_平均下的各物理参数
    //If T_pj > 273 Then T_pj = T_pj - 273.15

	double ppc=0.0f;
	double tpc=0.0f;
	double Pr=0.0f;
	double tr=0.0f;
    Pr = P_pj * pow(10,6);
    tr = T_pj;
    //在未知临界温度Tpc、临界压力Ppc条件下  //////
    //-------------对凝析气
    //-------------对干气
    if (Bz_g >= 0.7)
	{
        tpc = 92.2 + 176.6 * Bz_g; //K
        ppc = (4.881 - 0.3861 * Bz_g) * pow(10 , 6); //Pa
	}
    else
	{
        tpc = 92.2 + 176.7 * Bz_g; //K
        ppc = (4.778 - 0.248 * Bz_g) * pow(10 ,6);//Pa
	}
   
	///////////////
	double Tpr=0.0f;
	double PPr=0.0f;
	double Zg=0.0f;
	double Zg1=0.0f;
	double dspr=0.0f;
	double A=0.0f;
	double B=0.0f;
	double C=0.0f;
	double D=0.0f;
	double tr2=0.0f;
	double tr3=0.0f;
	double Dr=0.0f;
	double f_f=0.0f;
	double Fp=0.0f;
	double Bg=0.0f;
	double Bg1=0.0f;


	 Tpr = (T_pj + 273) / tpc;
	 PPr = Pr / ppc;
	 if (P_pj < 35)//Cranmer方法
	 {
		 Zg = 1;
	     do
		 {        
			 Zg1 = Zg;
              dspr = 0.27 * PPr / (Zg * Tpr);
              Zg = 1 + (0.31506 - 1.0467 / Tpr - 0.5783 / pow(Tpr, 3)) * dspr + (0.5353 - 0.6123 / pow(Tpr,2)) * pow(dspr, 2) + 0.6815 * pow(dspr , 2 )/ pow(Tpr,3);
		 }while (fabs(Zg - Zg1) > 0.001);
	 }
     
      else //Pr > 35 Then 'Hall-Yarbough方法
	  {
            tr = Tpr;
            tr2 = pow(tr,2);
			tr3 = pow(tr,3);
            A = -(14.76 / tr - 9.76 / tr2 + 4.58 / tr3);
            B = 90.7 / tr - 242.2 / tr2 + 42.4 / tr3;
            C = 1.18 + 2.82 / tr;
            D = 0.06125 * PPr * exp(-1.2 * pow((1 - 1 / tr),2)) / tr;
            Dr = D;
           do
		   {
                f_f = (1 + Dr +pow( Dr, 2) - pow(Dr ,3)) /pow( (1 - Dr), 3);
                f_f = f_f + A * Dr + B * pow(Dr,C) - D / Dr;
                Fp = (4 + 4 * Dr - 2 * pow(Dr, 2)) / pow((1 - Dr), 4);
                Fp = Fp + A + B * C *pow( Dr ,(C - 1)) + D / pow(Dr, 2);
                Dr = Dr - f_f / Fp;
		   }while(fabs(f_f / Fp) > 0.0001);
           Zg = D / Dr;
	  }
        Bg = 3.458 *pow( 10 ,-4) * Zg * (tr + 273) / (Pr * pow(10 ,-6));// 'm^3/m^3(标)   T_平均——K；P_平均——MPa
        Bg1 = 0.101 * pow(10 ,6) * Zg * (T_pj + 273) / (293 * P_pj);
        return Zg;

}

double B_water(double T,double P)//水体积系数
{
	double Ceta=0.0f;
    Ceta = 1.8 * T + 32;  //T——C
    double a1=0.0f;
	double a2=0.0f;
	double a3=0.0f;
    a1 = 0.9947 + 5.8 * pow(10,-6)* Ceta + 1.02 * pow(10,-6) *pow( Ceta , 2);
    a2 = -4.228 * pow(10,-6) + 1.8376 * pow(10,-8) * Ceta - 6.77 * pow(10,-11) *pow( Ceta, 2);
    a3 = 1.3 * pow(10,-10) - 1.3855 * pow(10,-12) * Ceta + 4.285 *pow( 10 ,-15) * pow(Ceta , 2);
	return a1 + a2 * (145.03 * P) + a3 *pow( 145.03 * P , 2);
}
double Rsg_w(double T,double P,double s)//溶解气水比，s,矿化度%
{
    double Ceta=0.0f;
    Ceta = 1.8 * T + 32;  //T——C
	double A=0.0f;
	double B=0.0f;
	double C=0.0f;
	double SC=0.0f;
	double Raw=0.0f;
	double Rsw=0.0f;
    A = 2.12 + 3.45 * pow(10,-3) * Ceta - 3.59 * pow(10,-15) * pow(Ceta, 2);
    B = 0.0107 - 5.26 * pow(10,-5) * Ceta + 1.48 * pow(10,-7) * pow(Ceta,2);
    C = -8.75 * pow(10,-7) + 3.9 * pow(10,-9) - 1.02 * pow(10,-11) * pow(Ceta,2);
    SC = 1 - (0.0753 - 0.000173 * Ceta) * s;
    Rsw = (A + B * (145.03 * P) + C * pow((145.03 * P), 2)) / 5.615;
    return  Rsw * SC;
}
double Rsg_o(double T,double P,double Md_o,double Bz_g)//Glaso法(1980) 溶解汽油比温度 压力 密度 气比重 有问题
{
	double Ceta=0.0f;
    Ceta = 1.8 * T + 32;  //T——C
    double D=0.01;
	double Ezj=0.01;
	double po=0.01;
	D = 141.5 / (Md_o / 1000) - 131.5;
    Ezj = 2.8869 - pow((14.1811 - 3.3093 * log10(145.03 * P) ) ,0.5);
    po =pow( 10, Ezj);
    return  Bz_g / 5.615 *pow( ((pow(D , 0.989) /pow( Ceta, 0.172)) * po) , 1.2255);

}

double B_oil(double T,double p,double Bz_g,double Md_o,double Rs )//Glaso法(1980)油体积系数
{
    double Ceta=0.0f;
    Ceta = 1.8 * T + 32;  //T——C
    double Bz_o=0.01;
	double aa=0.01;
	double bb=0.01;
    Bz_o = Md_o / 1000;
    aa = 5.615 * Rs *pow( (Bz_g / Bz_o) , 0.526 )+ 0.968 * Ceta;
    bb = -6.58511 + 2.91329 * log10(aa) - 0.27683 * pow((log10(aa)) , 2);
    return 1.0 + pow(10 , bb);
}

double Nd_oil(double T,double P,double Pb,double Md_o,double Bzg ,double Rs)//油粘度
{
    double Ceta=0.01;
    Ceta = 1.8 * T + 32;  //T——C
	double D=0.01;
	double aa=0.01;
	double Nd_od=0.01;
    D = 141.5 / (Md_o / 1000) - 131.5;
    aa = 10.313 * log10(Ceta) - 36.447;
    Nd_od = (3.141 * pow(10, 10)) *pow( Ceta ,-3.444) * pow((log10(D) ),aa);//Glaso法(1980)
	//Rs = Rsg_o(T, P, Md_o, Bzg)
	//Standing法(1981)
    double E1=0.0f;
	double D1=0.0f;
	double c1=0.0f;
	double B1=0.0f;
	double a1=0.0f;
	double Nd_ob=0.0f;
	double Nd_on=0.0f;
    E1 = 2.1 * pow(10,-2) * Rs;
    D1 = 6.18 * pow(10,-3) * Rs;
    c1 = 4.84 * pow(10,-4) * Rs;
    B1 = 0.68 * pow(10,-c1) + 0.25 * pow(10,-D1) + 0.062 * pow(10,-E1);
    a1 = (5.615 * Rs) * (1.2353 * pow(10,-6) * Rs - 7.4 * pow(10,-4));
    Nd_ob = pow(10,a1) * pow(Nd_od,B1);
    Nd_on = Nd_ob + 0.14503 * (P - Pb) * (0.024 * pow(Nd_ob,1.6) + 0.038 * pow(Nd_ob,0.56));
    if (P <= Pb)
	{
         return Nd_ob;
	}
    else
	{
         return Nd_on;
	}

}

double Nd_gas(double Bz_g,double P_yali,double Temperature,double Zg)//天然气粘度比重 压力 温度 偏差因子
{
	double Yinz_yas=0.01;
	double Mg=0.01;
	double T_t=0.01;
	double K_k=0.01;
	double X_x=0.01;
	double  Y_y =0.01;
	double  Midu_g =0.01;
	Mg = 29 * Bz_g;
    T_t = Temperature + 273;
    K_k = ((9.4 + 0.02 * Mg) * pow((1.8 * T_t) , 1.5)) / (209 + 19 * Mg + 1.8 * T_t);
    X_x = 3.5 + 986 / (1.8 * T_t) + 0.01 * Mg;
    Y_y = 2.4 - 0.2 * X_x;
    Midu_g = 3.4844 * Bz_g * P_yali / (Zg * T_t);
	return pow( 10 ,-4) * K_k * exp(X_x * pow(Midu_g, Y_y));  //mpa·s
}

double Md_gas(double Bz_g,double P_yali,double Temperature,double Zg)//天然气密度 比重 压力 温度 偏差因子
{
	double Yinz_yas=0.01;
	double T_t=0.01;
    T_t = Temperature; //+ 273;
    return  3484.4 * Bz_g * P_yali / (Zg * T_t);
}


double Md_water(double T,double Khd)//T——C；   Md_w——kg/m^3  水密度 温度 矿化度
{
    if (Khd = 0) 
	{
		return pow(10,3) * (0.996732 - 0.461464 * pow(10,-4) * T - 0.306254 * pow(10,-5) * pow(T,2));
	}
    else
	{
        return pow(10,3) * (1.083886 - 5.10546 * pow(10,-4) * T - 3.06254 * pow(10,-6) * pow(T,2));
	}
}

double fm(double Czd_xd, double Nrem, double Qg)//以下计算摩阻系数fm
{
   double Ccc1=0.01;
   double Ccc2=0.01;
   double f0=0.0f;
   double f1=0.0f;
   double f2=0.0f;
   double f3=0.0f;
   int k =0;

   Ccc1 = 80 / Czd_xd;
   Ccc2 = 4160 / pow((2 * Czd_xd) , 0.85);
   if (Nrem <= 100000)
   {
	   return 64 / Nrem;
   }
   else if(Nrem > 2000 && Nrem <= Ccc1)
   {
       if (Nrem <= 100000)
	   {
	       return 0.3164 / pow(Nrem,0.25);
	   }
	   else
	   {
	       f0 = 0.1;
	   }
	   do
	   {
	       f1 =pow( f0 ,-0.5);
           f2 = 2 * log10(Nrem * pow(f0,0.5))  - 0.8;
           f3 = fabs(f1 - f2) / f1;
           if (f3 <= 0.05)
		   {
			   return f0;
		   }
           f0 = f0 - 0.01;
	   }while(f3 > 0.05);
   }
   else if(Nrem > Ccc1 & Nrem <= Ccc2 )
   {
        k = 0;
        if (Qg < 10000) 
		{
            f0 = 0.12;
		}
        else
		{
			f0 = 0.08;
		}
		do
		{
			f1 =pow( f0 ,-0.5);
            f2 = 2 * log10(2.51 / (Nrem * pow(f0, 0.5)) + 1 / (3.7 * Czd_xd)) ;
            f3 = fabs(f1 - f2) / f1;
			if(f3 <= 0.01)
			{
				return f0;
			}
            f0 = f0 - 0.001;
            k = k + 1;

		}while( f3 > 0.01);
   }
   else
   {
	   return pow((2 * log10(1 / (2 * Czd_xd) + 1.74)) ,-2);

   }

 }

double Nd_water(double T)
{
	double Ceta=0.0f;
	Ceta = 1.8 * T + 32;// 'T——C
	return exp(1.003 - 1.479 * pow(10 ,-2) * Ceta + 1.982 * pow(10 ,-5) * pow(Ceta,2));

}

int __stdcall CalfinitDT(GJJS_Parameter* GJJSParameter)
{
	//动态数组分配
	double* pfPwd=(double *)malloc(sizeof(double)*GJJSParameter->InputDTNum);
	if(pfPwd==NULL)
	{
		return 0;
	}

	double* pfTD=(double *)malloc(sizeof(double)*GJJSParameter->InputDTNum);
	if(pfTD==NULL)
	{
		return 0;
	}


	//无因次化
//	double fa=GJJSParameter->InputDiaPP[0];
//	double fb=GJJSParameter->InputDiaTP[0];
	double InputDiaNum=GJJSParameter->InputDTNum;
//	double ftp=GJJSParameter->ftp;
	double fS=GJJSParameter->fResS;
//	double mintime=0.0f;
//	double maxtime=0.0f;
//	int    TheoryNum=0;
	double fCD=0.0f;
	int i=0;
    for(i=0;i<=InputDiaNum-2;i++)
	{
		if(GJJSParameter->innerboundary==0)
		{
        	pfTD[i]=GJJSParameter->InputDiaTP[i+1]-GJJSParameter->InputDiaTP[i];
        	pfTD[i]=0.0036*(GJJSParameter->fResK)*pfTD[i]/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
        		(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
    	}
    	if(GJJSParameter->innerboundary==4 ||GJJSParameter->innerboundary==5)//表皮叠加公式
		{
         	pfTD[i]=GJJSParameter->InputDiaTP[i+1]-GJJSParameter->InputDiaTP[i];
        	pfTD[i]=0.0036*(GJJSParameter->fResK)*pfTD[i]/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
    		(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
    	}
    	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
		{
	    	if(fS<0.0)
			{
		    	GJJSParameter->fResS=0.0f;
	    		fS=GJJSParameter->fResS;
			}
	    	pfTD[i]=GJJSParameter->InputDiaTP[i+1]-GJJSParameter->InputDiaTP[i];
        	pfTD[i]=0.0036*(GJJSParameter->fResK)*pfTD[i]/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
        		(GJJSParameter->BasicProCt)*(GJJSParameter->fResXf)*(GJJSParameter->fResXf));
        }

	}


	if(GJJSParameter->innerboundary==0)
	{
    	fCD=(GJJSParameter->fResC)/
	    	(2.0*3.1415926*(GJJSParameter->BasicProFai)*(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProH)*
    		(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
//	    ftp=(1.0/fCD)*0.0036*(GJJSParameter->fResK)*ftp/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
//	    	(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
	}
	if(GJJSParameter->innerboundary==4 ||GJJSParameter->innerboundary==5)//表皮叠加公式
	{
    	fCD=(GJJSParameter->fResC)/
	    	(2.0*3.1415926*(GJJSParameter->BasicProFai)*(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProH)*
    		(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
 //	    ftp=0.0036*(GJJSParameter->fResK)*ftp/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
//	    	(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProRw)*(GJJSParameter->BasicProRw));
	}
	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
	{
		if(fS<0.0)
		{
			GJJSParameter->fResS=0.0f;
			fS=GJJSParameter->fResS;
		}
		fCD=(GJJSParameter->fResC)/
	    	(2.0*3.1415926*(GJJSParameter->BasicProFai)*(GJJSParameter->BasicProCt)*(GJJSParameter->BasicProH)*
    		(GJJSParameter->fResXf)*(GJJSParameter->fResXf));
//        ftp=0.0036*(GJJSParameter->fResK)*ftp/((GJJSParameter->BasicProFai)*(GJJSParameter->BasicProMiuo)*
//	    	(GJJSParameter->BasicProCt)*(GJJSParameter->fResXf)*(GJJSParameter->fResXf));
	}
	////////////////////////////////////////////////////////////
	//计算无因次理论压力（压恢）
    int j;
	int pnumber=0;
    double ff=0.00f;
	int nmin;
	int nmax;
	double timestep;
	double fd1=0.0f,fd2=0.0f;
	double fd3=0.0f,fd4=0.0f;
	double fretype=0.0f;
	double flinelong=0.0f;
	if(GJJSParameter->innerboundary==0 )//无因此化
	{
		fretype=exp(fS)/GJJSParameter->BasicProRw;
	}
	if(GJJSParameter->innerboundary==4 || GJJSParameter->innerboundary==5)//无因此化
	{
		fretype=1.0/GJJSParameter->BasicProRw;
	}
	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)
	{
		fretype=1.0/GJJSParameter->fResXf;
	}

	if(GJJSParameter->mboundary==1 || GJJSParameter->mboundary==8)//一条边界
	{
	    fd1=GJJSParameter->fResBD1*fretype;
	}
	if(GJJSParameter->mboundary==2 || GJJSParameter->mboundary==3
			|| GJJSParameter->mboundary==4 || GJJSParameter->mboundary==5)//2条边界
	{
		fd1=GJJSParameter->fResBD1*fretype;
		fd2=GJJSParameter->fResBD2*fretype;
	}
	if(GJJSParameter->mboundary==6)//3条边界
	{
		fd1=GJJSParameter->fResBD1*fretype;
		fd2=GJJSParameter->fResBD2*fretype;
		fd3=GJJSParameter->fResBD3*fretype;
	}
	if(GJJSParameter->mboundary==7)//4条边界
	{
		fd1=GJJSParameter->fResBD1*fretype;
		fd2=GJJSParameter->fResBD2*fretype;
		fd3=GJJSParameter->fResBD3*fretype;
		fd4=GJJSParameter->fResBD4*fretype;
	}

	if(GJJSParameter->innerboundary==0 )//有效井井
	{
    	for(i=0;i<=InputDiaNum-2;i++)
		{
	 	          int num;
		          num=0;
		          double numpwd;
		          double s[8];
		          double nump;
		          nump=0;
		          numpwd=0.0f;
	              for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/(pfTD[i]*(1.0/fCD));
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
	    	      pfPwd[i]=log(2)*numpwd/(pfTD[i]*(1.0/fCD));
 
		    	 // pfPwd[i]=pfPwd[i]+log(2)*numpwd/((pfTD[i]*(1.0/fCD)));	

		  }
		
	}

	if(GJJSParameter->innerboundary==4 || GJJSParameter->innerboundary==5)//井
	{
    	for(i=0;i<=InputDiaNum-2;i++)
		{
	 	          int num;
		          num=0;
		          double numpwd;
		          double s[8];
		          double nump;
		          nump=0;
		          numpwd=0.0f;
	              for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/(pfTD[i]);
			           numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
				  }
	    	      pfPwd[i]=log(2)*numpwd/(pfTD[i]);	    	 		    	  
		}
	}
	if(GJJSParameter->innerboundary==1 || GJJSParameter->innerboundary==2)//裂缝
	{
		for(i=0;i<=InputDiaNum-2;i++)
		{
	     	      int num;
	    	      num=0;
		          double numpwd;
		          double s[8];
	    	      double nump;
	    	      nump=0;
		          numpwd=0.0f;
	              for(num=0;num<=7;num++)
				  {
		               s[num]=log(2)*(num+1)/(pfTD[i]);
	    		       numpwd+=V(num+1)*laplas(s[num],fCD,fS,GJJSParameter,1);
					   if(GJJSParameter->mboundary>=1)
					   {
				    	  numpwd+=BoundaryLap(s[num],fCD,fS,GJJSParameter,fd1,fd2,fd3,fd4,V(num+1));
					   }
					   
				  }
		          pfPwd[i]=log(2)*numpwd/(pfTD[i]);
		}
	}

	GJJSParameter->OutputDTP[0]=GJJSParameter->InputDTP[0];
	GJJSParameter->OutputDTQ[0]=GJJSParameter->InputDTQ[0];

	double fzi=0.0f;//ping jun pian cha ying zhi
	double fmiui=0.0f;//ping jun nian du
	double deltap=0.0f;
	double deltaq=0.0f;
    fzi=Z_yasyz(GJJSParameter->BasicProPi, GJJSParameter->BasicProT, GJJSParameter->BasicBZG);
	fmiui=Nd_gas( GJJSParameter->BasicBZG,GJJSParameter->BasicProPi,GJJSParameter->BasicProT,fzi);
	for(i=0;i<=InputDiaNum-2;i++)
	{
		deltap=	GJJSParameter->OutputDTPI[i+1]*GJJSParameter->OutputDTPI[i+1]-
			(pfPwd[i]*0.01273*fzi*(GJJSParameter->BasicProT+273)*GJJSParameter->InputDTQ[i+1]*fmiui)/
			(0.001*GJJSParameter->fResK*GJJSParameter->BasicProH);	
		deltaq=0.001*GJJSParameter->fResK*GJJSParameter->BasicProH*(GJJSParameter->OutputDTPI[i+1]*GJJSParameter->OutputDTPI[i+1]-GJJSParameter->InputDTP[i+1]*GJJSParameter->InputDTP[i+1])/
			(pfPwd[i]*0.01273*fzi*(GJJSParameter->BasicProT+273)*fmiui);

		if(deltap>0)
		{
			GJJSParameter->OutputDTP[i+1]=sqrt(deltap);
			GJJSParameter->OutputDTQ[i+1]=deltaq;
		}
		else
		{
			GJJSParameter->OutputDTP[i+1]=0;
			GJJSParameter->OutputDTQ[i+1]=0;

		}
	}

	free(pfTD);
	free(pfPwd);

	return 0;
}

int __stdcall PZCover(GJJS_Parameter* GJJSParameter)
{
/*	double* pfPiz=(double *)malloc(sizeof(double)*GJJSParameter->InputDTNum);
	if(pfPiz==NULL)
	{
		return 0;
	}*/
	double fzi=0.0f;
	fzi=Z_yasyz(GJJSParameter->BasicProPi, GJJSParameter->BasicProT, GJJSParameter->BasicBZG);
	double fpichazhi[1000];
	double fpizchazhi[1000];
	double fzg=0.0f;
	double fPiZmax=0.0f;
	fPiZmax=GJJSParameter->BasicProPi*2;
	double fdelta=0.0f;
	fdelta=fPiZmax/1000;
	int i=0;
	for(i=0;i<1000;i++)
	{
		fpichazhi[i]=0.2+fdelta*i;
		fzg= Z_yasyz(fpichazhi[i], GJJSParameter->BasicProT, GJJSParameter->BasicBZG);
		fpizchazhi[i]=fpichazhi[i]/fzg;
	}
	int j=0;
	GJJSParameter->OutputDTPI[0]=GJJSParameter->BasicProPi;
	double GP=0.0f;//leichanqiliang
	double GasG=0.0f;//chuliang
	GasG = GJJSParameter->BasicProFai*GJJSParameter->BasicProH*PI*GJJSParameter->fResBD1*GJJSParameter->fResBD1/GJJSParameter->BasicProBo;
	GJJSParameter->GasG = GasG;
	double pfPiz=0.0f;

	for(i=1;i<GJJSParameter->InputDTNum;i++)
	{
		GJJSParameter->OutputDTPI[i]=0.1;
		GP=GP+GJJSParameter->InputDTQ[i-1]*10000;
		pfPiz=(GJJSParameter->BasicProPi/fzi)*(1-GP/GasG);
		if (pfPiz <= 0)
		{
			GJJSParameter->OutputDTPI[i] = 0.1;
		} 
		else 
		{
			for(j=1;j<1000;j++)
			{
				if(fpizchazhi[j]>pfPiz)
				{
					GJJSParameter->OutputDTPI[i]=fpichazhi[j-1]+(fpichazhi[j]-fpichazhi[j-1])*(pfPiz-fpizchazhi[j-1])/(fpizchazhi[j]-fpizchazhi[j-1]);
					break;
				}
			}
		}
	}
//	free(pfPiz);
	return 0;
}
int __stdcall ChanLianYuCe(GJJS_Parameter* GJJSParameter)
{
	return 0;
}

double ffsy(double fu,GJJS_Parameter* GJJSParameter)
{
    double fw;
	double flamda;
	double fbeita;
    double felta;
	double fd;
	if(GJJSParameter->fResWomiga>0.1)
	{
		GJJSParameter->fResWomiga=0.1;
	}
	if(GJJSParameter->fResksxs<10000)
	{
		GJJSParameter->fResksxs=10000;

	}

	fw=GJJSParameter->fResWomiga;
	flamda=GJJSParameter->fResksxs;
	felta=0.0f;
	fbeita=1000;
	fd=0.0f;
	felta=1.0/(tanh(sqrt(flamda*fu)));
	fd=fw*fu+((1-fw)*fbeita*(sqrt(flamda*fu)*felta-1))/flamda;
	return fd;
}

