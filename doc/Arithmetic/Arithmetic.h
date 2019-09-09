typedef struct GJJS_Parameter{
	int InputDiaNum;
	double* InputDiaTP;
	double* InputDiaPP;	//导入的原始流量数据

	double* FixDcTP;
	double* FixDcPP;	//导入的原始流量数据

	// 动态拟合
	int InputDTNum;
	double* InputDTT;  //导入的动态数据-时间
	double* InputDTP;  //导入的动态数据-压力
	double* InputDTQ;  //导入的动态数据-流量

	double* OutputDTT;//输出动态数据-时间
	double* OutputDTP;//输出动态数据-时间
	double* OutputDTQ;//输出动态数据-时间
	double* OutputDTPI;//输出动态数据-原始压力

	//油藏参数
	double BasicProH;
	double BasicProFai;
	double BasicProRw;
	double BasicProMiuo;
	double BasicProBo;
	double BasicProCt;
	double BasicProQ;
	double BasicProPi;	// 原始压力
	double BasicProT;	// 原始温度
	double BasicBZG;	// 气比重

	//均质
	double fResK;//结果
	double fResS;
    
    int IsArray;         // 是否采用数组形式，0表示采用单个的表皮系数，1表示数组
    int RessCount;       // 采用数组形式的时候的数量
    double* TimeArray;   // 时间数组
    double* RessArray;   // 表皮系数数组

	double fResC;
	double fResksxs;     // 页岩扩散系数；//
	double fRescrxs;     // 页岩储容系数; //
    double fResjxxs;     // 页岩解析系数; //

	//双重介质
	double fResWomiga;
	double fResNamuda;

	//三重
	double fResWomiga1;
	double fResNamuda1;
	double fResWomiga2;
	double fResNamuda2;

	//部分射孔
	double fResHpt;
	double fResHpd;

	//两区复合
	double fResM12;
	double fResEta12;
	double fResRdi;
	double fResRdo;

	//压裂模型
	double fResFcd;
	double fResXf;
	
	//边界距离
	double fResBD1;
	double fResBD2;
	double fResBD3;
	double fResBD4;
	double fQidongyali;
	
	double ftp;
	
	int modeltype;
	int innerboundary;
	int mboundary;

	double GasG;

	// 预测模块
	int fyccount;           // 预测数组个数//
	double fycyccl[30];     // 预测产量//
	double fycjxyl[30];     // 预测极限压力//
	double fycycsj[30];     // 预测时间//

    int ycsjcount;			// 预测数据数量
	double ycsjycsj[50000]; // 预测时间//
	double fycdcyl[50000];  // 预测地层压力//
	double fycjdly[50000];  // 预测井底压力//
	double fycdjcl[50000];  // 预测单井产量//
	double fycljcl[50000];  // 预测累积产量//

};

double ffsy(double fu,GJJS_Parameter* GJJSParameter);

extern "C" int __declspec(dllexport) __stdcall Derivate(double *pfX, double *pfY, int iNum, int bLn, int iN, double *pfdY);
extern "C" double __declspec(dllexport) __stdcall I0(double x);
extern "C" double __declspec(dllexport) __stdcall K0(double x);
extern "C" double __declspec(dllexport) __stdcall I1(double x);
extern "C" double __declspec(dllexport) __stdcall K1(double x);
extern "C" double __declspec(dllexport) __stdcall IntegelK0x(double fx);
extern "C" double __declspec(dllexport) __stdcall multp(int nx);
extern "C" double __declspec(dllexport) __stdcall IntegelK0(double fx);
extern "C" double __declspec(dllexport) __stdcall Ki2(double fu);
extern "C" double __declspec(dllexport) __stdcall IntegelI0(double fx);
extern "C" double __declspec(dllexport) __stdcall Gam1(double x);
extern "C" double __declspec(dllexport) __stdcall Gamm(double x);
extern "C" double __declspec(dllexport) __stdcall Gam2To3(double x);
extern "C" double __declspec(dllexport) __stdcall Iv(double v, double z);
extern "C" double __declspec(dllexport) __stdcall Kv(double v, double z);
extern "C" double __declspec(dllexport) __stdcall V(int nx);

extern "C" double __declspec(dllexport) __stdcall Spole1(int jl1, int jr1, int jl2, int jr2, double *x, double *y);
extern "C" double __declspec(dllexport) __stdcall Spole(double x, double x1, double x2, double y1, double y2);
extern "C" int __declspec(dllexport) __stdcall DataDerivateSmooth(	double *td,      //输入时间数据序列点
																	double *pd,      //输入压力数据序列点
																	int    npoint,   //数据序列点数
																	int    Flag,     //0为对数求导,1为自然求导
																	double l,
																	double *pdd);    //输出压力导数数据序列点;        //磨光系数

extern "C" int __declspec(dllexport) __stdcall funJZZZ(double** fA, double** fB, int NumX, int NumY);
extern "C" int __declspec(dllexport) __stdcall funJZCJ(double** fA, double** fB, double** fC, int NumAX, int NumAY, int NumBX, int NumBY);
extern "C" int __declspec(dllexport) __stdcall funJZN(double** fA, double** fB, int NumX);

extern "C" double __declspec(dllexport) __stdcall ffs(double fu, GJJS_Parameter* GJJSParameter, double fCDe2S, double fS);
extern "C" double __declspec(dllexport) __stdcall ffs3(double fu, GJJS_Parameter* GJJSParameter, double fCDe2S, double fS);

extern "C" double __declspec(dllexport) __stdcall laplas(	double fu, 
															double fCD, 
															double fS,
															GJJS_Parameter* GJJSParameter,
															double frd);
extern "C" double __declspec(dllexport) __stdcall BoundaryLap(	double fu, 
																double fCD, 
																double fS, 
																GJJS_Parameter* GJJSParameter,
																double fd1,
																double fd2, 
																double fd3,
																double fd4,
																double fvnum);
extern "C" int __declspec(dllexport) __stdcall Calfinit(double* pfPwd, 
														double* pfTD, 
														int& number, 
														GJJS_Parameter* GJJSParameter);
extern "C" int __declspec(dllexport) __stdcall CalAutofix(GJJS_Parameter* GJJSParameter);

extern "C" double __declspec(dllexport) __stdcall Z_yasyz(double P_pj ,double T_pj,double Bz_g); //偏差系数计算
extern "C" double __declspec(dllexport) __stdcall B_water(double T,double P);//水体积系数
extern "C" double __declspec(dllexport) __stdcall Rsg_w(double T,double P,double s);//溶解气水比，s,矿化度%
extern "C" double __declspec(dllexport) __stdcall Rsg_o(double T,double P,double Md_o,double Bz_g);//Glaso法(1980) 溶解汽油比温度 压力 密度 气比重 有问题
extern "C" double __declspec(dllexport) __stdcall B_oil(double T,double p,double Bz_g,double Md_o,double Rs );//Glaso法(1980)油体积系数
extern "C" double __declspec(dllexport) __stdcall Nd_oil(double T,double P,double Pb,double Md_o,double Bzg ,double Rs);//油粘度
extern "C" double __declspec(dllexport) __stdcall Nd_gas(double Bz_g,double P_yali,double Temperature,double Zg);//天然气粘度比重 压力 温度 偏差因子
extern "C" double __declspec(dllexport) __stdcall Md_gas(double Bz_g,double P_yali,double Temperature,double Zg);//天然气密度 比重 压力 温度 偏差因子
extern "C" double __declspec(dllexport) __stdcall Md_water(double T,double Khd);//T——C；   Md_w——kg/m^3  水密度 温度 矿化度
extern "C" double __declspec(dllexport) __stdcall fm(double Czd_xd, double Nrem, double Qg);//以下计算摩阻系数fm
extern "C" double __declspec(dllexport) __stdcall Nd_water(double T);

// 井筒压降计算程序
extern "C" double __declspec(dllexport) __stdcall H_Hagedorn_Brown(double Panbie, double Qg_sc, double Qw_sc, double Qg_ing,
																   double D_youg_nj, double Tjingkou_sc, double Tjingdi_sc,
																   double Hyouc_sc, double bz_scg, double Bz_ing, double Ppr_sc,
																   double Ttr_sc, double s, double Md_w_sc, double Czd_xd, 
																   double P_qid, double Bmyz_w, double Md_o_sc, double Qo_sc);
extern "C" int __declspec(dllexport) __stdcall CalfinitDT(GJJS_Parameter* GJJSParameter);
extern "C" int __declspec(dllexport) __stdcall PZCover(GJJS_Parameter* GJJSParameter);
extern "C" int __declspec(dllexport) __stdcall ChanLianYuCe(GJJS_Parameter* GJJSParameter);

