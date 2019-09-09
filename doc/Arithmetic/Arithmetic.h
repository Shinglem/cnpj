typedef struct GJJS_Parameter{
	int InputDiaNum;
	double* InputDiaTP;
	double* InputDiaPP;	//�����ԭʼ��������

	double* FixDcTP;
	double* FixDcPP;	//�����ԭʼ��������

	// ��̬���
	int InputDTNum;
	double* InputDTT;  //����Ķ�̬����-ʱ��
	double* InputDTP;  //����Ķ�̬����-ѹ��
	double* InputDTQ;  //����Ķ�̬����-����

	double* OutputDTT;//�����̬����-ʱ��
	double* OutputDTP;//�����̬����-ʱ��
	double* OutputDTQ;//�����̬����-ʱ��
	double* OutputDTPI;//�����̬����-ԭʼѹ��

	//�Ͳز���
	double BasicProH;
	double BasicProFai;
	double BasicProRw;
	double BasicProMiuo;
	double BasicProBo;
	double BasicProCt;
	double BasicProQ;
	double BasicProPi;	// ԭʼѹ��
	double BasicProT;	// ԭʼ�¶�
	double BasicBZG;	// ������

	//����
	double fResK;//���
	double fResS;
    
    int IsArray;         // �Ƿ����������ʽ��0��ʾ���õ����ı�Ƥϵ����1��ʾ����
    int RessCount;       // ����������ʽ��ʱ�������
    double* TimeArray;   // ʱ������
    double* RessArray;   // ��Ƥϵ������

	double fResC;
	double fResksxs;     // ҳ����ɢϵ����//
	double fRescrxs;     // ҳ�Ҵ���ϵ��; //
    double fResjxxs;     // ҳ�ҽ���ϵ��; //

	//˫�ؽ���
	double fResWomiga;
	double fResNamuda;

	//����
	double fResWomiga1;
	double fResNamuda1;
	double fResWomiga2;
	double fResNamuda2;

	//�������
	double fResHpt;
	double fResHpd;

	//��������
	double fResM12;
	double fResEta12;
	double fResRdi;
	double fResRdo;

	//ѹ��ģ��
	double fResFcd;
	double fResXf;
	
	//�߽����
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

	// Ԥ��ģ��
	int fyccount;           // Ԥ���������//
	double fycyccl[30];     // Ԥ�����//
	double fycjxyl[30];     // Ԥ�⼫��ѹ��//
	double fycycsj[30];     // Ԥ��ʱ��//

    int ycsjcount;			// Ԥ����������
	double ycsjycsj[50000]; // Ԥ��ʱ��//
	double fycdcyl[50000];  // Ԥ��ز�ѹ��//
	double fycjdly[50000];  // Ԥ�⾮��ѹ��//
	double fycdjcl[50000];  // Ԥ�ⵥ������//
	double fycljcl[50000];  // Ԥ���ۻ�����//

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
extern "C" int __declspec(dllexport) __stdcall DataDerivateSmooth(	double *td,      //����ʱ���������е�
																	double *pd,      //����ѹ���������е�
																	int    npoint,   //�������е���
																	int    Flag,     //0Ϊ������,1Ϊ��Ȼ��
																	double l,
																	double *pdd);    //���ѹ�������������е�;        //ĥ��ϵ��

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

extern "C" double __declspec(dllexport) __stdcall Z_yasyz(double P_pj ,double T_pj,double Bz_g); //ƫ��ϵ������
extern "C" double __declspec(dllexport) __stdcall B_water(double T,double P);//ˮ���ϵ��
extern "C" double __declspec(dllexport) __stdcall Rsg_w(double T,double P,double s);//�ܽ���ˮ�ȣ�s,�󻯶�%
extern "C" double __declspec(dllexport) __stdcall Rsg_o(double T,double P,double Md_o,double Bz_g);//Glaso��(1980) �ܽ����ͱ��¶� ѹ�� �ܶ� ������ ������
extern "C" double __declspec(dllexport) __stdcall B_oil(double T,double p,double Bz_g,double Md_o,double Rs );//Glaso��(1980)�����ϵ��
extern "C" double __declspec(dllexport) __stdcall Nd_oil(double T,double P,double Pb,double Md_o,double Bzg ,double Rs);//��ճ��
extern "C" double __declspec(dllexport) __stdcall Nd_gas(double Bz_g,double P_yali,double Temperature,double Zg);//��Ȼ��ճ�ȱ��� ѹ�� �¶� ƫ������
extern "C" double __declspec(dllexport) __stdcall Md_gas(double Bz_g,double P_yali,double Temperature,double Zg);//��Ȼ���ܶ� ���� ѹ�� �¶� ƫ������
extern "C" double __declspec(dllexport) __stdcall Md_water(double T,double Khd);//T����C��   Md_w����kg/m^3  ˮ�ܶ� �¶� �󻯶�
extern "C" double __declspec(dllexport) __stdcall fm(double Czd_xd, double Nrem, double Qg);//���¼���Ħ��ϵ��fm
extern "C" double __declspec(dllexport) __stdcall Nd_water(double T);

// ��Ͳѹ���������
extern "C" double __declspec(dllexport) __stdcall H_Hagedorn_Brown(double Panbie, double Qg_sc, double Qw_sc, double Qg_ing,
																   double D_youg_nj, double Tjingkou_sc, double Tjingdi_sc,
																   double Hyouc_sc, double bz_scg, double Bz_ing, double Ppr_sc,
																   double Ttr_sc, double s, double Md_w_sc, double Czd_xd, 
																   double P_qid, double Bmyz_w, double Md_o_sc, double Qo_sc);
extern "C" int __declspec(dllexport) __stdcall CalfinitDT(GJJS_Parameter* GJJSParameter);
extern "C" int __declspec(dllexport) __stdcall PZCover(GJJS_Parameter* GJJSParameter);
extern "C" int __declspec(dllexport) __stdcall ChanLianYuCe(GJJS_Parameter* GJJSParameter);

