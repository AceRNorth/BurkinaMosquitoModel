#include "headDsxModel.h"

	/*--------------------------------------global variables----------------------------------*/
	/*----------------------------------------------------------------------------------------*/
	/*-------------random number seed-----------------*/
	int32 seed = (int32)time(0);
	CRandomMersenne rg(seed);
	/*------------------------------------------------*/
	/*----------------global structs------------------*/
	Pars pa;//parameters 
	initials in; //initial conditions
	totals to; // totals
	Times ti; // timekeeping parameters
	vector<Patch> Pat; // information on each population
	/*------------------------------------------------*/
	/*------------------Output files------------------*/
	ostringstream os1,os2;
	ofstream globalinfo,localinfo;
	/*------------------------------------------------*/
	/*------------Array to store rain input data--------------*/
	double rain[nx][ny][38][52];//nx is x grid cell,ny is y grid cell, 38 is number of years in data, 52 is weeks in year
	/*--------------------------------------------------------*/
	/*----------------------------------------------------------------------------------------*/
	int main(void)
		{
	/*---------input parameters---------------------*/
		cin>>pa.set; //enumeration for output files
	//	cin>>pa.offset;
		cin>>in.inputfile;// name of file with settlement data
		cin>>in.recPatFreq;//if outputting local data, how many sites to collect for (1 is all sites, 10 is 1 in 10 etc)
		cin>>ti.interval;//time interval for outputting global data
		cin>>ti.maxT;// maximum time to simulate until (days)
		cin>>ti.rec;// how often to collect local data (every rec days)
		cin>>ti.NumRuns;// how many simulation runs
		cin>>in.driver_time;// time to start releasing drive alleles; if negative, the release dates are random in year following in.driver_time, if positive, the releases are on the specific date
		cin>>in.NumDriver;// number of drive homozygous male mosquitoes per release
		cin>>in.NumDriverPat;// number of release sites per year
		cin>>pa.muJ;//juvenile density independent mortality per day
		cin>>pa.muA;//adult mortality per day
		cin>>pa.d;//adult dispersal rate
		cin>>pa.Fgamma;//rate of r2 allele formation from W/D meiosis in females
		cin>>pa.Mgamma;//rate of r2 allele formation from W/D meiosis in males
		cin>>pa.beta;//parameter that controls mating rate
		cin>>pa.theta;//egg laying rate of wildtype females (eggs per day)
		cin>>pa.Frho;//Maternal Cas9 deposition fitness cost
		cin>>pa.Mrho;//Paternal Cas9 deposition fitness cost
		cin>>pa.xi;// somatic Cas9 expression fitness cost
		cin>>pa.ef;// homing rate in females
		cin>>pa.em;//homing rate in males
		cin>>pa.LD;//maximum distance at which settlements are connected
		cin>>pa.dL;//Long-distance migration rate
		cin>>pa.muLD;//Long-distance migration mortality
		cin>>pa.t_disp1;//start day of NE->SW migration
		cin>>pa.t_disp2;//end day of NE->SW migration
		cin>>pa.t_disp3;//start day of SW->NE migration
		cin>>pa.t_disp4;//end day of SW->NE migration
		cin>>pa.psi;//Aestivation rate
		cin>>pa.muAES;//aestivation mortality
		cin>>pa.t_hide1;//start day of going into aestivation
		cin>>pa.t_hide2;//end day of going into aestivation
		cin>>pa.t_wake1;//start day of emerging from aestivation
		cin>>pa.t_wake2;//end day of emerging from aestivation
		cin>>pa.alpha0;//baseline contribution to carrying capacity
		cin>>pa.alpha1;//Maximum contribution from rain
		cin>>pa.alpha2;//maximum contribution from water bodies
		cin>>pa.phi;// Increase in carrying capacity/alpha_1 per mm rain per week (when rainfall low)
		cin>>pa.delta;//Increase in length of standing water from non-permanent waterways per km non-permanent waterways (within 5km)per mm rain per week (when rainfall low)
		cin>>pa.kappa;//Increase in carrying capacity/α2 per km standing water (within 5km; when water bodies rare)
		cin>>pa.al0var;// variance in baseline carrying capacity
		cin>>in.rainfile;// rain data input file
		cin>>pa.bias;// male proportion in progeney of males carrying drive allele
	/*----------------------------------------------------------------------------------------*/
	/*---------------------set mean and variance of local permanent water--------------------------------*/
		if(pa.alpha0>0)
			{
			pa.mu=log(pa.alpha0*pa.alpha0/sqrt(pa.al0var+pa.alpha0*pa.alpha0));
			pa.sig=sqrt(log(1+pa.al0var/(pa.alpha0*pa.alpha0)));	
			};
	/*---------------------------------------------------------------------------------------------------*/

	/*-------------------------------input rain data file-------------------------------------------------*/
		ostringstream fff;
		string line ="";
		string cell;
		fff.str(in.rainfile);
		ifstream raindata(fff.str().c_str()); 
		int xs,ys,week,year;
		double rr;
		while( getline(raindata, line ))
		{       
			stringstream lineStream(line);
			getline(lineStream,cell,',');xs=(int)strtod(cell.c_str(),NULL)-1;
			getline(lineStream,cell,',');ys=(int)strtod(cell.c_str(),NULL)-1;
			getline(lineStream,cell,',');year=(int)strtod(cell.c_str(),NULL)-1;
			getline(lineStream,cell,',');week=(int)strtod(cell.c_str(),NULL)-1;
			getline(lineStream,cell,',');rr=strtod(cell.c_str(),NULL);
			rain[xs][ys][year][week]=rr;
		};
	/*---------------------------------------------------------------------------------------------------*/
	/*-------------------------------numbers to use when initiating populations--------------------------*/
		in.NumAdultsWV=100000; in.NumAdultsWM=500000; in.NumAdultsWF=400000;
		in.NumAdultsDV=100000; in.NumAdultsDM=500000; in.NumAdultsDF=400000;
		in.NumAdultsRV=100000; in.NumAdultsRM=500000; in.NumAdultsRF=400000;
	/*---------------------------------------------------------------------------------------------------*/
	
	/*--------------------------------Set inheritance architecture---------------------------------------*/
		SetFertility();
	/*---------------------------------------------------------------------------------------------------*/

	/*--------------------------------run model NumRuns times--------------------------------------------*/
		RunNReps(ti.NumRuns);
	/*---------------------------------------------------------------------------------------------------*/
	return 0;};
	

	void RunNReps(int N)
	      {
		for(int i=pa.offset;i<N+pa.offset;i++)
		{
		os1<<"LocalData"<<pa.set<<"run"<<i+1<<".txt";// make a file for outputing local data
		os2<<"Totals"<<pa.set<<"run"<<i+1<<".txt";// make a file for outputing global data
		localinfo.open(os1.str().c_str());
		globalinfo.open(os2.str().c_str());
			initiate();
			RunMaxT();
		os2.str("");
		os1.str("");
		globalinfo.close();
		localinfo.close();
		};
		return;};

	void RunMaxT(void)
		{
		int TT=0;
		int TTT;
		int swit=0;
		int switR=0;
		int relpatches[in.NumDriverPat];
		int reltimes[in.NumDriverPat];
		for(int ii=0;ii<in.NumDriverPat;ii++)relpatches[ii]=500000;
		int uniquepat,relpat,year2;
		while (TT<ti.maxT+1)
			{
			if(TT%365==0)
				{
				ti.yearnow=rg.IRandom(0,37);
				if(TT==365)year2=ti.yearnow;
				if(TT==1825)ti.yearnow=year2;
				if(TT==3285)ti.yearnow=year2;
				if(TT==4745)ti.yearnow=year2;
				if(TT==6205)ti.yearnow=year2;
				if(TT==7665)ti.yearnow=year2;
					for(int jj=0;jj<in.NumDriverPat;jj++)
						{
						if(in.NumDriverPat>=Pat.size()-1){relpatches[jj]=jj;}
						else{
						uniquepat=0;
						while(uniquepat==0)
							{
							uniquepat=1;
							relpat=rg.IRandom(0,Pat.size()-1);
							for(int ii=0;ii<in.NumDriverPat;ii++){if(relpat==relpatches[ii])uniquepat=0;};
							};
						relpatches[jj]=relpat;
						};
						if(in.driver_time<0)reltimes[jj]=rg.IRandom(1,364); else reltimes[jj]=in.driver_time;
						};
				};
			if(TT%ti.interval==0)
				{
				globalinfo<<TT<<"   "<<ti.yearnow<<"    "<<to.Jww<<"     "<<to.Jwd<<"     "<<to.Jdd<<"     "<<to.Jwr<<"     "<<to.Jrr<<"     "<<to.Jdr<<"   "<<to.Mww<<"     "<<to.Mwd<<"     "<<to.Mdd<<"     "<<to.Mwr<<"     "<<to.Mrr<<"     "<<to.Mdr<<"  "<<to.Vww<<"     "<<to.Vwd<<"     "<<to.Vdd<<"     "<<to.Vwr<<"     "<<to.Vrr<<"     "<<to.Vdr<<"  "<<to.Fww<<"     "<<to.Fwd<<"     "<<to.Fdd<<"     "<<to.Fwr<<"     "<<to.Frr<<"     "<<to.Fdr<<"    "<<to.distW<<"    "<<to.distD<<"   "<<to.distR<<endl;
				};
			TT++;
			if(TT>=abs(in.driver_time))
			{
				for(int jj=0;jj<in.NumDriverPat;jj++)
					{
					if(TT%365==reltimes[jj])
						{
						PutDriverPat(relpatches[jj]);
						};
					};
			};
			TTT=TT-abs(in.driver_time)+365;
			if(TTT%ti.rec==292){record(TT);};
			OneStep(TT);
			};
        return;};

					
	void initiate(void){
		to.Jww=0;to.Jwd=0;to.Jdd=0;to.Jwr=0;to.Jrr=0;to.Jdr=0;
		to.Mww=0;to.Mwd=0;to.Mdd=0;to.Mwr=0;to.Mrr=0;to.Mdr=0;
		to.Vww=0;to.Vwd=0;to.Vdd=0;to.Vwr=0;to.Vrr=0;to.Vdr=0;
		to.Fww=0;to.Fwd=0;to.Fdd=0;to.Fwr=0;to.Frr=0;to.Fdr=0;
		to.JTot=0;to.MTot=0;to.VTot=0;to.FTot=0;
		to.distW=0;to.distD=0;to.distR=0;
		Pat.clear();
		Patch pp;
		/*----------------------------------------------------------------------------------*/
		/*------------------------------input the settlement data---------------------------*/
		ostringstream ddd;
		ddd.str(in.inputfile);
		ifstream indata(ddd.str().c_str()); 
		string line ="";
		string cell,name,type;
		while( getline(indata, line ))
		{       
			stringstream lineStream(line);
			getline(lineStream,cell,','); pp.x=strtod(cell.c_str(),NULL);
			getline(lineStream,cell,','); pp.y=strtod(cell.c_str(),NULL);
			getline(lineStream,cell,',');pp.WaterTemp=strtod(cell.c_str(),NULL);
			getline(lineStream,cell,',');pp.WaterPerm=strtod(cell.c_str(),NULL);
			getline(lineStream,cell,',');pp.sqx=strtod(cell.c_str(),NULL)-1;
			getline(lineStream,cell,',');pp.sqy=strtod(cell.c_str(),NULL)-1;
			getline(lineStream,cell,',');type=cell;
			pp.connecIND.clear();
			for(int i=0;i<TL;i++)
				{
				pp.Jww[i]=0; pp.fJww[i]=0; pp.mJww[i]=0; pp.bJww[i]=0; pp.fJwd[i]=0; pp.mJwd[i]=0; pp.bJwd[i]=0; pp.Jdd[i]=0; pp.Jwr[i]=0; pp.fJwr[i]=0; pp.mJwr[i]=0; pp.bJwr[i]=0; pp.Jrr[i]=0;pp.Jdr[i]=0;pp.mbJrr[i]=0;pp.mbJdr[i]=0;
				};
			pp.Mww=0;pp.Mwd=0;pp.Mdd=0;pp.Mwr=0;pp.Mrr=0;pp.Mdr=0;
			pp.Vww=0;pp.fVww=0;pp.mVww=0;pp.bVww=0;pp.Vdd=0;pp.fVwd=0;pp.mVwd=0;pp.bVwd=0;pp.Vwr=0;pp.fVwr=0;pp.mVwr=0;pp.bVwr=0;pp.Vrr=0;pp.Vdr=0;
			pp.Fwwww=0;pp.Fwwwd=0;pp.Fwwdd=0;pp.Fwwwr=0;pp.Fwwrr=0;pp.Fwwdr=0;
			pp.fFwwww=0;pp.fFwwwd=0;pp.fFwwdd=0;pp.fFwwwr=0;pp.fFwwrr=0;pp.fFwwdr=0;
			pp.mFwwww=0;pp.mFwwwd=0;pp.mFwwdd=0;pp.mFwwwr=0;pp.mFwwrr=0;pp.mFwwdr=0;
			pp.bFwwww=0;pp.bFwwwd=0;pp.bFwwdd=0;pp.bFwwwr=0;pp.bFwwrr=0;pp.bFwwdr=0;
			pp.fFwdww=0;pp.fFwdwd=0;pp.fFwddd=0;pp.fFwdwr=0;pp.fFwdrr=0;pp.fFwddr=0;
			pp.mFwdww=0;pp.mFwdwd=0;pp.mFwddd=0;pp.mFwdwr=0;pp.mFwdrr=0;pp.mFwddr=0;
			pp.bFwdww=0;pp.bFwdwd=0;pp.bFwddd=0;pp.bFwdwr=0;pp.bFwdrr=0;pp.bFwddr=0;
			pp.Fddww=0;pp.Fddwd=0;pp.Fdddd=0;pp.Fddwr=0;pp.Fddrr=0;pp.Fdddr=0;
			pp.Fwrww=0;pp.Fwrwd=0;pp.Fwrdd=0;pp.Fwrwr=0;pp.Fwrrr=0;pp.Fwrdr=0;
			pp.fFwrww=0;pp.fFwrwd=0;pp.fFwrdd=0;pp.fFwrwr=0;pp.fFwrrr=0;pp.fFwrdr=0;
			pp.mFwrww=0;pp.mFwrwd=0;pp.mFwrdd=0;pp.mFwrwr=0;pp.mFwrrr=0;pp.mFwrdr=0;
			pp.bFwrww=0;pp.bFwrwd=0;pp.bFwrdd=0;pp.bFwrwr=0;pp.bFwrrr=0;pp.bFwrdr=0;
			pp.Frrww=0;pp.Frrwd=0;pp.Frrdd=0;pp.Frrwr=0;pp.Frrrr=0;pp.Frrdr=0;
			pp.Fdrww=0;pp.Fdrwd=0;pp.Fdrdd=0;pp.Fdrwr=0;pp.Fdrrr=0;pp.Fdrdr=0;
			pp.AesFwwww=0;pp.AesFwwwd=0;pp.AesFwwdd=0;pp.AesFwwwr=0;pp.AesFwwrr=0;pp.AesFwwdr=0;
			pp.AesfFwwww=0;pp.AesfFwwwd=0;pp.AesfFwwdd=0;pp.AesfFwwwr=0;pp.AesfFwwrr=0;pp.AesfFwwdr=0;
			pp.AesmFwwww=0;pp.AesmFwwwd=0;pp.AesmFwwdd=0;pp.AesmFwwwr=0;pp.AesmFwwrr=0;pp.AesmFwwdr=0;
			pp.AesbFwwww=0;pp.AesbFwwwd=0;pp.AesbFwwdd=0;pp.AesbFwwwr=0;pp.AesbFwwrr=0;pp.AesbFwwdr=0;
			pp.AesfFwdww=0;pp.AesfFwdwd=0;pp.AesfFwddd=0;pp.AesfFwdwr=0;pp.AesfFwdrr=0;pp.AesfFwddr=0;
			pp.AesmFwdww=0;pp.AesmFwdwd=0;pp.AesmFwddd=0;pp.AesmFwdwr=0;pp.AesmFwdrr=0;pp.AesmFwddr=0;
			pp.AesbFwdww=0;pp.AesbFwdwd=0;pp.AesbFwddd=0;pp.AesbFwdwr=0;pp.AesbFwdrr=0;pp.AesbFwddr=0;
			pp.AesFddww=0;pp.AesFddwd=0;pp.AesFdddd=0;pp.AesFddwr=0;pp.AesFddrr=0;pp.AesFdddr=0;
			pp.AesFwrww=0;pp.AesFwrwd=0;pp.AesFwrdd=0;pp.AesFwrwr=0;pp.AesFwrrr=0;pp.AesFwrdr=0;
			pp.AesfFwrww=0;pp.AesfFwrwd=0;pp.AesfFwrdd=0;pp.AesfFwrwr=0;pp.AesfFwrrr=0;pp.AesfFwrdr=0;
			pp.AesmFwrww=0;pp.AesmFwrwd=0;pp.AesmFwrdd=0;pp.AesmFwrwr=0;pp.AesmFwrrr=0;pp.AesmFwrdr=0;
			pp.AesbFwrww=0;pp.AesbFwrwd=0;pp.AesbFwrdd=0;pp.AesbFwrwr=0;pp.AesbFwrrr=0;pp.AesbFwrdr=0;
			pp.AesFrrww=0;pp.AesFrrwd=0;pp.AesFrrdd=0;pp.AesFrrwr=0;pp.AesFrrrr=0;pp.AesFrrdr=0;
			pp.AesFdrww=0;pp.AesFdrwd=0;pp.AesFdrdd=0;pp.AesFdrwr=0;pp.AesFdrrr=0;pp.AesFdrdr=0;
			pp.LDMFwwww=0;pp.LDMFwwwd=0;pp.LDMFwwdd=0;pp.LDMFwwwr=0;pp.LDMFwwrr=0;pp.LDMFwwdr=0;
			pp.LDMfFwwww=0;pp.LDMfFwwwd=0;pp.LDMfFwwdd=0;pp.LDMfFwwwr=0;pp.LDMfFwwrr=0;pp.LDMfFwwdr=0;
			pp.LDMmFwwww=0;pp.LDMmFwwwd=0;pp.LDMmFwwdd=0;pp.LDMmFwwwr=0;pp.LDMmFwwrr=0;pp.LDMmFwwdr=0;
			pp.LDMbFwwww=0;pp.LDMbFwwwd=0;pp.LDMbFwwdd=0;pp.LDMbFwwwr=0;pp.LDMbFwwrr=0;pp.LDMbFwwdr=0;
			pp.LDMfFwdww=0;pp.LDMfFwdwd=0;pp.LDMfFwddd=0;pp.LDMfFwdwr=0;pp.LDMfFwdrr=0;pp.LDMfFwddr=0;
			pp.LDMmFwdww=0;pp.LDMmFwdwd=0;pp.LDMmFwddd=0;pp.LDMmFwdwr=0;pp.LDMmFwdrr=0;pp.LDMmFwddr=0;
			pp.LDMbFwdww=0;pp.LDMbFwdwd=0;pp.LDMbFwddd=0;pp.LDMbFwdwr=0;pp.LDMbFwdrr=0;pp.LDMbFwddr=0;
			pp.LDMFddww=0;pp.LDMFddwd=0;pp.LDMFdddd=0;pp.LDMFddwr=0;pp.LDMFddrr=0;pp.LDMFdddr=0;
			pp.LDMFwrww=0;pp.LDMFwrwd=0;pp.LDMFwrdd=0;pp.LDMFwrwr=0;pp.LDMFwrrr=0;pp.LDMFwrdr=0;
			pp.LDMfFwrww=0;pp.LDMfFwrwd=0;pp.LDMfFwrdd=0;pp.LDMfFwrwr=0;pp.LDMfFwrrr=0;pp.LDMfFwrdr=0;
			pp.LDMmFwrww=0;pp.LDMmFwrwd=0;pp.LDMmFwrdd=0;pp.LDMmFwrwr=0;pp.LDMmFwrrr=0;pp.LDMmFwrdr=0;
			pp.LDMbFwrww=0;pp.LDMbFwrwd=0;pp.LDMbFwrdd=0;pp.LDMbFwrwr=0;pp.LDMbFwrrr=0;pp.LDMbFwrdr=0;
			pp.LDMFrrww=0;pp.LDMFrrwd=0;pp.LDMFrrdd=0;pp.LDMFrrwr=0;pp.LDMFrrrr=0;pp.LDMFrrdr=0;
			pp.LDMFdrww=0;pp.LDMFdrwd=0;pp.LDMFdrdd=0;pp.LDMFdrwr=0;pp.LDMFdrrr=0;pp.LDMFdrdr=0;
			pp.comp=0;pp.mate_rate=0;pp.JTot=0;pp.MTot=0;
			pp.distOrigin=0;
			if(pa.alpha0>0)pp.alpha0=exp(random_normal(pa.mu,pa.sig));else pp.alpha0=0;
			Pat.push_back(pp);
			if(type=="W") { PatPopulate(Pat.size()-1,'W'); };
			if(type=="D") {  PatPopulate(Pat.size()-1,'D');};
			if(type=="R") {  PatPopulate(Pat.size()-1,'R');};
			SetsPerCell[pp.sqx][pp.sqy].push_back(Pat.size()-1);
		};
		/*----------------------------------------------------------------------------------*/
		for(int ind=0;ind<Pat.size();ind++) { UpdateConnec(ind); };
	return;};

	void UpdateConnec(int index)
	{
		double dd,ww;
		for(int ii=0;ii<Pat.size();ii++)Pat[index].connecIND.clear();
		for(int ii=0;ii<Pat.size();ii++)Pat[index].connecW.clear();
		Pat[index].TotW=0;
		for(int ii=0;ii<Pat.size();ii++)
			{
				dd=dist(Pat[index].x,Pat[index].y,Pat[ii].x,Pat[ii].y);
				if(dd<pa.LD)
				{
					Pat[index].connecIND.push_back(ii); 
					ww=1-dd/pa.LD;
					Pat[index].connecW.push_back(ww); 
					Pat[index].TotW+=ww;
				};
			};
		Pat[index].distOrigin=dist(Pat[index].x,Pat[index].y,Pat[0].x,Pat[0].y);
	return;};

	void PutDriver(void){
		int pat;
		int test=0;
		while(test==0)
		{
		pat=rg.IRandom(0,Pat.size()-1);
		if(Pat[pat].type=='W') test=1;
		};
		for(int ii=0;ii<in.NumDriver;ii++)
		{
			if(in.dist=='u')pat=rg.IRandom(0,Pat.size()-1);
			Pat[pat].Mdd++;
			Pat[pat].MTot++;
			to.Mdd++;
			to.MTot++;
		};
		UpdateComp(0);
		UpdateMate();
		return;};


	void PutDriverPat(int pat){
		for(int ii=0;ii<in.NumDriver;ii++)
		{
			Pat[pat].Mdd++;
			Pat[pat].MTot++;
			to.Mdd++;
			to.MTot++;
		};
		UpdateMate();
		return;};


	void PatPopulate(int pat,char type){
		if(type=='W')
		{
			for(int a=0;a<TL;a++){Pat[pat].Jww[a]+=in.NumJW[a];to.Jww+=in.NumJW[a];to.JTot+=in.NumJW[a];Pat[pat].JTot+=in.NumJW[a];};
			Pat[pat].Mww=int(in.NumAdultsWM);to.Mww+=int(in.NumAdultsWM);to.MTot+=int(in.NumAdultsWM);Pat[pat].MTot+=int(in.NumAdultsWM);
			Pat[pat].Vww=int(in.NumAdultsWV);to.Vww+=int(in.NumAdultsWV);to.VTot+=int(in.NumAdultsWV);
			Pat[pat].Fwwww=int(in.NumAdultsWF);to.Fww+=int(in.NumAdultsWF);to.FTot+=int(in.NumAdultsWF);
			Pat[pat].type='W';
			if(Pat[pat].distOrigin>to.distW)to.distW=Pat[pat].distOrigin;
		};
		if (type=='D')
		{
			for(int a=0;a<TL;a++){Pat[pat].Jdd[a]+=in.NumJD[a];to.Jdd+=in.NumJD[a];to.JTot+=in.NumJD[a];Pat[pat].JTot+=in.NumJD[a];};
			Pat[pat].Mdd=int(in.NumAdultsDM);to.Mdd+=int(in.NumAdultsDM);to.MTot+=int(in.NumAdultsDM);Pat[pat].MTot+=int(in.NumAdultsDM);
			Pat[pat].Vdd=int(in.NumAdultsDV);to.Vdd+=int(in.NumAdultsDV);to.VTot+=int(in.NumAdultsDV);
			Pat[pat].Fdddd=int(in.NumAdultsDF);to.Fdd+=int(in.NumAdultsDF);to.FTot+=int(in.NumAdultsDF);
			Pat[pat].type='D';
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};
		if(type=='R')
		{
			for(int a=0;a<TL;a++){Pat[pat].Jrr[a]+=in.NumJR[a];to.Jrr+=in.NumJR[a];to.JTot+=in.NumJR[a];Pat[pat].JTot+=in.NumJR[a];};
			Pat[pat].Mrr=int(in.NumAdultsRM);to.Mrr+=int(in.NumAdultsRM);to.MTot+=int(in.NumAdultsRM);Pat[pat].MTot+=int(in.NumAdultsRM);
			Pat[pat].Vrr=int(in.NumAdultsRV);to.Vrr+=int(in.NumAdultsRV);to.VTot+=int(in.NumAdultsRV);
			Pat[pat].Frrrr=int(in.NumAdultsRF);to.Frr+=int(in.NumAdultsRF);to.FTot+=int(in.NumAdultsRF);
			Pat[pat].type='R';
			if(Pat[pat].distOrigin>to.distR)to.distR=Pat[pat].distOrigin;
		};
		if(type=='E')Pat[pat].type='E';
			Pat[pat].comp=std::pow((double)pa.alpha1/(pa.alpha1+Pat[pat].JTot),1.0/double(TL));
			Pat[pat].mate_rate=Pat[pat].MTot/(pa.beta+Pat[pat].MTot);
	return;};
	void OneStep(int day){
		if(day%365 > pa.t_disp1 && day%365<=pa.t_disp2 && pa.dL>0.00001)LDM('S');
		if(day%365 > pa.t_disp3 && day%365<=pa.t_disp4 && pa.dL>0.00001)LDM('N');
		JuvGetOlder();
		AdultsDie();
		VirginsMate();
		AdultsMove();
		LayEggs();
		JuvEmerge();
		if(day%365 > pa.t_hide1 && day%365<=pa.t_hide2 && pa.psi>0.00001)Hide();
		if(day%365 > pa.t_wake1 && day%365<=pa.t_wake2 &&pa.psi>0.00001)Wake(day);
		UpdateComp(day);
		UpdateMate();
	return;};
	

	void JuvEmerge(void){
		int surv,survM;
		double pcomp;
		for(int pat=0;pat<Pat.size();pat++)
		{
		pcomp=(1-pa.muJ)*Pat[pat].comp;
		surv=random_binomial(Pat[pat].Jww[TL-1],pcomp); 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mww+=survM; to.Mww+=survM;to.MTot+=survM;Pat[pat].Vww+=surv-survM;to.Vww+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].fJww[TL-1],pcomp); 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mww+=survM; to.Mww+=survM;to.MTot+=survM;Pat[pat].fVww+=surv-survM;to.Vww+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].mJww[TL-1],pcomp); 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mww+=survM; to.Mww+=survM;to.MTot+=survM;Pat[pat].mVww+=surv-survM;to.Vww+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].bJww[TL-1],pcomp); 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mww+=survM; to.Mww+=survM;to.MTot+=survM;Pat[pat].bVww+=surv-survM;to.Vww+=surv-survM;to.VTot+=surv-survM;
		};

		surv=random_binomial(Pat[pat].fJwd[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mwd+=survM; to.Mwd+=survM;to.MTot+=survM;Pat[pat].fVwd+=surv-survM;to.Vwd+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};
		surv=random_binomial(Pat[pat].mJwd[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mwd+=survM; to.Mwd+=survM;to.MTot+=survM;Pat[pat].mVwd+=surv-survM;to.Vwd+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};

		surv=random_binomial(Pat[pat].bJwd[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mwd+=survM; to.Mwd+=survM;to.MTot+=survM;Pat[pat].bVwd+=surv-survM;to.Vwd+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};

		surv=random_binomial(Pat[pat].Jdd[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mdd+=survM; to.Mdd+=survM;to.MTot+=survM;Pat[pat].Vdd+=surv-survM;to.Vdd+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};
		surv=random_binomial(Pat[pat].Jwr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mwr+=survM; to.Mwr+=survM;to.MTot+=survM;Pat[pat].Vwr+=surv-survM;to.Vwr+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distW)to.distW=Pat[pat].distOrigin;
			if(Pat[pat].distOrigin>to.distR)to.distR=Pat[pat].distOrigin;
		};
		surv=random_binomial(Pat[pat].fJwr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mwr+=survM; to.Mwr+=survM;to.MTot+=survM;Pat[pat].fVwr+=surv-survM;to.Vwr+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].mJwr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mwr+=survM; to.Mwr+=survM;to.MTot+=survM;Pat[pat].mVwr+=surv-survM;to.Vwr+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].bJwr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mwr+=survM; to.Mwr+=survM;to.MTot+=survM;Pat[pat].bVwr+=surv-survM;to.Vwr+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].Jrr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mrr+=survM; to.Mrr+=survM;to.MTot+=survM;Pat[pat].Vrr+=surv-survM;to.Vrr+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].Jdr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,0.5); Pat[pat].MTot+=survM;Pat[pat].Mdr+=survM; to.Mdr+=survM;to.MTot+=survM;Pat[pat].Vdr+=surv-survM;to.Vdr+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};
		surv=random_binomial(Pat[pat].mbJrr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mrr+=survM; to.Mrr+=survM;to.MTot+=survM;Pat[pat].Vrr+=surv-survM;to.Vrr+=surv-survM;to.VTot+=surv-survM;
		};
		surv=random_binomial(Pat[pat].mbJdr[TL-1],pcomp); 
		if(surv>0)
		{		
		survM=random_binomial(surv,pa.bias); Pat[pat].MTot+=survM;Pat[pat].Mdr+=survM; to.Mdr+=survM;to.MTot+=survM;Pat[pat].Vdr+=surv-survM;to.Vdr+=surv-survM;to.VTot+=surv-survM;
			if(Pat[pat].distOrigin>to.distD)to.distD=Pat[pat].distOrigin;
		};
		};
		return;};

	void JuvGetOlder(void){
		long int jtPatch,jtAll,jtwwAll,jtwdAll,jtddAll,jtwrAll,jtrrAll,jtdrAll;
		jtAll=0; jtwwAll=0; jtwdAll=0; jtddAll=0; jtwrAll=0; jtrrAll=0; jtdrAll=0;
		double pcomp;
		for(int pat=0;pat<Pat.size();pat++)
		{
		pcomp=(1-pa.muJ)*Pat[pat].comp;
			jtPatch=0;
			for(int age=TL-1;age>0;age--)
			{
				Pat[pat].Jww[age]=random_binomial(Pat[pat].Jww[age-1],pcomp);jtPatch+=Pat[pat].Jww[age];jtwwAll+=Pat[pat].Jww[age];
				Pat[pat].fJww[age]=random_binomial(Pat[pat].fJww[age-1],pcomp);jtPatch+=Pat[pat].fJww[age];jtwwAll+=Pat[pat].fJww[age];
				Pat[pat].mJww[age]=random_binomial(Pat[pat].mJww[age-1],pcomp);jtPatch+=Pat[pat].mJww[age];jtwwAll+=Pat[pat].mJww[age];
				Pat[pat].bJww[age]=random_binomial(Pat[pat].bJww[age-1],pcomp);jtPatch+=Pat[pat].bJww[age];jtwwAll+=Pat[pat].bJww[age];
				Pat[pat].fJwd[age]=random_binomial(Pat[pat].fJwd[age-1],pcomp);jtPatch+=Pat[pat].fJwd[age];jtwdAll+=Pat[pat].fJwd[age];
				Pat[pat].mJwd[age]=random_binomial(Pat[pat].mJwd[age-1],pcomp);jtPatch+=Pat[pat].mJwd[age];jtwdAll+=Pat[pat].mJwd[age];
				Pat[pat].bJwd[age]=random_binomial(Pat[pat].bJwd[age-1],pcomp);jtPatch+=Pat[pat].bJwd[age];jtwdAll+=Pat[pat].bJwd[age];
				Pat[pat].Jdd[age]=random_binomial(Pat[pat].Jdd[age-1],pcomp);jtPatch+=Pat[pat].Jdd[age];jtddAll+=Pat[pat].Jdd[age];
				Pat[pat].Jwr[age]=random_binomial(Pat[pat].Jwr[age-1],pcomp);jtPatch+=Pat[pat].Jwr[age];jtwrAll+=Pat[pat].Jwr[age];
				Pat[pat].fJwr[age]=random_binomial(Pat[pat].fJwr[age-1],pcomp);jtPatch+=Pat[pat].fJwr[age];jtwrAll+=Pat[pat].fJwr[age];
				Pat[pat].mJwr[age]=random_binomial(Pat[pat].mJwr[age-1],pcomp);jtPatch+=Pat[pat].mJwr[age];jtwrAll+=Pat[pat].mJwr[age];
				Pat[pat].bJwr[age]=random_binomial(Pat[pat].bJwr[age-1],pcomp);jtPatch+=Pat[pat].bJwr[age];jtwrAll+=Pat[pat].bJwr[age];
				Pat[pat].Jrr[age]=random_binomial(Pat[pat].Jrr[age-1],pcomp);jtPatch+=Pat[pat].Jrr[age];jtrrAll+=Pat[pat].Jrr[age];
				Pat[pat].Jdr[age]=random_binomial(Pat[pat].Jdr[age-1],pcomp);jtPatch+=Pat[pat].Jdr[age];jtdrAll+=Pat[pat].Jdr[age];
				Pat[pat].mbJrr[age]=random_binomial(Pat[pat].mbJrr[age-1],pcomp);jtPatch+=Pat[pat].mbJrr[age];jtrrAll+=Pat[pat].mbJrr[age];
				Pat[pat].mbJdr[age]=random_binomial(Pat[pat].mbJdr[age-1],pcomp);jtPatch+=Pat[pat].mbJdr[age];jtdrAll+=Pat[pat].mbJdr[age];
			};
			Pat[pat].JTot=jtPatch;
			jtAll+=jtPatch;
		};
		to.Jww=jtwwAll;
		to.Jwd=jtwdAll;
		to.Jdd=jtddAll;
		to.Jwr=jtwrAll;
		to.Jrr=jtrrAll;
		to.Jdr=jtdrAll;
		to.JTot=jtAll;
	return;};

	void VirginsMate(){
		int vww=0;	int fvww=0;	int mvww=0;	int bvww=0;	
		int fvwd=0; int mvwd=0; int bvwd=0; int vdd=0;
		int vwr=0; int fvwr=0; int mvwr=0; int bvwr=0;
		int vrr=0;	int vdr=0;	
		int fww=0;	int ffww=0;	int mfww=0;	int bfww=0;	
		int ffwd=0; int mfwd=0; int bfwd=0; int fdd=0;
		int fwr=0; int ffwr=0; int mfwr=0; int bfwr=0;
		int frr=0;	int fdr=0;	
		int* mates;
		int num;
		for(int pat=0;pat<Pat.size();pat++)
				{
			long long int probs[6]={Pat[pat].Mww,Pat[pat].Mwd,Pat[pat].Mdd,Pat[pat].Mwr,Pat[pat].Mrr,Pat[pat].Mdr};
					vww=random_binomial(Pat[pat].Vww,Pat[pat].mate_rate);
					fvww=random_binomial(Pat[pat].fVww,Pat[pat].mate_rate);
					mvww=random_binomial(Pat[pat].mVww,Pat[pat].mate_rate);
					bvww=random_binomial(Pat[pat].bVww,Pat[pat].mate_rate);
					fvwd=random_binomial(Pat[pat].fVwd,Pat[pat].mate_rate);
					mvwd=random_binomial(Pat[pat].mVwd,Pat[pat].mate_rate);
					bvwd=random_binomial(Pat[pat].bVwd,Pat[pat].mate_rate);
					vdd=random_binomial(Pat[pat].Vdd,Pat[pat].mate_rate);
					vwr=random_binomial(Pat[pat].Vwr,Pat[pat].mate_rate);
					fvwr=random_binomial(Pat[pat].fVwr,Pat[pat].mate_rate);
					mvwr=random_binomial(Pat[pat].mVwr,Pat[pat].mate_rate);
					bvwr=random_binomial(Pat[pat].bVwr,Pat[pat].mate_rate);
					vrr=random_binomial(Pat[pat].Vrr,Pat[pat].mate_rate);
					vdr=random_binomial(Pat[pat].Vdr,Pat[pat].mate_rate);

					if(vww>0)
					{
					mates=random_multinom(vww,probs);
					Pat[pat].Fwwww+=*mates; Pat[pat].Fwwwd+=*(mates+1);Pat[pat].Fwwdd+=*(mates+2);Pat[pat].Fwwwr+=*(mates+3);Pat[pat].Fwwrr+=*(mates+4);Pat[pat].Fwwdr+=*(mates+5);
					to.Fww+=vww;
					delete[] mates;
					Pat[pat].Vww-=vww; to.Vww-=vww; to.VTot-=vww; to.FTot+=vww;
					};

					if(fvww>0)
					{
					mates=random_multinom(fvww,probs);
					Pat[pat].fFwwww+=*mates; Pat[pat].fFwwwd+=*(mates+1);Pat[pat].fFwwdd+=*(mates+2);Pat[pat].fFwwwr+=*(mates+3);Pat[pat].fFwwrr+=*(mates+4);Pat[pat].fFwwdr+=*(mates+5);
					to.Fww+=fvww;
					delete[] mates;
					Pat[pat].fVww-=fvww; to.Vww-=fvww; to.VTot-=fvww; to.FTot+=fvww;
					};

					if(mvww>0)
					{
					mates=random_multinom(mvww,probs);
					Pat[pat].mFwwww+=*mates; Pat[pat].mFwwwd+=*(mates+1);Pat[pat].mFwwdd+=*(mates+2);Pat[pat].mFwwwr+=*(mates+3);Pat[pat].mFwwrr+=*(mates+4);Pat[pat].mFwwdr+=*(mates+5);
					to.Fww+=mvww;
					delete[] mates;
					Pat[pat].mVww-=mvww; to.Vww-=mvww; to.VTot-=mvww; to.FTot+=mvww;
					};

					if(bvww>0)
					{
					mates=random_multinom(bvww,probs);
					Pat[pat].bFwwww+=*mates; Pat[pat].bFwwwd+=*(mates+1);Pat[pat].bFwwdd+=*(mates+2);Pat[pat].bFwwwr+=*(mates+3);Pat[pat].bFwwrr+=*(mates+4);Pat[pat].bFwwdr+=*(mates+5);
					to.Fww+=bvww;
					delete[] mates;
					Pat[pat].bVww-=bvww; to.Vww-=bvww; to.VTot-=bvww; to.FTot+=bvww;
					};
					if(fvwd>0)
					{
					mates=random_multinom(fvwd,probs);
					Pat[pat].fFwdww+=*mates; Pat[pat].fFwdwd+=*(mates+1);Pat[pat].fFwddd+=*(mates+2);Pat[pat].fFwdwr+=*(mates+3);Pat[pat].fFwdrr+=*(mates+4);Pat[pat].fFwddr+=*(mates+5);
					to.Fwd+=fvwd;
					delete[] mates;
					Pat[pat].fVwd-=fvwd; to.Vwd-=fvwd; to.VTot-=fvwd; to.FTot+=fvwd;
					};
					if(mvwd>0)
					{
					mates=random_multinom(mvwd,probs);
					Pat[pat].mFwdww+=*mates; Pat[pat].mFwdwd+=*(mates+1);Pat[pat].mFwddd+=*(mates+2);Pat[pat].mFwdwr+=*(mates+3);Pat[pat].mFwdrr+=*(mates+4);Pat[pat].mFwddr+=*(mates+5);
					to.Fwd+=mvwd;
					delete[] mates;
					Pat[pat].mVwd-=mvwd; to.Vwd-=mvwd; to.VTot-=mvwd; to.FTot+=mvwd;
					};
					if(bvwd>0)
					{
					mates=random_multinom(bvwd,probs);
					Pat[pat].bFwdww+=*mates; Pat[pat].bFwdwd+=*(mates+1);Pat[pat].bFwddd+=*(mates+2);Pat[pat].bFwdwr+=*(mates+3);Pat[pat].bFwdrr+=*(mates+4);Pat[pat].bFwddr+=*(mates+5);
					to.Fwd+=bvwd;
					delete[] mates;
					Pat[pat].bVwd-=bvwd; to.Vwd-=bvwd; to.VTot-=bvwd; to.FTot+=bvwd;
					};
					if(vdd>0)
					{
					mates=random_multinom(vdd,probs);
					Pat[pat].Fddww+=*mates; Pat[pat].Fddwd+=*(mates+1);Pat[pat].Fdddd+=*(mates+2);Pat[pat].Fddwr+=*(mates+3);Pat[pat].Fddrr+=*(mates+4);Pat[pat].Fdddr+=*(mates+5);
					to.Fdd+=vdd;
					delete[] mates;
					Pat[pat].Vdd-=vdd; to.Vdd-=vdd; to.VTot-=vdd; to.FTot+=vdd;
					};
					if(vwr>0)
					{
					mates=random_multinom(vwr,probs);
					Pat[pat].Fwrww+=*mates; Pat[pat].Fwrwd+=*(mates+1);Pat[pat].Fwrdd+=*(mates+2);Pat[pat].Fwrwr+=*(mates+3);Pat[pat].Fwrrr+=*(mates+4);Pat[pat].Fwrdr+=*(mates+5);
					to.Fwr+=vwr;
					delete[] mates;
					Pat[pat].Vwr-=vwr; to.Vwr-=vwr; to.VTot-=vwr; to.FTot+=vwr;
					};
					if(fvwr>0)
					{
					mates=random_multinom(fvwr,probs);
					Pat[pat].fFwrww+=*mates; Pat[pat].fFwrwd+=*(mates+1);Pat[pat].fFwrdd+=*(mates+2);Pat[pat].fFwrwr+=*(mates+3);Pat[pat].fFwrrr+=*(mates+4);Pat[pat].fFwrdr+=*(mates+5);
					to.Fwr+=fvwr;
					delete[] mates;
					Pat[pat].fVwr-=fvwr; to.Vwr-=fvwr; to.VTot-=fvwr; to.FTot+=fvwr;
					};
					if(mvwr>0)
					{
					mates=random_multinom(mvwr,probs);
					Pat[pat].mFwrww+=*mates; Pat[pat].mFwrwd+=*(mates+1);Pat[pat].mFwrdd+=*(mates+2);Pat[pat].mFwrwr+=*(mates+3);Pat[pat].mFwrrr+=*(mates+4);Pat[pat].mFwrdr+=*(mates+5);
					to.Fwr+=mvwr;
					delete[] mates;
					Pat[pat].mVwr-=mvwr; to.Vwr-=mvwr; to.VTot-=mvwr; to.FTot+=mvwr;
					};
					if(bvwr>0)
					{
					mates=random_multinom(bvwr,probs);
					Pat[pat].bFwrww+=*mates; Pat[pat].bFwrwd+=*(mates+1);Pat[pat].bFwrdd+=*(mates+2);Pat[pat].bFwrwr+=*(mates+3);Pat[pat].bFwrrr+=*(mates+4);Pat[pat].bFwrdr+=*(mates+5);
					to.Fwr+=bvwr;
					delete[] mates;
					Pat[pat].bVwr-=bvwr; to.Vwr-=bvwr; to.VTot-=bvwr; to.FTot+=bvwr;
					};
					if(vrr>0)
					{
					mates=random_multinom(vrr,probs);
					Pat[pat].Frrww+=*mates; Pat[pat].Frrwd+=*(mates+1);Pat[pat].Frrdd+=*(mates+2);Pat[pat].Frrwr+=*(mates+3);Pat[pat].Frrrr+=*(mates+4);Pat[pat].Frrdr+=*(mates+5);
					to.Frr+=vrr;
					delete[] mates;
					Pat[pat].Vrr-=vrr; to.Vrr-=vrr; to.VTot-=vrr; to.FTot+=vrr;
					};
					if(vdr>0)
					{
					mates=random_multinom(vdr,probs);
					Pat[pat].Fdrww+=*mates; Pat[pat].Fdrwd+=*(mates+1);Pat[pat].Fdrdd+=*(mates+2);Pat[pat].Fdrwr+=*(mates+3);Pat[pat].Fdrrr+=*(mates+4);Pat[pat].Fdrdr+=*(mates+5);
					to.Fdr+=vdr;
					delete[] mates;
					Pat[pat].Vdr-=vdr; to.Vdr-=vdr; to.VTot-=vdr; to.FTot+=vdr;
					};
				};
	return;};






	void AdultsMove(void){
		int pat,newpat;
		int* aww;
		int howmany;
		if(Pat.size()>1)
		{
		for(pat=0;pat<Pat.size();pat++)
				{
		Pat[pat].MoveMww=random_binomial(Pat[pat].Mww,pa.d);Pat[pat].Mww-=Pat[pat].MoveMww;Pat[pat].MTot-=Pat[pat].MoveMww;
		Pat[pat].MoveMwd=random_binomial(Pat[pat].Mwd,pa.d);Pat[pat].Mwd-=Pat[pat].MoveMwd;Pat[pat].MTot-=Pat[pat].MoveMwd;
		Pat[pat].MoveMdd=random_binomial(Pat[pat].Mdd,pa.d);Pat[pat].Mdd-=Pat[pat].MoveMdd;Pat[pat].MTot-=Pat[pat].MoveMdd;
		Pat[pat].MoveMwr=random_binomial(Pat[pat].Mwr,pa.d);Pat[pat].Mwr-=Pat[pat].MoveMwr;Pat[pat].MTot-=Pat[pat].MoveMwr;
		Pat[pat].MoveMrr=random_binomial(Pat[pat].Mrr,pa.d);Pat[pat].Mrr-=Pat[pat].MoveMrr;Pat[pat].MTot-=Pat[pat].MoveMrr;
		Pat[pat].MoveMdr=random_binomial(Pat[pat].Mdr,pa.d);Pat[pat].Mdr-=Pat[pat].MoveMdr;Pat[pat].MTot-=Pat[pat].MoveMdr;

		Pat[pat].MoveFwwww=random_binomial(Pat[pat].Fwwww,pa.d);Pat[pat].Fwwww-=Pat[pat].MoveFwwww; 
		Pat[pat].MoveFwwwd=random_binomial(Pat[pat].Fwwwd,pa.d);Pat[pat].Fwwwd-=Pat[pat].MoveFwwwd; 
		Pat[pat].MoveFwwdd=random_binomial(Pat[pat].Fwwdd,pa.d);Pat[pat].Fwwdd-=Pat[pat].MoveFwwdd; 
		Pat[pat].MoveFwwwr=random_binomial(Pat[pat].Fwwwr,pa.d);Pat[pat].Fwwwr-=Pat[pat].MoveFwwwr; 
		Pat[pat].MoveFwwrr=random_binomial(Pat[pat].Fwwrr,pa.d);Pat[pat].Fwwrr-=Pat[pat].MoveFwwrr; 
		Pat[pat].MoveFwwdr=random_binomial(Pat[pat].Fwwdr,pa.d);Pat[pat].Fwwdr-=Pat[pat].MoveFwwdr; 
		Pat[pat].MovefFwwww=random_binomial(Pat[pat].fFwwww,pa.d);Pat[pat].fFwwww-=Pat[pat].MovefFwwww; 
		Pat[pat].MovefFwwwd=random_binomial(Pat[pat].fFwwwd,pa.d);Pat[pat].fFwwwd-=Pat[pat].MovefFwwwd; 
		Pat[pat].MovefFwwdd=random_binomial(Pat[pat].fFwwdd,pa.d);Pat[pat].fFwwdd-=Pat[pat].MovefFwwdd; 
		Pat[pat].MovefFwwwr=random_binomial(Pat[pat].fFwwwr,pa.d);Pat[pat].fFwwwr-=Pat[pat].MovefFwwwr; 
		Pat[pat].MovefFwwrr=random_binomial(Pat[pat].fFwwrr,pa.d);Pat[pat].fFwwrr-=Pat[pat].MovefFwwrr; 
		Pat[pat].MovefFwwdr=random_binomial(Pat[pat].fFwwdr,pa.d);Pat[pat].fFwwdr-=Pat[pat].MovefFwwdr; 
		Pat[pat].MovemFwwww=random_binomial(Pat[pat].mFwwww,pa.d);Pat[pat].mFwwww-=Pat[pat].MovemFwwww; 
		Pat[pat].MovemFwwwd=random_binomial(Pat[pat].mFwwwd,pa.d);Pat[pat].mFwwwd-=Pat[pat].MovemFwwwd; 
		Pat[pat].MovemFwwdd=random_binomial(Pat[pat].mFwwdd,pa.d);Pat[pat].mFwwdd-=Pat[pat].MovemFwwdd; 
		Pat[pat].MovemFwwwr=random_binomial(Pat[pat].mFwwwr,pa.d);Pat[pat].mFwwwr-=Pat[pat].MovemFwwwr; 
		Pat[pat].MovemFwwrr=random_binomial(Pat[pat].mFwwrr,pa.d);Pat[pat].mFwwrr-=Pat[pat].MovemFwwrr; 
		Pat[pat].MovemFwwdr=random_binomial(Pat[pat].mFwwdr,pa.d);Pat[pat].mFwwdr-=Pat[pat].MovemFwwdr; 
		Pat[pat].MovebFwwww=random_binomial(Pat[pat].bFwwww,pa.d);Pat[pat].bFwwww-=Pat[pat].MovebFwwww; 
		Pat[pat].MovebFwwwd=random_binomial(Pat[pat].bFwwwd,pa.d);Pat[pat].bFwwwd-=Pat[pat].MovebFwwwd; 
		Pat[pat].MovebFwwdd=random_binomial(Pat[pat].bFwwdd,pa.d);Pat[pat].bFwwdd-=Pat[pat].MovebFwwdd; 
		Pat[pat].MovebFwwwr=random_binomial(Pat[pat].bFwwwr,pa.d);Pat[pat].bFwwwr-=Pat[pat].MovebFwwwr; 
		Pat[pat].MovebFwwrr=random_binomial(Pat[pat].bFwwrr,pa.d);Pat[pat].bFwwrr-=Pat[pat].MovebFwwrr; 
		Pat[pat].MovebFwwdr=random_binomial(Pat[pat].bFwwdr,pa.d);Pat[pat].bFwwdr-=Pat[pat].MovebFwwdr; 

		Pat[pat].MovefFwdww=random_binomial(Pat[pat].fFwdww,pa.d);Pat[pat].fFwdww-=Pat[pat].MovefFwdww; 
		Pat[pat].MovefFwdwd=random_binomial(Pat[pat].fFwdwd,pa.d);Pat[pat].fFwdwd-=Pat[pat].MovefFwdwd; 
		Pat[pat].MovefFwddd=random_binomial(Pat[pat].fFwddd,pa.d);Pat[pat].fFwddd-=Pat[pat].MovefFwddd; 
		Pat[pat].MovefFwdwr=random_binomial(Pat[pat].fFwdwr,pa.d);Pat[pat].fFwdwr-=Pat[pat].MovefFwdwr; 
		Pat[pat].MovefFwdrr=random_binomial(Pat[pat].fFwdrr,pa.d);Pat[pat].fFwdrr-=Pat[pat].MovefFwdrr; 
		Pat[pat].MovefFwddr=random_binomial(Pat[pat].fFwddr,pa.d);Pat[pat].fFwddr-=Pat[pat].MovefFwddr; 
		Pat[pat].MovemFwdww=random_binomial(Pat[pat].mFwdww,pa.d);Pat[pat].mFwdww-=Pat[pat].MovemFwdww; 
		Pat[pat].MovemFwdwd=random_binomial(Pat[pat].mFwdwd,pa.d);Pat[pat].mFwdwd-=Pat[pat].MovemFwdwd; 
		Pat[pat].MovemFwddd=random_binomial(Pat[pat].mFwddd,pa.d);Pat[pat].mFwddd-=Pat[pat].MovemFwddd; 
		Pat[pat].MovemFwdwr=random_binomial(Pat[pat].mFwdwr,pa.d);Pat[pat].mFwdwr-=Pat[pat].MovemFwdwr; 
		Pat[pat].MovemFwdrr=random_binomial(Pat[pat].mFwdrr,pa.d);Pat[pat].mFwdrr-=Pat[pat].MovemFwdrr; 
		Pat[pat].MovemFwddr=random_binomial(Pat[pat].mFwddr,pa.d);Pat[pat].mFwddr-=Pat[pat].MovemFwddr; 
		Pat[pat].MovebFwdww=random_binomial(Pat[pat].bFwdww,pa.d);Pat[pat].bFwdww-=Pat[pat].MovebFwdww; 
		Pat[pat].MovebFwdwd=random_binomial(Pat[pat].bFwdwd,pa.d);Pat[pat].bFwdwd-=Pat[pat].MovebFwdwd; 
		Pat[pat].MovebFwddd=random_binomial(Pat[pat].bFwddd,pa.d);Pat[pat].bFwddd-=Pat[pat].MovebFwddd; 
		Pat[pat].MovebFwdwr=random_binomial(Pat[pat].bFwdwr,pa.d);Pat[pat].bFwdwr-=Pat[pat].MovebFwdwr; 
		Pat[pat].MovebFwdrr=random_binomial(Pat[pat].bFwdrr,pa.d);Pat[pat].bFwdrr-=Pat[pat].MovebFwdrr; 
		Pat[pat].MovebFwddr=random_binomial(Pat[pat].bFwddr,pa.d);Pat[pat].bFwddr-=Pat[pat].MovebFwddr; 
		Pat[pat].MoveFddww=random_binomial(Pat[pat].Fddww,pa.d);Pat[pat].Fddww-=Pat[pat].MoveFddww; 
		Pat[pat].MoveFddwd=random_binomial(Pat[pat].Fddwd,pa.d);Pat[pat].Fddwd-=Pat[pat].MoveFddwd; 
		Pat[pat].MoveFdddd=random_binomial(Pat[pat].Fdddd,pa.d);Pat[pat].Fdddd-=Pat[pat].MoveFdddd; 
		Pat[pat].MoveFddwr=random_binomial(Pat[pat].Fddwr,pa.d);Pat[pat].Fddwr-=Pat[pat].MoveFddwr; 
		Pat[pat].MoveFddrr=random_binomial(Pat[pat].Fddrr,pa.d);Pat[pat].Fddrr-=Pat[pat].MoveFddrr; 
		Pat[pat].MoveFdddr=random_binomial(Pat[pat].Fdddr,pa.d);Pat[pat].Fdddr-=Pat[pat].MoveFdddr; 
		Pat[pat].MoveFwrww=random_binomial(Pat[pat].Fwrww,pa.d);Pat[pat].Fwrww-=Pat[pat].MoveFwrww; 
		Pat[pat].MoveFwrwd=random_binomial(Pat[pat].Fwrwd,pa.d);Pat[pat].Fwrwd-=Pat[pat].MoveFwrwd; 
		Pat[pat].MoveFwrdd=random_binomial(Pat[pat].Fwrdd,pa.d);Pat[pat].Fwrdd-=Pat[pat].MoveFwrdd; 
		Pat[pat].MoveFwrwr=random_binomial(Pat[pat].Fwrwr,pa.d);Pat[pat].Fwrwr-=Pat[pat].MoveFwrwr; 
		Pat[pat].MoveFwrrr=random_binomial(Pat[pat].Fwrrr,pa.d);Pat[pat].Fwrrr-=Pat[pat].MoveFwrrr; 
		Pat[pat].MoveFwrdr=random_binomial(Pat[pat].Fwrdr,pa.d);Pat[pat].Fwrdr-=Pat[pat].MoveFwrdr; 
		Pat[pat].MovefFwrww=random_binomial(Pat[pat].fFwrww,pa.d);Pat[pat].fFwrww-=Pat[pat].MovefFwrww; 
		Pat[pat].MovefFwrwd=random_binomial(Pat[pat].fFwrwd,pa.d);Pat[pat].fFwrwd-=Pat[pat].MovefFwrwd; 
		Pat[pat].MovefFwrdd=random_binomial(Pat[pat].fFwrdd,pa.d);Pat[pat].fFwrdd-=Pat[pat].MovefFwrdd; 
		Pat[pat].MovefFwrwr=random_binomial(Pat[pat].fFwrwr,pa.d);Pat[pat].fFwrwr-=Pat[pat].MovefFwrwr; 
		Pat[pat].MovefFwrrr=random_binomial(Pat[pat].fFwrrr,pa.d);Pat[pat].fFwrrr-=Pat[pat].MovefFwrrr; 
		Pat[pat].MovefFwrdr=random_binomial(Pat[pat].fFwrdr,pa.d);Pat[pat].fFwrdr-=Pat[pat].MovefFwrdr; 
		Pat[pat].MovemFwrww=random_binomial(Pat[pat].mFwrww,pa.d);Pat[pat].mFwrww-=Pat[pat].MovemFwrww; 
		Pat[pat].MovemFwrwd=random_binomial(Pat[pat].mFwrwd,pa.d);Pat[pat].mFwrwd-=Pat[pat].MovemFwrwd; 
		Pat[pat].MovemFwrdd=random_binomial(Pat[pat].mFwrdd,pa.d);Pat[pat].mFwrdd-=Pat[pat].MovemFwrdd; 
		Pat[pat].MovemFwrwr=random_binomial(Pat[pat].mFwrwr,pa.d);Pat[pat].mFwrwr-=Pat[pat].MovemFwrwr; 
		Pat[pat].MovemFwrrr=random_binomial(Pat[pat].mFwrrr,pa.d);Pat[pat].mFwrrr-=Pat[pat].MovemFwrrr; 
		Pat[pat].MovemFwrdr=random_binomial(Pat[pat].mFwrdr,pa.d);Pat[pat].mFwrdr-=Pat[pat].MovemFwrdr; 
		Pat[pat].MovebFwrww=random_binomial(Pat[pat].bFwrww,pa.d);Pat[pat].bFwrww-=Pat[pat].MovebFwrww; 
		Pat[pat].MovebFwrwd=random_binomial(Pat[pat].bFwrwd,pa.d);Pat[pat].bFwrwd-=Pat[pat].MovebFwrwd; 
		Pat[pat].MovebFwrdd=random_binomial(Pat[pat].bFwrdd,pa.d);Pat[pat].bFwrdd-=Pat[pat].MovebFwrdd; 
		Pat[pat].MovebFwrwr=random_binomial(Pat[pat].bFwrwr,pa.d);Pat[pat].bFwrwr-=Pat[pat].MovebFwrwr; 
		Pat[pat].MovebFwrrr=random_binomial(Pat[pat].bFwrrr,pa.d);Pat[pat].bFwrrr-=Pat[pat].MovebFwrrr; 
		Pat[pat].MovebFwrdr=random_binomial(Pat[pat].bFwrdr,pa.d);Pat[pat].bFwrdr-=Pat[pat].MovebFwrdr; 
		Pat[pat].MoveFrrww=random_binomial(Pat[pat].Frrww,pa.d);Pat[pat].Frrww-=Pat[pat].MoveFrrww; 
		Pat[pat].MoveFrrwd=random_binomial(Pat[pat].Frrwd,pa.d);Pat[pat].Frrwd-=Pat[pat].MoveFrrwd; 
		Pat[pat].MoveFrrdd=random_binomial(Pat[pat].Frrdd,pa.d);Pat[pat].Frrdd-=Pat[pat].MoveFrrdd; 
		Pat[pat].MoveFrrwr=random_binomial(Pat[pat].Frrwr,pa.d);Pat[pat].Frrwr-=Pat[pat].MoveFrrwr; 
		Pat[pat].MoveFrrrr=random_binomial(Pat[pat].Frrrr,pa.d);Pat[pat].Frrrr-=Pat[pat].MoveFrrrr; 
		Pat[pat].MoveFrrdr=random_binomial(Pat[pat].Frrdr,pa.d);Pat[pat].Frrdr-=Pat[pat].MoveFrrdr; 
		Pat[pat].MoveFdrww=random_binomial(Pat[pat].Fdrww,pa.d);Pat[pat].Fdrww-=Pat[pat].MoveFdrww; 
		Pat[pat].MoveFdrwd=random_binomial(Pat[pat].Fdrwd,pa.d);Pat[pat].Fdrwd-=Pat[pat].MoveFdrwd; 
		Pat[pat].MoveFdrdd=random_binomial(Pat[pat].Fdrdd,pa.d);Pat[pat].Fdrdd-=Pat[pat].MoveFdrdd; 
		Pat[pat].MoveFdrwr=random_binomial(Pat[pat].Fdrwr,pa.d);Pat[pat].Fdrwr-=Pat[pat].MoveFdrwr; 
		Pat[pat].MoveFdrrr=random_binomial(Pat[pat].Fdrrr,pa.d);Pat[pat].Fdrrr-=Pat[pat].MoveFdrrr; 
		Pat[pat].MoveFdrdr=random_binomial(Pat[pat].Fdrdr,pa.d);Pat[pat].Fdrdr-=Pat[pat].MoveFdrdr; 
				};
		for(pat=0;pat<Pat.size();pat++)
				{
				howmany=Pat[pat].connecIND.size();

				if(Pat[pat].MoveMww>0)
				{
				aww=random_multinom_var(Pat[pat].MoveMww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
				for(newpat=howmany-1;newpat>=0;newpat--) {Pat[Pat[pat].connecIND[newpat]].Mww+=*(aww+newpat);Pat[Pat[pat].connecIND[newpat]].MTot+=*(aww+newpat);};
				delete[] aww;
				};
				if(Pat[pat].MoveMwd>0)
				{
				aww=random_multinom_var(Pat[pat].MoveMwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
				for(newpat=howmany-1;newpat>=0;newpat--) {Pat[Pat[pat].connecIND[newpat]].Mwd+=*(aww+newpat);Pat[Pat[pat].connecIND[newpat]].MTot+=*(aww+newpat);};
				delete[] aww;
				};
				if(Pat[pat].MoveMdd>0)
				{
				aww=random_multinom_var(Pat[pat].MoveMdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
				for(newpat=howmany-1;newpat>=0;newpat--) {Pat[Pat[pat].connecIND[newpat]].Mdd+=*(aww+newpat);Pat[Pat[pat].connecIND[newpat]].MTot+=*(aww+newpat);};
				delete[] aww;
				};
				if(Pat[pat].MoveMwr>0)
				{
				aww=random_multinom_var(Pat[pat].MoveMwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
				for(newpat=howmany-1;newpat>=0;newpat--) {Pat[Pat[pat].connecIND[newpat]].Mwr+=*(aww+newpat);Pat[Pat[pat].connecIND[newpat]].MTot+=*(aww+newpat);};
				delete[] aww;
				};
				if(Pat[pat].MoveMrr>0)
				{
				aww=random_multinom_var(Pat[pat].MoveMrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
				for(newpat=howmany-1;newpat>=0;newpat--) {Pat[Pat[pat].connecIND[newpat]].Mrr+=*(aww+newpat);Pat[Pat[pat].connecIND[newpat]].MTot+=*(aww+newpat);};
				delete[] aww;
				};
				if(Pat[pat].MoveMdr>0)
				{
				aww=random_multinom_var(Pat[pat].MoveMdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
				for(newpat=howmany-1;newpat>=0;newpat--) {Pat[Pat[pat].connecIND[newpat]].Mdr+=*(aww+newpat);Pat[Pat[pat].connecIND[newpat]].MTot+=*(aww+newpat);};
				delete[] aww;
				};

		if(Pat[pat].MoveFwwww>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwwww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwwww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwwwd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwwwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwwwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwwdd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwwdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwwdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwwwr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwwwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwwwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwwrr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwwrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwwrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwwdr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwwdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwwdr+=*(aww+newpat);
		delete[] aww;
		};

		if(Pat[pat].MovefFwwww>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwwww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwwww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwwwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwwwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwwwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwwdd>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwwdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwwdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwwwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwwwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwwwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwwrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwwrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwwrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwwdr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwwdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwwdr+=*(aww+newpat);
		delete[] aww;
		};

		if(Pat[pat].MovemFwwww>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwwww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwwww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwwwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwwwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwwwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwwdd>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwwdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwwdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwwwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwwwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwwwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwwrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwwrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwwrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwwdr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwwdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwwdr+=*(aww+newpat);
		delete[] aww;
		};

		if(Pat[pat].MovebFwwww>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwwww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwwww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwwwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwwwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwwwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwwdd>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwwdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwwdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwwwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwwwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwwwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwwrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwwrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwwrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwwdr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwwdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwwdr+=*(aww+newpat);
		delete[] aww;
		};


		if(Pat[pat].MovefFwdww>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwdww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwdww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwdwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwdwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwdwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwddd>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwddd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwddd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwdwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwdwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwdwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwdrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwdrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwdrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwddr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwddr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwddr+=*(aww+newpat);
		delete[] aww;
		};

		if(Pat[pat].MovemFwdww>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwdww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwdww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwdwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwdwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwdwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwddd>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwddd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwddd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwdwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwdwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwdwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwdrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwdrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwdrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwddr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwddr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwddr+=*(aww+newpat);
		delete[] aww;
		};

		if(Pat[pat].MovebFwdww>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwdww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwdww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwdwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwdwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwdwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwddd>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwddd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwddd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwdwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwdwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwdwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwdrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwdrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwdrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwddr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwddr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwddr+=*(aww+newpat);
		delete[] aww;
		};



		if(Pat[pat].MoveFddww>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFddww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fddww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFddwd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFddwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fddwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdddd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdddd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdddd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFddwr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFddwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fddwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFddrr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFddrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fddrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdddr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdddr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdddr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwrww>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwrww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwrww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwrwd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwrwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwrwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwrdd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwrdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwrdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwrwr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwrwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwrwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwrrr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwrrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwrrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFwrdr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFwrdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fwrdr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwrww>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwrww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwrww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwrwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwrwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwrwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwrdd>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwrdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwrdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwrwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwrwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwrwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwrrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwrrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwrrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovefFwrdr>0)
		{
		aww=random_multinom_var(Pat[pat].MovefFwrdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].fFwrdr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwrww>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwrww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwrww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwrwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwrwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwrwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwrdd>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwrdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwrdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwrwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwrwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwrwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwrrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwrrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwrrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovemFwrdr>0)
		{
		aww=random_multinom_var(Pat[pat].MovemFwrdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].mFwrdr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwrww>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwrww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwrww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwrwd>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwrwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwrwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwrdd>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwrdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwrdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwrwr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwrwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwrwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwrrr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwrrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwrrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MovebFwrdr>0)
		{
		aww=random_multinom_var(Pat[pat].MovebFwrdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].bFwrdr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFrrww>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFrrww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Frrww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFrrwd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFrrwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Frrwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFrrdd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFrrdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Frrdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFrrwr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFrrwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Frrwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFrrrr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFrrrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Frrrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFrrdr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFrrdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Frrdr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdrww>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdrww,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdrww+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdrwd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdrwd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdrwd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdrdd>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdrdd,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdrdd+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdrwr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdrwr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdrwr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdrrr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdrrr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdrrr+=*(aww+newpat);
		delete[] aww;
		};
		if(Pat[pat].MoveFdrdr>0)
		{
		aww=random_multinom_var(Pat[pat].MoveFdrdr,howmany,&Pat[pat].connecW[0],Pat[pat].TotW);
		for(newpat=howmany-1;newpat>=0;newpat--) Pat[Pat[pat].connecIND[newpat]].Fdrdr+=*(aww+newpat);
		delete[] aww;
		};
		};
		};
	return;};




	void Hide(){
		int num;
		for(int pat=0;pat<Pat.size();pat++)
				{
	num=random_binomial(Pat[pat].Fwwww,pa.psi);Pat[pat].Fwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesFwwww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwwwd,pa.psi);Pat[pat].Fwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesFwwwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwwdd,pa.psi);Pat[pat].Fwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesFwwdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwwwr,pa.psi);Pat[pat].Fwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesFwwwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwwrr,pa.psi);Pat[pat].Fwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesFwwrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwwdr,pa.psi);Pat[pat].Fwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesFwwdr+=random_binomial(num,1-pa.muAES);	


	num=random_binomial(Pat[pat].fFwwww,pa.psi);Pat[pat].fFwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesfFwwww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwwwd,pa.psi);Pat[pat].fFwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesfFwwwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwwdd,pa.psi);Pat[pat].fFwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesfFwwdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwwwr,pa.psi);Pat[pat].fFwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesfFwwwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwwrr,pa.psi);Pat[pat].fFwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesfFwwrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwwdr,pa.psi);Pat[pat].fFwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesfFwwdr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].mFwwww,pa.psi);Pat[pat].mFwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesmFwwww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwwwd,pa.psi);Pat[pat].mFwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesmFwwwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwwdd,pa.psi);Pat[pat].mFwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesmFwwdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwwwr,pa.psi);Pat[pat].mFwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesmFwwwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwwrr,pa.psi);Pat[pat].mFwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesmFwwrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwwdr,pa.psi);Pat[pat].mFwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesmFwwdr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].bFwwww,pa.psi);Pat[pat].bFwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesbFwwww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwwwd,pa.psi);Pat[pat].bFwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesbFwwwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwwdd,pa.psi);Pat[pat].bFwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesbFwwdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwwwr,pa.psi);Pat[pat].bFwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesbFwwwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwwrr,pa.psi);Pat[pat].bFwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesbFwwrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwwdr,pa.psi);Pat[pat].bFwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].AesbFwwdr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].fFwdww,pa.psi);Pat[pat].fFwdww-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesfFwdww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwdwd,pa.psi);Pat[pat].fFwdwd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesfFwdwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwddd,pa.psi);Pat[pat].fFwddd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesfFwddd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwdwr,pa.psi);Pat[pat].fFwdwr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesfFwdwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwdrr,pa.psi);Pat[pat].fFwdrr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesfFwdrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwddr,pa.psi);Pat[pat].fFwddr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesfFwddr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].mFwdww,pa.psi);Pat[pat].mFwdww-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesmFwdww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwdwd,pa.psi);Pat[pat].mFwdwd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesmFwdwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwddd,pa.psi);Pat[pat].mFwddd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesmFwddd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwdwr,pa.psi);Pat[pat].mFwdwr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesmFwdwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwdrr,pa.psi);Pat[pat].mFwdrr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesmFwdrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwddr,pa.psi);Pat[pat].mFwddr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesmFwddr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].bFwdww,pa.psi);Pat[pat].bFwdww-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesbFwdww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwdwd,pa.psi);Pat[pat].bFwdwd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesbFwdwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwddd,pa.psi);Pat[pat].bFwddd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesbFwddd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwdwr,pa.psi);Pat[pat].bFwdwr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesbFwdwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwdrr,pa.psi);Pat[pat].bFwdrr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesbFwdrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwddr,pa.psi);Pat[pat].bFwddr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].AesbFwddr+=random_binomial(num,1-pa.muAES);	


	num=random_binomial(Pat[pat].Fddww,pa.psi);Pat[pat].Fddww-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].AesFddww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fddwd,pa.psi);Pat[pat].Fddwd-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].AesFddwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdddd,pa.psi);Pat[pat].Fdddd-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].AesFdddd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fddwr,pa.psi);Pat[pat].Fddwr-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].AesFddwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fddrr,pa.psi);Pat[pat].Fddrr-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].AesFddrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdddr,pa.psi);Pat[pat].Fdddr-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].AesFdddr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].Fwrww,pa.psi);Pat[pat].Fwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesFwrww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwrwd,pa.psi);Pat[pat].Fwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesFwrwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwrdd,pa.psi);Pat[pat].Fwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesFwrdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwrwr,pa.psi);Pat[pat].Fwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesFwrwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwrrr,pa.psi);Pat[pat].Fwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesFwrrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fwrdr,pa.psi);Pat[pat].Fwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesFwrdr+=random_binomial(num,1-pa.muAES);	


	num=random_binomial(Pat[pat].fFwrww,pa.psi);Pat[pat].fFwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesfFwrww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwrwd,pa.psi);Pat[pat].fFwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesfFwrwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwrdd,pa.psi);Pat[pat].fFwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesfFwrdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwrwr,pa.psi);Pat[pat].fFwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesfFwrwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwrrr,pa.psi);Pat[pat].fFwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesfFwrrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].fFwrdr,pa.psi);Pat[pat].fFwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesfFwrdr+=random_binomial(num,1-pa.muAES);	

	num=random_binomial(Pat[pat].mFwrww,pa.psi);Pat[pat].mFwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesmFwrww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwrwd,pa.psi);Pat[pat].mFwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesmFwrwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwrdd,pa.psi);Pat[pat].mFwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesmFwrdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwrwr,pa.psi);Pat[pat].mFwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesmFwrwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwrrr,pa.psi);Pat[pat].mFwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesmFwrrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].mFwrdr,pa.psi);Pat[pat].mFwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesmFwrdr+=random_binomial(num,1-pa.muAES);	


	num=random_binomial(Pat[pat].bFwrww,pa.psi);Pat[pat].bFwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesbFwrww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwrwd,pa.psi);Pat[pat].bFwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesbFwrwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwrdd,pa.psi);Pat[pat].bFwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesbFwrdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwrwr,pa.psi);Pat[pat].bFwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesbFwrwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwrrr,pa.psi);Pat[pat].bFwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesbFwrrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].bFwrdr,pa.psi);Pat[pat].bFwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].AesbFwrdr+=random_binomial(num,1-pa.muAES);	





	num=random_binomial(Pat[pat].Frrww,pa.psi);Pat[pat].Frrww-=num;to.Frr-=num;to.FTot-=num;Pat[pat].AesFrrww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Frrwd,pa.psi);Pat[pat].Frrwd-=num;to.Frr-=num;to.FTot-=num;Pat[pat].AesFrrwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Frrdd,pa.psi);Pat[pat].Frrdd-=num;to.Frr-=num;to.FTot-=num;Pat[pat].AesFrrdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Frrwr,pa.psi);Pat[pat].Frrwr-=num;to.Frr-=num;to.FTot-=num;Pat[pat].AesFrrwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Frrrr,pa.psi);Pat[pat].Frrrr-=num;to.Frr-=num;to.FTot-=num;Pat[pat].AesFrrrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Frrdr,pa.psi);Pat[pat].Frrdr-=num;to.Frr-=num;to.FTot-=num;Pat[pat].AesFrrdr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdrww,pa.psi);Pat[pat].Fdrww-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].AesFdrww+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdrwd,pa.psi);Pat[pat].Fdrwd-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].AesFdrwd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdrdd,pa.psi);Pat[pat].Fdrdd-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].AesFdrdd+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdrwr,pa.psi);Pat[pat].Fdrwr-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].AesFdrwr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdrrr,pa.psi);Pat[pat].Fdrrr-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].AesFdrrr+=random_binomial(num,1-pa.muAES);	
	num=random_binomial(Pat[pat].Fdrdr,pa.psi);Pat[pat].Fdrdr-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].AesFdrdr+=random_binomial(num,1-pa.muAES);	
				};
		return;};

	void Wake(int day){
		int num;
		double prob=1.0/(1.0+pa.t_wake2-(day%365));
		for(int pat=0;pat<Pat.size();pat++)
				{
	num=random_binomial(Pat[pat].AesFwwww,prob);Pat[pat].Fwwww+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesFwwww-=num;	
	num=random_binomial(Pat[pat].AesFwwwd,prob);Pat[pat].Fwwwd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesFwwwd-=num;	
	num=random_binomial(Pat[pat].AesFwwdd,prob);Pat[pat].Fwwdd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesFwwdd-=num;	
	num=random_binomial(Pat[pat].AesFwwwr,prob);Pat[pat].Fwwwr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesFwwwr-=num;	
	num=random_binomial(Pat[pat].AesFwwrr,prob);Pat[pat].Fwwrr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesFwwrr-=num;	
	num=random_binomial(Pat[pat].AesFwwdr,prob);Pat[pat].Fwwdr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesFwwdr-=num;	
	num=random_binomial(Pat[pat].AesfFwwww,prob);Pat[pat].fFwwww+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesfFwwww-=num;	
	num=random_binomial(Pat[pat].AesfFwwwd,prob);Pat[pat].fFwwwd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesfFwwwd-=num;	
	num=random_binomial(Pat[pat].AesfFwwdd,prob);Pat[pat].fFwwdd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesfFwwdd-=num;	
	num=random_binomial(Pat[pat].AesfFwwwr,prob);Pat[pat].fFwwwr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesfFwwwr-=num;	
	num=random_binomial(Pat[pat].AesfFwwrr,prob);Pat[pat].fFwwrr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesfFwwrr-=num;	
	num=random_binomial(Pat[pat].AesfFwwdr,prob);Pat[pat].fFwwdr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesfFwwdr-=num;	
		
	num=random_binomial(Pat[pat].AesmFwwww,prob);Pat[pat].mFwwww+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesmFwwww-=num;	
	num=random_binomial(Pat[pat].AesmFwwwd,prob);Pat[pat].mFwwwd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesmFwwwd-=num;	
	num=random_binomial(Pat[pat].AesmFwwdd,prob);Pat[pat].mFwwdd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesmFwwdd-=num;	
	num=random_binomial(Pat[pat].AesmFwwwr,prob);Pat[pat].mFwwwr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesmFwwwr-=num;	
	num=random_binomial(Pat[pat].AesmFwwrr,prob);Pat[pat].mFwwrr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesmFwwrr-=num;	
	num=random_binomial(Pat[pat].AesmFwwdr,prob);Pat[pat].mFwwdr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesmFwwdr-=num;	
	num=random_binomial(Pat[pat].AesbFwwww,prob);Pat[pat].bFwwww+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesbFwwww-=num;	
	num=random_binomial(Pat[pat].AesbFwwwd,prob);Pat[pat].bFwwwd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesbFwwwd-=num;	
	num=random_binomial(Pat[pat].AesbFwwdd,prob);Pat[pat].bFwwdd+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesbFwwdd-=num;	
	num=random_binomial(Pat[pat].AesbFwwwr,prob);Pat[pat].bFwwwr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesbFwwwr-=num;	
	num=random_binomial(Pat[pat].AesbFwwrr,prob);Pat[pat].bFwwrr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesbFwwrr-=num;	
	num=random_binomial(Pat[pat].AesbFwwdr,prob);Pat[pat].bFwwdr+=num;to.Fww+=num;to.FTot+=num;Pat[pat].AesbFwwdr-=num;	
	num=random_binomial(Pat[pat].AesfFwdww,prob);Pat[pat].fFwdww+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesfFwdww-=num;	
	num=random_binomial(Pat[pat].AesfFwdwd,prob);Pat[pat].fFwdwd+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesfFwdwd-=num;	
	num=random_binomial(Pat[pat].AesfFwddd,prob);Pat[pat].fFwddd+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesfFwddd-=num;	
	num=random_binomial(Pat[pat].AesfFwdwr,prob);Pat[pat].fFwdwr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesfFwdwr-=num;	
	num=random_binomial(Pat[pat].AesfFwdrr,prob);Pat[pat].fFwdrr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesfFwdrr-=num;	
	num=random_binomial(Pat[pat].AesfFwddr,prob);Pat[pat].fFwddr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesfFwddr-=num;	
	num=random_binomial(Pat[pat].AesmFwdww,prob);Pat[pat].mFwdww+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesmFwdww-=num;	
	num=random_binomial(Pat[pat].AesmFwdwd,prob);Pat[pat].mFwdwd+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesmFwdwd-=num;	
	num=random_binomial(Pat[pat].AesmFwddd,prob);Pat[pat].mFwddd+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesmFwddd-=num;	
	num=random_binomial(Pat[pat].AesmFwdwr,prob);Pat[pat].mFwdwr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesmFwdwr-=num;	
	num=random_binomial(Pat[pat].AesmFwdrr,prob);Pat[pat].mFwdrr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesmFwdrr-=num;	
	num=random_binomial(Pat[pat].AesmFwddr,prob);Pat[pat].mFwddr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesmFwddr-=num;	
					
	num=random_binomial(Pat[pat].AesbFwdww,prob);Pat[pat].bFwdww+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesbFwdww-=num;	
	num=random_binomial(Pat[pat].AesbFwdwd,prob);Pat[pat].bFwdwd+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesbFwdwd-=num;	
	num=random_binomial(Pat[pat].AesbFwddd,prob);Pat[pat].bFwddd+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesbFwddd-=num;	
	num=random_binomial(Pat[pat].AesbFwdwr,prob);Pat[pat].bFwdwr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesbFwdwr-=num;	
	num=random_binomial(Pat[pat].AesbFwdrr,prob);Pat[pat].bFwdrr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesbFwdrr-=num;	
	num=random_binomial(Pat[pat].AesbFwddr,prob);Pat[pat].bFwddr+=num;to.Fwd+=num;to.FTot+=num;Pat[pat].AesbFwddr-=num;	
	num=random_binomial(Pat[pat].AesFddww,prob);Pat[pat].Fddww+=num;to.Fdd+=num;to.FTot+=num;Pat[pat].AesFddww-=num;	
	num=random_binomial(Pat[pat].AesFddwd,prob);Pat[pat].Fddwd+=num;to.Fdd+=num;to.FTot+=num;Pat[pat].AesFddwd-=num;	
	num=random_binomial(Pat[pat].AesFdddd,prob);Pat[pat].Fdddd+=num;to.Fdd+=num;to.FTot+=num;Pat[pat].AesFdddd-=num;	
	num=random_binomial(Pat[pat].AesFddwr,prob);Pat[pat].Fddwr+=num;to.Fdd+=num;to.FTot+=num;Pat[pat].AesFddwr-=num;	
	num=random_binomial(Pat[pat].AesFddrr,prob);Pat[pat].Fddrr+=num;to.Fdd+=num;to.FTot+=num;Pat[pat].AesFddrr-=num;	
	num=random_binomial(Pat[pat].AesFdddr,prob);Pat[pat].Fdddr+=num;to.Fdd+=num;to.FTot+=num;Pat[pat].AesFdddr-=num;	

	num=random_binomial(Pat[pat].AesFwrww,prob);Pat[pat].Fwrww+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesFwrww-=num;	
	num=random_binomial(Pat[pat].AesFwrwd,prob);Pat[pat].Fwrwd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesFwrwd-=num;	
	num=random_binomial(Pat[pat].AesFwrdd,prob);Pat[pat].Fwrdd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesFwrdd-=num;	
	num=random_binomial(Pat[pat].AesFwrwr,prob);Pat[pat].Fwrwr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesFwrwr-=num;	
	num=random_binomial(Pat[pat].AesFwrrr,prob);Pat[pat].Fwrrr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesFwrrr-=num;	
	num=random_binomial(Pat[pat].AesFwrdr,prob);Pat[pat].Fwrdr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesFwrdr-=num;	

	num=random_binomial(Pat[pat].AesfFwrww,prob);Pat[pat].fFwrww+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesfFwrww-=num;	
	num=random_binomial(Pat[pat].AesfFwrwd,prob);Pat[pat].fFwrwd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesfFwrwd-=num;	
	num=random_binomial(Pat[pat].AesfFwrdd,prob);Pat[pat].fFwrdd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesfFwrdd-=num;	
	num=random_binomial(Pat[pat].AesfFwrwr,prob);Pat[pat].fFwrwr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesfFwrwr-=num;	
	num=random_binomial(Pat[pat].AesfFwrrr,prob);Pat[pat].fFwrrr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesfFwrrr-=num;	
	num=random_binomial(Pat[pat].AesfFwrdr,prob);Pat[pat].fFwrdr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesfFwrdr-=num;	
	num=random_binomial(Pat[pat].AesmFwrww,prob);Pat[pat].mFwrww+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesmFwrww-=num;	
	num=random_binomial(Pat[pat].AesmFwrwd,prob);Pat[pat].mFwrwd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesmFwrwd-=num;	
	num=random_binomial(Pat[pat].AesmFwrdd,prob);Pat[pat].mFwrdd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesmFwrdd-=num;	
	num=random_binomial(Pat[pat].AesmFwrwr,prob);Pat[pat].mFwrwr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesmFwrwr-=num;	
	num=random_binomial(Pat[pat].AesmFwrrr,prob);Pat[pat].mFwrrr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesmFwrrr-=num;	
	num=random_binomial(Pat[pat].AesmFwrdr,prob);Pat[pat].mFwrdr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesmFwrdr-=num;	
	num=random_binomial(Pat[pat].AesbFwrww,prob);Pat[pat].bFwrww+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesbFwrww-=num;	
	num=random_binomial(Pat[pat].AesbFwrwd,prob);Pat[pat].bFwrwd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesbFwrwd-=num;	
	num=random_binomial(Pat[pat].AesbFwrdd,prob);Pat[pat].bFwrdd+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesbFwrdd-=num;	
	num=random_binomial(Pat[pat].AesbFwrwr,prob);Pat[pat].bFwrwr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesbFwrwr-=num;	
	num=random_binomial(Pat[pat].AesbFwrrr,prob);Pat[pat].bFwrrr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesbFwrrr-=num;	
	num=random_binomial(Pat[pat].AesbFwrdr,prob);Pat[pat].bFwrdr+=num;to.Fwr+=num;to.FTot+=num;Pat[pat].AesbFwrdr-=num;	
	num=random_binomial(Pat[pat].AesFrrww,prob);Pat[pat].Frrww+=num;to.Frr+=num;to.FTot+=num;Pat[pat].AesFrrww-=num;	
	num=random_binomial(Pat[pat].AesFrrwd,prob);Pat[pat].Frrwd+=num;to.Frr+=num;to.FTot+=num;Pat[pat].AesFrrwd-=num;	
	num=random_binomial(Pat[pat].AesFrrdd,prob);Pat[pat].Frrdd+=num;to.Frr+=num;to.FTot+=num;Pat[pat].AesFrrdd-=num;	
	num=random_binomial(Pat[pat].AesFrrwr,prob);Pat[pat].Frrwr+=num;to.Frr+=num;to.FTot+=num;Pat[pat].AesFrrwr-=num;	
	num=random_binomial(Pat[pat].AesFrrrr,prob);Pat[pat].Frrrr+=num;to.Frr+=num;to.FTot+=num;Pat[pat].AesFrrrr-=num;	
	num=random_binomial(Pat[pat].AesFrrdr,prob);Pat[pat].Frrdr+=num;to.Frr+=num;to.FTot+=num;Pat[pat].AesFrrdr-=num;	

	num=random_binomial(Pat[pat].AesFdrww,prob);Pat[pat].Fdrww+=num;to.Fdr+=num;to.FTot+=num;Pat[pat].AesFdrww-=num;	
	num=random_binomial(Pat[pat].AesFdrwd,prob);Pat[pat].Fdrwd+=num;to.Fdr+=num;to.FTot+=num;Pat[pat].AesFdrwd-=num;	
	num=random_binomial(Pat[pat].AesFdrdd,prob);Pat[pat].Fdrdd+=num;to.Fdr+=num;to.FTot+=num;Pat[pat].AesFdrdd-=num;	
	num=random_binomial(Pat[pat].AesFdrwr,prob);Pat[pat].Fdrwr+=num;to.Fdr+=num;to.FTot+=num;Pat[pat].AesFdrwr-=num;	
	num=random_binomial(Pat[pat].AesFdrrr,prob);Pat[pat].Fdrrr+=num;to.Fdr+=num;to.FTot+=num;Pat[pat].AesFdrrr-=num;	
	num=random_binomial(Pat[pat].AesFdrdr,prob);Pat[pat].Fdrdr+=num;to.Fdr+=num;to.FTot+=num;Pat[pat].AesFdrdr-=num;	
				};
		return;};




	void LDM(char NorS){
		int ind,pat,newysq,newxsq,howmany,otherpat,surv,num;
		int* aww; 
		int mid=(nx+.5)/2;
		for(int pat=0;pat<Pat.size();pat++)
				{
	num=random_binomial(Pat[pat].Fwwww,pa.dL);Pat[pat].Fwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMFwwww=num;	
	num=random_binomial(Pat[pat].Fwwwd,pa.dL);Pat[pat].Fwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMFwwwd=num;	
	num=random_binomial(Pat[pat].Fwwdd,pa.dL);Pat[pat].Fwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMFwwdd=num;	
	num=random_binomial(Pat[pat].Fwwwr,pa.dL);Pat[pat].Fwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMFwwwr=num;	
	num=random_binomial(Pat[pat].Fwwrr,pa.dL);Pat[pat].Fwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMFwwrr=num;	
	num=random_binomial(Pat[pat].Fwwdr,pa.dL);Pat[pat].Fwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMFwwdr=num;	

	num=random_binomial(Pat[pat].fFwwww,pa.dL);Pat[pat].fFwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMfFwwww=num;	
	num=random_binomial(Pat[pat].fFwwwd,pa.dL);Pat[pat].fFwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMfFwwwd=num;	
	num=random_binomial(Pat[pat].fFwwdd,pa.dL);Pat[pat].fFwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMfFwwdd=num;	
	num=random_binomial(Pat[pat].fFwwwr,pa.dL);Pat[pat].fFwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMfFwwwr=num;	
	num=random_binomial(Pat[pat].fFwwrr,pa.dL);Pat[pat].fFwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMfFwwrr=num;	
	num=random_binomial(Pat[pat].fFwwdr,pa.dL);Pat[pat].fFwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMfFwwdr=num;	
	num=random_binomial(Pat[pat].mFwwww,pa.dL);Pat[pat].mFwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMmFwwww=num;	
	num=random_binomial(Pat[pat].mFwwwd,pa.dL);Pat[pat].mFwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMmFwwwd=num;	
	num=random_binomial(Pat[pat].mFwwdd,pa.dL);Pat[pat].mFwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMmFwwdd=num;	
	num=random_binomial(Pat[pat].mFwwwr,pa.dL);Pat[pat].mFwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMmFwwwr=num;	
	num=random_binomial(Pat[pat].mFwwrr,pa.dL);Pat[pat].mFwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMmFwwrr=num;	
	num=random_binomial(Pat[pat].mFwwdr,pa.dL);Pat[pat].mFwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMmFwwdr=num;	

	num=random_binomial(Pat[pat].bFwwww,pa.dL);Pat[pat].bFwwww-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMbFwwww=num;	
	num=random_binomial(Pat[pat].bFwwwd,pa.dL);Pat[pat].bFwwwd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMbFwwwd=num;	
	num=random_binomial(Pat[pat].bFwwdd,pa.dL);Pat[pat].bFwwdd-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMbFwwdd=num;	
	num=random_binomial(Pat[pat].bFwwwr,pa.dL);Pat[pat].bFwwwr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMbFwwwr=num;	
	num=random_binomial(Pat[pat].bFwwrr,pa.dL);Pat[pat].bFwwrr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMbFwwrr=num;	
	num=random_binomial(Pat[pat].bFwwdr,pa.dL);Pat[pat].bFwwdr-=num;to.Fww-=num;to.FTot-=num;Pat[pat].LDMbFwwdr=num;	

	num=random_binomial(Pat[pat].fFwdww,pa.dL);Pat[pat].fFwdww-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMfFwdww=num;	
	num=random_binomial(Pat[pat].fFwdwd,pa.dL);Pat[pat].fFwdwd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMfFwdwd=num;	
	num=random_binomial(Pat[pat].fFwddd,pa.dL);Pat[pat].fFwddd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMfFwddd=num;	
	num=random_binomial(Pat[pat].fFwdwr,pa.dL);Pat[pat].fFwdwr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMfFwdwr=num;	
	num=random_binomial(Pat[pat].fFwdrr,pa.dL);Pat[pat].fFwdrr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMfFwdrr=num;	
	num=random_binomial(Pat[pat].fFwddr,pa.dL);Pat[pat].fFwddr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMfFwddr=num;	

	num=random_binomial(Pat[pat].mFwdww,pa.dL);Pat[pat].mFwdww-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMmFwdww=num;	
	num=random_binomial(Pat[pat].mFwdwd,pa.dL);Pat[pat].mFwdwd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMmFwdwd=num;	
	num=random_binomial(Pat[pat].mFwddd,pa.dL);Pat[pat].mFwddd-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMmFwddd=num;	
	num=random_binomial(Pat[pat].mFwdwr,pa.dL);Pat[pat].mFwdwr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMmFwdwr=num;	
	num=random_binomial(Pat[pat].mFwdrr,pa.dL);Pat[pat].mFwdrr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMmFwdrr=num;	
	num=random_binomial(Pat[pat].mFwddr,pa.dL);Pat[pat].mFwddr-=num;to.Fwd-=num;to.FTot-=num;Pat[pat].LDMmFwddr=num;	

	num=random_binomial(Pat[pat].Fddww,pa.dL);Pat[pat].Fddww-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].LDMFddww=num;	
	num=random_binomial(Pat[pat].Fddwd,pa.dL);Pat[pat].Fddwd-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].LDMFddwd=num;	
	num=random_binomial(Pat[pat].Fdddd,pa.dL);Pat[pat].Fdddd-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].LDMFdddd=num;	
	num=random_binomial(Pat[pat].Fddwr,pa.dL);Pat[pat].Fddwr-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].LDMFddwr=num;	
	num=random_binomial(Pat[pat].Fddrr,pa.dL);Pat[pat].Fddrr-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].LDMFddrr=num;	
	num=random_binomial(Pat[pat].Fdddr,pa.dL);Pat[pat].Fdddr-=num;to.Fdd-=num;to.FTot-=num;Pat[pat].LDMFdddr=num;	

	num=random_binomial(Pat[pat].Fwrww,pa.dL);Pat[pat].Fwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMFwrww=num;	
	num=random_binomial(Pat[pat].Fwrwd,pa.dL);Pat[pat].Fwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMFwrwd=num;	
	num=random_binomial(Pat[pat].Fwrdd,pa.dL);Pat[pat].Fwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMFwrdd=num;	
	num=random_binomial(Pat[pat].Fwrwr,pa.dL);Pat[pat].Fwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMFwrwr=num;	
	num=random_binomial(Pat[pat].Fwrrr,pa.dL);Pat[pat].Fwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMFwrrr=num;	
	num=random_binomial(Pat[pat].Fwrdr,pa.dL);Pat[pat].Fwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMFwrdr=num;	


	num=random_binomial(Pat[pat].fFwrww,pa.dL);Pat[pat].fFwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMfFwrww=num;	
	num=random_binomial(Pat[pat].fFwrwd,pa.dL);Pat[pat].fFwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMfFwrwd=num;	
	num=random_binomial(Pat[pat].fFwrdd,pa.dL);Pat[pat].fFwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMfFwrdd=num;	
	num=random_binomial(Pat[pat].fFwrwr,pa.dL);Pat[pat].fFwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMfFwrwr=num;	
	num=random_binomial(Pat[pat].fFwrrr,pa.dL);Pat[pat].fFwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMfFwrrr=num;	
	num=random_binomial(Pat[pat].fFwrdr,pa.dL);Pat[pat].fFwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMfFwrdr=num;	

	num=random_binomial(Pat[pat].mFwrww,pa.dL);Pat[pat].mFwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMmFwrww=num;	
	num=random_binomial(Pat[pat].mFwrwd,pa.dL);Pat[pat].mFwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMmFwrwd=num;	
	num=random_binomial(Pat[pat].mFwrdd,pa.dL);Pat[pat].mFwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMmFwrdd=num;	
	num=random_binomial(Pat[pat].mFwrwr,pa.dL);Pat[pat].mFwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMmFwrwr=num;	
	num=random_binomial(Pat[pat].mFwrrr,pa.dL);Pat[pat].mFwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMmFwrrr=num;	
	num=random_binomial(Pat[pat].mFwrdr,pa.dL);Pat[pat].mFwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMmFwrdr=num;	

	num=random_binomial(Pat[pat].bFwrww,pa.dL);Pat[pat].bFwrww-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMbFwrww=num;	
	num=random_binomial(Pat[pat].bFwrwd,pa.dL);Pat[pat].bFwrwd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMbFwrwd=num;	
	num=random_binomial(Pat[pat].bFwrdd,pa.dL);Pat[pat].bFwrdd-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMbFwrdd=num;	
	num=random_binomial(Pat[pat].bFwrwr,pa.dL);Pat[pat].bFwrwr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMbFwrwr=num;	
	num=random_binomial(Pat[pat].bFwrrr,pa.dL);Pat[pat].bFwrrr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMbFwrrr=num;	
	num=random_binomial(Pat[pat].bFwrdr,pa.dL);Pat[pat].bFwrdr-=num;to.Fwr-=num;to.FTot-=num;Pat[pat].LDMbFwrdr=num;	

	num=random_binomial(Pat[pat].Frrww,pa.dL);Pat[pat].Frrww-=num;to.Frr-=num;to.FTot-=num;Pat[pat].LDMFrrww=num;	
	num=random_binomial(Pat[pat].Frrwd,pa.dL);Pat[pat].Frrwd-=num;to.Frr-=num;to.FTot-=num;Pat[pat].LDMFrrwd=num;	
	num=random_binomial(Pat[pat].Frrdd,pa.dL);Pat[pat].Frrdd-=num;to.Frr-=num;to.FTot-=num;Pat[pat].LDMFrrdd=num;	
	num=random_binomial(Pat[pat].Frrwr,pa.dL);Pat[pat].Frrwr-=num;to.Frr-=num;to.FTot-=num;Pat[pat].LDMFrrwr=num;	
	num=random_binomial(Pat[pat].Frrrr,pa.dL);Pat[pat].Frrrr-=num;to.Frr-=num;to.FTot-=num;Pat[pat].LDMFrrrr=num;	
	num=random_binomial(Pat[pat].Frrdr,pa.dL);Pat[pat].Frrdr-=num;to.Frr-=num;to.FTot-=num;Pat[pat].LDMFrrdr=num;	

	num=random_binomial(Pat[pat].Fdrww,pa.dL);Pat[pat].Fdrww-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].LDMFdrww=num;	
	num=random_binomial(Pat[pat].Fdrwd,pa.dL);Pat[pat].Fdrwd-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].LDMFdrwd=num;	
	num=random_binomial(Pat[pat].Fdrdd,pa.dL);Pat[pat].Fdrdd-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].LDMFdrdd=num;	
	num=random_binomial(Pat[pat].Fdrwr,pa.dL);Pat[pat].Fdrwr-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].LDMFdrwr=num;	
	num=random_binomial(Pat[pat].Fdrrr,pa.dL);Pat[pat].Fdrrr-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].LDMFdrrr=num;	
	num=random_binomial(Pat[pat].Fdrdr,pa.dL);Pat[pat].Fdrdr-=num;to.Fdr-=num;to.FTot-=num;Pat[pat].LDMFdrdr=num;	
				};
		if(NorS=='N')
		{
		for(int xx=0;xx<nx-1;xx++)
		{
			for(int yy=0;yy<ny-2;yy++)
			{
				for(ind=0;ind<SetsPerCell[xx][yy].size();ind++)
				{
				pat=SetsPerCell[xx][yy][ind];
				surv=random_binomial(Pat[pat].LDMFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};








				surv=random_binomial(Pat[pat].LDMfFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};


				surv=random_binomial(Pat[pat].LDMmFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};





				surv=random_binomial(Pat[pat].LDMbFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};





				surv=random_binomial(Pat[pat].LDMfFwdww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdww+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwdwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdwd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwddd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwdwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdwr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwdrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdrr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwddr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};




				surv=random_binomial(Pat[pat].LDMmFwdww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdww+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwdwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdwd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwddd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwdwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdwr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwdrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdrr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwddr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};


				surv=random_binomial(Pat[pat].LDMbFwdww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdww+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwdwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdwd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwddd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwdwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdwr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwdrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdrr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwddr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};



				surv=random_binomial(Pat[pat].LDMFddww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddww+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFddwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddwd+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdddd+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFddwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddwr+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFddrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddrr+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdddr+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};







				surv=random_binomial(Pat[pat].LDMfFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};






				surv=random_binomial(Pat[pat].LDMmFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};



				surv=random_binomial(Pat[pat].LDMbFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};





				surv=random_binomial(Pat[pat].LDMFrrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrww+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrwd+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrdd+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrwr+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrrr+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrdr+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrww+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrwd+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrdd+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrwr+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrrr+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy,yy+int((nx+mid-yy-xx)/2));
					newxsq=xx+(newysq-yy);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrdr+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};

				};
				};
				};
				};
		if(NorS=='S')
		{
		for(int xx=nx-1;xx>=0;xx--)
		{
			for(int yy=ny-1;yy>0;yy--)
			{
				for(ind=0;ind<SetsPerCell[xx][yy].size();ind++)
				{
				pat=SetsPerCell[xx][yy][ind];
				surv=random_binomial(Pat[pat].LDMFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};






				surv=random_binomial(Pat[pat].LDMfFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};


				surv=random_binomial(Pat[pat].LDMmFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};



				surv=random_binomial(Pat[pat].LDMbFwwww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwww+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwwd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwdd+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwwr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwrr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwwdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwwdr+=*(aww+otherpat);
						to.Fww+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};




				surv=random_binomial(Pat[pat].LDMfFwdww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdww+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwdwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdwd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwddd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwdwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdwr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwdrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwdrr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwddr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};


				surv=random_binomial(Pat[pat].LDMmFwdww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdww+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwdwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdwd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwddd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwdwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdwr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwdrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwdrr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwddr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};



				surv=random_binomial(Pat[pat].LDMbFwdww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdww+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwdwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdwd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwddd+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwdwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdwr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwdrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwdrr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwddr+=*(aww+otherpat);
						to.Fwd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};



				surv=random_binomial(Pat[pat].LDMFddww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddww+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFddwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddwd+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdddd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdddd+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFddwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddwr+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFddrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fddrr+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdddr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdddr+=*(aww+otherpat);
						to.Fdd+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};

				surv=random_binomial(Pat[pat].LDMFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};




				surv=random_binomial(Pat[pat].LDMfFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMfFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].fFwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};





				surv=random_binomial(Pat[pat].LDMmFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMmFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].mFwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};





				surv=random_binomial(Pat[pat].LDMbFwrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrww+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrwd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrdd+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrwr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrrr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMbFwrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].bFwrdr+=*(aww+otherpat);
						to.Fwr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};









				surv=random_binomial(Pat[pat].LDMFrrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrww+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrwd+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrdd+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrwr+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrrr+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFrrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Frrdr+=*(aww+otherpat);
						to.Frr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};

				surv=random_binomial(Pat[pat].LDMFdrww,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrww+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrwd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrwd+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrdd,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrdd+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrwr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrwr+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrrr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrrr+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				surv=random_binomial(Pat[pat].LDMFdrdr,1-pa.muLD);
				if(surv>0)
				{
					newysq=rg.IRandom(yy-int((yy+xx-mid)/2.0),yy);
					newxsq=xx-(yy-newysq);
					howmany=SetsPerCell[newxsq][newysq].size();
					if(howmany>0)
					{
					aww=random_multinomEqualProb(surv,howmany);
					for(otherpat=0;otherpat<SetsPerCell[newxsq][newysq].size();otherpat++) 
						{
						Pat[SetsPerCell[newxsq][newysq][otherpat]].Fdrdr+=*(aww+otherpat);
						to.Fdr+=*(aww+otherpat); to.FTot+=*(aww+otherpat);
						};
					delete[] aww;
					};
				};
				};
				};
		};
		};
	return;};



	void LayEggs (void){


		for(int pat=0;pat<Pat.size();pat++)
				{

//fraction of 1:ww,2:fww,3:mww,4:fwd,5:mwd,6:dd,7:wr,8:fwr,9:mwr,10:rr,11:dr,12:bww,13:bwd,14:bwr
	Pat[pat].Jww[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[0]+ Pat[pat].Fwwwd*pa.fwwwd[0]+ Pat[pat].Fwwdd*pa.fwwdd[0]+ Pat[pat].Fwwwr*pa.fwwwr[0]+ Pat[pat].Fwwrr*pa.fwwrr[0]+ Pat[pat].Fwwdr*pa.fwwdr[0]+ Pat[pat].fFwwww*pa.ffwwww[0]+ Pat[pat].fFwwwd*pa.ffwwwd[0]+ Pat[pat].fFwwdd*pa.ffwwdd[0]+ Pat[pat].fFwwwr*pa.ffwwwr[0]+ Pat[pat].fFwwrr*pa.ffwwrr[0]+ Pat[pat].fFwwdr*pa.ffwwdr[0]+ Pat[pat].mFwwww*pa.mfwwww[0]+ Pat[pat].mFwwwd*pa.mfwwwd[0]+ Pat[pat].mFwwdd*pa.mfwwdd[0]+ Pat[pat].mFwwwr*pa.mfwwwr[0]+ Pat[pat].mFwwrr*pa.mfwwrr[0]+ Pat[pat].mFwwdr*pa.mfwwdr[0]+ Pat[pat].bFwwww*pa.bfwwww[0]+ Pat[pat].bFwwwd*pa.bfwwwd[0]+ Pat[pat].bFwwdd*pa.bfwwdd[0]+ Pat[pat].bFwwwr*pa.bfwwwr[0]+ Pat[pat].bFwwrr*pa.bfwwrr[0]+ Pat[pat].bFwwdr*pa.bfwwdr[0] +Pat[pat].fFwdww*pa.ffwdww[0]+ Pat[pat].fFwdwd*pa.ffwdwd[0]+ Pat[pat].fFwddd*pa.ffwddd[0]+ Pat[pat].fFwdwr*pa.ffwdwr[0]+ Pat[pat].fFwdrr*pa.ffwdrr[0]+ Pat[pat].fFwddr*pa.ffwddr[0] +Pat[pat].mFwdww*pa.mfwdww[0]+ Pat[pat].mFwdwd*pa.mfwdwd[0]+ Pat[pat].mFwddd*pa.mfwddd[0]+ Pat[pat].mFwdwr*pa.mfwdwr[0]+ Pat[pat].mFwdrr*pa.mfwdrr[0]+ Pat[pat].mFwddr*pa.mfwddr[0] +Pat[pat].bFwdww*pa.bfwdww[0]+ Pat[pat].bFwdwd*pa.bfwdwd[0]+ Pat[pat].bFwddd*pa.bfwddd[0]+ Pat[pat].bFwdwr*pa.bfwdwr[0]+ Pat[pat].bFwdrr*pa.bfwdrr[0]+ Pat[pat].bFwddr*pa.bfwddr[0] + Pat[pat].Fddww*pa.fddww[0]+ Pat[pat].Fddwd*pa.fddwd[0]+ Pat[pat].Fdddd*pa.fdddd[0]+ Pat[pat].Fddwr*pa.fddwr[0]+ Pat[pat].Fddrr*pa.fddrr[0]+ Pat[pat].Fdddr*pa.fdddr[0] + Pat[pat].Fwrww*pa.fwrww[0]+ Pat[pat].Fwrwd*pa.fwrwd[0]+ Pat[pat].Fwrdd*pa.fwrdd[0]+ Pat[pat].Fwrwr*pa.fwrwr[0]+ Pat[pat].Fwrrr*pa.fwrrr[0]+ Pat[pat].Fwrdr*pa.fwrdr[0] + Pat[pat].fFwrww*pa.ffwrww[0]+ Pat[pat].fFwrwd*pa.ffwrwd[0]+ Pat[pat].fFwrdd*pa.ffwrdd[0]+ Pat[pat].fFwrwr*pa.ffwrwr[0]+ Pat[pat].fFwrrr*pa.ffwrrr[0]+ Pat[pat].fFwrdr*pa.ffwrdr[0] + Pat[pat].mFwrww*pa.mfwrww[0]+ Pat[pat].mFwrwd*pa.mfwrwd[0]+ Pat[pat].mFwrdd*pa.mfwrdd[0]+ Pat[pat].mFwrwr*pa.mfwrwr[0]+ Pat[pat].mFwrrr*pa.mfwrrr[0]+ Pat[pat].mFwrdr*pa.mfwrdr[0] + Pat[pat].bFwrww*pa.bfwrww[0]+ Pat[pat].bFwrwd*pa.bfwrwd[0]+ Pat[pat].bFwrdd*pa.bfwrdd[0]+ Pat[pat].bFwrwr*pa.bfwrwr[0]+ Pat[pat].bFwrrr*pa.bfwrrr[0]+ Pat[pat].bFwrdr*pa.bfwrdr[0] + Pat[pat].Frrww*pa.frrww[0]+ Pat[pat].Frrwd*pa.frrwd[0]+ Pat[pat].Frrdd*pa.frrdd[0]+ Pat[pat].Frrwr*pa.frrwr[0]+ Pat[pat].Frrrr*pa.frrrr[0]+ Pat[pat].Frrdr*pa.frrdr[0] + Pat[pat].Fdrww*pa.fdrww[0]+ Pat[pat].Fdrwd*pa.fdrwd[0]+ Pat[pat].Fdrdd*pa.fdrdd[0]+ Pat[pat].Fdrwr*pa.fdrwr[0]+ Pat[pat].Fdrrr*pa.fdrrr[0]+ Pat[pat].Fdrdr*pa.fdrdr[0]));
	Pat[pat].fJww[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[1]+ Pat[pat].Fwwwd*pa.fwwwd[1]+ Pat[pat].Fwwdd*pa.fwwdd[1]+ Pat[pat].Fwwwr*pa.fwwwr[1]+ Pat[pat].Fwwrr*pa.fwwrr[1]+ Pat[pat].Fwwdr*pa.fwwdr[1]+ Pat[pat].fFwwww*pa.ffwwww[1]+ Pat[pat].fFwwwd*pa.ffwwwd[1]+ Pat[pat].fFwwdd*pa.ffwwdd[1]+ Pat[pat].fFwwwr*pa.ffwwwr[1]+ Pat[pat].fFwwrr*pa.ffwwrr[1]+ Pat[pat].fFwwdr*pa.ffwwdr[1]+ Pat[pat].mFwwww*pa.mfwwww[1]+ Pat[pat].mFwwwd*pa.mfwwwd[1]+ Pat[pat].mFwwdd*pa.mfwwdd[1]+ Pat[pat].mFwwwr*pa.mfwwwr[1]+ Pat[pat].mFwwrr*pa.mfwwrr[1]+ Pat[pat].mFwwdr*pa.mfwwdr[1]+ Pat[pat].bFwwww*pa.bfwwww[1]+ Pat[pat].bFwwwd*pa.bfwwwd[1]+ Pat[pat].bFwwdd*pa.bfwwdd[1]+ Pat[pat].bFwwwr*pa.bfwwwr[1]+ Pat[pat].bFwwrr*pa.bfwwrr[1]+ Pat[pat].bFwwdr*pa.bfwwdr[1] +Pat[pat].fFwdww*pa.ffwdww[1]+ Pat[pat].fFwdwd*pa.ffwdwd[1]+ Pat[pat].fFwddd*pa.ffwddd[1]+ Pat[pat].fFwdwr*pa.ffwdwr[1]+ Pat[pat].fFwdrr*pa.ffwdrr[1]+ Pat[pat].fFwddr*pa.ffwddr[1] +Pat[pat].mFwdww*pa.mfwdww[1]+ Pat[pat].mFwdwd*pa.mfwdwd[1]+ Pat[pat].mFwddd*pa.mfwddd[1]+ Pat[pat].mFwdwr*pa.mfwdwr[1]+ Pat[pat].mFwdrr*pa.mfwdrr[1]+ Pat[pat].mFwddr*pa.mfwddr[1] +Pat[pat].bFwdww*pa.bfwdww[1]+ Pat[pat].bFwdwd*pa.bfwdwd[1]+ Pat[pat].bFwddd*pa.bfwddd[1]+ Pat[pat].bFwdwr*pa.bfwdwr[1]+ Pat[pat].bFwdrr*pa.bfwdrr[1]+ Pat[pat].bFwddr*pa.bfwddr[1] + Pat[pat].Fddww*pa.fddww[1]+ Pat[pat].Fddwd*pa.fddwd[1]+ Pat[pat].Fdddd*pa.fdddd[1]+ Pat[pat].Fddwr*pa.fddwr[1]+ Pat[pat].Fddrr*pa.fddrr[1]+ Pat[pat].Fdddr*pa.fdddr[1] + Pat[pat].Fwrww*pa.fwrww[1]+ Pat[pat].Fwrwd*pa.fwrwd[1]+ Pat[pat].Fwrdd*pa.fwrdd[1]+ Pat[pat].Fwrwr*pa.fwrwr[1]+ Pat[pat].Fwrrr*pa.fwrrr[1]+ Pat[pat].Fwrdr*pa.fwrdr[1] + Pat[pat].fFwrww*pa.ffwrww[1]+ Pat[pat].fFwrwd*pa.ffwrwd[1]+ Pat[pat].fFwrdd*pa.ffwrdd[1]+ Pat[pat].fFwrwr*pa.ffwrwr[1]+ Pat[pat].fFwrrr*pa.ffwrrr[1]+ Pat[pat].fFwrdr*pa.ffwrdr[1] + Pat[pat].mFwrww*pa.mfwrww[1]+ Pat[pat].mFwrwd*pa.mfwrwd[1]+ Pat[pat].mFwrdd*pa.mfwrdd[1]+ Pat[pat].mFwrwr*pa.mfwrwr[1]+ Pat[pat].mFwrrr*pa.mfwrrr[1]+ Pat[pat].mFwrdr*pa.mfwrdr[1] + Pat[pat].bFwrww*pa.bfwrww[1]+ Pat[pat].bFwrwd*pa.bfwrwd[1]+ Pat[pat].bFwrdd*pa.bfwrdd[1]+ Pat[pat].bFwrwr*pa.bfwrwr[1]+ Pat[pat].bFwrrr*pa.bfwrrr[1]+ Pat[pat].bFwrdr*pa.bfwrdr[1] + Pat[pat].Frrww*pa.frrww[1]+ Pat[pat].Frrwd*pa.frrwd[1]+ Pat[pat].Frrdd*pa.frrdd[1]+ Pat[pat].Frrwr*pa.frrwr[1]+ Pat[pat].Frrrr*pa.frrrr[1]+ Pat[pat].Frrdr*pa.frrdr[1] + Pat[pat].Fdrww*pa.fdrww[1]+ Pat[pat].Fdrwd*pa.fdrwd[1]+ Pat[pat].Fdrdd*pa.fdrdd[1]+ Pat[pat].Fdrwr*pa.fdrwr[1]+ Pat[pat].Fdrrr*pa.fdrrr[1]+ Pat[pat].Fdrdr*pa.fdrdr[1]));
	Pat[pat].mJww[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[2]+ Pat[pat].Fwwwd*pa.fwwwd[2]+ Pat[pat].Fwwdd*pa.fwwdd[2]+ Pat[pat].Fwwwr*pa.fwwwr[2]+ Pat[pat].Fwwrr*pa.fwwrr[2]+ Pat[pat].Fwwdr*pa.fwwdr[2]+ Pat[pat].fFwwww*pa.ffwwww[2]+ Pat[pat].fFwwwd*pa.ffwwwd[2]+ Pat[pat].fFwwdd*pa.ffwwdd[2]+ Pat[pat].fFwwwr*pa.ffwwwr[2]+ Pat[pat].fFwwrr*pa.ffwwrr[2]+ Pat[pat].fFwwdr*pa.ffwwdr[2]+ Pat[pat].mFwwww*pa.mfwwww[2]+ Pat[pat].mFwwwd*pa.mfwwwd[2]+ Pat[pat].mFwwdd*pa.mfwwdd[2]+ Pat[pat].mFwwwr*pa.mfwwwr[2]+ Pat[pat].mFwwrr*pa.mfwwrr[2]+ Pat[pat].mFwwdr*pa.mfwwdr[2]+ Pat[pat].bFwwww*pa.bfwwww[2]+ Pat[pat].bFwwwd*pa.bfwwwd[2]+ Pat[pat].bFwwdd*pa.bfwwdd[2]+ Pat[pat].bFwwwr*pa.bfwwwr[2]+ Pat[pat].bFwwrr*pa.bfwwrr[2]+ Pat[pat].bFwwdr*pa.bfwwdr[2] +Pat[pat].fFwdww*pa.ffwdww[2]+ Pat[pat].fFwdwd*pa.ffwdwd[2]+ Pat[pat].fFwddd*pa.ffwddd[2]+ Pat[pat].fFwdwr*pa.ffwdwr[2]+ Pat[pat].fFwdrr*pa.ffwdrr[2]+ Pat[pat].fFwddr*pa.ffwddr[2] +Pat[pat].mFwdww*pa.mfwdww[2]+ Pat[pat].mFwdwd*pa.mfwdwd[2]+ Pat[pat].mFwddd*pa.mfwddd[2]+ Pat[pat].mFwdwr*pa.mfwdwr[2]+ Pat[pat].mFwdrr*pa.mfwdrr[2]+ Pat[pat].mFwddr*pa.mfwddr[2] +Pat[pat].bFwdww*pa.bfwdww[2]+ Pat[pat].bFwdwd*pa.bfwdwd[2]+ Pat[pat].bFwddd*pa.bfwddd[2]+ Pat[pat].bFwdwr*pa.bfwdwr[2]+ Pat[pat].bFwdrr*pa.bfwdrr[2]+ Pat[pat].bFwddr*pa.bfwddr[2] + Pat[pat].Fddww*pa.fddww[2]+ Pat[pat].Fddwd*pa.fddwd[2]+ Pat[pat].Fdddd*pa.fdddd[2]+ Pat[pat].Fddwr*pa.fddwr[2]+ Pat[pat].Fddrr*pa.fddrr[2]+ Pat[pat].Fdddr*pa.fdddr[2] + Pat[pat].Fwrww*pa.fwrww[2]+ Pat[pat].Fwrwd*pa.fwrwd[2]+ Pat[pat].Fwrdd*pa.fwrdd[2]+ Pat[pat].Fwrwr*pa.fwrwr[2]+ Pat[pat].Fwrrr*pa.fwrrr[2]+ Pat[pat].Fwrdr*pa.fwrdr[2] + Pat[pat].fFwrww*pa.ffwrww[2]+ Pat[pat].fFwrwd*pa.ffwrwd[2]+ Pat[pat].fFwrdd*pa.ffwrdd[2]+ Pat[pat].fFwrwr*pa.ffwrwr[2]+ Pat[pat].fFwrrr*pa.ffwrrr[2]+ Pat[pat].fFwrdr*pa.ffwrdr[2] + Pat[pat].mFwrww*pa.mfwrww[2]+ Pat[pat].mFwrwd*pa.mfwrwd[2]+ Pat[pat].mFwrdd*pa.mfwrdd[2]+ Pat[pat].mFwrwr*pa.mfwrwr[2]+ Pat[pat].mFwrrr*pa.mfwrrr[2]+ Pat[pat].mFwrdr*pa.mfwrdr[2] + Pat[pat].bFwrww*pa.bfwrww[2]+ Pat[pat].bFwrwd*pa.bfwrwd[2]+ Pat[pat].bFwrdd*pa.bfwrdd[2]+ Pat[pat].bFwrwr*pa.bfwrwr[2]+ Pat[pat].bFwrrr*pa.bfwrrr[2]+ Pat[pat].bFwrdr*pa.bfwrdr[2] + Pat[pat].Frrww*pa.frrww[2]+ Pat[pat].Frrwd*pa.frrwd[2]+ Pat[pat].Frrdd*pa.frrdd[2]+ Pat[pat].Frrwr*pa.frrwr[2]+ Pat[pat].Frrrr*pa.frrrr[2]+ Pat[pat].Frrdr*pa.frrdr[2] + Pat[pat].Fdrww*pa.fdrww[2]+ Pat[pat].Fdrwd*pa.fdrwd[2]+ Pat[pat].Fdrdd*pa.fdrdd[2]+ Pat[pat].Fdrwr*pa.fdrwr[2]+ Pat[pat].Fdrrr*pa.fdrrr[2]+ Pat[pat].Fdrdr*pa.fdrdr[2]));
	Pat[pat].fJwd[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[3]+ Pat[pat].Fwwwd*pa.fwwwd[3]+ Pat[pat].Fwwdd*pa.fwwdd[3]+ Pat[pat].Fwwwr*pa.fwwwr[3]+ Pat[pat].Fwwrr*pa.fwwrr[3]+ Pat[pat].Fwwdr*pa.fwwdr[3]+ Pat[pat].fFwwww*pa.ffwwww[3]+ Pat[pat].fFwwwd*pa.ffwwwd[3]+ Pat[pat].fFwwdd*pa.ffwwdd[3]+ Pat[pat].fFwwwr*pa.ffwwwr[3]+ Pat[pat].fFwwrr*pa.ffwwrr[3]+ Pat[pat].fFwwdr*pa.ffwwdr[3]+ Pat[pat].mFwwww*pa.mfwwww[3]+ Pat[pat].mFwwwd*pa.mfwwwd[3]+ Pat[pat].mFwwdd*pa.mfwwdd[3]+ Pat[pat].mFwwwr*pa.mfwwwr[3]+ Pat[pat].mFwwrr*pa.mfwwrr[3]+ Pat[pat].mFwwdr*pa.mfwwdr[3]+ Pat[pat].bFwwww*pa.bfwwww[3]+ Pat[pat].bFwwwd*pa.bfwwwd[3]+ Pat[pat].bFwwdd*pa.bfwwdd[3]+ Pat[pat].bFwwwr*pa.bfwwwr[3]+ Pat[pat].bFwwrr*pa.bfwwrr[3]+ Pat[pat].bFwwdr*pa.bfwwdr[3] +Pat[pat].fFwdww*pa.ffwdww[3]+ Pat[pat].fFwdwd*pa.ffwdwd[3]+ Pat[pat].fFwddd*pa.ffwddd[3]+ Pat[pat].fFwdwr*pa.ffwdwr[3]+ Pat[pat].fFwdrr*pa.ffwdrr[3]+ Pat[pat].fFwddr*pa.ffwddr[3] +Pat[pat].mFwdww*pa.mfwdww[3]+ Pat[pat].mFwdwd*pa.mfwdwd[3]+ Pat[pat].mFwddd*pa.mfwddd[3]+ Pat[pat].mFwdwr*pa.mfwdwr[3]+ Pat[pat].mFwdrr*pa.mfwdrr[3]+ Pat[pat].mFwddr*pa.mfwddr[3] +Pat[pat].bFwdww*pa.bfwdww[3]+ Pat[pat].bFwdwd*pa.bfwdwd[3]+ Pat[pat].bFwddd*pa.bfwddd[3]+ Pat[pat].bFwdwr*pa.bfwdwr[3]+ Pat[pat].bFwdrr*pa.bfwdrr[3]+ Pat[pat].bFwddr*pa.bfwddr[3] + Pat[pat].Fddww*pa.fddww[3]+ Pat[pat].Fddwd*pa.fddwd[3]+ Pat[pat].Fdddd*pa.fdddd[3]+ Pat[pat].Fddwr*pa.fddwr[3]+ Pat[pat].Fddrr*pa.fddrr[3]+ Pat[pat].Fdddr*pa.fdddr[3] + Pat[pat].Fwrww*pa.fwrww[3]+ Pat[pat].Fwrwd*pa.fwrwd[3]+ Pat[pat].Fwrdd*pa.fwrdd[3]+ Pat[pat].Fwrwr*pa.fwrwr[3]+ Pat[pat].Fwrrr*pa.fwrrr[3]+ Pat[pat].Fwrdr*pa.fwrdr[3] + Pat[pat].fFwrww*pa.ffwrww[3]+ Pat[pat].fFwrwd*pa.ffwrwd[3]+ Pat[pat].fFwrdd*pa.ffwrdd[3]+ Pat[pat].fFwrwr*pa.ffwrwr[3]+ Pat[pat].fFwrrr*pa.ffwrrr[3]+ Pat[pat].fFwrdr*pa.ffwrdr[3] + Pat[pat].mFwrww*pa.mfwrww[3]+ Pat[pat].mFwrwd*pa.mfwrwd[3]+ Pat[pat].mFwrdd*pa.mfwrdd[3]+ Pat[pat].mFwrwr*pa.mfwrwr[3]+ Pat[pat].mFwrrr*pa.mfwrrr[3]+ Pat[pat].mFwrdr*pa.mfwrdr[3] + Pat[pat].bFwrww*pa.bfwrww[3]+ Pat[pat].bFwrwd*pa.bfwrwd[3]+ Pat[pat].bFwrdd*pa.bfwrdd[3]+ Pat[pat].bFwrwr*pa.bfwrwr[3]+ Pat[pat].bFwrrr*pa.bfwrrr[3]+ Pat[pat].bFwrdr*pa.bfwrdr[3] + Pat[pat].Frrww*pa.frrww[3]+ Pat[pat].Frrwd*pa.frrwd[3]+ Pat[pat].Frrdd*pa.frrdd[3]+ Pat[pat].Frrwr*pa.frrwr[3]+ Pat[pat].Frrrr*pa.frrrr[3]+ Pat[pat].Frrdr*pa.frrdr[3] + Pat[pat].Fdrww*pa.fdrww[3]+ Pat[pat].Fdrwd*pa.fdrwd[3]+ Pat[pat].Fdrdd*pa.fdrdd[3]+ Pat[pat].Fdrwr*pa.fdrwr[3]+ Pat[pat].Fdrrr*pa.fdrrr[3]+ Pat[pat].Fdrdr*pa.fdrdr[3]));
	Pat[pat].mJwd[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[4]+ Pat[pat].Fwwwd*pa.fwwwd[4]+ Pat[pat].Fwwdd*pa.fwwdd[4]+ Pat[pat].Fwwwr*pa.fwwwr[4]+ Pat[pat].Fwwrr*pa.fwwrr[4]+ Pat[pat].Fwwdr*pa.fwwdr[4]+ Pat[pat].fFwwww*pa.ffwwww[4]+ Pat[pat].fFwwwd*pa.ffwwwd[4]+ Pat[pat].fFwwdd*pa.ffwwdd[4]+ Pat[pat].fFwwwr*pa.ffwwwr[4]+ Pat[pat].fFwwrr*pa.ffwwrr[4]+ Pat[pat].fFwwdr*pa.ffwwdr[4]+ Pat[pat].mFwwww*pa.mfwwww[4]+ Pat[pat].mFwwwd*pa.mfwwwd[4]+ Pat[pat].mFwwdd*pa.mfwwdd[4]+ Pat[pat].mFwwwr*pa.mfwwwr[4]+ Pat[pat].mFwwrr*pa.mfwwrr[4]+ Pat[pat].mFwwdr*pa.mfwwdr[4]+ Pat[pat].bFwwww*pa.bfwwww[4]+ Pat[pat].bFwwwd*pa.bfwwwd[4]+ Pat[pat].bFwwdd*pa.bfwwdd[4]+ Pat[pat].bFwwwr*pa.bfwwwr[4]+ Pat[pat].bFwwrr*pa.bfwwrr[4]+ Pat[pat].bFwwdr*pa.bfwwdr[4] +Pat[pat].fFwdww*pa.ffwdww[4]+ Pat[pat].fFwdwd*pa.ffwdwd[4]+ Pat[pat].fFwddd*pa.ffwddd[4]+ Pat[pat].fFwdwr*pa.ffwdwr[4]+ Pat[pat].fFwdrr*pa.ffwdrr[4]+ Pat[pat].fFwddr*pa.ffwddr[4] +Pat[pat].mFwdww*pa.mfwdww[4]+ Pat[pat].mFwdwd*pa.mfwdwd[4]+ Pat[pat].mFwddd*pa.mfwddd[4]+ Pat[pat].mFwdwr*pa.mfwdwr[4]+ Pat[pat].mFwdrr*pa.mfwdrr[4]+ Pat[pat].mFwddr*pa.mfwddr[4] +Pat[pat].bFwdww*pa.bfwdww[4]+ Pat[pat].bFwdwd*pa.bfwdwd[4]+ Pat[pat].bFwddd*pa.bfwddd[4]+ Pat[pat].bFwdwr*pa.bfwdwr[4]+ Pat[pat].bFwdrr*pa.bfwdrr[4]+ Pat[pat].bFwddr*pa.bfwddr[4] + Pat[pat].Fddww*pa.fddww[4]+ Pat[pat].Fddwd*pa.fddwd[4]+ Pat[pat].Fdddd*pa.fdddd[4]+ Pat[pat].Fddwr*pa.fddwr[4]+ Pat[pat].Fddrr*pa.fddrr[4]+ Pat[pat].Fdddr*pa.fdddr[4] + Pat[pat].Fwrww*pa.fwrww[4]+ Pat[pat].Fwrwd*pa.fwrwd[4]+ Pat[pat].Fwrdd*pa.fwrdd[4]+ Pat[pat].Fwrwr*pa.fwrwr[4]+ Pat[pat].Fwrrr*pa.fwrrr[4]+ Pat[pat].Fwrdr*pa.fwrdr[4] + Pat[pat].fFwrww*pa.ffwrww[4]+ Pat[pat].fFwrwd*pa.ffwrwd[4]+ Pat[pat].fFwrdd*pa.ffwrdd[4]+ Pat[pat].fFwrwr*pa.ffwrwr[4]+ Pat[pat].fFwrrr*pa.ffwrrr[4]+ Pat[pat].fFwrdr*pa.ffwrdr[4] + Pat[pat].mFwrww*pa.mfwrww[4]+ Pat[pat].mFwrwd*pa.mfwrwd[4]+ Pat[pat].mFwrdd*pa.mfwrdd[4]+ Pat[pat].mFwrwr*pa.mfwrwr[4]+ Pat[pat].mFwrrr*pa.mfwrrr[4]+ Pat[pat].mFwrdr*pa.mfwrdr[4] + Pat[pat].bFwrww*pa.bfwrww[4]+ Pat[pat].bFwrwd*pa.bfwrwd[4]+ Pat[pat].bFwrdd*pa.bfwrdd[4]+ Pat[pat].bFwrwr*pa.bfwrwr[4]+ Pat[pat].bFwrrr*pa.bfwrrr[4]+ Pat[pat].bFwrdr*pa.bfwrdr[4] + Pat[pat].Frrww*pa.frrww[4]+ Pat[pat].Frrwd*pa.frrwd[4]+ Pat[pat].Frrdd*pa.frrdd[4]+ Pat[pat].Frrwr*pa.frrwr[4]+ Pat[pat].Frrrr*pa.frrrr[4]+ Pat[pat].Frrdr*pa.frrdr[4] + Pat[pat].Fdrww*pa.fdrww[4]+ Pat[pat].Fdrwd*pa.fdrwd[4]+ Pat[pat].Fdrdd*pa.fdrdd[4]+ Pat[pat].Fdrwr*pa.fdrwr[4]+ Pat[pat].Fdrrr*pa.fdrrr[4]+ Pat[pat].Fdrdr*pa.fdrdr[4]));
	Pat[pat].Jdd[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[5]+ Pat[pat].Fwwwd*pa.fwwwd[5]+ Pat[pat].Fwwdd*pa.fwwdd[5]+ Pat[pat].Fwwwr*pa.fwwwr[5]+ Pat[pat].Fwwrr*pa.fwwrr[5]+ Pat[pat].Fwwdr*pa.fwwdr[5]+ Pat[pat].fFwwww*pa.ffwwww[5]+ Pat[pat].fFwwwd*pa.ffwwwd[5]+ Pat[pat].fFwwdd*pa.ffwwdd[5]+ Pat[pat].fFwwwr*pa.ffwwwr[5]+ Pat[pat].fFwwrr*pa.ffwwrr[5]+ Pat[pat].fFwwdr*pa.ffwwdr[5]+ Pat[pat].mFwwww*pa.mfwwww[5]+ Pat[pat].mFwwwd*pa.mfwwwd[5]+ Pat[pat].mFwwdd*pa.mfwwdd[5]+ Pat[pat].mFwwwr*pa.mfwwwr[5]+ Pat[pat].mFwwrr*pa.mfwwrr[5]+ Pat[pat].mFwwdr*pa.mfwwdr[5]+ Pat[pat].bFwwww*pa.bfwwww[5]+ Pat[pat].bFwwwd*pa.bfwwwd[5]+ Pat[pat].bFwwdd*pa.bfwwdd[5]+ Pat[pat].bFwwwr*pa.bfwwwr[5]+ Pat[pat].bFwwrr*pa.bfwwrr[5]+ Pat[pat].bFwwdr*pa.bfwwdr[5] +Pat[pat].fFwdww*pa.ffwdww[5]+ Pat[pat].fFwdwd*pa.ffwdwd[5]+ Pat[pat].fFwddd*pa.ffwddd[5]+ Pat[pat].fFwdwr*pa.ffwdwr[5]+ Pat[pat].fFwdrr*pa.ffwdrr[5]+ Pat[pat].fFwddr*pa.ffwddr[5] +Pat[pat].mFwdww*pa.mfwdww[5]+ Pat[pat].mFwdwd*pa.mfwdwd[5]+ Pat[pat].mFwddd*pa.mfwddd[5]+ Pat[pat].mFwdwr*pa.mfwdwr[5]+ Pat[pat].mFwdrr*pa.mfwdrr[5]+ Pat[pat].mFwddr*pa.mfwddr[5] +Pat[pat].bFwdww*pa.bfwdww[5]+ Pat[pat].bFwdwd*pa.bfwdwd[5]+ Pat[pat].bFwddd*pa.bfwddd[5]+ Pat[pat].bFwdwr*pa.bfwdwr[5]+ Pat[pat].bFwdrr*pa.bfwdrr[5]+ Pat[pat].bFwddr*pa.bfwddr[5] + Pat[pat].Fddww*pa.fddww[5]+ Pat[pat].Fddwd*pa.fddwd[5]+ Pat[pat].Fdddd*pa.fdddd[5]+ Pat[pat].Fddwr*pa.fddwr[5]+ Pat[pat].Fddrr*pa.fddrr[5]+ Pat[pat].Fdddr*pa.fdddr[5] + Pat[pat].Fwrww*pa.fwrww[5]+ Pat[pat].Fwrwd*pa.fwrwd[5]+ Pat[pat].Fwrdd*pa.fwrdd[5]+ Pat[pat].Fwrwr*pa.fwrwr[5]+ Pat[pat].Fwrrr*pa.fwrrr[5]+ Pat[pat].Fwrdr*pa.fwrdr[5] + Pat[pat].fFwrww*pa.ffwrww[5]+ Pat[pat].fFwrwd*pa.ffwrwd[5]+ Pat[pat].fFwrdd*pa.ffwrdd[5]+ Pat[pat].fFwrwr*pa.ffwrwr[5]+ Pat[pat].fFwrrr*pa.ffwrrr[5]+ Pat[pat].fFwrdr*pa.ffwrdr[5] + Pat[pat].mFwrww*pa.mfwrww[5]+ Pat[pat].mFwrwd*pa.mfwrwd[5]+ Pat[pat].mFwrdd*pa.mfwrdd[5]+ Pat[pat].mFwrwr*pa.mfwrwr[5]+ Pat[pat].mFwrrr*pa.mfwrrr[5]+ Pat[pat].mFwrdr*pa.mfwrdr[5] + Pat[pat].bFwrww*pa.bfwrww[5]+ Pat[pat].bFwrwd*pa.bfwrwd[5]+ Pat[pat].bFwrdd*pa.bfwrdd[5]+ Pat[pat].bFwrwr*pa.bfwrwr[5]+ Pat[pat].bFwrrr*pa.bfwrrr[5]+ Pat[pat].bFwrdr*pa.bfwrdr[5] + Pat[pat].Frrww*pa.frrww[5]+ Pat[pat].Frrwd*pa.frrwd[5]+ Pat[pat].Frrdd*pa.frrdd[5]+ Pat[pat].Frrwr*pa.frrwr[5]+ Pat[pat].Frrrr*pa.frrrr[5]+ Pat[pat].Frrdr*pa.frrdr[5] + Pat[pat].Fdrww*pa.fdrww[5]+ Pat[pat].Fdrwd*pa.fdrwd[5]+ Pat[pat].Fdrdd*pa.fdrdd[5]+ Pat[pat].Fdrwr*pa.fdrwr[5]+ Pat[pat].Fdrrr*pa.fdrrr[5]+ Pat[pat].Fdrdr*pa.fdrdr[5]));
	Pat[pat].Jwr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[6]+ Pat[pat].Fwwwd*pa.fwwwd[6]+ Pat[pat].Fwwdd*pa.fwwdd[6]+ Pat[pat].Fwwwr*pa.fwwwr[6]+ Pat[pat].Fwwrr*pa.fwwrr[6]+ Pat[pat].Fwwdr*pa.fwwdr[6]+ Pat[pat].fFwwww*pa.ffwwww[6]+ Pat[pat].fFwwwd*pa.ffwwwd[6]+ Pat[pat].fFwwdd*pa.ffwwdd[6]+ Pat[pat].fFwwwr*pa.ffwwwr[6]+ Pat[pat].fFwwrr*pa.ffwwrr[6]+ Pat[pat].fFwwdr*pa.ffwwdr[6]+ Pat[pat].mFwwww*pa.mfwwww[6]+ Pat[pat].mFwwwd*pa.mfwwwd[6]+ Pat[pat].mFwwdd*pa.mfwwdd[6]+ Pat[pat].mFwwwr*pa.mfwwwr[6]+ Pat[pat].mFwwrr*pa.mfwwrr[6]+ Pat[pat].mFwwdr*pa.mfwwdr[6]+ Pat[pat].bFwwww*pa.bfwwww[6]+ Pat[pat].bFwwwd*pa.bfwwwd[6]+ Pat[pat].bFwwdd*pa.bfwwdd[6]+ Pat[pat].bFwwwr*pa.bfwwwr[6]+ Pat[pat].bFwwrr*pa.bfwwrr[6]+ Pat[pat].bFwwdr*pa.bfwwdr[6] +Pat[pat].fFwdww*pa.ffwdww[6]+ Pat[pat].fFwdwd*pa.ffwdwd[6]+ Pat[pat].fFwddd*pa.ffwddd[6]+ Pat[pat].fFwdwr*pa.ffwdwr[6]+ Pat[pat].fFwdrr*pa.ffwdrr[6]+ Pat[pat].fFwddr*pa.ffwddr[6] +Pat[pat].mFwdww*pa.mfwdww[6]+ Pat[pat].mFwdwd*pa.mfwdwd[6]+ Pat[pat].mFwddd*pa.mfwddd[6]+ Pat[pat].mFwdwr*pa.mfwdwr[6]+ Pat[pat].mFwdrr*pa.mfwdrr[6]+ Pat[pat].mFwddr*pa.mfwddr[6] +Pat[pat].bFwdww*pa.bfwdww[6]+ Pat[pat].bFwdwd*pa.bfwdwd[6]+ Pat[pat].bFwddd*pa.bfwddd[6]+ Pat[pat].bFwdwr*pa.bfwdwr[6]+ Pat[pat].bFwdrr*pa.bfwdrr[6]+ Pat[pat].bFwddr*pa.bfwddr[6] + Pat[pat].Fddww*pa.fddww[6]+ Pat[pat].Fddwd*pa.fddwd[6]+ Pat[pat].Fdddd*pa.fdddd[6]+ Pat[pat].Fddwr*pa.fddwr[6]+ Pat[pat].Fddrr*pa.fddrr[6]+ Pat[pat].Fdddr*pa.fdddr[6] + Pat[pat].Fwrww*pa.fwrww[6]+ Pat[pat].Fwrwd*pa.fwrwd[6]+ Pat[pat].Fwrdd*pa.fwrdd[6]+ Pat[pat].Fwrwr*pa.fwrwr[6]+ Pat[pat].Fwrrr*pa.fwrrr[6]+ Pat[pat].Fwrdr*pa.fwrdr[6] + Pat[pat].fFwrww*pa.ffwrww[6]+ Pat[pat].fFwrwd*pa.ffwrwd[6]+ Pat[pat].fFwrdd*pa.ffwrdd[6]+ Pat[pat].fFwrwr*pa.ffwrwr[6]+ Pat[pat].fFwrrr*pa.ffwrrr[6]+ Pat[pat].fFwrdr*pa.ffwrdr[6] + Pat[pat].mFwrww*pa.mfwrww[6]+ Pat[pat].mFwrwd*pa.mfwrwd[6]+ Pat[pat].mFwrdd*pa.mfwrdd[6]+ Pat[pat].mFwrwr*pa.mfwrwr[6]+ Pat[pat].mFwrrr*pa.mfwrrr[6]+ Pat[pat].mFwrdr*pa.mfwrdr[6] + Pat[pat].bFwrww*pa.bfwrww[6]+ Pat[pat].bFwrwd*pa.bfwrwd[6]+ Pat[pat].bFwrdd*pa.bfwrdd[6]+ Pat[pat].bFwrwr*pa.bfwrwr[6]+ Pat[pat].bFwrrr*pa.bfwrrr[6]+ Pat[pat].bFwrdr*pa.bfwrdr[6] + Pat[pat].Frrww*pa.frrww[6]+ Pat[pat].Frrwd*pa.frrwd[6]+ Pat[pat].Frrdd*pa.frrdd[6]+ Pat[pat].Frrwr*pa.frrwr[6]+ Pat[pat].Frrrr*pa.frrrr[6]+ Pat[pat].Frrdr*pa.frrdr[6] + Pat[pat].Fdrww*pa.fdrww[6]+ Pat[pat].Fdrwd*pa.fdrwd[6]+ Pat[pat].Fdrdd*pa.fdrdd[6]+ Pat[pat].Fdrwr*pa.fdrwr[6]+ Pat[pat].Fdrrr*pa.fdrrr[6]+ Pat[pat].Fdrdr*pa.fdrdr[6]));
	Pat[pat].fJwr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[7]+ Pat[pat].Fwwwd*pa.fwwwd[7]+ Pat[pat].Fwwdd*pa.fwwdd[7]+ Pat[pat].Fwwwr*pa.fwwwr[7]+ Pat[pat].Fwwrr*pa.fwwrr[7]+ Pat[pat].Fwwdr*pa.fwwdr[7]+ Pat[pat].fFwwww*pa.ffwwww[7]+ Pat[pat].fFwwwd*pa.ffwwwd[7]+ Pat[pat].fFwwdd*pa.ffwwdd[7]+ Pat[pat].fFwwwr*pa.ffwwwr[7]+ Pat[pat].fFwwrr*pa.ffwwrr[7]+ Pat[pat].fFwwdr*pa.ffwwdr[7]+ Pat[pat].mFwwww*pa.mfwwww[7]+ Pat[pat].mFwwwd*pa.mfwwwd[7]+ Pat[pat].mFwwdd*pa.mfwwdd[7]+ Pat[pat].mFwwwr*pa.mfwwwr[7]+ Pat[pat].mFwwrr*pa.mfwwrr[7]+ Pat[pat].mFwwdr*pa.mfwwdr[7]+ Pat[pat].bFwwww*pa.bfwwww[7]+ Pat[pat].bFwwwd*pa.bfwwwd[7]+ Pat[pat].bFwwdd*pa.bfwwdd[7]+ Pat[pat].bFwwwr*pa.bfwwwr[7]+ Pat[pat].bFwwrr*pa.bfwwrr[7]+ Pat[pat].bFwwdr*pa.bfwwdr[7] +Pat[pat].fFwdww*pa.ffwdww[7]+ Pat[pat].fFwdwd*pa.ffwdwd[7]+ Pat[pat].fFwddd*pa.ffwddd[7]+ Pat[pat].fFwdwr*pa.ffwdwr[7]+ Pat[pat].fFwdrr*pa.ffwdrr[7]+ Pat[pat].fFwddr*pa.ffwddr[7] +Pat[pat].mFwdww*pa.mfwdww[7]+ Pat[pat].mFwdwd*pa.mfwdwd[7]+ Pat[pat].mFwddd*pa.mfwddd[7]+ Pat[pat].mFwdwr*pa.mfwdwr[7]+ Pat[pat].mFwdrr*pa.mfwdrr[7]+ Pat[pat].mFwddr*pa.mfwddr[7] +Pat[pat].bFwdww*pa.bfwdww[7]+ Pat[pat].bFwdwd*pa.bfwdwd[7]+ Pat[pat].bFwddd*pa.bfwddd[7]+ Pat[pat].bFwdwr*pa.bfwdwr[7]+ Pat[pat].bFwdrr*pa.bfwdrr[7]+ Pat[pat].bFwddr*pa.bfwddr[7] + Pat[pat].Fddww*pa.fddww[7]+ Pat[pat].Fddwd*pa.fddwd[7]+ Pat[pat].Fdddd*pa.fdddd[7]+ Pat[pat].Fddwr*pa.fddwr[7]+ Pat[pat].Fddrr*pa.fddrr[7]+ Pat[pat].Fdddr*pa.fdddr[7] + Pat[pat].Fwrww*pa.fwrww[7]+ Pat[pat].Fwrwd*pa.fwrwd[7]+ Pat[pat].Fwrdd*pa.fwrdd[7]+ Pat[pat].Fwrwr*pa.fwrwr[7]+ Pat[pat].Fwrrr*pa.fwrrr[7]+ Pat[pat].Fwrdr*pa.fwrdr[7] + Pat[pat].fFwrww*pa.ffwrww[7]+ Pat[pat].fFwrwd*pa.ffwrwd[7]+ Pat[pat].fFwrdd*pa.ffwrdd[7]+ Pat[pat].fFwrwr*pa.ffwrwr[7]+ Pat[pat].fFwrrr*pa.ffwrrr[7]+ Pat[pat].fFwrdr*pa.ffwrdr[7] + Pat[pat].mFwrww*pa.mfwrww[7]+ Pat[pat].mFwrwd*pa.mfwrwd[7]+ Pat[pat].mFwrdd*pa.mfwrdd[7]+ Pat[pat].mFwrwr*pa.mfwrwr[7]+ Pat[pat].mFwrrr*pa.mfwrrr[7]+ Pat[pat].mFwrdr*pa.mfwrdr[7] + Pat[pat].bFwrww*pa.bfwrww[7]+ Pat[pat].bFwrwd*pa.bfwrwd[7]+ Pat[pat].bFwrdd*pa.bfwrdd[7]+ Pat[pat].bFwrwr*pa.bfwrwr[7]+ Pat[pat].bFwrrr*pa.bfwrrr[7]+ Pat[pat].bFwrdr*pa.bfwrdr[7] + Pat[pat].Frrww*pa.frrww[7]+ Pat[pat].Frrwd*pa.frrwd[7]+ Pat[pat].Frrdd*pa.frrdd[7]+ Pat[pat].Frrwr*pa.frrwr[7]+ Pat[pat].Frrrr*pa.frrrr[7]+ Pat[pat].Frrdr*pa.frrdr[7] + Pat[pat].Fdrww*pa.fdrww[7]+ Pat[pat].Fdrwd*pa.fdrwd[7]+ Pat[pat].Fdrdd*pa.fdrdd[7]+ Pat[pat].Fdrwr*pa.fdrwr[7]+ Pat[pat].Fdrrr*pa.fdrrr[7]+ Pat[pat].Fdrdr*pa.fdrdr[7]));
	Pat[pat].mJwr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[8]+ Pat[pat].Fwwwd*pa.fwwwd[8]+ Pat[pat].Fwwdd*pa.fwwdd[8]+ Pat[pat].Fwwwr*pa.fwwwr[8]+ Pat[pat].Fwwrr*pa.fwwrr[8]+ Pat[pat].Fwwdr*pa.fwwdr[8]+ Pat[pat].fFwwww*pa.ffwwww[8]+ Pat[pat].fFwwwd*pa.ffwwwd[8]+ Pat[pat].fFwwdd*pa.ffwwdd[8]+ Pat[pat].fFwwwr*pa.ffwwwr[8]+ Pat[pat].fFwwrr*pa.ffwwrr[8]+ Pat[pat].fFwwdr*pa.ffwwdr[8]+ Pat[pat].mFwwww*pa.mfwwww[8]+ Pat[pat].mFwwwd*pa.mfwwwd[8]+ Pat[pat].mFwwdd*pa.mfwwdd[8]+ Pat[pat].mFwwwr*pa.mfwwwr[8]+ Pat[pat].mFwwrr*pa.mfwwrr[8]+ Pat[pat].mFwwdr*pa.mfwwdr[8]+ Pat[pat].bFwwww*pa.bfwwww[8]+ Pat[pat].bFwwwd*pa.bfwwwd[8]+ Pat[pat].bFwwdd*pa.bfwwdd[8]+ Pat[pat].bFwwwr*pa.bfwwwr[8]+ Pat[pat].bFwwrr*pa.bfwwrr[8]+ Pat[pat].bFwwdr*pa.bfwwdr[8] +Pat[pat].fFwdww*pa.ffwdww[8]+ Pat[pat].fFwdwd*pa.ffwdwd[8]+ Pat[pat].fFwddd*pa.ffwddd[8]+ Pat[pat].fFwdwr*pa.ffwdwr[8]+ Pat[pat].fFwdrr*pa.ffwdrr[8]+ Pat[pat].fFwddr*pa.ffwddr[8] +Pat[pat].mFwdww*pa.mfwdww[8]+ Pat[pat].mFwdwd*pa.mfwdwd[8]+ Pat[pat].mFwddd*pa.mfwddd[8]+ Pat[pat].mFwdwr*pa.mfwdwr[8]+ Pat[pat].mFwdrr*pa.mfwdrr[8]+ Pat[pat].mFwddr*pa.mfwddr[8] +Pat[pat].bFwdww*pa.bfwdww[8]+ Pat[pat].bFwdwd*pa.bfwdwd[8]+ Pat[pat].bFwddd*pa.bfwddd[8]+ Pat[pat].bFwdwr*pa.bfwdwr[8]+ Pat[pat].bFwdrr*pa.bfwdrr[8]+ Pat[pat].bFwddr*pa.bfwddr[8] + Pat[pat].Fddww*pa.fddww[8]+ Pat[pat].Fddwd*pa.fddwd[8]+ Pat[pat].Fdddd*pa.fdddd[8]+ Pat[pat].Fddwr*pa.fddwr[8]+ Pat[pat].Fddrr*pa.fddrr[8]+ Pat[pat].Fdddr*pa.fdddr[8] + Pat[pat].Fwrww*pa.fwrww[8]+ Pat[pat].Fwrwd*pa.fwrwd[8]+ Pat[pat].Fwrdd*pa.fwrdd[8]+ Pat[pat].Fwrwr*pa.fwrwr[8]+ Pat[pat].Fwrrr*pa.fwrrr[8]+ Pat[pat].Fwrdr*pa.fwrdr[8] + Pat[pat].fFwrww*pa.ffwrww[8]+ Pat[pat].fFwrwd*pa.ffwrwd[8]+ Pat[pat].fFwrdd*pa.ffwrdd[8]+ Pat[pat].fFwrwr*pa.ffwrwr[8]+ Pat[pat].fFwrrr*pa.ffwrrr[8]+ Pat[pat].fFwrdr*pa.ffwrdr[8] + Pat[pat].mFwrww*pa.mfwrww[8]+ Pat[pat].mFwrwd*pa.mfwrwd[8]+ Pat[pat].mFwrdd*pa.mfwrdd[8]+ Pat[pat].mFwrwr*pa.mfwrwr[8]+ Pat[pat].mFwrrr*pa.mfwrrr[8]+ Pat[pat].mFwrdr*pa.mfwrdr[8] + Pat[pat].bFwrww*pa.bfwrww[8]+ Pat[pat].bFwrwd*pa.bfwrwd[8]+ Pat[pat].bFwrdd*pa.bfwrdd[8]+ Pat[pat].bFwrwr*pa.bfwrwr[8]+ Pat[pat].bFwrrr*pa.bfwrrr[8]+ Pat[pat].bFwrdr*pa.bfwrdr[8] + Pat[pat].Frrww*pa.frrww[8]+ Pat[pat].Frrwd*pa.frrwd[8]+ Pat[pat].Frrdd*pa.frrdd[8]+ Pat[pat].Frrwr*pa.frrwr[8]+ Pat[pat].Frrrr*pa.frrrr[8]+ Pat[pat].Frrdr*pa.frrdr[8] + Pat[pat].Fdrww*pa.fdrww[8]+ Pat[pat].Fdrwd*pa.fdrwd[8]+ Pat[pat].Fdrdd*pa.fdrdd[8]+ Pat[pat].Fdrwr*pa.fdrwr[8]+ Pat[pat].Fdrrr*pa.fdrrr[8]+ Pat[pat].Fdrdr*pa.fdrdr[8]));
	Pat[pat].Jrr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[9]+ Pat[pat].Fwwwd*pa.fwwwd[9]+ Pat[pat].Fwwdd*pa.fwwdd[9]+ Pat[pat].Fwwwr*pa.fwwwr[9]+ Pat[pat].Fwwrr*pa.fwwrr[9]+ Pat[pat].Fwwdr*pa.fwwdr[9]+ Pat[pat].fFwwww*pa.ffwwww[9]+ Pat[pat].fFwwwd*pa.ffwwwd[9]+ Pat[pat].fFwwdd*pa.ffwwdd[9]+ Pat[pat].fFwwwr*pa.ffwwwr[9]+ Pat[pat].fFwwrr*pa.ffwwrr[9]+ Pat[pat].fFwwdr*pa.ffwwdr[9]+ Pat[pat].mFwwww*pa.mfwwww[9]+ Pat[pat].mFwwwd*pa.mfwwwd[9]+ Pat[pat].mFwwdd*pa.mfwwdd[9]+ Pat[pat].mFwwwr*pa.mfwwwr[9]+ Pat[pat].mFwwrr*pa.mfwwrr[9]+ Pat[pat].mFwwdr*pa.mfwwdr[9]+ Pat[pat].bFwwww*pa.bfwwww[9]+ Pat[pat].bFwwwd*pa.bfwwwd[9]+ Pat[pat].bFwwdd*pa.bfwwdd[9]+ Pat[pat].bFwwwr*pa.bfwwwr[9]+ Pat[pat].bFwwrr*pa.bfwwrr[9]+ Pat[pat].bFwwdr*pa.bfwwdr[9] +Pat[pat].fFwdww*pa.ffwdww[9]+ Pat[pat].fFwdwd*pa.ffwdwd[9]+ Pat[pat].fFwddd*pa.ffwddd[9]+ Pat[pat].fFwdwr*pa.ffwdwr[9]+ Pat[pat].fFwdrr*pa.ffwdrr[9]+ Pat[pat].fFwddr*pa.ffwddr[9] +Pat[pat].mFwdww*pa.mfwdww[9]+ Pat[pat].mFwdwd*pa.mfwdwd[9]+ Pat[pat].mFwddd*pa.mfwddd[9]+ Pat[pat].mFwdwr*pa.mfwdwr[9]+ Pat[pat].mFwdrr*pa.mfwdrr[9]+ Pat[pat].mFwddr*pa.mfwddr[9] +Pat[pat].bFwdww*pa.bfwdww[9]+ Pat[pat].bFwdwd*pa.bfwdwd[9]+ Pat[pat].bFwddd*pa.bfwddd[9]+ Pat[pat].bFwdwr*pa.bfwdwr[9]+ Pat[pat].bFwdrr*pa.bfwdrr[9]+ Pat[pat].bFwddr*pa.bfwddr[9] + Pat[pat].Fddww*pa.fddww[9]+ Pat[pat].Fddwd*pa.fddwd[9]+ Pat[pat].Fdddd*pa.fdddd[9]+ Pat[pat].Fddwr*pa.fddwr[9]+ Pat[pat].Fddrr*pa.fddrr[9]+ Pat[pat].Fdddr*pa.fdddr[9] + Pat[pat].Fwrww*pa.fwrww[9]+ Pat[pat].Fwrwd*pa.fwrwd[9]+ Pat[pat].Fwrdd*pa.fwrdd[9]+ Pat[pat].Fwrwr*pa.fwrwr[9]+ Pat[pat].Fwrrr*pa.fwrrr[9]+ Pat[pat].Fwrdr*pa.fwrdr[9] + Pat[pat].fFwrww*pa.ffwrww[9]+ Pat[pat].fFwrwd*pa.ffwrwd[9]+ Pat[pat].fFwrdd*pa.ffwrdd[9]+ Pat[pat].fFwrwr*pa.ffwrwr[9]+ Pat[pat].fFwrrr*pa.ffwrrr[9]+ Pat[pat].fFwrdr*pa.ffwrdr[9] + Pat[pat].mFwrww*pa.mfwrww[9]+ Pat[pat].mFwrwd*pa.mfwrwd[9]+ Pat[pat].mFwrdd*pa.mfwrdd[9]+ Pat[pat].mFwrwr*pa.mfwrwr[9]+ Pat[pat].mFwrrr*pa.mfwrrr[9]+ Pat[pat].mFwrdr*pa.mfwrdr[9] + Pat[pat].bFwrww*pa.bfwrww[9]+ Pat[pat].bFwrwd*pa.bfwrwd[9]+ Pat[pat].bFwrdd*pa.bfwrdd[9]+ Pat[pat].bFwrwr*pa.bfwrwr[9]+ Pat[pat].bFwrrr*pa.bfwrrr[9]+ Pat[pat].bFwrdr*pa.bfwrdr[9] + Pat[pat].Frrww*pa.frrww[9]+ Pat[pat].Frrwd*pa.frrwd[9]+ Pat[pat].Frrdd*pa.frrdd[9]+ Pat[pat].Frrwr*pa.frrwr[9]+ Pat[pat].Frrrr*pa.frrrr[9]+ Pat[pat].Frrdr*pa.frrdr[9] + Pat[pat].Fdrww*pa.fdrww[9]+ Pat[pat].Fdrwd*pa.fdrwd[9]+ Pat[pat].Fdrdd*pa.fdrdd[9]+ Pat[pat].Fdrwr*pa.fdrwr[9]+ Pat[pat].Fdrrr*pa.fdrrr[9]+ Pat[pat].Fdrdr*pa.fdrdr[9]));
	Pat[pat].Jdr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[10]+ Pat[pat].Fwwwd*pa.fwwwd[10]+ Pat[pat].Fwwdd*pa.fwwdd[10]+ Pat[pat].Fwwwr*pa.fwwwr[10]+ Pat[pat].Fwwrr*pa.fwwrr[10]+ Pat[pat].Fwwdr*pa.fwwdr[10]+ Pat[pat].fFwwww*pa.ffwwww[10]+ Pat[pat].fFwwwd*pa.ffwwwd[10]+ Pat[pat].fFwwdd*pa.ffwwdd[10]+ Pat[pat].fFwwwr*pa.ffwwwr[10]+ Pat[pat].fFwwrr*pa.ffwwrr[10]+ Pat[pat].fFwwdr*pa.ffwwdr[10]+ Pat[pat].mFwwww*pa.mfwwww[10]+ Pat[pat].mFwwwd*pa.mfwwwd[10]+ Pat[pat].mFwwdd*pa.mfwwdd[10]+ Pat[pat].mFwwwr*pa.mfwwwr[10]+ Pat[pat].mFwwrr*pa.mfwwrr[10]+ Pat[pat].mFwwdr*pa.mfwwdr[10]+ Pat[pat].bFwwww*pa.bfwwww[10]+ Pat[pat].bFwwwd*pa.bfwwwd[10]+ Pat[pat].bFwwdd*pa.bfwwdd[10]+ Pat[pat].bFwwwr*pa.bfwwwr[10]+ Pat[pat].bFwwrr*pa.bfwwrr[10]+ Pat[pat].bFwwdr*pa.bfwwdr[10] +Pat[pat].fFwdww*pa.ffwdww[10]+ Pat[pat].fFwdwd*pa.ffwdwd[10]+ Pat[pat].fFwddd*pa.ffwddd[10]+ Pat[pat].fFwdwr*pa.ffwdwr[10]+ Pat[pat].fFwdrr*pa.ffwdrr[10]+ Pat[pat].fFwddr*pa.ffwddr[10] +Pat[pat].mFwdww*pa.mfwdww[10]+ Pat[pat].mFwdwd*pa.mfwdwd[10]+ Pat[pat].mFwddd*pa.mfwddd[10]+ Pat[pat].mFwdwr*pa.mfwdwr[10]+ Pat[pat].mFwdrr*pa.mfwdrr[10]+ Pat[pat].mFwddr*pa.mfwddr[10] +Pat[pat].bFwdww*pa.bfwdww[10]+ Pat[pat].bFwdwd*pa.bfwdwd[10]+ Pat[pat].bFwddd*pa.bfwddd[10]+ Pat[pat].bFwdwr*pa.bfwdwr[10]+ Pat[pat].bFwdrr*pa.bfwdrr[10]+ Pat[pat].bFwddr*pa.bfwddr[10] + Pat[pat].Fddww*pa.fddww[10]+ Pat[pat].Fddwd*pa.fddwd[10]+ Pat[pat].Fdddd*pa.fdddd[10]+ Pat[pat].Fddwr*pa.fddwr[10]+ Pat[pat].Fddrr*pa.fddrr[10]+ Pat[pat].Fdddr*pa.fdddr[10] + Pat[pat].Fwrww*pa.fwrww[10]+ Pat[pat].Fwrwd*pa.fwrwd[10]+ Pat[pat].Fwrdd*pa.fwrdd[10]+ Pat[pat].Fwrwr*pa.fwrwr[10]+ Pat[pat].Fwrrr*pa.fwrrr[10]+ Pat[pat].Fwrdr*pa.fwrdr[10] + Pat[pat].fFwrww*pa.ffwrww[10]+ Pat[pat].fFwrwd*pa.ffwrwd[10]+ Pat[pat].fFwrdd*pa.ffwrdd[10]+ Pat[pat].fFwrwr*pa.ffwrwr[10]+ Pat[pat].fFwrrr*pa.ffwrrr[10]+ Pat[pat].fFwrdr*pa.ffwrdr[10] + Pat[pat].mFwrww*pa.mfwrww[10]+ Pat[pat].mFwrwd*pa.mfwrwd[10]+ Pat[pat].mFwrdd*pa.mfwrdd[10]+ Pat[pat].mFwrwr*pa.mfwrwr[10]+ Pat[pat].mFwrrr*pa.mfwrrr[10]+ Pat[pat].mFwrdr*pa.mfwrdr[10] + Pat[pat].bFwrww*pa.bfwrww[10]+ Pat[pat].bFwrwd*pa.bfwrwd[10]+ Pat[pat].bFwrdd*pa.bfwrdd[10]+ Pat[pat].bFwrwr*pa.bfwrwr[10]+ Pat[pat].bFwrrr*pa.bfwrrr[10]+ Pat[pat].bFwrdr*pa.bfwrdr[10] + Pat[pat].Frrww*pa.frrww[10]+ Pat[pat].Frrwd*pa.frrwd[10]+ Pat[pat].Frrdd*pa.frrdd[10]+ Pat[pat].Frrwr*pa.frrwr[10]+ Pat[pat].Frrrr*pa.frrrr[10]+ Pat[pat].Frrdr*pa.frrdr[10] + Pat[pat].Fdrww*pa.fdrww[10]+ Pat[pat].Fdrwd*pa.fdrwd[10]+ Pat[pat].Fdrdd*pa.fdrdd[10]+ Pat[pat].Fdrwr*pa.fdrwr[10]+ Pat[pat].Fdrrr*pa.fdrrr[10]+ Pat[pat].Fdrdr*pa.fdrdr[10]));
	Pat[pat].bJww[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[11]+ Pat[pat].Fwwwd*pa.fwwwd[11]+ Pat[pat].Fwwdd*pa.fwwdd[11]+ Pat[pat].Fwwwr*pa.fwwwr[11]+ Pat[pat].Fwwrr*pa.fwwrr[11]+ Pat[pat].Fwwdr*pa.fwwdr[11]+ Pat[pat].fFwwww*pa.ffwwww[11]+ Pat[pat].fFwwwd*pa.ffwwwd[11]+ Pat[pat].fFwwdd*pa.ffwwdd[11]+ Pat[pat].fFwwwr*pa.ffwwwr[11]+ Pat[pat].fFwwrr*pa.ffwwrr[11]+ Pat[pat].fFwwdr*pa.ffwwdr[11]+ Pat[pat].mFwwww*pa.mfwwww[11]+ Pat[pat].mFwwwd*pa.mfwwwd[11]+ Pat[pat].mFwwdd*pa.mfwwdd[11]+ Pat[pat].mFwwwr*pa.mfwwwr[11]+ Pat[pat].mFwwrr*pa.mfwwrr[11]+ Pat[pat].mFwwdr*pa.mfwwdr[11]+ Pat[pat].bFwwww*pa.bfwwww[11]+ Pat[pat].bFwwwd*pa.bfwwwd[11]+ Pat[pat].bFwwdd*pa.bfwwdd[11]+ Pat[pat].bFwwwr*pa.bfwwwr[11]+ Pat[pat].bFwwrr*pa.bfwwrr[11]+ Pat[pat].bFwwdr*pa.bfwwdr[11] +Pat[pat].fFwdww*pa.ffwdww[11]+ Pat[pat].fFwdwd*pa.ffwdwd[11]+ Pat[pat].fFwddd*pa.ffwddd[11]+ Pat[pat].fFwdwr*pa.ffwdwr[11]+ Pat[pat].fFwdrr*pa.ffwdrr[11]+ Pat[pat].fFwddr*pa.ffwddr[11] +Pat[pat].mFwdww*pa.mfwdww[11]+ Pat[pat].mFwdwd*pa.mfwdwd[11]+ Pat[pat].mFwddd*pa.mfwddd[11]+ Pat[pat].mFwdwr*pa.mfwdwr[11]+ Pat[pat].mFwdrr*pa.mfwdrr[11]+ Pat[pat].mFwddr*pa.mfwddr[11] +Pat[pat].bFwdww*pa.bfwdww[11]+ Pat[pat].bFwdwd*pa.bfwdwd[11]+ Pat[pat].bFwddd*pa.bfwddd[11]+ Pat[pat].bFwdwr*pa.bfwdwr[11]+ Pat[pat].bFwdrr*pa.bfwdrr[11]+ Pat[pat].bFwddr*pa.bfwddr[11] + Pat[pat].Fddww*pa.fddww[11]+ Pat[pat].Fddwd*pa.fddwd[11]+ Pat[pat].Fdddd*pa.fdddd[11]+ Pat[pat].Fddwr*pa.fddwr[11]+ Pat[pat].Fddrr*pa.fddrr[11]+ Pat[pat].Fdddr*pa.fdddr[11] + Pat[pat].Fwrww*pa.fwrww[11]+ Pat[pat].Fwrwd*pa.fwrwd[11]+ Pat[pat].Fwrdd*pa.fwrdd[11]+ Pat[pat].Fwrwr*pa.fwrwr[11]+ Pat[pat].Fwrrr*pa.fwrrr[11]+ Pat[pat].Fwrdr*pa.fwrdr[11] + Pat[pat].fFwrww*pa.ffwrww[11]+ Pat[pat].fFwrwd*pa.ffwrwd[11]+ Pat[pat].fFwrdd*pa.ffwrdd[11]+ Pat[pat].fFwrwr*pa.ffwrwr[11]+ Pat[pat].fFwrrr*pa.ffwrrr[11]+ Pat[pat].fFwrdr*pa.ffwrdr[11] + Pat[pat].mFwrww*pa.mfwrww[11]+ Pat[pat].mFwrwd*pa.mfwrwd[11]+ Pat[pat].mFwrdd*pa.mfwrdd[11]+ Pat[pat].mFwrwr*pa.mfwrwr[11]+ Pat[pat].mFwrrr*pa.mfwrrr[11]+ Pat[pat].mFwrdr*pa.mfwrdr[11] + Pat[pat].bFwrww*pa.bfwrww[11]+ Pat[pat].bFwrwd*pa.bfwrwd[11]+ Pat[pat].bFwrdd*pa.bfwrdd[11]+ Pat[pat].bFwrwr*pa.bfwrwr[11]+ Pat[pat].bFwrrr*pa.bfwrrr[11]+ Pat[pat].bFwrdr*pa.bfwrdr[11] + Pat[pat].Frrww*pa.frrww[11]+ Pat[pat].Frrwd*pa.frrwd[11]+ Pat[pat].Frrdd*pa.frrdd[11]+ Pat[pat].Frrwr*pa.frrwr[11]+ Pat[pat].Frrrr*pa.frrrr[11]+ Pat[pat].Frrdr*pa.frrdr[11] + Pat[pat].Fdrww*pa.fdrww[11]+ Pat[pat].Fdrwd*pa.fdrwd[11]+ Pat[pat].Fdrdd*pa.fdrdd[11]+ Pat[pat].Fdrwr*pa.fdrwr[11]+ Pat[pat].Fdrrr*pa.fdrrr[11]+ Pat[pat].Fdrdr*pa.fdrdr[11]));
	Pat[pat].bJwd[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[12]+ Pat[pat].Fwwwd*pa.fwwwd[12]+ Pat[pat].Fwwdd*pa.fwwdd[12]+ Pat[pat].Fwwwr*pa.fwwwr[12]+ Pat[pat].Fwwrr*pa.fwwrr[12]+ Pat[pat].Fwwdr*pa.fwwdr[12]+ Pat[pat].fFwwww*pa.ffwwww[12]+ Pat[pat].fFwwwd*pa.ffwwwd[12]+ Pat[pat].fFwwdd*pa.ffwwdd[12]+ Pat[pat].fFwwwr*pa.ffwwwr[12]+ Pat[pat].fFwwrr*pa.ffwwrr[12]+ Pat[pat].fFwwdr*pa.ffwwdr[12]+ Pat[pat].mFwwww*pa.mfwwww[12]+ Pat[pat].mFwwwd*pa.mfwwwd[12]+ Pat[pat].mFwwdd*pa.mfwwdd[12]+ Pat[pat].mFwwwr*pa.mfwwwr[12]+ Pat[pat].mFwwrr*pa.mfwwrr[12]+ Pat[pat].mFwwdr*pa.mfwwdr[12]+ Pat[pat].bFwwww*pa.bfwwww[12]+ Pat[pat].bFwwwd*pa.bfwwwd[12]+ Pat[pat].bFwwdd*pa.bfwwdd[12]+ Pat[pat].bFwwwr*pa.bfwwwr[12]+ Pat[pat].bFwwrr*pa.bfwwrr[12]+ Pat[pat].bFwwdr*pa.bfwwdr[12] +Pat[pat].fFwdww*pa.ffwdww[12]+ Pat[pat].fFwdwd*pa.ffwdwd[12]+ Pat[pat].fFwddd*pa.ffwddd[12]+ Pat[pat].fFwdwr*pa.ffwdwr[12]+ Pat[pat].fFwdrr*pa.ffwdrr[12]+ Pat[pat].fFwddr*pa.ffwddr[12] +Pat[pat].mFwdww*pa.mfwdww[12]+ Pat[pat].mFwdwd*pa.mfwdwd[12]+ Pat[pat].mFwddd*pa.mfwddd[12]+ Pat[pat].mFwdwr*pa.mfwdwr[12]+ Pat[pat].mFwdrr*pa.mfwdrr[12]+ Pat[pat].mFwddr*pa.mfwddr[12] +Pat[pat].bFwdww*pa.bfwdww[12]+ Pat[pat].bFwdwd*pa.bfwdwd[12]+ Pat[pat].bFwddd*pa.bfwddd[12]+ Pat[pat].bFwdwr*pa.bfwdwr[12]+ Pat[pat].bFwdrr*pa.bfwdrr[12]+ Pat[pat].bFwddr*pa.bfwddr[12] + Pat[pat].Fddww*pa.fddww[12]+ Pat[pat].Fddwd*pa.fddwd[12]+ Pat[pat].Fdddd*pa.fdddd[12]+ Pat[pat].Fddwr*pa.fddwr[12]+ Pat[pat].Fddrr*pa.fddrr[12]+ Pat[pat].Fdddr*pa.fdddr[12] + Pat[pat].Fwrww*pa.fwrww[12]+ Pat[pat].Fwrwd*pa.fwrwd[12]+ Pat[pat].Fwrdd*pa.fwrdd[12]+ Pat[pat].Fwrwr*pa.fwrwr[12]+ Pat[pat].Fwrrr*pa.fwrrr[12]+ Pat[pat].Fwrdr*pa.fwrdr[12] + Pat[pat].fFwrww*pa.ffwrww[12]+ Pat[pat].fFwrwd*pa.ffwrwd[12]+ Pat[pat].fFwrdd*pa.ffwrdd[12]+ Pat[pat].fFwrwr*pa.ffwrwr[12]+ Pat[pat].fFwrrr*pa.ffwrrr[12]+ Pat[pat].fFwrdr*pa.ffwrdr[12] + Pat[pat].mFwrww*pa.mfwrww[12]+ Pat[pat].mFwrwd*pa.mfwrwd[12]+ Pat[pat].mFwrdd*pa.mfwrdd[12]+ Pat[pat].mFwrwr*pa.mfwrwr[12]+ Pat[pat].mFwrrr*pa.mfwrrr[12]+ Pat[pat].mFwrdr*pa.mfwrdr[12] + Pat[pat].bFwrww*pa.bfwrww[12]+ Pat[pat].bFwrwd*pa.bfwrwd[12]+ Pat[pat].bFwrdd*pa.bfwrdd[12]+ Pat[pat].bFwrwr*pa.bfwrwr[12]+ Pat[pat].bFwrrr*pa.bfwrrr[12]+ Pat[pat].bFwrdr*pa.bfwrdr[12] + Pat[pat].Frrww*pa.frrww[12]+ Pat[pat].Frrwd*pa.frrwd[12]+ Pat[pat].Frrdd*pa.frrdd[12]+ Pat[pat].Frrwr*pa.frrwr[12]+ Pat[pat].Frrrr*pa.frrrr[12]+ Pat[pat].Frrdr*pa.frrdr[12] + Pat[pat].Fdrww*pa.fdrww[12]+ Pat[pat].Fdrwd*pa.fdrwd[12]+ Pat[pat].Fdrdd*pa.fdrdd[12]+ Pat[pat].Fdrwr*pa.fdrwr[12]+ Pat[pat].Fdrrr*pa.fdrrr[12]+ Pat[pat].Fdrdr*pa.fdrdr[12]));
	Pat[pat].bJwr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[13]+ Pat[pat].Fwwwd*pa.fwwwd[13]+ Pat[pat].Fwwdd*pa.fwwdd[13]+ Pat[pat].Fwwwr*pa.fwwwr[13]+ Pat[pat].Fwwrr*pa.fwwrr[13]+ Pat[pat].Fwwdr*pa.fwwdr[13]+ Pat[pat].fFwwww*pa.ffwwww[13]+ Pat[pat].fFwwwd*pa.ffwwwd[13]+ Pat[pat].fFwwdd*pa.ffwwdd[13]+ Pat[pat].fFwwwr*pa.ffwwwr[13]+ Pat[pat].fFwwrr*pa.ffwwrr[13]+ Pat[pat].fFwwdr*pa.ffwwdr[13]+ Pat[pat].mFwwww*pa.mfwwww[13]+ Pat[pat].mFwwwd*pa.mfwwwd[13]+ Pat[pat].mFwwdd*pa.mfwwdd[13]+ Pat[pat].mFwwwr*pa.mfwwwr[13]+ Pat[pat].mFwwrr*pa.mfwwrr[13]+ Pat[pat].mFwwdr*pa.mfwwdr[13]+ Pat[pat].bFwwww*pa.bfwwww[13]+ Pat[pat].bFwwwd*pa.bfwwwd[13]+ Pat[pat].bFwwdd*pa.bfwwdd[13]+ Pat[pat].bFwwwr*pa.bfwwwr[13]+ Pat[pat].bFwwrr*pa.bfwwrr[13]+ Pat[pat].bFwwdr*pa.bfwwdr[13] +Pat[pat].fFwdww*pa.ffwdww[13]+ Pat[pat].fFwdwd*pa.ffwdwd[13]+ Pat[pat].fFwddd*pa.ffwddd[13]+ Pat[pat].fFwdwr*pa.ffwdwr[13]+ Pat[pat].fFwdrr*pa.ffwdrr[13]+ Pat[pat].fFwddr*pa.ffwddr[13] +Pat[pat].mFwdww*pa.mfwdww[13]+ Pat[pat].mFwdwd*pa.mfwdwd[13]+ Pat[pat].mFwddd*pa.mfwddd[13]+ Pat[pat].mFwdwr*pa.mfwdwr[13]+ Pat[pat].mFwdrr*pa.mfwdrr[13]+ Pat[pat].mFwddr*pa.mfwddr[13] +Pat[pat].bFwdww*pa.bfwdww[13]+ Pat[pat].bFwdwd*pa.bfwdwd[13]+ Pat[pat].bFwddd*pa.bfwddd[13]+ Pat[pat].bFwdwr*pa.bfwdwr[13]+ Pat[pat].bFwdrr*pa.bfwdrr[13]+ Pat[pat].bFwddr*pa.bfwddr[13] + Pat[pat].Fddww*pa.fddww[13]+ Pat[pat].Fddwd*pa.fddwd[13]+ Pat[pat].Fdddd*pa.fdddd[13]+ Pat[pat].Fddwr*pa.fddwr[13]+ Pat[pat].Fddrr*pa.fddrr[13]+ Pat[pat].Fdddr*pa.fdddr[13] + Pat[pat].Fwrww*pa.fwrww[13]+ Pat[pat].Fwrwd*pa.fwrwd[13]+ Pat[pat].Fwrdd*pa.fwrdd[13]+ Pat[pat].Fwrwr*pa.fwrwr[13]+ Pat[pat].Fwrrr*pa.fwrrr[13]+ Pat[pat].Fwrdr*pa.fwrdr[13] + Pat[pat].fFwrww*pa.ffwrww[13]+ Pat[pat].fFwrwd*pa.ffwrwd[13]+ Pat[pat].fFwrdd*pa.ffwrdd[13]+ Pat[pat].fFwrwr*pa.ffwrwr[13]+ Pat[pat].fFwrrr*pa.ffwrrr[13]+ Pat[pat].fFwrdr*pa.ffwrdr[13] + Pat[pat].mFwrww*pa.mfwrww[13]+ Pat[pat].mFwrwd*pa.mfwrwd[13]+ Pat[pat].mFwrdd*pa.mfwrdd[13]+ Pat[pat].mFwrwr*pa.mfwrwr[13]+ Pat[pat].mFwrrr*pa.mfwrrr[13]+ Pat[pat].mFwrdr*pa.mfwrdr[13] + Pat[pat].bFwrww*pa.bfwrww[13]+ Pat[pat].bFwrwd*pa.bfwrwd[13]+ Pat[pat].bFwrdd*pa.bfwrdd[13]+ Pat[pat].bFwrwr*pa.bfwrwr[13]+ Pat[pat].bFwrrr*pa.bfwrrr[13]+ Pat[pat].bFwrdr*pa.bfwrdr[13] + Pat[pat].Frrww*pa.frrww[13]+ Pat[pat].Frrwd*pa.frrwd[13]+ Pat[pat].Frrdd*pa.frrdd[13]+ Pat[pat].Frrwr*pa.frrwr[13]+ Pat[pat].Frrrr*pa.frrrr[13]+ Pat[pat].Frrdr*pa.frrdr[13] + Pat[pat].Fdrww*pa.fdrww[13]+ Pat[pat].Fdrwd*pa.fdrwd[13]+ Pat[pat].Fdrdd*pa.fdrdd[13]+ Pat[pat].Fdrwr*pa.fdrwr[13]+ Pat[pat].Fdrrr*pa.fdrrr[13]+ Pat[pat].Fdrdr*pa.fdrdr[13]));
	Pat[pat].mbJrr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[14]+ Pat[pat].Fwwwd*pa.fwwwd[14]+ Pat[pat].Fwwdd*pa.fwwdd[14]+ Pat[pat].Fwwwr*pa.fwwwr[14]+ Pat[pat].Fwwrr*pa.fwwrr[14]+ Pat[pat].Fwwdr*pa.fwwdr[14]+ Pat[pat].fFwwww*pa.ffwwww[14]+ Pat[pat].fFwwwd*pa.ffwwwd[14]+ Pat[pat].fFwwdd*pa.ffwwdd[14]+ Pat[pat].fFwwwr*pa.ffwwwr[14]+ Pat[pat].fFwwrr*pa.ffwwrr[14]+ Pat[pat].fFwwdr*pa.ffwwdr[14]+ Pat[pat].mFwwww*pa.mfwwww[14]+ Pat[pat].mFwwwd*pa.mfwwwd[14]+ Pat[pat].mFwwdd*pa.mfwwdd[14]+ Pat[pat].mFwwwr*pa.mfwwwr[14]+ Pat[pat].mFwwrr*pa.mfwwrr[14]+ Pat[pat].mFwwdr*pa.mfwwdr[14]+ Pat[pat].bFwwww*pa.bfwwww[14]+ Pat[pat].bFwwwd*pa.bfwwwd[14]+ Pat[pat].bFwwdd*pa.bfwwdd[14]+ Pat[pat].bFwwwr*pa.bfwwwr[14]+ Pat[pat].bFwwrr*pa.bfwwrr[14]+ Pat[pat].bFwwdr*pa.bfwwdr[14] +Pat[pat].fFwdww*pa.ffwdww[14]+ Pat[pat].fFwdwd*pa.ffwdwd[14]+ Pat[pat].fFwddd*pa.ffwddd[14]+ Pat[pat].fFwdwr*pa.ffwdwr[14]+ Pat[pat].fFwdrr*pa.ffwdrr[14]+ Pat[pat].fFwddr*pa.ffwddr[14] +Pat[pat].mFwdww*pa.mfwdww[14]+ Pat[pat].mFwdwd*pa.mfwdwd[14]+ Pat[pat].mFwddd*pa.mfwddd[14]+ Pat[pat].mFwdwr*pa.mfwdwr[14]+ Pat[pat].mFwdrr*pa.mfwdrr[14]+ Pat[pat].mFwddr*pa.mfwddr[14] +Pat[pat].bFwdww*pa.bfwdww[14]+ Pat[pat].bFwdwd*pa.bfwdwd[14]+ Pat[pat].bFwddd*pa.bfwddd[14]+ Pat[pat].bFwdwr*pa.bfwdwr[14]+ Pat[pat].bFwdrr*pa.bfwdrr[14]+ Pat[pat].bFwddr*pa.bfwddr[14] + Pat[pat].Fddww*pa.fddww[14]+ Pat[pat].Fddwd*pa.fddwd[14]+ Pat[pat].Fdddd*pa.fdddd[14]+ Pat[pat].Fddwr*pa.fddwr[14]+ Pat[pat].Fddrr*pa.fddrr[14]+ Pat[pat].Fdddr*pa.fdddr[14] + Pat[pat].Fwrww*pa.fwrww[14]+ Pat[pat].Fwrwd*pa.fwrwd[14]+ Pat[pat].Fwrdd*pa.fwrdd[14]+ Pat[pat].Fwrwr*pa.fwrwr[14]+ Pat[pat].Fwrrr*pa.fwrrr[14]+ Pat[pat].Fwrdr*pa.fwrdr[14] + Pat[pat].fFwrww*pa.ffwrww[14]+ Pat[pat].fFwrwd*pa.ffwrwd[14]+ Pat[pat].fFwrdd*pa.ffwrdd[14]+ Pat[pat].fFwrwr*pa.ffwrwr[14]+ Pat[pat].fFwrrr*pa.ffwrrr[14]+ Pat[pat].fFwrdr*pa.ffwrdr[14] + Pat[pat].mFwrww*pa.mfwrww[14]+ Pat[pat].mFwrwd*pa.mfwrwd[14]+ Pat[pat].mFwrdd*pa.mfwrdd[14]+ Pat[pat].mFwrwr*pa.mfwrwr[14]+ Pat[pat].mFwrrr*pa.mfwrrr[14]+ Pat[pat].mFwrdr*pa.mfwrdr[14] + Pat[pat].bFwrww*pa.bfwrww[14]+ Pat[pat].bFwrwd*pa.bfwrwd[14]+ Pat[pat].bFwrdd*pa.bfwrdd[14]+ Pat[pat].bFwrwr*pa.bfwrwr[14]+ Pat[pat].bFwrrr*pa.bfwrrr[14]+ Pat[pat].bFwrdr*pa.bfwrdr[14] + Pat[pat].Frrww*pa.frrww[14]+ Pat[pat].Frrwd*pa.frrwd[14]+ Pat[pat].Frrdd*pa.frrdd[14]+ Pat[pat].Frrwr*pa.frrwr[14]+ Pat[pat].Frrrr*pa.frrrr[14]+ Pat[pat].Frrdr*pa.frrdr[14] + Pat[pat].Fdrww*pa.fdrww[14]+ Pat[pat].Fdrwd*pa.fdrwd[14]+ Pat[pat].Fdrdd*pa.fdrdd[14]+ Pat[pat].Fdrwr*pa.fdrwr[14]+ Pat[pat].Fdrrr*pa.fdrrr[14]+ Pat[pat].Fdrdr*pa.fdrdr[14]));
	Pat[pat].mbJdr[0]=random_poisson(pa.theta*( Pat[pat].Fwwww*pa.fwwww[15]+ Pat[pat].Fwwwd*pa.fwwwd[15]+ Pat[pat].Fwwdd*pa.fwwdd[15]+ Pat[pat].Fwwwr*pa.fwwwr[15]+ Pat[pat].Fwwrr*pa.fwwrr[15]+ Pat[pat].Fwwdr*pa.fwwdr[15]+ Pat[pat].fFwwww*pa.ffwwww[15]+ Pat[pat].fFwwwd*pa.ffwwwd[15]+ Pat[pat].fFwwdd*pa.ffwwdd[15]+ Pat[pat].fFwwwr*pa.ffwwwr[15]+ Pat[pat].fFwwrr*pa.ffwwrr[15]+ Pat[pat].fFwwdr*pa.ffwwdr[15]+ Pat[pat].mFwwww*pa.mfwwww[15]+ Pat[pat].mFwwwd*pa.mfwwwd[15]+ Pat[pat].mFwwdd*pa.mfwwdd[15]+ Pat[pat].mFwwwr*pa.mfwwwr[15]+ Pat[pat].mFwwrr*pa.mfwwrr[15]+ Pat[pat].mFwwdr*pa.mfwwdr[15]+ Pat[pat].bFwwww*pa.bfwwww[15]+ Pat[pat].bFwwwd*pa.bfwwwd[15]+ Pat[pat].bFwwdd*pa.bfwwdd[15]+ Pat[pat].bFwwwr*pa.bfwwwr[15]+ Pat[pat].bFwwrr*pa.bfwwrr[15]+ Pat[pat].bFwwdr*pa.bfwwdr[15] +Pat[pat].fFwdww*pa.ffwdww[15]+ Pat[pat].fFwdwd*pa.ffwdwd[15]+ Pat[pat].fFwddd*pa.ffwddd[15]+ Pat[pat].fFwdwr*pa.ffwdwr[15]+ Pat[pat].fFwdrr*pa.ffwdrr[15]+ Pat[pat].fFwddr*pa.ffwddr[15] +Pat[pat].mFwdww*pa.mfwdww[15]+ Pat[pat].mFwdwd*pa.mfwdwd[15]+ Pat[pat].mFwddd*pa.mfwddd[15]+ Pat[pat].mFwdwr*pa.mfwdwr[15]+ Pat[pat].mFwdrr*pa.mfwdrr[15]+ Pat[pat].mFwddr*pa.mfwddr[15] +Pat[pat].bFwdww*pa.bfwdww[15]+ Pat[pat].bFwdwd*pa.bfwdwd[15]+ Pat[pat].bFwddd*pa.bfwddd[15]+ Pat[pat].bFwdwr*pa.bfwdwr[15]+ Pat[pat].bFwdrr*pa.bfwdrr[15]+ Pat[pat].bFwddr*pa.bfwddr[15] + Pat[pat].Fddww*pa.fddww[15]+ Pat[pat].Fddwd*pa.fddwd[15]+ Pat[pat].Fdddd*pa.fdddd[15]+ Pat[pat].Fddwr*pa.fddwr[15]+ Pat[pat].Fddrr*pa.fddrr[15]+ Pat[pat].Fdddr*pa.fdddr[15] + Pat[pat].Fwrww*pa.fwrww[15]+ Pat[pat].Fwrwd*pa.fwrwd[15]+ Pat[pat].Fwrdd*pa.fwrdd[15]+ Pat[pat].Fwrwr*pa.fwrwr[15]+ Pat[pat].Fwrrr*pa.fwrrr[15]+ Pat[pat].Fwrdr*pa.fwrdr[15] + Pat[pat].fFwrww*pa.ffwrww[15]+ Pat[pat].fFwrwd*pa.ffwrwd[15]+ Pat[pat].fFwrdd*pa.ffwrdd[15]+ Pat[pat].fFwrwr*pa.ffwrwr[15]+ Pat[pat].fFwrrr*pa.ffwrrr[15]+ Pat[pat].fFwrdr*pa.ffwrdr[15] + Pat[pat].mFwrww*pa.mfwrww[15]+ Pat[pat].mFwrwd*pa.mfwrwd[15]+ Pat[pat].mFwrdd*pa.mfwrdd[15]+ Pat[pat].mFwrwr*pa.mfwrwr[15]+ Pat[pat].mFwrrr*pa.mfwrrr[15]+ Pat[pat].mFwrdr*pa.mfwrdr[15] + Pat[pat].bFwrww*pa.bfwrww[15]+ Pat[pat].bFwrwd*pa.bfwrwd[15]+ Pat[pat].bFwrdd*pa.bfwrdd[15]+ Pat[pat].bFwrwr*pa.bfwrwr[15]+ Pat[pat].bFwrrr*pa.bfwrrr[15]+ Pat[pat].bFwrdr*pa.bfwrdr[15] + Pat[pat].Frrww*pa.frrww[15]+ Pat[pat].Frrwd*pa.frrwd[15]+ Pat[pat].Frrdd*pa.frrdd[15]+ Pat[pat].Frrwr*pa.frrwr[15]+ Pat[pat].Frrrr*pa.frrrr[15]+ Pat[pat].Frrdr*pa.frrdr[15] + Pat[pat].Fdrww*pa.fdrww[15]+ Pat[pat].Fdrwd*pa.fdrwd[15]+ Pat[pat].Fdrdd*pa.fdrdd[15]+ Pat[pat].Fdrwr*pa.fdrwr[15]+ Pat[pat].Fdrrr*pa.fdrrr[15]+ Pat[pat].Fdrdr*pa.fdrdr[15]));



		to.Jww+=Pat[pat].Jww[0]+Pat[pat].fJww[0]+Pat[pat].mJww[0]+Pat[pat].bJww[0];
		to.Jwd+=Pat[pat].fJwd[0]+Pat[pat].mJwd[0]+Pat[pat].bJwd[0];
		to.Jdd+=Pat[pat].Jdd[0];
		to.Jwr+=Pat[pat].Jwr[0]+Pat[pat].fJwr[0]+Pat[pat].mJwr[0]+Pat[pat].bJwr[0];
		to.Jrr+=Pat[pat].Jrr[0]+Pat[pat].mbJrr[0];
		to.Jdr+=Pat[pat].Jdr[0]+Pat[pat].mbJdr[0];
		to.JTot+=Pat[pat].Jww[0]+Pat[pat].fJww[0]+Pat[pat].mJww[0]+Pat[pat].bJww[0]+ Pat[pat].fJwd[0]+Pat[pat].mJwd[0]+Pat[pat].bJwd[0]+Pat[pat].Jdd[0]+Pat[pat].Jwr[0]+Pat[pat].fJwr[0]+Pat[pat].mJwr[0]+Pat[pat].bJwr[0]+Pat[pat].Jrr[0]+Pat[pat].Jdr[0]+Pat[pat].mbJrr[0]+Pat[pat].mbJdr[0];
				};
				return;};


	void AdultsDie(void){
		int num;
		for(int pat=0;pat<Pat.size();pat++)
				{
num=random_binomial(Pat[pat].Mww,pa.muA);Pat[pat].Mww-=num;Pat[pat].MTot-=num;to.Mww-=num;to.MTot-=num;	
num=random_binomial(Pat[pat].Mwd,pa.muA);Pat[pat].Mwd-=num;Pat[pat].MTot-=num;to.Mwd-=num;to.MTot-=num;	
num=random_binomial(Pat[pat].Mdd,pa.muA);Pat[pat].Mdd-=num;Pat[pat].MTot-=num;to.Mdd-=num;to.MTot-=num;	
num=random_binomial(Pat[pat].Mwr,pa.muA);Pat[pat].Mwr-=num;Pat[pat].MTot-=num;to.Mwr-=num;to.MTot-=num;	
num=random_binomial(Pat[pat].Mrr,pa.muA);Pat[pat].Mrr-=num;Pat[pat].MTot-=num;to.Mrr-=num;to.MTot-=num;	
num=random_binomial(Pat[pat].Mdr,pa.muA);Pat[pat].Mdr-=num;Pat[pat].MTot-=num;to.Mdr-=num;to.MTot-=num;	

num=random_binomial(Pat[pat].Vww,pa.muA);Pat[pat].Vww-=num;to.Vww-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].fVww,pa.muA);Pat[pat].fVww-=num;to.Vww-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].mVww,pa.muA);Pat[pat].mVww-=num;to.Vww-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].bVww,pa.muA);Pat[pat].bVww-=num;to.Vww-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].fVwd,pa.muA);Pat[pat].fVwd-=num;to.Vwd-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].mVwd,pa.muA);Pat[pat].mVwd-=num;to.Vwd-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].bVwd,pa.muA);Pat[pat].bVwd-=num;to.Vwd-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].Vdd,pa.muA);Pat[pat].Vdd-=num;to.Vdd-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].Vwr,pa.muA);Pat[pat].Vwr-=num;to.Vwr-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].fVwr,pa.muA);Pat[pat].fVwr-=num;to.Vwr-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].mVwr,pa.muA);Pat[pat].mVwr-=num;to.Vwr-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].bVwr,pa.muA);Pat[pat].bVwr-=num;to.Vwr-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].Vrr,pa.muA);Pat[pat].Vrr-=num;to.Vrr-=num;to.VTot-=num;	
num=random_binomial(Pat[pat].Vdr,pa.muA);Pat[pat].Vdr-=num;to.Vdr-=num;to.VTot-=num;	

num=random_binomial(Pat[pat].Fwwww,pa.muA);Pat[pat].Fwwww-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwwwd,pa.muA);Pat[pat].Fwwwd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwwdd,pa.muA);Pat[pat].Fwwdd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwwwr,pa.muA);Pat[pat].Fwwwr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwwrr,pa.muA);Pat[pat].Fwwrr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwwdr,pa.muA);Pat[pat].Fwwdr-=num;to.Fww-=num;to.FTot-=num;	

num=random_binomial(Pat[pat].fFwwww,pa.muA);Pat[pat].fFwwww-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwwwd,pa.muA);Pat[pat].fFwwwd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwwdd,pa.muA);Pat[pat].fFwwdd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwwwr,pa.muA);Pat[pat].fFwwwr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwwrr,pa.muA);Pat[pat].fFwwrr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwwdr,pa.muA);Pat[pat].fFwwdr-=num;to.Fww-=num;to.FTot-=num;	

num=random_binomial(Pat[pat].mFwwww,pa.muA);Pat[pat].mFwwww-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwwwd,pa.muA);Pat[pat].mFwwwd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwwdd,pa.muA);Pat[pat].mFwwdd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwwwr,pa.muA);Pat[pat].mFwwwr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwwrr,pa.muA);Pat[pat].mFwwrr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwwdr,pa.muA);Pat[pat].mFwwdr-=num;to.Fww-=num;to.FTot-=num;	

num=random_binomial(Pat[pat].bFwwww,pa.muA);Pat[pat].bFwwww-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwwwd,pa.muA);Pat[pat].bFwwwd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwwdd,pa.muA);Pat[pat].bFwwdd-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwwwr,pa.muA);Pat[pat].bFwwwr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwwrr,pa.muA);Pat[pat].bFwwrr-=num;to.Fww-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwwdr,pa.muA);Pat[pat].bFwwdr-=num;to.Fww-=num;to.FTot-=num;	


num=random_binomial(Pat[pat].fFwdww,pa.muA);Pat[pat].fFwdww-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwdwd,pa.muA);Pat[pat].fFwdwd-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwddd,pa.muA);Pat[pat].fFwddd-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwdwr,pa.muA);Pat[pat].fFwdwr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwdrr,pa.muA);Pat[pat].fFwdrr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwddr,pa.muA);Pat[pat].fFwddr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwdww,pa.muA);Pat[pat].mFwdww-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwdwd,pa.muA);Pat[pat].mFwdwd-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwddd,pa.muA);Pat[pat].mFwddd-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwdwr,pa.muA);Pat[pat].mFwdwr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwdrr,pa.muA);Pat[pat].mFwdrr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwddr,pa.muA);Pat[pat].mFwddr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwdww,pa.muA);Pat[pat].bFwdww-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwdwd,pa.muA);Pat[pat].bFwdwd-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwddd,pa.muA);Pat[pat].bFwddd-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwdwr,pa.muA);Pat[pat].bFwdwr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwdrr,pa.muA);Pat[pat].bFwdrr-=num;to.Fwd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwddr,pa.muA);Pat[pat].bFwddr-=num;to.Fwd-=num;to.FTot-=num;	

num=random_binomial(Pat[pat].Fddww,pa.muA);Pat[pat].Fddww-=num;to.Fdd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fddwd,pa.muA);Pat[pat].Fddwd-=num;to.Fdd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdddd,pa.muA);Pat[pat].Fdddd-=num;to.Fdd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fddwr,pa.muA);Pat[pat].Fddwr-=num;to.Fdd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fddrr,pa.muA);Pat[pat].Fddrr-=num;to.Fdd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdddr,pa.muA);Pat[pat].Fdddr-=num;to.Fdd-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwrww,pa.muA);Pat[pat].Fwrww-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwrwd,pa.muA);Pat[pat].Fwrwd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwrdd,pa.muA);Pat[pat].Fwrdd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwrwr,pa.muA);Pat[pat].Fwrwr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwrrr,pa.muA);Pat[pat].Fwrrr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fwrdr,pa.muA);Pat[pat].Fwrdr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwrww,pa.muA);Pat[pat].fFwrww-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwrwd,pa.muA);Pat[pat].fFwrwd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwrdd,pa.muA);Pat[pat].fFwrdd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwrwr,pa.muA);Pat[pat].fFwrwr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwrrr,pa.muA);Pat[pat].fFwrrr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].fFwrdr,pa.muA);Pat[pat].fFwrdr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwrww,pa.muA);Pat[pat].mFwrww-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwrwd,pa.muA);Pat[pat].mFwrwd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwrdd,pa.muA);Pat[pat].mFwrdd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwrwr,pa.muA);Pat[pat].mFwrwr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwrrr,pa.muA);Pat[pat].mFwrrr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].mFwrdr,pa.muA);Pat[pat].mFwrdr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwrww,pa.muA);Pat[pat].bFwrww-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwrwd,pa.muA);Pat[pat].bFwrwd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwrdd,pa.muA);Pat[pat].bFwrdd-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwrwr,pa.muA);Pat[pat].bFwrwr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwrrr,pa.muA);Pat[pat].bFwrrr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].bFwrdr,pa.muA);Pat[pat].bFwrdr-=num;to.Fwr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Frrww,pa.muA);Pat[pat].Frrww-=num;to.Frr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Frrwd,pa.muA);Pat[pat].Frrwd-=num;to.Frr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Frrdd,pa.muA);Pat[pat].Frrdd-=num;to.Frr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Frrwr,pa.muA);Pat[pat].Frrwr-=num;to.Frr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Frrrr,pa.muA);Pat[pat].Frrrr-=num;to.Frr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Frrdr,pa.muA);Pat[pat].Frrdr-=num;to.Frr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdrww,pa.muA);Pat[pat].Fdrww-=num;to.Fdr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdrwd,pa.muA);Pat[pat].Fdrwd-=num;to.Fdr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdrdd,pa.muA);Pat[pat].Fdrdd-=num;to.Fdr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdrwr,pa.muA);Pat[pat].Fdrwr-=num;to.Fdr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdrrr,pa.muA);Pat[pat].Fdrrr-=num;to.Fdr-=num;to.FTot-=num;	
num=random_binomial(Pat[pat].Fdrdr,pa.muA);Pat[pat].Fdrdr-=num;to.Fdr-=num;to.FTot-=num;	
				};
		return;};

	void UpdateComp(int day){
		int week=(int)(day % 365)/7;
		if(week==52)week=51;
		double al;
		for(int pat=0;pat<Pat.size();pat++)
				{
	al=Pat[pat].alpha0+pa.alpha1*(1-exp(-1*pa.phi*rain[Pat[pat].sqx][Pat[pat].sqy][ti.yearnow][week])) +pa.alpha2*(1-exp(-1*pa.kappa*(Pat[pat].WaterPerm +(1-exp(-1*pa.delta*rain[Pat[pat].sqx][Pat[pat].sqy][ti.yearnow][week]))*Pat[pat].WaterTemp)));
				if(day>0){Pat[pat].comp=std::pow(al/(al+Pat[pat].JTot+0.0001),1/double(TL));}
				else {Pat[pat].comp=std::pow(pa.alpha1/(pa.alpha1+Pat[pat].JTot),1/double(TL));};
				if(!(Pat[pat].comp+1>0)){cout<<"uc  "<<Pat[pat].comp<<"     "<<pat<<"   rain  "<<rain[Pat[pat].sqx][Pat[pat].sqy][ti.yearnow][week]<<endl;exit(0);};
				};
		return;};


	void UpdateMate(void){
		for(int pat=0;pat<Pat.size();pat++)
				{
				Pat[pat].mate_rate=Pat[pat].MTot/(pa.beta+Pat[pat].MTot);
				};
		return;};

	void SetFertility()
{
//fraction of 1:ww,2:fww,3:mww,4:fwd,5:mwd,6:dd,7:wr,8:fwr,9:mwr,10:fnrr,11:fndr,12:bww,13:bwd,14:bwr,15:mbrr,16:mbdr
	double Fwwww[16]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double Fwwwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0,0,0};
	double Fwwdd[16]={0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
	double Fwwwr[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double Fwwrr[16]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double Fwwdr[16]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0};

	double fFwwww[16]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double fFwwwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0,0,0};
	double fFwwdd[16]={0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
	double fFwwwr[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double fFwwrr[16]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double fFwwdr[16]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0};
	double mFwwww[16]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double mFwwwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0,0,0};
	double mFwwdd[16]={0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
	double mFwwwr[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double mFwwrr[16]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double mFwwdr[16]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0};

	double bFwwww[16]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double bFwwwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0,0,0};
	double bFwwdd[16]={0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
	double bFwwwr[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double bFwwrr[16]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double bFwwdr[16]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0};


	double fFwdww[16]={0,(1-pa.ef-pa.Fgamma)*0.5,0,(1+pa.ef)*0.5,0,0,0,pa.Fgamma*0.5,0,0,0,0,0,0,0,0};
	double fFwdwd[16]={0,0,0,0,0, (1+pa.ef)*(1+pa.em)*0.25,0,0,0,0,0 ,(1-pa.ef-pa.Fgamma)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*(1+pa.em)*0.25+ (1+pa.ef)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*pa.Mgamma*0.25+ (1-pa.em-pa.Mgamma)*pa.Fgamma*0.25,pa.Mgamma*pa.Fgamma*0.25,(1+pa.ef)*pa.Mgamma*0.25+(1+pa.em)*pa.Fgamma*0.25};

	double fFwddd[16]={0,0,0,0,0,(1+pa.ef)*0.5,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,0,0,pa.Fgamma*0.5};



	double fFwdwr[16]={0,(1-pa.ef-pa.Fgamma)*0.25,0,(1+pa.ef)*0.25,0,0,0,(1-pa.ef-pa.Fgamma)*0.25+pa.Fgamma*0.25,0,pa.Fgamma*0.25,(1+pa.ef)*0.25,0,0,0,0,0};
	double fFwdrr[16]={0,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,0,pa.Fgamma*0.5,(1+pa.ef)*0.5,0,0,0,0,0};

	double fFwddr[16]={0,0,0,0,0,(1+pa.ef)*0.25,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.25,(1-pa.ef-pa.Fgamma)*0.25,pa.Fgamma*0.25,((1+pa.ef)*0.25+pa.Fgamma*0.25)};

	double mFwdww[16]={0,(1-pa.ef-pa.Fgamma)*0.5,0,(1+pa.ef)*0.5,0,0,0,pa.Fgamma*0.5,0,0,0,0,0,0,0,0};

	double mFwdwd[16]={0,0,0,0,0,(1+pa.ef)*(1+pa.em)*0.25,0,0,0 ,0,0 ,(1-pa.ef-pa.Fgamma)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*(1+pa.em)*0.25+ (1+pa.ef)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*pa.Mgamma*0.25+ (1-pa.em-pa.Mgamma)*pa.Fgamma*0.25,pa.Mgamma*pa.Fgamma*0.25,(1+pa.em)*pa.Fgamma*0.25+(1+pa.ef)*pa.Mgamma*0.25};

	double mFwddd[16]={0,0,0,0,0,(1+pa.ef)*0.5,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,0,0,pa.Fgamma*0.5};

	double mFwdwr[16]={0,(1-pa.ef-pa.Fgamma)*0.25,0,(1+pa.ef)*0.25,0,0,0,((1-pa.ef-pa.Fgamma)*0.25+pa.Fgamma*0.25),0,pa.Fgamma*0.25,(1+pa.ef)*0.25,0,0,0,0,0};
	double mFwdrr[16]={0,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,pa.Fgamma*0.5,0,(1+pa.ef)*0.5,0,0,0,0,0};
	double mFwddr[16]={0,0,0,0,0,(1+pa.ef)*0.25,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.25,(1-pa.ef-pa.Fgamma)*0.25,pa.Fgamma*0.25,((1+pa.ef)*0.25+pa.Fgamma*0.25)};
	double bFwdww[16]={0,(1-pa.ef-pa.Fgamma)*0.5,0,(1+pa.ef)*0.5,0,0,0,pa.Fgamma*0.5,0,0,0,0,0,0,0,0};
	double bFwdwd[16]={0,0,0,0,0,(1+pa.ef)*(1+pa.em)*0.25,0,0,0 ,0,0,(1-pa.ef-pa.Fgamma)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*(1+pa.em)*0.25+ (1+pa.ef)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*pa.Mgamma*0.25+ (1-pa.em-pa.Mgamma)*pa.Fgamma*0.25, pa.Mgamma*pa.Fgamma*0.25,(1+pa.em)*pa.Fgamma*0.25+(1+pa.ef)*pa.Mgamma*0.25};
	double bFwddd[16]={0,0,0,0,0,(1+pa.ef)*0.5,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,0,0,pa.Fgamma*0.5};
	double bFwdwr[16]={0,(1-pa.ef-pa.Fgamma)*0.25,0,(1+pa.ef)*0.25,0,0,0,((1-pa.ef-pa.Fgamma)*0.25+pa.Fgamma*0.25),0,pa.Fgamma*0.25,(1+pa.ef)*0.25,0,0,0,0,0};
	double bFwdrr[16]={0,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,pa.Fgamma*0.5,0,(1+pa.ef)*0.5,0,0,0,0,0};
	double bFwddr[16]={0,0,0,0,0,(1+pa.ef)*0.25,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.25,(1-pa.ef-pa.Fgamma)*0.25,pa.Fgamma*0.25,((1+pa.ef)*0.25+pa.Fgamma*0.25)};
	double Fddww[16]={0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
	double Fddwd[16]={0,0,0,0,0,(1+pa.em)*0.5,0,0,0,0,0,0,(1-pa.em-pa.Mgamma)*0.5,0,0,pa.Mgamma*0.5};
	double Fdddd[16]={0,0,0,0,0,1,0,0,0,0,0,0,0,0};
	double Fddwr[16]={0,0,0,0.5,0,0,0,0,0,0,0.5,0,0,0,0,0};
	double Fddrr[16]={0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0};
	double Fdddr[16]={0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0.5};
	double Fwrww[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double Fwrwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),0,0,0,0,0,pa.Mgamma*0.25,(1+pa.em)*0.25};
	double Fwrdd[16]={0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0.5};
	double Fwrwr[16]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0,0,0};
	double Fwrrr[16]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0};
	double Fwrdr[16]={0,0,0,0,0.25,0,0,0,0.25,0,0,0,0,0,0.25,0.25};
	double fFwrww[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double fFwrwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),0,0,0,0,0,pa.Mgamma*0.25,(1+pa.em)*0.25};
	double fFwrdd[16]={0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0.5};
	double fFwrwr[16]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0,0,0};
	double fFwrrr[16]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0};
	double fFwrdr[16]={0,0,0,0,0.25,0,0,0,0.25,0,0,0,0,0,0.25,0.25};
	double mFwrww[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double mFwrwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),0,0,0,0,0,pa.Mgamma*0.25,(1+pa.em)*0.25};
	double mFwrdd[16]={0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0.5};
	double mFwrwr[16]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0,0,0};
	double mFwrrr[16]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0};
	double mFwrdr[16]={0,0,0,0,0.25,0,0,0,0.25,0,0,0,0,0,0.25,0.25};
	double bFwrww[16]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0};
	double bFwrwd[16]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),0,0,0,0,0,pa.Mgamma*0.25,(1+pa.em)*0.25};
	double bFwrdd[16]={0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0.5};
	double bFwrwr[16]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0,0,0};
	double bFwrrr[16]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0};
	double bFwrdr[16]={0,0,0,0,0.25,0,0,0,0.25,0,0,0,0,0,0.25,0.25};
	double Frrww[16]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double Frrwd[16]={0,0,0,0,0,0,0,0,(1-pa.em-pa.Mgamma)*0.5,0,0,0,0,0,0.5*pa.Mgamma,(1+pa.em)*0.5};
	double Frrdd[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	double Frrwr[16]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0};
	double Frrrr[16]={0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
	double Frrdr[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0.5};
	double Fdrww[16]={0,0,0,0.5,0,0,0,0.5,0,0,0,0,0,0,0,0};
	double Fdrwd[16]={0,0,0,0,(1-pa.em-pa.Mgamma)*0.25,(1+pa.em)*0.25,0,0,(1-pa.em-pa.Mgamma)*0.25,0,0,0,0,0,pa.Mgamma*0.25,((1+pa.em)*0.25+pa.Mgamma*0.25)};
	double Fdrdd[16]={0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0.5};
	double Fdrwr[16]={0,0,0,0.25,0,0,0,0.25,0,0.25,0.25,0,0,0,0,0};
	double Fdrrr[16]={0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0,0,0};
	double Fdrdr[16]={0,0,0,0,0,0.25,0,0,0,0,0,0,0,0,0.25,0.5};


	double omegafww=(1-pa.Frho); double omegamww=(1-pa.Mrho); double omegabww=(1-pa.Frho)*(1-pa.Mrho);
	double omegafwd=(1-pa.Frho)*(1-pa.xi); double omegamwd=(1-pa.Mrho)*(1-pa.xi); double omegabwd=(1-pa.Frho)*(1-pa.Mrho)*(1-pa.xi);
	double omegafwr=(1-pa.Frho); double omegamwr=(1-pa.Mrho); double omegabwr=(1-pa.Frho)*(1-pa.Mrho);
	double omegadd=0; double omegadr=0; double omegarr=0; double omegawr=1;


	for(int i=0;i<16;i++)
	{
	fFwwww[i]*=omegafww; fFwwwd[i]*=omegafww; fFwwdd[i]*=omegafww; fFwwwr[i]*=omegafww; fFwwrr[i]*=omegafww; fFwwdr[i]*=omegafww;
	mFwwww[i]*=omegamww; mFwwwd[i]*=omegamww; mFwwdd[i]*=omegamww; mFwwwr[i]*=omegamww; mFwwrr[i]*=omegamww; mFwwdr[i]*=omegamww;
	bFwwww[i]*=omegabww; bFwwwd[i]*=omegabww; bFwwdd[i]*=omegabww; bFwwwr[i]*=omegabww; bFwwrr[i]*=omegabww; bFwwdr[i]*=omegabww;
	fFwdww[i]*=omegafwd; fFwdwd[i]*=omegafwd; fFwddd[i]*=omegafwd; fFwdwr[i]*=omegafwd; fFwdrr[i]*=omegafwd; fFwddr[i]*=omegafwd;
	mFwdww[i]*=omegamwd; mFwdwd[i]*=omegamwd; mFwddd[i]*=omegamwd; mFwdwr[i]*=omegamwd; mFwdrr[i]*=omegamwd; mFwddr[i]*=omegamwd;
	bFwdww[i]*=omegabwd; bFwdwd[i]*=omegabwd; bFwddd[i]*=omegabwd; bFwdwr[i]*=omegabwd; bFwdrr[i]*=omegabwd; bFwddr[i]*=omegabwd;
	Fddww[i]*=omegadd; Fddwd[i]*=omegadd; Fdddd[i]*=omegadd; Fddwr[i]*=omegadd; Fddrr[i]*=omegadd; Fdddr[i]*=omegadd;
	Fwrww[i]*=omegawr; Fwrwd[i]*=omegawr; Fwrdd[i]*=omegawr; Fwrwr[i]*=omegawr; Fwrrr[i]*=omegawr; Fwrdr[i]*=omegawr;
	fFwrww[i]*=omegafwr; fFwrwd[i]*=omegafwr; fFwrdd[i]*=omegafwr; fFwrwr[i]*=omegafwr; fFwrrr[i]*=omegafwr; fFwrdr[i]*=omegafwr;
	mFwrww[i]*=omegamwr; mFwrwd[i]*=omegamwr; mFwrdd[i]*=omegamwr; mFwrwr[i]*=omegamwr; mFwrrr[i]*=omegamwr; mFwrdr[i]*=omegamwr;
	bFwrww[i]*=omegabwr; bFwrwd[i]*=omegabwr; bFwrdd[i]*=omegabwr; bFwrwr[i]*=omegabwr; bFwrrr[i]*=omegabwr; bFwrdr[i]*=omegabwr;
	Frrww[i]*=omegarr; Frrwd[i]*=omegarr; Frrdd[i]*=omegarr; Frrwr[i]*=omegarr; Frrrr[i]*=omegarr; Frrdr[i]*=omegarr;
	Fdrww[i]*=omegadr; Fdrwd[i]*=omegadr; Fdrdd[i]*=omegadr; Fdrwr[i]*=omegadr; Fdrrr[i]*=omegadr; Fdrdr[i]*=omegadr;
	};




	for(int k=0;k<16;k++)
	{
	pa.fwwww[k]=Fwwww[k]; pa.fwwwd[k]=Fwwwd[k]; pa.fwwdd[k]=Fwwdd[k]; pa.fwwwr[k]=Fwwwr[k]; pa.fwwrr[k]=Fwwrr[k]; pa.fwwdr[k]=Fwwdr[k]; 
	pa.ffwwww[k]=fFwwww[k]; pa.ffwwwd[k]=fFwwwd[k]; pa.ffwwdd[k]=fFwwdd[k]; pa.ffwwwr[k]=fFwwwr[k]; pa.ffwwrr[k]=fFwwrr[k]; pa.ffwwdr[k]=fFwwdr[k]; 
	pa.mfwwww[k]=mFwwww[k]; pa.mfwwwd[k]=mFwwwd[k]; pa.mfwwdd[k]=mFwwdd[k]; pa.mfwwwr[k]=mFwwwr[k]; pa.mfwwrr[k]=mFwwrr[k]; pa.mfwwdr[k]=mFwwdr[k]; 
	pa.bfwwww[k]=bFwwww[k]; pa.bfwwwd[k]=bFwwwd[k]; pa.bfwwdd[k]=bFwwdd[k]; pa.bfwwwr[k]=bFwwwr[k]; pa.bfwwrr[k]=bFwwrr[k]; pa.bfwwdr[k]=bFwwdr[k]; 
	pa.ffwdww[k]=fFwdww[k]; pa.ffwdwd[k]=fFwdwd[k]; pa.ffwddd[k]=fFwddd[k]; pa.ffwdwr[k]=fFwdwr[k]; pa.ffwdrr[k]=fFwdrr[k]; pa.ffwddr[k]=fFwddr[k]; 
	pa.mfwdww[k]=mFwdww[k]; pa.mfwdwd[k]=mFwdwd[k]; pa.mfwddd[k]=mFwddd[k]; pa.mfwdwr[k]=mFwdwr[k]; pa.mfwdrr[k]=mFwdrr[k]; pa.mfwddr[k]=mFwddr[k]; 
	pa.bfwdww[k]=bFwdww[k]; pa.bfwdwd[k]=bFwdwd[k]; pa.bfwddd[k]=bFwddd[k]; pa.bfwdwr[k]=bFwdwr[k]; pa.bfwdrr[k]=bFwdrr[k]; pa.bfwddr[k]=bFwddr[k]; 
	pa.fddww[k]=Fddww[k]; pa.fddwd[k]=Fddwd[k]; pa.fdddd[k]=Fdddd[k]; pa.fddwr[k]=Fddwr[k]; pa.fddrr[k]=Fddrr[k]; pa.fdddr[k]=Fdddr[k]; 
	pa.fwrww[k]=Fwrww[k]; pa.fwrwd[k]=Fwrwd[k]; pa.fwrdd[k]=Fwrdd[k]; pa.fwrwr[k]=Fwrwr[k]; pa.fwrrr[k]=Fwrrr[k]; pa.fwrdr[k]=Fwrdr[k]; 
	pa.ffwrww[k]=fFwrww[k]; pa.ffwrwd[k]=fFwrwd[k]; pa.ffwrdd[k]=fFwrdd[k]; pa.ffwrwr[k]=fFwrwr[k]; pa.ffwrrr[k]=fFwrrr[k]; pa.ffwrdr[k]=fFwrdr[k]; 
	pa.mfwrww[k]=mFwrww[k]; pa.mfwrwd[k]=mFwrwd[k]; pa.mfwrdd[k]=mFwrdd[k]; pa.mfwrwr[k]=mFwrwr[k]; pa.mfwrrr[k]=mFwrrr[k]; pa.mfwrdr[k]=mFwrdr[k]; 
	pa.bfwrww[k]=bFwrww[k]; pa.bfwrwd[k]=bFwrwd[k]; pa.bfwrdd[k]=bFwrdd[k]; pa.bfwrwr[k]=bFwrwr[k]; pa.bfwrrr[k]=bFwrrr[k]; pa.bfwrdr[k]=bFwrdr[k]; 
	pa.frrww[k]=Frrww[k]; pa.frrwd[k]=Frrwd[k]; pa.frrdd[k]=Frrdd[k]; pa.frrwr[k]=Frrwr[k]; pa.frrrr[k]=Frrrr[k]; pa.frrdr[k]=Frrdr[k]; 
	pa.fdrww[k]=Fdrww[k]; pa.fdrwd[k]=Fdrwd[k]; pa.fdrdd[k]=Fdrdd[k]; pa.fdrwr[k]=Fdrwr[k]; pa.fdrrr[k]=Fdrrr[k]; pa.fdrdr[k]=Fdrdr[k]; 
	};


	return;};

	int random_binomial(int N,double p){
		int ran;
		if(N==0){ran=0;}
		else if(p>0.999999){ran=N;}
		else if(p<0.000001){ran=0;}
		else if(N*p>10 && N*(1-p)>10) {ran=max(0,min(N,int(random_normal(N*p,sqrt(N*p*(1-p))))));}
		else if((N>20 && p<0.05) || (N>100 && N*p<10)){ran=random_poisson(N*p);}
		else if((N>20 && p>0.95) || (N>100 && N*(1-p)<10)){ran=N-random_poisson(N*(1-p));}
		else { variate_generator<mt19937, binomial_distribution<> > b(mt19937(int(rg.Random()*time(NULL))), binomial_distribution<>(N,p)); ran=b();};
		return ran;};


       double random_normal(double mu, double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rg.Random();
	   u2 = rg.Random();
	 }
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(TWOPI* u2);
	return z0 * sigma + mu;
};
double dist (double x1, double y1, double x2, double y2)
	{
		return double(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
	};


	int random_poisson(double landa)
	{
		int k;
		if(landa<1e-5){k=0;}
		else if(landa>30) {k=max(0,(int)random_normal(landa,sqrt(landa)));}
		else
			{
			double p=exp(-landa);
			double g=p;
			double u=rg.Random()*0.999999999;
			k=0;
			while (u>g)
			    {
				p*=(landa/(double)(++k));
				g+=p;
			};
			};
	return k;
	};

	int* random_multinom(int N,long long int probs[6])
{
	int *ran=new int[6];
	double sum=double(probs[0]+probs[1]+probs[2]+probs[3]+probs[4]+probs[5]);
	ran[0]=random_binomial(N,probs[0]/sum);
	ran[1]=random_binomial(N-ran[0],probs[1]/(sum-probs[0]));
	ran[2]=random_binomial(N-ran[0]-ran[1],probs[2]/(sum-probs[0]-probs[1]));
	ran[3]=random_binomial(N-ran[0]-ran[1]-ran[2],probs[3]/(sum-probs[0]-probs[1]-probs[2]));
	ran[4]=random_binomial(N-ran[0]-ran[1]-ran[2]-ran[3],probs[4]/(sum-probs[0]-probs[1]-probs[2]-probs[3]));
	ran[5]=N-ran[0]-ran[1]-ran[2]-ran[3]-ran[4];
	return ran;};

int* random_multinom_var(int N,int howmany,double *relprobs,double tot)
{
	int *ran=new int[howmany];
	double sum=tot;
	int Nused=N;
	for(int pat=0;pat<howmany;pat++)
	{
		if(N>0)
		{
		ran[pat]=random_binomial(Nused,*(relprobs+pat)/sum);
		sum-=*(relprobs+pat);
		Nused-=ran[pat];
		}
		else ran[pat]=0;
	};
	return ran;};

int* random_multinomEqualProb(int N,int howmany)
{
	double prob=1/(1.0*howmany);
	int *ran=new int[howmany];
	double sum=1;
	int Nused=N;
	for(int pat=0;pat<howmany;pat++)
	{
		if(N>0)
		{
		ran[pat]=random_binomial(Nused,prob/sum);
		sum-=prob;
		Nused-=ran[pat];
		}
		else ran[pat]=0;
	};
	return ran;};


void record(int index)
	{	

	for(int pat=0;pat<Pat.size();pat+=in.recPatFreq)
				{
	localinfo<<index<<"   "<<Pat[pat].Mww<<"   "<<Pat[pat].Mwd<<"   "<<Pat[pat].Mdd<<"   "<<Pat[pat].Mwr<<"   "<<Pat[pat].Mrr<<"   "<<Pat[pat].Mdr<<endl;
				};
	return;};




//Below is code for the random number generator---------------------------------------------------------------------------------------------------------------------------------



void CRandomMersenne::Init0(uint32 seed) {
   // Detect computer architecture
   union {float f; uint32 i[2];} convert;
   convert.f = 1.0;
   if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
   else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
   else Architecture = NONIEEE;

   // Seed generator
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(uint32 seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}

uint32 CRandomMersenne::BRandom() {// Generate 32 random bits

   uint32 y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32 LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32 UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32 mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }

   y = mt[mti++];

#if 1
   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;
#endif

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   union {float f; uint32 i[2];} convert; //Union allows one portion of the memory to be accessed as different data types.
   double ra=1.1;
	while(ra>=1)
	{
   uint32 r = BRandom();               // Get 32 random bits
   // The fastest way to convert random bits to floating point is as follows:
   // Set the binary exponent of a floating point number to 1+bias and set
   // the mantissa to random bits. This will give a random number in the
   // interval [1,2). Then subtract 1.0 to get a random number in the interval
   // [0,1). This procedure requires that we know how floating point numbers
   // are stored. The storing method is tested in function RandomInit and saved
   // in the variable Architecture.

   // This shortcut allows the compiler to optimize away the following switch
   // statement for the most common architectures:
#if defined(_M_IX86) || defined(_M_X64) || defined(__LITTLE_ENDIAN__)
   Architecture = LITTLE_ENDIAN1;
#elif defined(__BIG_ENDIAN__)
   Architecture = BIG_ENDIAN1;
#endif

   switch (Architecture) {
   case LITTLE_ENDIAN1:
      convert.i[0] =  r << 20;
      convert.i[1] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case BIG_ENDIAN1:
      convert.i[1] =  r << 20;
      convert.i[0] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case NONIEEE: default: ;
   }
   // This somewhat slower method works for all architectures, including
   // non-IEEE floating point representation:
   //return 0.0000005 +0.999999*(double)r * (1./((double)(uint32)(-1L)+1.));
   ra= (double)r * (1./((double)(uint32)(-1L)+1.));
	};
   return ra;
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((max - min + 1) * Random()) + min;
   if (r > max) r = max;
   return r;
}


