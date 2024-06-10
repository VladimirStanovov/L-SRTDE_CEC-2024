#include <iostream>
#include <time.h>
#include <fstream>
#include <random>
#include <chrono>

const int Benchmark = 17; // 17 (CEC 2017 and CEC 2024) or 22 (CEC 2022)
//#include "cec22_test_func.cpp"
#include "cec17_test_func.cpp"
//#include "cec17_test_fast_pow.cpp" // Should work faster

//Ctrl+F "CHANGE" in the file and replace 17 with 22
const int ResTsize1 = 29; // number of functions //12 for CEC 2022, 29 for CEC 2024(2017)
const int ResTsize2 = 1001; // number of records per function //1000+1 = 1001 for CEC 2024

using namespace std;
/*typedef std::chrono::high_resolution_clock myclock; // for random seed
myclock::time_point beginning = myclock::now();
myclock::duration d1 = myclock::now() - beginning;
#ifdef __linux__
    unsigned globalseed = d1.count();
#elif _WIN32
    unsigned globalseed = unsigned(time(NULL));
#else

#endif*/
unsigned globalseed = 2024;
unsigned seed1 = globalseed;
unsigned seed2 = globalseed+100;
unsigned seed3 = globalseed+200;
unsigned seed4 = globalseed+300;
std::mt19937 generator_uni_i(seed1);
std::mt19937 generator_uni_r(seed2);
std::mt19937 generator_norm(seed3);
std::mt19937 generator_uni_i_3(seed4);
std::uniform_int_distribution<int> uni_int(0,32768);
std::uniform_real_distribution<double> uni_real(0.0,1.0);
std::normal_distribution<double> norm_dist(0.0,1.0);

int IntRandom(int target) {if(target == 0) return 0; return uni_int(generator_uni_i)%target;}
double Random(double minimal, double maximal){return uni_real(generator_uni_r)*(maximal-minimal)+minimal;}
double NormRand(double mu, double sigma){return norm_dist(generator_norm)*sigma + mu;}

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;
int stepsFEval[ResTsize2-1];
double ResultsArray[ResTsize2];
int LastFEcount;
int NFEval = 0;
int MaxFEval = 0;
int GNVars;
double tempF[1];
double fopt;
char buffer[500];
double globalbest;
bool globalbestinit;
bool TimeComplexity = true;

void qSort2int(double* Mass, int* Mass2, int low, int high)
{
    int i=low;
    int j=high;
    double x=Mass[(low+high)>>1];
    do
    {
        while(Mass[i]<x)    ++i;
        while(Mass[j]>x)    --j;
        if(i<=j)
        {
            double temp=Mass[i];
            Mass[i]=Mass[j];
            Mass[j]=temp;
            int temp2=Mass2[i];
            Mass2[i]=Mass2[j];
            Mass2[j]=temp2;
            i++;    j--;
        }
    } while(i<=j);
    if(low<j)   qSort2int(Mass,Mass2,low,j);
    if(i<high)  qSort2int(Mass,Mass2,i,high);
}

double GetOptimum(int func_num, double* xopt)
{
    FILE *fpt;
    char FileName[50];
    int res = 0;
    //CHANGE!
    sprintf(FileName, "input_data/shift_data_%d.txt",func_num);
    //sprintf(FileName, "input_data_22/shift_data_%d.txt",func_num);
    fpt = fopen(FileName,"r");
    if (fpt==NULL)
        printf("\n Error: Cannot open input file for reading 1 \n");
    for(int k=0;k<GNVars;k++)
        res = fscanf(fpt,"%lf",&xopt[k]);
    fclose(fpt);
    //CHANGE!
    cec17_test_func(xopt, tempF, GNVars, 1, func_num);
    //cec22_test_func(xopt, tempF, GNVars, 1, func_num);
    return tempF[0];
}
void SaveBestValues(int func_index)
{
    double temp = globalbest - fopt;
    if(temp <= 1E-8 && ResultsArray[ResTsize2-1] == MaxFEval)
        ResultsArray[ResTsize2-1] = NFEval;
    for(int stepFEcount=LastFEcount;stepFEcount<ResTsize2-1;stepFEcount++)
    {
        if(NFEval == stepsFEval[stepFEcount])
        {
            if(temp <= 1E-8)
                temp = 0;
            ResultsArray[stepFEcount] = temp;
            LastFEcount = stepFEcount;
        }
    }
}
class Optimizer
{
public:
    int MemorySize;
    int MemoryIter;
    int SuccessFilled;
    int MemoryCurrentIndex;
    int NVars;			    // размерность пространства
    int NIndsCurrent;
    int NIndsFront;
    int NIndsFrontMax;
    int newNIndsFront;
    int PopulSize;
    int func_num;
    int func_index;
    int TheChosenOne;
    int Generation;
    int PFIndex;

    double bestfit;
    double SuccessRate;
    double F;       /*параметры*/
    double Cr;
    double Right;		    // верхняя граница
    double Left;		    // нижняя граница

    double** Popul;	        // массив для частиц
    double** PopulFront;
    double** PopulTemp;
    double* FitArr;		// значения функции пригодности
    double* FitArrCopy;
    double* FitArrFront;
    double* Trial;
    double* tempSuccessCr;
    double* MemoryCr;
    double* FitDelta;
    double* Weights;

    int* Indices;
    int* Indices2;

    void Initialize(int _newNInds, int _newNVars, int _newfunc_num, int _newfunc_index);
    void Clean();
    void MainCycle();
    void FindNSaveBest(bool init, int IndIter);
    void UpdateMemoryCr();
    double MeanWL(double* Vector, double* TempWeights);
    void RemoveWorst(int NInds, int NewNInds);
};
double cec_24(double* HostVector,int func_num)
{
    //CHANGE!
    cec17_test_func(HostVector, tempF, GNVars, 1, func_num);
    //cec22_test_func(HostVector, tempF, GNVars, 1, func_num);	
    NFEval++;
    return tempF[0];
}
void Optimizer::Initialize(int _newNInds, int _newNVars, int _newfunc_num, int _newfunc_index)
{
    NVars = _newNVars;
    NIndsCurrent = _newNInds;
    NIndsFront = _newNInds;
    NIndsFrontMax = _newNInds;
    PopulSize = _newNInds*2;
    Left = -100;
    Right = 100;
    Generation = 0;
    TheChosenOne = 0;
    MemorySize = 5;
    MemoryIter = 0;
    SuccessFilled = 0;
    SuccessRate = 0.5;
    func_num = _newfunc_num;
    func_index = _newfunc_index;
    for(int steps_k=0;steps_k!=ResTsize2-1;steps_k++)
        stepsFEval[steps_k] = 10000.0/double(ResTsize2-1)*GNVars*(steps_k+1);

    Popul = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        Popul[i] = new double[NVars];
    PopulFront = new double*[NIndsFront];
    for(int i=0;i!=NIndsFront;i++)
        PopulFront[i] = new double[NVars];
    PopulTemp = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulTemp[i] = new double[NVars];
    FitArr = new double[PopulSize];
    FitArrCopy = new double[PopulSize];
    FitArrFront = new double[NIndsFront];
    Weights = new double[PopulSize];
    tempSuccessCr = new double[PopulSize];
    FitDelta = new double[PopulSize];
    MemoryCr = new double[MemorySize];
    Trial = new double[NVars];
    Indices = new int[PopulSize];
    Indices2 = new int[PopulSize];

	for (int i = 0; i<PopulSize; i++)
		for (int j = 0; j<NVars; j++)
			Popul[i][j] = Random(Left,Right);
    for(int i=0;i!=PopulSize;i++)
        tempSuccessCr[i] = 0;
    for(int i=0;i!=MemorySize;i++)
        MemoryCr[i] = 1.0;
}
void Optimizer::UpdateMemoryCr()
{
    if(SuccessFilled != 0)
    {
        MemoryCr[MemoryIter] = 0.5*(MeanWL(tempSuccessCr,FitDelta) + MemoryCr[MemoryIter]);
        MemoryIter = (MemoryIter+1)%MemorySize;
    }
}
double Optimizer::MeanWL(double* Vector, double* TempWeights)
{
    double SumWeight = 0;
    double SumSquare = 0;
    double Sum = 0;
    for(int i=0;i!=SuccessFilled;i++)
        SumWeight += TempWeights[i];
    for(int i=0;i!=SuccessFilled;i++)
        Weights[i] = TempWeights[i]/SumWeight;
    for(int i=0;i!=SuccessFilled;i++)
        SumSquare += Weights[i]*Vector[i]*Vector[i];
    for(int i=0;i!=SuccessFilled;i++)
        Sum += Weights[i]*Vector[i];
    if(fabs(Sum) > 1e-8)
        return SumSquare/Sum;
    else
        return 1.0;
}
void Optimizer::FindNSaveBest(bool init, int IndIter)
{
    if(FitArr[IndIter] <= bestfit || init)
        bestfit = FitArr[IndIter];
    if(bestfit < globalbest || init)
	{
		globalbest = bestfit;		
	}
}
void Optimizer::RemoveWorst(int _NIndsFront, int _newNIndsFront)
{
    int PointsToRemove = _NIndsFront - _newNIndsFront;
    for(int L=0;L!=PointsToRemove;L++)
    {
        double WorstFit = FitArrFront[0];
        int WorstNum = 0;
        for(int i=1;i!=_NIndsFront;i++)
        {
            if(FitArrFront[i] > WorstFit)
            {
                WorstFit = FitArrFront[i];
                WorstNum = i;
            }
        }
        for(int i=WorstNum;i!=_NIndsFront-1;i++)
        {
            for(int j=0;j!=NVars;j++)
                PopulFront[i][j] = PopulFront[i+1][j];
            FitArrFront[i] = FitArrFront[i+1];
        }
    }
}
void Optimizer::MainCycle()
{
    vector<double> FitTemp2;
    for(int IndIter=0;IndIter<NIndsFront;IndIter++)
    {
        FitArr[IndIter] = cec_24(Popul[IndIter],func_num);
        FindNSaveBest(IndIter == 0,IndIter);
        if(!globalbestinit || bestfit < globalbest)
        {
            globalbest = bestfit;
            globalbestinit = true;
        }
        SaveBestValues(func_index);
    }
    double minfit = FitArr[0];
    double maxfit = FitArr[0];
    for(int i=0;i!=NIndsFront;i++)
    {
        FitArrCopy[i] = FitArr[i];
        Indices[i] = i;
        maxfit = max(maxfit,FitArr[i]);
        minfit = min(minfit,FitArr[i]);
    }
    if(minfit != maxfit)
        qSort2int(FitArrCopy,Indices,0,NIndsFront-1);
    for(int i=0;i!=NIndsFront;i++)
    {
        for(int j=0;j!=NVars;j++)
            PopulFront[i][j] = Popul[Indices[i]][j];
        FitArrFront[i] = FitArrCopy[i];
    }
    PFIndex = 0;
    while(NFEval < MaxFEval)
    {
        double meanF = 0.4+tanh(SuccessRate*5)*0.25;
        double sigmaF = 0.02;
        minfit = FitArr[0];
        maxfit = FitArr[0];
        for(int i=0;i!=NIndsFront;i++)
        {
            FitArrCopy[i] = FitArr[i];
            Indices[i] = i;
            maxfit = max(maxfit,FitArr[i]);
            minfit = min(minfit,FitArr[i]);
        }
        if(minfit != maxfit)
            qSort2int(FitArrCopy,Indices,0,NIndsFront-1);
        minfit = FitArrFront[0];
        maxfit = FitArrFront[0];
        for(int i=0;i!=NIndsFront;i++)
        {
            FitArrCopy[i] = FitArrFront[i];
            Indices2[i] = i;
            maxfit = max(maxfit,FitArrFront[i]);
            minfit = min(minfit,FitArrFront[i]);
        }
        if(minfit != maxfit)
            qSort2int(FitArrCopy,Indices2,0,NIndsFront-1);
        FitTemp2.resize(NIndsFront);
        for(int i=0;i!=NIndsFront;i++)
            FitTemp2[i] = exp(-double(i)/double(NIndsFront)*3);
        std::discrete_distribution<int> ComponentSelectorFront (FitTemp2.begin(),FitTemp2.end());
        int prand = 0;
        int Rand1 = 0;
        int Rand2 = 0;
        int psizeval = max(2,int(NIndsFront*0.7*exp(-SuccessRate*7)));//int(0.3*NIndsFront));//
        for(int IndIter=0;IndIter<NIndsFront;IndIter++)
        {
            TheChosenOne = IntRandom(NIndsFront);
            MemoryCurrentIndex = IntRandom(MemorySize);
            do
                prand = Indices[IntRandom(psizeval)];
            while(prand == TheChosenOne);
            do
                Rand1 = Indices2[ComponentSelectorFront(generator_uni_i_3)];
            while(Rand1 == prand);
            do
                Rand2 = Indices[IntRandom(NIndsFront)];
            while(Rand2 == prand || Rand2 == Rand1);
            do
                F = NormRand(meanF,sigmaF);
            while(F < 0.0 || F > 1.0);
            Cr = NormRand(MemoryCr[MemoryCurrentIndex],0.05);
            Cr = min(max(Cr,0.0),1.0);
            double ActualCr = 0;
            int WillCrossover = IntRandom(NVars);
            for(int j=0;j!=NVars;j++)
            {
                if(Random(0,1) < Cr || WillCrossover == j)
                {
                    Trial[j] = PopulFront[TheChosenOne][j] + F*(Popul[prand][j] - PopulFront[TheChosenOne][j]) + F*(PopulFront[Rand1][j] - Popul[Rand2][j]);
                    if(Trial[j] < Left)
                        Trial[j] = Random(Left,Right);
                    if(Trial[j] > Right)
                        Trial[j] = Random(Left,Right);
                    ActualCr++;
                }
                else
                    Trial[j] = PopulFront[TheChosenOne][j];
            }
            ActualCr = ActualCr / double(NVars);
            double TempFit = cec_24(Trial,func_num);
            if(TempFit <= FitArrFront[TheChosenOne])
            {
                for(int j=0;j!=NVars;j++)
                {
                    Popul[NIndsCurrent+SuccessFilled][j] = Trial[j];
                    PopulFront[PFIndex][j] = Trial[j];
                }
                FitArr[NIndsCurrent+SuccessFilled] = TempFit;
                FitArrFront[PFIndex] = TempFit;
                FindNSaveBest(false,NIndsCurrent+SuccessFilled);
                tempSuccessCr[SuccessFilled] = ActualCr;//Cr;
                FitDelta[SuccessFilled] = fabs(FitArrFront[TheChosenOne]-TempFit);
                SuccessFilled++;
                PFIndex = (PFIndex + 1)%NIndsFront;
            }
            SaveBestValues(func_index);
        }
        SuccessRate = double(SuccessFilled)/double(NIndsFront);
        newNIndsFront = int(double(4-NIndsFrontMax)/double(MaxFEval)*NFEval + NIndsFrontMax);
        RemoveWorst(NIndsFront,newNIndsFront);
        NIndsFront = newNIndsFront;
        UpdateMemoryCr();
        NIndsCurrent = NIndsFront + SuccessFilled;
        SuccessFilled = 0;
        Generation++;		
        if(NIndsCurrent > NIndsFront)
        {
            minfit = FitArr[0];
            maxfit = FitArr[0];
            for(int i=0;i!=NIndsCurrent;i++)
            {
                Indices[i] = i;
                maxfit = max(maxfit,FitArr[i]);
                minfit = min(minfit,FitArr[i]);
            }
            if(minfit != maxfit)
                qSort2int(FitArr,Indices,0,NIndsCurrent-1);
            NIndsCurrent = NIndsFront;
            for(int i=0;i!=NIndsCurrent;i++)
                for(int j=0;j!=NVars;j++)
                    PopulTemp[i][j] = Popul[Indices[i]][j];
            for(int i=0;i!=NIndsCurrent;i++)
                for(int j=0;j!=NVars;j++)
                    Popul[i][j] = PopulTemp[i][j];
        }
    }
}
void Optimizer::Clean()
{
    delete Trial;
    for(int i=0;i!=PopulSize;i++)
        delete Popul[i];
    for(int i=0;i!=NIndsFrontMax;i++)
        delete PopulFront[i];
    for(int i=0;i!=PopulSize;i++)
        delete PopulTemp[i];
    delete PopulTemp;
    delete Popul;
    delete PopulFront;
    delete FitArr;
    delete FitArrCopy;
    delete FitArrFront;
    delete Indices;
    delete Indices2;
    delete tempSuccessCr;
    delete FitDelta;
    delete MemoryCr;
    delete Weights;
}
int main(int argc, char** argv)
{
    unsigned t0g=clock(),t1g;
    int TotalNRuns;
    if(Benchmark == 17)
        TotalNRuns = 25;//25 for CEC 2024, 51 for CEC 2017
    else
        TotalNRuns = 30;
	
	if(TimeComplexity) //Works for CEC 2017/2024
	{
		ofstream fout_t("time_complexity.txt");
		cout<<"Running time complexity code"<<endl;
		double T1, T2;
		unsigned t1=clock(),t0;
		GNVars = 30;
		MaxFEval = 10000;
		double* xtmp = new double[GNVars];
		for(int j=0;j!=GNVars;j++)
			xtmp[j] = 0;
					
		t0=clock();
		for(int func_num=1;func_num!=31;func_num++)
		{
			if(func_num == 2)
				continue;
			NFEval = 0;
			for(int j=0;j!=MaxFEval;j++)
			{
				// CHANGE !
				cec17_test_func(xtmp, tempF, GNVars, 1, func_num);
				//cec22_test_func(xtmp, tempF, GNVars, 1, func_num);
			}	
			//cout<<func_num<<endl;
		}
		t1=clock()-t0;
		T1 = double(t1)/29.0;			
		cout<<"T1 = "<<T1<<endl;
		fout_t<<"T1 = "<<T1<<endl;						
		
		t0=clock();
		for(int func_num=1;func_num!=31;func_num++)
		{
			if(func_num == 2)
				continue;				
			globalbestinit = false;
			LastFEcount = 0;
			NFEval = 0;
			int PopSize = 20;
			fopt = 100*func_num;				
			Optimizer OptZ;
			OptZ.Initialize(PopSize*GNVars, GNVars, 1, func_num);
			OptZ.MainCycle();
			OptZ.Clean();				
			//cout<<func_num<<endl;
		}
		t1=clock()-t0;
		T2 = double(t1)/29.0;	
		cout<<"T2 = " << T2 << endl;
		fout_t<<"T2 = " << T2 << endl;     
	
		delete xtmp;
	}	
	//Should be 0:4 for CEC 2017, 0:1 for CEC 2022, 1:2 for CEC 2024
    for(int GNVarsIter = 0;GNVarsIter!=2;GNVarsIter++) 
    {
        int maxNFunc;
        if(Benchmark == 17)
        {
            if(GNVarsIter == 0)
            {
                GNVars = 10;
                MaxFEval = 100000;
            }
            if(GNVarsIter == 1)
            {
                GNVars = 30;
                MaxFEval = 300000;
            }
            if(GNVarsIter == 2)
            {
                GNVars = 50;
                MaxFEval = 500000;
            }
            if(GNVarsIter == 3)
            {
                GNVars = 100;
                MaxFEval = 1000000;
            }
            maxNFunc = 30;
        }
        else
        {
            if(GNVarsIter == 0)
            {
                GNVars = 10;
                MaxFEval = 200000;
            }
            if(GNVarsIter == 1)
            {
                GNVars = 20;
                MaxFEval = 1000000;
            }
            maxNFunc = 12;
        }        

        int func_index = 1;
        for(int func_num = 1; func_num < maxNFunc+1; func_num++)
        {
            if(func_num == 2 && Benchmark == 17) // Skip F2
                continue;
            double* xopt = new double[GNVars];
            fopt = 0;
            if(Benchmark == 17)
                //fopt = GetOptimum(func_num,xopt);
				fopt = 100*func_num;
            else
            {
                switch(func_num)
                {
                    case 1:{fopt = 300;  break;}
                    case 2:{fopt = 400;  break;}
                    case 3:{fopt = 600;  break;}
                    case 4:{fopt = 800;  break;}
                    case 5:{fopt = 900;  break;}
                    case 6:{fopt = 1800; break;}
                    case 7:{fopt = 2000; break;}
                    case 8:{fopt = 2200; break;}
                    case 9:{fopt = 2300; break;}
                    case 10:{fopt = 2400;break;}
                    case 11:{fopt = 2600;break;}
                    case 12:{fopt = 2700;break;}
                    default:{fopt = 0;break;}
                }
            }
            sprintf(buffer, "L-SRTDE_F%d_D%d.txt",func_index,GNVars);
            ofstream fout(buffer);
            int PopSize=20;
            for(int run=0;run!=TotalNRuns;run++)
            {
                ResultsArray[ResTsize2-1] = MaxFEval;
                cout<<"func\t"<<func_num<<"\trun\t"<<run<<endl;
                globalbestinit = false;
                LastFEcount = 0;
                NFEval = 0;
                Optimizer OptZ;
                OptZ.Initialize(PopSize*GNVars, GNVars, func_num, func_index);
                OptZ.MainCycle();
                OptZ.Clean();
                for(int j=0;j!=ResTsize2;j++)
                    fout<<ResultsArray[j]<<"\t";
                fout<<endl;
            }
            delete xopt;
            func_index++;
        }
    }
    t1g=clock()-t0g;
    double T0g = t1g/double(CLOCKS_PER_SEC);
    cout << "Time spent: " << T0g << endl;
	return 0;
}
