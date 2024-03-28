#include "LinearInterpolation.hh"

LinearInterpolation::LinearInterpolation(Double_t MXX[], Double_t MYY[], Int_t nbins) : MASSXX(new Double_t[nbins]), MASSYY(new Double_t[nbins]), NUMBINS(nbins)
{
    YMIN = MYY[0]; YMAX = YMIN;
    for(int i = 0; i < NUMBINS; i++)
    {
       MASSXX[i] = MXX[i];
       MASSYY[i] = MYY[i];
       YMIN = (MYY[i] < YMIN) ? MYY[i] : YMIN;
       YMAX = (MYY[i] > YMAX) ? MYY[i] : YMAX;
    }
}

LinearInterpolation::LinearInterpolation(std::vector<Double_t> vx, std::vector<Double_t> vy) : MASSXX(new Double_t[vx.size()]), MASSYY(new Double_t[vx.size()]), NUMBINS(vx.size())
{
    YMIN = vy[0]; YMAX = YMIN;
    for(int i = 0; i < NUMBINS; i++)
    {
       MASSXX[i] = vx[i];
       MASSYY[i] = vy[i];
       YMIN = (vy[i] < YMIN) ? vy[i] : YMIN;
       YMAX = (vy[i] > YMAX) ? vy[i] : YMAX;
    }
}

LinearInterpolation* LinearInterpolation::Initialization(TString innam, Double_t unit1, Double_t unit2)
{
    FILE *file = fopen(innam, "r");
    if (!file)
    {
        std::cout << "\n" << "LinearInterpolation: Wrong name of the input file of data for interpolation! \n" << "File name: " << innam << "\n\n";

        return nullptr;
    }
    else
    {
        fseek(file, 0, SEEK_END);
        if(ftell(file) == 0) 
        {
            std::cout << "LinearInterpolation: Empty file - " << innam << "\n"; 
            fclose(file); 
            return nullptr;
        }
        {
            fseek(file, 0, SEEK_SET);
            return new LinearInterpolation(file, unit1, unit2);
        }
    }
}

LinearInterpolation::LinearInterpolation(FILE* file, Double_t unit1, Double_t unit2)
{
    std::vector<Double_t> list_x, list_y;

    double xx, yy;

    char STR[256];
    while (STR == fgets(STR, 255, file))
    {     
        int nit = sscanf(STR, "%lf %lf", &xx, &yy);
        if(nit > 1)
        {
            list_x.push_back(xx * unit1);
            list_y.push_back(yy * unit2);
        }
        
        list_x.push_back(xx * unit1);
        list_y.push_back(yy * unit2);
    }

    NUMBINS = list_x.size();
    MASSXX = new Double_t[NUMBINS];
    MASSYY = new Double_t[NUMBINS];

    YMIN = list_y[0]; YMAX = YMIN;
    for(int i = 0; i < NUMBINS; i++)
    {
       MASSXX[i] = list_x[i];
       MASSYY[i] = list_y[i];
       YMIN = (list_y[i] < YMIN) ? list_y[i] : YMIN;
       YMAX = (list_y[i] > YMAX) ? list_y[i] : YMAX;
    }

    fclose(file); 
}

LinearInterpolation::~LinearInterpolation()
{
    delete [] MASSXX;
    delete [] MASSYY;
}

Double_t LinearInterpolation::GetValue(Double_t xx)
{
    if(NUMBINS == 0) {std::cout << "LinearInterpolation: Empty  arrays!\n"; return -1;}
    if(NUMBINS == 1) return MASSYY[0];

    if(xx >= MASSXX[NUMBINS-1]) return MASSYY[NUMBINS-1];
    if(xx <= MASSXX[0]) return MASSYY[0];
    
    int k = 0;
    for (int i = 0; i < NUMBINS-1; i++)
    {
        if(xx < MASSXX[i+1]) {k = i; break;}
    }

    return MASSYY[k] + (xx - MASSXX[k]) * ((MASSYY[k+1] - MASSYY[k])/(MASSXX[k+1] - MASSXX[k]));
}

Double_t LinearInterpolation::GetXmin()
{
    return MASSXX[0];
}

Double_t LinearInterpolation::GetXmax()
{
    return MASSXX[NUMBINS-1];
}

Double_t LinearInterpolation::GetYmin()
{
    return YMIN;
}

Double_t LinearInterpolation::GetYmax()
{
    return YMAX;
}

Int_t LinearInterpolation::GetNBins()
{
    return NUMBINS;
}

Double_t* LinearInterpolation::GetXMass()
{
    return MASSXX;
}

Double_t* LinearInterpolation::GetYMass()
{
    return MASSYY;
}
