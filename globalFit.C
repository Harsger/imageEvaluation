#include <stdlib.h>
#include <stdio.h>
#include <algorithm> 
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <initializer_list>

#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TBranch.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TArc.h>

using namespace std;

unsigned int numberOFpoints = 6;
double slopeEstimate = 0.;
bool freeFit = true;

map< string , unsigned int > maskOrder = {
    { "RASFORK_BLOCK_0504" , 0 },
    { "RASFORK_BLOCK_0506" , 1 },
    { "RASFORK_BLOCK_0306" , 2 },
    { "RASFORK_BLOCK_0304" , 3 },
    { "RASFORK_BLOCK_0104" , 4 },
    { "RASFORK_BLOCK_0106" , 5 }
};

vector< vector<double> > nominal = {
    { -793.074 ,  435.2 },
    {  793.074 ,  435.2 },
    {  712.414 ,    0.  },
    { -712.414 ,    0.  },
    { -631.755 , -435.2 },
    {  631.755 , -435.2 }
};

vector< vector<double> > measured;

vector<double> fitResults;

double estimate[2] = { 0. , 0. };
double residuals[6][2];

void SM2function( Int_t & , Double_t * , Double_t& f , Double_t* par , Int_t ) {
   
    f = 0;

    Double_t tangens = par[2];
    Double_t sinus = - tangens / TMath::Sqrt( 1. + tangens * tangens );
    Double_t cosin = - 1. / TMath::Sqrt( 1. + tangens * tangens );

    for ( Int_t i=0; i<numberOFpoints; i++){

        Double_t mX = measured.at(i).at(0) - par[0] ;
        Double_t mY = measured.at(i).at(1) - par[1];
        
        Double_t nX = nominal.at(i).at(0);
        Double_t nY = nominal.at(i).at(1);
        
        Double_t dX = cosin * mX - sinus * mY - nX;
        Double_t dY = sinus * mX + cosin * mY - nY;
        
        residuals[i][0] = dX;
        residuals[i][1] = dY;
        
        Double_t dr = TMath::Sqrt( dX * dX + dY * dY );
        
        f += dr * dr;
        
    }
   
}

void doSM2fit(){
    
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    Double_t arglist[1] = {0};
    arglist[0] = -1;
    fitter->ExecuteCommand("SET PRINT", arglist, 1);
    fitter->SetFCN( SM2function );
    
    fitter->SetParameter(0, "X" , estimate[0] , 0.01 , estimate[0] - 100 , estimate[0] + 100 );
    fitter->SetParameter(1, "Y" , estimate[1] , 0.01 , estimate[1] - 100 , estimate[1] + 100 );
    fitter->SetParameter(2, "P" , 0. , 0.001 , -0.2 , 0.2 );
        
    if( !freeFit ){ 
//         fitter->SetParameter(0, "X" ,   estimate[0] , 0. ,   estimate[0] ,   estimate[0] );
        fitter->SetParameter(1, "Y" ,   estimate[1] , 0. ,   estimate[1] ,   estimate[1] );
        fitter->SetParameter(2, "P" , slopeEstimate , 0. , slopeEstimate , slopeEstimate );
//         fitter->FixParameter(0);
        fitter->FixParameter(1);
        fitter->FixParameter(2);
    }

    fitter->ExecuteCommand("MIGRAD", arglist, 0);
    
    fitResults.clear();
    
    fitResults.push_back( fitter->GetParameter(0) );
    fitResults.push_back( fitter->GetParameter(1) );
    fitResults.push_back( fitter->GetParameter(2) );
    
    fitResults.push_back( fitter->GetParError(0) );
    fitResults.push_back( fitter->GetParError(1) );
    fitResults.push_back( fitter->GetParError(2) );
    
    double amin, edm, errdef;
    int nvpar, nparx;
    fitter->GetStats( amin, edm, errdef, nvpar, nparx);
    
    fitResults.push_back( amin / ( numberOFpoints - nparx ) ); // Chi^2 / NDF
    
}

vector< vector<string> > getInput( string filename ){
    
    vector< vector<string> > input;
    
    if( filename.compare("") == 0 ){ 
        cout << " WARNING : no input to read from " << endl;
        return input;
    }
    
    ifstream ifile(filename.c_str());
    if( !( ifile ) ){ 
        cout << " WARNING : could not read input file " << filename << endl;
        return input;
    }
    
    string line = "";
    string word = "";
    vector<string> dummy;
    
    while( getline( ifile, line) ){
        
        stringstream sline(line);
        
        while( !( sline.eof() ) ){ 
            
            sline >> skipws >> word;
            if( word != "" ) dummy.push_back(word);
            word="";
            
        }
        
        if( dummy.size() > 0 ) input.push_back(dummy);
        dummy.clear();
    }
    
    ifile.close();
    
    return input;
    
}

vector<unsigned int> getSortedIndices(vector<double> order){
   
    vector<unsigned int> sorted;
    unsigned int nPoints = order.size();
    double lower = 0.;
    double lowest = 0.;
    unsigned int index = 0;
    
    if( nPoints < 1 ) return sorted;

    for(unsigned int l=0; l<nPoints; l++){

        if( sorted.size() < 1 ) lowest = order.at(0);
        else{ 
            lower = order.at( sorted.at(l-1) );
            unsigned int newone = 0;
            bool found = false;
            for(unsigned int p=0; p<nPoints; p++){
                bool inlist = false;
                for(unsigned int s=0; s<sorted.size(); s++){
                    if( sorted.at(s) == p ){ 
                        inlist = true;
                        break;
                    }
                }
                if( !inlist){ 
                    newone = p;
                    found = true;
                    break;
                }
            }
            if( !found ){
                cout << " WARNING : no index found " << endl;
                break;
            }
            else{ 
                lowest = order.at(newone);
                index = newone;
            }
        }

        for(unsigned int p=0; p<nPoints; p++){

            if( sorted.size() < 1 && order.at(p) < lowest ){
                lowest = order.at(p);
                index = p;
            }
            if( order.at(p) < lowest && order.at(p) > lower  ){ 
                index = p;
                lowest = order.at(p);
            }

        }

        sorted.push_back( index );

    }
    
    return sorted;
  
}

void globalFit(
    TString dataName,
    bool doFreeFit = false,
    TString destination = "/project/etp4/mherrmann/imageEvaluation/residuals"
){
    freeFit = doFreeFit;
    
    vector< vector<string> > stData = getInput( dataName.Data() );
    if( stData.size() < 1 ){ 
        cout << " ERROR : data not found in " << dataName << " => abort " << endl;
        return 1;
    }
    
    vector<double> devecdummy = { 0. , 0. };
    for(unsigned int p=0; p<nominal.size(); p++) measured.push_back( devecdummy );
    
    for( auto line : stData ){
        measured.at( maskOrder[ line.at(0) ] ).at(0) = atof( line.at(1).c_str() );
        measured.at( maskOrder[ line.at(0) ] ).at(1) = atof( line.at(2).c_str() );
    }
    
    for( auto p : measured ){
        for(unsigned int c=0; c<2; c++) estimate[c] += p.at(c);
    }
    for(unsigned int c=0; c<2; c++) estimate[c] /= (double)numberOFpoints;
    
    if( !freeFit ){
        
        estimate[0] = ( measured.at(2).at(0) + measured.at(3).at(0) ) * 0.5;
        estimate[1] = ( measured.at(2).at(1) + measured.at(3).at(1) ) * 0.5;
        
        double dX = measured.at(2).at(0) - measured.at(3).at(0);
        double dY = measured.at(3).at(1) - measured.at(2).at(1);        // turn by 180 degree => index swapped
        slopeEstimate = dY / dX;
        
    }
    
    doSM2fit();
    
    cout << " fitResults " << endl;
    for(unsigned int c=0; c<3; c++) cout << " \t " << fitResults.at(c) << " +/- " << fitResults.at(c+3) << endl;
    cout << " \t " << fitResults.at( fitResults.size()-1 ) << endl;
    
    TString outname = destination;
    if( !destination.EndsWith( ".txt" ) ){
        
        TString strDummy = dataName;
        strDummy.ReplaceAll( ".txt" , "_residuals.txt" );
        strDummy = strDummy( strDummy.Last('/')+1 , strDummy.Sizeof() );
        
        outname += "/";
        outname += strDummy;
        
    }
    
    ofstream output( outname.Data() ); 
    
    for( auto o : maskOrder ) output << o.first << " " << residuals[ o.second ][0] << " " << residuals[ o.second ][1] << endl;
    
    output.close();
    
    cout << " data written to : " << outname.Data() << endl;
    
}