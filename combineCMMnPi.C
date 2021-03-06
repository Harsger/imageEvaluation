#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TBranch.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TProfile.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

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

void combineCMMnPi(
    TString dataName,
    TString refName,
    TString destination = "/project/etp4/mherrmann/imageEvaluation/cmmPoints",
    bool flip = true,
    double cmmScale = 0.001
){
    
    vector< vector<string> > stData = getInput( dataName.Data() );
    if( stData.size() < 1 ) return 1;
    
    vector< vector<string> > stReference = getInput( refName.Data() );
    if( stReference.size() < 1 ) return 1;
    
    map< unsigned int , string > imageCode;
    imageCode[2] = "RASFORK_BLOCK_0504";
    imageCode[3] = "RASFORK_BLOCK_0506";
    imageCode[4] = "RASFORK_BLOCK_0306";
    imageCode[5] = "RASFORK_BLOCK_0304";
    imageCode[6] = "RASFORK_BLOCK_0104";
    imageCode[7] = "RASFORK_BLOCK_0106";
    
    map< string , pair< double, double > > data;
    map< string , pair< double, double > > reference;
    
    for( auto line : stData ){
        
        TString imageNumber = line.at(0);
        
        data[ imageNumber.Data() ].first = atof( line.at(1).c_str() );
        data[ imageNumber.Data() ].second = atof( line.at(2).c_str() );
        
    }
    
    for(unsigned int i=2; i<stReference.size(); i++){
        
        reference[ imageCode[i] ].first = atof( stReference.at(i).at(1).c_str() ) * cmmScale;
        reference[ imageCode[i] ].second = atof( stReference.at(i).at(2).c_str() ) * cmmScale;
        
    }
    
    map< string , pair< double, double > > result;
    
    double flipper = 1.;
    if( flip ) flipper = -1.;
    
    for( auto d : data ){
        
        result[ d.first ].first = reference[ d.first ].first + flipper * d.second.first;
        result[ d.first ].second = reference[ d.first ].second + flipper * d.second.second;
        
    }
    
    TString outname = destination;
    if( !destination.EndsWith( ".txt" ) ){
        
        TString strDummy = dataName;
        strDummy.ReplaceAll( ".txt" , "_cmm.txt" );
        strDummy = strDummy( strDummy.Last('/')+1 , strDummy.Sizeof() );
        
        outname += "/";
        outname += strDummy;
        
    }
    
    ofstream output( outname.Data() ); 
    
    for( auto r : result ) output << r.first << " " << r.second.first << " " << r.second.second << endl;
    
    output.close();
    
    cout << " data written to : " << outname.Data() << endl;
    
}

