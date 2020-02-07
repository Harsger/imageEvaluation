#include <stdlib.h>
#include <stdio.h>
#include <algorithm> 
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <initializer_list>

#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

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

using namespace cv;
using namespace std;

Mat globalMat;
TGraph * pointsTOfit;
TFile * outputfile; 

string inname = "";
string templatename = "";

bool debug = false;
bool showing = false;
bool toAbort = false;

unsigned int bestEdge[2];
unsigned int templateSize[2];

RNG rng(12345);
unsigned int blurSize = 3;

double alpha = 1.;
int beta = 0;
unsigned int excludeCorners = 25;

unsigned int waittime = 0;

namespace patch{
    
    template < typename T > string to_string( const T& n ){
        
        ostringstream stm;
        stm << n;
        return stm.str();
        
    }
    
}

void showImage( Mat src , TString saveName="" );

Mat blurNconvert( Mat src );

void estimateCenterViaTemplate();
    
int main(int argc, char* argv[]){
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "   clip [options]\n"
        "\n"
        " -i\tinput picture name  \t(default:  \"" << inname << "\" , MANDATORY )\n"
        " -t\ttemplate filename   \t(default:  \"" << templatename << "\" , MANDATORY )\n"
        " -a\tscalefactor         \t(default:  \"" << alpha << "\")\n"
        " -b\toffsetfactor        \t(default:  \"" << beta << "\")\n"
        " -w\twaittime to show    \t(default:  \"" << waittime << "\")\n"
        " -D\tdebugging mode      \t(default:  \"" << debug << "\")\n"
        "\n";
        return 0;
    }
    
    char c;
    while( ( c = getopt( argc , argv , "i:t:a:b:w:D" ) ) != -1 ){
        switch (c){
            case 'i':
                inname = optarg;
                break;
            case 't':
                templatename = optarg;
                break;
            case 'a':
                alpha = atof( optarg );
                break;
            case 'b':
                beta = atoi( optarg );
                break;
            case 'w':
                waittime = atoi( optarg );
                showing = true;
                break;
            case 'D':
                debug = true;
                break;
            case '?':
                if( isprint(optopt) ) fprintf( stderr , " Unknown option `-%c'.\n" , optopt );
                else fprintf( stderr , " Unknown option character `\\x%x'.\n" , optopt );
                return 1;
            default:
                abort();
        }
    }
    
    if( inname == "" || templatename == "" ){
        cout << " inputname and templatename required " << endl;
        return 1;
    }
  
    cout << " inputname       : " << inname << "\n";
    cout << " templatename    : " << templatename << "\n";
    if( alpha != 1. || beta != 0 ){ 
        cout << " picture will be rescaled with : new = " << alpha << " x original";
        if( beta < 0. ) cout << " - " << abs(beta) << endl;
        else cout << " + "<< beta << endl;
    }
    if(showing) cout << " pictures are shown " << waittime << " seconds " << "\n";
    if(debug) cout << " debugging mode enabled " << "\n";
    
    Mat original = imread(inname, IMREAD_COLOR);
    
    if( !original.data ){
        cout << " ERROR : image not found \"" << inname << "\"" << endl;
        return 1;
    }
    
    templateSize[0] = original.cols;
    templateSize[1] = original.rows;
    
    showImage( original , "original.bmp" );
    
    Mat scaled = original.clone();
    if( alpha != 1. || beta != 0 ){ 
        scaled.convertTo(scaled,-1,alpha,beta);
        showImage( scaled , "scaled.bmp" );
    }
        
    globalMat = scaled.clone();
    estimateCenterViaTemplate();
    if(toAbort) return 1;
    
    cout << " estimated center : \t" << bestEdge[0] << "\t" << bestEdge[1] << endl;
    
    ofstream outfile;
    outfile.open("matchResults.txt", std::ios_base::app);
    unsigned int waited = 0;
    while( !( outfile.is_open() ) ){
        if( waited > 500 ){
            cout << " can not find file " << endl;
            return 1;
        }
        usleep(100000);
        waited++;
        outfile.open("matchResults.txt", std::ios_base::app);
    }
    
    outfile << inname;
    for(unsigned int p=0; p<2; p++) outfile << " " << bestEdge[p];
    for(unsigned int p=0; p<2; p++) outfile << " " << templateSize[p];
    outfile << endl;
    
    Rect roi;
    roi.x = bestEdge[0];
    roi.y = bestEdge[1];
    roi.width  = templateSize[0];
    roi.height = templateSize[1];
    
    Mat crop = original(roi);
    
    showImage( crop , "cropped.bmp" );
    
    return 0;
    
}

void showImage( Mat src , TString saveName ){
        
    if( debug && ( saveName.EndsWith(".bmp") || saveName.EndsWith(".png") ) ) imwrite( saveName.Data() , src );
    
    if( showing ){
        
        TString windowName = saveName;
        windowName = windowName.ReplaceAll( ".bmp" , "" );
        windowName = windowName.ReplaceAll( ".png" , "" );
        
        if(debug) cout << " SHOWING " << windowName << endl;
        
        namedWindow( windowName.Data() , WINDOW_NORMAL );
        resizeWindow( windowName.Data() , 1000 , 1000 );
        imshow( windowName.Data() , src );
        waitKey(100);
        sleep( waittime );
        destroyWindow( windowName.Data() );
        
    }
    
}

Mat blurNconvert( Mat src ){
    
    Mat destiny = src.clone();
    
    if(debug){ 
        cout << " size : \t x " << destiny.cols << " \t y " << destiny.rows;
        cout << " \t depth : " << destiny.depth() << "\t channels : " << destiny.channels() << endl;
    }
    
    GaussianBlur( destiny , destiny , Size( blurSize , blurSize ) , 0 , 0 , BORDER_DEFAULT );
    cvtColor( destiny , destiny , CV_BGR2GRAY );
    
    return destiny;
    
}

void estimateCenterViaTemplate(){
    
    Mat tempImage = imread( templatename );
    
    if( tempImage.empty() ){ 
        cout << " ERROR: template file not found for : " << templatename << endl;
        bestEdge[0] = 0;
        bestEdge[1] = 0;
        toAbort = true;
        return;
    }
    
    Mat source = blurNconvert( globalMat );
    Mat tempTOuse = blurNconvert( tempImage );
    
    templateSize[0] = tempTOuse.cols;
    templateSize[1] = tempTOuse.rows;
    
    Mat matching;
    
    matchTemplate( source , tempTOuse , matching , CV_TM_CCOEFF_NORMED );
    
    double minVal , maxVal;
    Point minLoc , maxLoc;
    
    minMaxLoc( matching , &minVal , &maxVal , &minLoc , &maxLoc , Mat() );
    
    matching.convertTo( matching , CV_8UC3 , 255.0 ); 
    showImage( matching , "matching.bmp" );
    
    bestEdge[0] = maxLoc.x;
    bestEdge[1] = maxLoc.y;
        
}
