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

unsigned int estimate[2];
vector<double> fitResults;

RNG rng(12345);
unsigned int blurSize = 3;

int cannyThreshold = 50;
int minRadius = 450;
int maxRadius = 500;
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

void circleFunction( Int_t & , Double_t * , Double_t& f , Double_t* par , Int_t );

void doCircleFit( TGraph * graphTOfit );

void fitContour();
    
int main(int argc, char* argv[]){
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "   matchNfit [options]\n"
        "\n"
        " -i\tinput picture name  \t(default:  \"" << inname << "\" , MANDATORY )\n"
        " -t\ttemplate filename   \t(default:  \"" << templatename << "\" , MANDATORY )\n"
        " -c\tcanny threshold     \t(default:  \"" << cannyThreshold << "\")\n"
        " -l\tminimal fit radius  \t(default:  \"" << minRadius << "\")\n"
        " -u\tmaximal fit radius  \t(default:  \"" << maxRadius << "\")\n"
        " -a\tscalefactor         \t(default:  \"" << alpha << "\")\n"
        " -b\toffsetfactor        \t(default:  \"" << beta << "\")\n"
        " -w\twaittime to show    \t(default:  \"" << waittime << "\")\n"
        " -D\tdebugging mode      \t(default:  \"" << debug << "\")\n"
        "\n";
        return 0;
    }
    
    char c;
    while( ( c = getopt( argc , argv , "i:t:c:l:u:a:b:w:D" ) ) != -1 ){
        switch (c){
            case 'i':
                inname = optarg;
                break;
            case 't':
                templatename = optarg;
                break;
            case 'c':
                cannyThreshold = atof( optarg );
                break;
            case 'l':
                minRadius = atof( optarg );
                break;
            case 'u':
                maxRadius = atof( optarg );
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
    cout << " canny threshold : " << cannyThreshold << endl;
    cout << " radius range    : [ " << minRadius << " , " << maxRadius << " ] "<< endl;
    if( alpha != 1. || beta != 0 ){ 
        cout << " picture will be rescaled with : new = " << alpha << " x original";
        if( beta < 0. ) cout << " - " << abs(beta) << endl;
        else cout << " + "<< beta << endl;
    }
    if(showing) cout << " pictures are shown " << waittime << " seconds " << "\n";
    if(debug) cout << " debugging mode enabled " << "\n";
    
    if(debug){
        outputfile = new TFile("graphsNhistsForPictures.root","RECREATE");
        outputfile->cd();
    }
    
    Mat original = imread(inname);
    
    if( !original.data ){
        cout << " ERROR : image not found \"" << inname << "\"" << endl;
        return 1;
    }
    
    showImage( original , "original.bmp" );
    
    Mat scaled = original.clone();
    if( alpha != 1. || beta != 0 ){ 
        scaled.convertTo(scaled,-1,alpha,beta);
        showImage( scaled , "scaled.bmp" );
    }
        
    globalMat = scaled.clone();
    estimateCenterViaTemplate();
    if(toAbort) return 1;
    
    cout << " estimated center : \t" << estimate[0] << "\t" << estimate[1] << endl;
    
//     Mat scaled = original.clone();
//     if( alpha != 1. || beta != 0 ){ 
//         scaled.convertTo(scaled,-1,alpha,beta);
//         showImage( scaled , "scaled.bmp" );
//     }
//     globalMat = scaled.clone();
    
    fitContour();
    
    if(debug) outputfile->Close();
    
    ofstream outfile;
    outfile.open("fitResults.txt", std::ios_base::app);
    unsigned int waited = 0;
    while( !( outfile.is_open() ) ){
        if( waited > 500 ){
            cout << " can not find file " << endl;
            return 1;
        }
        usleep(100000);
        waited++;
        outfile.open("fitResults.txt", std::ios_base::app);
    }
    
    cout << " mean results ";
    for(unsigned int p=0; p<fitResults.size(); p++) cout << " " << fitResults.at(p);
    cout << endl;
    
    outfile << inname;
    for(unsigned int p=0; p<fitResults.size(); p++) outfile << " " << fitResults.at(p);
    outfile << endl;
    
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
        estimate[0] = 0;
        estimate[1] = 0;
        toAbort = true;
        return;
    }
    
    Mat source = blurNconvert( globalMat );
    Mat tempTOuse = blurNconvert( tempImage );
    
    Mat matching;
    
    matchTemplate( source , tempTOuse , matching , CV_TM_CCOEFF_NORMED );
    
    double minVal , maxVal;
    Point minLoc , maxLoc;
    
    minMaxLoc( matching , &minVal , &maxVal , &minLoc , &maxLoc , Mat() );
    
    matching.convertTo( matching , CV_8UC3 , 255.0 ); 
    showImage( matching , "matching.bmp" );
    
    estimate[0] = maxLoc.x + tempImage.cols/2.;
    estimate[1] = maxLoc.y + tempImage.rows/2.;
        
}

void circleFunction( Int_t & , Double_t * , Double_t& f , Double_t* par , Int_t ) {
    
   Int_t np = pointsTOfit->GetN();
   Double_t *x = pointsTOfit->GetX();
   Double_t *y = pointsTOfit->GetY();
   
   f = 0;
   
   for ( Int_t i=0; i<np; i++){
       
      Double_t u = x[i] - par[0];
      Double_t v = y[i] - par[1];
      
      Double_t dr = par[2] - TMath::Sqrt( u * u + v * v );
      
      f += dr * dr;
      
   }
   
}

void doCircleFit(){
    
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    Double_t arglist[1] = {0};
    arglist[0] = -1;
    fitter->ExecuteCommand("SET PRINT", arglist, 1);
    fitter->SetFCN(circleFunction);
    
    fitter->SetParameter(0, "X" , estimate[0] , 0.01 , estimate[0] - 100 , estimate[0] + 100 );
    fitter->SetParameter(1, "Y" , estimate[1] , 0.01 , estimate[1] - 100 , estimate[1] + 100 );
    fitter->SetParameter(2, "R" , 0.5 * ( minRadius + maxRadius ) , 0.01 , minRadius , maxRadius );

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
    
    fitResults.push_back( amin / ( pointsTOfit->GetN() - nparx ) ); // Chi^2 / NDF
    
}

void fitContour(){
    
    Mat workMat = blurNconvert( globalMat );
    
    vector< vector<Point> > contours;
    vector<Vec4i> hierarchy;
    
    Canny( workMat , workMat , cannyThreshold , 3 * cannyThreshold , 3 );
    
    findContours( workMat , contours , hierarchy , CV_RETR_TREE , CV_CHAIN_APPROX_NONE );
    
    Mat withContours = Mat::zeros( workMat.size(), CV_8UC3 );
    for( int i=0; i<contours.size(); i++){
        Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
        drawContours( withContours, contours, i, color, 2, 8, hierarchy, 0, Point() );
    }
    
    showImage( withContours , "contours.bmp" );
        
    TGraph *** halfCircle = new TGraph**[2];
    TArc *** fitCircle = new TArc**[2];
    
    for(unsigned int s=0; s<2; s++){
        halfCircle[s] = new TGraph*[2];
        fitCircle[s] = new TArc*[2];
        for(unsigned int v=0; v<2; v++){
            halfCircle[s][v] = new TGraph();
        }
    }
        
    double mean[2][2] = { { 0. , 0. } , { 0. , 0. } };
    TString graphName[2][2] = { { "right" , "left" } , { "top" , "bot"} };
    unsigned int maxDistance = (unsigned int)( sqrt( pow( workMat.rows , 2 ) + pow( workMat.cols , 2 ) ) );
    TH1I * distances = new TH1I( "distances" , "distances" , maxDistance , 0 , maxDistance );
    
    vector<Point> consideredPoints;
    
    for(unsigned int c=0; c<contours.size(); c++){
        for(unsigned int p=0; p<contours.at(c).size(); p++){
        
            double value[2] = { (double)contours.at(c).at(p).x , (double)contours.at(c).at(p).y };
            double dif[2] = { value[0] - estimate[0] , value[1] - estimate[1] };
            double distance = sqrt( pow(dif[0],2) + pow(dif[1],2) );
            
            distances->Fill( distance );
            
            if( 
                distance > minRadius && distance < maxRadius 
                && abs( dif[0] ) > excludeCorners && abs( dif[1] ) > excludeCorners
            ){
                
                for(unsigned int s=0; s<2; s++){
                    for(unsigned int v=0; v<2; v++){
                        
                        double parity = 1.;
                        if( v > 0 ) parity = -1.;
                        
                        unsigned int other = 1;
                        if( s > 0 ) other = 0;
                        
                        if( parity * dif[s] > 0 ){
                            
                            mean[s][v] += value[s];
                            halfCircle[s][v]->SetPoint( halfCircle[s][v]->GetN() , value[other] , value[s] );
                            
                        }
                    
                    }
                }
                
                consideredPoints.push_back( Point( contours.at(c).at(p).x , contours.at(c).at(p).y ) );
                
            }
        
        }
    }
    
    if(debug){
        outputfile->cd();
        distances->Write();
        cout << " found # contour points : " << consideredPoints.size();
        for(unsigned int s=0; s<2; s++){
            for(unsigned int v=0; v<2; v++){
                cout << "\t " << graphName[s][v] << " " << halfCircle[s][v]->GetN();
            }
        }
        cout << endl;
    }
    else distances->Delete();
    
    unsigned int index[2][2];
    unsigned int counter=0;
    
    vector< vector<double> > results;
                                
    for(unsigned int s=0; s<2; s++){
        
        unsigned int oldEstimate[2] = { estimate[0] , estimate[1] };
        
        if( s == 0 ){
            estimate[0] = oldEstimate[1];
            estimate[1] = oldEstimate[0];
        }
        
        for(unsigned int v=0; v<2; v++){
            
            index[s][v] = counter;
            counter++;
            
            pointsTOfit = halfCircle[s][v];
            
            doCircleFit();
            results.push_back( fitResults );
            
            if(debug){ 
                
                cout << " " << graphName[s][v];
                for(unsigned int p=0; p<fitResults.size(); p++) cout << "\t" << fitResults.at(p);
                cout << endl;
                
                outputfile->cd();
                
                unsigned int angle[2] = { 0 , 180 };
                if( v == 1 ){
                    angle[0] = 180;
                    angle[1] = 360;
                }
                
                fitCircle[s][v] = new TArc( fitResults.at(0) , fitResults.at(1) , fitResults.at(2) , angle[0] , angle[1] );
                
                TString name = graphName[s][v];
                name += "Points";
                halfCircle[s][v]->SetTitle( graphName[s][v] );
                halfCircle[s][v]->SetName( graphName[s][v] );
                
                name = graphName[s][v];
                name += "Fit";
                
                halfCircle[s][v]->Write();
                fitCircle[s][v]->Write(name);
                
            }
            
        }
        
        if( s == 0 ){ 
            estimate[0] = oldEstimate[0];
            estimate[1] = oldEstimate[1];
        }
        
    }
    
    double meanResults[2][2];
    double meanRadius[2] = { 0. , 0. };
                                
    for(unsigned int s=0; s<2; s++){
        meanResults[s][0] = 0.5 * ( results.at( index[s][0] ).at(0) + results.at( index[s][1] ).at(0) );
        meanResults[s][1] = 0.5 * ( results.at( index[s][0] ).at(3) + results.at( index[s][1] ).at(3) );
        meanRadius[0] += results.at( index[s][0] ).at(2) + results.at( index[s][1] ).at(2);
        meanRadius[1] += results.at( index[s][0] ).at(5) + results.at( index[s][1] ).at(5);
    }
    
    meanRadius[0] /= 4.;
    meanRadius[1] /= 4.;
    
    fitResults.clear();
    
    fitResults.push_back( meanResults[0][0] );
    fitResults.push_back( meanResults[1][0] );
    fitResults.push_back( meanRadius[0] );
    fitResults.push_back( meanResults[0][1] );
    fitResults.push_back( meanResults[1][1] );
    fitResults.push_back( meanRadius[1] );
        
    double lineLength = 30.;
    double lineWidth = 2;
    Mat showMat = globalMat.clone();
    
    line( showMat , Point( estimate[0] , estimate[1] - lineLength ) , Point( estimate[0] , estimate[1] + lineLength ) , Scalar(0,0,255) , lineWidth , 8 , 0 );
    line( showMat , Point( estimate[0] - lineLength , estimate[1] ) , Point( estimate[0] + lineLength , estimate[1] ) , Scalar(0,0,255) , lineWidth , 8 , 0 );
    
//     line( showMat , Point( fitResults.at(0) , fitResults.at(1) - lineLength ) , Point( fitResults.at(0) , fitResults.at(1) + lineLength ) , Scalar(255,0,0) , lineWidth , 8 , 0 ); 
//     line( showMat , Point( fitResults.at(0) - lineLength , fitResults.at(1) ) , Point( fitResults.at(0) + lineLength , fitResults.at(1) ) , Scalar(255,0,0) , lineWidth , 8 , 0 ); 
//     ellipse(showMat, Point( fitResults.at(0) , fitResults.at(1) ), Size( fitResults.at(2) , fitResults.at(2) ) , 0 , 0 , 360 , Scalar(255,0,0) , lineWidth , 8 , 0 ); 
    
    line( showMat , Point( fitResults.at(1) , fitResults.at(0) - lineLength ) , Point( fitResults.at(1) , fitResults.at(0) + lineLength ) , Scalar(255,0,0) , lineWidth , 8 , 0 ); 
    line( showMat , Point( fitResults.at(1) - lineLength , fitResults.at(0) ) , Point( fitResults.at(1) + lineLength , fitResults.at(0) ) , Scalar(255,0,0) , lineWidth , 8 , 0 ); 
    ellipse(showMat, Point( fitResults.at(1) , fitResults.at(0) ), Size( fitResults.at(2) , fitResults.at(2) ) , 0 , 0 , 360 , Scalar(255,0,0) , lineWidth , 8 , 0 ); 
        
    cvtColor( showMat , showMat , CV_BGR2RGB );
    showImage( showMat , "meanFitted.bmp" );
    
    if(debug){
        
        unsigned int angle[2][2][2] = {
            { { 270 , 450 } , {  90 , 270 } },
            { {   0 , 180 } , { 180 , 360 } }
        };
        
        for(unsigned int s=0; s<2; s++){
            
            showMat = globalMat.clone();
            
            for(unsigned int v=0; v<2; v++){
                
                unsigned int color[2] = { 0 , 255 };
                if( v == 1 ) swap( color[0] , color[1] );
                
                unsigned int xNy[2] = { 1 , 0 };
                if( s == 1 ) swap( xNy[0] , xNy[1] );
                
                unsigned int touse[2] = { 1 , 0 };
                if( s == 1 ) swap( touse[0] , touse[1] );
                
                line( showMat , 
                      Point( results.at(index[s][v]).at(xNy[0]) - touse[0] * lineLength , results.at(index[s][v]).at(xNy[1]) - touse[1] * lineLength ) , 
                      Point( results.at(index[s][v]).at(xNy[0]) + touse[0] * lineLength , results.at(index[s][v]).at(xNy[1]) + touse[1] * lineLength ) , 
                      Scalar(color[0],0,color[1]) , lineWidth , 8 , 0 ); 
                    
                ellipse(showMat, 
                        Point( results.at(index[s][v]).at(xNy[0]) , results.at(index[s][v]).at(xNy[1]) ), 
                        Size( results.at(index[s][v]).at(2) , results.at(index[s][v]).at(2) ) , 
                        0 , angle[s][v][0] , angle[s][v][1] , Scalar(color[0],0,color[1]) , lineWidth , 8 , 0 );
                
            }
            
            TString pictureName = "leftNrightFit.bmp";
            if( s == 1 ) pictureName = "upNdownFit.bmp";
            
            showImage( showMat , pictureName );
            
        }
        
    }
    
}
