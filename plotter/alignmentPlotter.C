#include "PanelPlotHelper.h"
#include "AtlasUtils.h"
#include <boost/format.hpp>
#include <map>
#include <fstream>

// #ifndef __CLING__
#include "json.hpp"
// #endif

// Include external code here
#include "PanelPlotHelper.cxx"
#include "AtlasUtils.C"

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
    unsigned int lower = 0;
    unsigned int lowest = 0;
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

using xyMap_t = std::map<std::string, std::pair<double,double>>;

xyMap_t read_mask_residuals(std::string filename, string layer="") {
  xyMap_t xvalues;
  std::ifstream in(filename);

// #ifndef __CLING__
  using json = nlohmann::json;
  using parse_error = nlohmann::detail::parse_error;
  json j;
  while (in) {
    try {
      in >> j; // Throws a parse_error at end of file...
    } catch (parse_error&) {
      break;
    }
    std::string blockId = j["blockId"];
    if (blockId.find("CALIRASFORK")!=std::string::npos) continue; // Do not remember old CALIRASFORK data
    if (!(bool)(j["fitStatus"])) continue;

    string measurmentName = "meas";
    measurmentName += layer;
    measurmentName += "_x";
    double xval = j["results"][measurmentName.c_str()];
    xvalues[blockId].first = xval;
    measurmentName = "meas";
    measurmentName += layer;
    measurmentName += "_y";
    double yval = j["results"][measurmentName.c_str()];
    xvalues[blockId].second = yval;
  }
// #endif

  return xvalues;
}

struct RasforkParser : NswAsap::PanelPlotHelper::maskParser_t {
  double direction = 1.; // 1 for bottom mask fit, -1 for top mask fit
  std::string tag = "LE"; // "SE" for SM1 and SM2
  int start_pcb = 1; // 6 for SM2
//   std::string prefix="RASFORK_BLOCK";
  string prefix="RASFORK_BLOCK";
  virtual void parseMaskId(std::string _mask_id) {
    //std::cout << "parseMaskId: " << _mask_id << std::endl;
    // 012345678901234567
    // RASFORK_BLOCK_0306
    // PCB_LE1_B_RI_MASK3
    int ioff = prefix.length();
    int ipcb = std::stoi(_mask_id.substr(ioff+1,2));
    ipcb = (ipcb+1)/2 + (start_pcb-1);
    int imask = std::stoi(_mask_id.substr(ioff+3,2));
    std::string side = ((imask%4)==0) ? "RI" : "RE";
    imask = imask/4;
    switch (imask) {
      case 0:
	imask = 4;
	break;
      case 1:
	imask = 1;
	break;
      case 2:
	imask = 3;
	break;
    }
    mask_id = (boost::format("PCB_%s%d_B_%s_MASK%d")%tag%ipcb%side%imask).str();
    //std::cout << "Returning " << mask_id << std::endl;
  }
  virtual double arrowDirection() const {
    return direction*(isRI()?-1.:+1.);
  }
};

void alignmentPlotter(
    TString filename,
    bool toScale = false,
    bool singleArrow = true
){
    
    TString panelID = filename;
    panelID = panelID( panelID.Last('/')+1 , panelID.Sizeof() );
    panelID.ReplaceAll( "camOldPinsNmarker_" , "" );
    panelID.ReplaceAll( "camLeftPinsNmarker_" , "" );
    panelID.ReplaceAll( "camRightPinsNmarker_" , "" );
    panelID.ReplaceAll( "_fits" , "" );
    panelID.ReplaceAll( "_dif" , "" );
    panelID.ReplaceAll( "_cmm" , "" );
    panelID.ReplaceAll( "_residuals" , "" );
    panelID.ReplaceAll( ".txt" , "" );
    
    vector< vector<string> > data = getInput( filename.Data() );
    
    xyMap_t xyvalues;
    
    for( auto d : data ){
        xyvalues[ d.at(0) ].first = atof( d.at(1).c_str() );
        xyvalues[ d.at(0) ].second = atof( d.at(2).c_str() );
        cout << d.at(0) << " \t " << xyvalues[ d.at(0) ].first << " \t " << xyvalues[ d.at(0) ].second << endl;
        if( xyvalues[ d.at(0) ].first == 0. ) xyvalues[ d.at(0) ].first = 1e-6;
        if( xyvalues[ d.at(0) ].second == 0. ) xyvalues[ d.at(0) ].second = 1e-6;
    }
    
    RasforkParser parser;
    parser.direction = -1.; 
    parser.tag = "SE"; // for SM1 or SM2
    parser.toScale = toScale;
    parser.singleArrow = singleArrow;
        
    NswAsap::PanelPlotHelper::pcbPlot( xyvalues , "SM2", panelID.Data() , parser );
        
    string forPad = panelID.Data();
    savePad( forPad.c_str() );

}

