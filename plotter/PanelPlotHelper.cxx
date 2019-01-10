#include "PanelPlotHelper.h"
#include "TLine.h"
#include "TArrow.h"
#include "TMathText.h"
#include "TLatex.h"
#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <boost/format.hpp>
#include <set>
#include <iostream>
#include <cmath>
using namespace NswAsap;

void PanelPlotHelper::pcbPlot (xyMap_t xyvalues, const std::string& module, const std::string& title, maskParser_t& parser) {
  std::map<int,board_layout_t> bmap;
  std::set<int> pcb_numbers;
  for(auto p : xyvalues) {
    std::string id = p.first;
    parser.parseMaskId(id);
    pcb_numbers.insert(parser.ipcb());
  }
  for (int ipcb : pcb_numbers) {
    bmap[ipcb] = makeBoardLayout(ipcb, module);
  }
  TH2* h2 = makeHistoLayout(bmap, title);

  histoMaskResiduals(h2, bmap, xyvalues, parser);

  plotHisto(h2);
}

void PanelPlotHelper::maskPlot (yMap_t yvalues, const std::string& module, const std::string& title, maskParser_t& parser) {
  std::map<int,board_layout_t> bmap;
  std::set<int> pcb_numbers;
  for(auto p : yvalues) {
    std::string id = p.first;
    parser.parseMaskId(id);
    pcb_numbers.insert(parser.ipcb());
  }
  for (int ipcb : pcb_numbers) {
    bmap[ipcb] = makeBoardLayout(ipcb, module);
  }
  TH2* h2 = makeHistoLayout(bmap, title);

  histoMaskResiduals(h2, bmap, yvalues, parser);

  plotHisto(h2);
}

PanelPlotHelper::board_layout_t PanelPlotHelper::makeBoardLayout (int ipcb, std::string module) {
  board_layout_t ret;

  double x0 = 963.7; // halfway through the sole
  double y0 = 2021.7; // pins of PCB5

  if (module=="LM1")
    y0 -= (5-ipcb) * 460.8; // pins of PCB ipcb
  else
    y0 -= (5-ipcb) * 435.2; // pins of PCB ipcb for SM1 and SM2

  double dx; // distance in x between middle of sole and pins of PCB ipb
  if (module=="LM1")
    dx = 860.025-(5-ipcb)*460.8*tan(16.5*M_PI/180.);
  else if (module=="SM1")
    dx = 542.890-(5-ipcb)*435.2*tan(10.5*M_PI/180.);
  else if (module=="SM2")
    dx = 1586.147*0.5-(8-ipcb)*435.2*tan(10.5*M_PI/180.);
  else throw std::runtime_error("Not implemented");

  std::vector<std::tuple<std::string,double,double>> masks; // coordinates of the 4 masks relative to the pin
  if (module=="LM1")
    masks
      = {
        std::make_tuple("MASK2",-19.,0.),
        std::make_tuple("MASK1",-64.,0.),
        std::make_tuple("MASK3",-110.862,158.205),
        std::make_tuple("MASK4",-17.138, -158.205)
      };
  else if (module=="SM1")
    masks
      = {
	std::make_tuple("MASK2",-29.,0.),
	std::make_tuple("MASK1",-64.,0.),
	std::make_tuple("MASK3",-94.069,162.237),
	std::make_tuple("MASK4",-33.231, -162.237)
      };
  else if (module=="SM2")
    masks
      = {
	std::make_tuple("MASK2",-29.,0.),
	std::make_tuple("MASK1",-64.,0.),
	std::make_tuple("MASK3",-94.069,162.237),
	std::make_tuple("MASK4",-33.931, -162.237)
      };

  int imask = 0;
  for (int ori : {1, -1}) {
    std::string side=(ori==1)?"RI":"RE";
    std::string pcbtag=(module=="LM1")?"LE":"SE";
    for (auto p : masks) {
      ++imask;
      std::string maskid;
      double x, y;
      std::tie(maskid,x,y) = p;
      x = x0 + ori*(-dx) + ori*x;
      y = y0 + y;
      std::string tagString = (boost::format("PCB_%s%1d_B_%s_%s")%pcbtag%ipcb%side%maskid).str();
      ret.masks[tagString] = TVector2(x,y);
    }
  }
  {
    std::vector<std::tuple<double,double>> d_pos;
    std::string pcbtag=(module=="LM1")?"LE":"SE";
    if (module=="LM1") {
      //d_pos = { std::make_tuple(0., 0.), std::make_tuple(0., 158.205), std::make_tuple(0., -158.205) };
      // TMP DEBUG Strip window test: valid for LE4 only.
      d_pos = {
	std::make_tuple(   0.   , 0.   ), std::make_tuple(   0.   , 158.205), std::make_tuple(   0.   , -158.205),
	std::make_tuple(-250.   , 0.   ), std::make_tuple(-250.   , 158.205), std::make_tuple(-250.   , -158.205),
	std::make_tuple( 250.   , 0.   ), std::make_tuple( 250.   , 158.205), std::make_tuple( 250.   , -158.205),
	std::make_tuple(-500.   , 0.   ), std::make_tuple(-500.   , 158.205), std::make_tuple(-500.   , -158.205),
	std::make_tuple( 500.   , 0.   ), std::make_tuple( 500.   , 158.205), std::make_tuple( 500.   , -158.205),
	std::make_tuple(-655.   , 0.   ), std::make_tuple(-705.   , 158.205), std::make_tuple(-610.   , -158.205),
	std::make_tuple( 655.   , 0.   ), std::make_tuple( 705.   , 158.205), std::make_tuple( 610.   , -158.205),
      };
    } else if (module=="SM1")
      d_pos = { std::make_tuple(0., 0.), std::make_tuple(0., 162.237), std::make_tuple(0., -162.237) };
    else if (module=="SM2")
      d_pos = { std::make_tuple(0., 0.), std::make_tuple(0., 162.237), std::make_tuple(0., -162.237) };
    int istrip=0;
    for (auto p: d_pos) {
      double dx, dy;
      std::tie(dx,dy) = p;
      ++istrip;
      double x = x0 + dx;
      double y = y0 + dy;
      std::string tagString = (boost::format("PCB_%s%1d_B_CN_STRI%1d")%pcbtag%ipcb%istrip).str();
      ret.masks[tagString] = TVector2(x,y);
    }
  }

  // Corners of PCB very approximate (to guide the eye)
  double corner_dx, corner_dy;
  if (module=="LM1") {
    corner_dx = 0.5*460.8*tan(16.5*M_PI/180.);
    corner_dy = 0.5*460.8;
  } else {
    corner_dx = 0.5*435.2*tan(10.5*M_PI/180.);
    corner_dy = 0.5*435.2;
  }
  ret.corners[0] = TVector2(x0-dx-75.+corner_dx,y0-corner_dy);
  ret.corners[1] = TVector2(x0+dx+75.-corner_dx,y0-corner_dy);
  ret.corners[2] = TVector2(x0+dx+75.+corner_dx,y0+corner_dy);
  ret.corners[3] = TVector2(x0-dx-75.-corner_dx,y0+corner_dy);
  return ret;
}

TH2* PanelPlotHelper::makeHistoLayout (const std::map<int,board_layout_t>& bmap, std::string title)  {
  double xmin(std::numeric_limits<double>::max()), xmax(std::numeric_limits<double>::min());
  double ymin(std::numeric_limits<double>::max()), ymax(std::numeric_limits<double>::min());
  for (const auto& p : bmap) {
    const board_layout_t& l = p.second;
    for (int i=0; i<4; ++i) {
      xmin = std::min(xmin, l.corners[i].X());
      xmax = std::max(xmax, l.corners[i].X());
      ymin = std::min(ymin, l.corners[i].Y());
      ymax = std::max(ymax, l.corners[i].Y());
    }
  }
  double rx = std::min(0.2*(xmax-xmin), 200.);
  double ry = std::min(0.2*(ymax-ymin), 200.);
  xmin -= rx;
  xmax += rx;
  ymin -= ry;
  ymax += ry;

  TH2F* h = new TH2F("h", title.c_str(), 100, xmin, xmax, 100, ymin, ymax);
  h->SetStats(0);
  h->GetXaxis()->SetTitle("X [mm]");
  h->GetYaxis()->SetTitle("Y [mm]");
  h->SetDirectory(0);
  for (const auto& p : bmap) {
    const board_layout_t& l = p.second;
    for (int i=0; i<4; ++i) {
      const TVector2& a = l.corners[i];
      const TVector2& b = l.corners[(i+1)%4];
      TLine* l = new TLine(a.X(), a.Y(), b.X(), b.Y());
      h->GetListOfFunctions()->Add(l);
    }
  }
  return h;
}

void PanelPlotHelper::plotHisto (TH2* h2) {
  double xmin = h2->GetXaxis()->GetXmin();
  double xmax = h2->GetXaxis()->GetXmax();
  double ymin = h2->GetYaxis()->GetXmin();
  double ymax = h2->GetYaxis()->GetXmax();

  double w, h; // canvas size
  double l_1000 = 250; // 1 meter should appear as 400 pixels
  double rm = 0.05*700.;
  double lm = 0.15*700.;
  double tm = 0.10*500.;
  double bm = 0.15*500.;
  double l_xaxis = (xmax-xmin)/1000.*l_1000;
  double l_yaxis = (ymax-ymin)/1000.*l_1000;
  w = l_xaxis+rm+lm;
  h = l_yaxis+tm+bm;

  double scale = 1.;
//   double scale = 700./double(std::min(w,h));
//   std::cout << "scale is " << scale << std::endl;

//   gStyle->SetTitleFontSize(0.05*scale);
//   h2->GetXaxis()->SetTitleSize(gStyle->GetTitleXSize()*scale);
//   h2->GetXaxis()->SetLabelSize(gStyle->GetLabelSize()*scale);
//   h2->GetYaxis()->SetTitleSize(gStyle->GetTitleYSize()*scale);
//   h2->GetYaxis()->SetLabelSize(gStyle->GetLabelSize()*scale);

  TIter next(h2->GetListOfFunctions());
  TObject* obj;
  while ((obj = next())) {
    TLatex* text = dynamic_cast<TLatex*>(obj);
    if (text) {
//       text->SetTextSize(text->GetTextSize()*scale);
//       text->SetTextSize(gStyle->GetTitleYSize());
      text->SetTextSize(0.03);
    }
  }

//   std::cout << w << " " << h << std::endl;
  TCanvas * c1 = new TCanvas("c", "c", w, h);
  if (gROOT->IsBatch()) {
    c1->SetCanvasSize(w,h);
  } else {
    c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));
  }
  c1->SetRightMargin(rm/w);
  c1->SetLeftMargin(lm/w);
  c1->SetTopMargin(tm/h);
  c1->SetBottomMargin(bm/h);
  c1->SetFixedAspectRatio();
  
  gStyle->SetOptTitle(1);

  h2->Draw();
}

void PanelPlotHelper::histoMaskResiduals (TH2* h2, const std::map<int, board_layout_t>& bmap, const yMap_t& yvalues, maskParser_t& parser, double scale) {
  for (const auto& p : yvalues) {
    parser.parseMaskId(p.first);
    double deltay = p.second;
    double deltay_arrow = parser.arrowDirection()*deltay;

    int col = parser.color();

    const board_layout_t& l = bmap.find(parser.ipcb())->second;
    TVector2 vmask = l.masks.find(parser.asapId())->second;
    TArrow* arr = new TArrow(vmask.X(), vmask.Y(), vmask.X(), vmask.Y()+scale*deltay_arrow, 0.01, "|>");
    arr->SetLineColor(col);
    arr->SetFillColor(col);
    h2->GetListOfFunctions()->Add(arr);

    bool do_right_align = parser.isRI();
    if (parser.imask() == 2) do_right_align = !do_right_align;

    double offset = (do_right_align ? -10. : + 10.);
    TLatex* text = new TLatex(vmask.X()+offset, vmask.Y(), (boost::format("%.0f#mum")%(deltay*1000.)).str().c_str());
    text->SetTextFont(42);
    text->SetTextSize(0.02);
    if (do_right_align)
      text->SetTextAlign(32); // right adjusted
    else
      text->SetTextAlign(12); // left adjusted
    text->SetTextColor(col);
    h2->GetListOfFunctions()->Add(text);
  }
}

void PanelPlotHelper::histoMaskResiduals (TH2* h2, const std::map<int, board_layout_t>& bmap, const xyMap_t& xyvalues, maskParser_t& parser, double scale) {
  for (const auto& p : xyvalues) {
    parser.parseMaskId(p.first);
    double deltax = p.second.first;
    double deltay = p.second.second;
    double deltax_arrow = parser.arrowDirection()*deltax;
    double deltay_arrow = parser.arrowDirection()*deltay;
    
    if( !parser.toScale ){
        deltax_arrow = parser.arrowDirection()*0.02*deltax_arrow/abs(deltax_arrow);
        deltay_arrow = parser.arrowDirection()*0.02*deltay_arrow/abs(deltay_arrow);
    }

    int col = parser.color();

    const board_layout_t& l = bmap.find(parser.ipcb())->second;
    TVector2 vmask = l.masks.find(parser.asapId())->second;
    TArrow* arr;
//     if (std::isnan(deltax))
//       arr = new TArrow(vmask.X(), vmask.Y(), vmask.X(), vmask.Y()+scale*deltay_arrow, 0.01, "|>");
//     else
//       arr = new TArrow(vmask.X(), vmask.Y(), vmask.X()+scale*deltax_arrow, vmask.Y()+scale*deltay_arrow, 0.01, "|>");
    if( parser.singleArrow ){
        arr = new TArrow(vmask.X(), vmask.Y(), vmask.X()+scale*deltax_arrow, vmask.Y()+scale*deltay_arrow, 0.01, "|>");
        arr->SetLineColor(col);
        arr->SetFillColor(col);
        h2->GetListOfFunctions()->Add(arr);
    }
    else{
        arr = new TArrow(vmask.X(), vmask.Y(), vmask.X()+scale*deltax_arrow, vmask.Y(), 0.01, "|>");
        arr->SetLineColor(col);
        arr->SetFillColor(col);
        h2->GetListOfFunctions()->Add(arr);
        
        arr = new TArrow(vmask.X(), vmask.Y(), vmask.X(), vmask.Y()+scale*deltay_arrow, 0.01, "|>");
        arr->SetLineColor(col);
        arr->SetFillColor(col);
        h2->GetListOfFunctions()->Add(arr);
    }

    bool do_right_align = parser.isRI();
    if (parser.imask() == 2) do_right_align = !do_right_align;

    double offset = (do_right_align ? -10. : + 10.);
    TLatex* text;
    if (std::isnan(deltax))
      text = new TLatex(vmask.X()+offset, vmask.Y(), (boost::format("%.0f")%(deltay*1000.)).str().c_str());
    else
      text = new TLatex(vmask.X()+offset, vmask.Y(), (boost::format("(%.0f,%.0f)")%(deltax*1000.)%(deltay*1000.)).str().c_str());
    text->SetTextFont(42);
    text->SetTextSize(0.02);
    if (do_right_align)
      text->SetTextAlign(32); // right adjusted
    else
      text->SetTextAlign(12); // left adjusted
    text->SetTextColor(col);
    h2->GetListOfFunctions()->Add(text);
  }
}

