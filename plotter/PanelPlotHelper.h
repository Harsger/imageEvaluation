#ifndef NSWASAP_PANELPLOTHELPER_H
#define NSWASAP_PANELPLOTHELPER_H

#include "TH2.h"
#include "TVector2.h"
#include <map>

namespace NswAsap {
  class PanelPlotHelper {
    public:
      using yMap_t = std::map<std::string, double>;
      using xyMap_t = std::map<std::string, std::pair<double,double>>;

      struct maskParser_t {
	// Default implementation for the gantry scan
	std::string mask_id;
	bool is_strip = false;
        bool toScale = false;
        bool singleArrow = true;
	int color1 = kBlack;
	int color2 = kGray+2;
	virtual void parseMaskId(std::string _mask_id) { mask_id = _mask_id; is_strip=mask_id[13]=='S'; }
	virtual std::string asapId() const { return mask_id; }
	virtual int ipcb() const { return (int)(asapId()[6]-'0'); } // unchecked...
	virtual int imask() const { return (int)(asapId()[17]-'0'); } // unchecked...
	virtual double arrowDirection() const { if (is_strip) return true; return isRI()?-1.:+1.; }
	virtual int color() const { return (imask()<=2) ? color1 : color2; }
	virtual bool isRI() const { return asapId()[11] == 'I'; }
      };

      static void maskPlot(yMap_t yvalues, const std::string& module, const std::string& title, maskParser_t& parser);
      static void pcbPlot(xyMap_t xyvalues, const std::string& module, const std::string& title, maskParser_t& parser);

    private:
      struct board_layout_t {
	std::map<std::string, TVector2> masks;
	std::array<TVector2,4> corners;
      };

      static board_layout_t makeBoardLayout(int ipcb, std::string module);
      static TH2* makeHistoLayout(const std::map<int, board_layout_t>& bmap, std::string title="") ;
      static void plotHisto(TH2* h2);
      static void histoMaskResiduals(TH2* h2, const std::map<int, board_layout_t>& lvec, const yMap_t& yvalues, maskParser_t& parser, double scale=2000.);
      static void histoMaskResiduals(TH2* h2, const std::map<int, board_layout_t>& lvec, const xyMap_t& xyvalues, maskParser_t& parser, double scale=2000.);
  };
}

#endif

