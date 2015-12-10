///////////////////////////////////////////////////////////////////////////////
#ifndef __CONTIG_SPECTRUM_H__
#define __CONTIG_SPECTRUM_H__
///////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "RendererBase.h"
#include "ReportDefines.h"
#include "PlotBase.h"
#include "abruijn.h"
///////////////////////////////////////////////////////////////////////////////
using namespace std;
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
//
struct abMassItem {
  double  mass;
  int     idxInContig;
};
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
class ContigSpectrum : public PlotBase {

   // spectrum to draw
  specnets::Spectrum *m_spectrum;

  // ABinfo data set
  //specnets::abinfo_t *m_abinfo;
  abContigData_t  *abData;

  // vector if masses from the abruijn for the star
  vector<abMassItem> m_abMasses;

  // sps_seqs data, properly reversed, if the case
  vector<double> *m_deNovoIntervals;

  // Contig to draw
  int m_contigIndex;

  // spectrum index
  int m_spectrumIndex;


  // build mass value list from abruijn graph masses
  void buildVector(void);
  // find mass value in mass value list
  int find(double &value);

  // red horizontal line
  void addRedArrow(double x1, double y1, double x2, double y2);
  // gren horizontal line
  void addGreenArrow(double x1, double y1, double x2, double y2);
  // add the red label with the mass shift value
  void addRedLabel(double x1, double y1, string &label);



 public:

  ContigSpectrum(RendererBase *renderer, int ci, abContigData_t *a, specnets::Spectrum *s, vector<double> *seqs, int idx) :
    m_contigIndex(ci), m_spectrum(s), abData(a), m_spectrumIndex(idx), m_deNovoIntervals(seqs)
    {
      m_rendererObject              = renderer;
      m_ownRenderer                 = false;
      m_rendererImage.lineThickness = 4;
    };

  ~ContigSpectrum() {};

  virtual void draw(double limit, double factor = 1.0);

  virtual void drawExec(void) {};

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
