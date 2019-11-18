#include "TH1F.h"

Double_t typPdf(Double_t *x, Double_t *par){
  Double_t f=par[0]+par[1]*x[0];
  return f;
}
