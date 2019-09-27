#include "stdio.h"
#include <iostream>
#include "dsigma.h"

using namespace std;

int read_spectrum(char *fname, float &low, float &bwidth, double *spectrum){
    FILE *file = fopen(fname, "r");
    char sdum[1000];
    double x[1000];
    fscanf(file, "%[^\n]\n", sdum);
    int count = 0;
    while(!feof(file)){
        fscanf(file, "%lf %lf\n", &(x[count]), &(spectrum[count]));
        count++;
    }
    low = x[0];
    bwidth = x[1] - x[0];
    fclose(file);
    return 0;
}

int main(){
    char *fname = "/disk02/usr6/yashida/relic/mc/relic.flux/nakazato/nakazato_nuebar_ref_nh.txt";
    double spectrum[1000];
    float low, bwidth;
    read_spectrum(fname, low, bwidth, spectrum);
    cout << weight(6, low, bwidth, spectrum) << endl;
}
