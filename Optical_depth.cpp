#include <iostream>
#include <fstream>
#include <cmath>
#include <Vector3d.h>
#include "Date.h"
#include "Geoparameter.h"
#include "opt_vec.h"
#include "pmc.h"
#include "pmc_simulation.h"

double Optical_depth(const double lambda, Date date, AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1, const int num_pmc, PMC *pmc){  
   
    /*任意の2点の光学的深さ*/
    double Lambda = lambda;
    int Day = date.doy();
    double delta = opt_vec(r2, r1, Lambda, Day);

    /*Mie散乱の寄与*/
    // Region pmc_region = search_pmc_region(r2, r1);
    // for(int i = 0; i < pmc_region.num; i++){
    //     delta += Mie_opt_depth(pmc_region.r_from[i], pmc_region.r_to[i], num_pmc, pmc);
    // }

    return delta;

}