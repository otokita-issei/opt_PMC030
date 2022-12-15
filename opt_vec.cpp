#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector3d.h"
#include "opt_vec.h"
#include "Geoparameter.h"

double opt_vec(AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1){  
    // AndoLab::Vector3d <double> r1;
    // AndoLab::Vector3d <double> r2;
    AndoLab::Vector3d <double> d;
    AndoLab::Vector3d <double> tmp;

    // r1.set(r_1, th1, phi1, AndoLab::coordinate::Spherical);
    // r2.set(r_2, th2, phi2, AndoLab::coordinate::Spherical);

    if(r1.r() > r2.r()){    /*2つのベクトルを条件によって交換する*/
        tmp = r1;
        r1 = r2;
        r2 = tmp;
    }

    d = (r2 - r1) / abs(r2 - r1);

    double r1_d = r1 % d;

    /*経路1の計算(delta1)*/
    int z0, z_max, N_z;
    double xi0, r_1, t1, delta1;

    /*r1が自然数であった場合の対処*/
    if(r1.r() / 1e3 - int(r1.r() / 1e3) == 0 ){     
        z0 = (r1.r() - RADIUS_OF_EARTH)*MtoKM;
        t1 = 0.0;
    } else {
        z0 = (r1.r() - RADIUS_OF_EARTH)*MtoKM + 1;
    }

    z_max = (r2.r() - RADIUS_OF_EARTH)*MtoKM;
    r_1 = r1.r();

    /*z0 > z_maxのとき (小さい高度差 (xx.y - xx.z[km]))*/
    if(z0 > z_max){
        N_z = 0;
        xi0 = r2.r(); 
        t1 = -1.0*r1_d + std::sqrt(r1_d*r1_d + xi0*xi0 - r_1*r_1) - (-1.0*r1_d + std::sqrt(r1_d*r1_d)); /*[m]*/
        delta1 = opt_vertical(z0-1, (r1.latitude() + r2.latitude())/2, (r1.longitude() + r2.longitude())/2) * t1 / 1e3; /*opt_vertical(z, Lat, Lon, date)*/
    } else {
        N_z = z_max - z0;
        xi0 = z0 * KMtoM + RADIUS_OF_EARTH;      /* xi0 = |r1 + t1d|*/
        t1 = -1.0*r1_d + std::sqrt(r1_d*r1_d + xi0*xi0 - r_1*r_1) - (-1.0*r1_d + std::sqrt(r1_d*r1_d)); /*[m]*/
        delta1 = opt_vertical(z0-1, r1.latitude(), r1.longitude()) * t1 / 1e3;
    }

    /*高度が同じ場合*/
    if(r2.r() - r1.r() == 0){
        N_z = 0;
        t1 = std::sqrt((r1.x()-r2.x())*(r1.x()-r2.x()) + (r1.y()-r2.y())*(r1.y()-r2.y()) + (r1.z()-r2.z())*(r1.z()-r2.z())); 
        AndoLab::Vector3d <double> p = r1 + t1*d/2.0;
        delta1 = opt_vertical(int((p.r()-RADIUS_OF_EARTH)*MtoKM), (r1.latitude() + r2.latitude())/2, (r1.longitude() + r2.longitude())/2) * t1 / 1e3;
    } else {
        delta1 = delta1;
    }  
 
    /*経路2の計算(delta2)*/
    int z_i;
    double xi_i, xi_i_1, t_i, t_i_1;
    double delta2 = 0.0;

    for(int i = 0; i < N_z; i++){
        z_i = z0 + i; 
        xi_i = (z_i + 1) * KMtoM + RADIUS_OF_EARTH;
        xi_i_1 = z_i * KMtoM + RADIUS_OF_EARTH;
        t_i = -1.0*r1_d + std::sqrt(r1_d*r1_d + xi_i*xi_i - r_1*r_1);
        t_i_1 = -1.0*r1_d + std::sqrt(r1_d*r1_d + xi_i_1*xi_i_1 - r_1*r_1);

        delta2 = delta2 + opt_vertical(z_i, (r1 + t_i*d).latitude(), (r1 + t_i*d).longitude()) * (t_i - t_i_1) / 1e3;
    }

    /*経路3の計算(delta3)*/
    double xi_max, t_max, t_R, delta3;

    xi_max = z_max * KMtoM + RADIUS_OF_EARTH;
    t_max = -1.0*r1_d + std::sqrt(r1_d*r1_d + xi_max*xi_max - r_1*r_1);
    
    /*r2が自然数であった場合の対処*/
    if(r2.r() / 1e3 - int(r2.r() / 1e3) == 0 ){    
        t_R = 0.0;
    } else {
        t_R = abs(r2 - r1) - t_max;
    }

    /*z0 > z_maxのとき*/
    if(z0 > z_max){
        t_R = 0.0; 
    }

    /*高度が同じ場合*/
    if(r2.r() - r1.r() == 0){
        t_R = 0.0;
     } 
    
    delta3 = opt_vertical(z_max, r2.latitude(), r2.longitude()) * t_R / 1e3;

    /*任意の2点の光学的深さ*/
    double Delta = delta1 + delta2 + delta3;

    return Delta;
}