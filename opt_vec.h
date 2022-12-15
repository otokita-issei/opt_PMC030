/*固定するパラメタ*/
constexpr double MtoKM {1.0e-3};    /*mをkmに変換*/
constexpr double KMtoM {1.0e3};     /*kmをmに変換*/

/*固定の値*/
constexpr double z_km {100};    /*高度*/
constexpr double t0 = 100;      /*基準の日にち*/
constexpr double L0 = 40;       /*基準の緯度*/
constexpr double inc_t = 10;    /*日にちの増加量*/
constexpr double inc_L = 10;    /*緯度の増加量*/

/*任意の2点間の光学的深さ*/
double opt_vec(AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1);

/*垂直に伝搬する2点間の光学的深さ*/
double opt_vertical(int z, double Lat, double Lon);






/*変化させるパラメタ*/
/*1つめの任意の点*/
// constexpr double Lat1 {42};
// constexpr double Lon1 {0};
// constexpr double Alt1 {54.9e3}; /*高度[m]*/

// /*2つめの任意の点*/
// constexpr double Lat2 {85};
// constexpr double Lon2 {0};
// constexpr double Alt2 {60.1e3}; /*Alt1 < Alt2*/

//constexpr double day_t {172};     /*任意の日にち*/

/*固定するパラメタ*/

//constexpr double LONGITUDE_HIMAWARI {140.7}; /*気象衛生ひまわりの経度*/
//constexpr double Deg2Rad {M_PI/180.0};  /*DegをRadに変換*/
//constexpr double Rad2Deg {180.0/M_PI};  /*RadをDegに変換*/
// constexpr double MtoKM {1.0e-3};    /*mをkmに変換*/
// constexpr double KMtoM {1.0e3};     /*kmをmに変換*/
//constexpr double R0 {6370.0e3};     /*地球の半径[m]*/

/*1つめの任意の点の球座標*/
// constexpr double r_1 {R0 + Alt1};
// constexpr double th1 {(90 - Lat1) * Deg2Rad};
// constexpr double phi1 {(Lon1 - LONGITUDE_HIMAWARI) * Deg2Rad};

/*2つめの任意の点の球座標*/
// constexpr double r_2 {R0 + Alt2};
// constexpr double th2 {(90 - Lat2) * Deg2Rad};
// constexpr double phi2 {(Lon2 - LONGITUDE_HIMAWARI) * Deg2Rad};

//constexpr int z0 = Alt1*MtoKM + 1;
//constexpr int z_max = Alt2*MtoKM;
//constexpr int N_z = z_max - z0 + 1;

// /*opt_vertical.cpp*/
// double opt_vertical(int z, double Lat, double Lon); /*opt_vertical(int z, double Lat, double Lon, Date date)*/

// /*固定の値*/
// constexpr double z_km {100};    /*高度*/
// constexpr double t0 = 100;      /*基準の日にち*/
// constexpr double L0 = 40;       /*基準の緯度*/
// constexpr double inc_t = 10;    /*日にちの増加量*/
// constexpr double inc_L = 10;    /*緯度の増加量*/

// constexpr double t = day_t;
// constexpr int t_n = day_t/10; 
// constexpr int tn = t_n*10;  /*切り捨てたday_t*/
// constexpr int tn_1 = tn + 10;
