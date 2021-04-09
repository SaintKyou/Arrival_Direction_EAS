#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include "structs.h"
#include "funcs.h"
#include <sstream>
#include <string>
#include <QFileDialog>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "matrix.h"

//Instalation parameters
const double c = 0.2998; // m/nsec
const double sigma = 3; //nsec  погрешность задержки срабатывания ДС
//const double q_st_max = 250*2.56; // 2.56 - sq m площадь ДС

//double q_st_max[36] = {1000,1000,1259,1000,1000,1000,1000,1000,100,100,100,100,
//                       251.2,316.2,251.2,199.5,1000,1000,1000,1000,1000,1000,1000,1000,
//                       1000,1000,1000,1000,1585,1000,1000,1259,1000,1000,1000,1000}; // m^-2
double q_st_max[36] = {125.893,100,158.489,125.893,158.489,158.489,125.893,100,125.893,158.489,100,
                       125.893,125.893,158.489,199.526,158.489,100,158.489,79.433,100,125.893,125.893,
                       100,79.433,100,125.893,125.893,100,158.489,100,100,158.489,100,100,158.489,199.526};
const double q_add_max = 5000./13.;  //  5000 пКл  не надо делить на 2.56

double x_conf[36] = {-25.36,-37.61,-37.61,-25.36,-25.36,-37.61,-37.61,-25.36,6.48,-6.87,-6.87,6.48,37.37,22.47,
                     22.47,37.37,37.37,22.98,22.98,37.37,6.91,-2.79,-2.79,6.91,-9.54,-20.56,-20.56,-9.54,42.14,26.94,
                     26.94,42.14,-25.36,-37.25,-37.25,-25.36};
double y_conf[36] = {5.89,5.89,-7.32,-7.32,37.34,37.34,24.14,24.14,12.65,12.65,-12.63,-12.63,10.96,10.96,-3.95,-3.95,
                     45.57,45.57,28.48,28.48,-56.16,-56.16,-70.54,-70.54,53.62,64.71,53.62,42.54,-16.80,-16.80,-32.10,
                     -32.10,-25.75,-25.75,-39.55,-39.55};
double z_conf[36] = {-6.68,-6.68,-6.68,-6.68,-6.68,-6.68,-6.68,-6.68,0.00,0.00,0.00,0.00,-14.95,-15.27,-15.35,-15.08,
                     -16.46,-16.18,-16.52,-16.52,-16.15,-15.91,-15.93,-16.10,-17.03,-17.80,-17.50,-16.93,-15.14,-15.46,
                     -15.46,-14.69,-7.33,-7.36,-7.34,-7.28};
double m_delays[54] = {-1.21,-2.89,-2.98,-1.88,-1.90,-0.07,-3.36,0.26,-1.24,3.43,1.96,-1.42,1.75,2.77,3.82,1.01,2.05,1.48,
                        1.22,-2.19,-5.03,-4.03,-6.85,-2.57,-0.35,1.96,0.41,2.14,0.32,-1.32,-0.01,0.22,2.55,0.11,2.7,2.67,
                       -2.76,-0.29,0.45,2.75,4.25,0.89,-2.46,-5.34,-2.64,-2.58,0.09,2.77,1.29,-1.86,-2.22,-3.29,-3.25,-0.19};  // задержки м/у станциями
double m_crlnk[36] = {0};
std::vector<double> mdel; // задержки м/у дс
std::vector<double> crlnk;
//std::vector<short> odel; //оптические задержки
short o_del[9] = {};
std::vector<int> g_size;


extern struct TDateTimeKadr dt;
extern struct SSETUP set;
extern struct SCOORD sc;
extern struct STATION ST;
extern struct DATA DataEAS;
extern struct CLUSTER CL;
extern struct HEADER hdr;
extern struct SCROSSLINK scnk;
extern struct SDELAYS del;
struct SModelEAS
{
                float N_event;                  //[0,100] разбивание по 10 событий
                float NRUN, NEVENT, PART0, E0; // NRUN - имя файла, NEVENT - [0,10], PART0 = 14 id частицы
                float Teta, Fi;
                float XAxisShift, YAxisShift;
                float H1INT;
                float NGAM, NEL, NHADR, NMU;
                float NeNKGlong, sNKGlong;
                float NVD_edep;
                float NVD_npe;
                float MuBundle;
                float MuTrackLenNVD;
                float nMuNVD, eMuNVD, eMuNVD1, nMuDCR, nHitSM, nMuDCRwat, nHitSMwat;
                float AmplKSM[7][4][4][6], EdepCntSCT[9][5][2];
                float EdepDetNE[9][4][4]; //cluster, station, detector
                float TimDetNE[9][4][4][4]; //cluster, station, detector, time (max, th = 8.2 MeV, min, mean) пороговое время
                float EdepStNE[9][4], TimStNE[9][4][4]; //### new cluster, station, time (max, th, min, mean)
                float marker; //useless
};

std::string out_folder{};
std::string in_folder{};
const int M = 5000;
double *E = new double[M];
double *tetha = new double[M];
double *phi = new double[M];
double *s = new double[M];
double *X = new double[M];
double *Y = new double[M];
double *Ne = new double[M];

double *t = new double[M*36];   //50000*36
double *t1 = new double[M*36];
double *Q = new double[M*36];
double *Q_dumb = new double[M*36*4];
double *Q_add = new double[M*36];

long *t_dop = new long[4];
double *t_cl = new double[M*9];
double *t_tr = new double[M*9];
long *t_add = new long[M*36];
long *t_trig = new long[M*36];

//exp
//const int M1 = 4000;
//double *t_exp = new double[M1*36];
//double *t_cl_exp = new double[M1*9];
//double *t1_exp = new double[M1*36];

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->spinBox_krat->setRange(7,9);
    //ui->spinBox_krat->value();
    setFixedSize ( 400, 374 );
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_toolButton_input_clicked()
{
   // if(ui->checkBox_model->isChecked()){        // модельный формат
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),"/home",
                  QFileDialog::ShowDirsOnly|QFileDialog::DontResolveSymlinks);
    QDir directory(dir);
    ui->lineEdit_file->setText(dir);
    QStringList files = directory.entryList(QStringList() << "*.dat" << "*.DAT",QDir::Files);
    int ii{}, jj{}, count{};
    double q_s{}, n_i{};          //суммарный заряд
    //long t_d{};
    //int number_of_sob = M;
    std::ifstream file;
    SModelEAS vg;

    foreach (QString filename, files) {
        std::string current_locale_text = filename.toLocal8Bit().constData();
        std::string curent_folder = dir.toStdString();
        std::fstream f(curent_folder+"\\" + current_locale_text);
        std::string name = curent_folder+"\\" + current_locale_text;
        file.open(name.c_str(), std::fstream::binary | std::fstream::out);
        if (file){
            file.seekg(0); // To make sure that the data is read from the starting position of the file.
            while (file.read((char *)&vg, sizeof(vg))){
            E[ii] = vg.E0;
            phi[ii] = vg.Fi;
            tetha[ii] = vg.Teta;
            s[ii] = vg.sNKGlong;
            X[ii] = vg.XAxisShift;
            Y[ii] = vg.YAxisShift;
            Ne[ii] = vg.NeNKGlong;
            ++ii;
            for (int j = 0; j < 9; j++) {
                for (int m = 0; m < 4; m++) {
                    t[jj] = vg.TimStNE[j][m][1];
                    //++jj;
                    for (int l = 0; l < 4 ; ++l) {
                        q_s += (vg.EdepDetNE[j][m][l]); // Мэв
                        ++count;
                        if (count % 4 == 0)
                        {
                            n_i = q_s/8.2;      //0.75 - эксперимент
//                            if(n_i < 1)       // делал для нейронки
//                            {
//                                n_i = 0;
//                                q_s = 0;
//                            }
                            Q[jj] = n_i;    //n_i {q_s}
                            Q_add[jj] = (vg.EdepDetNE[j][m][3]/8.2)/500.0;   //чтобы 13пКл получилось
                            q_s = 0;
                            ++jj;
                        }
                    }
                }
             }
           }
        }
      }
  //  }

    //qDebug() << ui->lineEdit_out->text();
    std::ofstream fo3("C:\\test\\q_norm.txt");
    int countn = 0;
    double max_q{};
    if(ui->checkBox_norm->isChecked()){
        for (long long ii = 0; ii < 36*M; ++ii) {
            ++countn;
            if(countn%36==0){
                for (int kk = countn-36; kk < countn; ++kk) {
                  if(Q[kk] > max_q) max_q = Q[kk];
                }
                for (int kk = countn-36; kk < countn; ++kk) {
                    if(max_q!=0) Q[kk]/= max_q;
                }
                fo3 << Q[countn-36] << ',' << Q[countn-35] << ',' << Q[countn-34] << ',' << Q[countn-33] << ','
                << Q[countn-32] << ',' << Q[countn-31] << ',' << Q[countn-30] << ',' << Q[countn-29] << ','
                << Q[countn-28] << ',' << Q[countn-27] << ',' << Q[countn-26] << ',' << Q[countn-25] << ','
                << Q[countn-24] << ',' << Q[countn-23] << ',' << Q[countn-22] << ',' << Q[countn-21] << ','
                << Q[countn-20] << ',' << Q[countn-19] << ',' << Q[countn-18] << ',' << Q[countn-17] << ','
                << Q[countn-16] << ',' << Q[countn-15] << ',' << Q[countn-14] << ',' << Q[countn-13] << ','
                << Q[countn-12] << ',' << Q[countn-11] << ',' << Q[countn-10] << ',' << Q[countn-9] << ','
                << Q[countn-8] << ',' << Q[countn-7] << ',' << Q[countn-6] << ',' << Q[countn-5] << ','
                << Q[countn-4] << ',' << Q[countn-3] << ',' << Q[countn-2] << ',' << Q[countn-1] << "\n";
                max_q = 0;
            }
        }
    }

/*    foreach(QString filename, files) {
    std::string current_locale_text = filename.toLocal8Bit().constData();
    std::string curent_folder = dir.toStdString();
    std::fstream f(curent_folder+"//" + current_locale_text);
    std::string name = curent_folder+"//" + current_locale_text;
    FILE *fp = fopen(name.c_str(), "rb");
    SModelEAS *records;
    records = new SModelEAS[number_of_sob];
    fread((SModelEAS *)records, sizeof(SModelEAS), number_of_sob, fp);
    fclose(fp);

    for(int i = 0; i < number_of_sob; i++){
    if (records[i].E0 == 0) break;
    E[ii] = records[i].E0;
    phi[ii] = records[i].Fi;
    tetha[ii] = records[i].Teta;
    s[ii] = records[i].sNKGlong;
    X[ii] = records[i].XAxisShift;
    Y[ii] = records[i].YAxisShift;
    Ne[ii] = records[i].NeNKGlong;
    ++ii;
    for (int j = 0; j < 9; j++) {
        for (int m = 0; m < 4; m++) {
            for (int l = 0; l < 4; l++) {
              t_dop[l] = records[i].TimDetNE[j][m][l][1];
              q_s += (records[i].EdepDetNE[j][m][l]*1.5854);
              ++count;
              if(count%4==0)
              {
                if(t_dop[0]!=-1)t_d = t_dop[0];
                else if(t_dop[1]!=-1)t_d = t_dop[1];
                else if(t_dop[2]!=-1)t_d = t_dop[2];
                else t_d = t_dop[3];
                for (int ll = 0; l < 4; l++) {
                    if(t_dop[ll] < t_d && t_dop[ll]!=-1)t_d = t_dop[ll];
                }
                Q[jj] = q_s;
                t[jj] = records[i].TimStNE[j][m][1];
                Q_add[jj] = (records[i].EdepDetNE[j][m][3]*1.5854)/500.0;   //чтобы 13пКл получилось
                q_s = 0;
              }
            }
            ++jj;
        }
      }
     }
    }*/

//    std::ofstream ft("C:\\test\\new_f\\t.txt");
//    std::ofstream fq("C:\\test\\new_f\\q.txt");
//    for (int i = 0; i < M*36;++i) {
//       ft << t[i] << '\n';
//       fq << Q[i] << '\n';
//    }

//    for (int i = 0; i < M*36; ++i) {
//    if(t[i]>0)
//    {
//        t_add[i]= t[i] - 2320;
//    }
//     else {
//        t_add[i] = -1;
//    }
//  }
}

void MainWindow::on_toolButton_exp_clicked()    // эксперимент
{
    std::ifstream file;
    std::ifstream crl;  // поток для файлов конфигурации установки (crosslink, delays, ssetup)
    std::ofstream input("C:\\Data\\NEAS_DATA\\2020-08-01\\2020-08-01.txt");
    std::ofstream input_size("C:\\Data\\NEAS_DATA\\2020-08-01\\size.txt");
    int number_of_file{};
    std::ofstream out("C:\\test\\exp3.txt");
    std::ofstream out1("C:\\test\\tcl.txt");
    std::ofstream out2("C:\\test\\exp2.txt");
    std::ofstream output("C:\\test\\exp.txt");
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),"C:\\Data\\datas",
                  QFileDialog::ShowDirsOnly|QFileDialog::DontResolveSymlinks);
    QDir directory(dir);
    ui->lineEdit_exp->setText(dir);
    in_folder = dir.toStdString();

    QStringList files = directory.entryList(QStringList() << "*.ne" << "*.NE",QDir::Files);
    int chislo{};
    std::ofstream big_q("C:\\test\\Big_q.txt");

    foreach (QString filename, files) {      
        std::string current_locale_text = filename.toLocal8Bit().constData();
        std::string curent_folder = dir.toStdString();
        std::fstream f(curent_folder + "\\" + current_locale_text);
        std::string name = curent_folder + "\\" + current_locale_text;

        //std::string namee = "C:\\Data\\NEAS_DATA\\2020-08-01\\2020-08-01.cl";
        std::string namee = name;
        std::string namcl, namdel, namset;
        namee.push_back('.');
        for(int i = 0; i < namee.size(); ++i){
            if(namee[i]!=namee[namee.size()-1]) namcl.push_back(name[i]);
            else break;
        }
        namdel = namcl;
        namset = namcl;
        namcl = namcl + ".cl";
        namdel = namdel + ".del";
        namset = namset + ".set";
        ++chislo;

        qDebug() << namcl.c_str() << ' ' << namdel.c_str() << ' ' << name.c_str() << static_cast<double>(chislo)/files.length();
        crl.open(namcl.c_str(), std::fstream::binary | std::fstream::out);
        int k{};
        if(crl){
            crl.seekg(0);
            crl.read((char *)&scnk, sizeof (SCROSSLINK));
            for (int i = 0; i < 9 ; ++i) {
                input << '\n';
                for (int j = 0; j < 4 ; ++j ) {
                    input << scnk.CrossLink[i][j] << ",";
                    crlnk.push_back(scnk.CrossLink[i][j]);
                    m_crlnk[k] = scnk.CrossLink[i][j];
                    ++k;
                }
            }
        }
        else{
            qDebug() << "Error while opening " << namcl.c_str() << ".";
        }
        crl.close();
        k = 0;
        //namee = "C:\\Data\\NEAS_DATA\\2020-08-01\\2020-08-01.del";
        crl.open(namdel.c_str(), std::fstream::binary | std::fstream::out);
        if(crl){
            crl.seekg(0);
            crl.read((char *)&del, sizeof (SDELAYS));
            for (int i = 0; i < 9 ; ++i) {
                input << '\n';
                for (int j = 0; j < 6 ; ++j ) {
                    input << del.Delay[i][j] << ",";
                    mdel.push_back(del.Delay[i][j]);
                    //crlnk[k] = del.Delay[i][j];
                    //++k;
                }
            }
        }
        else{
            qDebug() << "Error while opening " << namdel.c_str() << ".";
        }
        crl.close();

/*        k = 0;
        crl.open(namset.c_str(), std::fstream::binary | std::fstream::out);
        if(crl){
            crl.seekg(0);
            crl.read((char *)&set, sizeof (SSETUP));
            for (int i = 0; i < 9 ; ++i) {
                input << '\n';
                    input << set.OpticDelay[i] << ",";
                    o_del[k] = set.OpticDelay[i];
                    ++k;
                    //odel.push_back(set.OpticDelay[i]);
            }
        }
        else{
            qDebug() << "Error while opening " << namset.c_str() << ".";
        }
        crl.close();
*/

        file.open(name.c_str(), std::fstream::binary | std::fstream::out);
        std::vector<unsigned int> step;
        long long size{}, length{};
        int i{};
        if (file){
            ++number_of_file;
            file.seekg (0, file.end);
            length = file.tellg();
            file.seekg(0);
            while (file.read((char *)&hdr, sizeof (HEADER))) {
                size+=48;
                // out << hdr.lendata << endl;
                //out << hdr.start[0] << hdr.start[1]<< hdr.start[2]<< hdr.start[3]<< hdr.start[4] << endl;
                step.push_back(hdr.lendata);
                size+=hdr.lendata;
            //    out /*<< hdr.id[0] << hdr.id[1]<< hdr.id[2]<< hdr.id[3]<< hdr.id[4]<< hdr.id[5]<< hdr.id[6]<< hdr.id[7]
              //          << hdr.id[8]<< hdr.id[9]<< hdr.id[10]<< hdr.id[11]<< hdr.id[12]<< hdr.id[13]*/<< hdr.id[14]<< hdr.id[15]
                        /*<< hdr.id[16]<< hdr.id[17]<< hdr.id[18]<< hdr.id[19]<< hdr.id[20] << hdr.id[21]<< hdr.id[22]<< hdr.id[23]
                        << hdr.id[24]<< hdr.id[25]<< hdr.id[26]<< hdr.id[27]<< hdr.id[28]<< hdr.id[29]<< hdr.id[30]<< hdr.id[31]
                       << hdr.id[32]<< hdr.id[33]<< hdr.id[34]<< hdr.id[35]<< hdr.id[36]<< hdr.id[37]*/ // <<std::endl;
                //out << hdr.lendata << ' ' << size << ' ' << length << endl;

                out2 << hdr.id[0] << hdr.id[1]<< hdr.id[2]<< hdr.id[3]<< hdr.id[4]<< hdr.id[5]<< hdr.id[6]<< hdr.id[7]
                     << hdr.id[8]<< hdr.id[9]<< hdr.id[10]<< hdr.id[11]<< hdr.id[12]<< hdr.id[13]<< hdr.id[14]<< hdr.id[15]
                     << hdr.id[16]<< hdr.id[17]<< hdr.id[18]<< hdr.id[19]<< hdr.id[20] << hdr.id[21]<< hdr.id[22]<< hdr.id[23]
                     << hdr.id[24]<< hdr.id[25]<< hdr.id[26]<< hdr.id[27]<< hdr.id[28]<< hdr.id[29]<< hdr.id[30]<< hdr.id[31]
                     << hdr.id[32]<< hdr.id[33] /*<< hdr.id[34]<< hdr.id[35]<< hdr.id[36]<< hdr.id[37]*/ << std::endl;

                file.seekg(size);
            }

            file.close();
            size = 48;
            std::vector<unsigned int> v;
            std::vector<unsigned int> v1;
            std::vector<unsigned int> v2;
            std::vector<double> buff;           //buffer for Q
            std::vector<double> buff_add;       //buffer for Q_add
            std::vector<double> buff_t;         //buffer for time (кластерное время)
            std::vector<double> buff_t2;        //время в событии

            long long metka = 48;
            int ff{};

            std::ifstream is (name.c_str(), std::ifstream::binary);
            if (is)
            {
                while(metka < length)
                {
                    is.seekg(metka);
                    is.read ((char *)&DataEAS, step[i]);
                    if(get_dig(dec2bit(DataEAS.Hit)) < 7 /*ui->spinBox_krat->value()*/ || dec2bit(DataEAS.maska)!=0 ){
                        metka+=48;
                        metka+=step[i];
                        ++i;
                        v = make_v(dec2bit(DataEAS.Hit));
                        //out2 << dec2bit(DataEAS.Hit) << ' ' << dec2bit(DataEAS.maska) << std::endl;
                        //for(auto i:v) out2 << i;
                        //out2  << ' ' << get_dig(dec2bit(DataEAS.Hit))  << endl;
//                        for(int i = 0; i < 9; ++i){
//                            if(v[i]==1){
//                               // ++schet[i];
//                            }
//                        }
//                        for(int i = 0; i < 9; ++i){
//                         //   output << schet[i] << '|';
//                        }
                        //output << std::endl;
                        //print_vel(dec2bit(DataEAS.Hit));
                        //out2 << dec2bit(DataEAS.Hit) << " " << get_dig(dec2bit(DataEAS.Hit)) << endl;
                    }
                    else {
                        //cout << "Computing " << ff << " sob " << endl;
                        ++ff;
                        is.seekg(metka);
                        is.read ((char *)&DataEAS, step[i]);
                        v = make_v(dec2bit(DataEAS.Hit));
                        reverse(v.begin(),v.end());
                        for(auto i:v) out1 << i;
                        auto rz = std::accumulate( v.begin(), v.end(), 0, []( int l, int r ) {
                            return l * 10 + r;
                        } );
//                        for(int i = 0; i < 9; ++i){
//                            if(v[i]==1){
//                             //   ++schet[i];
//                            }
//                        }
//                        for(int i = 0; i < 9; ++i){
//                        //    output << schet[i] << '|';
//                        }
                        //output << std::endl;
                        v1.push_back(rz);
                        //v1.push_back(dec2bit(DataEAS.Hit));
                        out1 << " " << get_dig(dec2bit(DataEAS.Hit)) <<  " " << rz  << " " << dec2bit(DataEAS.Hit) << std::endl;
                        //out2 << dec2bit(DataEAS.Hit) << ' ' << dec2bit(DataEAS.maska) << std::endl;
                        //print_vel(dec2bit(DataEAS.Hit));
                        metka+=4;
                        std::vector<long> buff_t1;        // buffer for время события (classic)
                        for (int ii = 0; ii < get_dig(dec2bit(DataEAS.Hit)); ++ii) {
                            is.seekg(metka);
                            is.read ((char *)&CL, 150);
                            is.seekg(metka);
                            is.read ((char *)&CL, 150-(4-get_dig(dec2bit(CL.hit)))*26);
                            v=make_v2(dec2bit(CL.hit));
                            reverse(v.begin(),v.end());
                            for(auto i:v) out1 << i;
                            auto rz = std::accumulate( v.begin(), v.end(), 0, []( int l, int r ) {
                                return l * 10 + r;
                            } );
                            v2.push_back(rz);
                            buff_t1.push_back(CL.time);
                            //output << CL.time << ' ' << get_dig(dec2bit(DataEAS.Hit)) << std::endl;
                            //v2.push_back(dec2bit(CL.hit));
                            out1 << " " << get_dig(dec2bit(CL.hit)) <<  " " << rz << " " << CL.hit << std::endl;
                            metka+=46;
                            if(ii==get_dig(dec2bit(DataEAS.Hit))-1){ //norm_vec(buff_t1);  // нормируем времена кластеров
                            output << '\n';
                            for (int i = 0 ; i < buff_t1.size() ; ++i ){
                            buff_t2.push_back(buff_t1[i]);
                            output << buff_t1[i] << ',';
                                }
                            }
                            for (int jj = 0; jj < get_dig(dec2bit(CL.hit)); ++jj) {
                                is.seekg(metka);
                                is.read ((char *)&ST, 26);
                                //out1 << ST.Q/8.2 << endl;
                                out1 << ST.time << std::endl;
                                //out2 << ST.time << std::endl;
                                buff.push_back(ST.Q/13.); //+ koef sc
                                buff_add.push_back(ST.Q_add/13.);
                                buff_t.push_back(ST.time);
                                metka+=26;
                             // out2 << ST.Q/8.2 << " " << ST.Q_add/8.2 << " " << ST.time << std::endl;
                            }
                        }
                        buff_t1.clear();
                        metka+=48;
                        ++i;
                    }
                }
            }
            std::string outfile = curent_folder+"\\" + "Q" + "\\" + "q" + std::to_string(number_of_file) + ".txt";
            std::ofstream fo1(outfile);
            std::string outfile2 = curent_folder+"\\" + "Q_add" + "\\" + "q_add" + std::to_string(number_of_file) + ".txt";
            std::ofstream fo2(outfile2);
            std::string outfile3 = curent_folder+"\\" + "Q_norm" + "\\" + "q_add" + std::to_string(number_of_file) + ".txt";
            std::ofstream fo3(outfile3);
            //        fo1  << "s1" <<  '\t'  << "s2" <<  '\t'  << "s3" <<  '\t'  << "s4" <<  '\t'  << "s5" <<  '\t'  << "s6" <<  '\t'  << "s7" <<  '\t'  << "s8" <<  '\t'
            //        << "s9" <<  '\t'  << "s10" <<  '\t'  << "s11" <<  '\t'  << "s12" <<  '\t'   << "s13" <<  '\t'  << "s14" <<  '\t'  << "s15" <<  '\t'  << "s16" <<  '\t'
            //        << "s17" <<  '\t'  << "s18" <<  '\t'  << "s19" <<  '\t'  << "s20" <<  '\t'  << "s21" <<  '\t'  << "s22" <<  '\t'  << "s23" <<  '\t'  << "s24" <<  '\t'
            //        << "s25" <<  '\t'  << "s26" <<  '\t'  << "s27" <<  '\t'  << "s28" <<  '\t'  << "s29" <<  '\t'  << "s30" <<  '\t'  << "s31" <<  '\t'  << "s32" <<  '\t'
            //        << "s33" <<  '\t'  << "s34" <<  '\t'  << "s35" <<  '\t'  << "s36" << "\n";

            std::vector <double> Q;
            std::vector <double> Q_add;
            std::vector <double> t;     // cluster time
            std::vector <double> t_dop;
            //std::vector <double> t1;    // classic time
            bool flag{};
            long long m{}, l{}, n{};
            for(int i = 0; i < v1.size(); ++i){
                for(m = 0; m < 9; ++m){
                    for(int j = 0; j < 4; ++j){
                        if(make_v(v1[i])[m]==0){
                            //out << '0' << '\n' << '0' << '\n'<< '0' << '\n'<< '0' << '\n';
                            Q.push_back(0);
                            Q.push_back(0);
                            Q.push_back(0);
                            Q.push_back(0);
                            Q_add.push_back(0);
                            Q_add.push_back(0);
                            Q_add.push_back(0);
                            Q_add.push_back(0);
                            t.push_back(-1);
                            t.push_back(-1);
                            t.push_back(-1);
                            t.push_back(-1);
                            //++m;
                            j+=4;
                            flag = 0;
                        }
                        else {
                            if(make_v2(v2[n])[j]==0){
                                Q.push_back(0);
                                Q_add.push_back(0);
                                t.push_back(-1);
                                //out << "0" << endl;
                            }
                            else{
                                //out << buff[l] << endl;
                                Q.push_back(buff[l]);
                                Q_add.push_back(buff_add[l]);
                                t.push_back(buff_t[l]);
                                ++l;
                            }
                            flag = 1;
                        }
                    }
                    if(flag) ++n;
                    //++m;
                }
            }


            long ll{};
            for(int i = 0; i < v1.size(); ++i){
                for(m = 0; m < 9; ++m){
                   if(make_v(v1[i])[m]==0){
                    t_dop.push_back(-1);
                   }
                   else {
                       t_dop.push_back(buff_t2[ll]);
                       ++ll;
                   }
                   //if(t_dop[i%m]!=-1) t_dop[i%m]-=o_del[m];
                  }
                }

//            for (int i = 0; i < t_dop.size() ; ++i ) {
//                out << t_dop[i] << '\n';
//            }

            long long countn{};
            for (long long ii = 0; ii < Q.size(); ++ii) {
                //        if(t[i] < 0){
                //            Q[i] = 0;
                //        }
                ++countn;
                if(countn%36==0){
                    for (int kk = countn-36; kk < countn ; ++kk) {  //смотреть на ампл. доп канала
                        if(Q[kk] > (q_st_max[kk%36]*2.56) && Q_add[kk]>0){      //не надо делить на 13
                            if(kk%36==12 || kk%36 == 13) Q[kk] = m_crlnk[kk%36]*1.2*Q_add[kk];
                            else if (kk%36==14 || kk%36==25 || kk%36==26) Q[kk] = m_crlnk[kk%36]*1.1*Q_add[kk];
                            else if (kk%36==32) Q[kk] = m_crlnk[kk%36]*0.8*Q_add[kk];
                            else Q[kk] = m_crlnk[kk%36]*Q_add[kk];
                        }
                        out << m_crlnk[kk%36] << ',';
                        if(Q_add[kk] > q_add_max){
                            if(kk%36==12 || kk%36 == 13) Q[kk] = m_crlnk[kk%36]*1.2*q_add_max;
                            else if (kk%36==14 || kk%36==25 || kk%36==26) Q[kk] = m_crlnk[kk%36]*1.1*q_add_max;
                            else if (kk%36==32) Q[kk] = m_crlnk[kk%36]*0.8*q_add_max;
                            else Q[kk] = m_crlnk[kk%36]*q_add_max;
                        }
                      //  if(Q[kk]!=0) out2 << Q[kk] << " " << Q_add[kk] << std::endl;

                    }
                    out << '\n';
                fo1 << Q[countn-36] << ',' << Q[countn-35] << ',' << Q[countn-34] << ',' << Q[countn-33] << ','
                << Q[countn-32] << ',' << Q[countn-31] << ',' << Q[countn-30] << ',' << Q[countn-29] << ','
                << Q[countn-28] << ',' << Q[countn-27] << ',' << Q[countn-26] << ',' << Q[countn-25] << ','
                << Q[countn-24] << ',' << Q[countn-23] << ',' << Q[countn-22] << ',' << Q[countn-21] << ','
                << Q[countn-20] << ',' << Q[countn-19] << ',' << Q[countn-18] << ',' << Q[countn-17] << ','
                << Q[countn-16] << ',' << Q[countn-15] << ',' << Q[countn-14] << ',' << Q[countn-13] << ','
                << Q[countn-12] << ',' << Q[countn-11] << ',' << Q[countn-10] << ',' << Q[countn-9] << ','
                << Q[countn-8] << ',' << Q[countn-7] << ',' << Q[countn-6] << ',' << Q[countn-5] << ','
                << Q[countn-4] << ',' << Q[countn-3] << ',' << Q[countn-2] << ',' << Q[countn-1] << "\n";

                fo2 << Q_add[countn-36] << ',' << Q_add[countn-35] << ',' << Q_add[countn-34] << ',' << Q_add[countn-33] << ','
                << Q_add[countn-32] << ',' << Q_add[countn-31] << ',' << Q_add[countn-30] << ',' << Q_add[countn-29] << ','
                << Q_add[countn-28] << ',' << Q_add[countn-27] << ',' << Q_add[countn-26] << ',' << Q_add[countn-25] << ','
                << Q_add[countn-24] << ',' << Q_add[countn-23] << ',' << Q_add[countn-22] << ',' << Q_add[countn-21] << ','
                << Q_add[countn-20] << ',' << Q_add[countn-19] << ',' << Q_add[countn-18] << ',' << Q_add[countn-17] << ','
                << Q_add[countn-16] << ',' << Q_add[countn-15] << ',' << Q_add[countn-14] << ',' << Q_add[countn-13] << ','
                << Q_add[countn-12] << ',' << Q_add[countn-11] << ',' << Q_add[countn-10] << ',' << Q_add[countn-9] << ','
                << Q_add[countn-8] << ',' << Q_add[countn-7] << ',' << Q_add[countn-6] << ',' << Q_add[countn-5] << ','
                << Q_add[countn-4] << ',' << Q_add[countn-3] << ',' << Q_add[countn-2] << ',' << Q_add[countn-1] << "\n";

                big_q << Q[countn-36] << ',' << Q[countn-35] << ',' << Q[countn-34] << ',' << Q[countn-33] << ','
                << Q[countn-32] << ',' << Q[countn-31] << ',' << Q[countn-30] << ',' << Q[countn-29] << ','
                << Q[countn-28] << ',' << Q[countn-27] << ',' << Q[countn-26] << ',' << Q[countn-25] << ','
                << Q[countn-24] << ',' << Q[countn-23] << ',' << Q[countn-22] << ',' << Q[countn-21] << ','
                << Q[countn-20] << ',' << Q[countn-19] << ',' << Q[countn-18] << ',' << Q[countn-17] << ','
                << Q[countn-16] << ',' << Q[countn-15] << ',' << Q[countn-14] << ',' << Q[countn-13] << ','
                << Q[countn-12] << ',' << Q[countn-11] << ',' << Q[countn-10] << ',' << Q[countn-9] << ','
                << Q[countn-8] << ',' << Q[countn-7] << ',' << Q[countn-6] << ',' << Q[countn-5] << ','
                << Q[countn-4] << ',' << Q[countn-3] << ',' << Q[countn-2] << ',' << Q[countn-1] << "\n";
                }
            }

            countn = 0;
            double max_q{};
            if(ui->checkBox_norm->isChecked()){
                for (long long ii = 0; ii < Q.size(); ++ii) {
                    ++countn;
                    if(countn%36==0){
                        for (int kk = countn-36; kk < countn; ++kk) {
                          if(Q[kk] > max_q) max_q = Q[kk];
                        }
                        for (int kk = countn-36; kk < countn; ++kk) {
                            if(max_q!=0) Q[kk]/= max_q;
                        }
                        fo3 << Q[countn-36] << ',' << Q[countn-35] << ',' << Q[countn-34] << ',' << Q[countn-33] << ','
                        << Q[countn-32] << ',' << Q[countn-31] << ',' << Q[countn-30] << ',' << Q[countn-29] << ','
                        << Q[countn-28] << ',' << Q[countn-27] << ',' << Q[countn-26] << ',' << Q[countn-25] << ','
                        << Q[countn-24] << ',' << Q[countn-23] << ',' << Q[countn-22] << ',' << Q[countn-21] << ','
                        << Q[countn-20] << ',' << Q[countn-19] << ',' << Q[countn-18] << ',' << Q[countn-17] << ','
                        << Q[countn-16] << ',' << Q[countn-15] << ',' << Q[countn-14] << ',' << Q[countn-13] << ','
                        << Q[countn-12] << ',' << Q[countn-11] << ',' << Q[countn-10] << ',' << Q[countn-9] << ','
                        << Q[countn-8] << ',' << Q[countn-7] << ',' << Q[countn-6] << ',' << Q[countn-5] << ','
                        << Q[countn-4] << ',' << Q[countn-3] << ',' << Q[countn-2] << ',' << Q[countn-1] << "\n";
                        max_q = 0;
                    }
                }
            }


            std::string help;
            if(number_of_file < 10) help = "00" + std::to_string(number_of_file);
            else if(number_of_file < 100 and number_of_file >= 10) help = "0" + std::to_string(number_of_file);
            else help = std::to_string(number_of_file);
            std::string outfile_t = curent_folder+ "\\" + "t" + help + ".txt";
//            std::string outfile_t1 = curent_folder+ "\\" + "t_cl" + help + ".txt";       //classic
//            std::ofstream tt1(outfile_t1);
//            for(auto ii:t1) tt1 << ii  << std::endl;
            std::ofstream tt(outfile_t);
            for(auto ii:t) tt << ii  << std::endl;
            input_size << t.size()/36. << '\n';
            g_size.push_back(t.size()/36.);


            Q.clear();
            Q_add.clear();
            t.clear();
            //t1.clear();
            //t_dop.clear();
            buff.clear();
//            buff_t1.clear();
            buff_add.clear();
            v.clear();
            v1.clear();
            v2.clear();
            buff_t.clear();

        file.close();
        }
        else {
            qDebug() << "Error while opening file.";
        }
            step.clear();
            //file.close();
    }
}

void MainWindow::on_toolButton_output_clicked()
{
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),"/home",
                  QFileDialog::ShowDirsOnly|QFileDialog::DontResolveSymlinks);
    QDir directory(dir);
    ui->lineEdit_out->setText(dir);
    out_folder = dir.toStdString();
}

void MainWindow::on_pushButton_clicked()
{
    if(ui->checkBox_mdl->isChecked()){ //model

    this->setStatusTip("Идет преобразование...");
    QString pat = ui->lineEdit_out->text();
    std::string str = pat.toStdString();
    str.c_str();
    std::ofstream fo(str +"\\" + "output_new.csv");
    if(ui->checkBox_write->isChecked()){
    std::ofstream model_data(str + "\\" + "model.csv");
    model_data << "E0" << ", " <<  "phi" << ", " << "tetha" << ", "  << "s" << ", " << "X" << ", " << "Y" << ", " << "Ne" << std::endl;
        for (int i = 0; i < 5000; ++i) {
            model_data << E[i] << ", " <<  phi[i] << ", " << tetha[i] << ", "  << s[i] << ", " << X[i] << ", " << Y[i] << ", " << Ne[i] << std::endl;
        }
    }

    double *x1 = new double[3600];
    double *y1 = new double[3600];
    double *z1 = new double[3600];

    std::ofstream out("C:\\test\\vremya.txt");
    int count{}, count_cl{};
    double min{};
    int kk{}, jj{};

    // Задержки сигнала между ДС
    count = 0;
/*    for (jj = 0; jj < M*36;++jj) {   // вычитать задержку из первого сработавшего ДС
        ++count;
        if(count_cl==9) count_cl = 0;
        if(count%4 ==0)
        {
            if(t[jj-3]<0) t[jj-3] = 40000000.0;
            if(t[jj-2]<0) t[jj-2] = 40000000.0;
            if(t[jj-1]<0) t[jj-1] = 40000000.0;
            if(t[jj]<0)   t[jj] = 40000000.0;

            if(t[jj-3] < t[jj-2] && t[jj-3] < t[jj-1] && t[jj-3] < t[jj])   //если 1 дс самый первый
            {
                t[jj-2] += m_delays[count_cl*6];                            // 1-2
                t[jj-1] += m_delays[count_cl*6+1];                          // 1-3
                t[jj] += m_delays[count_cl*6+2];                            // 1-4
            }
            else if (t[jj-2] < t[jj-3] && t[jj-2] < t[jj-1] && t[jj-2] < t[jj]) // если 2 самый первый
            {
                t[jj-3] += m_delays[count_cl*6];                                // 2-1
                t[jj-1] += m_delays[count_cl*6+3];                              // 2-3
                t[jj] += m_delays[count_cl*6+4];                                // 2-4
            }
            else if (t[jj-1] < t[jj-3] && t[jj-1] < t[jj-2] && t[jj-1] < t[jj]) // если 3 самый первый
            {
                t[jj-3] += m_delays[count_cl*6+1];                              // 3-1
                t[jj-2] += m_delays[count_cl*6+3];                              // 3-2
                t[jj] += m_delays[count_cl*6+5];                                // 3-4
            }
            else if (t[jj] < t[jj-3] && t[jj] < t[jj-2] && t[jj] < t[jj-1]) {   // если 4 самый первый
                t[jj-3] += m_delays[count_cl*6+2];                              // 4-1
                t[jj-2] += m_delays[count_cl*6+4];                              // 4-2
                t[jj-1] += m_delays[count_cl*6+5];                              // 4-3
            }

            if(t[jj-3]>20000000) t[jj-3] = -1;
            if(t[jj-2]>20000000) t[jj-2] = -1;
            if(t[jj-1]>20000000) t[jj-1] = -1;
            if(t[jj]>20000000) t[jj] = -1;

            count = 0;
            ++count_cl;
        }
    }*/

    // Задержки сигнала синхронизации [аддитивная величина на все кластеры 5-10нс]
    for (jj = 0; jj < M*36;++jj) {
        if(t[jj] > 0){             // временные ворота
        //    t[jj]+=5.0;
        }
    }

    //НормировОЧКА
    for (jj = 0; jj < M*36;) {
    ++jj;
    ++kk;
    if(kk==36)
    {
        for (int i = jj-36; i < jj; ++i) {
            if(t[i]!=-1){
                min = t[i];
                break;
            }
        }
            if(min!=0)
            {
                for(int j = jj-36; j < jj; ++j){
                    if(t[j] > 0 && t[j] < min){
                        min = t[j];
                    }
                }
            }
         for (int i = jj-36; i < jj ; ++i) {
               if(t[i]>0) t[i]-=min;
        }
        kk = 0;
    }
  }

    //  Устранение выбросов по времени
        for (jj = 0; jj < M*36;++jj) {
            if(t[jj]>490){             // временные ворота
                  t[jj] = -1;
            }
        }

    count = 0;
    //Срабатывание мин. 2 ДС
    for (jj = 0; jj < M*36;++jj) {
        ++count;
        if(count % 4 == 0)
        {
            if(t[jj-3] < 0 && t[jj-2] < 0 && t[jj-1] < 0 || t[jj-3] < 0 && t[jj-2] < 0 && t[jj] < 0 || t[jj-3] < 0 && t[jj-1] < 0 && t[jj] < 0 || t[jj-2] < 0 && t[jj-1] < 0 && t[jj] < 0 )
            {
                t[jj] = t[jj-1] = t[jj-2] = t[jj-3] = -1;
            }
//            count = 0;
        }
    }

    for (int i = 0; i < M*36; i++) {
         //out << t[i] << '\n';
     }

//    kk = 0;

    count = 0;
    int flag1{};
    std::vector<double> tr_time;
    // Времена кластеров
    jj = 0;
    int pp{};
    for (int i = 0; i < M*36;) {
        ++i;
        ++count;
        if(count==4)
        {
            tr_time.push_back(t[i-4]);
            tr_time.push_back(t[i-3]);
            tr_time.push_back(t[i-2]);
            tr_time.push_back(t[i-1]);
            for (int j = i-4; j < i; ++j) {
                if(t[j]!=-1 && flag1 == 0){
                    min = t[j];
                    ++flag1;
                    break;
                }
            }
            for(int i = 0; i < 4; ++i) if(tr_time[i]==-1) tr_time[i] = 500;
            std::sort(tr_time.begin(), tr_time.begin()+4);
                if(flag1!=0)
                {
                    for(int j = i-4; j < i; ++j){
                        if(t[j] != -1 && t[j] < min){
                            min = t[j];
                        }
                    }
                }
                else if (flag1 == 0) {
                    min = -1;
                }
        t_cl[jj] = min;
        if(tr_time[1] == 500) t_tr[jj] = -1;
        else t_tr[jj] = tr_time[1];
        tr_time.push_back(t[i-4]);
        tr_time.push_back(t[i-3]);
        tr_time.push_back(t[i-2]);
        tr_time.push_back(t[i-1]);
        for(int i = 0; i < 4; ++i){
        if(tr_time[i]!=t_tr[jj]) t_trig[pp] = -1;
        else t_trig[pp] = t_tr[jj];
        ++pp;
        }
        ++jj;
        count = 0;
        min = 0;
        tr_time.clear();
        }
        flag1 = 0;
    }
        std::ofstream out1("C:\\Data\\datas\\l1l.txt");
        for (int i = 0; i < M*36;++i) {
             out1 << t_trig[i] << '\n';
        }

    //Времена для кластерного времени
    count = jj = 0;
    for (int i = 0; i < M*36; ) {
        ++i;
        ++count;
        if(count==4){
            for (int j = i-4; j < i; ++j) {
                if(t[j]>0)
                {
                    t1[j]= t[j] - t_cl[jj];         //раскоментить если надо
                }
                else {
                    t1[j] = t[j];
                }
            }
            ++jj;
            count = 0;
        }
    }

    // Временные ворота для кластеров 130 нс
    for (jj = 0; jj < M*36;++jj) {
        if(t1[jj]>110){             // временные ворота 130 нс -> 110
              t1[jj] = -1;
        }
    }

    // Устранение выбросов для всех по кластерам
    for (jj = 0; jj < M*36;++jj) {
        if(t1[jj]<0 && t[jj]>0){             // временные ворота
            t[jj] = -1;
        }
//        if(Q[jj] < 2){
//            t[jj] = -1;
//        }
    }

    //срабатывание мин 2 (2 раз)
    for (jj = 0; jj < M*36;++jj) {
        ++count;
        if(count % 4 == 0)
        {
            if(t[jj-3] < 0 && t[jj-2] < 0 && t[jj-1] < 0 || t[jj-3] < 0 && t[jj-2] < 0 && t[jj] < 0 || t[jj-3] < 0 && t[jj-1] < 0 && t[jj] < 0 || t[jj-2] < 0 && t[jj-1] < 0 && t[jj] < 0 )
            {
                t[jj] = t[jj-1] = t[jj-2] = t[jj-3] = -1;
            }
//            count = 0;
        }
    }

    std::ofstream out2("C:\\test\\l2l.txt");
    for (int i = 0; i < M*36;++i) {
         out2 << t[i] << '\n';
    }
    int countn{};   // запись количества частиц
    fo  << "s1" << "," << "s2" << "," << "s3" << "," << "s4" << "," << "s5" << "," << "s6" << "," << "s7" << "," << "s8" << ","
    << "s9" << "," << "s10" << "," << "s11" << "," << "s12" << ","  << "s13" << "," << "s14" << "," << "s15" << "," << "s16" << ","
    << "s17" << "," << "s18" << "," << "s19" << "," << "s20" << "," << "s21" << "," << "s22" << "," << "s23" << "," << "s24" << ","
    << "s25" << "," << "s26" << "," << "s27" << "," << "s28" << "," << "s29" << "," << "s30" << "," << "s31" << "," << "s32" << ","
    << "s33" << "," << "s34" << "," << "s35" << "," << "s36" << "\n";
        for (int i = 0; i < 5000*36; ++i) {
            if(t[i] < 0){
                Q[i] = 0;
            }
            ++countn;
            if(countn%36==0){
            fo << Q[countn-36] << "," << Q[countn-35] << "," << Q[countn-34] << "," << Q[countn-33] << ","
            << Q[countn-32] << "," << Q[countn-31] << "," << Q[countn-30] << "," << Q[countn-29] << ","
            << Q[countn-28] << "," << Q[countn-27] << "," << Q[countn-26] << "," << Q[countn-25] << ","
            << Q[countn-24] << "," << Q[countn-23] << "," << Q[countn-22] << "," << Q[countn-21] << ","
            << Q[countn-20] << "," << Q[countn-19] << "," << Q[countn-18] << "," << Q[countn-17] << ","
            << Q[countn-16] << "," << Q[countn-15] << "," << Q[countn-14] << "," << Q[countn-13] << ","
            << Q[countn-12] << "," << Q[countn-11] << "," << Q[countn-10] << "," << Q[countn-9] << ","
            << Q[countn-8] << "," << Q[countn-7] << "," << Q[countn-6] << "," << Q[countn-5] << ","
            << Q[countn-4] << "," << Q[countn-3] << "," << Q[countn-2] << "," << Q[countn-1] << "\n";
                }
        }

    // Расчет начальных приближений
  /*      std::ofstream r(str+"r0.csv");
        r << "X0" << "," << "Y0" << std::endl;
        count = 0;
        double Q_all{}, Summ_x{}, Summ_y{};
        int k = 0;
        double *x0 = new double[M];
        double *y0 = new double[M];
        for (int i = 0; i < M*36; ++i) {
            ++count;
            Summ_x+=Q[i]*x_conf[k];
            Summ_y+=Q[i]*y_conf[k];
            Q_all+=Q[i];
            ++k;
            if(count%36==0)
            {
                x0[i/36] = Summ_x/Q_all;
                y0[i/36] = Summ_y/Q_all;
                r << x0[i/36] << "," << y0[i/36] << std::endl;
                Summ_x = 0;
                Summ_y = 0;
                Q_all = 0;
                k = 0;
            }
        }*/

    // Координата самого красного ДС (или кластера??)
/*        std::ofstream r(str+"rr.csv");
        count = countn = 0;
        int k = 0;
        double qmax{},xx{},yy{};
        double *xr = new double[M];
        double *yr = new double[M];
        r << "xr" << "," << "yr" << std::endl;
        for (int i = 0; i < M*36; ++i) {
            ++count;
            if(Q[i] > qmax)
            {
                qmax = Q[i];
                xx = x_conf[i-countn*36];
                yy = y_conf[i-countn*36];
            }
            if(count%36==0){
                ++countn;
                r << xx << "," << yy << std::endl;
                qmax = 0;
            }
        }*/


    int j{};
    for (int i = 0; i < 3600; i++) {
        x1[i] = x_conf[j];
        y1[i] = y_conf[j];
        z1[i] = z_conf[j];
        if(j==35)j=0;
        else j++;
    }

    if(ui->radioButton_classic->isChecked())        //классический метод
    {        
        std::ofstream vosst_data(str+ "\\" +"classic.csv");
        std::ofstream cl_time(str+ "\\" +"times.csv");
        vosst_data << "tetha_v" << "," << "phi_v" << "," << "krat" << std::endl;
        double phi_v{}, tetha_v{};
        double *min = new double[3600];

        double *x_total = new double[36];
        double *y_total = new double[36];
        double *z_total = new double[36];
        double *t_total = new double[36];

        double x_first{},y_first{};
        double alpha{}, betta{}, C{};

        int p{};
        j = 0;
        unsigned int count_st{}, count_av{}, count_cl{}, u{};

        Matrix M(4,4), Omega(4,1), X(4,1);
        int metka{};
        count = j = 0;
        for (int w = 1; w < 51; w++) {

        int k{};

        // модифицированный классический метод

        for (int i = metka; i < 3600+metka; i++) {
//              if(t1[i]==0){
//            min[j] = t_cl[k];
//            ++k;
//            }
//            else min[j] = -1;
              min[j] = t[i];
      //      min[j] = t_trig[i];
            ++j;
        }
        j = 0;

        for (int i = 0; i < 3600; ++i) {
            count++;

            x_total[p] = x1[i];
            y_total[p] = y1[i];
            z_total[p] = z1[i];
            t_total[p] = min[i];

            if(p == 35){
            cl_time << t_total[0] << "," << t_total[1] << "," << t_total[2] << "," << t_total[3] << ","
            << t_total[4] << "," << t_total[5] << "," << t_total[6] << "," << t_total[7] << ","
            << t_total[8] << "," << t_total[9] << "," << t_total[10] << "," << t_total[11] << ","
            << t_total[12] << "," << t_total[13] << "," << t_total[14] << "," << t_total[15] << ","
            << t_total[16] << "," << t_total[17] << "," << t_total[18] << "," << t_total[19] << ","
            << t_total[20] << "," << t_total[21] << "," << t_total[22] << "," << t_total[23] << ","
            << t_total[24] << "," << t_total[25] << "," << t_total[26] << "," << t_total[27] << ","
            << t_total[28] << "," << t_total[29] << "," << t_total[30] << "," << t_total[31] << ","
            << t_total[32] << "," << t_total[33] << "," << t_total[34] << "," << t_total[35] << "\n";
            p = 0;
            for (int i = 0; i < 36;++i) {
                if(t_total[i]==0){
                    x_first = x_total[i];
                    y_first = y_total[i];
                }
            }
            }
            else p++;

            if(count%36 == 0)
     {
            for (int j = 0; j < 36; j++) {               // заполнение матриц
                ++count_st;
                if(count_st%4==0 && (t_total[j] !=-1 || t_total[j-1] !=-1 || t_total[j-2] !=-1 || t_total[j-3] !=-1 )){count_cl++;}
                if(t_total[j]!=-1){
                ++count_av;
                ++u;

                M.Add_to_el(0,0,x_total[j]*x_total[j]);
                M.Add_to_el(0,1,x_total[j]*y_total[j]);
                M.Add_to_el(0,2,x_total[j]*z_total[j]);
                M.Add_to_el(0,3,x_total[j]);

                M.Add_to_el(1,0,x_total[j]*y_total[j]);
                M.Add_to_el(1,1,y_total[j]*y_total[j]);
                M.Add_to_el(1,2,y_total[j]*z_total[j]);
                M.Add_to_el(1,3,y_total[j]);

                M.Add_to_el(2,0,x_total[j]*z_total[j]);
                M.Add_to_el(2,1,y_total[j]*z_total[j]);
                M.Add_to_el(2,2,z_total[j]*z_total[j]);
                M.Add_to_el(2,3,z_total[j]);

                M.Add_to_el(3,0,x_total[j]);
                M.Add_to_el(3,1,y_total[j]);
                M.Add_to_el(3,2,z_total[j]);
                M.Set_el(3,3,u);

                Omega.Add_to_el(0,0,c*t_total[j]*x_total[j]);
                Omega.Add_to_el(1,0,c*t_total[j]*y_total[j]);
                Omega.Add_to_el(2,0,c*t_total[j]*z_total[j]);
                Omega.Add_to_el(3,0,c*t_total[j]);
                }
            }
              u = 0;

              M.reverse();
              X.mult_arr(M,Omega);

              alpha = X.Get_el(0,0);
              betta = X.Get_el(1,0);
              C = X.Get_el(2,0);

              alpha = -alpha/(sqrt(alpha*alpha+betta*betta+C*C));
              betta = -betta/(sqrt(alpha*alpha+betta*betta+C*C));
              C = -C/(sqrt(alpha*alpha+betta*betta+C*C));

              phi_v = atan2(betta,alpha)*180.0/M_PI;
              tetha_v = acos(C)*180.0/M_PI;
              if(phi_v < 0) phi_v+=360.0;

              vosst_data << tetha_v << "," << phi_v << "," << count_cl  << std::endl;

              count_av = count_st = count_cl = 0;

              for (int i = 0; i < 4; ++i) {
                  Omega.Set_el(i,0,0);
                  for (int j = 0; j < 4; ++j) {
                      M.Set_el(i,j,0);
                  }
              }
        }
      }
        metka+=3600;
    }
 }

    if(ui->radioButton_clust3->isChecked() || ui->radioButton_clust4->isChecked()){     //кластерный метод

        std::ofstream vosst_data(str+"\\" +"clustern.csv");
        std::ofstream angels(str+"\\" +"angels.csv");
        vosst_data << "tetha_v" << "," << "phi_v" << std::endl;
        double *min = new double[3600];
        double *tcl = new double[4];
        double *xcl = new double[4];
        double *ycl = new double[4];
        int metka{}, p{};
        unsigned int count_st{}, count{},u{}, method_ch{}, count_av{};
        double A{}, B{}, alpha{},betta{}, C{}, tethacl{}, phicl{};

        double *phi_m = new double[9*M];
        double *tetha_m = new double[9*M];
        double *A_m = new double[9*M];
        double *B_m = new double[9*M];

        if(ui->radioButton_clust3->isChecked()) method_ch = 1;
        else method_ch = 0;

        Matrix M(3,3), Omega(3,1), X(3,1);
        j = 0;
        for (int w = 1; w < 51; w++) {
        for (int i = metka; i < 3600+metka; i++) {
            min[j] = t1[i];
            ++j;
        }
        j = 0;

        for (int i = 0; i < 3600; ++i) {
            ++count;
            if(min[i]==-1)++count_st;       // счетчик несработавших станций
            if(count%4 == 0)
     {
                tcl[3] = min[i];
                tcl[2] = min[i-1];
                tcl[1] = min[i-2];
                tcl[0] = min[i-3];

                xcl[3] = x1[i];
                xcl[2] = x1[i-1];
                xcl[1] = x1[i-2];
                xcl[0] = x1[i-3];

                ycl[3] = y1[i];
                ycl[2] = y1[i-1];
                ycl[1] = y1[i-2];
                ycl[0] = y1[i-3];

                if(count_st > method_ch){       // если не сработало больше 1 станции угол восстановить нельзя   [для 4 > 0]
                    tetha_m[p] = -1;
                    phi_m[p] = -1;
                    A_m[p] = 0;
                    B_m[p] = 0;
                    ++p;
                    count_st = 0;
                }
                else {
                    for (int j = 0; j < 4; ++j) {               // заполнение матриц
                        if(tcl[j]!=-1){
                            ++u;
                            M.Add_to_el(0,0,xcl[j]*xcl[j]);
                            M.Add_to_el(0,1,xcl[j]*ycl[j]);
                            M.Add_to_el(0,2,xcl[j]);
                            M.Add_to_el(1,1,ycl[j]*ycl[j]);
                            M.Add_to_el(1,2,ycl[j]);
                            M.Add_to_el(1,0,xcl[j]*ycl[j]);
                            M.Add_to_el(2,0,xcl[j]);
                            M.Add_to_el(2,1,ycl[j]);
                            M.Set_el(2,2,u);
                            Omega.Add_to_el(0,0,c*tcl[j]*xcl[j]);
                            Omega.Add_to_el(1,0,c*tcl[j]*ycl[j]);
                            Omega.Add_to_el(2,0,c*tcl[j]);
                        }
                    }
                    u = 0;

                    M.reverse();
                    X.mult_arr(M,Omega);

                    alpha = X.Get_el(0,0);
                    betta = X.Get_el(1,0);
                    C = sqrt(1-alpha*alpha-betta*betta);

                    if(C < 0){
                        A_m[p] = 0;
                        B_m[p] = 0;
                    }
                    else{
                    alpha = alpha*(-1);
                    betta = betta*(-1);
                    A_m[p] = alpha;
                    B_m[p] = betta;
                    ++count_av;
                    }

                    ++p;
                    count_st = 0;
                }
                for (int i = 0; i < 3; ++i) {
                    Omega.Set_el(i,0,0);
                    for (int j = 0; j < 3; ++j) {
                        M.Set_el(i,j,0);
                    }
                }
        }
        }
        metka+=3600;
    }

        //обработка полученных углов (кластерный метод)

        unsigned int count_av_phi{};
        count = 0;

        for (int i = 0; i < 900*50; ++i) {
            ++count;                                        //счетчик для обработки кластера
            if(A_m[i]!= 0 && A_m[i] > -1.0 && A_m[i] < 1.0 && B_m[i]!= 0 && B_m[i] < 1.0 && B_m[i] > -1.0){
            A+=A_m[i];
            B+=B_m[i];
//          test_data << acos(1-A_m[i]*A_m[i] - B_m[i]*B_m[i])*180/M_PI << ' ' << atan2(B_m[i],A_m[i])*180/M_PI << std::endl;
            ++count_av_phi;
            }
            else {
                A_m[i] = 0;
                B_m[i] = 0;
            }
            if (count % 9 == 0){
                    A = A/count_av_phi;
                    B = B/count_av_phi;
                    C = sqrt(1-A*A-B*B);
            angels << acos(sqrt(1-A_m[count-9]*A_m[count-9]-B_m[count-9]*B_m[count-9]))*180/M_PI << "," << acos(sqrt(1-A_m[count-8]*A_m[count-8]-B_m[count-8]*B_m[count-8]))*180/M_PI << ","
            << acos(sqrt(1-A_m[count-7]*A_m[count-7]-B_m[count-7]*B_m[count-7]))*180/M_PI << "," << acos(sqrt(1-A_m[count-6]*A_m[count-6]-B_m[count-6]*B_m[count-6]))*180/M_PI << ","
            << acos(sqrt(1-A_m[count-5]*A_m[count-5]-B_m[count-5]*B_m[count-5]))*180/M_PI << "," << acos(sqrt(1-A_m[count-4]*A_m[count-4]-B_m[count-4]*B_m[count-4]))*180/M_PI << ","
            << acos(sqrt(1-A_m[count-3]*A_m[count-3]-B_m[count-3]*B_m[count-3]))*180/M_PI << "," << acos(sqrt(1-A_m[count-2]*A_m[count-2]-B_m[count-2]*B_m[count-2]))*180/M_PI << ","
            << acos(sqrt(1-A_m[count-1]*A_m[count-1]-B_m[count-1]*B_m[count-1]))*180/M_PI << '\n'; // tetha
            angels << atan2(B_m[count-9],A_m[count-9])*180/M_PI << "," << atan2(B_m[count-8],A_m[count-8])*180/M_PI << ","
            << atan2(B_m[count-7],A_m[count-7])*180/M_PI << "," << atan2(B_m[count-6],A_m[count-6])*180/M_PI << ","
            << atan2(B_m[count-5],A_m[count-5])*180/M_PI << "," << atan2(B_m[count-4],A_m[count-4])*180/M_PI << ","
            << atan2(B_m[count-3],A_m[count-3])*180/M_PI << "," << atan2(B_m[count-2],A_m[count-2])*180/M_PI << ","
            << atan2(B_m[count-1],A_m[count-1])*180/M_PI << '\n'; // phi
                    tethacl = acos(C)*180/M_PI;
                    phicl = atan2(B,A)*180/M_PI;
                    if(phicl < 0) phicl+=360.0;

                vosst_data << tethacl << "," << phicl << "," << count_av_phi << std::endl;
                A = B = count_av = count_av_phi = 0;
            }
        }
       }
    }

    if(ui->checkBox_exp->isChecked()){   // exp
        std::ifstream input_size("C:\\Data\\NEAS_DATA\\2020-08-01\\size.txt");
        QString dir = "C:\\Data\\dutas";
        std::ifstream file;
        QDir directory(dir);
        QStringList files = directory.entryList(QStringList() << "*.txt" << "*.TXT",QDir::Files);
        long chislo{};
        //int *in_size = new int[1];
        int in_size{};
        std::ofstream big_a("C:\\test\\Big_a.txt");

        foreach (QString filename, files) {
            std::string current_locale_text = filename.toLocal8Bit().constData();
            std::string curent_folder = dir.toStdString();
            std::fstream f(curent_folder + "\\" + current_locale_text);
            std::string name = curent_folder + "\\" + current_locale_text;
            file.open(name.c_str(), std::fstream::binary | std::fstream::out);
            if(file){

//          for (i = chislo ; i < g_size.size() ; ++i) {
//              in_size = g_size[i];
//              if(in_size[i]<=0)break;
//              //}
                //int M1 = in_size[0];

                in_size = g_size[chislo];

                int M1 = 5500;  //15500 - 5
                double *t_exp = new double[M1*36];
                double *t_cl_exp = new double[M1*9];
                double *t1_exp = new double[M1*36];

                qDebug() << "Opened successfuly: " << name.c_str() << " Events: " << std::to_string(in_size).c_str();
                ++chislo;
                for (int i = 0; i < M1*36; ++i) {
                    file >> t_exp[i];
                }

                for (int i = 0; i < 54; ++i) {
                    if(i < 36) m_crlnk[i] = crlnk[i+36*(chislo-1)];     // нигде не используется
                    m_delays[i] = mdel[i+54*(chislo-1)];
                }

                int count{}, count_cl{};
                int kk{}, jj{};
                double min{};

                // Задержки сигнала между ДС
                count = 0;
                for (jj = 0; jj < M1*36;++jj) {   // вычитать задержку из первого сработавшего ДС
                    ++count;
                    if(count_cl==9) count_cl = 0;
                    if(count%4 ==0)
                    {
                        if(t_exp[jj-3]<0) t_exp[jj-3] = 40000000.0;
                        if(t_exp[jj-2]<0) t_exp[jj-2] = 40000000.0;
                        if(t_exp[jj-1]<0) t_exp[jj-1] = 40000000.0;
                        if(t_exp[jj]<0)   t_exp[jj] = 40000000.0;

                        if(t_exp[jj-3] < t_exp[jj-2] && t_exp[jj-3] < t_exp[jj-1] && t_exp[jj-3] < t_exp[jj])   //если 1 дс самый первый
                        {
                            t_exp[jj-2] -= m_delays[count_cl*6];                            // 1-2
                            t_exp[jj-1] -= m_delays[count_cl*6+1];                          // 1-3
                            t_exp[jj] -= m_delays[count_cl*6+2];                            // 1-4
                        }
                        else if (t_exp[jj-2] < t_exp[jj-3] && t_exp[jj-2] < t_exp[jj-1] && t_exp[jj-2] < t_exp[jj]) // если 2 самый первый
                        {
                            t_exp[jj-3] -= m_delays[count_cl*6];                                // 2-1
                            t_exp[jj-1] -= m_delays[count_cl*6+3];                              // 2-3
                            t_exp[jj] -= m_delays[count_cl*6+4];                                // 2-4
                        }
                        else if (t_exp[jj-1] < t_exp[jj-3] && t_exp[jj-1] < t_exp[jj-2] && t_exp[jj-1] < t_exp[jj]) // если 3 самый первый
                        {
                            t_exp[jj-3] -= m_delays[count_cl*6+1];                              // 3-1
                            t_exp[jj-2] -= m_delays[count_cl*6+3];                              // 3-2
                            t_exp[jj] -= m_delays[count_cl*6+5];                                // 3-4
                        }
                        else if (t_exp[jj] < t_exp[jj-3] && t_exp[jj] < t_exp[jj-2] && t_exp[jj] < t_exp[jj-1]) {   // если 4 самый первый
                            t_exp[jj-3] -= m_delays[count_cl*6+2];                              // 4-1
                            t_exp[jj-2] -= m_delays[count_cl*6+4];                              // 4-2
                            t_exp[jj-1] -= m_delays[count_cl*6+5];                              // 4-3
                        }

                        if(t_exp[jj-3]>20000000) t_exp[jj-3] = -1;
                        if(t_exp[jj-2]>20000000) t_exp[jj-2] = -1;
                        if(t_exp[jj-1]>20000000) t_exp[jj-1] = -1;
                        if(t_exp[jj]>20000000) t_exp[jj] = -1;

                        count = 0;
                        ++count_cl;
                    }
                }

                //НормировОЧКА
                for (jj = 0; jj < M1*36;) {
                    ++jj;
                    ++kk;
                    if(kk==36)
                    {
                        for (int i = jj-36; i < jj; ++i) {
                            if(t_exp[i]!=-1){
                                min = t_exp[i];
                                break;
                            }
                        }
                        if(min!=0)
                        {
                            for(int j = jj-36; j < jj; ++j){
                                if(t_exp[j] > 0 && t_exp[j] < min){
                                    min = t_exp[j];
                                }
                            }
                        }
                        for (int i = jj-36; i < jj ; ++i) {
                            if(t_exp[i]>0) t_exp[i]-=min;
                        }
                        kk = 0;
                    }
                }

                //  Устранение выбросов по времени
                    for (jj = 0; jj < M1*36;++jj) {
                        if(t_exp[jj]>490){             // временные ворота
                              t_exp[jj] = -1;
                        }
                    }

                    count = 0;
                    //Срабатывание мин. 2 ДС
                    for (jj = 0; jj < M1*36;++jj) {
                        ++count;
                        if(count % 4 == 0)
                        {
                            if(t_exp[jj-3] < 0 && t_exp[jj-2] < 0 && t_exp[jj-1] < 0 || t_exp[jj-3] < 0 && t_exp[jj-2] < 0 && t_exp[jj] < 0 || t_exp[jj-3] < 0 && t_exp[jj-1] < 0 && t_exp[jj] < 0 || t_exp[jj-2] < 0 && t_exp[jj-1] < 0 && t_exp[jj] < 0 )
                            {
                                t_exp[jj] = t_exp[jj-1] = t_exp[jj-2] = t_exp[jj-3] = -1;
                            }
                //            count = 0;
                        }
                    }
                    count = 0;

                    int flag1{};
                    // Времена кластеров
                    jj = 0;
                    for (int i = 0; i < M1*36;) {
                        ++i;
                        ++count;
                        if(count==4)
                        {
                            for (int j = i-4; j < i; ++j) {
                                if(t_exp[j]!=-1 && flag1 == 0){
                                    min = t_exp[j];
                                    ++flag1;
                                    break;
                                }
                            }
                                if(flag1!=0)
                                {
                                    for(int j = i-4; j < i; ++j){
                                        if(t_exp[j] != -1 && t_exp[j] < min){
                                            min = t_exp[j];
                                        }
                                    }
                                }
                                else if (flag1 == 0) {
                                    min = -1;
                                }
                        t_cl_exp[jj] = min;
                        ++jj;
                        count = 0;
                        min = 0;
                        }
                        flag1 = 0;
                    }

                    //Времена для кластерного времени
                    count = jj = 0;
                    for (int i = 0; i < M1*36; ) {
                        ++i;
                        ++count;
                        if(count==4){
                            for (int j = i-4; j < i; ++j) {
                                if(t_exp[j]>0)
                                {
                                    t1_exp[j]= t_exp[j] - t_cl_exp[jj];         //раскоментить если надо
                                }
                                else {
                                    t1_exp[j] = t_exp[j];
                                }
                            }
                            ++jj;
                            count = 0;
                        }
                    }

                    // Временные ворота для кластеров 130 нс ?? 120 нс
                    for (jj = 0; jj < M1*36;++jj) {
                        if(t1_exp[jj]>110){             // временные ворота 130 нс -> 110
                              t1_exp[jj] = -1;
                        }
                    }

                    // Устранение выбросов для всех по кластерам
                    for (jj = 0; jj < M1*36;++jj) {
                        if(t1_exp[jj]<0 && t_exp[jj]>0){             // временные ворота
                            t_exp[jj] = -1;
                        }
                    }

                    //срабатывание мин 2 (2 раз)
                    for (jj = 0; jj < M1*36;++jj) {
                        ++count;
                        if(count % 4 == 0)
                        {
                            if(t_exp[jj-3] < 0 && t_exp[jj-2] < 0 && t_exp[jj-1] < 0 || t_exp[jj-3] < 0 && t_exp[jj-2] < 0 && t_exp[jj] < 0 || t_exp[jj-3] < 0 && t_exp[jj-1] < 0 && t_exp[jj] < 0 || t_exp[jj-2] < 0 && t_exp[jj-1] < 0 && t_exp[jj] < 0 )
                            {
                                t_exp[jj] = t_exp[jj-1] = t_exp[jj-2] = t_exp[jj-3] = -1;
                            }
                //            count = 0;
                        }
                    }

                    std::ofstream f2("C:\\test\\t_expa.txt");
                    for (int i = 0; i < M1*36; ++i) {
                        f2 << t_exp[i] << std::endl;
                    }

                    double *x1 = new double[3600];
                    double *y1 = new double[3600];
                    double *z1 = new double[3600];

                    int j{};
                    for (int i = 0; i < 3600; i++) {
                        x1[i] = x_conf[j];
                        y1[i] = y_conf[j];
                        z1[i] = z_conf[j];
                        if(j==35)j=0;
                        else j++;
                    }


                    if(ui->radioButton_classic->isChecked()) {       //классический метод
                    std::ofstream vosst_data("C:\\test\\classic.txt");
                    vosst_data << "tetha_v" << '\t' << "phi_v" << '\t' << "krat" << std::endl;
                    double phi_v{}, tetha_v{};
                    double *min1 = new double[3600];

                    double *x_total = new double[36];
                    double *y_total = new double[36];
                    double *z_total = new double[36];
                    double *t_total = new double[36];

                    double x_first{},y_first{};
                    double alpha{}, betta{}, C{};

                    int p{};
                    j = 0;
                    unsigned int count_st{}, count_av{}, count_cl1{}, u{};

                    Matrix M(4,4), Omega(4,1), X(4,1);
                    int metka{};
                    count = j = 0;
                    for (int w = 1; w < M1/100; ++w) {

                    for (int i = metka; i < 3600+metka; ++i) {
                        min1[j] = t_exp[i];
                        ++j;
                    }
                    j = 0;

                    for (int i = 0; i < 3600; ++i) {
                        count++;

                        x_total[p] = x1[i];
                        y_total[p] = y1[i];
                        z_total[p] = z1[i];
                        t_total[p] = min1[i];

                        if(p == 35){
                        p = 0;
                        for (int i = 0; i < 36;++i) {
                            if(t_total[i]==0){
                                x_first = x_total[i];
                                y_first = y_total[i];
                            }
                        }
                        }
                        else p++;

                        if(count%36 == 0)
                 {
                        for (int j = 0; j < 36; j++) {               // заполнение матриц
                            ++count_st;
                            if(count_st%4==0 && (t_total[j] !=-1 || t_total[j-1] !=-1 || t_total[j-2] !=-1 || t_total[j-3] !=-1 )){count_cl1++;}
                            if(t_total[j]!=-1){
                            ++count_av;
                            ++u;

                            M.Add_to_el(0,0,x_total[j]*x_total[j]);
                            M.Add_to_el(0,1,x_total[j]*y_total[j]);
                            M.Add_to_el(0,2,x_total[j]*z_total[j]);
                            M.Add_to_el(0,3,x_total[j]);

                            M.Add_to_el(1,0,x_total[j]*y_total[j]);
                            M.Add_to_el(1,1,y_total[j]*y_total[j]);
                            M.Add_to_el(1,2,y_total[j]*z_total[j]);
                            M.Add_to_el(1,3,y_total[j]);

                            M.Add_to_el(2,0,x_total[j]*z_total[j]);
                            M.Add_to_el(2,1,y_total[j]*z_total[j]);
                            M.Add_to_el(2,2,z_total[j]*z_total[j]);
                            M.Add_to_el(2,3,z_total[j]);

                            M.Add_to_el(3,0,x_total[j]);
                            M.Add_to_el(3,1,y_total[j]);
                            M.Add_to_el(3,2,z_total[j]);
                            M.Set_el(3,3,u);

                            Omega.Add_to_el(0,0,c*t_total[j]*x_total[j]);
                            Omega.Add_to_el(1,0,c*t_total[j]*y_total[j]);
                            Omega.Add_to_el(2,0,c*t_total[j]*z_total[j]);
                            Omega.Add_to_el(3,0,c*t_total[j]);
                            }
                        }
                          u = 0;

                          M.reverse();
                          X.mult_arr(M,Omega);

                          alpha = X.Get_el(0,0);
                          betta = X.Get_el(1,0);
                          C = X.Get_el(2,0);

                          alpha = -alpha/(sqrt(alpha*alpha+betta*betta+C*C));
                          betta = -betta/(sqrt(alpha*alpha+betta*betta+C*C));
                          C = -C/(sqrt(alpha*alpha+betta*betta+C*C));

                          phi_v = atan2(betta,alpha)*180.0/M_PI;
                          tetha_v = acos(C)*180.0/M_PI;
                          if(phi_v < 0) phi_v+=360.0;

                          vosst_data << tetha_v << '\t' << phi_v << '\t' << count_cl1 << std::endl;

                          count_av = count_st = count_cl1 = 0;

                          for (int i = 0; i < 4; ++i) {
                              Omega.Set_el(i,0,0);
                              for (int j = 0; j < 4; ++j) {
                                  M.Set_el(i,j,0);
                              }
                          }
                    }
                  }
                    metka+=3600;
                }
             }

                    if(ui->radioButton_clust3->isChecked() || ui->radioButton_clust4->isChecked()){
                    std::ofstream vosst_data(curent_folder + "\\" + "angels\\" + std::to_string(chislo) + "clustern.txt");
                    //std::ofstream vosst_data1("C:\\test\\classic.txt");
                    //std::ofstream vosst_data2("C:\\test\\classic2.txt");
                    //std::ofstream angels("C:\\test\\classic.xls");
                    //vosst_data << "tetha_v" << '\t' << "phi_v" << std::endl;
                    double *min1 = new double[3600];
                    double *tcl = new double[4];
                    double *xcl = new double[4];
                    double *ycl = new double[4];
                    int metka{}, p{};
                    unsigned int count_st{}, count1{},u{}, method_ch{}, count_av{};
                    double A{}, B{}, alpha{},betta{}, C{}, tethacl{}, phicl{};

                    double *phi_m = new double[9*M1]; //2900
                    double *tetha_m = new double[9*M1];
                    double *A_m = new double[9*M1];
                    double *B_m = new double[9*M1];

                    if(ui->radioButton_clust3->isChecked()) method_ch = 1;
                    else method_ch = 0;

                    Matrix M(3,3), Omega(3,1), X(3,1);
                    j = 0;
                    for (int w = 1; w < M1/100; ++w) {              //n_sob = w*100
                    for (int i = metka; i < 3600+metka; ++i) {
                        min1[j] = t1_exp[i];
                        ++j;
                    }
                    j = 0;

                    for (int i = 0; i < 3600; ++i) {
                        ++count1;
                        if(min1[i]==-1)++count_st;       // счетчик несработавших станций
                        if(count1%4 == 0)
                 {
                            tcl[3] = min1[i];
                            tcl[2] = min1[i-1];
                            tcl[1] = min1[i-2];
                            tcl[0] = min1[i-3];

                            xcl[3] = x1[i];
                            xcl[2] = x1[i-1];
                            xcl[1] = x1[i-2];
                            xcl[0] = x1[i-3];

                            ycl[3] = y1[i];
                            ycl[2] = y1[i-1];
                            ycl[1] = y1[i-2];
                            ycl[0] = y1[i-3];

                            if(count_st > method_ch){       // если не сработало больше 1 станции угол восстановить нельзя   [для 4 > 0]
                                tetha_m[p] = -1;
                                phi_m[p] = -1;
                                A_m[p] = 0;
                                B_m[p] = 0;
                                ++p;
                                count_st = 0;
                            }
                            else {
                                for (int j = 0; j < 4; ++j) {               // заполнение матриц
                                    if(tcl[j]!=-1){
                                        ++u;
                                        M.Add_to_el(0,0,xcl[j]*xcl[j]);
                                        M.Add_to_el(0,1,xcl[j]*ycl[j]);
                                        M.Add_to_el(0,2,xcl[j]);
                                        M.Add_to_el(1,1,ycl[j]*ycl[j]);
                                        M.Add_to_el(1,2,ycl[j]);
                                        M.Add_to_el(1,0,xcl[j]*ycl[j]);
                                        M.Add_to_el(2,0,xcl[j]);
                                        M.Add_to_el(2,1,ycl[j]);
                                        M.Set_el(2,2,u);
                                        Omega.Add_to_el(0,0,c*tcl[j]*xcl[j]);
                                        Omega.Add_to_el(1,0,c*tcl[j]*ycl[j]);
                                        Omega.Add_to_el(2,0,c*tcl[j]);
                                    }
                                }
                                u = 0;

                                M.reverse();
                                X.mult_arr(M,Omega);

                                alpha = X.Get_el(0,0);
                                betta = X.Get_el(1,0);
                                C = sqrt(1-alpha*alpha-betta*betta);

                                if(C < 0){
                                    A_m[p] = 0;
                                    B_m[p] = 0;
                                }
                                else{
                                alpha = alpha*(-1);
                                betta = betta*(-1);
                                A_m[p] = alpha;
                                B_m[p] = betta;
                                ++count_av;
                                }

                                ++p;
                                count_st = 0;
                            }
                            for (int i = 0; i < 3; ++i) {
                                Omega.Set_el(i,0,0);
                                for (int j = 0; j < 3; ++j) {
                                    M.Set_el(i,j,0);
                                }
                            }
                    }
                    }
                    metka+=3600;
                }

                    //обработка полученных углов (кластерный метод)

                    unsigned int count_av_phi{};
                    count1 = 0;

                    for (int i = 0; i < 9*M1; ++i) {
                        ++count1;                                        //счетчик для обработки кластера
                        if(A_m[i]!= 0 && A_m[i] > -1.0 && A_m[i] < 1.0 && B_m[i]!= 0 && B_m[i] < 1.0 && B_m[i] > -1.0){
                        A+=A_m[i];
                        B+=B_m[i];
            //          test_data << acos(1-A_m[i]*A_m[i] - B_m[i]*B_m[i])*180/M_PI << ' ' << atan2(B_m[i],A_m[i])*180/M_PI << std::endl;
                        ++count_av_phi;
                        }
                        else {
                            A_m[i] = 0;
                            B_m[i] = 0;
                        }
                        if (count1 % 9 == 0){
                                A = A/count_av_phi;
                                B = B/count_av_phi;
                                C = sqrt(1-A*A-B*B);
        //                angels << acos(sqrt(1-A_m[count1-9]*A_m[count1-9]-B_m[count1-9]*B_m[count1-9]))*180/M_PI << '\t' << acos(sqrt(1-A_m[count1-8]*A_m[count1-8]-B_m[count1-8]*B_m[count1-8]))*180/M_PI << '\t'
        //                << acos(sqrt(1-A_m[count1-7]*A_m[count1-7]-B_m[count1-7]*B_m[count1-7]))*180/M_PI << '\t' << acos(sqrt(1-A_m[count1-6]*A_m[count1-6]-B_m[count1-6]*B_m[count1-6]))*180/M_PI << '\t'
        //                << acos(sqrt(1-A_m[count1-5]*A_m[count1-5]-B_m[count1-5]*B_m[count1-5]))*180/M_PI << '\t' << acos(sqrt(1-A_m[count1-4]*A_m[count1-4]-B_m[count1-4]*B_m[count1-4]))*180/M_PI << '\t'
        //                << acos(sqrt(1-A_m[count1-3]*A_m[count1-3]-B_m[count1-3]*B_m[count1-3]))*180/M_PI << '\t' << acos(sqrt(1-A_m[count1-2]*A_m[count1-2]-B_m[count1-2]*B_m[count1-2]))*180/M_PI << '\t'
        //                << acos(sqrt(1-A_m[count1-1]*A_m[count1-1]-B_m[count1-1]*B_m[count1-1]))*180/M_PI << '\n'; // tetha
        //                angels << atan2(B_m[count1-9],A_m[count1-9])*180/M_PI << "," << atan2(B_m[count1-8],A_m[count1-8])*180/M_PI << '\t'
        //                << atan2(B_m[count1-7],A_m[count1-7])*180/M_PI << '\t' << atan2(B_m[count1-6],A_m[count1-6])*180/M_PI << '\t'
        //                << atan2(B_m[count1-5],A_m[count1-5])*180/M_PI << '\t' << atan2(B_m[count1-4],A_m[count1-4])*180/M_PI << '\t'
        //                << atan2(B_m[count1-3],A_m[count1-3])*180/M_PI << '\t' << atan2(B_m[count1-2],A_m[count1-2])*180/M_PI << '\t'
        //                << atan2(B_m[count1-1],A_m[count1-1])*180/M_PI << '\n'; // phi
                                tethacl = acos(C)*180/M_PI;
                                phicl = atan2(B,A)*180/M_PI;
                                if(phicl < 0) phicl+=360.0;

                            if(i<=in_size*9){
                            vosst_data << tethacl << ',' << phicl << std::endl;
                            big_a << tethacl << ',' << phicl << std::endl;
                            }
                            //vosst_data1 << tethacl << std::endl;
                            //vosst_data2 << phicl << std::endl;
                            A = B = count_av = count_av_phi = 0;
                        }
                    }
                }

            file.close();
            }

            else{
                qDebug() << "Error while opening " << name.c_str();
            }
        }


    }

    this->setStatusTip("Готово.");
}







