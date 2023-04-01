// Abkowitz-RK4.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <stdlib.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;


using namespace std;
void write(string filename, vector<vector<double> >& data)
{
    //filename是待写入的文件名，data是保存待写文件的数据，保证filename一定存在
    ofstream outFile;
    outFile.open(filename, ios::out); //打开这个文件
    for (auto d : data) //遍历待写入的数据
    {
        for (int i = 0; i < d.size(); i++)
        {
            if (i == d.size() - 1) outFile << d[i] << endl;
            else outFile << d[i] << ",";
        }
    }
    cout << "OK" << endl;
}

double DegToRad(double deg)
{
    return deg * 3.1415926 / 180;
}

class Abkowitz
{
    double L = 171.8;
    double m = 798e-5;
    double xg = -3.7014/L;
    double Izz = 39.2e-5;
    double U = 15;
    int U_ = 1;

    double X0 = 0;
    double Xu = -184e-5;
    double Xuu = 0;
    double Xuuu = 0;
    double Xvv = 0;
    double Xvvu = 0;
    double Xrr = 0;
    double Xrru = 0;
    double Xoou = 0;
    double Xvr = 0;
    double Xvru = 0;
    double Xvo = 0;
    double Xvou = 0;
    double Xro = 0;
    double Xrou = 0;
    double Xudot = -42e-5;

    double Y0 = 0;
    double Y0u = 0;
    double Y0uu = 0;
    double Yv = -1160e-5;
    double Yvu = 0;
    double Yvuu = 0;
    double Yvvv = 0;
    double Yvrr = 0;
    double Yvoo = 0;
    double Yr = -499e-5;
    double Yru = 0;
    double Yruu = 0;
    double Yrrr = 0;
    double Yroo = 0;
    double Yo = 278e-5;
    double You = 0;
    double Youu = 0;
    double Yooo = 0;
    double Yvvr = 0;
    double Yvvo = 0;
    double Yrro = 0;
    double Yvro = 0;
    double Yvdot = -748e-5;
    double Yrdot = -9.354e-5;

    double N0 = 0;
    double N0u = 0;
    double N0uu = 0;
    double Nv = -264e-5;
    double Nvu = 0;
    double Nvuu = 0;
    double Nvvv = 0;
    double Nvrr = 0;
    double Nvoo = 0;
    double Nr = -166e-5;
    double Nru = 0;
    double Nruu = 0;
    double Nrrr = 0;
    double Nroo = 0;
    double No = -139e-5;
    double Nou = 0;
    double Nouu = 0;
    double Nooo = 0;
    double Nvvr = 0;
    double Nvvo = 0;
    double Nrro = 0;
    double Nvro = 0;
    double Nvdot = 4.646e-5;
    double Nrdot = -43.8e-5;

    public:
        double u0 = 0;
        double v0 = 0;
        double r0 = 0;
        double t0 = 0;
        double tn = 0;
        double h = 0;
        vector<double> r_data;
        vector<double> v_data;
        vector<double> u_data;
        vector<double> x_data;
        vector<double> y_data;
        vector<double> theta_data;
        vector<double> o_data;
        vector<vector<double> >data;
        MatrixXd A_;
        Abkowitz(double _u0, double _v0,double _r0,double _tn,double _h)
        {
            u0 = _u0/U;
            v0 = _v0/U;
            r0 = _r0*L/U;
            tn = _tn;
            h = _h;
            A_ = MatrixXd::Random(2, 2);
        }

        double* solve_realtime(double dt, double ut, double vt, double rt, double yt, double ot, double o_aim, double Xt, double Yt, double theta_t)
        {
            double ydot, rdot, odot, vdot;
            double C, T1T2, T1PlusT2, K, T3;
            C = Yv * Nr + Nv * ((m + 0) * U_ - Yr);
            T1T2 = (m + Yvdot) * (Izz + Nrdot) / C;
            T1PlusT2 = ((m + Yvdot) * Nr + (Izz + Nrdot) * Yv) / C;
            K = (Nv * Yo - Yv * No) / C;
            T3 = (m + Yvdot) * No / (Nv * Yo - Yv * No);
            rdot = yt;
            odot = -(ot - o_aim) / 2.5;
            vdot = (Yo * ot - (m * xg - Yrdot) * rdot - (m * U_ - Yr) * rt + Yv * vt) / (m - Yvdot);
            ydot = (K * ot + K * T3 * odot - (T1PlusT2)*yt - rt) / T1T2;

            //RK4
            double K1v = vdot;
            double K1y = ydot;
            double K1r = rdot;
            double K1o = odot;

            double K2r = yt + dt / 2 * K1y;
            double K2v = (Yo * (ot + dt / 2 * K1o) - (m * xg - Yrdot) * K2r - (m * U_ - Yr) * (rt + dt / 2 * K1r) + Yv * (vt + dt / 2 * K1v)) / (m - Yvdot);
            double K2o = -((ot + dt / 2 * K1o) - o_aim) / 2.5;
            double K2y = (K * (ot + dt / 2 * K1o) + K * T3 * K2o - (T1PlusT2) * (yt + dt / 2 * K1y) - (rt + dt / 2 * K1r)) / T1T2;

            double K3r = yt + dt / 2 * K2y;
            double K3v = (Yo * (ot + dt / 2 * K1o) - (m * xg - Yrdot) * K3r - (m * U_ - Yr) * (rt + dt / 2 * K1r) + Yv * (vt + dt / 2 * K1v)) / (m - Yvdot);
            double K3o = -((ot + dt / 2 * K2o) - o_aim) / 2.5;
            double K3y = (K * (ot + dt / 2 * K2o) + K * T3 * K3o - (T1PlusT2) * (yt + dt / 2 * K2y) - (rt + dt / 2 * K2r)) / T1T2;

            double K4r = yt + dt / 2 * K3y;
            double K4v = (Yo * (ot + dt / 2 * K1o) - (m * xg - Yrdot) * K4r - (m * U_ - Yr) * (rt + dt / 2 * K1r) + Yv * (vt + dt / 2 * K1v)) / (m - Yvdot);
            double K4o = -((ot + dt / 2 * K3o) - o_aim) / 2.5;
            double K4y = (K * (ot + dt / 2 * K3o) + K * T3 * K4o - (T1PlusT2) * (yt + dt / 2 * K3y) - (rt + dt / 2 * K3r)) / T1T2;

            yt += dt / 6 * (K1y + 2 * K2y + 2 * K3y + K4y);
            vt += dt / 6 * (K1v + 2 * K2v + 2 * K3v + K4v);
            rt += dt / 6 * (K1r + 2 * K2r + 2 * K3r + K4r);
            ot += dt / 6 * (K1o + 2 * K2o + 2 * K3o + K4o);

            ut = sqrt(1 - vt * vt);
            Xt += dt * (ut * cos(theta_t) - vt * sin(theta_t));
            Yt += dt * (ut * sin(theta_t) + vt * cos(theta_t));
            theta_t += dt * rt;
            
            double data_t[9] = { ut,vt,rt,yt,ot,o_aim,Xt,Yt,theta_t};
            return data_t;
        }

        void solve()
        {
            double y,ydot,rdot,odot,vdot;
            double u, v, r, o, theta, X=0, Y=0;
            double C, T1T2, T1PlusT2, K, T3;
            for (double t = t0; t < tn; t += h)
            {
                //double o_aim = t < 20 ? 20 * 3.14 / 180 : t < 40 ? -20 * 3.14 / 180 : t < 60 ? 20 * 3.14 / 180 : t < 70 ? -20 * 3.14 / 180 : 0;
                double o_aim = t < 5 ? 0 : DegToRad(20);
                C = Yv * Nr + Nv * ((m + 0) * U_-Yr);
                T1T2 = (m + Yvdot) * (Izz + Nrdot) / C;
                T1PlusT2 = ((m + Yvdot) * Nr + (Izz + Nrdot) * Yv) / C;
                K = (Nv * Yo - Yv * No) / C;
                T3 = (m + Yvdot) * No / (Nv * Yo - Yv * No);
                if (t == t0)
                {
                    u = u0;
                    v = v0;
                    r = r0;
                    o = 0;
                    y = 0;
                    odot = 0;
                    theta = 0;
                    rdot = y;
                    vdot = (Yo * o - (m * xg - Yrdot) * rdot - (m * U_ - Yr) * r + Yv * v) / (m - Yvdot);
                    ydot = (K * o + K * T3 * odot - (T1PlusT2)*y - r) / T1T2;
                    odot = -(o - o_aim) / 2.5;
                }
                else
                {
                    rdot = y;
                    vdot = (Yo * o - (m * xg - Yrdot) * rdot - (m * U_ - Yr) * r + Yv * v) / (m - Yvdot);
                    ydot = (K * o + K * T3 * odot - (T1PlusT2)*y - r) / T1T2;
                    odot = -(o - o_aim) / 2.5;
                }
                
                //RK4
                double K1v = vdot;
                double K1y = ydot;
                double K1r = rdot;
                double K1o = odot;

                double K2r = y + h / 2 * K1y;
                double K2v = (Yo * (o + h / 2 * K1o) - (m * xg - Yrdot) * K2r - (m * U_ - Yr) * (r + h / 2 * K1r) + Yv * (v + h / 2 * K1v)) / (m - Yvdot);
                double K2o = -((o + h / 2 * K1o) - o_aim) / 2.5;
                double K2y = (K * (o + h / 2 * K1o) + K * T3 * K2o - (T1PlusT2) * (y + h / 2 * K1y) - (r + h / 2 * K1r)) / T1T2;
                
                double K3r = y + h / 2 * K2y;
                double K3v = (Yo * (o + h / 2 * K1o) - (m * xg - Yrdot) * K3r - (m * U_ - Yr) * (r + h / 2 * K1r) + Yv * (v + h / 2 * K1v)) / (m - Yvdot);
                double K3o = -((o + h / 2 * K2o) - o_aim) / 2.5;
                double K3y = (K * (o + h / 2 * K2o) + K * T3 * K3o - (T1PlusT2) * (y + h / 2 * K2y) - (r + h / 2 * K2r)) / T1T2;
                
                double K4r = y + h / 2 * K3y;
                double K4v = (Yo * (o + h / 2 * K1o) - (m * xg - Yrdot) * K4r - (m * U_ - Yr) * (r + h / 2 * K1r) + Yv * (v + h / 2 * K1v)) / (m - Yvdot);
                double K4o = -((o + h / 2 * K3o) - o_aim) / 2.5;
                double K4y = (K * (o + h / 2 * K3o) + K * T3 * K4o - (T1PlusT2) * (y + h / 2 * K3y) - (r + h / 2 * K3r)) / T1T2;
                
                y += h / 6 * (K1y + 2 * K2y + 2 * K3y + K4y);
                v += h / 6 * (K1v + 2 * K2v + 2 * K3v + K4v);
                r += h / 6 * (K1r + 2 * K2r + 2 * K3r + K4r);
                o += h / 6 * (K1o + 2 * K2o + 2 * K3o + K4o);

                u = sqrt(1 - v * v);
                theta += h * r;
                X += h * (u * cos(theta) - v * sin(theta));
                Y += h * (u * sin(theta) + v * cos(theta));
                r_data.push_back(r);
                v_data.push_back(v);
                u_data.push_back(u);
                theta_data.push_back(theta);
                x_data.push_back(X);
                y_data.push_back(Y);
                o_data.push_back(o);
                /*
                    std::cout << "t=" << t << std::endl;
                    std::cout << "u=" << u << std::endl;
                    std::cout << "o=" << o << std::endl;
                    std::cout << "v=" << v << std::endl;
                    std::cout << "odot=" << odot << std::endl;
                    std::cout << "r=" << r << std::endl;
                    std::cout << "theta=" << theta << std::endl;
                    std::cout << std::endl;
                    */
            }
            data.push_back(u_data);
            data.push_back(v_data);
            data.push_back(r_data);
            data.push_back(o_data);
            data.push_back(x_data);
            data.push_back(y_data);
            data.push_back(theta_data);
        }
};

int main()
{
    Abkowitz model(15,0,0.0,0,0);
    /*
    string filename = "H:/code/MSVCprojects/ConsoleApplication1Abko/x64/Debug/test.csv";
    model.solve();
    write(filename, model.data);
    */
    double dt = 0.01;
    double ut = 1;
    double vt = 0;
    double rt = 0;
    double yt = 0;
    double ot = 0;
    double Xt = 0;
    double Yt = 0;
    double theta_t = 0;
    double* data_t = (double*)malloc(sizeof(double) * 9);
    if (data_t != NULL)
    {
        memset(data_t, 0, sizeof(double) * 9);
    }
    for (double t = 0; t <= 100; t += dt)
    {
        double o_aim = t < 5 ? 0 : DegToRad(20);
        memcpy_s(data_t, sizeof(double) * 9, model.solve_realtime(dt, ut, vt, rt, yt, ot, o_aim, Xt, Yt, theta_t), sizeof(double) * 9);
        if (data_t != NULL)
        {
            ut = data_t[0];
            vt = data_t[1];
            rt = data_t[2];
            yt = data_t[3];
            ot = data_t[4];
            o_aim = data_t[5];
            Xt = data_t[6];
            Yt = data_t[7];
            theta_t = data_t[8];
        }
        std::cout << "t=" << t << std::endl;
        std::cout << "u=" << data_t[0] << std::endl;
        std::cout << "v=" << data_t[1] << std::endl;
        std::cout << "r=" << data_t[2] << std::endl;
        std::cout << "ot=" << data_t[3] << std::endl;
        std::cout << std::endl;
    }
    
    //std::cout << model.A_ << std::endl;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
