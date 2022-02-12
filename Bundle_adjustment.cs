using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;
using System.IO;
using System.Windows.Forms;
using System.Data.SqlClient;
using System.Collections;
using Accord.Math;
using Accord.Math.Decompositions;
using System.IO.MemoryMappedFiles;
using System.Diagnostics;
using inv2;
using MathWorks.MATLAB.NET.Arrays;
using MathWorks.MATLAB.NET.Utility;
using outlier_bundle;
using System.Threading;
namespace theses
{
    class Bundle_adjustment
    {
       public static object obj = new object();
       public static void thread_bundle(int ii, int ie, string[] image_name_Eterior, List<string> point_id, List<int> count_camera, double[,] observe_value_points_weight, double[,] cam_observe, List<double> camera_weight, double[,] Exterior_Orientation_Eterior, double[] interior, double[,] Initial_value, double[,] image_m, double[,] image_sin_cos, double[,] observe_value_camera, double[,] observe_value_points, out double[,] S_out, out double[,] eawveb_out, out double[,] llln_out)
        {
            int shomwhile = 0;
            double[,] S = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
            double[,] eawveb = new double[6 * (image_name_Eterior.Length), 1];
            double[,] llln = new double[6 * (image_name_Eterior.Length), 1];
            for (int i1 = 0; i1 < ii; i1++)
            {
                shomwhile += count_camera[i1];
            }
            for (int i = ii; i < ie; i++)
            {


                double[,] v = new double[3, 3];
                List<double> w1 = new List<double>();
                List<double> w2 = new List<double>();
                List<double> w3 = new List<double>();
                List<int> wj = new List<int>();
                int ij2 = shomwhile;
                shomwhile += count_camera[i];

                if (point_id[i] != "non")
                {
                    double[,] pp2nB = new double[3, 3];

                    pp2nB[0, 0] = observe_value_points_weight[i, 0];
                    pp2nB[1, 1] = observe_value_points_weight[i, 1];
                    pp2nB[2, 2] = observe_value_points_weight[i, 2];

                    double[,] ebii = new double[3, 1];
                    double[,] lllnb = new double[3, 1];
                    double[,] pp2 = new double[2, 2];
                    double[,] pp2nA = new double[6, 6];
                    for (int ij = ij2; ij < shomwhile; ij++)
                    {

                        int j2 = Convert.ToInt32(cam_observe[ij, 0]);

                        pp2[0, 0] = camera_weight[8 * j2];
                        pp2[1, 1] = camera_weight[8 * j2 + 1];
                        pp2nA[0, 0] = camera_weight[8 * j2 + 2];
                        pp2nA[1, 1] = camera_weight[8 * j2 + 3];
                        pp2nA[2, 2] = camera_weight[8 * j2 + 4];
                        pp2nA[3, 3] = camera_weight[8 * j2 + 5];
                        pp2nA[4, 4] = camera_weight[8 * j2 + 6];
                        pp2nA[5, 5] = camera_weight[8 * j2 + 7];

                        double XL = Exterior_Orientation_Eterior[j2, 0];
                        double YL = Exterior_Orientation_Eterior[j2, 1];
                        double ZL = Exterior_Orientation_Eterior[j2, 2];
                        double om = Exterior_Orientation_Eterior[j2, 3];
                        double phi = Exterior_Orientation_Eterior[j2, 4];
                        double k = Exterior_Orientation_Eterior[j2, 5];
                        double ff = interior[0];
                        double x0 = interior[1];
                        double y0 = interior[2];
                        double K1 = interior[3];
                        double K2 = interior[4];
                        double K3 = interior[5];
                        double P1 = interior[6];
                        double P2 = interior[7];
                        double B1 = interior[8];
                        double B2 = interior[9];
                        double X = Initial_value[i, 0];
                        double Y = Initial_value[i, 1];
                        double Z = Initial_value[i, 2];
                        //      MessageBox.Show("" + X+"   "+Y+"    "+Z);
                        double m11 = image_m[j2, 0];
                        double m12 = image_m[j2, 1];
                        //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                        double m13 = image_m[j2, 2];
                        //m21=-cos(phi)*sin(k)
                        double m21 = image_m[j2, 3];
                        //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                        double m22 = image_m[j2, 4];
                        //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                        double m23 = image_m[j2, 5];
                        //m31=sin(phi)
                        double m31 = image_m[j2, 6];
                        //m32=-sin(om)*cos(phi)
                        double m32 = image_m[j2, 7];
                        //m33=cos(om)*cos(phi)
                        double m33 = image_m[j2, 8];

                        double sinom = image_sin_cos[j2, 0]; double sinphi = image_sin_cos[j2, 1]; double sink = image_sin_cos[j2, 2];
                        double cosom = image_sin_cos[j2, 3]; double cosphi = image_sin_cos[j2, 4]; double cosk = image_sin_cos[j2, 5];


                        double xx = cam_observe[ij, 1];
                        double yy = cam_observe[ij, 2];
                        double sinphiXXL = (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL));
                        double sinksinom = (sink * sinom - cosk * cosom * sinphi);
                        double cosomsink = (cosom * sink + cosk * sinom * sinphi);
                        double coskcosom = (cosk * cosom - sink * sinom * sinphi);
                        double cosksinom = (cosk * sinom + cosom * sink * sinphi);
                        double xom = ((ff * (sinksinom * (Y - YL) - cosomsink * (Z - ZL))) / sinphiXXL - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));

                        double xphi = ((ff * (cosk * sinphi * (X - XL) + cosk * cosom * cosphi * (Z - ZL) - cosk * cosphi * sinom * (Y - YL))) / sinphiXXL + (ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double xk = (-(ff * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / sinphiXXL);
                        double xXL = ((ff * cosk * cosphi) / sinphiXXL - (ff * sinphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double xYL = ((ff * cosomsink) / sinphiXXL + (ff * cosphi * sinom * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double xZL = (-((ff * cosom * cosphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * sinksinom) / sinphiXXL));
                        //   (f * cos(om) * cos(phi) * ((cos(om) * sin(k) + cos(k) * sin(om) * sin(phi)) * (Y - YL) + (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi)) * (Z - ZL) + cos(k) * cos(phi) * (X - XL))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL)) ^ 2 - (f * (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL))];
                        double yom = ((ff * (cosksinom * (Y - YL) - coskcosom * (Z - ZL))) / sinphiXXL - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double yphi = ((ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * (sink * sinphi * (X - XL) + cosom * cosphi * sink * (Z - ZL) - cosphi * sink * sinom * (Y - YL))) / sinphiXXL);
                        double yk = ((ff * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / sinphiXXL);
                        double yXL = (-(ff * sinphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosphi * sink) / sinphiXXL);
                        double yYL = ((ff * coskcosom) / sinphiXXL + (ff * cosphi * sinom * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double yZL = ((ff * cosksinom) / sinphiXXL - (ff * cosom * cosphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double xX = ((ff * sinphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosk * cosphi) / sinphiXXL);
                        double xY = (-(ff * cosomsink) / sinphiXXL - (ff * cosphi * sinom * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double xZ = ((ff * cosom * cosphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * sinksinom) / sinphiXXL);
                        double yX = ((ff * sinphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) + (ff * cosphi * sink) / sinphiXXL);
                        double yY = (-(ff * coskcosom) / sinphiXXL - (ff * cosphi * sinom * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                        double yZ = ((ff * cosom * cosphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosksinom) / sinphiXXL);



                        double xpi = (xx - x0);
                        double ypi = (yy - y0);

                        double ri = Math.Sqrt(((xpi * xpi) + (ypi * ypi)));
                        double deltax = (xpi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P1 * (((ri * ri)) + (2 * ((xpi * xpi))))) + (2 * P2 * xpi * ypi) + (B1 * xpi) + (B2 * ypi);
                        double deltay = (ypi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P2 * (((ri * ri)) + (2 * ((ypi * ypi))))) + (2 * P1 * xpi * ypi);
                        double x = -(ff * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltax;
                        double y = -(ff * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltay;
                        //    MessageBox.Show("" + deltax + "   " + deltay);
                        // MessageBox.Show("" + (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)))+"   "+ Math.Pow((x0 - xx), 2));

                        //   double[,] B = new double[2, 2];


                        double a00 = xom;

                        double a01 = xphi;

                        double a02 = xk;

                        double a03 = xXL;

                        double a04 = xYL;

                        double a05 = xZL;

                        double a10 = yom;

                        double a11 = yphi;

                        double a12 = yk;

                        double a13 = yXL;

                        double a14 = yYL;

                        double a15 = yZL;

                        double bb00 = xX;

                        double bb01 = xY;

                        double bb02 = xZ;

                        double bb10 = yX;

                        double bb11 = yY;

                        double bb12 = yZ;


                        //  MessageBox.Show("" + c[0, 0] + "   " + c[0, 1]+ "   " + c[0, 2] + "   " + c[0, 3] + "   " + c[0, 4] + "   " + c[0, 5] + "   " + c[0, 6] + "   " + c[0, 7] + "   " + c[0, 8] + "   " + c[0, 9]);
                        // MessageBox.Show("" + c[1, 0] + "   " + c[1, 1] + "   " + c[1, 2] + "   " + c[1, 3] + "   " + c[1, 4] + "   " + c[1, 5] + "   " + c[1, 6] + "   " + c[1, 7] + "   " + c[1, 8] + "   " + c[1, 9]);

                        llln[6 * j2, 0] = observe_value_camera[j2, 0] - om; llln[6 * j2 + 1, 0] = observe_value_camera[j2, 1] - phi; llln[6 * j2 + 2, 0] = observe_value_camera[j2, 2] - k;
                        llln[6 * j2 + 3, 0] = observe_value_camera[j2, 3] - XL; llln[6 * j2 + 4, 0] = observe_value_camera[j2, 4] - YL; llln[6 * j2 + 5, 0] = observe_value_camera[j2, 5] - ZL; lllnb[0, 0] = observe_value_points[i, 0] - X; lllnb[1, 0] = observe_value_points[i, 1] - Y; lllnb[2, 0] = observe_value_points[i, 2] - Z;
                        double ll0 = xx - x; ; double ll1 = yy - y;
                        //if (iteration == 10)
                        // {
                        //     MessageBox.Show("" + lll[0, 0] + "   " + lll[1, 0]);
                        // }
                        double pb00 = pp2[0, 0]; double pb01 = 0; double pb10 = 0; double pb11 = pp2[1, 1];
                        double[,] eajj = new double[6, 1];
                        eajj[0, 0] = ll0 * (pb00 * (a00) + pb10 * (a10)) + ll1 * (pb01 * (a00) + pb11 * (a10));
                        eajj[1, 0] = ll0 * (pb00 * (a01) + pb10 * (a11)) + ll1 * (pb01 * (a01) + pb11 * (a11));
                        eajj[2, 0] = ll0 * (pb00 * (a02) + pb10 * (a12)) + ll1 * (pb01 * (a02) + pb11 * (a12));
                        eajj[3, 0] = ll0 * (pb00 * (a03) + pb10 * (a13)) + ll1 * (pb01 * (a03) + pb11 * (a13));
                        eajj[4, 0] = ll0 * (pb00 * (a04) + pb10 * (a14)) + ll1 * (pb01 * (a04) + pb11 * (a14));
                        eajj[5, 0] = ll0 * (pb00 * (a05) + pb10 * (a15)) + ll1 * (pb01 * (a05) + pb11 * (a15));

                        //      MessageBox.Show("" + lll[1, 0] + "   " + lll[1, 0] + "   " + lll[2, 0] + "   " + lll[3, 0] + "   " + lll[4, 0] + "   " + lll[5, 0] + "   " + lll[6, 0] + "   " + lll[7, 0] + "   " + lll[8, 0] + "   " + lll[9, 0] + "   " + lll[10, 0] + "   " + lll[11, 0] + "   " + lll[12, 0] + "   " + lll[13, 0] + "   " + lll[14, 0] + "   " + lll[15, 0] + "   " + lll[16, 0] + "   " + lll[17, 0] + "   " + lll[18, 0] + "   " + lll[19, 0] + "   " + lll[20, 0]);

                        for (int pp = 0; pp < 6; pp++)
                        {
                            eawveb[6 * j2 + pp, 0] = eawveb[6 * j2 + pp, 0] + eajj[pp, 0];

                        }
                        //     MessageBox.Show("1sina" + eawveb[6 * 0 + 3, 0]);

                        double[,] ebii2 = new double[3, 1];
                        ebii2[0, 0] = ll0 * (pb00 * (bb00) + pb10 * (bb10)) + ll1 * (pb01 * (bb00) + pb11 * (bb10));
                        ebii2[1, 0] = ll0 * (pb00 * (bb01) + pb10 * (bb11)) + ll1 * (pb01 * (bb01) + pb11 * (bb11));
                        ebii2[2, 0] = ll0 * (pb00 * (bb02) + pb10 * (bb12)) + ll1 * (pb01 * (bb02) + pb11 * (bb12));
                        ebii[0, 0] = ebii[0, 0] + ebii2[0, 0]; ebii[1, 0] = ebii[1, 0] + ebii2[1, 0];
                        ebii[2, 0] = ebii[2, 0] + ebii2[2, 0];
                        //  MessageBox.Show("" + ebii[0, 0] + "   " + ebii2[0, 0]);
                        double[,] v2 = new double[3, 3];
                        v2[0, 0] = bb00 * (pb00 * (bb00) + pb10 * (bb10)) + bb10 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 1] = bb01 * (pb00 * (bb00) + pb10 * (bb10)) + bb11 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 2] = bb02 * (pb00 * (bb00) + pb10 * (bb10)) + bb12 * (pb01 * (bb00) + pb11 * (bb10));
                        v2[1, 0] = bb00 * (pb00 * (bb01) + pb10 * (bb11)) + bb10 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 1] = bb01 * (pb00 * (bb01) + pb10 * (bb11)) + bb11 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 2] = bb02 * (pb00 * (bb01) + pb10 * (bb11)) + bb12 * (pb01 * (bb01) + pb11 * (bb11));
                        v2[2, 0] = bb00 * (pb00 * (bb02) + pb10 * (bb12)) + bb10 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 1] = bb01 * (pb00 * (bb02) + pb10 * (bb12)) + bb11 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 2] = bb02 * (pb00 * (bb02) + pb10 * (bb12)) + bb12 * (pb01 * (bb02) + pb11 * (bb12));
                        double[,] ww = new double[6, 3];

                        ww[0, 0] = bb00 * (pb00 * (a00) + pb10 * (a10)) + bb10 * (pb01 * (a00) + pb11 * (a10)); ww[0, 1] = bb01 * (pb00 * (a00) + pb10 * (a10)) + bb11 * (pb01 * (a00) + pb11 * (a10)); ww[0, 2] = bb02 * (pb00 * (a00) + pb10 * (a10)) + bb12 * (pb01 * (a00) + pb11 * (a10));
                        ww[1, 0] = bb00 * (pb00 * (a01) + pb10 * (a11)) + bb10 * (pb01 * (a01) + pb11 * (a11)); ww[1, 1] = bb01 * (pb00 * (a01) + pb10 * (a11)) + bb11 * (pb01 * (a01) + pb11 * (a11)); ww[1, 2] = bb02 * (pb00 * (a01) + pb10 * (a11)) + bb12 * (pb01 * (a01) + pb11 * (a11));
                        ww[2, 0] = bb00 * (pb00 * (a02) + pb10 * (a12)) + bb10 * (pb01 * (a02) + pb11 * (a12)); ww[2, 1] = bb01 * (pb00 * (a02) + pb10 * (a12)) + bb11 * (pb01 * (a02) + pb11 * (a12)); ww[2, 2] = bb02 * (pb00 * (a02) + pb10 * (a12)) + bb12 * (pb01 * (a02) + pb11 * (a12));
                        ww[3, 0] = bb00 * (pb00 * (a03) + pb10 * (a13)) + bb10 * (pb01 * (a03) + pb11 * (a13)); ww[3, 1] = bb01 * (pb00 * (a03) + pb10 * (a13)) + bb11 * (pb01 * (a03) + pb11 * (a13)); ww[3, 2] = bb02 * (pb00 * (a03) + pb10 * (a13)) + bb12 * (pb01 * (a03) + pb11 * (a13));
                        ww[4, 0] = bb00 * (pb00 * (a04) + pb10 * (a14)) + bb10 * (pb01 * (a04) + pb11 * (a14)); ww[4, 1] = bb01 * (pb00 * (a04) + pb10 * (a14)) + bb11 * (pb01 * (a04) + pb11 * (a14)); ww[4, 2] = bb02 * (pb00 * (a04) + pb10 * (a14)) + bb12 * (pb01 * (a04) + pb11 * (a14));
                        ww[5, 0] = bb00 * (pb00 * (a05) + pb10 * (a15)) + bb10 * (pb01 * (a05) + pb11 * (a15)); ww[5, 1] = bb01 * (pb00 * (a05) + pb10 * (a15)) + bb11 * (pb01 * (a05) + pb11 * (a15)); ww[5, 2] = bb02 * (pb00 * (a05) + pb10 * (a15)) + bb12 * (pb01 * (a05) + pb11 * (a15));

                        w1.Add(ww[0, 0]);
                        w1.Add(ww[1, 0]);
                        w1.Add(ww[2, 0]);
                        w1.Add(ww[3, 0]);
                        w1.Add(ww[4, 0]);
                        w1.Add(ww[5, 0]);

                        w2.Add(ww[0, 1]);
                        w2.Add(ww[1, 1]);
                        w2.Add(ww[2, 1]);
                        w2.Add(ww[3, 1]);
                        w2.Add(ww[4, 1]);
                        w2.Add(ww[5, 1]);

                        w3.Add(ww[0, 2]);
                        w3.Add(ww[1, 2]);
                        w3.Add(ww[2, 2]);
                        w3.Add(ww[3, 2]);
                        w3.Add(ww[4, 2]);
                        w3.Add(ww[5, 2]);

                        wj.Add(j2);
                        //WW[6 * j2 + 0, 0] = ww[0, 0]; WW[6 * j2 + 1, 0] = ww[1, 0]; WW[6 * j2 + 2, 0] = ww[2, 0]; WW[6 * j2 + 3, 0] = ww[3, 0]; WW[6 * j2 + 4, 0] = ww[4, 0]; WW[6 * j2 + 5, 0] = ww[5, 0];
                        //WW[6 * j2 + 0, 1] = ww[0, 1]; WW[6 * j2 + 1, 1] = ww[1,1]; WW[6 * j2 + 2, 1] = ww[2, 1]; WW[6 * j2 + 3, 1] = ww[3, 1]; WW[6 * j2 + 4, 1] = ww[4, 1]; WW[6 * j2 + 5, 1] = ww[5, 1];
                        //WW[6 * j2 + 0, 2] = ww[0, 2]; WW[6 * j2 + 1, 2] = ww[1, 2]; WW[6 * j2 + 2, 2] = ww[2, 2]; WW[6 * j2 + 3, 2] = ww[3,2]; WW[6 * j2 + 4,2] = ww[4, 2]; WW[6 * j2 + 5, 2] = ww[5, 2];
                        v[0, 0] = v[0, 0] + v2[0, 0]; v[1, 0] = v[1, 0] + v2[1, 0]; v[2, 0] = v[2, 0] + v2[2, 0];
                        v[0, 1] = v[0, 1] + v2[0, 1]; v[1, 1] = v[1, 1] + v2[1, 1]; v[2, 1] = v[2, 1] + v2[2, 1];
                        v[0, 2] = v[0, 2] + v2[0, 2]; v[1, 2] = v[1, 2] + v2[1, 2]; v[2, 2] = v[2, 2] + v2[2, 2];


                        double[,] u = new double[6, 6];

                        //  MessageBox.Show("" + F[0, 0]+"   "+ FF[0, 0]+"    "+ BB[0,0]+"  "+j2+"  "+i);
                        u[0, 0] = a00 * (pb00 * (a00) + pb10 * (a10)) + a10 * (pb01 * (a00) + pb11 * (a10)); u[0, 1] = a01 * (pb00 * (a00) + pb10 * (a10)) + a11 * (pb01 * (a00) + pb11 * (a10)); u[0, 2] = a02 * (pb00 * (a00) + pb10 * (a10)) + a12 * (pb01 * (a00) + pb11 * (a10)); u[0, 3] = a03 * (pb00 * (a00) + pb10 * (a10)) + a13 * (pb01 * (a00) + pb11 * (a10)); u[0, 4] = a04 * (pb00 * (a00) + pb10 * (a10)) + a14 * (pb01 * (a00) + pb11 * (a10)); u[0, 5] = a05 * (pb00 * (a00) + pb10 * (a10)) + a15 * (pb01 * (a00) + pb11 * (a10));
                        u[1, 0] = a00 * (pb00 * (a01) + pb10 * (a11)) + a10 * (pb01 * (a01) + pb11 * (a11)); u[1, 1] = a01 * (pb00 * (a01) + pb10 * (a11)) + a11 * (pb01 * (a01) + pb11 * (a11)); u[1, 2] = a02 * (pb00 * (a01) + pb10 * (a11)) + a12 * (pb01 * (a01) + pb11 * (a11)); u[1, 3] = a03 * (pb00 * (a01) + pb10 * (a11)) + a13 * (pb01 * (a01) + pb11 * (a11)); u[1, 4] = a04 * (pb00 * (a01) + pb10 * (a11)) + a14 * (pb01 * (a01) + pb11 * (a11)); u[1, 5] = a05 * (pb00 * (a01) + pb10 * (a11)) + a15 * (pb01 * (a01) + pb11 * (a11));
                        u[2, 0] = a00 * (pb00 * (a02) + pb10 * (a12)) + a10 * (pb01 * (a02) + pb11 * (a12)); u[2, 1] = a01 * (pb00 * (a02) + pb10 * (a12)) + a11 * (pb01 * (a02) + pb11 * (a12)); u[2, 2] = a02 * (pb00 * (a02) + pb10 * (a12)) + a12 * (pb01 * (a02) + pb11 * (a12)); u[2, 3] = a03 * (pb00 * (a02) + pb10 * (a12)) + a13 * (pb01 * (a02) + pb11 * (a12)); u[2, 4] = a04 * (pb00 * (a02) + pb10 * (a12)) + a14 * (pb01 * (a02) + pb11 * (a12)); u[2, 5] = a05 * (pb00 * (a02) + pb10 * (a12)) + a15 * (pb01 * (a02) + pb11 * (a12));
                        u[3, 0] = a00 * (pb00 * (a03) + pb10 * (a13)) + a10 * (pb01 * (a03) + pb11 * (a13)); u[3, 1] = a01 * (pb00 * (a03) + pb10 * (a13)) + a11 * (pb01 * (a03) + pb11 * (a13)); u[3, 2] = a02 * (pb00 * (a03) + pb10 * (a13)) + a12 * (pb01 * (a03) + pb11 * (a13)); u[3, 3] = a03 * (pb00 * (a03) + pb10 * (a13)) + a13 * (pb01 * (a03) + pb11 * (a13)); u[3, 4] = a04 * (pb00 * (a03) + pb10 * (a13)) + a14 * (pb01 * (a03) + pb11 * (a13)); u[3, 5] = a05 * (pb00 * (a03) + pb10 * (a13)) + a15 * (pb01 * (a03) + pb11 * (a13));
                        u[4, 0] = a00 * (pb00 * (a04) + pb10 * (a14)) + a10 * (pb01 * (a04) + pb11 * (a14)); u[4, 1] = a01 * (pb00 * (a04) + pb10 * (a14)) + a11 * (pb01 * (a04) + pb11 * (a14)); u[4, 2] = a02 * (pb00 * (a04) + pb10 * (a14)) + a12 * (pb01 * (a04) + pb11 * (a14)); u[4, 3] = a03 * (pb00 * (a04) + pb10 * (a14)) + a13 * (pb01 * (a04) + pb11 * (a14)); u[4, 4] = a04 * (pb00 * (a04) + pb10 * (a14)) + a14 * (pb01 * (a04) + pb11 * (a14)); u[4, 5] = a05 * (pb00 * (a04) + pb10 * (a14)) + a15 * (pb01 * (a04) + pb11 * (a14));
                        u[5, 0] = a00 * (pb00 * (a05) + pb10 * (a15)) + a10 * (pb01 * (a05) + pb11 * (a15)); u[5, 1] = a01 * (pb00 * (a05) + pb10 * (a15)) + a11 * (pb01 * (a05) + pb11 * (a15)); u[5, 2] = a02 * (pb00 * (a05) + pb10 * (a15)) + a12 * (pb01 * (a05) + pb11 * (a15)); u[5, 3] = a03 * (pb00 * (a05) + pb10 * (a15)) + a13 * (pb01 * (a05) + pb11 * (a15)); u[5, 4] = a04 * (pb00 * (a05) + pb10 * (a15)) + a14 * (pb01 * (a05) + pb11 * (a15)); u[5, 5] = a05 * (pb00 * (a05) + pb10 * (a15)) + a15 * (pb01 * (a05) + pb11 * (a15));


                        for (int pp = 0; pp < 6; pp++)
                        {
                            for (int ppp = 0; ppp < 6; ppp++)
                            {

                                S[6 * j2 + pp, 6 * j2 + ppp] = S[6 * j2 + pp, 6 * j2 + ppp] + u[pp, ppp];

                            }
                        }

                    }


                    //   MessageBox.Show("" + eaj[0, 0] + "   " + eaj[6 * image_name_Eterior.Length + 1, 0]);
                    v[0, 0] = v[0, 0] + pp2nB[0, 0]; v[1, 1] = v[1, 1] + pp2nB[1, 1]; v[2, 2] = v[2, 2] + pp2nB[2, 2];
                    ebii[0, 0] = ebii[0, 0] + (pp2nB[0, 0] * lllnb[0, 0]);
                    ebii[1, 0] = ebii[1, 0] + (pp2nB[1, 1] * lllnb[1, 0]);
                    ebii[2, 0] = ebii[2, 0] + (pp2nB[2, 2] * lllnb[2, 0]);
                    double ebii00 = ebii[0, 0] - (pp2nB[0, 0] * lllnb[0, 0]);
                    double ebii10 = ebii[1, 0] - (pp2nB[1, 1] * lllnb[1, 0]);
                    double ebii20 = ebii[2, 0] - (pp2nB[2, 2] * lllnb[2, 0]);
                    double[,] v3 = new double[3, 3];
                    double detv = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                    v3[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                    v3[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                    v3[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                    v3[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                    v3[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                    v3[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                    v3[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                    v3[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                    v3[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                    double v300 = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                    double v301 = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                    double v302 = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                    double v310 = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                    double v311 = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                    double v312 = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                    double v320 = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                    double v321 = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                    double v322 = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;


                    for (int j2 = 0; j2 < wj.Count; j2++)
                    {
                        double[,] wcalj2 = new double[6, 3];
                        double wcalj200, wcalj201, wcalj202, wcalj210, wcalj211, wcalj212, wcalj220, wcalj221, wcalj222, wcalj230, wcalj231, wcalj232, wcalj240, wcalj241, wcalj242, wcalj250, wcalj251, wcalj252;
                        for (int iop = 0; iop < 6; iop++)
                        {

                            wcalj2[iop, 0] = w1[6 * j2 + iop];
                            wcalj2[iop, 1] = w2[6 * j2 + iop];
                            wcalj2[iop, 2] = w3[6 * j2 + iop];

                        }
                        wcalj200 = wcalj2[0, 0]; wcalj201 = wcalj2[0, 1]; wcalj202 = wcalj2[0, 2]; wcalj210 = wcalj2[1, 0]; wcalj211 = wcalj2[1, 1]; wcalj212 = wcalj2[1, 2]; wcalj220 = wcalj2[2, 0]; wcalj221 = wcalj2[2, 1]; wcalj222 = wcalj2[2, 2]; wcalj230 = wcalj2[3, 0]; wcalj231 = wcalj2[3, 1]; wcalj232 = wcalj2[3, 2]; wcalj240 = wcalj2[4, 0]; wcalj241 = wcalj2[4, 1]; wcalj242 = wcalj2[4, 2]; wcalj250 = wcalj2[5, 0]; wcalj251 = wcalj2[5, 1]; wcalj252 = wcalj2[5, 2];
                        //  double[,] eawveb2 = Matrix.Multiply((Matrix.Multiply(wcalj2, v3)), ebii);
                        double[,] eawveb2 = new double[6, 1];
                        eawveb2[0, 0] = ebii00 * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + ebii10 * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + ebii20 * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202);
                        eawveb2[1, 0] = ebii00 * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + ebii10 * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + ebii20 * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212);
                        eawveb2[2, 0] = ebii00 * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + ebii10 * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + ebii20 * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222);
                        eawveb2[3, 0] = ebii00 * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + ebii10 * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + ebii20 * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232);
                        eawveb2[4, 0] = ebii00 * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + ebii10 * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + ebii20 * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242);
                        eawveb2[5, 0] = ebii00 * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + ebii10 * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + ebii20 * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252);


                        for (int iop = 0; iop < 6; iop++)
                        {
                            eawveb[6 * wj[j2] + iop, 0] = eawveb[6 * wj[j2] + iop, 0] - eawveb2[iop, 0];
                        }
                        //   MessageBox.Show("" + eawveb2[0, 0]);
                        for (int j = 0; j < wj.Count; j++)
                        {
                            double[,] wcalj = new double[6, 3];
                            for (int iop = 0; iop < 6; iop++)
                            {
                                wcalj[iop, 0] = w1[6 * j + iop];
                                wcalj[iop, 1] = w2[6 * j + iop];
                                wcalj[iop, 2] = w3[6 * j + iop];
                            }
                            double wcalj00 = wcalj[0, 0]; double wcalj01 = wcalj[0, 1]; double wcalj02 = wcalj[0, 2]; double wcalj10 = wcalj[1, 0]; double wcalj11 = wcalj[1, 1]; double wcalj12 = wcalj[1, 2]; double wcalj20 = wcalj[2, 0]; double wcalj21 = wcalj[2, 1]; double wcalj22 = wcalj[2, 2]; double wcalj30 = wcalj[3, 0]; double wcalj31 = wcalj[3, 1]; double wcalj32 = wcalj[3, 2]; double wcalj40 = wcalj[4, 0]; double wcalj41 = wcalj[4, 1]; double wcalj42 = wcalj[4, 2]; double wcalj50 = wcalj[5, 0]; double wcalj51 = wcalj[5, 1]; double wcalj52 = wcalj[5, 2];

                            double[,] wcalj2wcalj = new double[6, 6];
                            wcalj2wcalj[0, 0] = (wcalj00) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj01) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj02) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 1] = (wcalj10) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj11) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj12) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 2] = (wcalj20) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj21) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj22) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 3] = (wcalj30) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj31) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj32) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 4] = (wcalj40) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj41) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj42) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 5] = (wcalj50) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj51) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj52) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202);
                            wcalj2wcalj[1, 0] = (wcalj00) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj01) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj02) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 1] = (wcalj10) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj11) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj12) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 2] = (wcalj20) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj21) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj22) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 3] = (wcalj30) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj31) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj32) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 4] = (wcalj40) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj41) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj42) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 5] = (wcalj50) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj51) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj52) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212);
                            wcalj2wcalj[2, 0] = (wcalj00) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj01) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj02) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 1] = (wcalj10) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj11) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj12) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 2] = (wcalj20) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj21) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj22) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 3] = (wcalj30) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj31) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj32) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 4] = (wcalj40) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj41) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj42) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 5] = (wcalj50) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj51) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj52) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222);
                            wcalj2wcalj[3, 0] = (wcalj00) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj01) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj02) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 1] = (wcalj10) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj11) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj12) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 2] = (wcalj20) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj21) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj22) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 3] = (wcalj30) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj31) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj32) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 4] = (wcalj40) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj41) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj42) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 5] = (wcalj50) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj51) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj52) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232);
                            wcalj2wcalj[4, 0] = (wcalj00) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj01) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj02) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 1] = (wcalj10) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj11) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj12) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 2] = (wcalj20) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj21) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj22) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 3] = (wcalj30) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj31) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj32) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 4] = (wcalj40) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj41) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj42) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 5] = (wcalj50) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj51) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj52) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242);
                            wcalj2wcalj[5, 0] = (wcalj00) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj01) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj02) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 1] = (wcalj10) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj11) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj12) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 2] = (wcalj20) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj21) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj22) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 3] = (wcalj30) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj31) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj32) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 4] = (wcalj40) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj41) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj42) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 5] = (wcalj50) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj51) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj52) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252);

                            for (int pp = 0; pp < 6; pp++)
                            {
                                for (int ppp = 0; ppp < 6; ppp++)
                                {
                                    S[6 * wj[j2] + pp, 6 * wj[j] + ppp] = S[6 * wj[j2] + pp, 6 * wj[j] + ppp] - wcalj2wcalj[pp, ppp];
                                }
                            }

                        }
                    }
                }

            }
            eawveb_out = eawveb; S_out = S; llln_out = llln;
        }
        public static void thread_bundle_interior_exterior(int ii, int ie, double[,] IOP , string[] image_name_Eterior, List<string> point_id, List<int> count_camera, double[,] observe_value_points_weight, double[,] cam_observe, List<double> camera_weight, double[,] Exterior_Orientation_Eterior, double[] interior, double[,] Initial_value, double[,] image_m, double[,] image_sin_cos, double[,] observe_value_camera, double[,] observe_value_points, out double[,] S_out, out double[,] eawveb_out, out double[,] llln_out, out double[,] lllnc_out)
        {
            int shomwhile = 0;
            double[,] S = new double[(6 * (image_name_Eterior.Length)) + 10, (6 * (image_name_Eterior.Length)) + 10];
            double[,] llln = new double[6 * (image_name_Eterior.Length), 1];
            double[,] lllnc = new double[10, 1];
            double[,] eawveb = new double[6 * (image_name_Eterior.Length) + 10, 1];
            double[,] pp2nC = new double[10, 10];
            pp2nC[0, 0] = 1 / 10000000000000000000;
            pp2nC[1, 1] = 1 / 10000000000000000000;
            pp2nC[2, 2] = 1 / 10000000000000000000;
            pp2nC[3, 3] = 1 / 10000000000000000000;
            pp2nC[4, 4] = 1 / 10000000000000000000;
            pp2nC[5, 5] = 1 / 10000000000000000000;
            pp2nC[6, 6] = 1 / 10000000000000000000;
            pp2nC[7, 7] = 1 / 10000000000000000000;
            pp2nC[8, 8] = 1 / 10000000000000000000;
            pp2nC[9, 9] = 1 / 10000000000000000000;
            for (int i1 = 0; i1 < ii; i1++)
            {
                shomwhile += count_camera[i1];
            }
            for (int i = ii; i < ie; i++)
            {

                double[,] v = new double[3, 3];
                double[,] v2 = new double[3, 3];
                double[,] ww = new double[6, 3];
                double[,] wwt = new double[3, 3];
                double[,] u = new double[6, 6];
                double[,] Fc = new double[10, 3];
                double Fc00, Fc01, Fc02, Fc10, Fc11, Fc12, Fc20, Fc21, Fc22, Fc30, Fc31, Fc32, Fc40, Fc41, Fc42, Fc50, Fc51, Fc52, Fc60, Fc61, Fc62, Fc70, Fc71, Fc72, Fc80, Fc81, Fc82, Fc90, Fc91, Fc92;
                List<double> w1 = new List<double>();
                List<double> w2 = new List<double>();
                List<double> w3 = new List<double>();
                List<int> wj = new List<int>();
                int ij2 = shomwhile;
                shomwhile += count_camera[i];

                double ff = IOP[0, 0];
                double x0 = IOP[1, 0];
                double y0 = IOP[2, 0];
                double K1 = IOP[3, 0];
                double K2 = IOP[4, 0];
                double K3 = IOP[5, 0];
                double P1 = IOP[6, 0];
                double P2 = IOP[7, 0];
                double B1 = IOP[8, 0];
                double B2 = IOP[9, 0];

             
                lllnc[0, 0] = ff - interior[0];
                lllnc[1, 0] = x0 - interior[1];
                lllnc[2, 0] = y0 - interior[2];
                lllnc[3, 0] = K1 - interior[3];
                lllnc[4, 0] = K2 - interior[4];
                lllnc[5, 0] = K3 - interior[5];
                lllnc[6, 0] = P1 - interior[6];
                lllnc[7, 0] = P2 - interior[7];
                lllnc[8, 0] = B1 - interior[8];
                lllnc[9, 0] = B2 - interior[9];
                if (point_id[i] != "non")
                {
                    double[,] pp2nB = new double[3, 3];

                    pp2nB[0, 0] = observe_value_points_weight[i, 0];
                    pp2nB[1, 1] = observe_value_points_weight[i, 1];
                    pp2nB[2, 2] = observe_value_points_weight[i, 2];

                    double[,] ebii = new double[3, 1];
                    double[,] lllnb = new double[3, 1];
                    double[,] pp2 = new double[2, 2];
                    double[,] pp2nA = new double[6, 6];




                    for (int ij = ij2; ij < shomwhile; ij++)
                    {

                        int j2 = Convert.ToInt32(cam_observe[ij, 0]);

                        pp2[0, 0] = camera_weight[8 * j2];
                        pp2[1, 1] = camera_weight[8 * j2 + 1];
                        pp2nA[0, 0] = camera_weight[8 * j2 + 2];
                        pp2nA[1, 1] = camera_weight[8 * j2 + 3];
                        pp2nA[2, 2] = camera_weight[8 * j2 + 4];
                        pp2nA[3, 3] = camera_weight[8 * j2 + 5];
                        pp2nA[4, 4] = camera_weight[8 * j2 + 6];
                        pp2nA[5, 5] = camera_weight[8 * j2 + 7];



                            double XL = Exterior_Orientation_Eterior[j2, 0];
                            double YL = Exterior_Orientation_Eterior[j2, 1];
                            double ZL = Exterior_Orientation_Eterior[j2, 2];
                            double om = Exterior_Orientation_Eterior[j2, 3];
                            double phi = Exterior_Orientation_Eterior[j2, 4];
                            double k = Exterior_Orientation_Eterior[j2, 5];

                            double X = Initial_value[i, 0];
                            double Y = Initial_value[i, 1];
                            double Z = Initial_value[i, 2];
                        double m11 = image_m[j2, 0];
                        double m12 = image_m[j2, 1];
                        //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                        double m13 = image_m[j2, 2];
                        //m21=-cos(phi)*sin(k)
                        double m21 = image_m[j2, 3];
                        //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                        double m22 = image_m[j2, 4];
                        //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                        double m23 = image_m[j2, 5];
                        //m31=sin(phi)
                        double m31 = image_m[j2, 6];
                        //m32=-sin(om)*cos(phi)
                        double m32 = image_m[j2, 7];
                        //m33=cos(om)*cos(phi)
                        double m33 = image_m[j2, 8];

                        double sinom = image_sin_cos[j2, 0]; double sinphi = image_sin_cos[j2, 1]; double sink = image_sin_cos[j2, 2];
                        double cosom = image_sin_cos[j2, 3]; double cosphi = image_sin_cos[j2, 4]; double cosk = image_sin_cos[j2, 5];


                        double xx = cam_observe[ij, 1];
                        double yy = cam_observe[ij, 2];
                            // Fx_xx = 2 * P2 * (y0 - yy) - K1 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2)                             - B1 - K2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2                                                                                   - K3 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 3 + P1 * (6 * x0 - 6 * xx) - (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) + 3 * K3 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2) + 1
                        double Fx_xx = 2 * P2 * (y0 - yy) - K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - B1 - K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + P1 * (6 * x0 - 6 * xx) - (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) + 1;
                            double Fy_xx = 2 * P1 * (y0 - yy) + P2 * (2 * x0 - 2 * xx) - (y0 - yy) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))));
                            double Fx_yy = 2 * P2 * (x0 - xx) - B2 + P1 * (2 * y0 - 2 * yy) - (x0 - xx) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))));
                            double Fy_yy = 2 * P1 * (x0 - xx) - K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + P2 * (6 * y0 - 6 * yy) - (y0 - yy) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) + 1;



                            double xom = -((ff * ((sink * sinom - cosk * cosom * sinphi) * (Y - YL) - (cosom * sink + cosk * sinom * sinphi) * (Z - ZL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));

                            double xphi = -((ff * (cosk * sinphi * (X - XL) + cosk * cosom * cosphi * (Z - ZL) - cosk * cosphi * sinom * (Y - YL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) + (ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                            double xk = -(-(ff * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double xXL = -((ff * cosk * cosphi) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * sinphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                            double xYL = -((ff * (cosom * sink + cosk * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) + (ff * cosphi * sinom * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                            double xZL = -(-((ff * cosom * cosphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * (sink * sinom - cosk * cosom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))));
                            //   (f * cos(om) * cos(phi) * ((cos(om) * sin(k) + cos(k) * sin(om) * sin(phi)) * (Y - YL) + (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi)) * (Z - ZL) + cos(k) * cos(phi) * (X - XL))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL)) ^ 2 - (f * (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL))];
                            double yom = -((ff * ((cosk * sinom + cosom * sink * sinphi) * (Y - YL) - (cosk * cosom - sink * sinom * sinphi) * (Z - ZL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                            double yphi = -((ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2) - (ff * (sink * sinphi * (X - XL) + cosom * cosphi * sink * (Z - ZL) - cosphi * sink * sinom * (Y - YL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double yk = -((ff * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double yXL = -(-(ff * sinphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2) - (ff * cosphi * sink) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double yYL = -((ff * (cosk * cosom - sink * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) + (ff * cosphi * sinom * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                            double yZL = -((ff * (cosk * sinom + cosom * sink * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * cosom * cosphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                            double xX = -((ff * sinphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * cosk * cosphi) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double xY = -(-(ff * (cosom * sink + cosk * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * cosphi * sinom * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))));
                            double xZ = -((ff * cosom * cosphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * (sink * sinom - cosk * cosom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double yX = -((ff * sinphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) + (ff * cosphi * sink) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                            double yY = -(-(ff * (cosk * cosom - sink * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * cosphi * sinom * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))));
                            double yZ = -((ff * cosom * cosphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * (cosk * sinom + cosom * sink * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));


                            double xpi = (xx - x0);
                            double ypi = (yy - y0);

                            double ri = Math.Sqrt(((xpi * xpi) + (ypi * ypi)));
                            double deltax = (xpi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P1 * (((ri * ri)) + (2 * ((xpi * xpi))))) + (2 * P2 * xpi * ypi) + (B1 * xpi) + (B2 * ypi);
                            double deltay = (ypi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P2 * (((ri * ri)) + (2 * ((ypi * ypi))))) + (2 * P1 * xpi * ypi);
                            if (i == 1 && j2 == 0)
                            {
                                //       MessageBox.Show("" + deltax + "  " + deltay);
                            }
                            double x = -(ff * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltax;
                            double y = -(ff * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltay;
                            double[,] lll = new double[2, 1];
                            double Fx = (xpi) + (ff * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) - (deltax);
                            double Fy = (ypi) + (ff * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) - (deltay);
                            //    MessageBox.Show("" + deltax + "   " + deltay);
                            double Fxf = (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                            double Fyf = (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                        //Fxx0= B1 + K1 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) - 2 * P2 * (y0 - yy) + K2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2 + K3 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 3 -                                                                                                                                                                                                                                                                      P1 * (6 * x0 - 6 * xx) + (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2)                              + 3 * K3 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2) - 1
        
                        double Fxx0 = B1 + K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - 2 * P2 * (y0 - yy) + K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - P1 * (6 * x0 - 6 * xx) + (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) - 1;

                        //      Fxy0= B2 - 2 * P2 * (x0 - xx) - P1 * (2 * y0 - 2 * yy) + (x0 - xx) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2)                             + 3 * K3 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2)
                        double Fxy0 = B2 - 2 * P2 * (x0 - xx) - P1 * (2 * y0 - 2 * yy) + (x0 - xx) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))));
                        //FxK1=  ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) * (x0 - xx)
                        double FxK1 = ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))* (x0 - xx);
                        //FxK2 =   ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2 * (x0 - xx)
                        double FxK2 = ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * (x0 - xx);
                        //FxK3 =((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 3 * (x0 - xx)
                        double FxK3 = ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))* ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))* ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * (x0 - xx);
                       // FxP1 = -3 * (x0 - xx) ^ 2 - (y0 - yy) ^ 2
                        double FxP1 = ((x0 - xx) * (x0 - xx))- ((y0 - yy) * (y0 - yy));
                        //FxP2 = -2 * (x0 - xx) * (y0 - yy)
                            double FxP2 = -2 * xpi * ypi;

                        double FxB1 = -xpi;

                        double FxB2 = -ypi;
                     //   Fyx0 =          (y0 - yy) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2)                             + 3 * K3 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2)                                                                                   - P2 * (2 * x0 - 2 * xx) - 2 * P1 * (y0 - yy)
                            double Fyx0 = (y0 - yy) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) - P2 * (2 * x0 - 2 * xx) - 2 * P1 * (y0 - yy);
                        //Fyy0 =      K1 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2)                             - 2 * P1 * (x0 - xx) + K2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2 +                                                                                   K3 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 3                                                                                                                                             - P2 * (6 * y0 - 6 * yy) + (y0 - yy) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2)                             + 3 * K3 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2)                                                                                   - 1
                        double Fyy0 = K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - 2 * P1 * (x0 - xx) + K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - P2 * (6 * y0 - 6 * yy) + (y0 - yy) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) - 1;

                        // double FyK1 = ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) * (y0 - yy)
                        double FyK1 = -ypi * (ri * ri);
                        //FyK2 =   ((x0 - xx)^2 + (y0 - yy)^2)^2*(y0 - yy)
                        double FyK2 = -ypi * (ri * ri) * (ri * ri);
                        // double FyK3 = ((x0 - xx)^2 + (y0 - yy)^2)^3*(y0 - yy)
                        double FyK3 = -ypi * (ri * ri) * (ri * ri) * (ri * ri);
                        // double FyP1 = -2*(x0 - xx)*(y0 - yy)
                        double FyP1 = -2 * ypi * xpi;
                        // double FyP2 = - (x0 - xx)^2 - 3*(y0 - yy)^2
                        double FyP2 = -(((ri * ri)) + (2 * (ypi * ypi)));
                            double[,] AA = new double[2, 6];
                            double[,] BB = new double[2, 3];
                            double[,] c = new double[2, 10];


                            double b00 = Fx_xx; double b01 = Fx_yy; double b10 = Fy_xx; double b11 = Fy_yy;

                            double p00 = 1.0 / pp2[0, 0]; double p11 = 1.0 / pp2[1, 1];
                            double detbpb = (b00 * b11 * p00 * p11 * (b00) * (b11) - b00 * b11 * p00 * p11 * (b01) * (b10) - b01 * b10 * p00 * p11 * (b00) * (b11) + b01 * b10 * p00 * p11 * (b01) * (b10));
                        //  double[,] PB = Matrix.Inverse(Matrix.Multiply(Matrix.Multiply(B, PP), Matrix.Transpose(B)));
                  
                        double[,] PB = new double[2, 2];
                            PB[0, 0] = (b10 * p00 * (b10) + b11 * p11 * (b11)) / detbpb;
                            double pb00 = (b10 * p00 * (b10) + b11 * p11 * (b11)) / detbpb;
                            PB[0, 1] = -(b00 * p00 * (b10) + b01 * p11 * (b11)) / detbpb;
                            double pb01 = -(b00 * p00 * (b10) + b01 * p11 * (b11)) / detbpb;
                            PB[1, 0] = -(b10 * p00 * (b00) + b11 * p11 * (b01)) / detbpb;
                            double pb10 = -(b10 * p00 * (b00) + b11 * p11 * (b01)) / detbpb;
                            PB[1, 1] = (b00 * p00 * (b00) + b01 * p11 * (b01)) / detbpb;
                            double pb11 = (b00 * p00 * (b00) + b01 * p11 * (b01)) / detbpb;
                            //MessageBox.Show("" +PB[0, 0] + "  " + PB[1, 1] + "  " + PP[1, 0] + "  " + PB[0, 1]);

                            AA[0, 0] = xom;
                            double a00 = xom;
                            AA[0, 1] = xphi;
                            double a01 = xphi;
                            AA[0, 2] = xk;
                            double a02 = xk;
                            AA[0, 3] = xXL;
                            double a03 = xXL;
                            AA[0, 4] = xYL;
                            double a04 = xYL;
                            AA[0, 5] = xZL;
                            double a05 = xZL;
                            AA[1, 0] = yom;
                            double a10 = yom;
                            AA[1, 1] = yphi;
                            double a11 = yphi;
                            AA[1, 2] = yk;
                            double a12 = yk;
                            AA[1, 3] = yXL;
                            double a13 = yXL;
                            AA[1, 4] = yYL;
                            double a14 = yYL;
                            AA[1, 5] = yZL;
                            double a15 = yZL;
                            BB[0, 0] = xX;
                            double bb00 = xX;
                            BB[0, 1] = xY;
                            double bb01 = xY;
                            BB[0, 2] = xZ;
                            double bb02 = xZ;
                            BB[1, 0] = yX;
                            double bb10 = yX;
                            BB[1, 1] = yY;
                            double bb11 = yY;
                            BB[1, 2] = yZ;
                            double bb12 = yZ;
                            c[0, 0] = Fxx0;
                            double c00 = Fxx0;
                            c[0, 1] = Fxy0;
                            double c01 = Fxy0;
                            c[0, 2] = FxK1;
                            double c02 = FxK1;
                            c[0, 3] = FxK2;
                            double c03 = FxK2;
                            c[0, 4] = FxK3;
                            double c04 = FxK3;
                            c[0, 5] = FxP1;
                            double c05 = FxP1;
                            c[0, 6] = FxP2;
                            double c06 = FxP2;
                            c[0, 7] = FxB1;
                            double c07 = FxB1;
                            c[0, 8] = FxB2;
                            double c08 = FxB2;
                            c[0, 9] = Fxf;
                            double c09 = Fxf;
                            c[1, 0] = Fyx0;
                            double c10 = Fyx0;
                            c[1, 1] = Fyy0;
                            double c11 = Fyy0;
                            c[1, 2] = FyK1;
                            double c12 = FyK1;
                            c[1, 3] = FyK2;
                            double c13 = FyK2;
                            c[1, 4] = FyK3;
                            double c14 = FyK3;
                            c[1, 5] = FyP1;
                            double c15 = FyP1;
                            c[1, 6] = FyP2;
                            double c16 = FyP2;
                            c[1, 7] = 0;
                            double c17 = 0;
                            c[1, 8] = 0;
                            double c18 = 0;
                            c[1, 9] = Fyf;
                            double c19 = Fyf;
                              lll[0, 0] = Fx; lll[1, 0] = Fy;
                            double ll0 = Fx; double ll1 = Fy;
                     //    a21*(((a21)*(b12*p11*(b12) + b11*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b12*p11*(b22) + b11*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))) - a11*(((a21)*(b22*p11*(b12) + b21*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b22*p11*(b22) + b21*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))), a22*(((a21)*(b12*p11*(b12) + b11*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b12*p11*(b22) + b11*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))) - a12*(((a21)*(b22*p11*(b12) + b21*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b22*p11*(b22) + b21*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))), a23*(((a21)*(b12*p11*(b12) + b11*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b12*p11*(b22) + b11*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))) - a13*(((a21)*(b22*p11*(b12) + b21*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b22*p11*(b22) + b21*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))), a24*(((a21)*(b12*p11*(b12) + b11*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b12*p11*(b22) + b11*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))) - a14*(((a21)*(b22*p11*(b12) + b21*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b22*p11*(b22) + b21*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))), a25*(((a21)*(b12*p11*(b12) + b11*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b12*p11*(b22) + b11*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))) - a15*(((a21)*(b22*p11*(b12) + b21*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b22*p11*(b22) + b21*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))), a26*(((a21)*(b12*p11*(b12) + b11*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b12*p11*(b22) + b11*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21))) - a16*(((a21)*(b22*p11*(b12) + b21*p22*(b11)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)) - ((a11)*(b22*p11*(b22) + b21*p22*(b21)))/(b11*b22*(b11)*(b22) - b11*b22*(b12)*(b21) - b12*b21*(b11)*(b22) + b12*b21*(b12)*(b21)))
                        llln[6 * j2, 0] = observe_value_camera[j2, 0] - om; llln[6 * j2 + 1, 0] = observe_value_camera[j2, 1] - phi; llln[6 * j2 + 2, 0] = observe_value_camera[j2, 2] - k;
                            llln[6 * j2 + 3, 0] = observe_value_camera[j2, 3] - XL; llln[6 * j2 + 4, 0] = observe_value_camera[j2, 4] - YL; llln[6 * j2 + 5, 0] = observe_value_camera[j2, 5] - ZL; lllnb[0, 0] = observe_value_points[i, 0] - X; lllnb[1, 0] = observe_value_points[i, 1] - Y; lllnb[2, 0] = observe_value_points[i, 2] - Z;
                             double[,] Ec = new double[6, 10];
                            Ec[0, 0] = c00 * (pb00 * (a00) + pb10 * (a10)) + c10 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 1] = c01 * (pb00 * (a00) + pb10 * (a10)) + c11 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 2] = c02 * (pb00 * (a00) + pb10 * (a10)) + c12 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 3] = c03 * (pb00 * (a00) + pb10 * (a10)) + c13 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 4] = c04 * (pb00 * (a00) + pb10 * (a10)) + c14 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 5] = c05 * (pb00 * (a00) + pb10 * (a10)) + c15 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 6] = c06 * (pb00 * (a00) + pb10 * (a10)) + c16 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 7] = c07 * (pb00 * (a00) + pb10 * (a10)) + c17 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 8] = c08 * (pb00 * (a00) + pb10 * (a10)) + c18 * (pb01 * (a00) + pb11 * (a10)); Ec[0, 9] = c09 * (pb00 * (a00) + pb10 * (a10)) + c19 * (pb01 * (a00) + pb11 * (a10));
                            Ec[1, 0] = c00 * (pb00 * (a01) + pb10 * (a11)) + c10 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 1] = c01 * (pb00 * (a01) + pb10 * (a11)) + c11 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 2] = c02 * (pb00 * (a01) + pb10 * (a11)) + c12 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 3] = c03 * (pb00 * (a01) + pb10 * (a11)) + c13 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 4] = c04 * (pb00 * (a01) + pb10 * (a11)) + c14 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 5] = c05 * (pb00 * (a01) + pb10 * (a11)) + c15 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 6] = c06 * (pb00 * (a01) + pb10 * (a11)) + c16 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 7] = c07 * (pb00 * (a01) + pb10 * (a11)) + c17 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 8] = c08 * (pb00 * (a01) + pb10 * (a11)) + c18 * (pb01 * (a01) + pb11 * (a11)); Ec[1, 9] = c09 * (pb00 * (a01) + pb10 * (a11)) + c19 * (pb01 * (a01) + pb11 * (a11));
                            Ec[2, 0] = c00 * (pb00 * (a02) + pb10 * (a12)) + c10 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 1] = c01 * (pb00 * (a02) + pb10 * (a12)) + c11 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 2] = c02 * (pb00 * (a02) + pb10 * (a12)) + c12 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 3] = c03 * (pb00 * (a02) + pb10 * (a12)) + c13 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 4] = c04 * (pb00 * (a02) + pb10 * (a12)) + c14 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 5] = c05 * (pb00 * (a02) + pb10 * (a12)) + c15 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 6] = c06 * (pb00 * (a02) + pb10 * (a12)) + c16 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 7] = c07 * (pb00 * (a02) + pb10 * (a12)) + c17 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 8] = c08 * (pb00 * (a02) + pb10 * (a12)) + c18 * (pb01 * (a02) + pb11 * (a12)); Ec[2, 9] = c09 * (pb00 * (a02) + pb10 * (a12)) + c19 * (pb01 * (a02) + pb11 * (a12));
                            Ec[3, 0] = c00 * (pb00 * (a03) + pb10 * (a13)) + c10 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 1] = c01 * (pb00 * (a03) + pb10 * (a13)) + c11 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 2] = c02 * (pb00 * (a03) + pb10 * (a13)) + c12 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 3] = c03 * (pb00 * (a03) + pb10 * (a13)) + c13 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 4] = c04 * (pb00 * (a03) + pb10 * (a13)) + c14 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 5] = c05 * (pb00 * (a03) + pb10 * (a13)) + c15 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 6] = c06 * (pb00 * (a03) + pb10 * (a13)) + c16 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 7] = c07 * (pb00 * (a03) + pb10 * (a13)) + c17 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 8] = c08 * (pb00 * (a03) + pb10 * (a13)) + c18 * (pb01 * (a03) + pb11 * (a13)); Ec[3, 9] = c09 * (pb00 * (a03) + pb10 * (a13)) + c19 * (pb01 * (a03) + pb11 * (a13));
                            Ec[4, 0] = c00 * (pb00 * (a04) + pb10 * (a14)) + c10 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 1] = c01 * (pb00 * (a04) + pb10 * (a14)) + c11 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 2] = c02 * (pb00 * (a04) + pb10 * (a14)) + c12 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 3] = c03 * (pb00 * (a04) + pb10 * (a14)) + c13 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 4] = c04 * (pb00 * (a04) + pb10 * (a14)) + c14 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 5] = c05 * (pb00 * (a04) + pb10 * (a14)) + c15 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 6] = c06 * (pb00 * (a04) + pb10 * (a14)) + c16 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 7] = c07 * (pb00 * (a04) + pb10 * (a14)) + c17 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 8] = c08 * (pb00 * (a04) + pb10 * (a14)) + c18 * (pb01 * (a04) + pb11 * (a14)); Ec[4, 9] = c09 * (pb00 * (a04) + pb10 * (a14)) + c19 * (pb01 * (a04) + pb11 * (a14));
                            Ec[5, 0] = c00 * (pb00 * (a05) + pb10 * (a15)) + c10 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 1] = c01 * (pb00 * (a05) + pb10 * (a15)) + c11 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 2] = c02 * (pb00 * (a05) + pb10 * (a15)) + c12 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 3] = c03 * (pb00 * (a05) + pb10 * (a15)) + c13 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 4] = c04 * (pb00 * (a05) + pb10 * (a15)) + c14 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 5] = c05 * (pb00 * (a05) + pb10 * (a15)) + c15 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 6] = c06 * (pb00 * (a05) + pb10 * (a15)) + c16 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 7] = c07 * (pb00 * (a05) + pb10 * (a15)) + c17 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 8] = c08 * (pb00 * (a05) + pb10 * (a15)) + c18 * (pb01 * (a05) + pb11 * (a15)); Ec[5, 9] = c09 * (pb00 * (a05) + pb10 * (a15)) + c19 * (pb01 * (a05) + pb11 * (a15));
                        

                     
                        

                        double[,] Ect = new double[10, 6];
                            Ect[0, 0] = Ec[0, 0]; Ect[1, 0] = Ec[0, 1]; Ect[2, 0] = Ec[0, 2]; Ect[3, 0] = Ec[0, 3]; Ect[4, 0] = Ec[0, 4]; Ect[5, 0] = Ec[0, 5]; Ect[6, 0] = Ec[0, 6]; Ect[7, 0] = Ec[0, 7]; Ect[8, 0] = Ec[0, 8]; Ect[9, 0] = Ec[0, 9];
                            Ect[0, 1] = Ec[1, 0]; Ect[1, 1] = Ec[1, 1]; Ect[2, 1] = Ec[1, 2]; Ect[3, 1] = Ec[1, 3]; Ect[4, 1] = Ec[1, 4]; Ect[5, 1] = Ec[1, 5]; Ect[6, 1] = Ec[1, 6]; Ect[7, 1] = Ec[1, 7]; Ect[8, 1] = Ec[1, 8]; Ect[9, 1] = Ec[1, 9];
                            Ect[0, 2] = Ec[2, 0]; Ect[1, 2] = Ec[2, 1]; Ect[2, 2] = Ec[2, 2]; Ect[3, 2] = Ec[2, 3]; Ect[4, 2] = Ec[2, 4]; Ect[5, 2] = Ec[2, 5]; Ect[6, 2] = Ec[2, 6]; Ect[7, 2] = Ec[2, 7]; Ect[8, 2] = Ec[2, 8]; Ect[9, 2] = Ec[2, 9];
                            Ect[0, 3] = Ec[3, 0]; Ect[1, 3] = Ec[3, 1]; Ect[2, 3] = Ec[3, 2]; Ect[3, 3] = Ec[3, 3]; Ect[4, 3] = Ec[3, 4]; Ect[5, 3] = Ec[3, 5]; Ect[6, 3] = Ec[3, 6]; Ect[7, 3] = Ec[3, 7]; Ect[8, 3] = Ec[3, 8]; Ect[9, 3] = Ec[3, 9];
                            Ect[0, 4] = Ec[4, 0]; Ect[1, 4] = Ec[4, 1]; Ect[2, 4] = Ec[4, 2]; Ect[3, 4] = Ec[4, 3]; Ect[4, 4] = Ec[4, 4]; Ect[5, 4] = Ec[4, 5]; Ect[6, 4] = Ec[4, 6]; Ect[7, 4] = Ec[4, 7]; Ect[8, 4] = Ec[4, 8]; Ect[9, 4] = Ec[4, 9];
                            Ect[0, 5] = Ec[5, 0]; Ect[1, 5] = Ec[5, 1]; Ect[2, 5] = Ec[5, 2]; Ect[3, 5] = Ec[5, 3]; Ect[4, 5] = Ec[5, 4]; Ect[5, 5] = Ec[5, 5]; Ect[6, 5] = Ec[5, 6]; Ect[7, 5] = Ec[5, 7]; Ect[8, 5] = Ec[5, 8]; Ect[9, 5] = Ec[5, 9];

                            double[,] Fc2 = new double[10, 3];
                            Fc2[0, 0] = bb00 * (pb00 * (c00) + pb10 * (c10)) + bb10 * (pb01 * (c00) + pb11 * (c10)); Fc2[0, 1] = bb01 * (pb00 * (c00) + pb10 * (c10)) + bb11 * (pb01 * (c00) + pb11 * (c10)); Fc2[0, 2] = bb02 * (pb00 * (c00) + pb10 * (c10)) + bb12 * (pb01 * (c00) + pb11 * (c10));
                            Fc2[1, 0] = bb00 * (pb00 * (c01) + pb10 * (c11)) + bb10 * (pb01 * (c01) + pb11 * (c11)); Fc2[1, 1] = bb01 * (pb00 * (c01) + pb10 * (c11)) + bb11 * (pb01 * (c01) + pb11 * (c11)); Fc2[1, 2] = bb02 * (pb00 * (c01) + pb10 * (c11)) + bb12 * (pb01 * (c01) + pb11 * (c11));
                            Fc2[2, 0] = bb00 * (pb00 * (c02) + pb10 * (c12)) + bb10 * (pb01 * (c02) + pb11 * (c12)); Fc2[2, 1] = bb01 * (pb00 * (c02) + pb10 * (c12)) + bb11 * (pb01 * (c02) + pb11 * (c12)); Fc2[2, 2] = bb02 * (pb00 * (c02) + pb10 * (c12)) + bb12 * (pb01 * (c02) + pb11 * (c12));
                            Fc2[3, 0] = bb00 * (pb00 * (c03) + pb10 * (c13)) + bb10 * (pb01 * (c03) + pb11 * (c13)); Fc2[3, 1] = bb01 * (pb00 * (c03) + pb10 * (c13)) + bb11 * (pb01 * (c03) + pb11 * (c13)); Fc2[3, 2] = bb02 * (pb00 * (c03) + pb10 * (c13)) + bb12 * (pb01 * (c03) + pb11 * (c13));
                            Fc2[4, 0] = bb00 * (pb00 * (c04) + pb10 * (c14)) + bb10 * (pb01 * (c04) + pb11 * (c14)); Fc2[4, 1] = bb01 * (pb00 * (c04) + pb10 * (c14)) + bb11 * (pb01 * (c04) + pb11 * (c14)); Fc2[4, 2] = bb02 * (pb00 * (c04) + pb10 * (c14)) + bb12 * (pb01 * (c04) + pb11 * (c14));
                            Fc2[5, 0] = bb00 * (pb00 * (c05) + pb10 * (c15)) + bb10 * (pb01 * (c05) + pb11 * (c15)); Fc2[5, 1] = bb01 * (pb00 * (c05) + pb10 * (c15)) + bb11 * (pb01 * (c05) + pb11 * (c15)); Fc2[5, 2] = bb02 * (pb00 * (c05) + pb10 * (c15)) + bb12 * (pb01 * (c05) + pb11 * (c15));
                            Fc2[6, 0] = bb00 * (pb00 * (c06) + pb10 * (c16)) + bb10 * (pb01 * (c06) + pb11 * (c16)); Fc2[6, 1] = bb01 * (pb00 * (c06) + pb10 * (c16)) + bb11 * (pb01 * (c06) + pb11 * (c16)); Fc2[6, 2] = bb02 * (pb00 * (c06) + pb10 * (c16)) + bb12 * (pb01 * (c06) + pb11 * (c16));
                            Fc2[7, 0] = bb00 * (pb00 * (c07) + pb10 * (c17)) + bb10 * (pb01 * (c07) + pb11 * (c17)); Fc2[7, 1] = bb01 * (pb00 * (c07) + pb10 * (c17)) + bb11 * (pb01 * (c07) + pb11 * (c17)); Fc2[7, 2] = bb02 * (pb00 * (c07) + pb10 * (c17)) + bb12 * (pb01 * (c07) + pb11 * (c17));
                            Fc2[8, 0] = bb00 * (pb00 * (c08) + pb10 * (c18)) + bb10 * (pb01 * (c08) + pb11 * (c18)); Fc2[8, 1] = bb01 * (pb00 * (c08) + pb10 * (c18)) + bb11 * (pb01 * (c08) + pb11 * (c18)); Fc2[8, 2] = bb02 * (pb00 * (c08) + pb10 * (c18)) + bb12 * (pb01 * (c08) + pb11 * (c18));
                            Fc2[9, 0] = bb00 * (pb00 * (c09) + pb10 * (c19)) + bb10 * (pb01 * (c09) + pb11 * (c19)); Fc2[9, 1] = bb01 * (pb00 * (c09) + pb10 * (c19)) + bb11 * (pb01 * (c09) + pb11 * (c19)); Fc2[9, 2] = bb02 * (pb00 * (c09) + pb10 * (c19)) + bb12 * (pb01 * (c09) + pb11 * (c19));
                     
                        double[,] FcV2 = new double[3, 10];
                            FcV2[0, 0] = c00 * (pb00 * (bb00) + pb10 * (bb10)) + c10 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 1] = c01 * (pb00 * (bb00) + pb10 * (bb10)) + c11 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 2] = c02 * (pb00 * (bb00) + pb10 * (bb10)) + c12 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 3] = c03 * (pb00 * (bb00) + pb10 * (bb10)) + c13 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 4] = c04 * (pb00 * (bb00) + pb10 * (bb10)) + c14 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 5] = c05 * (pb00 * (bb00) + pb10 * (bb10)) + c15 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 6] = c06 * (pb00 * (bb00) + pb10 * (bb10)) + c16 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 7] = c07 * (pb00 * (bb00) + pb10 * (bb10)) + c17 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 8] = c08 * (pb00 * (bb00) + pb10 * (bb10)) + c18 * (pb01 * (bb00) + pb11 * (bb10)); FcV2[0, 9] = c09 * (pb00 * (bb00) + pb10 * (bb10)) + c19 * (pb01 * (bb00) + pb11 * (bb10));
                            FcV2[1, 0] = c00 * (pb00 * (bb01) + pb10 * (bb11)) + c10 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 1] = c01 * (pb00 * (bb01) + pb10 * (bb11)) + c11 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 2] = c02 * (pb00 * (bb01) + pb10 * (bb11)) + c12 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 3] = c03 * (pb00 * (bb01) + pb10 * (bb11)) + c13 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 4] = c04 * (pb00 * (bb01) + pb10 * (bb11)) + c14 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 5] = c05 * (pb00 * (bb01) + pb10 * (bb11)) + c15 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 6] = c06 * (pb00 * (bb01) + pb10 * (bb11)) + c16 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 7] = c07 * (pb00 * (bb01) + pb10 * (bb11)) + c17 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 8] = c08 * (pb00 * (bb01) + pb10 * (bb11)) + c18 * (pb01 * (bb01) + pb11 * (bb11)); FcV2[1, 9] = c09 * (pb00 * (bb01) + pb10 * (bb11)) + c19 * (pb01 * (bb01) + pb11 * (bb11));
                            FcV2[2, 0] = c00 * (pb00 * (bb02) + pb10 * (bb12)) + c10 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 1] = c01 * (pb00 * (bb02) + pb10 * (bb12)) + c11 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 2] = c02 * (pb00 * (bb02) + pb10 * (bb12)) + c12 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 3] = c03 * (pb00 * (bb02) + pb10 * (bb12)) + c13 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 4] = c04 * (pb00 * (bb02) + pb10 * (bb12)) + c14 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 5] = c05 * (pb00 * (bb02) + pb10 * (bb12)) + c15 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 6] = c06 * (pb00 * (bb02) + pb10 * (bb12)) + c16 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 7] = c07 * (pb00 * (bb02) + pb10 * (bb12)) + c17 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 8] = c08 * (pb00 * (bb02) + pb10 * (bb12)) + c18 * (pb01 * (bb02) + pb11 * (bb12)); FcV2[2, 9] = c09 * (pb00 * (bb02) + pb10 * (bb12)) + c19 * (pb01 * (bb02) + pb11 * (bb12));

                         
                            double[,] Pc = new double[10, 10];
                            Pc[0, 0] = c00 * (pb00 * (c00) + pb10 * (c10)) + c10 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 1] = c01 * (pb00 * (c00) + pb10 * (c10)) + c11 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 2] = c02 * (pb00 * (c00) + pb10 * (c10)) + c12 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 3] = c03 * (pb00 * (c00) + pb10 * (c10)) + c13 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 4] = c04 * (pb00 * (c00) + pb10 * (c10)) + c14 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 5] = c05 * (pb00 * (c00) + pb10 * (c10)) + c15 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 6] = c06 * (pb00 * (c00) + pb10 * (c10)) + c16 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 7] = c07 * (pb00 * (c00) + pb10 * (c10)) + c17 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 8] = c08 * (pb00 * (c00) + pb10 * (c10)) + c18 * (pb01 * (c00) + pb11 * (c10)); Pc[0, 9] = c09 * (pb00 * (c00) + pb10 * (c10)) + c19 * (pb01 * (c00) + pb11 * (c10));
                            Pc[1, 0] = c00 * (pb00 * (c01) + pb10 * (c11)) + c10 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 1] = c01 * (pb00 * (c01) + pb10 * (c11)) + c11 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 2] = c02 * (pb00 * (c01) + pb10 * (c11)) + c12 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 3] = c03 * (pb00 * (c01) + pb10 * (c11)) + c13 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 4] = c04 * (pb00 * (c01) + pb10 * (c11)) + c14 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 5] = c05 * (pb00 * (c01) + pb10 * (c11)) + c15 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 6] = c06 * (pb00 * (c01) + pb10 * (c11)) + c16 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 7] = c07 * (pb00 * (c01) + pb10 * (c11)) + c17 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 8] = c08 * (pb00 * (c01) + pb10 * (c11)) + c18 * (pb01 * (c01) + pb11 * (c11)); Pc[1, 9] = c09 * (pb00 * (c01) + pb10 * (c11)) + c19 * (pb01 * (c01) + pb11 * (c11));
                            Pc[2, 0] = c00 * (pb00 * (c02) + pb10 * (c12)) + c10 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 1] = c01 * (pb00 * (c02) + pb10 * (c12)) + c11 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 2] = c02 * (pb00 * (c02) + pb10 * (c12)) + c12 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 3] = c03 * (pb00 * (c02) + pb10 * (c12)) + c13 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 4] = c04 * (pb00 * (c02) + pb10 * (c12)) + c14 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 5] = c05 * (pb00 * (c02) + pb10 * (c12)) + c15 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 6] = c06 * (pb00 * (c02) + pb10 * (c12)) + c16 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 7] = c07 * (pb00 * (c02) + pb10 * (c12)) + c17 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 8] = c08 * (pb00 * (c02) + pb10 * (c12)) + c18 * (pb01 * (c02) + pb11 * (c12)); Pc[2, 9] = c09 * (pb00 * (c02) + pb10 * (c12)) + c19 * (pb01 * (c02) + pb11 * (c12));
                            Pc[3, 0] = c00 * (pb00 * (c03) + pb10 * (c13)) + c10 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 1] = c01 * (pb00 * (c03) + pb10 * (c13)) + c11 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 2] = c02 * (pb00 * (c03) + pb10 * (c13)) + c12 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 3] = c03 * (pb00 * (c03) + pb10 * (c13)) + c13 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 4] = c04 * (pb00 * (c03) + pb10 * (c13)) + c14 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 5] = c05 * (pb00 * (c03) + pb10 * (c13)) + c15 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 6] = c06 * (pb00 * (c03) + pb10 * (c13)) + c16 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 7] = c07 * (pb00 * (c03) + pb10 * (c13)) + c17 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 8] = c08 * (pb00 * (c03) + pb10 * (c13)) + c18 * (pb01 * (c03) + pb11 * (c13)); Pc[3, 9] = c09 * (pb00 * (c03) + pb10 * (c13)) + c19 * (pb01 * (c03) + pb11 * (c13));
                            Pc[4, 0] = c00 * (pb00 * (c04) + pb10 * (c14)) + c10 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 1] = c01 * (pb00 * (c04) + pb10 * (c14)) + c11 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 2] = c02 * (pb00 * (c04) + pb10 * (c14)) + c12 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 3] = c03 * (pb00 * (c04) + pb10 * (c14)) + c13 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 4] = c04 * (pb00 * (c04) + pb10 * (c14)) + c14 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 5] = c05 * (pb00 * (c04) + pb10 * (c14)) + c15 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 6] = c06 * (pb00 * (c04) + pb10 * (c14)) + c16 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 7] = c07 * (pb00 * (c04) + pb10 * (c14)) + c17 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 8] = c08 * (pb00 * (c04) + pb10 * (c14)) + c18 * (pb01 * (c04) + pb11 * (c14)); Pc[4, 9] = c09 * (pb00 * (c04) + pb10 * (c14)) + c19 * (pb01 * (c04) + pb11 * (c14));
                            Pc[5, 0] = c00 * (pb00 * (c05) + pb10 * (c15)) + c10 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 1] = c01 * (pb00 * (c05) + pb10 * (c15)) + c11 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 2] = c02 * (pb00 * (c05) + pb10 * (c15)) + c12 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 3] = c03 * (pb00 * (c05) + pb10 * (c15)) + c13 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 4] = c04 * (pb00 * (c05) + pb10 * (c15)) + c14 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 5] = c05 * (pb00 * (c05) + pb10 * (c15)) + c15 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 6] = c06 * (pb00 * (c05) + pb10 * (c15)) + c16 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 7] = c07 * (pb00 * (c05) + pb10 * (c15)) + c17 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 8] = c08 * (pb00 * (c05) + pb10 * (c15)) + c18 * (pb01 * (c05) + pb11 * (c15)); Pc[5, 9] = c09 * (pb00 * (c05) + pb10 * (c15)) + c19 * (pb01 * (c05) + pb11 * (c15));
                            Pc[6, 0] = c00 * (pb00 * (c06) + pb10 * (c16)) + c10 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 1] = c01 * (pb00 * (c06) + pb10 * (c16)) + c11 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 2] = c02 * (pb00 * (c06) + pb10 * (c16)) + c12 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 3] = c03 * (pb00 * (c06) + pb10 * (c16)) + c13 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 4] = c04 * (pb00 * (c06) + pb10 * (c16)) + c14 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 5] = c05 * (pb00 * (c06) + pb10 * (c16)) + c15 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 6] = c06 * (pb00 * (c06) + pb10 * (c16)) + c16 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 7] = c07 * (pb00 * (c06) + pb10 * (c16)) + c17 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 8] = c08 * (pb00 * (c06) + pb10 * (c16)) + c18 * (pb01 * (c06) + pb11 * (c16)); Pc[6, 9] = c09 * (pb00 * (c06) + pb10 * (c16)) + c19 * (pb01 * (c06) + pb11 * (c16));
                            Pc[7, 0] = c00 * (pb00 * (c07) + pb10 * (c17)) + c10 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 1] = c01 * (pb00 * (c07) + pb10 * (c17)) + c11 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 2] = c02 * (pb00 * (c07) + pb10 * (c17)) + c12 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 3] = c03 * (pb00 * (c07) + pb10 * (c17)) + c13 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 4] = c04 * (pb00 * (c07) + pb10 * (c17)) + c14 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 5] = c05 * (pb00 * (c07) + pb10 * (c17)) + c15 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 6] = c06 * (pb00 * (c07) + pb10 * (c17)) + c16 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 7] = c07 * (pb00 * (c07) + pb10 * (c17)) + c17 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 8] = c08 * (pb00 * (c07) + pb10 * (c17)) + c18 * (pb01 * (c07) + pb11 * (c17)); Pc[7, 9] = c09 * (pb00 * (c07) + pb10 * (c17)) + c19 * (pb01 * (c07) + pb11 * (c17));
                            Pc[8, 0] = c00 * (pb00 * (c08) + pb10 * (c18)) + c10 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 1] = c01 * (pb00 * (c08) + pb10 * (c18)) + c11 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 2] = c02 * (pb00 * (c08) + pb10 * (c18)) + c12 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 3] = c03 * (pb00 * (c08) + pb10 * (c18)) + c13 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 4] = c04 * (pb00 * (c08) + pb10 * (c18)) + c14 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 5] = c05 * (pb00 * (c08) + pb10 * (c18)) + c15 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 6] = c06 * (pb00 * (c08) + pb10 * (c18)) + c16 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 7] = c07 * (pb00 * (c08) + pb10 * (c18)) + c17 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 8] = c08 * (pb00 * (c08) + pb10 * (c18)) + c18 * (pb01 * (c08) + pb11 * (c18)); Pc[8, 9] = c09 * (pb00 * (c08) + pb10 * (c18)) + c19 * (pb01 * (c08) + pb11 * (c18));
                            Pc[9, 0] = c00 * (pb00 * (c09) + pb10 * (c19)) + c10 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 1] = c01 * (pb00 * (c09) + pb10 * (c19)) + c11 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 2] = c02 * (pb00 * (c09) + pb10 * (c19)) + c12 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 3] = c03 * (pb00 * (c09) + pb10 * (c19)) + c13 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 4] = c04 * (pb00 * (c09) + pb10 * (c19)) + c14 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 5] = c05 * (pb00 * (c09) + pb10 * (c19)) + c15 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 6] = c06 * (pb00 * (c09) + pb10 * (c19)) + c16 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 7] = c07 * (pb00 * (c09) + pb10 * (c19)) + c17 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 8] = c08 * (pb00 * (c09) + pb10 * (c19)) + c18 * (pb01 * (c09) + pb11 * (c19)); Pc[9, 9] = c09 * (pb00 * (c09) + pb10 * (c19)) + c19 * (pb01 * (c09) + pb11 * (c19));
                      //  double[,] invA = Matrix.Multiply((Matrix.Multiply(Matrix.Transpose(c), PB)), c);

                        double[,] eajj = new double[6, 1];
                            eajj[0, 0] = ll0 * (pb00 * (a00) + pb10 * (a10)) + ll1 * (pb01 * (a00) + pb11 * (a10));
                            eajj[1, 0] = ll0 * (pb00 * (a01) + pb10 * (a11)) + ll1 * (pb01 * (a01) + pb11 * (a11));
                            eajj[2, 0] = ll0 * (pb00 * (a02) + pb10 * (a12)) + ll1 * (pb01 * (a02) + pb11 * (a12));
                            eajj[3, 0] = ll0 * (pb00 * (a03) + pb10 * (a13)) + ll1 * (pb01 * (a03) + pb11 * (a13));
                            eajj[4, 0] = ll0 * (pb00 * (a04) + pb10 * (a14)) + ll1 * (pb01 * (a04) + pb11 * (a14));
                            eajj[5, 0] = ll0 * (pb00 * (a05) + pb10 * (a15)) + ll1 * (pb01 * (a05) + pb11 * (a15));
                         
                            double[,] eajjc = new double[10, 1];
                            eajjc[0, 0] = ll0 * (pb00 * (c00) + pb10 * (c10)) + ll1 * (pb01 * (c00) + pb11 * (c10));
                            eajjc[1, 0] = ll0 * (pb00 * (c01) + pb10 * (c11)) + ll1 * (pb01 * (c01) + pb11 * (c11));
                            eajjc[2, 0] = ll0 * (pb00 * (c02) + pb10 * (c12)) + ll1 * (pb01 * (c02) + pb11 * (c12));
                            eajjc[3, 0] = ll0 * (pb00 * (c03) + pb10 * (c13)) + ll1 * (pb01 * (c03) + pb11 * (c13));
                            eajjc[4, 0] = ll0 * (pb00 * (c04) + pb10 * (c14)) + ll1 * (pb01 * (c04) + pb11 * (c14));
                            eajjc[5, 0] = ll0 * (pb00 * (c05) + pb10 * (c15)) + ll1 * (pb01 * (c05) + pb11 * (c15));
                            eajjc[6, 0] = ll0 * (pb00 * (c06) + pb10 * (c16)) + ll1 * (pb01 * (c06) + pb11 * (c16));
                            eajjc[7, 0] = ll0 * (pb00 * (c07) + pb10 * (c17)) + ll1 * (pb01 * (c07) + pb11 * (c17));
                            eajjc[8, 0] = ll0 * (pb00 * (c08) + pb10 * (c18)) + ll1 * (pb01 * (c08) + pb11 * (c18));
                            eajjc[9, 0] = ll0 * (pb00 * (c09) + pb10 * (c19)) + ll1 * (pb01 * (c09) + pb11 * (c19));
                       
                            Fc[0, 0] = Fc[0, 0] + Fc2[0, 0]; Fc[0, 1] = Fc[0, 1] + Fc2[0, 1]; Fc[0, 2] = Fc[0, 2] + Fc2[0, 2];
                            Fc[1, 0] = Fc[1, 0] + Fc2[1, 0]; Fc[1, 1] = Fc[1, 1] + Fc2[1, 1]; Fc[1, 2] = Fc[1, 2] + Fc2[1, 2];
                            Fc[2, 0] = Fc[2, 0] + Fc2[2, 0]; Fc[2, 1] = Fc[2, 1] + Fc2[2, 1]; Fc[2, 2] = Fc[2, 2] + Fc2[2, 2];
                            Fc[3, 0] = Fc[3, 0] + Fc2[3, 0]; Fc[3, 1] = Fc[3, 1] + Fc2[3, 1]; Fc[3, 2] = Fc[3, 2] + Fc2[3, 2];
                            Fc[4, 0] = Fc[4, 0] + Fc2[4, 0]; Fc[4, 1] = Fc[4, 1] + Fc2[4, 1]; Fc[4, 2] = Fc[4, 2] + Fc2[4, 2];
                            Fc[5, 0] = Fc[5, 0] + Fc2[5, 0]; Fc[5, 1] = Fc[5, 1] + Fc2[5, 1]; Fc[5, 2] = Fc[5, 2] + Fc2[5, 2];
                            Fc[6, 0] = Fc[6, 0] + Fc2[6, 0]; Fc[6, 1] = Fc[6, 1] + Fc2[6, 1]; Fc[6, 2] = Fc[6, 2] + Fc2[6, 2];
                            Fc[7, 0] = Fc[7, 0] + Fc2[7, 0]; Fc[7, 1] = Fc[7, 1] + Fc2[7, 1]; Fc[7, 2] = Fc[7, 2] + Fc2[7, 2];
                            Fc[8, 0] = Fc[8, 0] + Fc2[8, 0]; Fc[8, 1] = Fc[8, 1] + Fc2[8, 1]; Fc[8, 2] = Fc[8, 2] + Fc2[8, 2];
                            Fc[9, 0] = Fc[9, 0] + Fc2[9, 0]; Fc[9, 1] = Fc[9, 1] + Fc2[9, 1]; Fc[9, 2] = Fc[9, 2] + Fc2[9, 2];

                          
                            for (int pp = 0; pp < 6; pp++)
                            {
                                eawveb[6 * j2 + pp, 0] = eawveb[6 * j2 + pp, 0] + eajj[pp, 0];
                                for (int ppp = 0; ppp < 10; ppp++)
                                {

                                    S[6 * j2 + pp, 6 * image_name_Eterior.Length + ppp] = S[6 * j2 + pp, 6 * image_name_Eterior.Length + ppp] + Ec[pp, ppp];
                                    S[6 * image_name_Eterior.Length + ppp, 6 * j2 + pp] = S[6 * image_name_Eterior.Length + ppp, 6 * j2 + pp] + Ect[ppp, pp];
                                }
                            }

                            for (int pp = 0; pp < 10; pp++)
                            {
                                eawveb[6 * image_name_Eterior.Length + pp, 0] = eawveb[6 * image_name_Eterior.Length + pp, 0] + eajjc[pp, 0];
                                for (int ppp = 0; ppp < 10; ppp++)
                                {
                                    S[6 * image_name_Eterior.Length + pp, 6 * image_name_Eterior.Length + ppp] = S[6 * image_name_Eterior.Length + pp, 6 * image_name_Eterior.Length + ppp] + Pc[pp, ppp];
                                }
                            }



                            //     MessageBox.Show("1sina" + eawveb[6 * 0 + 3, 0]);

                            //double[,] ebii2 = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(BB), PB), lll);
                            double[,] ebii2 = new double[3, 1];
                            ebii2[0, 0] = ll0 * (pb00 * (bb00) + pb10 * (bb10)) + ll1 * (pb01 * (bb00) + pb11 * (bb10));
                            ebii2[1, 0] = ll0 * (pb00 * (bb01) + pb10 * (bb11)) + ll1 * (pb01 * (bb01) + pb11 * (bb11));
                            ebii2[2, 0] = ll0 * (pb00 * (bb02) + pb10 * (bb12)) + ll1 * (pb01 * (bb02) + pb11 * (bb12));
                            ebii[0, 0] = ebii[0, 0] + ebii2[0, 0]; ebii[1, 0] = ebii[1, 0] + ebii2[1, 0];
                            ebii[2, 0] = ebii[2, 0] + ebii2[2, 0];

                            //  MessageBox.Show("" + ebii[0, 0] + "   " + ebii2[0, 0]);
                            //  v2 = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(BB), PB), BB);
                            v2[0, 0] = bb00 * (pb00 * (bb00) + pb10 * (bb10)) + bb10 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 1] = bb01 * (pb00 * (bb00) + pb10 * (bb10)) + bb11 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 2] = bb02 * (pb00 * (bb00) + pb10 * (bb10)) + bb12 * (pb01 * (bb00) + pb11 * (bb10));
                            v2[1, 0] = bb00 * (pb00 * (bb01) + pb10 * (bb11)) + bb10 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 1] = bb01 * (pb00 * (bb01) + pb10 * (bb11)) + bb11 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 2] = bb02 * (pb00 * (bb01) + pb10 * (bb11)) + bb12 * (pb01 * (bb01) + pb11 * (bb11));
                            v2[2, 0] = bb00 * (pb00 * (bb02) + pb10 * (bb12)) + bb10 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 1] = bb01 * (pb00 * (bb02) + pb10 * (bb12)) + bb11 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 2] = bb02 * (pb00 * (bb02) + pb10 * (bb12)) + bb12 * (pb01 * (bb02) + pb11 * (bb12));
                            // ww = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(AA), PB), BB);
                            ww[0, 0] = bb00 * (pb00 * (a00) + pb10 * (a10)) + bb10 * (pb01 * (a00) + pb11 * (a10)); ww[0, 1] = bb01 * (pb00 * (a00) + pb10 * (a10)) + bb11 * (pb01 * (a00) + pb11 * (a10)); ww[0, 2] = bb02 * (pb00 * (a00) + pb10 * (a10)) + bb12 * (pb01 * (a00) + pb11 * (a10));
                            ww[1, 0] = bb00 * (pb00 * (a01) + pb10 * (a11)) + bb10 * (pb01 * (a01) + pb11 * (a11)); ww[1, 1] = bb01 * (pb00 * (a01) + pb10 * (a11)) + bb11 * (pb01 * (a01) + pb11 * (a11)); ww[1, 2] = bb02 * (pb00 * (a01) + pb10 * (a11)) + bb12 * (pb01 * (a01) + pb11 * (a11));
                            ww[2, 0] = bb00 * (pb00 * (a02) + pb10 * (a12)) + bb10 * (pb01 * (a02) + pb11 * (a12)); ww[2, 1] = bb01 * (pb00 * (a02) + pb10 * (a12)) + bb11 * (pb01 * (a02) + pb11 * (a12)); ww[2, 2] = bb02 * (pb00 * (a02) + pb10 * (a12)) + bb12 * (pb01 * (a02) + pb11 * (a12));
                            ww[3, 0] = bb00 * (pb00 * (a03) + pb10 * (a13)) + bb10 * (pb01 * (a03) + pb11 * (a13)); ww[3, 1] = bb01 * (pb00 * (a03) + pb10 * (a13)) + bb11 * (pb01 * (a03) + pb11 * (a13)); ww[3, 2] = bb02 * (pb00 * (a03) + pb10 * (a13)) + bb12 * (pb01 * (a03) + pb11 * (a13));
                            ww[4, 0] = bb00 * (pb00 * (a04) + pb10 * (a14)) + bb10 * (pb01 * (a04) + pb11 * (a14)); ww[4, 1] = bb01 * (pb00 * (a04) + pb10 * (a14)) + bb11 * (pb01 * (a04) + pb11 * (a14)); ww[4, 2] = bb02 * (pb00 * (a04) + pb10 * (a14)) + bb12 * (pb01 * (a04) + pb11 * (a14));
                            ww[5, 0] = bb00 * (pb00 * (a05) + pb10 * (a15)) + bb10 * (pb01 * (a05) + pb11 * (a15)); ww[5, 1] = bb01 * (pb00 * (a05) + pb10 * (a15)) + bb11 * (pb01 * (a05) + pb11 * (a15)); ww[5, 2] = bb02 * (pb00 * (a05) + pb10 * (a15)) + bb12 * (pb01 * (a05) + pb11 * (a15));
                            w1.Add(ww[0, 0]);
                            w1.Add(ww[1, 0]);
                            w1.Add(ww[2, 0]);
                            w1.Add(ww[3, 0]);
                            w1.Add(ww[4, 0]);
                            w1.Add(ww[5, 0]);

                            w2.Add(ww[0, 1]);
                            w2.Add(ww[1, 1]);
                            w2.Add(ww[2, 1]);
                            w2.Add(ww[3, 1]);
                            w2.Add(ww[4, 1]);
                            w2.Add(ww[5, 1]);

                            w3.Add(ww[0, 2]);
                            w3.Add(ww[1, 2]);
                            w3.Add(ww[2, 2]);
                            w3.Add(ww[3, 2]);
                            w3.Add(ww[4, 2]);
                            w3.Add(ww[5, 2]);

                            wj.Add(j2);
                            //WW[6 * j2 + 0, 0] = ww[0, 0]; WW[6 * j2 + 1, 0] = ww[1, 0]; WW[6 * j2 + 2, 0] = ww[2, 0]; WW[6 * j2 + 3, 0] = ww[3, 0]; WW[6 * j2 + 4, 0] = ww[4, 0]; WW[6 * j2 + 5, 0] = ww[5, 0];
                            //WW[6 * j2 + 0, 1] = ww[0, 1]; WW[6 * j2 + 1, 1] = ww[1,1]; WW[6 * j2 + 2, 1] = ww[2, 1]; WW[6 * j2 + 3, 1] = ww[3, 1]; WW[6 * j2 + 4, 1] = ww[4, 1]; WW[6 * j2 + 5, 1] = ww[5, 1];
                            //WW[6 * j2 + 0, 2] = ww[0, 2]; WW[6 * j2 + 1, 2] = ww[1, 2]; WW[6 * j2 + 2, 2] = ww[2, 2]; WW[6 * j2 + 3, 2] = ww[3,2]; WW[6 * j2 + 4,2] = ww[4, 2]; WW[6 * j2 + 5, 2] = ww[5, 2];
                            v[0, 0] = v[0, 0] + v2[0, 0]; v[1, 0] = v[1, 0] + v2[1, 0]; v[2, 0] = v[2, 0] + v2[2, 0];
                            v[0, 1] = v[0, 1] + v2[0, 1]; v[1, 1] = v[1, 1] + v2[1, 1]; v[2, 1] = v[2, 1] + v2[2, 1];
                            v[0, 2] = v[0, 2] + v2[0, 2]; v[1, 2] = v[1, 2] + v2[1, 2]; v[2, 2] = v[2, 2] + v2[2, 2];




                            //  MessageBox.Show("" + F[0, 0]+"   "+ FF[0, 0]+"    "+ BB[0,0]+"  "+j2+"  "+i);
                            //   u = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(AA), PB), AA);
                            u[0, 0] = a00 * (pb00 * (a00) + pb10 * (a10)) + a10 * (pb01 * (a00) + pb11 * (a10)); u[0, 1] = a01 * (pb00 * (a00) + pb10 * (a10)) + a11 * (pb01 * (a00) + pb11 * (a10)); u[0, 2] = a02 * (pb00 * (a00) + pb10 * (a10)) + a12 * (pb01 * (a00) + pb11 * (a10)); u[0, 3] = a03 * (pb00 * (a00) + pb10 * (a10)) + a13 * (pb01 * (a00) + pb11 * (a10)); u[0, 4] = a04 * (pb00 * (a00) + pb10 * (a10)) + a14 * (pb01 * (a00) + pb11 * (a10)); u[0, 5] = a05 * (pb00 * (a00) + pb10 * (a10)) + a15 * (pb01 * (a00) + pb11 * (a10));
                            u[1, 0] = a00 * (pb00 * (a01) + pb10 * (a11)) + a10 * (pb01 * (a01) + pb11 * (a11)); u[1, 1] = a01 * (pb00 * (a01) + pb10 * (a11)) + a11 * (pb01 * (a01) + pb11 * (a11)); u[1, 2] = a02 * (pb00 * (a01) + pb10 * (a11)) + a12 * (pb01 * (a01) + pb11 * (a11)); u[1, 3] = a03 * (pb00 * (a01) + pb10 * (a11)) + a13 * (pb01 * (a01) + pb11 * (a11)); u[1, 4] = a04 * (pb00 * (a01) + pb10 * (a11)) + a14 * (pb01 * (a01) + pb11 * (a11)); u[1, 5] = a05 * (pb00 * (a01) + pb10 * (a11)) + a15 * (pb01 * (a01) + pb11 * (a11));
                            u[2, 0] = a00 * (pb00 * (a02) + pb10 * (a12)) + a10 * (pb01 * (a02) + pb11 * (a12)); u[2, 1] = a01 * (pb00 * (a02) + pb10 * (a12)) + a11 * (pb01 * (a02) + pb11 * (a12)); u[2, 2] = a02 * (pb00 * (a02) + pb10 * (a12)) + a12 * (pb01 * (a02) + pb11 * (a12)); u[2, 3] = a03 * (pb00 * (a02) + pb10 * (a12)) + a13 * (pb01 * (a02) + pb11 * (a12)); u[2, 4] = a04 * (pb00 * (a02) + pb10 * (a12)) + a14 * (pb01 * (a02) + pb11 * (a12)); u[2, 5] = a05 * (pb00 * (a02) + pb10 * (a12)) + a15 * (pb01 * (a02) + pb11 * (a12));
                            u[3, 0] = a00 * (pb00 * (a03) + pb10 * (a13)) + a10 * (pb01 * (a03) + pb11 * (a13)); u[3, 1] = a01 * (pb00 * (a03) + pb10 * (a13)) + a11 * (pb01 * (a03) + pb11 * (a13)); u[3, 2] = a02 * (pb00 * (a03) + pb10 * (a13)) + a12 * (pb01 * (a03) + pb11 * (a13)); u[3, 3] = a03 * (pb00 * (a03) + pb10 * (a13)) + a13 * (pb01 * (a03) + pb11 * (a13)); u[3, 4] = a04 * (pb00 * (a03) + pb10 * (a13)) + a14 * (pb01 * (a03) + pb11 * (a13)); u[3, 5] = a05 * (pb00 * (a03) + pb10 * (a13)) + a15 * (pb01 * (a03) + pb11 * (a13));
                            u[4, 0] = a00 * (pb00 * (a04) + pb10 * (a14)) + a10 * (pb01 * (a04) + pb11 * (a14)); u[4, 1] = a01 * (pb00 * (a04) + pb10 * (a14)) + a11 * (pb01 * (a04) + pb11 * (a14)); u[4, 2] = a02 * (pb00 * (a04) + pb10 * (a14)) + a12 * (pb01 * (a04) + pb11 * (a14)); u[4, 3] = a03 * (pb00 * (a04) + pb10 * (a14)) + a13 * (pb01 * (a04) + pb11 * (a14)); u[4, 4] = a04 * (pb00 * (a04) + pb10 * (a14)) + a14 * (pb01 * (a04) + pb11 * (a14)); u[4, 5] = a05 * (pb00 * (a04) + pb10 * (a14)) + a15 * (pb01 * (a04) + pb11 * (a14));
                            u[5, 0] = a00 * (pb00 * (a05) + pb10 * (a15)) + a10 * (pb01 * (a05) + pb11 * (a15)); u[5, 1] = a01 * (pb00 * (a05) + pb10 * (a15)) + a11 * (pb01 * (a05) + pb11 * (a15)); u[5, 2] = a02 * (pb00 * (a05) + pb10 * (a15)) + a12 * (pb01 * (a05) + pb11 * (a15)); u[5, 3] = a03 * (pb00 * (a05) + pb10 * (a15)) + a13 * (pb01 * (a05) + pb11 * (a15)); u[5, 4] = a04 * (pb00 * (a05) + pb10 * (a15)) + a14 * (pb01 * (a05) + pb11 * (a15)); u[5, 5] = a05 * (pb00 * (a05) + pb10 * (a15)) + a15 * (pb01 * (a05) + pb11 * (a15));
                            for (int pp = 0; pp < 6; pp++)
                            {
                                for (int ppp = 0; ppp < 6; ppp++)
                                {

                                    S[6 * j2 + pp, 6 * j2 + ppp] = S[6 * j2 + pp, 6 * j2 + ppp] + u[pp, ppp];

                                }
                            }
                            

                       
                    }

                    //   MessageBox.Show("" + eaj[0, 0] + "   " + eaj[6 * image_name_Eterior.Length + 1, 0]);
                    v[0, 0] = v[0, 0] + pp2nB[0, 0]; v[1, 1] = v[1, 1] + pp2nB[1, 1]; v[2, 2] = v[2, 2] + pp2nB[2, 2];
                    ebii[0, 0] = ebii[0, 0] - (pp2nB[0, 0] * lllnb[0, 0]);
                    ebii[1, 0] = ebii[1, 0] - (pp2nB[1, 1] * lllnb[1, 0]);
                    ebii[2, 0] = ebii[2, 0] - (pp2nB[2, 2] * lllnb[2, 0]);
                    double ebii00 = ebii[0, 0] - (pp2nB[0, 0] * lllnb[0, 0]);
                    double ebii10 = ebii[1, 0] - (pp2nB[1, 1] * lllnb[1, 0]);
                    double ebii20 = ebii[2, 0] - (pp2nB[2, 2] * lllnb[2, 0]);
                    double[,] v3 = new double[3, 3];
                    double detv = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                    if (detv == 0)
                    {

                        MessageBox.Show("detv");
                    }
                    v3[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                    v3[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                    v3[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                    v3[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                    v3[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                    v3[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                    v3[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                    v3[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                    v3[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                    double v300 = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                    double v301 = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                    double v302 = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                    double v310 = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                    double v311 = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                    double v312 = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                    double v320 = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                    double v321 = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                    double v322 = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;

                    Fc00 = Fc[0, 0]; Fc01 = Fc[0, 1]; Fc02 = Fc[0, 2]; Fc10 = Fc[1, 0]; Fc11 = Fc[1, 1]; Fc12 = Fc[1, 2];
                    Fc20 = Fc[2, 0]; Fc21 = Fc[2, 1]; Fc22 = Fc[2, 2]; Fc30 = Fc[3, 0]; Fc31 = Fc[3, 1]; Fc32 = Fc[3, 2];
                    Fc40 = Fc[4, 0]; Fc41 = Fc[4, 1]; Fc42 = Fc[4, 2]; Fc50 = Fc[5, 0]; Fc51 = Fc[5, 1]; Fc52 = Fc[5, 2];
                    Fc60 = Fc[6, 0]; Fc61 = Fc[6, 1]; Fc62 = Fc[6, 2]; Fc70 = Fc[7, 0]; Fc71 = Fc[7, 1]; Fc72 = Fc[7, 2];
                    Fc80 = Fc[8, 0]; Fc81 = Fc[8, 1]; Fc82 = Fc[8, 2]; Fc90 = Fc[9, 0]; Fc91 = Fc[9, 1]; Fc92 = Fc[9, 2];

                    // double[,] Fveb = Matrix.Multiply((Matrix.Multiply(Fc, v3)), ebii);
                    double[,] Fveb = new double[10, 1];
                    Fveb[0, 0] = ebii00 * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + ebii10 * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + ebii20 * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322);
                    Fveb[1, 0] = ebii00 * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + ebii10 * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + ebii20 * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322);
                    Fveb[2, 0] = ebii00 * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + ebii10 * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + ebii20 * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322);
                    Fveb[3, 0] = ebii00 * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + ebii10 * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + ebii20 * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322);
                    Fveb[4, 0] = ebii00 * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + ebii10 * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + ebii20 * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322);
                    Fveb[5, 0] = ebii00 * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + ebii10 * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + ebii20 * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322);
                    Fveb[6, 0] = ebii00 * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + ebii10 * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + ebii20 * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322);
                    Fveb[7, 0] = ebii00 * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + ebii10 * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + ebii20 * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322);
                    Fveb[8, 0] = ebii00 * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + ebii10 * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + ebii20 * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322);
                    Fveb[9, 0] = ebii00 * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + ebii10 * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + ebii20 * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322);
                    //var s1 = Stopwatch.StartNew();
                    for (int iop = 0; iop < 10; iop++)
                    {
                        eawveb[6 * image_name_Eterior.Length + iop, 0] = eawveb[6 * image_name_Eterior.Length + iop, 0] - Fveb[iop, 0];
                    }




                    for (int j2 = 0; j2 < wj.Count; j2++)
                    {
                        double[,] wcalj2 = new double[6, 3];
                        double wcalj200, wcalj201, wcalj202, wcalj210, wcalj211, wcalj212, wcalj220, wcalj221, wcalj222, wcalj230, wcalj231, wcalj232, wcalj240, wcalj241, wcalj242, wcalj250, wcalj251, wcalj252;
                        for (int iop = 0; iop < 6; iop++)
                        {

                            wcalj2[iop, 0] = w1[6 * j2 + iop];
                            wcalj2[iop, 1] = w2[6 * j2 + iop];
                            wcalj2[iop, 2] = w3[6 * j2 + iop];

                        }
                        wcalj200 = wcalj2[0, 0]; wcalj201 = wcalj2[0, 1]; wcalj202 = wcalj2[0, 2]; wcalj210 = wcalj2[1, 0]; wcalj211 = wcalj2[1, 1]; wcalj212 = wcalj2[1, 2]; wcalj220 = wcalj2[2, 0]; wcalj221 = wcalj2[2, 1]; wcalj222 = wcalj2[2, 2]; wcalj230 = wcalj2[3, 0]; wcalj231 = wcalj2[3, 1]; wcalj232 = wcalj2[3, 2]; wcalj240 = wcalj2[4, 0]; wcalj241 = wcalj2[4, 1]; wcalj242 = wcalj2[4, 2]; wcalj250 = wcalj2[5, 0]; wcalj251 = wcalj2[5, 1]; wcalj252 = wcalj2[5, 2];
                        //  double[,] eawveb2 = Matrix.Multiply((Matrix.Multiply(wcalj2, v3)), ebii);
                        double[,] eawveb2 = new double[6, 1];
                        eawveb2[0, 0] = ebii00 * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + ebii10 * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + ebii20 * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202);
                        eawveb2[1, 0] = ebii00 * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + ebii10 * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + ebii20 * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212);
                        eawveb2[2, 0] = ebii00 * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + ebii10 * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + ebii20 * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222);
                        eawveb2[3, 0] = ebii00 * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + ebii10 * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + ebii20 * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232);
                        eawveb2[4, 0] = ebii00 * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + ebii10 * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + ebii20 * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242);
                        eawveb2[5, 0] = ebii00 * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + ebii10 * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + ebii20 * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252);


                        for (int iop = 0; iop < 6; iop++)
                        {
                            eawveb[6 * wj[j2] + iop, 0] = eawveb[6 * wj[j2] + iop, 0] - eawveb2[iop, 0];
                        }
                        //   MessageBox.Show("" + eawveb2[0, 0]);
                        for (int j = 0; j < wj.Count; j++)
                        {

                            double[,] wcalj = new double[6, 3];
                            for (int iop = 0; iop < 6; iop++)
                            {
                                wcalj[iop, 0] = w1[6 * j + iop];
                                wcalj[iop, 1] = w2[6 * j + iop];
                                wcalj[iop, 2] = w3[6 * j + iop];
                            }
                            double wcalj00 = wcalj[0, 0]; double wcalj01 = wcalj[0, 1]; double wcalj02 = wcalj[0, 2]; double wcalj10 = wcalj[1, 0]; double wcalj11 = wcalj[1, 1]; double wcalj12 = wcalj[1, 2]; double wcalj20 = wcalj[2, 0]; double wcalj21 = wcalj[2, 1]; double wcalj22 = wcalj[2, 2]; double wcalj30 = wcalj[3, 0]; double wcalj31 = wcalj[3, 1]; double wcalj32 = wcalj[3, 2]; double wcalj40 = wcalj[4, 0]; double wcalj41 = wcalj[4, 1]; double wcalj42 = wcalj[4, 2]; double wcalj50 = wcalj[5, 0]; double wcalj51 = wcalj[5, 1]; double wcalj52 = wcalj[5, 2];
                            //   double[,] wcalj2wcalj = Matrix.Multiply((Matrix.Multiply(wcalj2, v3)), Matrix.Transpose(wcalj));
                            double[,] wcalj2wcalj = new double[6, 6];
                            wcalj2wcalj[0, 0] = (wcalj00) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj01) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj02) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 1] = (wcalj10) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj11) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj12) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 2] = (wcalj20) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj21) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj22) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 3] = (wcalj30) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj31) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj32) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 4] = (wcalj40) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj41) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj42) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2wcalj[0, 5] = (wcalj50) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (wcalj51) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (wcalj52) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202);
                            wcalj2wcalj[1, 0] = (wcalj00) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj01) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj02) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 1] = (wcalj10) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj11) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj12) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 2] = (wcalj20) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj21) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj22) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 3] = (wcalj30) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj31) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj32) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 4] = (wcalj40) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj41) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj42) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2wcalj[1, 5] = (wcalj50) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (wcalj51) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (wcalj52) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212);
                            wcalj2wcalj[2, 0] = (wcalj00) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj01) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj02) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 1] = (wcalj10) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj11) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj12) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 2] = (wcalj20) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj21) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj22) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 3] = (wcalj30) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj31) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj32) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 4] = (wcalj40) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj41) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj42) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2wcalj[2, 5] = (wcalj50) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (wcalj51) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (wcalj52) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222);
                            wcalj2wcalj[3, 0] = (wcalj00) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj01) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj02) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 1] = (wcalj10) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj11) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj12) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 2] = (wcalj20) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj21) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj22) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 3] = (wcalj30) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj31) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj32) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 4] = (wcalj40) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj41) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj42) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2wcalj[3, 5] = (wcalj50) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (wcalj51) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (wcalj52) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232);
                            wcalj2wcalj[4, 0] = (wcalj00) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj01) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj02) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 1] = (wcalj10) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj11) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj12) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 2] = (wcalj20) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj21) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj22) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 3] = (wcalj30) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj31) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj32) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 4] = (wcalj40) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj41) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj42) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2wcalj[4, 5] = (wcalj50) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (wcalj51) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (wcalj52) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242);
                            wcalj2wcalj[5, 0] = (wcalj00) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj01) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj02) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 1] = (wcalj10) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj11) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj12) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 2] = (wcalj20) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj21) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj22) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 3] = (wcalj30) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj31) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj32) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 4] = (wcalj40) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj41) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj42) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2wcalj[5, 5] = (wcalj50) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (wcalj51) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (wcalj52) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252);
                            for (int pp = 0; pp < 6; pp++)
                            {
                                for (int ppp = 0; ppp < 6; ppp++)
                                {
                                    S[6 * wj[j2] + pp, 6 * wj[j] + ppp] = S[6 * wj[j2] + pp, 6 * wj[j] + ppp] - wcalj2wcalj[pp, ppp];
                                }
                            }
                        }


                        // double[,] Fcvwcalj = Matrix.Multiply((Matrix.Multiply(Fc, v3)), Matrix.Transpose(wcalj2));
                        double[,] Fcvwcalj = new double[10, 6];
                        Fcvwcalj[0, 0] = (wcalj200) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (wcalj201) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (wcalj202) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); Fcvwcalj[0, 1] = (wcalj210) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (wcalj211) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (wcalj212) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); Fcvwcalj[0, 2] = (wcalj220) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (wcalj221) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (wcalj222) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); Fcvwcalj[0, 3] = (wcalj230) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (wcalj231) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (wcalj232) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); Fcvwcalj[0, 4] = (wcalj240) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (wcalj241) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (wcalj242) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); Fcvwcalj[0, 5] = (wcalj250) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (wcalj251) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (wcalj252) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322);
                        Fcvwcalj[1, 0] = (wcalj200) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (wcalj201) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (wcalj202) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); Fcvwcalj[1, 1] = (wcalj210) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (wcalj211) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (wcalj212) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); Fcvwcalj[1, 2] = (wcalj220) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (wcalj221) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (wcalj222) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); Fcvwcalj[1, 3] = (wcalj230) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (wcalj231) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (wcalj232) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); Fcvwcalj[1, 4] = (wcalj240) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (wcalj241) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (wcalj242) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); Fcvwcalj[1, 5] = (wcalj250) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (wcalj251) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (wcalj252) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322);
                        Fcvwcalj[2, 0] = (wcalj200) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (wcalj201) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (wcalj202) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); Fcvwcalj[2, 1] = (wcalj210) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (wcalj211) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (wcalj212) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); Fcvwcalj[2, 2] = (wcalj220) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (wcalj221) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (wcalj222) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); Fcvwcalj[2, 3] = (wcalj230) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (wcalj231) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (wcalj232) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); Fcvwcalj[2, 4] = (wcalj240) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (wcalj241) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (wcalj242) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); Fcvwcalj[2, 5] = (wcalj250) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (wcalj251) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (wcalj252) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322);
                        Fcvwcalj[3, 0] = (wcalj200) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (wcalj201) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (wcalj202) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); Fcvwcalj[3, 1] = (wcalj210) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (wcalj211) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (wcalj212) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); Fcvwcalj[3, 2] = (wcalj220) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (wcalj221) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (wcalj222) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); Fcvwcalj[3, 3] = (wcalj230) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (wcalj231) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (wcalj232) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); Fcvwcalj[3, 4] = (wcalj240) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (wcalj241) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (wcalj242) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); Fcvwcalj[3, 5] = (wcalj250) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (wcalj251) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (wcalj252) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322);
                        Fcvwcalj[4, 0] = (wcalj200) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (wcalj201) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (wcalj202) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); Fcvwcalj[4, 1] = (wcalj210) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (wcalj211) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (wcalj212) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); Fcvwcalj[4, 2] = (wcalj220) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (wcalj221) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (wcalj222) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); Fcvwcalj[4, 3] = (wcalj230) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (wcalj231) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (wcalj232) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); Fcvwcalj[4, 4] = (wcalj240) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (wcalj241) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (wcalj242) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); Fcvwcalj[4, 5] = (wcalj250) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (wcalj251) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (wcalj252) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322);
                        Fcvwcalj[5, 0] = (wcalj200) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (wcalj201) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (wcalj202) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); Fcvwcalj[5, 1] = (wcalj210) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (wcalj211) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (wcalj212) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); Fcvwcalj[5, 2] = (wcalj220) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (wcalj221) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (wcalj222) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); Fcvwcalj[5, 3] = (wcalj230) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (wcalj231) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (wcalj232) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); Fcvwcalj[5, 4] = (wcalj240) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (wcalj241) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (wcalj242) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); Fcvwcalj[5, 5] = (wcalj250) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (wcalj251) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (wcalj252) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322);
                        Fcvwcalj[6, 0] = (wcalj200) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (wcalj201) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (wcalj202) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); Fcvwcalj[6, 1] = (wcalj210) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (wcalj211) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (wcalj212) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); Fcvwcalj[6, 2] = (wcalj220) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (wcalj221) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (wcalj222) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); Fcvwcalj[6, 3] = (wcalj230) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (wcalj231) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (wcalj232) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); Fcvwcalj[6, 4] = (wcalj240) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (wcalj241) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (wcalj242) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); Fcvwcalj[6, 5] = (wcalj250) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (wcalj251) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (wcalj252) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322);
                        Fcvwcalj[7, 0] = (wcalj200) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (wcalj201) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (wcalj202) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); Fcvwcalj[7, 1] = (wcalj210) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (wcalj211) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (wcalj212) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); Fcvwcalj[7, 2] = (wcalj220) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (wcalj221) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (wcalj222) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); Fcvwcalj[7, 3] = (wcalj230) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (wcalj231) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (wcalj232) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); Fcvwcalj[7, 4] = (wcalj240) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (wcalj241) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (wcalj242) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); Fcvwcalj[7, 5] = (wcalj250) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (wcalj251) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (wcalj252) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322);
                        Fcvwcalj[8, 0] = (wcalj200) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (wcalj201) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (wcalj202) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); Fcvwcalj[8, 1] = (wcalj210) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (wcalj211) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (wcalj212) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); Fcvwcalj[8, 2] = (wcalj220) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (wcalj221) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (wcalj222) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); Fcvwcalj[8, 3] = (wcalj230) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (wcalj231) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (wcalj232) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); Fcvwcalj[8, 4] = (wcalj240) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (wcalj241) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (wcalj242) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); Fcvwcalj[8, 5] = (wcalj250) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (wcalj251) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (wcalj252) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322);
                        Fcvwcalj[9, 0] = (wcalj200) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (wcalj201) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (wcalj202) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); Fcvwcalj[9, 1] = (wcalj210) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (wcalj211) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (wcalj212) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); Fcvwcalj[9, 2] = (wcalj220) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (wcalj221) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (wcalj222) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); Fcvwcalj[9, 3] = (wcalj230) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (wcalj231) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (wcalj232) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); Fcvwcalj[9, 4] = (wcalj240) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (wcalj241) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (wcalj242) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); Fcvwcalj[9, 5] = (wcalj250) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (wcalj251) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (wcalj252) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322);


                        // double[,] wcalj2vFc = Matrix.Multiply((Matrix.Multiply(wcalj2, v3)), Matrix.Transpose(Fc));
                        double[,] wcalj2vFc = new double[6, 10];
                        wcalj2vFc[0, 0] = (Fc00) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc01) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc02) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 1] = (Fc10) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc11) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc12) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 2] = (Fc20) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc21) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc22) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 3] = (Fc30) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc31) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc32) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 4] = (Fc40) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc41) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc42) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 5] = (Fc50) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc51) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc52) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 6] = (Fc60) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc61) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc62) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 7] = (Fc70) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc71) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc72) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 8] = (Fc80) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc81) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc82) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202); wcalj2vFc[0, 9] = (Fc90) * (v300 * wcalj200 + v310 * wcalj201 + v320 * wcalj202) + (Fc91) * (v301 * wcalj200 + v311 * wcalj201 + v321 * wcalj202) + (Fc92) * (v302 * wcalj200 + v312 * wcalj201 + v322 * wcalj202);
                        wcalj2vFc[1, 0] = (Fc00) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc01) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc02) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 1] = (Fc10) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc11) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc12) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 2] = (Fc20) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc21) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc22) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 3] = (Fc30) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc31) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc32) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 4] = (Fc40) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc41) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc42) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 5] = (Fc50) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc51) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc52) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 6] = (Fc60) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc61) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc62) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 7] = (Fc70) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc71) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc72) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 8] = (Fc80) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc81) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc82) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212); wcalj2vFc[1, 9] = (Fc90) * (v300 * wcalj210 + v310 * wcalj211 + v320 * wcalj212) + (Fc91) * (v301 * wcalj210 + v311 * wcalj211 + v321 * wcalj212) + (Fc92) * (v302 * wcalj210 + v312 * wcalj211 + v322 * wcalj212);
                        wcalj2vFc[2, 0] = (Fc00) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc01) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc02) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 1] = (Fc10) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc11) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc12) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 2] = (Fc20) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc21) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc22) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 3] = (Fc30) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc31) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc32) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 4] = (Fc40) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc41) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc42) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 5] = (Fc50) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc51) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc52) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 6] = (Fc60) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc61) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc62) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 7] = (Fc70) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc71) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc72) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 8] = (Fc80) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc81) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc82) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222); wcalj2vFc[2, 9] = (Fc90) * (v300 * wcalj220 + v310 * wcalj221 + v320 * wcalj222) + (Fc91) * (v301 * wcalj220 + v311 * wcalj221 + v321 * wcalj222) + (Fc92) * (v302 * wcalj220 + v312 * wcalj221 + v322 * wcalj222);
                        wcalj2vFc[3, 0] = (Fc00) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc01) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc02) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 1] = (Fc10) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc11) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc12) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 2] = (Fc20) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc21) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc22) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 3] = (Fc30) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc31) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc32) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 4] = (Fc40) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc41) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc42) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 5] = (Fc50) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc51) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc52) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 6] = (Fc60) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc61) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc62) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 7] = (Fc70) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc71) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc72) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 8] = (Fc80) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc81) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc82) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232); wcalj2vFc[3, 9] = (Fc90) * (v300 * wcalj230 + v310 * wcalj231 + v320 * wcalj232) + (Fc91) * (v301 * wcalj230 + v311 * wcalj231 + v321 * wcalj232) + (Fc92) * (v302 * wcalj230 + v312 * wcalj231 + v322 * wcalj232);
                        wcalj2vFc[4, 0] = (Fc00) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc01) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc02) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 1] = (Fc10) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc11) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc12) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 2] = (Fc20) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc21) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc22) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 3] = (Fc30) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc31) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc32) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 4] = (Fc40) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc41) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc42) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 5] = (Fc50) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc51) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc52) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 6] = (Fc60) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc61) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc62) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 7] = (Fc70) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc71) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc72) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 8] = (Fc80) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc81) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc82) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242); wcalj2vFc[4, 9] = (Fc90) * (v300 * wcalj240 + v310 * wcalj241 + v320 * wcalj242) + (Fc91) * (v301 * wcalj240 + v311 * wcalj241 + v321 * wcalj242) + (Fc92) * (v302 * wcalj240 + v312 * wcalj241 + v322 * wcalj242);
                        wcalj2vFc[5, 0] = (Fc00) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc01) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc02) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 1] = (Fc10) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc11) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc12) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 2] = (Fc20) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc21) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc22) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 3] = (Fc30) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc31) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc32) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 4] = (Fc40) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc41) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc42) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 5] = (Fc50) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc51) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc52) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 6] = (Fc60) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc61) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc62) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 7] = (Fc70) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc71) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc72) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 8] = (Fc80) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc81) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc82) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252); wcalj2vFc[5, 9] = (Fc90) * (v300 * wcalj250 + v310 * wcalj251 + v320 * wcalj252) + (Fc91) * (v301 * wcalj250 + v311 * wcalj251 + v321 * wcalj252) + (Fc92) * (v302 * wcalj250 + v312 * wcalj251 + v322 * wcalj252);

                        for (int pp = 0; pp < 6; pp++)
                        {
                            for (int ppp = 0; ppp < 10; ppp++)
                            {
                                S[6 * wj[j2] + pp, 6 * image_name_Eterior.Length + ppp] = S[6 * wj[j2] + pp, 6 * image_name_Eterior.Length + ppp] - wcalj2vFc[pp, ppp];
                                S[6 * image_name_Eterior.Length + ppp, 6 * wj[j2] + pp] = S[6 * image_name_Eterior.Length + ppp, 6 * wj[j2] + pp] - Fcvwcalj[ppp, pp];
                            }
                        }

                    }
                    //  double[,] FcvFc = Matrix.Multiply((Matrix.Multiply(Fc, v3)), Matrix.Transpose(Fc));
                    double[,] FcvFc = new double[10, 10];
                    FcvFc[0, 0] = (Fc00) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc01) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc02) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 1] = (Fc10) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc11) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc12) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 2] = (Fc20) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc21) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc22) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 3] = (Fc30) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc31) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc32) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 4] = (Fc40) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc41) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc42) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 5] = (Fc50) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc51) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc52) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 6] = (Fc60) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc61) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc62) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 7] = (Fc70) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc71) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc72) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 8] = (Fc80) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc81) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc82) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322); FcvFc[0, 9] = (Fc90) * (Fc00 * v300 + Fc01 * v310 + Fc02 * v320) + (Fc91) * (Fc00 * v301 + Fc01 * v311 + Fc02 * v321) + (Fc92) * (Fc00 * v302 + Fc01 * v312 + Fc02 * v322);
                    FcvFc[1, 0] = (Fc00) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc01) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc02) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 1] = (Fc10) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc11) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc12) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 2] = (Fc20) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc21) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc22) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 3] = (Fc30) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc31) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc32) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 4] = (Fc40) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc41) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc42) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 5] = (Fc50) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc51) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc52) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 6] = (Fc60) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc61) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc62) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 7] = (Fc70) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc71) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc72) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 8] = (Fc80) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc81) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc82) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322); FcvFc[1, 9] = (Fc90) * (Fc10 * v300 + Fc11 * v310 + Fc12 * v320) + (Fc91) * (Fc10 * v301 + Fc11 * v311 + Fc12 * v321) + (Fc92) * (Fc10 * v302 + Fc11 * v312 + Fc12 * v322);
                    FcvFc[2, 0] = (Fc00) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc01) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc02) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 1] = (Fc10) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc11) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc12) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 2] = (Fc20) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc21) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc22) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 3] = (Fc30) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc31) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc32) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 4] = (Fc40) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc41) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc42) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 5] = (Fc50) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc51) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc52) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 6] = (Fc60) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc61) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc62) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 7] = (Fc70) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc71) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc72) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 8] = (Fc80) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc81) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc82) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322); FcvFc[2, 9] = (Fc90) * (Fc20 * v300 + Fc21 * v310 + Fc22 * v320) + (Fc91) * (Fc20 * v301 + Fc21 * v311 + Fc22 * v321) + (Fc92) * (Fc20 * v302 + Fc21 * v312 + Fc22 * v322);
                    FcvFc[3, 0] = (Fc00) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc01) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc02) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 1] = (Fc10) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc11) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc12) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 2] = (Fc20) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc21) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc22) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 3] = (Fc30) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc31) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc32) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 4] = (Fc40) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc41) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc42) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 5] = (Fc50) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc51) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc52) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 6] = (Fc60) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc61) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc62) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 7] = (Fc70) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc71) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc72) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 8] = (Fc80) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc81) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc82) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322); FcvFc[3, 9] = (Fc90) * (Fc30 * v300 + Fc31 * v310 + Fc32 * v320) + (Fc91) * (Fc30 * v301 + Fc31 * v311 + Fc32 * v321) + (Fc92) * (Fc30 * v302 + Fc31 * v312 + Fc32 * v322);
                    FcvFc[4, 0] = (Fc00) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc01) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc02) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 1] = (Fc10) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc11) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc12) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 2] = (Fc20) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc21) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc22) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 3] = (Fc30) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc31) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc32) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 4] = (Fc40) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc41) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc42) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 5] = (Fc50) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc51) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc52) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 6] = (Fc60) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc61) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc62) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 7] = (Fc70) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc71) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc72) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 8] = (Fc80) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc81) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc82) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322); FcvFc[4, 9] = (Fc90) * (Fc40 * v300 + Fc41 * v310 + Fc42 * v320) + (Fc91) * (Fc40 * v301 + Fc41 * v311 + Fc42 * v321) + (Fc92) * (Fc40 * v302 + Fc41 * v312 + Fc42 * v322);
                    FcvFc[5, 0] = (Fc00) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc01) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc02) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 1] = (Fc10) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc11) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc12) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 2] = (Fc20) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc21) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc22) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 3] = (Fc30) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc31) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc32) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 4] = (Fc40) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc41) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc42) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 5] = (Fc50) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc51) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc52) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 6] = (Fc60) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc61) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc62) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 7] = (Fc70) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc71) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc72) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 8] = (Fc80) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc81) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc82) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322); FcvFc[5, 9] = (Fc90) * (Fc50 * v300 + Fc51 * v310 + Fc52 * v320) + (Fc91) * (Fc50 * v301 + Fc51 * v311 + Fc52 * v321) + (Fc92) * (Fc50 * v302 + Fc51 * v312 + Fc52 * v322);
                    FcvFc[6, 0] = (Fc00) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc01) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc02) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 1] = (Fc10) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc11) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc12) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 2] = (Fc20) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc21) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc22) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 3] = (Fc30) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc31) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc32) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 4] = (Fc40) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc41) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc42) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 5] = (Fc50) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc51) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc52) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 6] = (Fc60) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc61) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc62) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 7] = (Fc70) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc71) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc72) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 8] = (Fc80) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc81) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc82) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322); FcvFc[6, 9] = (Fc90) * (Fc60 * v300 + Fc61 * v310 + Fc62 * v320) + (Fc91) * (Fc60 * v301 + Fc61 * v311 + Fc62 * v321) + (Fc92) * (Fc60 * v302 + Fc61 * v312 + Fc62 * v322);
                    FcvFc[7, 0] = (Fc00) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc01) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc02) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 1] = (Fc10) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc11) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc12) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 2] = (Fc20) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc21) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc22) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 3] = (Fc30) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc31) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc32) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 4] = (Fc40) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc41) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc42) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 5] = (Fc50) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc51) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc52) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 6] = (Fc60) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc61) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc62) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 7] = (Fc70) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc71) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc72) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 8] = (Fc80) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc81) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc82) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322); FcvFc[7, 9] = (Fc90) * (Fc70 * v300 + Fc71 * v310 + Fc72 * v320) + (Fc91) * (Fc70 * v301 + Fc71 * v311 + Fc72 * v321) + (Fc92) * (Fc70 * v302 + Fc71 * v312 + Fc72 * v322);
                    FcvFc[8, 0] = (Fc00) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc01) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc02) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 1] = (Fc10) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc11) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc12) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 2] = (Fc20) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc21) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc22) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 3] = (Fc30) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc31) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc32) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 4] = (Fc40) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc41) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc42) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 5] = (Fc50) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc51) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc52) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 6] = (Fc60) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc61) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc62) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 7] = (Fc70) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc71) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc72) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 8] = (Fc80) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc81) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc82) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322); FcvFc[8, 9] = (Fc90) * (Fc80 * v300 + Fc81 * v310 + Fc82 * v320) + (Fc91) * (Fc80 * v301 + Fc81 * v311 + Fc82 * v321) + (Fc92) * (Fc80 * v302 + Fc81 * v312 + Fc82 * v322);
                    FcvFc[9, 0] = (Fc00) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc01) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc02) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 1] = (Fc10) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc11) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc12) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 2] = (Fc20) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc21) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc22) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 3] = (Fc30) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc31) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc32) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 4] = (Fc40) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc41) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc42) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 5] = (Fc50) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc51) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc52) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 6] = (Fc60) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc61) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc62) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 7] = (Fc70) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc71) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc72) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 8] = (Fc80) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc81) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc82) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322); FcvFc[9, 9] = (Fc90) * (Fc90 * v300 + Fc91 * v310 + Fc92 * v320) + (Fc91) * (Fc90 * v301 + Fc91 * v311 + Fc92 * v321) + (Fc92) * (Fc90 * v302 + Fc91 * v312 + Fc92 * v322);

                    for (int pp = 0; pp < 10; pp++)
                    {
                        for (int ppp = 0; ppp < 10; ppp++)
                        {
                            S[6 * image_name_Eterior.Length + pp, 6 * image_name_Eterior.Length + ppp] = S[6 * image_name_Eterior.Length + pp, 6 * image_name_Eterior.Length + ppp] - FcvFc[pp, ppp];
                        }
                    }
                }
            }
            S_out = S; eawveb_out = eawveb; lllnc_out = lllnc; llln_out = llln;
        }
        public void Value_tie_point2(string[] image_name_Eterior, double[,] Exterior_Orientation_Eterior, string[] inv_observation_array, List<string> point_id, double[] interior, List<double> X_export, List<double> Y_export, List<double> Z_export, int camnumber, out double[,] cam_observe, out List<string> point_id2, out double[,] tie_point, out List<int> count_camera)
        {
            int shomtest = 0;
            int point_number = point_id.Count;
            List<string> pointid = new List<string>();
            string[] im_name_observe = new string[2];
            tie_point = new double[point_number, 3];
            if (camnumber < 2)
            {
                camnumber = 2;
            }
            double max_X = X_export.Max();
            double max_Y = Y_export.Max();
            double max_Z = Z_export.Max();
            double min_X = X_export.Min();
            double min_Y = Y_export.Min();
            double min_Z = Z_export.Min();
            List<double> ZZ = new List<double>();
           // MessageBox.Show("" + max_X + "    " + max_Y + "    " + max_Z + "    " + min_X + "    " + min_Y + "    " + min_Z);
            List<string> image_name_Eterior_list = image_name_Eterior.ToList();
            int znum = 0;
            double zmean = 0;
            for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
            {
                zmean = zmean + Exterior_Orientation_Eterior[j2, 2];
                znum = znum + 1;
            }
            zmean = zmean / znum;
            int znum2 = 0;
            double zmean2 = 0;
            string[] tokens = inv_observation_array[0].Split(new[] { "  " }, StringSplitOptions.None);
            string[] tokens2 = tokens[0].Split(new[] { " " }, StringSplitOptions.None);
            string point0 = tokens2[0];
            List<double> observe = new List<double>();
            List<double> Exterior = new List<double>();
            int count = 0;
            int pointnumber = 0;
            cam_observe = new double[inv_observation_array.Length, 3];
             count_camera = new List<int>();
            
            for (int i = 0; i < inv_observation_array.Length; i++)
            {
                string point1 = point0;
                tokens = inv_observation_array[i].Split(new[] { "  " }, StringSplitOptions.None);
                 tokens2 = tokens[0].Split(new[] { " " }, StringSplitOptions.None);
              string point = tokens2[0];
           
                string[] tokens3 = tokens[1].Split(new[] { " " }, StringSplitOptions.None);

                int j2 =Convert.ToInt32( tokens2[1]);
                // MessageBox.Show("" + Convert.ToDouble(observation_x_y[image_name_Eterior[1] + " " + point_id[48] + " " + "x"])+"  "+ point_id[47]);
                if (point0 == point)
                {
                    if (count == 0)
                    {
                       
                        observe.Add(Convert.ToDouble(tokens3[0]) / 1000);
                        observe.Add(Convert.ToDouble(tokens3[1]) / 1000);
                        cam_observe[i, 0] = j2;
                        cam_observe[i, 1] = Convert.ToDouble(tokens3[0]) / 1000;
                        cam_observe[i, 2] = Convert.ToDouble(tokens3[1]) / 1000;
                        for (int j = 0; j < 15; j++)
                        {

                            Exterior.Add(Exterior_Orientation_Eterior[j2, j]);


                        }
                        count = count + 1;
                    }
                    else
                    {
                        int shomar = 0;
                        for (int ij = 0; ij < count; ij++)
                        {
                            double difx = Exterior[15 * ij] - (Exterior_Orientation_Eterior[j2, 0]);
                            double dify = Exterior[15 * ij + 1] - (Exterior_Orientation_Eterior[j2, 1]);
                            double difz = Exterior[15 * ij + 2] - (Exterior_Orientation_Eterior[j2, 2]);
                            if ((Math.Pow((difx* difx + dify* dify + difz* difz), 1.0 / 2)) < 0)
                            {
                                shomar = shomar + 1;
                            }
                        }

                        if (shomar == 0)
                        {
                            observe.Add(Convert.ToDouble(tokens3[0]) / 1000);
                            observe.Add(Convert.ToDouble(tokens3[1]) / 1000);
                            for (int j = 0; j < 15; j++)
                            {
                                Exterior.Add(Exterior_Orientation_Eterior[j2, j]);
                            }
                            cam_observe[i, 0] = j2;
                            cam_observe[i, 1] = Convert.ToDouble(tokens3[0]) / 1000;
                            cam_observe[i, 2] = Convert.ToDouble(tokens3[1]) / 1000;
                            count = count + 1;
                        }
                    }
                    if(i== inv_observation_array.Length - 1)
                    {
                      
                        if (count >= camnumber)
                        {



                            double[,] XYZ1 = new double[3, 1];
                            XYZ_tie(observe, Exterior, interior[0],out XYZ1);
                            tie_point[pointnumber, 0] = XYZ1[0, 0];
                            tie_point[pointnumber, 1] = XYZ1[1, 0];
                            tie_point[pointnumber, 2] = XYZ1[2, 0];
                            ZZ.Add(XYZ1[2, 0]);
                            znum2 = znum2 + 1;
                            zmean2 = zmean2 + XYZ1[2, 0];

                            pointid.Add( point1);
                            pointnumber += 1;
                            count_camera.Add(count);

                        }
                        else
                        {
                            tie_point[pointnumber, 0] = 0; tie_point[pointnumber, 1] = 0; tie_point[pointnumber, 2] = 0;
                            pointid.Add("non");
                            pointnumber += 1;
                            shomtest += 1;
                            count_camera.Add(count);
                        }
                    }
                    continue;
                }
                else
                {
                    point0 = point;
                    i = i - 1;
                }
                
                
                if (count >= camnumber)
                {

                    double[,] XYZ1 = new double[3, 1];
                    XYZ_tie(observe, Exterior, interior[0],out XYZ1);
                    //Thread tr = new Thread(new ThreadStart(() => XYZ_tie(observe, Exterior, interior[0], out XYZ1)));
                    //tr.Start();
                    //tr.Join();
                    //   XYZ1 = XYZ_tie(observe, Exterior, interior[0]);
                    tie_point[pointnumber, 0] = XYZ1[0, 0];
                        tie_point[pointnumber, 1] = XYZ1[1, 0];
                        tie_point[pointnumber, 2] = XYZ1[2, 0];
                    ZZ.Add(XYZ1[2, 0]);
                    znum2 = znum2 + 1;
                    zmean2 = zmean2 + XYZ1[2, 0];
                    pointnumber += 1;
                    pointid.Add( point1);
                    count_camera.Add(count);

                }
                else
                {
                    tie_point[pointnumber, 0] = 0; tie_point[pointnumber, 1] = 0; tie_point[pointnumber, 2] = 0;
                    pointid.Add( "non");
                    pointnumber += 1;
                    shomtest += 1;
                    count_camera.Add(count);
                }
                 observe = new List<double>();
                 Exterior = new List<double>();
                count = 0;
            }
            
            zmean2 = zmean2 / znum2;
            point_id2 = pointid;
         //   MessageBox.Show("" + ZZ.Max() + "     " + ZZ.Min());
            // MessageBox.Show("" + initial_X_Y_Z[Array.FindIndex(point_id2_Observation, s => s.Equals("14326")),0]+"    "+ initial_X_Y_Z[Array.FindIndex(point_id2_Observation, s => s.Equals("14326")), 1]+"     "+ initial_X_Y_Z[Array.FindIndex(point_id2_Observation, s => s.Equals("14326")), 2]);
        }
        private void XYZ_tie(List<double> observe, List<double> Exterior, double f,out double [,] out_XYZ)
        {
      
            int num = observe.Count / 2;
            double[,] A = new double[2 * num, 3];
            double[,] l = new double[2 * num, 1];
            for (int i = 0; i < num; i++)
            {
                double XL = Exterior[15 * i];
                double YL = Exterior[15 * i + 1];
                double ZL = Exterior[15 * i + 2];
                //..........................
                //double m11 =Math.Cos(Exterior[15 * i + 4])*Math.Cos(Exterior[15 * i + 5]);
                //double m12 = (Math.Sin(Exterior[15 * i + 3]) * Math.Sin(Exterior[15 * i + 4])*Math.Cos(Exterior[15 * i + 5]))+(Math.Cos(Exterior[15 * i + 3]) * Math.Sin(Exterior[15 * i + 5]));
                //////m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                //double m13 = (Math.Sin(Exterior[15 * i + 3]) * Math.Sin(Exterior[15 * i + 5])) - (Math.Sin(Exterior[15 * i + 4]) * Math.Cos(Exterior[15 * i + 3]) * Math.Cos(Exterior[15 * i + 5]));
                //////m21=-cos(phi)*sin(k)
                //double m21 = -(Math.Cos(Exterior[15 * i + 4]) * Math.Sin(Exterior[15 * i + 5]));
                //////m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
             //   double aam22 = (Math.Cos(Exterior[15 * i + 3]) * Math.Cos(Exterior[15 * i + 5])) - (Math.Sin(Exterior[15 * i + 4]) * Math.Sin(Exterior[15 * i + 3]) * Math.Sin(Exterior[15 * i + 5]));
               
                //////m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))

                //double m23 = (Math.Sin(Exterior[15 * i + 3]) * Math.Cos(Exterior[15 * i + 5])) + (Math.Sin(Exterior[15 * i + 4]) * Math.Cos(Exterior[15 * i + 3]) * Math.Sin(Exterior[15 * i + 5]));
                //////m31=sin(phi)
                //double m31 = Math.Sin(Exterior[15 * i + 4]);
                //////m32=-sin(om)*cos(phi)
                //double m32 = -Math.Sin(Exterior[15 * i + 3]) * Math.Cos(Exterior[15 * i + 4]);
                //////m33=cos(om)*cos(phi)
                //double m33 = Math.Cos(Exterior[15 * i + 3]) * Math.Cos(Exterior[15 * i + 4]);
                //.........................

                double m11 = Exterior[15 * i + 6];
                double m12 = Exterior[15 * i + 7];
                double m13 = Exterior[15 * i + 8];
                double m21 = Exterior[15 * i + 9];
                double m22 = Exterior[15 * i + 10];
                double m23 = Exterior[15 * i + 11];
                double m31 = Exterior[15 * i + 12];
                double m32 = Exterior[15 * i + 13];
                double m33 = Exterior[15 * i + 14];
                //  MessageBox.Show("" + XL+"  "+YL + "  " + ZL + "  " + m11 + "  " + m12 + "  " +m13 + "  " + m21 + "  " + m22 + "  " + m23 + "  " + m31 + "  " + m32 + "  " +m33);
                A[2 * i, 0] = (observe[2 * i] * m31) + (f * m11);
                A[2 * i, 1] = (observe[2 * i] * m32) + (f * m12);
                A[2 * i, 2] = (observe[2 * i] * m33) + (f * m13);
                A[2 * i + 1, 0] = (observe[2 * i+1] * m31) + (f * m21);
                A[2 * i + 1, 1] = (observe[2 * i+1] * m32) + (f * m22);
                A[2 * i + 1, 2] = (observe[2 * i+1] * m33) + (f * m23);
                l[2 * i, 0] = (observe[2 * i] * ((XL * m31) + (YL * m32) + (ZL * m33))) + (f * ((XL * m11) + (YL * m12) + (ZL * m13)));
                l[2 * i + 1, 0] = (observe[2 * i + 1] * ((XL * m31) + (YL * m32) + (ZL * m33))) + (f * ((XL * m21) + (YL * m22) + (ZL * m23)));
            }
           
            // MessageBox.Show("" + AAr[0, 0]+"            "+ AAr[0, 1] + "            " + AAr[0, 2] + "            " + AAr[1, 0] + "            " + AAr[1, 1] + "            " + AAr[1, 2] + AAr[2, 0] + "            " + AAr[2, 1] + "            " + AAr[2, 2]);
            double[,] XYZ = new double[3, 1];
            double[,] invA = Matrix.Inverse(Matrix.Multiply(Matrix.Transpose(A), A));

            if (num < 3 & invA[2, 2] > 7000000)
            {
                int ppo = 1;
            }
            XYZ = Matrix.Multiply(Matrix.Multiply(invA, Matrix.Transpose(A)), l);

            out_XYZ = new double[3, 1];
            out_XYZ = XYZ;
            ///.......................
            ///
            //int ii = 0;
            //double am11 = Math.Cos(Exterior[15 * ii + 4]) * Math.Cos(Exterior[15 * ii + 5]);
            //double am12 = (Math.Sin(Exterior[15 * ii + 3]) * Math.Sin(Exterior[15 * ii + 4]) * Math.Cos(Exterior[15 * ii + 5])) + (Math.Cos(Exterior[15 * ii + 3]) * Math.Sin(Exterior[15 * ii + 5]));
            //////m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
            //double am13 = (Math.Sin(Exterior[15 * ii + 3]) * Math.Sin(Exterior[15 * ii + 5])) - (Math.Sin(Exterior[15 * ii + 4]) * Math.Cos(Exterior[15 * ii + 3]) * Math.Cos(Exterior[15 * ii + 5]));
            //////m21=-cos(phi)*sin(k)
            //double am21 = -(Math.Cos(Exterior[15 * ii + 4]) * Math.Sin(Exterior[15 * ii + 5]));
            //////m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
            //double am22 = (Math.Cos(Exterior[15 * ii + 4]) * Math.Cos(Exterior[15 * ii + 5])) - (Math.Sin(Exterior[15 * ii + 4]) * Math.Sin(Exterior[15 * ii + 3]) * Math.Sin(Exterior[15 * ii + 5]));
            //////m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
            //double am23 = (Math.Sin(Exterior[15 * ii + 3]) * Math.Cos(Exterior[15 * ii + 5])) + (Math.Sin(Exterior[15 * ii + 4]) * Math.Cos(Exterior[15 * ii + 3]) * Math.Sin(Exterior[15 * ii + 5]));
            //////m31=sin(phi)
            //double am31 = Math.Sin(Exterior[15 * ii + 4]);
            //////m32=-sin(om)*cos(phi)
            //double am32 = -Math.Sin(Exterior[15 * ii + 3]) * Math.Cos(Exterior[15 * ii + 4]);
            //////m33=cos(om)*cos(phi)
            //double am33 = Math.Cos(Exterior[15 * ii + 3]) * Math.Cos(Exterior[15 * ii + 4]);
            //double xx = -f * (((am11 * (XYZ[0,0] - Exterior[15 * ii])) + (am12 * (XYZ[1, 0] - Exterior[15 * ii + 1])) + (am13 * (XYZ[2, 0] - Exterior[15 * ii + 2]))) / ((am31 * (XYZ[0, 0] - Exterior[15 * ii])) + (am32 * (XYZ[1, 0] - Exterior[15 * ii + 1])) + (am33 * (XYZ[2, 0] - Exterior[15 * ii + 2]))));
            //double yy = -f * (((am21 * (XYZ[0, 0] - Exterior[15 * ii])) + (am22 * (XYZ[1, 0] - Exterior[15 * ii + 1])) + (am23 * (XYZ[2, 0] - Exterior[15 * ii + 2]))) / ((am31 * (XYZ[0, 0] - Exterior[15 * ii])) + (am32 * (XYZ[1, 0] - Exterior[15 * ii + 1])) + (am33 * (XYZ[2, 0] - Exterior[15 * ii + 2]))));

        }
        public static double[,] cov_Exterior_orientation_Image;
        public void Bundle_adjustment_fast(double[,] Initial_value, string[] GCP_point_name, double[,] GCP_observation, string[] GCP_point_name_weight, double[,] GCP_weight, string[] image_name_Eterior, double[,] Exterior_Orientation_Eterior, Hashtable GPS, double[,] cam_observe, List<string> point_id, double[] interior, double[,] GPS_weight, List<int> count_camera, out double[,] tie_points_coordinate, out double[,] Exterior_orientation, out double[,] cov_Exterior_orientation, out double[,] tie_out, out double[] Average_sigma)
        {
          
            Average_sigma = new double[6];
            List<int> GCP = new List<int>();
            List<string> GCP_in = new List<string>();
            List<string> GCP_we = new List<string>();
            Matrix_github S_g = new Matrix_github(6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length));
            Matrix_github eawveb_g = new Matrix_github(6 * (image_name_Eterior.Length), 1);
            if (GCP_point_name.Length > 1)
            {


                for (int i = 0; i < GCP_point_name.Length; i++)
                {
                    GCP_in.Add(GCP_point_name[i]);

                }
                for (int i = 0; i < GCP_point_name_weight.Length; i++)
                {
                    GCP_we.Add(GCP_point_name_weight[i]);
                }
            }
            double[,] LL = new double[6 * image_name_Eterior.Length, 1];

            int point_in = 0;

            List<double> camera_RMSE_x = new List<double>();
            List<double> camera_RMSE_y = new List<double>();
            List<double> camera_weight = new List<double>();

            Hashtable observe_value = new Hashtable();
            double[,] observe_value_camera = new double[image_name_Eterior.Length, 6];
            double[,] observe_value_points = new double[point_id.Count, 3];
            double[,] observe_value_points_weight = new double[point_id.Count, 3];
            double[,] cov_Ex = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
            for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
            {
                observe_value_camera[j2, 0] = Exterior_Orientation_Eterior[j2, 3];
                observe_value_camera[j2, 1] = Exterior_Orientation_Eterior[j2, 4];
                observe_value_camera[j2, 2] = Exterior_Orientation_Eterior[j2, 5];
                if (Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "X"]) != 0)
                {
                    observe_value_camera[j2, 3] = Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "X"]);
                }
                else
                {
                    observe_value_camera[j2, 3] = Exterior_Orientation_Eterior[j2, 0];
                }
                if (Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Y"]) != 0)
                {
                    observe_value_camera[j2, 4] = Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Y"]);
                }
                else
                {
                    observe_value_camera[j2, 4] = Exterior_Orientation_Eterior[j2, 1];
                }
                if (Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Z"]) != 0)
                {
                    observe_value_camera[j2, 5] = Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Z"]);
                }
                else
                {
                    observe_value_camera[j2, 5] = Exterior_Orientation_Eterior[j2, 2];
                }
                camera_weight.Add(1.0 / (interior[11] * interior[11]));
                camera_weight.Add(1.0 / (interior[11] * interior[11]));
                camera_weight.Add(1.0 / (GPS_weight[0, 0] * GPS_weight[0, 0]));
                camera_weight.Add(1.0 / (GPS_weight[0, 1] * GPS_weight[0, 1]));
                camera_weight.Add(1.0 / (GPS_weight[0, 2] * GPS_weight[0, 2]));
                camera_weight.Add(1.0 / (GPS_weight[0, 3] * GPS_weight[0, 3]));
                camera_weight.Add(1.0 / (GPS_weight[0, 4] * GPS_weight[0, 4]));
                camera_weight.Add(1.0 / (GPS_weight[0, 5] * GPS_weight[0, 5]));
            }

            for (int i = 0; i < point_id.Count; i++)
            {
                if (GCP_in.FindIndex(yy => yy == point_id[i]) != -1 && GCP_we.FindIndex(yy => yy == point_id[i]) != -1)
                {
                    observe_value_points_weight[i, 0] = 1.0 / (GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 0] * GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 0]);
                    observe_value_points_weight[i, 1] = 1.0 / (GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 1] * GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 1]);
                    observe_value_points_weight[i, 2] = 1.0 / (GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 2] * GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 2]);
                    observe_value_points[i, 0] = GCP_observation[(GCP_in.FindIndex(yy => yy == point_id[i])), 0];
                    observe_value_points[i, 1] = GCP_observation[(GCP_in.FindIndex(yy => yy == point_id[i])), 1];
                    observe_value_points[i, 2] = GCP_observation[(GCP_in.FindIndex(yy => yy == point_id[i])), 2];
                }
                else if (point_id[i] != "non")
                {
                    observe_value_points[i, 0] = Initial_value[i, 0];
                    observe_value_points[i, 1] = Initial_value[i, 1];
                    observe_value_points[i, 2] = Initial_value[i, 2];
                    observe_value_points_weight[i, 0] = 0.00000001;
                    observe_value_points_weight[i, 1] = 0.00000001;
                    observe_value_points_weight[i, 2] = 0.00000001;
                }
            }
            for (int iteration = 0; iteration < 1; iteration++)
            {


                //int iteration = 0;



                double[,] image_m = new double[image_name_Eterior.Length, 9];
                double[,] image_sin_cos = new double[image_name_Eterior.Length, 6];
                for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
                {
                    double om = Exterior_Orientation_Eterior[j2, 3];
                    double phi = Exterior_Orientation_Eterior[j2, 4];
                    double k = Exterior_Orientation_Eterior[j2, 5];
                    image_sin_cos[j2, 0] = Math.Sin(om); image_sin_cos[j2, 1] = Math.Sin(phi); image_sin_cos[j2, 2] = Math.Sin(k);
                    image_sin_cos[j2, 3] = Math.Cos(om); image_sin_cos[j2, 4] = Math.Cos(phi); image_sin_cos[j2, 5] = Math.Cos(k);
                    double sinom = Math.Sin(om); double sinphi = Math.Sin(phi); double sink = Math.Sin(k);
                    double cosom = Math.Cos(om); double cosphi = Math.Cos(phi); double cosk = Math.Cos(k);
                    //      MessageBox.Show("" + X+"   "+Y+"    "+Z);
                    image_m[j2, 0] = (cosphi) * (cosk);
                    image_m[j2, 1] = ((sinom) * (sinphi) * (cosk)) + ((cosom) * (sink));
                    //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                    image_m[j2, 2] = ((sinom) * (sink)) - ((sinphi) * (cosom) * (cosk));
                    //m21=-cos(phi)*sin(k)
                    image_m[j2, 3] = -(cosphi) * (sink);
                    //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                    image_m[j2, 4] = (cosom * cosk) - (sinom * sinphi * sink);
                    //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                    image_m[j2, 5] = (sinom * cosk) + (cosom * sinphi * sink);
                    //m31=sin(phi)
                    image_m[j2, 6] = (sinphi);
                    //m32=-sin(om)*cos(phi)
                    image_m[j2, 7] = -sinom * cosphi;
                    //m33=cos(om)*cos(phi)
                    image_m[j2, 8] = cosom * cosphi;
                }
                int shomwhile = 0;
                double[,] eawveb_out = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out2 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out2 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out2 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out3 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out3 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out3 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out4 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out4 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out4 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out5 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out5 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out5 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out6 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out6 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out6 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out7 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out7 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out7 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out8 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out8 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out8 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out9 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out9 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out9 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out10 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out10 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out10 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];

                // double[,] storeiv = new double[inv_observation_array.Length, 2];
                int odd = Convert.ToInt32(point_id.Count / 8);
                //  var s1 = System.Diagnostics.Stopwatch.StartNew();
                //    thread_bundle(0, point_id.Count, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out, out eawveb_out, out llln_out);

                //thread_bundle2(odd, odd * 2, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out2, out eawveb_out2, out llln_out2);

                //thread_bundle3(odd * 2, odd * 3, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out3, out eawveb_out3, out llln_out3);

                //thread_bundle4(odd * 3, point_id.Count, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out4, out eawveb_out4, out llln_out4);
                //   s1.Stop();
                //  MessageBox.Show("one " + s1.ElapsedMilliseconds);
                var s = System.Diagnostics.Stopwatch.StartNew();
                Thread t1 = new Thread(() => thread_bundle(0, odd, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out, out eawveb_out, out llln_out));

                Thread t2 = new Thread(() => thread_bundle(odd, odd * 2, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out2, out eawveb_out2, out llln_out2));

                Thread t3 = new Thread(() => thread_bundle(odd * 2, odd * 3, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out3, out eawveb_out3, out llln_out3));

                Thread t4 = new Thread(() => thread_bundle(odd * 3, odd * 4, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out4, out eawveb_out4, out llln_out4));
                Thread t5 = new Thread(() => thread_bundle(odd * 4, odd * 5, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out5, out eawveb_out5, out llln_out5));
                Thread t6 = new Thread(() => thread_bundle(odd * 5, odd * 6, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out6, out eawveb_out6, out llln_out6));
                Thread t7 = new Thread(() => thread_bundle(odd * 6, odd * 7, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out7, out eawveb_out7, out llln_out7));
                Thread t8 = new Thread(() => thread_bundle(odd * 7, point_id.Count, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out8, out eawveb_out8, out llln_out8));
                //    Thread t9 = new Thread(() => thread_bundle2(odd * 8, odd * 9, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out9, out eawveb_out9, out llln_out9));
                //     Thread t10 = new Thread(() => thread_bundle2(odd * 9, point_id.Count, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out10, out eawveb_out10, out llln_out10));

                t1.Start(); t2.Start(); t3.Start(); t4.Start(); t5.Start(); t6.Start(); t7.Start(); t8.Start();// t9.Start(); t10.Start();
                t1.Join(); t2.Join(); t3.Join(); t4.Join(); t5.Join(); t6.Join(); t7.Join(); t8.Join(); //t9.Join(); t10.Join();
                s.Stop();
             //   MessageBox.Show("two " + s.ElapsedMilliseconds);
                for (int i = 0; i < 6 * image_name_Eterior.Length; i++)
                {
                    for (int i2 = 0; i2 < 6 * image_name_Eterior.Length; i2++)
                    {
                        S_out[i, i2] = S_out[i, i2] + S_out2[i, i2] + S_out3[i, i2] + S_out4[i, i2] + S_out5[i, i2] + S_out6[i, i2] + S_out7[i, i2] + S_out8[i, i2] + S_out9[i, i2] + S_out10[i, i2];
                        S_g[i, i2] = S_out[i, i2];
                    }
                    eawveb_out[i, 0] = eawveb_out[i, 0] + eawveb_out2[i, 0] + eawveb_out3[i, 0] + eawveb_out4[i, 0] + eawveb_out5[i, 0] + eawveb_out6[i, 0] + eawveb_out7[i, 0] + eawveb_out8[i, 0] + eawveb_out9[i, 0] + eawveb_out10[i, 0];
                    //  llln_out[i, 0] = llln_out[i, 0] + llln_out2[i, 0] + llln_out3[i, 0] + llln_out4[i, 0]+ llln_out5[i, 0] + llln_out6[i, 0] + llln_out7[i, 0] + llln_out8[i, 0]+ llln_out9[i, 0] + llln_out10[i, 0];
                    eawveb_g[i, 0] = eawveb_out[i, 0];
                }
                for (int i = 0; i < image_name_Eterior.Length; i++)
                {
                    //eawveb_out[6 * i, 0] = eawveb_out[6 * i, 0] + eawveb_out2[6 * i, 0] + eawveb_out3[6 * i, 0] + eawveb_out4[6 * i, 0] + (camera_weight[8 * i + 2] * (llln_out[6 * i, 0]+ llln_out2[6 * i, 0]+ llln_out3[6 * i, 0]+ llln_out4[6 * i, 0]));
                    //eawveb_out[6 * i + 1, 0] = eawveb_out[6 * i + 1, 0] + eawveb_out2[6 * i + 1, 0] + eawveb_out3[6 * i + 1, 0] + eawveb_out4[6 * i + 1, 0] + (camera_weight[8 * i + 3] * (llln_out[6 * i + 1, 0] + llln_out2[6 * i + 1, 0] + llln_out3[6 * i + 1, 0] + llln_out4[6 * i + 1, 0]));
                    //eawveb_out[6 * i + 2, 0] = eawveb_out[6 * i + 2, 0] + eawveb_out2[6 * i + 2, 0] + eawveb_out3[6 * i + 2, 0] + eawveb_out4[6 * i + 2, 0] + (camera_weight[8 * i + 4] * (llln_out[6 * i + 2, 0] + llln_out2[6 * i + 2, 0] + llln_out3[6 * i + 2, 0] + llln_out4[6 * i + 2, 0]));
                    //eawveb_out[6 * i + 3, 0] = eawveb_out[6 * i + 3, 0] + eawveb_out2[6 * i + 3, 0] + eawveb_out3[6 * i + 3, 0] + eawveb_out4[6 * i + 3, 0] + (camera_weight[8 * i + 5] * (llln_out[6 * i + 3, 0] + llln_out2[6 * i + 3, 0] + llln_out3[6 * i + 3, 0] + llln_out4[6 * i + 3, 0]));
                    //eawveb_out[6 * i + 4, 0] = eawveb_out[6 * i + 4, 0] + eawveb_out2[6 * i + 4, 0] + eawveb_out3[6 * i + 4, 0] + eawveb_out4[6 * i + 4, 0] + (camera_weight[8 * i + 6] * (llln_out[6 * i + 4, 0] + llln_out2[6 * i + 4, 0] + llln_out3[6 * i + 4, 0] + llln_out4[6 * i + 4, 0]));
                    //eawveb_out[6 * i + 5, 0] = eawveb_out[6 * i + 5, 0] + eawveb_out2[6 * i + 5, 0] + eawveb_out3[6 * i + 5, 0] + eawveb_out4[6 * i + 5, 0] + (camera_weight[8 * i + 7] * (llln_out[6 * i + 5, 0] + llln_out2[6 * i + 5, 0] + llln_out3[6 * i + 5, 0] + llln_out4[6 * i + 5, 0]));

                    //S_out[6 * i, 6 * i] = S_out[6 * i, 6 * i] + S_out2[6 * i, 6 * i] + S_out3[6 * i, 6 * i] + S_out4[6 * i, 6 * i] + camera_weight[8 * i + 2];
                    //S_out[6 * i + 1, 6 * i + 1] = S_out[6 * i + 1, 6 * i + 1] + S_out2[6 * i + 1, 6 * i + 1] + S_out3[6 * i + 1, 6 * i + 1] + S_out4[6 * i + 1, 6 * i + 1] + camera_weight[8 * i + 3];
                    //S_out[6 * i + 2, 6 * i + 2] = S_out[6 * i + 2, 6 * i + 2] + S_out2[6 * i + 2, 6 * i + 2] + S_out3[6 * i + 2, 6 * i + 2] + S_out4[6 * i + 2, 6 * i + 2] + camera_weight[8 * i + 4];
                    //S_out[6 * i + 3, 6 * i + 3] = S_out[6 * i + 3, 6 * i + 3] + S_out2[6 * i + 3, 6 * i + 3] + S_out3[6 * i + 3, 6 * i + 3] + S_out4[6 * i + 3, 6 * i + 3] + camera_weight[8 * i + 5];
                    //S_out[6 * i + 4, 6 * i + 4] = S_out[6 * i + 4, 6 * i + 4] + S_out2[6 * i + 4, 6 * i + 4] + S_out3[6 * i + 4, 6 * i + 4] + S_out4[6 * i + 4, 6 * i + 4] + camera_weight[8 * i + 6];
                    //S_out[6 * i + 5, 6 * i + 5] = S_out[6 * i + 5, 6 * i + 5] + S_out2[6 * i + 5, 6 * i + 5] + S_out3[6 * i + 5, 6 * i + 5] + S_out4[6 * i + 5, 6 * i + 5] + camera_weight[8 * i + 7];
                    eawveb_out[6 * i, 0] = eawveb_out[6 * i, 0] + (camera_weight[8 * i + 2] * llln_out[6 * i, 0]);
                    eawveb_out[6 * i + 1, 0] = eawveb_out[6 * i + 1, 0] + (camera_weight[8 * i + 3] * llln_out[6 * i + 1, 0]);
                    eawveb_out[6 * i + 2, 0] = eawveb_out[6 * i + 2, 0] + (camera_weight[8 * i + 4] * llln_out[6 * i + 2, 0]);
                    eawveb_out[6 * i + 3, 0] = eawveb_out[6 * i + 3, 0] + (camera_weight[8 * i + 5] * llln_out[6 * i + 3, 0]);
                    eawveb_out[6 * i + 4, 0] = eawveb_out[6 * i + 4, 0] + (camera_weight[8 * i + 6] * llln_out[6 * i + 4, 0]);
                    eawveb_out[6 * i + 5, 0] = eawveb_out[6 * i + 5, 0] + (camera_weight[8 * i + 7] * llln_out[6 * i + 5, 0]);
                    //   MessageBox.Show("si" + S[6 * i+5, 6 * i+5]+"    "+ camera_weight[8 * i + 5]);
                    S_out[6 * i, 6 * i] = S_out[6 * i, 6 * i] + camera_weight[8 * i + 2];
                    S_out[6 * i + 1, 6 * i + 1] = S_out[6 * i + 1, 6 * i + 1] + camera_weight[8 * i + 3];
                    S_out[6 * i + 2, 6 * i + 2] = S_out[6 * i + 2, 6 * i + 2] + camera_weight[8 * i + 4];
                    S_out[6 * i + 3, 6 * i + 3] = S_out[6 * i + 3, 6 * i + 3] + camera_weight[8 * i + 5];
                    S_out[6 * i + 4, 6 * i + 4] = S_out[6 * i + 4, 6 * i + 4] + camera_weight[8 * i + 6];
                    S_out[6 * i + 5, 6 * i + 5] = S_out[6 * i + 5, 6 * i + 5] + camera_weight[8 * i + 7];
                 
                }

           
                ////inv matinv = new inv();
                ////dx_calculator.Class1 sgh = new dx_calculator.Class1();
                ////MWArray[] SMW = new MWNumericArray[] { S_out };
                ////MWNumericArray sdf = new MWNumericArray(S_out);
                ////MWNumericArray sdfl = new MWNumericArray(eawveb_out);
                ////MWArray invS2 = matinv.inv2(sdf);
                ////MWArray dx2 = sgh.dx_calculator(sdf, sdfl);
                ////double[,] invS = (double[,])((MWNumericArray)invS2).ToArray(MWArrayComponent.Real);
                ////double[,] dx = (double[,])((MWNumericArray)dx2).ToArray(MWArrayComponent.Real);
                double[,] invS = Matrix.Inverse(S_out);
                double[,] dx = Matrix.Multiply(invS, eawveb_out);
               
            //    System.IO.StreamWriter ll1l = new System.IO.StreamWriter(@"E:\test\LL.txt");
            //    System.IO.StreamWriter ll2l22 = new System.IO.StreamWriter(@"E:\test\invS.txt");
               
             //   for (int i = 0; i < 6 * image_name_Eterior.Length; i++)
             //   {
             //       ll1l.WriteLine("" + dx[i, 0]);
             //   }
             //   ll1l.Close();
                double s_omega = 0, s_phi = 0, s_k = 0, s_XL = 0, s_YL = 0, s_ZL = 0;
                for (int i = 0; i < image_name_Eterior.Length; i++)
                {
                    s_omega = s_omega + Math.Sqrt(invS[6 * i, 6 * i]);
                    s_phi = s_phi + Math.Sqrt(invS[6 * i + 1, 6 * i + 1]);
                    s_k = s_k + Math.Sqrt(invS[6 * i + 2, 6 * i + 2]);
                    s_XL = s_XL + Math.Sqrt(invS[6 * i + 3, 6 * i + 3]);
                    s_YL = s_YL + Math.Sqrt(invS[6 * i + 4, 6 * i + 4]);
                    s_ZL = s_ZL + Math.Sqrt(invS[6 * i + 5, 6 * i + 5]);
               //     ll2l22.WriteLine("" + image_name_Eterior[i] + " " + Math.Sqrt(invS[6 * i, 6 * i]) + " " + Math.Sqrt(invS[6 * i + 1, 6 * i + 1]) + " " + Math.Sqrt(invS[6 * i + 2, 6 * i + 2]) + " " + Math.Sqrt(invS[6 * i + 3, 6 * i + 3]) + " " + Math.Sqrt(invS[6 * i + 4, 6 * i + 4]) + " " + Math.Sqrt(invS[6 * i + 5, 6 * i + 5]));
                }

                Average_sigma[0] = s_omega / image_name_Eterior.Length;
                Average_sigma[1] = s_phi / image_name_Eterior.Length;
                Average_sigma[2] = s_k / image_name_Eterior.Length;
                Average_sigma[3] = s_XL / image_name_Eterior.Length;
                Average_sigma[4] = s_YL / image_name_Eterior.Length;
                Average_sigma[5] = s_ZL / image_name_Eterior.Length;
            //    ll2l22.Close();
                cov_Ex = invS;
         // ........................................................................................................................................
         //from there is diactivated for Error propagation
                //shomwhile = 0;
                //for (int i = 0; i < point_id.Count; i++)
                //{


                //    double[,] v = new double[3, 3];
                //    double[,] v2 = new double[3, 3];
                //    double[,] ww = new double[6, 3];
                //    List<double> w1 = new List<double>();
                //    List<double> w2 = new List<double>();
                //    List<double> w3 = new List<double>();

                //    double[,] ebwdx = new double[3, 1];
                //    List<int> wj = new List<int>();
                //    //  double[,] pp2 = new double[21, 21];
                //    int ij2 = shomwhile;
                //    shomwhile += count_camera[i];
                //    if (point_id[i] != "non")
                //    {

                //        double[,] pp2 = new double[2, 2];
                //        double[,] ebii = new double[3, 1];
                //        double[,] lllnb = new double[3, 1];

                //        for (int ij = ij2; ij < shomwhile; ij++)
                //        {

                //            int j2 = Convert.ToInt32(cam_observe[ij, 0]);

                //            pp2[0, 0] = camera_weight[8 * j2];
                //            pp2[1, 1] = camera_weight[8 * j2 + 1];

                //            double XL = Exterior_Orientation_Eterior[j2, 0];
                //            double YL = Exterior_Orientation_Eterior[j2, 1];
                //            double ZL = Exterior_Orientation_Eterior[j2, 2];
                //            double om = Exterior_Orientation_Eterior[j2, 3];
                //            double phi = Exterior_Orientation_Eterior[j2, 4];
                //            double k = Exterior_Orientation_Eterior[j2, 5];
                //            double ff = interior[0];
                //            double x0 = interior[1];
                //            double y0 = interior[2];
                //            double K1 = interior[3];
                //            double K2 = interior[4];
                //            double K3 = interior[5];
                //            double P1 = interior[6];
                //            double P2 = interior[7];
                //            double B1 = interior[8];
                //            double B2 = interior[9];
                //            double X = Initial_value[i, 0];
                //            double Y = Initial_value[i, 1];
                //            double Z = Initial_value[i, 2];
                //            double m11 = image_m[j2, 0];
                //            double m12 = image_m[j2, 1];
                //            //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                //            double m13 = image_m[j2, 2];
                //            //m21=-cos(phi)*sin(k)
                //            double m21 = image_m[j2, 3];
                //            //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                //            double m22 = image_m[j2, 4];
                //            //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                //            double m23 = image_m[j2, 5];
                //            //m31=sin(phi)
                //            double m31 = image_m[j2, 6];
                //            //m32=-sin(om)*cos(phi)
                //            double m32 = image_m[j2, 7];
                //            //m33=cos(om)*cos(phi)
                //            double m33 = image_m[j2, 8];

                //            double sinom = image_sin_cos[j2, 0]; double sinphi = image_sin_cos[j2, 1]; double sink = image_sin_cos[j2, 2];
                //            double cosom = image_sin_cos[j2, 3]; double cosphi = image_sin_cos[j2, 4]; double cosk = image_sin_cos[j2, 5];

                //            double xx = cam_observe[ij, 1];
                //            double yy = cam_observe[ij, 2];
                //            double sinphiXXL = (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL));
                //            double sinksinom = (sink * sinom - cosk * cosom * sinphi);
                //            double cosomsink = (cosom * sink + cosk * sinom * sinphi);
                //            double coskcosom = (cosk * cosom - sink * sinom * sinphi);
                //            double cosksinom = (cosk * sinom + cosom * sink * sinphi);
                //            double xom = ((ff * (sinksinom * (Y - YL) - cosomsink * (Z - ZL))) / sinphiXXL - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));

                //            double xphi = ((ff * (cosk * sinphi * (X - XL) + cosk * cosom * cosphi * (Z - ZL) - cosk * cosphi * sinom * (Y - YL))) / sinphiXXL + (ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double xk = (-(ff * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / sinphiXXL);
                //            double xXL = ((ff * cosk * cosphi) / sinphiXXL - (ff * sinphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double xYL = ((ff * cosomsink) / sinphiXXL + (ff * cosphi * sinom * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double xZL = (-((ff * cosom * cosphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * sinksinom) / sinphiXXL));
                //            //   (f * cos(om) * cos(phi) * ((cos(om) * sin(k) + cos(k) * sin(om) * sin(phi)) * (Y - YL) + (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi)) * (Z - ZL) + cos(k) * cos(phi) * (X - XL))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL)) ^ 2 - (f * (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL))];
                //            double yom = ((ff * (cosksinom * (Y - YL) - coskcosom * (Z - ZL))) / sinphiXXL - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double yphi = ((ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * (sink * sinphi * (X - XL) + cosom * cosphi * sink * (Z - ZL) - cosphi * sink * sinom * (Y - YL))) / sinphiXXL);
                //            double yk = ((ff * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / sinphiXXL);
                //            double yXL = (-(ff * sinphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosphi * sink) / sinphiXXL);
                //            double yYL = ((ff * coskcosom) / sinphiXXL + (ff * cosphi * sinom * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double yZL = ((ff * cosksinom) / sinphiXXL - (ff * cosom * cosphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double xX = ((ff * sinphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosk * cosphi) / sinphiXXL);
                //            double xY = (-(ff * cosomsink) / sinphiXXL - (ff * cosphi * sinom * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double xZ = ((ff * cosom * cosphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * sinksinom) / sinphiXXL);
                //            double yX = ((ff * sinphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) + (ff * cosphi * sink) / sinphiXXL);
                //            double yY = (-(ff * coskcosom) / sinphiXXL - (ff * cosphi * sinom * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
                //            double yZ = ((ff * cosom * cosphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosksinom) / sinphiXXL);


                //            double xpi = (xx - x0);
                //            double ypi = (yy - y0);

                //            double ri = Math.Sqrt(((xpi * xpi) + (ypi * ypi)));
                //            double deltax = (xpi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P1 * (((ri * ri)) + (2 * ((xpi * xpi))))) + (2 * P2 * xpi * ypi) + (B1 * xpi) + (B2 * ypi);
                //            double deltay = (ypi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P2 * (((ri * ri)) + (2 * ((ypi * ypi))))) + (2 * P1 * xpi * ypi);
                //            double x = -(ff * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltax;
                //            double y = -(ff * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltay;







                //            double a00 = xom;

                //            double a01 = xphi;

                //            double a02 = xk;

                //            double a03 = xXL;

                //            double a04 = xYL;

                //            double a05 = xZL;

                //            double a10 = yom;

                //            double a11 = yphi;

                //            double a12 = yk;

                //            double a13 = yXL;

                //            double a14 = yYL;

                //            double a15 = yZL;
                //            double bb00 = xX;

                //            double bb01 = xY;

                //            double bb02 = xZ;

                //            double bb10 = yX;

                //            double bb11 = yY;

                //            double bb12 = yZ;
                //            double pb00 = pp2[0, 0]; double pb01 = 0; double pb10 = 0; double pb11 = pp2[1, 1];

                //            //  MessageBox.Show("" + c[0, 0] + "   " + c[0, 1]+ "   " + c[0, 2] + "   " + c[0, 3] + "   " + c[0, 4] + "   " + c[0, 5] + "   " + c[0, 6] + "   " + c[0, 7] + "   " + c[0, 8] + "   " + c[0, 9]);
                //            // MessageBox.Show("" + c[1, 0] + "   " + c[1, 1] + "   " + c[1, 2] + "   " + c[1, 3] + "   " + c[1, 4] + "   " + c[1, 5] + "   " + c[1, 6] + "   " + c[1, 7] + "   " + c[1, 8] + "   " + c[1, 9]);

                //            lllnb[0, 0] = observe_value_points[i, 0] - X; lllnb[1, 0] = observe_value_points[i, 1] - Y; lllnb[2, 0] = observe_value_points[i, 2] - Z;
                //            // MessageBox.Show("" + lll[0, 0]+"   " + lll[1, 0]);
                //            double ll0 = xx - x; ; double ll1 = yy - y;
                //            double[,] ebii2 = new double[3, 1];
                //            ebii2[0, 0] = ll0 * (pb00 * (bb00) + pb10 * (bb10)) + ll1 * (pb01 * (bb00) + pb11 * (bb10));
                //            ebii2[1, 0] = ll0 * (pb00 * (bb01) + pb10 * (bb11)) + ll1 * (pb01 * (bb01) + pb11 * (bb11));
                //            ebii2[2, 0] = ll0 * (pb00 * (bb02) + pb10 * (bb12)) + ll1 * (pb01 * (bb02) + pb11 * (bb12));
                //            ebii[0, 0] = ebii[0, 0] + ebii2[0, 0]; ebii[1, 0] = ebii[1, 0] + ebii2[1, 0];
                //            ebii[2, 0] = ebii[2, 0] + ebii2[2, 0];
                //            //  MessageBox.Show("" + ebii[0, 0] + "   " + ebii2[0, 0]);
                //            v2[0, 0] = bb00 * (pb00 * (bb00) + pb10 * (bb10)) + bb10 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 1] = bb01 * (pb00 * (bb00) + pb10 * (bb10)) + bb11 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 2] = bb02 * (pb00 * (bb00) + pb10 * (bb10)) + bb12 * (pb01 * (bb00) + pb11 * (bb10));
                //            v2[1, 0] = bb00 * (pb00 * (bb01) + pb10 * (bb11)) + bb10 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 1] = bb01 * (pb00 * (bb01) + pb10 * (bb11)) + bb11 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 2] = bb02 * (pb00 * (bb01) + pb10 * (bb11)) + bb12 * (pb01 * (bb01) + pb11 * (bb11));
                //            v2[2, 0] = bb00 * (pb00 * (bb02) + pb10 * (bb12)) + bb10 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 1] = bb01 * (pb00 * (bb02) + pb10 * (bb12)) + bb11 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 2] = bb02 * (pb00 * (bb02) + pb10 * (bb12)) + bb12 * (pb01 * (bb02) + pb11 * (bb12));


                //            ww[0, 0] = bb00 * (pb00 * (a00) + pb10 * (a10)) + bb10 * (pb01 * (a00) + pb11 * (a10)); ww[0, 1] = bb01 * (pb00 * (a00) + pb10 * (a10)) + bb11 * (pb01 * (a00) + pb11 * (a10)); ww[0, 2] = bb02 * (pb00 * (a00) + pb10 * (a10)) + bb12 * (pb01 * (a00) + pb11 * (a10));
                //            ww[1, 0] = bb00 * (pb00 * (a01) + pb10 * (a11)) + bb10 * (pb01 * (a01) + pb11 * (a11)); ww[1, 1] = bb01 * (pb00 * (a01) + pb10 * (a11)) + bb11 * (pb01 * (a01) + pb11 * (a11)); ww[1, 2] = bb02 * (pb00 * (a01) + pb10 * (a11)) + bb12 * (pb01 * (a01) + pb11 * (a11));
                //            ww[2, 0] = bb00 * (pb00 * (a02) + pb10 * (a12)) + bb10 * (pb01 * (a02) + pb11 * (a12)); ww[2, 1] = bb01 * (pb00 * (a02) + pb10 * (a12)) + bb11 * (pb01 * (a02) + pb11 * (a12)); ww[2, 2] = bb02 * (pb00 * (a02) + pb10 * (a12)) + bb12 * (pb01 * (a02) + pb11 * (a12));
                //            ww[3, 0] = bb00 * (pb00 * (a03) + pb10 * (a13)) + bb10 * (pb01 * (a03) + pb11 * (a13)); ww[3, 1] = bb01 * (pb00 * (a03) + pb10 * (a13)) + bb11 * (pb01 * (a03) + pb11 * (a13)); ww[3, 2] = bb02 * (pb00 * (a03) + pb10 * (a13)) + bb12 * (pb01 * (a03) + pb11 * (a13));
                //            ww[4, 0] = bb00 * (pb00 * (a04) + pb10 * (a14)) + bb10 * (pb01 * (a04) + pb11 * (a14)); ww[4, 1] = bb01 * (pb00 * (a04) + pb10 * (a14)) + bb11 * (pb01 * (a04) + pb11 * (a14)); ww[4, 2] = bb02 * (pb00 * (a04) + pb10 * (a14)) + bb12 * (pb01 * (a04) + pb11 * (a14));
                //            ww[5, 0] = bb00 * (pb00 * (a05) + pb10 * (a15)) + bb10 * (pb01 * (a05) + pb11 * (a15)); ww[5, 1] = bb01 * (pb00 * (a05) + pb10 * (a15)) + bb11 * (pb01 * (a05) + pb11 * (a15)); ww[5, 2] = bb02 * (pb00 * (a05) + pb10 * (a15)) + bb12 * (pb01 * (a05) + pb11 * (a15));
                //            double ww00 = ww[0, 0]; double ww01 = ww[0, 1]; double ww02 = ww[0, 2];
                //            double ww10 = ww[1, 0]; double ww11 = ww[1, 1]; double ww12 = ww[1, 2];
                //            double ww20 = ww[2, 0]; double ww21 = ww[2, 1]; double ww22 = ww[2, 2];
                //            double ww30 = ww[3, 0]; double ww31 = ww[3, 1]; double ww32 = ww[3, 2];
                //            double ww40 = ww[4, 0]; double ww41 = ww[4, 1]; double ww42 = ww[4, 2];
                //            double ww50 = ww[5, 0]; double ww51 = ww[5, 1]; double ww52 = ww[5, 2];
                //            double[,] dxx2 = new double[6, 1];
                //            for (int pp = 0; pp < 6; pp++)
                //            {
                //                dxx2[0, 0] = dx[6 * j2, 0];
                //                dxx2[1, 0] = dx[6 * j2 + 1, 0];
                //                dxx2[2, 0] = dx[6 * j2 + 2, 0];
                //                dxx2[3, 0] = dx[6 * j2 + 3, 0];
                //                dxx2[4, 0] = dx[6 * j2 + 4, 0];
                //                dxx2[5, 0] = dx[6 * j2 + 5, 0];
                //            }
                //            double dxx200 = dxx2[0, 0]; double dxx210 = dxx2[1, 0]; double dxx220 = dxx2[2, 0];
                //            double dxx230 = dxx2[3, 0]; double dxx240 = dxx2[4, 0]; double dxx250 = dxx2[5, 0];
                //            //  double[,] dxw = Matrix.Multiply(Matrix.Transpose(ww), dxx2);
                //            double[,] dxw = new double[3, 1];
                //            dxw[0, 0] = dxx200 * (ww00) + dxx210 * (ww10) + dxx220 * (ww20) + dxx230 * (ww30) + dxx240 * (ww40) + dxx250 * (ww50);
                //            dxw[1, 0] = dxx200 * (ww01) + dxx210 * (ww11) + dxx220 * (ww21) + dxx230 * (ww31) + dxx240 * (ww41) + dxx250 * (ww51);
                //            dxw[2, 0] = dxx200 * (ww02) + dxx210 * (ww12) + dxx220 * (ww22) + dxx230 * (ww32) + dxx240 * (ww42) + dxx250 * (ww52);


                //            for (int ppp = 0; ppp < 3; ppp++)
                //            {
                //                ebwdx[ppp, 0] = ebwdx[ppp, 0] - dxw[ppp, 0];
                //            }

                //            w1.Add(ww[0, 0]);
                //            w1.Add(ww[1, 0]);
                //            w1.Add(ww[2, 0]);
                //            w1.Add(ww[3, 0]);
                //            w1.Add(ww[4, 0]);
                //            w1.Add(ww[5, 0]);

                //            w2.Add(ww[0, 1]);
                //            w2.Add(ww[1, 1]);
                //            w2.Add(ww[2, 1]);
                //            w2.Add(ww[3, 1]);
                //            w2.Add(ww[4, 1]);
                //            w2.Add(ww[5, 1]);

                //            w3.Add(ww[0, 2]);
                //            w3.Add(ww[1, 2]);
                //            w3.Add(ww[2, 2]);
                //            w3.Add(ww[3, 2]);
                //            w3.Add(ww[4, 2]);
                //            w3.Add(ww[5, 2]);

                //            wj.Add(j2);
                //            //WW[6 * j2 + 0, 0] = ww[0, 0]; WW[6 * j2 + 1, 0] = ww[1, 0]; WW[6 * j2 + 2, 0] = ww[2, 0]; WW[6 * j2 + 3, 0] = ww[3, 0]; WW[6 * j2 + 4, 0] = ww[4, 0]; WW[6 * j2 + 5, 0] = ww[5, 0];
                //            //WW[6 * j2 + 0, 1] = ww[0, 1]; WW[6 * j2 + 1, 1] = ww[1,1]; WW[6 * j2 + 2, 1] = ww[2, 1]; WW[6 * j2 + 3, 1] = ww[3, 1]; WW[6 * j2 + 4, 1] = ww[4, 1]; WW[6 * j2 + 5, 1] = ww[5, 1];
                //            //WW[6 * j2 + 0, 2] = ww[0, 2]; WW[6 * j2 + 1, 2] = ww[1, 2]; WW[6 * j2 + 2, 2] = ww[2, 2]; WW[6 * j2 + 3, 2] = ww[3,2]; WW[6 * j2 + 4,2] = ww[4, 2]; WW[6 * j2 + 5, 2] = ww[5, 2];
                //            v[0, 0] = v[0, 0] + v2[0, 0]; v[1, 0] = v[1, 0] + v2[1, 0]; v[2, 0] = v[2, 0] + v2[2, 0];
                //            v[0, 1] = v[0, 1] + v2[0, 1]; v[1, 1] = v[1, 1] + v2[1, 1]; v[2, 1] = v[2, 1] + v2[2, 1];
                //            v[0, 2] = v[0, 2] + v2[0, 2]; v[1, 2] = v[1, 2] + v2[1, 2]; v[2, 2] = v[2, 2] + v2[2, 2];
                //        }

                //        double[,] pp2nB = new double[3, 3];
                //        //    if (i == 1||i==2||i==3||i==4||i==5||i==6||i==7)
                //        //  {
                //        //   pp2nB[0, 0] = 10000000000000000;
                //        //  pp2nB[1, 1] = 100000000000000000;
                //        //  pp2nB[2, 2] = 10000000000000000;
                //        //  }
                //        //  else
                //        //   {
                //        pp2nB[0, 0] = observe_value_points_weight[i, 0];
                //        pp2nB[1, 1] = observe_value_points_weight[i, 1];
                //        pp2nB[2, 2] = observe_value_points_weight[i, 2];
                //        //   }
                //        v[0, 0] = v[0, 0] + pp2nB[0, 0]; v[1, 1] = v[1, 1] + pp2nB[1, 1]; v[2, 2] = v[2, 2] + pp2nB[2, 2];

                //        ebii[0, 0] = ebii[0, 0] + (pp2nB[0, 0] * lllnb[0, 0]);
                //        ebii[1, 0] = ebii[1, 0] + (pp2nB[1, 1] * lllnb[1, 0]);
                //        ebii[2, 0] = ebii[2, 0] + (pp2nB[2, 2] * lllnb[2, 0]);
                //        for (int ppp = 0; ppp < 3; ppp++)
                //        {
                //            ebwdx[ppp, 0] = ebwdx[ppp, 0] + ebii[ppp, 0];
                //        }
                //        double ebwdx00 = ebwdx[0, 0]; double ebwdx10 = ebwdx[1, 0]; double ebwdx20 = ebwdx[2, 0];

                //        //     MessageBox.Show("" + dxx2[0, 0]);
                //        //  MessageBox.Show("   " + Fdx[0, 0] + "    " +Fdx[1, 0] + "    " + Fdx[2, 0] );


                //        //  MessageBox.Show("" + "   " + Fdx[0, 0] + "   ");
                //        //var s1 = Stopwatch.StartNew();
                //        double[,] v3 = new double[3, 3];
                //        double detv = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                //        v3[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                //        v3[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                //        v3[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                //        v3[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                //        v3[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                //        v3[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                //        v3[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                //        v3[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                //        v3[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                //        double v300 = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                //        double v301 = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                //        double v302 = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                //        double v310 = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                //        double v311 = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                //        double v312 = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                //        double v320 = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                //        double v321 = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                //        double v322 = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;

                //        double[,] dxb = new double[3, 1];
                //        dxb[0, 0] = ebwdx00 * v300 + ebwdx10 * v301 + ebwdx20 * v302;
                //        dxb[1, 0] = ebwdx00 * v310 + ebwdx10 * v311 + ebwdx20 * v312;
                //        dxb[2, 0] = ebwdx00 * v320 + ebwdx10 * v321 + ebwdx20 * v322;
                //        //     if (dxb[0, 0]>100|| dxb[1, 0] > 100|| dxb[2, 0] > 10)
                //        //   {
                //        //     MessageBox.Show("   " + dxb[0, 0] + "    " + dxb[1, 0] + "    " + dxb[2, 0]+"   "+point_id[i] + "   " + v[0,0] + "   " + v[0, 1] + "   " + v[0, 2] + "   " + v[1, 0] + "   " + v[1, 1] + "   " + v[1, 2] + "   " + v[2, 0] + "   " + v[2, 1] + "   " + v[2, 2]);
                //        // point_id[i] = "non";
                //        //   }
                //        // s1.Stop();
                //        // MessageBox.Show("sina" + s1.ElapsedMilliseconds);
                //        point_in = point_in + 1;
                //        Initial_value[i, 0] = Initial_value[i, 0] + dxb[0, 0];
                //        Initial_value[i, 1] = Initial_value[i, 1] + dxb[1, 0];
                //        Initial_value[i, 2] = Initial_value[i, 2] + dxb[2, 0];
                //    }

                //}
                for (int i = 0; i < image_name_Eterior.Length; i++)
                {

                    Exterior_Orientation_Eterior[i, 0] = Exterior_Orientation_Eterior[i, 0] + dx[6 * i + 3, 0];
                    Exterior_Orientation_Eterior[i, 1] = Exterior_Orientation_Eterior[i, 1] + dx[6 * i + 4, 0];
                    Exterior_Orientation_Eterior[i, 2] = Exterior_Orientation_Eterior[i, 2] + dx[6 * i + 5, 0];
                    Exterior_Orientation_Eterior[i, 3] = Exterior_Orientation_Eterior[i, 3] + (dx[6 * i, 0]);
                    Exterior_Orientation_Eterior[i, 4] = Exterior_Orientation_Eterior[i, 4] + (dx[6 * i + 1, 0]);
                    Exterior_Orientation_Eterior[i, 5] = Exterior_Orientation_Eterior[i, 5] + (dx[6 * i + 2, 0]);
                }


                //     MessageBox.Show("" + dx[0, 0]);


            }

            tie_out = Initial_value;
            point_in = point_in - 1;
            tie_points_coordinate = Initial_value; Exterior_orientation = Exterior_Orientation_Eterior; cov_Exterior_orientation = cov_Ex;
            cov_Exterior_orientation_Image = new double[cov_Exterior_orientation.GetLength(0), cov_Exterior_orientation.GetLength(1)];
            cov_Exterior_orientation_Image = cov_Exterior_orientation;
        }


        public void Bundle_adjust_exterior_and_interior(double[,] Initial_value, string[] GCP_point_name, double[,] GCP_observation, string[] GCP_point_name_weight, double[,] GCP_weight, string[] image_name_Eterior, double[,] Exterior_Orientation_Eterior, Hashtable GPS, double[,] cam_observe, List<string> point_id, double[] interior, double[,] GPS_weight, List<int> count_camera, out double[,] tie_points_coordinate, out double[,] Exterior_orientation, out double[,] cov_Exterior_orientation, out double[] Average_sigma)
        {
            Average_sigma = new double[6];

            double[,] L = new double[1, 1];

            List<int> GCP = new List<int>();
            List<string> GCP_in = new List<string>();
            List<string> GCP_we = new List<string>();
            if (GCP_point_name.Length > 1)
            {
                for (int i = 0; i < GCP_point_name.Length; i++)
                {
                    GCP_in.Add(GCP_point_name[i]);

                }
                for (int i = 0; i < GCP_point_name_weight.Length; i++)
                {

                    GCP_we.Add(GCP_point_name_weight[i]);
                }
            }
            double[,] observe_value_camera = new double[image_name_Eterior.Length, 6];
            double[,] observe_value_points = new double[point_id.Count, 3];
            double[,] observe_value_points_weight = new double[point_id.Count, 3];
            List<double> camera_weight = new List<double>();
            for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
            {
                observe_value_camera[j2, 0] = Exterior_Orientation_Eterior[j2, 3];
                observe_value_camera[j2, 1] = Exterior_Orientation_Eterior[j2, 4];
                observe_value_camera[j2, 2] = Exterior_Orientation_Eterior[j2, 5];
                if (Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "X"]) != 0)
                {
                    observe_value_camera[j2, 3] = Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "X"]);
                }
                else
                {
                    observe_value_camera[j2, 3] = Exterior_Orientation_Eterior[j2, 0];
                }
                if (Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Y"]) != 0)
                {
                    observe_value_camera[j2, 4] = Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Y"]);
                }
                else
                {
                    observe_value_camera[j2, 4] = Exterior_Orientation_Eterior[j2, 1];
                }
                if (Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Z"]) != 0)
                {
                    observe_value_camera[j2, 5] = Convert.ToDouble(GPS[image_name_Eterior[j2] + " " + "Z"]);
                }
                else
                {
                    observe_value_camera[j2, 5] = Exterior_Orientation_Eterior[j2, 2];
                }
                camera_weight.Add(1.0 / (interior[11] * interior[11]));
                camera_weight.Add(1.0 / (interior[11] * interior[11]));
                camera_weight.Add(1.0 / (GPS_weight[0, 0] * GPS_weight[0, 0]));
                camera_weight.Add(1.0 / (GPS_weight[0, 1] * GPS_weight[0, 1]));
                camera_weight.Add(1.0 / (GPS_weight[0, 2] * GPS_weight[0, 2]));
                camera_weight.Add(1.0 / (GPS_weight[0, 3] * GPS_weight[0, 3]));
                camera_weight.Add(1.0 / (GPS_weight[0, 4] * GPS_weight[0, 4]));
                camera_weight.Add(1.0 / (GPS_weight[0, 5] * GPS_weight[0, 5]));
            }

            for (int i = 0; i < point_id.Count; i++)
            {
                if (GCP_in.FindIndex(yy => yy == point_id[i]) != -1 && GCP_we.FindIndex(yy => yy == point_id[i]) != -1)
                {
                    observe_value_points_weight[i, 0] = 1.0 / (GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 0] * GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 0]);
                    observe_value_points_weight[i, 1] = 1.0 / (GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 1] * GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 1]);
                    observe_value_points_weight[i, 2] = 1.0 / (GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 2] * GCP_weight[(GCP_we.FindIndex(yy => yy == point_id[i])), 2]);
                    observe_value_points[i, 0] = GCP_observation[(GCP_in.FindIndex(yy => yy == point_id[i])), 0];
                    observe_value_points[i, 1] = GCP_observation[(GCP_in.FindIndex(yy => yy == point_id[i])), 1];
                    observe_value_points[i, 2] = GCP_observation[(GCP_in.FindIndex(yy => yy == point_id[i])), 2];
                }
                else if (point_id[i] != "non")
                {
                    observe_value_points[i, 0] = Initial_value[i, 0];
                    observe_value_points[i, 1] = Initial_value[i, 1];
                    observe_value_points[i, 2] = Initial_value[i, 2];
                    observe_value_points_weight[i, 0] = 0.00000001;
                    observe_value_points_weight[i, 1] = 0.00000001;
                    observe_value_points_weight[i, 2] = 0.00000001;
                }
            }
            double[,] LL = new double[6 * image_name_Eterior.Length, 1];

            int point_in = 0;

            List<double> camera_RMSE_x = new List<double>();
            List<double> camera_RMSE_y = new List<double>();
            double[,] IOP = new double[10, 1];
            IOP[0, 0] = interior[0];
           
     

           

            double[,] cov_Ex = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
            double ff0 = interior[0];
            double x00 = interior[1];
            double y00 = interior[2];
            double K10 = interior[3];
            double K20 = interior[4];
            double K30 = interior[5];
            double P10 = interior[6];
            double P20 = interior[7];
            double B10 = interior[8];
            double B20 = interior[9];
            
   
           
            for (int iteration = 0; iteration < 20; iteration++)
                {
                double[,] image_m = new double[image_name_Eterior.Length, 9];
                double[,] image_sin_cos = new double[image_name_Eterior.Length, 6];
                for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
                {
                    double om = Exterior_Orientation_Eterior[j2, 3];
                    double phi = Exterior_Orientation_Eterior[j2, 4];
                    double k = Exterior_Orientation_Eterior[j2, 5];
                    image_sin_cos[j2, 0] = Math.Sin(om); image_sin_cos[j2, 1] = Math.Sin(phi); image_sin_cos[j2, 2] = Math.Sin(k);
                    image_sin_cos[j2, 3] = Math.Cos(om); image_sin_cos[j2, 4] = Math.Cos(phi); image_sin_cos[j2, 5] = Math.Cos(k);
                    double sinom = Math.Sin(om); double sinphi = Math.Sin(phi); double sink = Math.Sin(k);
                    double cosom = Math.Cos(om); double cosphi = Math.Cos(phi); double cosk = Math.Cos(k);
                    //      MessageBox.Show("" + X+"   "+Y+"    "+Z);
                    image_m[j2, 0] = (cosphi) * (cosk);
                    image_m[j2, 1] = ((sinom) * (sinphi) * (cosk)) + ((cosom) * (sink));
                    //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                    image_m[j2, 2] = ((sinom) * (sink)) - ((sinphi) * (cosom) * (cosk));
                    //m21=-cos(phi)*sin(k)
                    image_m[j2, 3] = -(cosphi) * (sink);
                    //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                    image_m[j2, 4] = (cosom * cosk) - (sinom * sinphi * sink);
                    //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                    image_m[j2, 5] = (sinom * cosk) + (cosom * sinphi * sink);
                    //m31=sin(phi)
                    image_m[j2, 6] = (sinphi);
                    //m32=-sin(om)*cos(phi)
                    image_m[j2, 7] = -sinom * cosphi;
                    //m33=cos(om)*cos(phi)
                    image_m[j2, 8] = cosom * cosphi;
                }

                double[,] llln_out = new double[6 * (image_name_Eterior.Length), 1];
                double[,] lllnc_out = new double[10, 1];

                double[,] pp2nC = new double[10, 10];
                pp2nC[0, 0] =  1/10000000000000000000;
                pp2nC[1, 1] = 1/10000000000000000000;
                pp2nC[2, 2] = 1/10000000000000000000;
                pp2nC[3, 3] = 1/10000000000000000000;
                pp2nC[4, 4] = 1/10000000000000000000;
                pp2nC[5, 5] = 1/10000000000000000000;
                pp2nC[6, 6] = 1/10000000000000000000;
                pp2nC[7, 7] = 1/10000000000000000000;
                pp2nC[8, 8] = 1/10000000000000000000;
                pp2nC[9, 9] = 1/ 10000000000000000000;
                double[,] eawveb_out = new double[6 * (image_name_Eterior.Length), 1];

                double[,] S_out = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out2 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out2 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out2 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out3 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out3 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out3 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out4 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out4 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out4 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out5 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out5 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out5 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out6 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out6 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out6 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out7 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out7 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out7 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out8 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out8 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out8 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out9 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out9 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out9 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];
                double[,] eawveb_out10 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] llln_out10 = new double[6 * (image_name_Eterior.Length), 1];
                double[,] S_out10 = new double[6 * (image_name_Eterior.Length), 6 * (image_name_Eterior.Length)];

                // double[,] storeiv = new double[inv_observation_array.Length, 2];
                int odd = Convert.ToInt32(point_id.Count / 8);
             

                
                //Thread t1 = new Thread(() => thread_bundle_interior_exterior(0, odd, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out, out eawveb_out, out llln_out, out lllnc_out));
                //Thread t2 = new Thread(() => thread_bundle_interior_exterior(odd, odd * 2, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out2, out eawveb_out2, out llln_out, out lllnc_out));
                //Thread t3 = new Thread(() => thread_bundle_interior_exterior(odd * 2, odd * 3, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out3, out eawveb_out3, out llln_out, out lllnc_out));
                //Thread t4 = new Thread(() => thread_bundle_interior_exterior(odd * 3, odd * 4, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out4, out eawveb_out4, out llln_out, out lllnc_out));
                //Thread t5 = new Thread(() => thread_bundle_interior_exterior(odd * 4, odd * 5, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out5, out eawveb_out5, out llln_out, out lllnc_out));
                //Thread t6 = new Thread(() => thread_bundle_interior_exterior(odd * 5, odd * 6, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out6, out eawveb_out6, out llln_out, out lllnc_out));
                //Thread t7 = new Thread(() => thread_bundle_interior_exterior(odd * 6, odd * 7, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out7, out eawveb_out7, out llln_out, out lllnc_out));
                //Thread t8 = new Thread(() => thread_bundle_interior_exterior(odd * 7, point_id.Count, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out8, out eawveb_out8, out llln_out, out lllnc_out));
           
                //t1.Start(); t2.Start(); t3.Start(); t4.Start(); t5.Start(); t6.Start(); t7.Start(); t8.Start();// t9.Start(); t10.Start();
                //t1.Join(); t2.Join(); t3.Join(); t4.Join(); t5.Join(); t6.Join(); t7.Join(); t8.Join(); //t9.Join(); t10.Join();
              
           
                  thread_bundle_interior_exterior(0, point_id.Count, IOP, image_name_Eterior, point_id, count_camera, observe_value_points_weight, cam_observe, camera_weight, Exterior_Orientation_Eterior, interior, Initial_value, image_m, image_sin_cos, observe_value_camera, observe_value_points, out S_out, out eawveb_out, out llln_out, out lllnc_out);
                for (int i = 0; i < 6 * image_name_Eterior.Length; i++)
                {
                    for (int i2 = 0; i2 < 6 * image_name_Eterior.Length; i2++)
                    {
                        S_out[i, i2] = S_out[i, i2] + S_out2[i, i2] + S_out3[i, i2] + S_out4[i, i2] + S_out5[i, i2] + S_out6[i, i2] + S_out7[i, i2] + S_out8[i, i2] + S_out9[i, i2] + S_out10[i, i2];
                    }
                    eawveb_out[i, 0] = eawveb_out[i, 0] + eawveb_out2[i, 0] + eawveb_out3[i, 0] + eawveb_out4[i, 0] + eawveb_out5[i, 0] + eawveb_out6[i, 0] + eawveb_out7[i, 0] + eawveb_out8[i, 0] + eawveb_out9[i, 0] + eawveb_out10[i, 0];
                    //  llln_out[i, 0] = llln_out[i, 0] + llln_out2[i, 0] + llln_out3[i, 0] + llln_out4[i, 0]+ llln_out5[i, 0] + llln_out6[i, 0] + llln_out7[i, 0] + llln_out8[i, 0]+ llln_out9[i, 0] + llln_out10[i, 0];
                }

                for (int i = 0; i < image_name_Eterior.Length; i++)
                {
                    eawveb_out[6 * i, 0] = eawveb_out[6 * i, 0] - (camera_weight[8 * i + 2] * llln_out[6 * i, 0]);
                    eawveb_out[6 * i + 1, 0] = eawveb_out[6 * i + 1, 0] - (camera_weight[8 * i + 3] * llln_out[6 * i + 1, 0]);
                    eawveb_out[6 * i + 2, 0] = eawveb_out[6 * i + 2, 0] - (camera_weight[8 * i + 4] * llln_out[6 * i + 2, 0]);
                    eawveb_out[6 * i + 3, 0] = eawveb_out[6 * i + 3, 0] - (camera_weight[8 * i + 5] * llln_out[6 * i + 3, 0]);
                    eawveb_out[6 * i + 4, 0] = eawveb_out[6 * i + 4, 0] - (camera_weight[8 * i + 6] * llln_out[6 * i + 4, 0]);
                    eawveb_out[6 * i + 5, 0] = eawveb_out[6 * i + 5, 0] - (camera_weight[8 * i + 7] * llln_out[6 * i + 5, 0]);
                    //   MessageBox.Show("si" + S[6 * i+5, 6 * i+5]+"    "+ camera_weight[8 * i + 5]);
                    S_out[6 * i, 6 * i] = S_out[6 * i, 6 * i] + camera_weight[8 * i + 2];
                    S_out[6 * i + 1, 6 * i + 1] = S_out[6 * i + 1, 6 * i + 1] + camera_weight[8 * i + 3];
                    S_out[6 * i + 2, 6 * i + 2] = S_out[6 * i + 2, 6 * i + 2] + camera_weight[8 * i + 4];
                    S_out[6 * i + 3, 6 * i + 3] = S_out[6 * i + 3, 6 * i + 3] + camera_weight[8 * i + 5];
                    S_out[6 * i + 4, 6 * i + 4] = S_out[6 * i + 4, 6 * i + 4] + camera_weight[8 * i + 6];
                    S_out[6 * i + 5, 6 * i + 5] = S_out[6 * i + 5, 6 * i + 5] + camera_weight[8 * i + 7];

                }
                eawveb_out[6 * image_name_Eterior.Length, 0] = eawveb_out[6 * image_name_Eterior.Length, 0] - (pp2nC[1, 1] * lllnc_out[1, 0]);
                eawveb_out[6 * image_name_Eterior.Length+1, 0] = eawveb_out[6 * image_name_Eterior.Length+1, 0] - (pp2nC[2, 2] * lllnc_out[2, 0]);
                eawveb_out[6 * image_name_Eterior.Length+2, 0] = eawveb_out[6 * image_name_Eterior.Length+2, 0] - (pp2nC[3, 3] * lllnc_out[3, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 3, 0] = eawveb_out[6 * image_name_Eterior.Length + 3, 0] - (pp2nC[4, 4] * lllnc_out[4, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 4, 0] = eawveb_out[6 * image_name_Eterior.Length + 4, 0] - (pp2nC[5, 5] * lllnc_out[5, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 5, 0] = eawveb_out[6 * image_name_Eterior.Length + 5, 0] - (pp2nC[6, 6] * lllnc_out[6, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 6, 0] = eawveb_out[6 * image_name_Eterior.Length + 6, 0] - (pp2nC[7, 7] * lllnc_out[7, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 7, 0] = eawveb_out[6 * image_name_Eterior.Length + 7, 0] - (pp2nC[8, 8] * lllnc_out[8, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 8, 0] = eawveb_out[6 * image_name_Eterior.Length + 8, 0] - (pp2nC[9, 9] * lllnc_out[9, 0]);
                eawveb_out[6 * image_name_Eterior.Length + 9, 0] = eawveb_out[6 * image_name_Eterior.Length + 9, 0] - (pp2nC[0, 0] * lllnc_out[0, 0]);
                S_out[6 * image_name_Eterior.Length + 0, 6 * image_name_Eterior.Length + 0] = S_out[6 * image_name_Eterior.Length + 0, 6 * image_name_Eterior.Length + 0] + pp2nC[1, 1];
                S_out[6 * image_name_Eterior.Length + 1, 6 * image_name_Eterior.Length + 1] = S_out[6 * image_name_Eterior.Length + 1, 6 * image_name_Eterior.Length + 1] + pp2nC[2, 2];
                S_out[6 * image_name_Eterior.Length + 2, 6 * image_name_Eterior.Length + 2] = S_out[6 * image_name_Eterior.Length + 2, 6 * image_name_Eterior.Length + 2] + pp2nC[3, 3];
                S_out[6 * image_name_Eterior.Length + 3, 6 * image_name_Eterior.Length + 3] = S_out[6 * image_name_Eterior.Length + 3, 6 * image_name_Eterior.Length + 3] + pp2nC[4, 4];
                S_out[6 * image_name_Eterior.Length + 4, 6 * image_name_Eterior.Length + 4] = S_out[6 * image_name_Eterior.Length + 4, 6 * image_name_Eterior.Length + 4] + pp2nC[5, 5];
                S_out[6 * image_name_Eterior.Length + 5, 6 * image_name_Eterior.Length + 5] = S_out[6 * image_name_Eterior.Length + 5, 6 * image_name_Eterior.Length + 5] + pp2nC[6, 6];
                S_out[6 * image_name_Eterior.Length + 6, 6 * image_name_Eterior.Length + 6] = S_out[6 * image_name_Eterior.Length + 6, 6 * image_name_Eterior.Length + 6] + pp2nC[7, 7];
                S_out[6 * image_name_Eterior.Length + 7, 6 * image_name_Eterior.Length + 7] = S_out[6 * image_name_Eterior.Length + 7, 6 * image_name_Eterior.Length + 7] + pp2nC[8, 8];
                S_out[6 * image_name_Eterior.Length + 8, 6 * image_name_Eterior.Length + 8] = S_out[6 * image_name_Eterior.Length + 8, 6 * image_name_Eterior.Length + 8] + pp2nC[9, 9];
                S_out[6 * image_name_Eterior.Length + 9, 6 * image_name_Eterior.Length + 9] = S_out[6 * image_name_Eterior.Length + 9, 6 * image_name_Eterior.Length + 9] + pp2nC[0, 0];
                System.IO.StreamWriter lS = new System.IO.StreamWriter(@"E:\test\S.txt");
                System.IO.StreamWriter leawveb = new System.IO.StreamWriter(@"E:\test\eawveb.txt");
                for (int i = 6 * image_name_Eterior.Length; i < 6 * image_name_Eterior.Length+10; i++)
                {
                    string line = "" + S_out[i, 0];
                    for (int j = 6 * image_name_Eterior.Length+1; j < 6 * image_name_Eterior.Length+10; j++)
                    {
                        line = line + " " + S_out[i, j];
                    }
                    lS.WriteLine(line);
                    leawveb.WriteLine("" + eawveb_out[i, 0]);

                }
                lS.Close();
                leawveb.Close();
                inv matinv = new inv();
                dx_calculator.Class1 sgh = new dx_calculator.Class1();
                MWArray[] SMW = new MWNumericArray[] { S_out };
                MWNumericArray sdf = new MWNumericArray(S_out);
                MWNumericArray sdfl = new MWNumericArray(eawveb_out);
                MWArray invS2 = matinv.inv2(sdf);
                MWArray dx2 = sgh.dx_calculator(sdf, sdfl);
                double[,] invS = (double[,])((MWNumericArray)invS2).ToArray(MWArrayComponent.Real);
                double[,] dx = (double[,])((MWNumericArray)dx2).ToArray(MWArrayComponent.Real);
                //  MessageBox.Show("" + dx[0, 0]+ "   "+eawveb[0,0]);
                //double[,] invS = Matrix.Inverse(S_out);
                //double[,] dx = Matrix.Multiply(invS, eawveb_out);
                System.IO.StreamWriter ll1l = new System.IO.StreamWriter(@"E:\test\LL.txt");
                System.IO.StreamWriter ll2l22 = new System.IO.StreamWriter(@"E:\test\invS.txt");
                //   MessageBox.Show("" + dx[0, 0]);
                for (int i = 0; i < 6 * image_name_Eterior.Length+10; i++)
                {
                    ll1l.WriteLine("" + dx[i, 0]);
                }
                
                ll1l.Close();
                double s_omega = 0, s_phi = 0, s_k = 0, s_XL = 0, s_YL = 0, s_ZL = 0;
                for (int i = 0; i < image_name_Eterior.Length; i++)
                {
                    s_omega = s_omega + Math.Sqrt(invS[6 * i, 6 * i]);
                    s_phi = s_phi + Math.Sqrt(invS[6 * i + 1, 6 * i + 1]);
                    s_k = s_k + Math.Sqrt(invS[6 * i + 2, 6 * i + 2]);
                    s_XL = s_XL + Math.Sqrt(invS[6 * i + 3, 6 * i + 3]);
                    s_YL = s_YL + Math.Sqrt(invS[6 * i + 4, 6 * i + 4]);
                    s_ZL = s_ZL + Math.Sqrt(invS[6 * i + 5, 6 * i + 5]);
                    ll2l22.WriteLine("" + image_name_Eterior[i] + " " + Math.Sqrt(invS[6 * i, 6 * i]) + " " + Math.Sqrt(invS[6 * i + 1, 6 * i + 1]) + " " + Math.Sqrt(invS[6 * i + 2, 6 * i + 2]) + " " + Math.Sqrt(invS[6 * i + 3, 6 * i + 3]) + " " + Math.Sqrt(invS[6 * i + 4, 6 * i + 4]) + " " + Math.Sqrt(invS[6 * i + 5, 6 * i + 5]));
                }

                Average_sigma[0] = s_omega / image_name_Eterior.Length;
                Average_sigma[1] = s_phi / image_name_Eterior.Length;
                Average_sigma[2] = s_k / image_name_Eterior.Length;
                Average_sigma[3] = s_XL / image_name_Eterior.Length;
                Average_sigma[4] = s_YL / image_name_Eterior.Length;
                Average_sigma[5] = s_ZL / image_name_Eterior.Length;

                ll2l22.Close();
                cov_Ex = invS;
   
                int shomwhile = 0;
                for (int i = 0; i < point_id.Count; i++)
                {
                  

                    double[,] v = new double[3, 3];
                    double[,] v2 = new double[3, 3];
                    double[,] ww = new double[6, 3];
                    double[,] wwt = new double[3, 3];
                    double[,] Fc = new double[10, 3];


                    List<double> w1 = new List<double>();
                    List<double> w2 = new List<double>();
                    List<double> w3 = new List<double>();

                    double[,] ebwdx = new double[3, 1];
                    List<int> wj = new List<int>();
                    //  double[,] pp2 = new double[21, 21];
                    int ij2 = shomwhile;
                    shomwhile += count_camera[i];
                    if (point_id[i] != "non")
                    {

                        double[,] pp2 = new double[2, 2];
                        double[,] ebii = new double[3, 1];
                        double[,] lllnb = new double[3, 1];

                           
                        

                                for (int ij = ij2; ij < shomwhile; ij++)
                                {

                                    int j2 = Convert.ToInt32(cam_observe[ij, 0]);
                            pp2[0, 0] = camera_weight[8 * j2];
                            pp2[1, 1] = camera_weight[8 * j2 + 1];
                            double XL = Exterior_Orientation_Eterior[j2, 0];
                                double YL = Exterior_Orientation_Eterior[j2, 1];
                                double ZL = Exterior_Orientation_Eterior[j2, 2];
                                double om = Exterior_Orientation_Eterior[j2, 3];
                                double phi = Exterior_Orientation_Eterior[j2, 4];
                                double k = Exterior_Orientation_Eterior[j2, 5];
                                double ff = IOP[0, 0];
                                double x0 = IOP[1, 0] ;
                                double y0 = IOP[2, 0] ;
                                double K1 = IOP[3, 0];
                                double K2 = IOP[4, 0];
                                double K3 = IOP[5, 0];
                                double P1 = IOP[6, 0];
                                double P2 = IOP[7, 0];
                                double B1 = IOP[8, 0];
                                double B2 = IOP[9, 0];
                                double X = Initial_value[i, 0];
                                double Y = Initial_value[i, 1];
                                double Z = Initial_value[i, 2];
                                double sinom = Math.Sin(om); double sinphi = Math.Sin(phi); double sink = Math.Sin(k);
                                double cosom = Math.Cos(om); double cosphi = Math.Cos(phi); double cosk = Math.Cos(k);
                                //      MessageBox.Show("" + X+"   "+Y+"    "+Z);
                                double m11 = (cosphi) * (cosk);
                                double m12 = ((sinom) * (sinphi) * (cosk)) + ((cosom) * (sink));
                                //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                                double m13 = ((sinom) * (sink)) - ((sinphi) * (cosom) * (cosk));
                                //m21=-cos(phi)*sin(k)
                                double m21 = -(cosphi) * (sink);
                                //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                                double m22 = (cosom * cosk) - (sinom * sinphi * sink);
                                //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                                double m23 = (sinom * cosk) + (cosom * sinphi * sink);
                                //m31=sin(phi)
                                double m31 = (sinphi);
                                //m32=-sin(om)*cos(phi)
                                double m32 = -sinom * cosphi;
                                //m33=cos(om)*cos(phi)
                                double m33 = cosom * cosphi;

                            double xx = cam_observe[ij, 1];
                            double yy = cam_observe[ij, 2];
                            // 2 * P2 * (y0 - yy) - K1 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) - B1 - K2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2 - K3 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 3 + P1 * (6 * x0 - 6 * xx) - (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) + 3 * K3 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ 2) + 1
                            double Fx_xx = 2 * P2 * (y0 - yy) - K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - B1 - K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + P1 * (6 * x0 - 6 * xx) - (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) + 1;
                                double Fy_xx = 2 * P1 * (y0 - yy) + P2 * (2 * x0 - 2 * xx) - (y0 - yy) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))));
                                double Fx_yy = 2 * P2 * (x0 - xx) - B2 + P1 * (2 * y0 - 2 * yy) - (x0 - xx) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))));
                                double Fy_yy = 2 * P1 * (x0 - xx) - K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + P2 * (6 * y0 - 6 * yy) - (y0 - yy) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) + 1;



                                double xom = -((ff * ((sink * sinom - cosk * cosom * sinphi) * (Y - YL) - (cosom * sink + cosk * sinom * sinphi) * (Z - ZL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));

                                double xphi = -((ff * (cosk * sinphi * (X - XL) + cosk * cosom * cosphi * (Z - ZL) - cosk * cosphi * sinom * (Y - YL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) + (ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                                double xk = -(-(ff * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double xXL = -((ff * cosk * cosphi) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * sinphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                                double xYL = -((ff * (cosom * sink + cosk * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) + (ff * cosphi * sinom * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                                double xZL = -(-((ff * cosom * cosphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * (sink * sinom - cosk * cosom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))));
                                //   (f * cos(om) * cos(phi) * ((cos(om) * sin(k) + cos(k) * sin(om) * sin(phi)) * (Y - YL) + (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi)) * (Z - ZL) + cos(k) * cos(phi) * (X - XL))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL)) ^ 2 - (f * (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL))];
                                double yom = -((ff * ((cosk * sinom + cosom * sink * sinphi) * (Y - YL) - (cosk * cosom - sink * sinom * sinphi) * (Z - ZL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                                double yphi = -((ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2) - (ff * (sink * sinphi * (X - XL) + cosom * cosphi * sink * (Z - ZL) - cosphi * sink * sinom * (Y - YL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double yk = -((ff * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double yXL = -(-(ff * sinphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2) - (ff * cosphi * sink) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double yYL = -((ff * (cosk * cosom - sink * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) + (ff * cosphi * sinom * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                                double yZL = -((ff * (cosk * sinom + cosom * sink * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * cosom * cosphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / Math.Pow((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)), 2));
                                double xX = -((ff * sinphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * cosk * cosphi) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double xY = -(-(ff * (cosom * sink + cosk * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * cosphi * sinom * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))));
                                double xZ = -((ff * cosom * cosphi * ((cosom * sink + cosk * sinom * sinphi) * (Y - YL) + (sink * sinom - cosk * cosom * sinphi) * (Z - ZL) + cosk * cosphi * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * (sink * sinom - cosk * cosom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double yX = -((ff * sinphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) + (ff * cosphi * sink) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));
                                double yY = -(-(ff * (cosk * cosom - sink * sinom * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) - (ff * cosphi * sinom * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))));
                                double yZ = -((ff * cosom * cosphi * ((cosk * cosom - sink * sinom * sinphi) * (Y - YL) + (cosk * sinom + cosom * sink * sinphi) * (Z - ZL) - cosphi * sink * (X - XL))) / ((sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)) * (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL))) - (ff * (cosk * sinom + cosom * sink * sinphi)) / (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL)));


                                double xpi = (xx - x0);
                                double ypi = (yy - y0);

                                double ri = Math.Sqrt(((xpi * xpi) + (ypi * ypi)));
                                double deltax = (xpi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P1 * (((ri * ri)) + (2 * ((xpi * xpi))))) + (2 * P2 * xpi * ypi) + (B1 * xpi) + (B2 * ypi);
                                double deltay = (ypi * ((K1 * ((ri * ri))) + (K2 * ((ri * ri) * (ri * ri))) + (K3 * ((ri * ri) * (ri * ri) * (ri * ri))))) + (P2 * (((ri * ri)) + (2 * ((ypi * ypi))))) + (2 * P1 * xpi * ypi);
                                if (i == 1 && j2 == 0)
                                {
                                    //       MessageBox.Show("" + deltax + "  " + deltay);
                                }
                                double x = -(ff * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltax;
                                double y = -(ff * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) + deltay;
                                double[,] lll = new double[2, 1];
                                double Fx = (xpi) + (ff * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) - (deltax);
                                double Fy = (ypi) + (ff * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))))) - (deltay);
                                //    MessageBox.Show("" + deltax + "   " + deltay);
                                double Fxf = (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                                double Fyf = (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                                // 2*P2*(y0 - yy) - K1*((x0 - xx)^2 + (y0 - yy)^2) - B1 - K2*((x0 - xx)^2 + (y0 - yy)^2)^2 - K3*((x0 - xx)^2 + (y0 - yy)^2)^3 + P1*(6*x0 - 6*xx) - (x0 - xx)*(K1*(2*x0 - 2*xx) + 2*K2*(2*x0 - 2*xx)*((x0 - xx)^2 + (y0 - yy)^2) + 3*K3*(2*x0 - 2*xx)*((x0 - xx)^2 + (y0 - yy)^2)^2) - 1

                                double Fxx0 = B1 + K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - 2 * P2 * (y0 - yy) + K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - P1 * (6 * x0 - 6 * xx) + (x0 - xx) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) - 1;
                                double Fxy0 = B2 - 2 * P2 * (x0 - xx) - P1 * (2 * y0 - 2 * yy) + (x0 - xx) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))));

                                //           2 * P2 * (x0 - xx) - B2 + P1 * (2 * y0 - 2 * yy) - ((x0 - xx) * ((3 * K1 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (1 / 2)) / 2 + (5 * K2 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (3 / 2)) / 2 + (7 * K3 * (2 * y0 - 2 * yy) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (5 / 2)) / 2)) / ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (1 / 2) + ((2 * y0 - 2 * yy) * (x0 - xx) * (K1 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (3 / 2) + K2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (5 / 2) + K3 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (7 / 2))) / (2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (3 / 2))
                                //  double FxK1 = -(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2)))) * (x0 - xx);
                                double FxK1 = -xpi * (ri * ri);
                                // MessageBox.Show("" + FxK1 + "   " +xx+"    "+yy+"    "+ri+"      ");
                                // double FxK2 = -(Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (2))) * (x0 - xx);
                                double FxK2 = -xpi * (ri * ri) * (ri * ri);
                                //  double FxK3 = -(Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (3))) * (x0 - xx);
                                double FxK3 = -xpi * (ri * ri) * (ri * ri) * (ri * ri);
                                //  double FxP1 = -((x0 - xx) * (3 * (((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2)))))) / (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)));
                                double FxP1 = -(((ri * ri)) + (2 * (xpi * xpi)));
                                // double FxP2 = -(2 * Math.Pow((x0 - xx), 2) * (y0 - yy)) / (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)));
                                double FxP2 = -2 * xpi * ypi;
                                // double FxB1 = Math.Pow((x0 - xx), 2) / (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)));
                                double FxB1 = -xpi;
                                //  double FxB2 = ((x0 - xx) * (y0 - yy)) / (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)));
                                double FxB2 = -ypi;
                                // 2 * P1 * (y0 - yy) + P2 * (2 * x0 - 2 * xx) - ((y0 - yy) * ((3 * K1 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (1 / 2)) / 2 + (5 * K2 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (3 / 2)) / 2 + (7 * K3 * (2 * x0 - 2 * xx) * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (5 / 2)) / 2)) / ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (1 / 2) + ((2 * x0 - 2 * xx) * (y0 - yy) * (K1 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (3 / 2) + K2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (5 / 2) + K3 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (7 / 2))) / (2 * ((x0 - xx) ^ 2 + (y0 - yy) ^ 2) ^ (3 / 2))
                                double Fyx0 = (y0 - yy) * (K1 * (2 * x0 - 2 * xx) + 2 * K2 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * x0 - 2 * xx) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) - P2 * (2 * x0 - 2 * xx) - 2 * P1 * (y0 - yy);
                                double Fyy0 = K1 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - 2 * P1 * (x0 - xx) + K2 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + K3 * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) - P2 * (6 * y0 - 6 * yy) + (y0 - yy) * (K1 * (2 * y0 - 2 * yy) + 2 * K2 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))))) + 3 * K3 * (2 * y0 - 2 * yy) * ((((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy)))) * ((((x0 - xx) * (x0 - xx))) + (((y0 - yy) * (y0 - yy))))))) - 1;

                                // double FyK1 = -(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2)))) * (y0 - yy);
                                double FyK1 = -ypi * (ri * ri);
                                //  double FyK2 = -(Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (2))) * (y0 - yy);
                                double FyK2 = -ypi * (ri * ri) * (ri * ri);
                                // double FyK3 = -(Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (3))) * (y0 - yy);
                                double FyK3 = -ypi * (ri * ri) * (ri * ri) * (ri * ri);
                                // double FyP1 = -(2 * (x0 - xx) * Math.Pow((y0 - yy), 2)) / (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)));
                                double FyP1 = -2 * ypi * xpi;
                                // double FyP2 = -((y0 - yy) * (Math.Pow((x0 - xx), 2) + 3 * Math.Pow((y0 - yy), 2))) / (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)));
                                double FyP2 = -(((ri * ri)) + (2 * (ypi * ypi)));
                                double[,] AA = new double[2, 6];
                                double[,] BB = new double[2, 3];
                                double[,] c = new double[2, 10];


                                // MessageBox.Show("" + (Math.Pow(((Math.Pow((x0 - xx), 2)) + (Math.Pow((y0 - yy), 2))), (1.0 / 2)))+"   "+ Math.Pow((x0 - xx), 2));

                                double[,] B = new double[2, 2];
                                // MessageBox.Show("" + xx0 + " " + xy0 + " " + xK1 + " " + xK2 + " " + xK3 + " " + xP1 + " " + xP2 + " " + xB1 + " " + xB2);
                                //  MessageBox.Show("" + yx0 + " " + yy0 + " " + yK1 + " " + yK2 + " " + yK3 + " " + yP1 + " " + yP2 );
                                //    System.IO.StreamWriter l2 = new System.IO.StreamWriter(@"E:\test\B.txt");
                                //   MessageBox.Show("" + dx[0, 0]);


                              //  B[0, 0] = Fx_xx; B[0, 1] = Fx_yy;
                               // B[1, 0] = Fy_xx; B[1, 1] = Fy_yy;
                                double b00 = Fx_xx; double b01 = Fx_yy; double b10 = Fy_xx; double b11 = Fy_yy;
                                // for (int i3 = 0; i3 < 2; i3++)
                                //  {
                                //    string line2 = "" + B[i3, 0];
                                //      for (int i4 = 1; i4 < 2; i4++)
                                //    {
                                //        line2 = line2 + " " + B[i3, i4];

                                //   }
                                //l2.WriteLine(line2);
                                // }
                                //  l2.Close();
                                //   double[,] PP = new double[2, 2];
                                //  PP[0, 0] = 1.0 / pp2[0, 0]; PP[1, 1] = 1.0 / pp2[1, 1];
                                double p00 = 1.0 / pp2[0, 0]; double p11 = 1.0 / pp2[1, 1];
                                double detbpb = (b00 * b11 * p00 * p11 * (b00) * (b11) - b00 * b11 * p00 * p11 * (b01) * (b10) - b01 * b10 * p00 * p11 * (b00) * (b11) + b01 * b10 * p00 * p11 * (b01) * (b10));
                                //  double[,] PB = Matrix.Inverse(Matrix.Multiply(Matrix.Multiply(B, PP), Matrix.Transpose(B)));
                                double[,] PB = new double[2, 2];
                                PB[0, 0] = (b10 * p00 * (b10) + b11 * p11 * (b11)) / detbpb;
                                double pb00 = (b10 * p00 * (b10) + b11 * p11 * (b11)) / detbpb;
                                PB[0, 1] = -(b00 * p00 * (b10) + b01 * p11 * (b11)) / detbpb;
                                double pb01 = -(b00 * p00 * (b10) + b01 * p11 * (b11)) / detbpb;
                                PB[1, 0] = -(b10 * p00 * (b00) + b11 * p11 * (b01)) / detbpb;
                                double pb10 = -(b10 * p00 * (b00) + b11 * p11 * (b01)) / detbpb;
                                PB[1, 1] = (b00 * p00 * (b00) + b01 * p11 * (b01)) / detbpb;
                                double pb11 = (b00 * p00 * (b00) + b01 * p11 * (b01)) / detbpb;

                                AA[0, 0] = xom;
                                double a00 = xom;
                                AA[0, 1] = xphi;
                                double a01 = xphi;
                                AA[0, 2] = xk;
                                double a02 = xk;
                                AA[0, 3] = xXL;
                                double a03 = xXL;
                                AA[0, 4] = xYL;
                                double a04 = xYL;
                                AA[0, 5] = xZL;
                                double a05 = xZL;
                                AA[1, 0] = yom;
                                double a10 = yom;
                                AA[1, 1] = yphi;
                                double a11 = yphi;
                                AA[1, 2] = yk;
                                double a12 = yk;
                                AA[1, 3] = yXL;
                                double a13 = yXL;
                                AA[1, 4] = yYL;
                                double a14 = yYL;
                                AA[1, 5] = yZL;
                                double a15 = yZL;
                                BB[0, 0] = xX;
                                double bb00 = xX;
                                BB[0, 1] = xY;
                                double bb01 = xY;
                                BB[0, 2] = xZ;
                                double bb02 = xZ;
                                BB[1, 0] = yX;
                                double bb10 = yX;
                                BB[1, 1] = yY;
                                double bb11 = yY;
                                BB[1, 2] = yZ;
                                double bb12 = yZ;
                                c[0, 0] = Fxx0;
                                double c00 = Fxx0;
                                c[0, 1] = Fxy0;
                                double c01 = Fxy0;
                                c[0, 2] = FxK1;
                                double c02 = FxK1;
                                c[0, 3] = FxK2;
                                double c03 = FxK2;
                                c[0, 4] = FxK3;
                                double c04 = FxK3;
                                c[0, 5] = FxP1;
                                double c05 = FxP1;
                                c[0, 6] = FxP2;
                                double c06 = FxP2;
                                c[0, 7] = FxB1;
                                double c07 = FxB1;
                                c[0, 8] = FxB2;
                                double c08 = FxB2;
                                c[0, 9] = Fxf;
                                double c09 = Fxf;
                                c[1, 0] = Fyx0;
                                double c10 = Fyx0;
                                c[1, 1] = Fyy0;
                                double c11 = Fyy0;
                                c[1, 2] = FyK1;
                                double c12 = FyK1;
                                c[1, 3] = FyK2;
                                double c13 = FyK2;
                                c[1, 4] = FyK3;
                                double c14 = FyK3;
                                c[1, 5] = FyP1;
                                double c15 = FyP1;
                                c[1, 6] = FyP2;
                                double c16 = FyP2;
                                c[1, 7] = 0;
                                double c17 = 0;
                                c[1, 8] = 0;
                                double c18 = 0;
                                c[1, 9] = Fyf;
                                double c19 = Fyf;
                                lll[0, 0] = Fx; lll[1, 0] = Fy;
                                double ll0 = Fx; double ll1 = Fy;
                              
                                lllnb[0, 0] = observe_value_points[i, 0] - X; lllnb[1, 0] = observe_value_points[i, 1] - Y; lllnb[2, 0] = observe_value_points[i, 2] - Z;
                                // MessageBox.Show("" + lll[0, 0]+"   " + lll[1, 0]);
                                //    double[,] Fc2 = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(c), PB), BB);
                                double[,] Fc2 = new double[10, 3];
                                Fc2[0, 0] = bb00 * (pb00 * (c00) + pb10 * (c10)) + bb10 * (pb01 * (c00) + pb11 * (c10)); Fc2[0, 1] = bb01 * (pb00 * (c00) + pb10 * (c10)) + bb11 * (pb01 * (c00) + pb11 * (c10)); Fc2[0, 2] = bb02 * (pb00 * (c00) + pb10 * (c10)) + bb12 * (pb01 * (c00) + pb11 * (c10));
                                Fc2[1, 0] = bb00 * (pb00 * (c01) + pb10 * (c11)) + bb10 * (pb01 * (c01) + pb11 * (c11)); Fc2[1, 1] = bb01 * (pb00 * (c01) + pb10 * (c11)) + bb11 * (pb01 * (c01) + pb11 * (c11)); Fc2[1, 2] = bb02 * (pb00 * (c01) + pb10 * (c11)) + bb12 * (pb01 * (c01) + pb11 * (c11));
                                Fc2[2, 0] = bb00 * (pb00 * (c02) + pb10 * (c12)) + bb10 * (pb01 * (c02) + pb11 * (c12)); Fc2[2, 1] = bb01 * (pb00 * (c02) + pb10 * (c12)) + bb11 * (pb01 * (c02) + pb11 * (c12)); Fc2[2, 2] = bb02 * (pb00 * (c02) + pb10 * (c12)) + bb12 * (pb01 * (c02) + pb11 * (c12));
                                Fc2[3, 0] = bb00 * (pb00 * (c03) + pb10 * (c13)) + bb10 * (pb01 * (c03) + pb11 * (c13)); Fc2[3, 1] = bb01 * (pb00 * (c03) + pb10 * (c13)) + bb11 * (pb01 * (c03) + pb11 * (c13)); Fc2[3, 2] = bb02 * (pb00 * (c03) + pb10 * (c13)) + bb12 * (pb01 * (c03) + pb11 * (c13));
                                Fc2[4, 0] = bb00 * (pb00 * (c04) + pb10 * (c14)) + bb10 * (pb01 * (c04) + pb11 * (c14)); Fc2[4, 1] = bb01 * (pb00 * (c04) + pb10 * (c14)) + bb11 * (pb01 * (c04) + pb11 * (c14)); Fc2[4, 2] = bb02 * (pb00 * (c04) + pb10 * (c14)) + bb12 * (pb01 * (c04) + pb11 * (c14));
                                Fc2[5, 0] = bb00 * (pb00 * (c05) + pb10 * (c15)) + bb10 * (pb01 * (c05) + pb11 * (c15)); Fc2[5, 1] = bb01 * (pb00 * (c05) + pb10 * (c15)) + bb11 * (pb01 * (c05) + pb11 * (c15)); Fc2[5, 2] = bb02 * (pb00 * (c05) + pb10 * (c15)) + bb12 * (pb01 * (c05) + pb11 * (c15));
                                Fc2[6, 0] = bb00 * (pb00 * (c06) + pb10 * (c16)) + bb10 * (pb01 * (c06) + pb11 * (c16)); Fc2[6, 1] = bb01 * (pb00 * (c06) + pb10 * (c16)) + bb11 * (pb01 * (c06) + pb11 * (c16)); Fc2[6, 2] = bb02 * (pb00 * (c06) + pb10 * (c16)) + bb12 * (pb01 * (c06) + pb11 * (c16));
                                Fc2[7, 0] = bb00 * (pb00 * (c07) + pb10 * (c17)) + bb10 * (pb01 * (c07) + pb11 * (c17)); Fc2[7, 1] = bb01 * (pb00 * (c07) + pb10 * (c17)) + bb11 * (pb01 * (c07) + pb11 * (c17)); Fc2[7, 2] = bb02 * (pb00 * (c07) + pb10 * (c17)) + bb12 * (pb01 * (c07) + pb11 * (c17));
                                Fc2[8, 0] = bb00 * (pb00 * (c08) + pb10 * (c18)) + bb10 * (pb01 * (c08) + pb11 * (c18)); Fc2[8, 1] = bb01 * (pb00 * (c08) + pb10 * (c18)) + bb11 * (pb01 * (c08) + pb11 * (c18)); Fc2[8, 2] = bb02 * (pb00 * (c08) + pb10 * (c18)) + bb12 * (pb01 * (c08) + pb11 * (c18));
                                Fc2[9, 0] = bb00 * (pb00 * (c09) + pb10 * (c19)) + bb10 * (pb01 * (c09) + pb11 * (c19)); Fc2[9, 1] = bb01 * (pb00 * (c09) + pb10 * (c19)) + bb11 * (pb01 * (c09) + pb11 * (c19)); Fc2[9, 2] = bb02 * (pb00 * (c09) + pb10 * (c19)) + bb12 * (pb01 * (c09) + pb11 * (c19));

                                //   MessageBox.Show("" + Fc[0, 0] + "    " + i + "    " + j2);
                                for (int pp = 0; pp < 10; pp++)
                                {
                                    for (int ppp = 0; ppp < 3; ppp++)
                                    {
                                        Fc[pp, ppp] = Fc[pp, ppp] + Fc2[pp, ppp];
                                    }
                                }
                                //  double[,] ebii2 = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(BB), PB), lll);
                                double[,] ebii2 = new double[3, 1];
                                ebii2[0, 0] = ll0 * (pb00 * (bb00) + pb10 * (bb10)) + ll1 * (pb01 * (bb00) + pb11 * (bb10));
                                ebii2[1, 0] = ll0 * (pb00 * (bb01) + pb10 * (bb11)) + ll1 * (pb01 * (bb01) + pb11 * (bb11));
                                ebii2[2, 0] = ll0 * (pb00 * (bb02) + pb10 * (bb12)) + ll1 * (pb01 * (bb02) + pb11 * (bb12));

                                for (int pp = 0; pp < 3; pp++)
                                {
                                    ebii[pp, 0] = ebii[pp, 0] + ebii2[pp, 0];
                                }
                                // v2 = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(BB), PB), BB);
                                v2[0, 0] = bb00 * (pb00 * (bb00) + pb10 * (bb10)) + bb10 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 1] = bb01 * (pb00 * (bb00) + pb10 * (bb10)) + bb11 * (pb01 * (bb00) + pb11 * (bb10)); v2[0, 2] = bb02 * (pb00 * (bb00) + pb10 * (bb10)) + bb12 * (pb01 * (bb00) + pb11 * (bb10));
                                v2[1, 0] = bb00 * (pb00 * (bb01) + pb10 * (bb11)) + bb10 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 1] = bb01 * (pb00 * (bb01) + pb10 * (bb11)) + bb11 * (pb01 * (bb01) + pb11 * (bb11)); v2[1, 2] = bb02 * (pb00 * (bb01) + pb10 * (bb11)) + bb12 * (pb01 * (bb01) + pb11 * (bb11));
                                v2[2, 0] = bb00 * (pb00 * (bb02) + pb10 * (bb12)) + bb10 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 1] = bb01 * (pb00 * (bb02) + pb10 * (bb12)) + bb11 * (pb01 * (bb02) + pb11 * (bb12)); v2[2, 2] = bb02 * (pb00 * (bb02) + pb10 * (bb12)) + bb12 * (pb01 * (bb02) + pb11 * (bb12));

                                // ww = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(AA), PB), BB);
                                ww[0, 0] = bb00 * (pb00 * (a00) + pb10 * (a10)) + bb10 * (pb01 * (a00) + pb11 * (a10)); ww[0, 1] = bb01 * (pb00 * (a00) + pb10 * (a10)) + bb11 * (pb01 * (a00) + pb11 * (a10)); ww[0, 2] = bb02 * (pb00 * (a00) + pb10 * (a10)) + bb12 * (pb01 * (a00) + pb11 * (a10));
                                ww[1, 0] = bb00 * (pb00 * (a01) + pb10 * (a11)) + bb10 * (pb01 * (a01) + pb11 * (a11)); ww[1, 1] = bb01 * (pb00 * (a01) + pb10 * (a11)) + bb11 * (pb01 * (a01) + pb11 * (a11)); ww[1, 2] = bb02 * (pb00 * (a01) + pb10 * (a11)) + bb12 * (pb01 * (a01) + pb11 * (a11));
                                ww[2, 0] = bb00 * (pb00 * (a02) + pb10 * (a12)) + bb10 * (pb01 * (a02) + pb11 * (a12)); ww[2, 1] = bb01 * (pb00 * (a02) + pb10 * (a12)) + bb11 * (pb01 * (a02) + pb11 * (a12)); ww[2, 2] = bb02 * (pb00 * (a02) + pb10 * (a12)) + bb12 * (pb01 * (a02) + pb11 * (a12));
                                ww[3, 0] = bb00 * (pb00 * (a03) + pb10 * (a13)) + bb10 * (pb01 * (a03) + pb11 * (a13)); ww[3, 1] = bb01 * (pb00 * (a03) + pb10 * (a13)) + bb11 * (pb01 * (a03) + pb11 * (a13)); ww[3, 2] = bb02 * (pb00 * (a03) + pb10 * (a13)) + bb12 * (pb01 * (a03) + pb11 * (a13));
                                ww[4, 0] = bb00 * (pb00 * (a04) + pb10 * (a14)) + bb10 * (pb01 * (a04) + pb11 * (a14)); ww[4, 1] = bb01 * (pb00 * (a04) + pb10 * (a14)) + bb11 * (pb01 * (a04) + pb11 * (a14)); ww[4, 2] = bb02 * (pb00 * (a04) + pb10 * (a14)) + bb12 * (pb01 * (a04) + pb11 * (a14));
                                ww[5, 0] = bb00 * (pb00 * (a05) + pb10 * (a15)) + bb10 * (pb01 * (a05) + pb11 * (a15)); ww[5, 1] = bb01 * (pb00 * (a05) + pb10 * (a15)) + bb11 * (pb01 * (a05) + pb11 * (a15)); ww[5, 2] = bb02 * (pb00 * (a05) + pb10 * (a15)) + bb12 * (pb01 * (a05) + pb11 * (a15));
                                double ww00 = ww[0, 0]; double ww01 = ww[0, 1]; double ww02 = ww[0, 2];
                                double ww10 = ww[1, 0]; double ww11 = ww[1, 1]; double ww12 = ww[1, 2];
                                double ww20 = ww[2, 0]; double ww21 = ww[2, 1]; double ww22 = ww[2, 2];
                                double ww30 = ww[3, 0]; double ww31 = ww[3, 1]; double ww32 = ww[3, 2];
                                double ww40 = ww[4, 0]; double ww41 = ww[4, 1]; double ww42 = ww[4, 2];
                                double ww50 = ww[5, 0]; double ww51 = ww[5, 1]; double ww52 = ww[5, 2];
                                double[,] dxx2 = new double[6, 1];
                                
                                for (int pp = 0; pp < 6; pp++)
                                {
                                    dxx2[0, 0] = dx[6 * j2, 0];
                                    dxx2[1, 0] = dx[6 * j2 + 1, 0];
                                    dxx2[2, 0] = dx[6 * j2 + 2, 0];
                                    dxx2[3, 0] = dx[6 * j2 + 3, 0];
                                    dxx2[4, 0] = dx[6 * j2 + 4, 0];
                                    dxx2[5, 0] = dx[6 * j2 + 5, 0];
                                }
                                double dxx200 = dxx2[0, 0]; double dxx210 = dxx2[1, 0]; double dxx220 = dxx2[2, 0];
                                double dxx230 = dxx2[3, 0]; double dxx240 = dxx2[4, 0]; double dxx250 = dxx2[5, 0];
                               //  double[,] dxw = Matrix.Multiply(Matrix.Transpose(ww), dxx2);
                                double[,] dxw = new double[3, 1];
                                dxw[0, 0] = dxx200 * (ww00) + dxx210 * (ww10) + dxx220 * (ww20) + dxx230 * (ww30) + dxx240 * (ww40) + dxx250 * (ww50);
                                dxw[1, 0] = dxx200 * (ww01) + dxx210 * (ww11) + dxx220 * (ww21) + dxx230 * (ww31) + dxx240 * (ww41) + dxx250 * (ww51);
                                dxw[2, 0] = dxx200 * (ww02) + dxx210 * (ww12) + dxx220 * (ww22) + dxx230 * (ww32) + dxx240 * (ww42) + dxx250 * (ww52);


                                for (int ppp = 0; ppp < 3; ppp++)
                                {
                                    ebwdx[ppp, 0] = ebwdx[ppp, 0] - dxw[ppp, 0];
                                }

                                w1.Add(ww[0, 0]);
                                w1.Add(ww[1, 0]);
                                w1.Add(ww[2, 0]);
                                w1.Add(ww[3, 0]);
                                w1.Add(ww[4, 0]);
                                w1.Add(ww[5, 0]);

                                w2.Add(ww[0, 1]);
                                w2.Add(ww[1, 1]);
                                w2.Add(ww[2, 1]);
                                w2.Add(ww[3, 1]);
                                w2.Add(ww[4, 1]);
                                w2.Add(ww[5, 1]);

                                w3.Add(ww[0, 2]);
                                w3.Add(ww[1, 2]);
                                w3.Add(ww[2, 2]);
                                w3.Add(ww[3, 2]);
                                w3.Add(ww[4, 2]);
                                w3.Add(ww[5, 2]);

                                wj.Add(j2);
                                //WW[6 * j2 + 0, 0] = ww[0, 0]; WW[6 * j2 + 1, 0] = ww[1, 0]; WW[6 * j2 + 2, 0] = ww[2, 0]; WW[6 * j2 + 3, 0] = ww[3, 0]; WW[6 * j2 + 4, 0] = ww[4, 0]; WW[6 * j2 + 5, 0] = ww[5, 0];
                                //WW[6 * j2 + 0, 1] = ww[0, 1]; WW[6 * j2 + 1, 1] = ww[1,1]; WW[6 * j2 + 2, 1] = ww[2, 1]; WW[6 * j2 + 3, 1] = ww[3, 1]; WW[6 * j2 + 4, 1] = ww[4, 1]; WW[6 * j2 + 5, 1] = ww[5, 1];
                                //WW[6 * j2 + 0, 2] = ww[0, 2]; WW[6 * j2 + 1, 2] = ww[1, 2]; WW[6 * j2 + 2, 2] = ww[2, 2]; WW[6 * j2 + 3, 2] = ww[3,2]; WW[6 * j2 + 4,2] = ww[4, 2]; WW[6 * j2 + 5, 2] = ww[5, 2];
                                v[0, 0] = v[0, 0] + v2[0, 0]; v[1, 0] = v[1, 0] + v2[1, 0]; v[2, 0] = v[2, 0] + v2[2, 0];
                                v[0, 1] = v[0, 1] + v2[0, 1]; v[1, 1] = v[1, 1] + v2[1, 1]; v[2, 1] = v[2, 1] + v2[2, 1];
                                v[0, 2] = v[0, 2] + v2[0, 2]; v[1, 2] = v[1, 2] + v2[1, 2]; v[2, 2] = v[2, 2] + v2[2, 2];

                        }
                        double Fc00, Fc01, Fc02, Fc10, Fc11, Fc12, Fc20, Fc21, Fc22, Fc30, Fc31, Fc32, Fc40, Fc41, Fc42, Fc50, Fc51, Fc52, Fc60, Fc61, Fc62, Fc70, Fc71, Fc72, Fc80, Fc81, Fc82, Fc90, Fc91, Fc92;

                        Fc00 = Fc[0, 0]; Fc01 = Fc[0, 1]; Fc02 = Fc[0, 2]; Fc10 = Fc[1, 0]; Fc11 = Fc[1, 1]; Fc12 = Fc[1, 2];
                        Fc20 = Fc[2, 0]; Fc21 = Fc[2, 1]; Fc22 = Fc[2, 2]; Fc30 = Fc[3, 0]; Fc31 = Fc[3, 1]; Fc32 = Fc[3, 2];
                        Fc40 = Fc[4, 0]; Fc41 = Fc[4, 1]; Fc42 = Fc[4, 2]; Fc50 = Fc[5, 0]; Fc51 = Fc[5, 1]; Fc52 = Fc[5, 2];
                        Fc60 = Fc[6, 0]; Fc61 = Fc[6, 1]; Fc62 = Fc[6, 2]; Fc70 = Fc[7, 0]; Fc71 = Fc[7, 1]; Fc72 = Fc[7, 2];
                        Fc80 = Fc[8, 0]; Fc81 = Fc[8, 1]; Fc82 = Fc[8, 2]; Fc90 = Fc[9, 0]; Fc91 = Fc[9, 1]; Fc92 = Fc[9, 2];

                        double[,] dxx22 = new double[10, 1];
                        for (int pp = 0; pp < 10; pp++)
                        {

                            dxx22[pp, 0] = dx[6 * image_name_Eterior.Length + pp, 0];

                        }
                        double dxx2200 = dxx22[0, 0]; double dxx2210 = dxx22[1, 0]; double dxx2220 = dxx22[2, 0]; double dxx2230 = dxx22[3, 0]; double dxx2240 = dxx22[4, 0];
                        double dxx2250 = dxx22[5, 0]; double dxx2260 = dxx22[6, 0]; double dxx2270 = dxx22[7, 0]; double dxx2280 = dxx22[8, 0]; double dxx2290 = dxx22[9, 0];
                       // double[,] Fdx = Matrix.Multiply(Matrix.Transpose(Fc), dxx22);
                        double[,] Fdx = new double[3, 1];
                        Fdx[0, 0] = dxx2200 * (Fc00) + dxx2210 * (Fc10) + dxx2220 * (Fc20) + dxx2230 * (Fc30) + dxx2240 * (Fc40) + dxx2250 * (Fc50) + dxx2260 * (Fc60) + dxx2270 * (Fc70) + dxx2280 * (Fc80) + dxx2290 * (Fc90);
                        Fdx[1, 0] = dxx2200 * (Fc01) + dxx2210 * (Fc11) + dxx2220 * (Fc21) + dxx2230 * (Fc31) + dxx2240 * (Fc41) + dxx2250 * (Fc51) + dxx2260 * (Fc61) + dxx2270 * (Fc71) + dxx2280 * (Fc81) + dxx2290 * (Fc91);
                        Fdx[2, 0] = dxx2200 * (Fc02) + dxx2210 * (Fc12) + dxx2220 * (Fc22) + dxx2230 * (Fc32) + dxx2240 * (Fc42) + dxx2250 * (Fc52) + dxx2260 * (Fc62) + dxx2270 * (Fc72) + dxx2280 * (Fc82) + dxx2290 * (Fc92);


                        double[,] pp2nB = new double[3, 3];
                      
                        pp2nB[0, 0] = observe_value_points_weight[i, 0];
                        pp2nB[1, 1] = observe_value_points_weight[i, 1];
                        pp2nB[2, 2] = observe_value_points_weight[i, 2];
                       
                        v[0, 0] = v[0, 0] + pp2nB[0, 0]; v[1, 1] = v[1, 1] + pp2nB[1, 1]; v[2, 2] = v[2, 2] + pp2nB[2, 2];

                        ebii[0, 0] = ebii[0, 0] - (pp2nB[0, 0] * lllnb[0, 0]);
                        ebii[1, 0] = ebii[1, 0] - (pp2nB[1, 1] * lllnb[1, 0]);
                        ebii[2, 0] = ebii[2, 0] - (pp2nB[2, 2] * lllnb[2, 0]);
                        for (int ppp = 0; ppp < 3; ppp++)
                        {
                            ebwdx[ppp, 0] = ebwdx[ppp, 0] + ebii[ppp, 0] - Fdx[ppp,0];
                        }
                        double ebwdx00 = ebwdx[0, 0]; double ebwdx10 = ebwdx[1, 0]; double ebwdx20 = ebwdx[2, 0];

                        double[,] v3 = new double[3, 3];
                        double detv = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                        v3[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                        v3[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                        v3[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                        v3[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                        v3[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                        v3[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                        v3[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                        v3[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                        v3[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                        double v300 = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                        double v301 = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                        double v302 = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                        double v310 = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                        double v311 = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                        double v312 = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                        double v320 = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                        double v321 = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                        double v322 = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                       // double[,] dxb = Matrix.Multiply(v3, ebwdx);
                        double[,] dxb = new double[3, 1];
                        dxb[0, 0] = ebwdx00 * v300 + ebwdx10 * v301 + ebwdx20 * v302;
                        dxb[1, 0] = ebwdx00 * v310 + ebwdx10 * v311 + ebwdx20 * v312;
                        dxb[2, 0] = ebwdx00 * v320 + ebwdx10 * v321 + ebwdx20 * v322;
                       // s1.Stop();
                       // MessageBox.Show("sina" + s1.ElapsedMilliseconds);
                        point_in = point_in + 1;
                        Initial_value[i, 0] = Initial_value[i, 0] - dxb[0, 0];
                        Initial_value[i, 1] = Initial_value[i, 1] - dxb[1, 0];
                        Initial_value[i, 2] = Initial_value[i, 2] - dxb[2, 0];
                      
                    }
                }
                for (int i = 0; i < image_name_Eterior.Length; i++)
                {

                    Exterior_Orientation_Eterior[i, 0] = Exterior_Orientation_Eterior[i, 0] - dx[6 * i + 3, 0];
                    Exterior_Orientation_Eterior[i, 1] = Exterior_Orientation_Eterior[i, 1] - dx[6 * i + 4, 0];
                    Exterior_Orientation_Eterior[i, 2] = Exterior_Orientation_Eterior[i, 2] - dx[6 * i + 5, 0];
                    Exterior_Orientation_Eterior[i, 3] = Exterior_Orientation_Eterior[i, 3] - (dx[6 * i, 0]);
                    Exterior_Orientation_Eterior[i, 4] = Exterior_Orientation_Eterior[i, 4] - (dx[6 * i + 1, 0]);
                    Exterior_Orientation_Eterior[i, 5] = Exterior_Orientation_Eterior[i, 5] - (dx[6 * i + 2, 0]);
                    //   MessageBox.Show("" + Exterior_Orientation_Eterior[i, 0]+"    "+ Exterior_Orientation_Eterior[i, 1] + "    " + Exterior_Orientation_Eterior[i, 2] + "    " + Exterior_Orientation_Eterior[i, 3] + "    " + Exterior_Orientation_Eterior[i, 4] + "    " + Exterior_Orientation_Eterior[i, 5]);
                }

                IOP[0, 0] = IOP[0, 0] - dx[6 * image_name_Eterior.Length + 9, 0];
                IOP[1, 0] = IOP[1, 0] - dx[6 * image_name_Eterior.Length , 0];
                IOP[2, 0] = IOP[2, 0] - dx[6 * image_name_Eterior.Length + 1, 0];
                IOP[3, 0] = IOP[3, 0] - dx[6 * image_name_Eterior.Length + 2, 0];
                IOP[4, 0] = IOP[4, 0] - dx[6 * image_name_Eterior.Length + 3, 0];
                IOP[5, 0] = IOP[5, 0] - dx[6 * image_name_Eterior.Length + 4, 0];
                IOP[6, 0] = IOP[6, 0] - dx[6 * image_name_Eterior.Length + 5, 0];
                IOP[7, 0] = IOP[7, 0] - dx[6 * image_name_Eterior.Length + 6, 0];
                IOP[8, 0] = IOP[8, 0] - dx[6 * image_name_Eterior.Length + 7, 0];
                IOP[9, 0] = IOP[9, 0] - dx[6 * image_name_Eterior.Length + 8, 0];
                //  if (iteration == 50 || iteration == 60 || iteration == 700 || iteration == 99)
                //   {
               //  MessageBox.Show("" + dx[0, 0]);
                //   }
                // camera_weight[6 * image_name_Eterior.Length] = invS[6 * image_name_Eterior.Length, 6 * image_name_Eterior.Length];
                // camera_weight[6 * image_name_Eterior.Length + 1] = invS[6 * image_name_Eterior.Length + 1, 6 * image_name_Eterior.Length + 1]; camera_weight[6 * image_name_Eterior.Length + 2] = invS[6 * image_name_Eterior.Length + 2, 6 * image_name_Eterior.Length + 2]; camera_weight[6 * image_name_Eterior.Length + 3] = invS[6 * image_name_Eterior.Length + 3, 6 * image_name_Eterior.Length + 3]; camera_weight[6 * image_name_Eterior.Length + 4] = invS[6 * image_name_Eterior.Length + 4, 6 * image_name_Eterior.Length + 4];
                // camera_weight[6 * image_name_Eterior.Length + 5] = invS[6 * image_name_Eterior.Length + 5, 6 * image_name_Eterior.Length + 5]; camera_weight[6 * image_name_Eterior.Length + 6] = invS[6 * image_name_Eterior.Length + 6, 6 * image_name_Eterior.Length + 6]; camera_weight[6 * image_name_Eterior.Length + 7] = invS[6 * image_name_Eterior.Length + 7, 6 * image_name_Eterior.Length + 7]; camera_weight[6 * image_name_Eterior.Length + 8] = invS[6 * image_name_Eterior.Length + 8, 6 * image_name_Eterior.Length + 8]; camera_weight[6 * image_name_Eterior.Length + 9] = invS[6 * image_name_Eterior.Length + 9, 6 * image_name_Eterior.Length + 9];

            //    s2.Stop();
              //  MessageBox.Show("sina" + s2.ElapsedMilliseconds);
            }
           

            point_in = point_in - 1;
            tie_points_coordinate = Initial_value; Exterior_orientation = Exterior_Orientation_Eterior; cov_Exterior_orientation = cov_Ex;
        }
    
        public void v_calculator(int i, int j2, List<double> camera_weight, double[,] Exterior_Orientation_Eterior, double[] interior, double[,] tie_com, double[,] image_m, double[,] image_sin_cos, double[,] im_ob_sig,out double[,] v2)
        {
            double sx = 1.0 / camera_weight[8 * j2];
            double sy = 1.0 / camera_weight[8 * j2 + 1];
            double XL = Exterior_Orientation_Eterior[j2, 0];
            double YL = Exterior_Orientation_Eterior[j2, 1];
            double ZL = Exterior_Orientation_Eterior[j2, 2];
            double om = Exterior_Orientation_Eterior[j2, 3];
            double phi = Exterior_Orientation_Eterior[j2, 4];
            double k = Exterior_Orientation_Eterior[j2, 5];
            double ff = interior[0];

            double X = tie_com[i, 0];
            double Y = tie_com[i, 1];
            double Z = tie_com[i, 2];
            double m11 = image_m[j2, 0];
            double m12 = image_m[j2, 1];
            //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
            double m13 = image_m[j2, 2];
            //m21=-cos(phi)*sin(k)
            double m21 = image_m[j2, 3];
            //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
            double m22 = image_m[j2, 4];
            //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
            double m23 = image_m[j2, 5];
            //m31=sin(phi)
            double m31 = image_m[j2, 6];
            //m32=-sin(om)*cos(phi)
            double m32 = image_m[j2, 7];
            //m33=cos(om)*cos(phi)
            double m33 = image_m[j2, 8];

            double sinom = image_sin_cos[j2, 0]; double sinphi = image_sin_cos[j2, 1]; double sink = image_sin_cos[j2, 2];
            double cosom = image_sin_cos[j2, 3]; double cosphi = image_sin_cos[j2, 4]; double cosk = image_sin_cos[j2, 5];





            double sinphiXXL = (sinphi * (X - XL) + cosom * cosphi * (Z - ZL) - cosphi * sinom * (Y - YL));
            double sinksinom = (sink * sinom - cosk * cosom * sinphi);
            double cosomsink = (cosom * sink + cosk * sinom * sinphi);
            double coskcosom = (cosk * cosom - sink * sinom * sinphi);
            double cosksinom = (cosk * sinom + cosom * sink * sinphi);
            double xom = -((ff * (sinksinom * (Y - YL) - cosomsink * (Z - ZL))) / sinphiXXL - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));

            double xphi = -((ff * (cosk * sinphi * (X - XL) + cosk * cosom * cosphi * (Z - ZL) - cosk * cosphi * sinom * (Y - YL))) / sinphiXXL + (ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
            double xk = -(-(ff * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / sinphiXXL);
            double xXL = -((ff * cosk * cosphi) / sinphiXXL - (ff * sinphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
            double xYL = -((ff * cosomsink) / sinphiXXL + (ff * cosphi * sinom * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
            double xZL = -(-((ff * cosom * cosphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * sinksinom) / sinphiXXL));
            //   (f * cos(om) * cos(phi) * ((cos(om) * sin(k) + cos(k) * sin(om) * sin(phi)) * (Y - YL) + (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi)) * (Z - ZL) + cos(k) * cos(phi) * (X - XL))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL)) ^ 2 - (f * (sin(k) * sin(om) - cos(k) * cos(om) * sin(phi))) / (sin(phi) * (X - XL) + cos(om) * cos(phi) * (Z - ZL) - cos(phi) * sin(om) * (Y - YL))];
            double yom = -((ff * (cosksinom * (Y - YL) - coskcosom * (Z - ZL))) / sinphiXXL - (ff * (cosom * cosphi * (Y - YL) + cosphi * sinom * (Z - ZL)) * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
            double yphi = -((ff * (cosphi * (X - XL) - cosom * sinphi * (Z - ZL) + sinom * sinphi * (Y - YL)) * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * (sink * sinphi * (X - XL) + cosom * cosphi * sink * (Z - ZL) - cosphi * sink * sinom * (Y - YL))) / sinphiXXL);
            double yk = -((ff * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / sinphiXXL);
            double yXL = -(-(ff * sinphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosphi * sink) / sinphiXXL);
            double yYL = -((ff * coskcosom) / sinphiXXL + (ff * cosphi * sinom * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
            double yZL = -((ff * cosksinom) / sinphiXXL - (ff * cosom * cosphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
            double xX = ((ff * sinphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosk * cosphi) / sinphiXXL);
            double xY = (-(ff * cosomsink) / sinphiXXL - (ff * cosphi * sinom * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL));
            double xZ = ((ff * cosom * cosphi * (cosomsink * (Y - YL) + sinksinom * (Z - ZL) + cosk * cosphi * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * sinksinom) / sinphiXXL);
            double yX = ((ff * sinphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) + (ff * cosphi * sink) / sinphiXXL);
            double yY = (-(ff * coskcosom) / sinphiXXL - (ff * cosphi * sinom * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL));
            double yZ = ((ff * cosom * cosphi * (coskcosom * (Y - YL) + cosksinom * (Z - ZL) - cosphi * sink * (X - XL))) / (sinphiXXL * sinphiXXL) - (ff * cosksinom) / sinphiXXL);



            double[,] BB = new double[2, 3];
            double[,] B = new double[2, 8];
            B[0, 0] = 1; B[0, 1] = 0; B[0, 2] = xom; B[0, 3] = xphi; B[0, 4] = xk; B[0, 5] = xXL; B[0, 6] = xYL; B[0, 7] = xZL;
            B[1, 0] = 0; B[1, 1] = 1; B[1, 2] = yom; B[1, 3] = yphi; B[1, 4] = yk; B[1, 5] = yXL; B[1, 6] = yYL; B[1, 7] = yZL;
            double som = im_ob_sig[6 * j2 + 0, 6 * j2 + 0]; double sphi = im_ob_sig[6 * j2 + 1, 6 * j2 + 1]; double sk = im_ob_sig[6 * j2 + 2, 6 * j2 + 2];
            double sxl = im_ob_sig[6 * j2 + 3, 6 * j2 + 3]; double syl = im_ob_sig[6 * j2 + 4, 6 * j2 + 4]; double szl = im_ob_sig[6 * j2 + 5, 6 * j2 + 5];
            double[,] pp = new double[2, 2];
            pp[0, 0] = sx + sxl * xXL * (xXL) + syl * xYL * (xYL) + szl * xZL * (xZL) + sk * xk * (xk) + som * xom * (xom) + sphi * xphi * (xphi);
            pp[0, 1] = sxl * xXL * (yXL) + syl * xYL * (yYL) + szl * xZL * (yZL) + sk * xk * (yk) + som * xom * (yom) + sphi * xphi * (yphi);
            pp[1, 0] = sxl * yXL * (xXL) + syl * yYL * (xYL) + szl * yZL * (xZL) + sk * yk * (xk) + som * yom * (xom) + sphi * yphi * (xphi);
            pp[1, 1] = sy + sxl * yXL * (yXL) + syl * yYL * (yYL) + szl * yZL * (yZL) + sk * yk * (yk) + som * yom * (yom) + sphi * yphi * (yphi);
            double[,] p = new double[2, 2];
            double detp = (pp[0, 0] * pp[1, 1]) - (pp[0, 1] * pp[1, 0]);
            p[0, 0] = pp[1, 1] / detp;
            p[0, 1] = -pp[0, 1] / detp;
            p[1, 0] = -pp[1, 0] / detp;
            p[1, 1] = pp[0, 0] / detp;
            BB[0, 0] = xX;
            BB[0, 1] = xY;
            BB[0, 2] = xZ;
            BB[1, 0] = yX;
            BB[1, 1] = yY;
            BB[1, 2] = yZ;
            v2 = new double[3, 3];
            // double [,] v211 = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(BB), pp2), BB);
            v2[0, 0] = BB[0, 0] * p[0, 0] * BB[0, 0] + BB[1, 0] * p[1, 1] * BB[1, 0];
            v2[0, 1] = BB[0, 1] * p[0, 0] * BB[0, 0] + BB[1, 1] * p[1, 1] * BB[1, 0];
            v2[0, 2] = BB[0, 2] * p[0, 0] * BB[0, 0] + BB[1, 2] * p[1, 1] * BB[1, 0];
            v2[1, 0] = BB[0, 0] * p[0, 0] * BB[0, 1] + BB[1, 0] * p[1, 1] * BB[1, 1];
            v2[1, 1] = BB[0, 1] * p[0, 0] * BB[0, 1] + BB[1, 1] * p[1, 1] * BB[1, 1];
            v2[1, 2] = BB[0, 2] * p[0, 0] * BB[0, 1] + BB[1, 2] * p[1, 1] * BB[1, 1];
            v2[2, 0] = BB[0, 0] * p[0, 0] * BB[0, 2] + BB[1, 0] * p[1, 1] * BB[1, 2];
            v2[2, 1] = BB[0, 1] * p[0, 0] * BB[0, 2] + BB[1, 1] * p[1, 1] * BB[1, 2];
            v2[2, 2] = BB[0, 2] * p[0, 0] * BB[0, 2] + BB[1, 2] * p[1, 1] * BB[1, 2];
        }
        public void Pre_analysis_image_mark(double[,] tie_com,double [,] normal, List<int>[] visibility_point_camera, List<int>[] visibility_camera_point, string[] image_name_Eterior, double[,] Exterior_Orientation_Eterior,   double[] interior, double[,] im_ob_sig, double[,] GPS_weight, out double[,] point_accuracy_complete, out List<string> List_of_removed_images, out List<double> List_of_Increasing_Errors)
        {
            double pi = Math.PI;
            int camera_point_taboo = 4;
            List<double> camera_weight = new List<double>();
            point_accuracy_complete = new double[tie_com.GetLength(0), tie_com.GetLength(1)];
            List_of_removed_images = new List<string>();
            List_of_Increasing_Errors = new List<double>();
            double[,] image_score_taboo = new double[image_name_Eterior.Length, 2];
            double[,] points_v = new double[tie_com.GetLength(0),9];
            double[,] points_v_complete = new double[tie_com.GetLength(0), 9];
            for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
            {
                camera_weight.Add(1.0 / (interior[11] * interior[11]));
                camera_weight.Add(1.0 / (interior[11] * interior[11]));
                camera_weight.Add(1.0 / (GPS_weight[0, 0] * GPS_weight[0, 0]));
                camera_weight.Add(1.0 / (GPS_weight[0, 1] * GPS_weight[0, 1]));
                camera_weight.Add(1.0 / (GPS_weight[0, 2] * GPS_weight[0, 2]));
                camera_weight.Add(1.0 / (GPS_weight[0, 3] * GPS_weight[0, 3]));
                camera_weight.Add(1.0 / (GPS_weight[0, 4] * GPS_weight[0, 4]));
                camera_weight.Add(1.0 / (GPS_weight[0, 5] * GPS_weight[0, 5]));
            }
            double[,] image_m = new double[image_name_Eterior.Length, 9];
            double[,] image_sin_cos = new double[image_name_Eterior.Length, 6];
            for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
            {
              double om = Exterior_Orientation_Eterior[j2, 3];
                double phi = Exterior_Orientation_Eterior[j2, 4];
                double k = Exterior_Orientation_Eterior[j2, 5];
                image_sin_cos[j2, 0] = Math.Sin(om); image_sin_cos[j2, 1] = Math.Sin(phi); image_sin_cos[j2, 2] = Math.Sin(k);
                image_sin_cos[j2, 3] = Math.Cos(om); image_sin_cos[j2, 4] = Math.Cos(phi); image_sin_cos[j2, 5] = Math.Cos(k);
                double sinom = Math.Sin(om); double sinphi = Math.Sin(phi); double sink = Math.Sin(k);
                double cosom = Math.Cos(om); double cosphi = Math.Cos(phi); double cosk = Math.Cos(k);
                //      MessageBox.Show("" + X+"   "+Y+"    "+Z);
                image_m[j2, 0] = (cosphi) * (cosk);
                image_m[j2, 1] = ((sinom) * (sinphi) * (cosk)) + ((cosom) * (sink));
                //m13=(sin(om)*sin(k))-(sin(phi)*cos(om)*cos(k));
                image_m[j2, 2] = ((sinom) * (sink)) - ((sinphi) * (cosom) * (cosk));
                //m21=-cos(phi)*sin(k)
                image_m[j2, 3] = -(cosphi) * (sink);
                //m22=(cos(om)*cos(k))-(sin(om)*sin(phi)*sin(k))
                image_m[j2, 4] = (cosom * cosk) - (sinom * sinphi * sink);
                //m23=(sin(om)*cos(k))+(cos(om)*sin(phi)*sin(k))
                image_m[j2, 5] = (sinom * cosk) + (cosom * sinphi * sink);
                //m31=sin(phi)
                image_m[j2, 6] = (sinphi);
                //m32=-sin(om)*cos(phi)
                image_m[j2, 7] = -sinom * cosphi;
                //m33=cos(om)*cos(phi)
                image_m[j2, 8] = cosom * cosphi;
            }
            //..................compare

            double point_complete_accuracy = 0;
            for (int ij = 0; ij < visibility_point_camera.GetLength(0); ij++)
            {
                if (Math.Acos(normal[ij, 2]) > (pi / 3) && visibility_point_camera[ij].Count < camera_point_taboo+3)
                {
                    for (int ji = 0; ji < visibility_point_camera[ij].Count; ji++)
                    {
                        image_score_taboo[visibility_point_camera[ij][ji], 1] = 1;
                    }
                }
                   
                    else if (visibility_point_camera[ij].Count < camera_point_taboo)
                    {
                        for(int ji=0;ji< visibility_point_camera[ij].Count; ji++)
                        {
                            image_score_taboo[visibility_point_camera[ij][ji], 1] = 1;
                        }
                    }
                    else
                    {

                        double[,] v = new double[3, 3];
                        for (int ji = 0; ji < visibility_point_camera[ij].Count; ji++)
                        {
                            double[,] v22 = new double[3, 3];
                            v_calculator(ij, visibility_point_camera[ij][ji], camera_weight, Exterior_Orientation_Eterior, interior, tie_com, image_m, image_sin_cos, im_ob_sig, out v22);
                            v[0, 0] = v[0, 0] + v22[0, 0]; v[1, 0] = v[1, 0] + v22[1, 0]; v[2, 0] = v[2, 0] + v22[2, 0];
                            v[0, 1] = v[0, 1] + v22[0, 1]; v[1, 1] = v[1, 1] + v22[1, 1]; v[2, 1] = v[2, 1] + v22[2, 1];
                            v[0, 2] = v[0, 2] + v22[0, 2]; v[1, 2] = v[1, 2] + v22[1, 2]; v[2, 2] = v[2, 2] + v22[2, 2];
                            points_v[ij, 0] = points_v[ij, 0] + v22[0, 0]; points_v[ij, 1] = points_v[ij, 1] + v22[0, 1]; points_v[ij, 2] = points_v[ij, 2] + v22[0, 2];
                            points_v[ij, 3] = points_v[ij, 3] + v22[1, 0]; points_v[ij, 4] = points_v[ij, 4] + v22[1, 1]; points_v[ij, 5] = points_v[ij, 5] + v22[1, 2];
                            points_v[ij, 6] = points_v[ij, 6] + v22[2, 0]; points_v[ij, 7] = points_v[ij, 7] + v22[2, 1]; points_v[ij, 8] = points_v[ij, 8] + v22[2, 2];
                            points_v_complete[ij, 0] = points_v_complete[ij, 0] + v22[0, 0]; points_v_complete[ij, 1] = points_v_complete[ij, 1] + v22[0, 1]; points_v_complete[ij, 2] = points_v_complete[ij, 2] + v22[0, 2];
                            points_v_complete[ij, 3] = points_v_complete[ij, 3] + v22[1, 0]; points_v_complete[ij, 4] = points_v_complete[ij, 4] + v22[1, 1]; points_v_complete[ij, 5] = points_v_complete[ij, 5] + v22[1, 2];
                            points_v_complete[ij, 6] = points_v_complete[ij, 6] + v22[2, 0]; points_v_complete[ij, 7] = points_v_complete[ij, 7] + v22[2, 1]; points_v_complete[ij, 8] = points_v_complete[ij, 8] + v22[2, 2];

                        }

                        double[,] invv = new double[3, 3];
                        double detv = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                        invv[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                        invv[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                        invv[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                        invv[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                        invv[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                        invv[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                        invv[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                        invv[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                        invv[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                        point_accuracy_complete[ij, 0] = Math.Sqrt(invv[0, 0]);
                        point_accuracy_complete[ij, 1] = Math.Sqrt(invv[1, 1]);
                        point_accuracy_complete[ij, 2] = Math.Sqrt(invv[2, 2]);
                    point_complete_accuracy = point_complete_accuracy + ((point_accuracy_complete[ij, 0] + point_accuracy_complete[ij, 1] + point_accuracy_complete[ij, 2]) / 3);
                    for (int ji = 0; ji < visibility_point_camera[ij].Count; ji++)
                        {
                            double[,] vp = new double[3, 3];
                            double[,] v22 = new double[3, 3];
                            v_calculator(ij, visibility_point_camera[ij][ji], camera_weight, Exterior_Orientation_Eterior, interior, tie_com, image_m, image_sin_cos, im_ob_sig, out v22);
                            vp[0, 0] = v[0, 0] - v22[0, 0]; vp[1, 0] = v[1, 0] - v22[1, 0]; vp[2, 0] = v[2, 0] - v22[2, 0];
                            vp[0, 1] = v[0, 1] - v22[0, 1]; vp[1, 1] = v[1, 1] - v22[1, 1]; vp[2, 1] = v[2, 1] - v22[2, 1];
                            vp[0, 2] = v[0, 2] - v22[0, 2]; vp[1, 2] = v[1, 2] - v22[1, 2]; vp[2, 2] = v[2, 2] - v22[2, 2];
                            double[,] invv2 = new double[3, 3];
                            double detv2 = (vp[0, 0] * vp[1, 1] * vp[2, 2] - vp[0, 0] * vp[1, 2] * vp[2, 1] - vp[0, 1] * vp[1, 0] * vp[2, 2] + vp[0, 1] * vp[1, 2] * vp[2, 0] + vp[0, 2] * vp[1, 0] * vp[2, 1] - vp[0, 2] * vp[1, 1] * vp[2, 0]);
                            invv2[0, 0] = (vp[1, 1] * vp[2, 2] - vp[1, 2] * vp[2, 1]) / detv2;
                            invv2[0, 1] = -(vp[0, 1] * vp[2, 2] - vp[0, 2] * vp[2, 1]) / detv2;
                            invv2[0, 2] = (vp[0, 1] * vp[1, 2] - vp[0, 2] * vp[1, 1]) / detv2;
                            invv2[1, 0] = -(vp[1, 0] * vp[2, 2] - vp[1, 2] * vp[2, 0]) / detv2;
                            invv2[1, 1] = (vp[0, 0] * vp[2, 2] - vp[0, 2] * vp[2, 0]) / detv2;
                            invv2[1, 2] = -(vp[0, 0] * vp[1, 2] - vp[0, 2] * vp[1, 0]) / detv2;
                            invv2[2, 0] = (vp[1, 0] * vp[2, 1] - vp[1, 1] * vp[2, 0]) / detv2;
                            invv2[2, 1] = -(vp[0, 0] * vp[2, 1] - vp[0, 1] * vp[2, 0]) / detv2;
                            invv2[2, 2] = (vp[0, 0] * vp[1, 1] - vp[0, 1] * vp[1, 0]) / detv2;
                            image_score_taboo[visibility_point_camera[ij][ji], 0] = image_score_taboo[visibility_point_camera[ij][ji], 0] + (((Math.Sqrt(invv2[0, 0]) - Math.Sqrt(invv[0, 0])) + (Math.Sqrt(invv2[1, 1]) - Math.Sqrt(invv[1, 1])) + (Math.Sqrt(invv2[2, 2]) - Math.Sqrt(invv[2, 2]))) / 3);
                        }
                    }
            }
            point_complete_accuracy = point_complete_accuracy / visibility_point_camera.GetLength(0);
            List<int> removed_images = new List<int>();
            List<string> removed_images_str = new List<string>();
            List<double>[] point_accuracy = new List<double>[visibility_point_camera.GetLength(0)];
            for (int ij = 0; ij < visibility_point_camera.GetLength(0); ij++)
            {
                point_accuracy[ij] = new List<double>();
            }
            List<int>[] point_involve_iteration = new List<int>[image_name_Eterior.Length];
            for (int ij = 0; ij < image_name_Eterior.Length; ij++)
            {
                point_involve_iteration[ij] = new List<int>();
            }
            for (int iteration = 0; iteration < image_name_Eterior.Length; iteration++)
            {
                double min_score = 1000000000000000000;
                int min_image = -1;
                for (int j2 = 0; j2 < image_name_Eterior.Length; j2++)
                {
                    if (image_score_taboo[j2, 1] == 0 && image_score_taboo[j2, 0] < min_score)
                    {
                        min_score = image_score_taboo[j2, 0];
                        min_image = j2;
                    }
                }
                if (min_image == -1)
                {
                  //  MessageBox.Show("");
                    break;
                }
             
                image_score_taboo[min_image, 1] = 1;
                removed_images.Add(min_image);
                removed_images_str.Add(image_name_Eterior[min_image]);
                for (int ij = 0; ij < visibility_camera_point[min_image].Count; ij++)
                {
                    int point0 = visibility_camera_point[min_image][ij];
                    
                    if (Math.Acos(normal[point0, 2]) > (pi / 3) && visibility_point_camera[point0].Count < camera_point_taboo + 3)
                    {
                        for (int ji = 0; ji < visibility_point_camera[point0].Count; ji++)
                        {
                            image_score_taboo[visibility_point_camera[point0][ji], 1] = 1;
                        }
                    }
                    if (visibility_point_camera[point0].Count < camera_point_taboo)
                    {
                        for (int ji = 0; ji < visibility_point_camera[point0].Count; ji++)
                        {
                            image_score_taboo[visibility_point_camera[point0][ji], 1] = 1;
                        }
                    }
                    else
                    {
                        double[,] v = new double[3, 3];
                        v[0, 0] = points_v[point0, 0]; v[0, 1] = points_v[point0, 1]; v[0, 2] = points_v[point0, 2];
                        v[1, 0] = points_v[point0, 3]; v[1, 1] = points_v[point0, 4]; v[1, 2] = points_v[point0, 5];
                        v[2, 0] = points_v[point0, 6]; v[2, 1] = points_v[point0, 7]; v[2, 2] = points_v[point0, 8];
                        double[,] invv = new double[3, 3];
                        double detv = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                        invv[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv;
                        invv[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv;
                        invv[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv;
                        invv[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv;
                        invv[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv;
                        invv[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv;
                        invv[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv;
                        invv[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv;
                        invv[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv;
                        for (int ji = 0; ji < visibility_point_camera[point0].Count; ji++)
                        {
                            double[,] vp = new double[3, 3];
                            double[,] v22 = new double[3, 3];
                            v_calculator(point0, visibility_point_camera[point0][ji], camera_weight, Exterior_Orientation_Eterior, interior, tie_com, image_m, image_sin_cos, im_ob_sig, out v22);
                            vp[0, 0] = v[0, 0] - v22[0, 0]; vp[1, 0] = v[1, 0] - v22[1, 0]; vp[2, 0] = v[2, 0] - v22[2, 0];
                            vp[0, 1] = v[0, 1] - v22[0, 1]; vp[1, 1] = v[1, 1] - v22[1, 1]; vp[2, 1] = v[2, 1] - v22[2, 1];
                            vp[0, 2] = v[0, 2] - v22[0, 2]; vp[1, 2] = v[1, 2] - v22[1, 2]; vp[2, 2] = v[2, 2] - v22[2, 2];
                            double[,] invv2 = new double[3, 3];
                            double detv2 = (vp[0, 0] * vp[1, 1] * vp[2, 2] - vp[0, 0] * vp[1, 2] * vp[2, 1] - vp[0, 1] * vp[1, 0] * vp[2, 2] + vp[0, 1] * vp[1, 2] * vp[2, 0] + vp[0, 2] * vp[1, 0] * vp[2, 1] - vp[0, 2] * vp[1, 1] * vp[2, 0]);
                            invv2[0, 0] = (vp[1, 1] * vp[2, 2] - vp[1, 2] * vp[2, 1]) / detv2;
                            invv2[0, 1] = -(vp[0, 1] * vp[2, 2] - vp[0, 2] * vp[2, 1]) / detv2;
                            invv2[0, 2] = (vp[0, 1] * vp[1, 2] - vp[0, 2] * vp[1, 1]) / detv2;
                            invv2[1, 0] = -(vp[1, 0] * vp[2, 2] - vp[1, 2] * vp[2, 0]) / detv2;
                            invv2[1, 1] = (vp[0, 0] * vp[2, 2] - vp[0, 2] * vp[2, 0]) / detv2;
                            invv2[1, 2] = -(vp[0, 0] * vp[1, 2] - vp[0, 2] * vp[1, 0]) / detv2;
                            invv2[2, 0] = (vp[1, 0] * vp[2, 1] - vp[1, 1] * vp[2, 0]) / detv2;
                            invv2[2, 1] = -(vp[0, 0] * vp[2, 1] - vp[0, 1] * vp[2, 0]) / detv2;
                            invv2[2, 2] = (vp[0, 0] * vp[1, 1] - vp[0, 1] * vp[1, 0]) / detv2;
                            image_score_taboo[visibility_point_camera[point0][ji], 0] = image_score_taboo[visibility_point_camera[point0][ji], 0] - (((Math.Sqrt(invv2[0, 0]) - Math.Sqrt(invv[0, 0])) + (Math.Sqrt(invv2[1, 1]) - Math.Sqrt(invv[1, 1])) + (Math.Sqrt(invv2[2, 2]) - Math.Sqrt(invv[2, 2]))) / 3);
                        }
                        visibility_point_camera[point0].Remove(min_image);
                        double[,] v222 = new double[3, 3];
                        v_calculator(point0, min_image, camera_weight, Exterior_Orientation_Eterior, interior, tie_com, image_m, image_sin_cos, im_ob_sig, out v222);
                        v[0, 0] = v[0, 0] - v222[0, 0]; v[1, 0] = v[1, 0] - v222[1, 0]; v[2, 0] = v[2, 0] - v222[2, 0];
                        v[0, 1] = v[0, 1] - v222[0, 1]; v[1, 1] = v[1, 1] - v222[1, 1]; v[2, 1] = v[2, 1] - v222[2, 1];
                        v[0, 2] = v[0, 2] - v222[0, 2]; v[1, 2] = v[1, 2] - v222[1, 2]; v[2, 2] = v[2, 2] - v222[2, 2];
                        points_v[point0, 0] = points_v[point0, 0] - v222[0, 0]; points_v[point0, 1] = points_v[point0, 1] - v222[0, 1]; points_v[point0, 2] = points_v[point0, 2] - v222[0, 2];
                        points_v[point0, 3] = points_v[point0, 3] - v222[1, 0]; points_v[point0, 4] = points_v[point0, 4] - v222[1, 1]; points_v[point0, 5] = points_v[point0, 5] - v222[1, 2];
                        points_v[point0, 6] = points_v[point0, 6] - v222[2, 0]; points_v[point0, 7] = points_v[point0, 7] - v222[2, 1]; points_v[point0, 8] = points_v[point0, 8] - v222[2, 2];
                        double[,] invv22 = new double[3, 3];
                        double detv22 = (v[0, 0] * v[1, 1] * v[2, 2] - v[0, 0] * v[1, 2] * v[2, 1] - v[0, 1] * v[1, 0] * v[2, 2] + v[0, 1] * v[1, 2] * v[2, 0] + v[0, 2] * v[1, 0] * v[2, 1] - v[0, 2] * v[1, 1] * v[2, 0]);
                        invv22[0, 0] = (v[1, 1] * v[2, 2] - v[1, 2] * v[2, 1]) / detv22;
                        invv22[0, 1] = -(v[0, 1] * v[2, 2] - v[0, 2] * v[2, 1]) / detv22;
                        invv22[0, 2] = (v[0, 1] * v[1, 2] - v[0, 2] * v[1, 1]) / detv22;
                        invv22[1, 0] = -(v[1, 0] * v[2, 2] - v[1, 2] * v[2, 0]) / detv22;
                        invv22[1, 1] = (v[0, 0] * v[2, 2] - v[0, 2] * v[2, 0]) / detv22;
                        invv22[1, 2] = -(v[0, 0] * v[1, 2] - v[0, 2] * v[1, 0]) / detv22;
                        invv22[2, 0] = (v[1, 0] * v[2, 1] - v[1, 1] * v[2, 0]) / detv22;
                        invv22[2, 1] = -(v[0, 0] * v[2, 1] - v[0, 1] * v[2, 0]) / detv22;
                        invv22[2, 2] = (v[0, 0] * v[1, 1] - v[0, 1] * v[1, 0]) / detv22;
                        double point_accuracy_X = Math.Sqrt(invv22[0, 0]);
                        double point_accuracy_Y = Math.Sqrt(invv22[1, 1]);
                        double point_accuracy_Z = Math.Sqrt(invv22[2, 2]);
                        point_accuracy[point0].Add(((point_accuracy_X - point_accuracy_complete[point0, 0]) + (point_accuracy_Y - point_accuracy_complete[point0, 1]) + (point_accuracy_Z - point_accuracy_complete[point0, 2])) / 3);
                        point_involve_iteration[iteration].Add(point0);
                        for (int ji = 0; ji < visibility_point_camera[point0].Count; ji++)
                        {
                            double[,] vp = new double[3, 3];
                            double[,] v22 = new double[3, 3];
                            v_calculator(point0, visibility_point_camera[point0][ji], camera_weight, Exterior_Orientation_Eterior, interior, tie_com, image_m, image_sin_cos, im_ob_sig, out v22);
                            vp[0, 0] = v[0, 0] - v22[0, 0]; vp[1, 0] = v[1, 0] - v22[1, 0]; vp[2, 0] = v[2, 0] - v22[2, 0];
                            vp[0, 1] = v[0, 1] - v22[0, 1]; vp[1, 1] = v[1, 1] - v22[1, 1]; vp[2, 1] = v[2, 1] - v22[2, 1];
                            vp[0, 2] = v[0, 2] - v22[0, 2]; vp[1, 2] = v[1, 2] - v22[1, 2]; vp[2, 2] = v[2, 2] - v22[2, 2];
                            double[,] invv2 = new double[3, 3];
                            double detv2 = (vp[0, 0] * vp[1, 1] * vp[2, 2] - vp[0, 0] * vp[1, 2] * vp[2, 1] - vp[0, 1] * vp[1, 0] * vp[2, 2] + vp[0, 1] * vp[1, 2] * vp[2, 0] + vp[0, 2] * vp[1, 0] * vp[2, 1] - vp[0, 2] * vp[1, 1] * vp[2, 0]);
                            invv2[0, 0] = (vp[1, 1] * vp[2, 2] - vp[1, 2] * vp[2, 1]) / detv2;
                            invv2[0, 1] = -(vp[0, 1] * vp[2, 2] - vp[0, 2] * vp[2, 1]) / detv2;
                            invv2[0, 2] = (vp[0, 1] * vp[1, 2] - vp[0, 2] * vp[1, 1]) / detv2;
                            invv2[1, 0] = -(vp[1, 0] * vp[2, 2] - vp[1, 2] * vp[2, 0]) / detv2;
                            invv2[1, 1] = (vp[0, 0] * vp[2, 2] - vp[0, 2] * vp[2, 0]) / detv2;
                            invv2[1, 2] = -(vp[0, 0] * vp[1, 2] - vp[0, 2] * vp[1, 0]) / detv2;
                            invv2[2, 0] = (vp[1, 0] * vp[2, 1] - vp[1, 1] * vp[2, 0]) / detv2;
                            invv2[2, 1] = -(vp[0, 0] * vp[2, 1] - vp[0, 1] * vp[2, 0]) / detv2;
                            invv2[2, 2] = (vp[0, 0] * vp[1, 1] - vp[0, 1] * vp[1, 0]) / detv2;
                            image_score_taboo[visibility_point_camera[point0][ji], 0] = image_score_taboo[visibility_point_camera[point0][ji], 0] + (((Math.Sqrt(invv2[0, 0]) - Math.Sqrt(invv[0, 0])) + (Math.Sqrt(invv2[1, 1]) - Math.Sqrt(invv[1, 1])) + (Math.Sqrt(invv2[2, 2]) - Math.Sqrt(invv[2, 2]))) / 3);
                        }
                    }
                }
            }
            Hashtable involve_points = new Hashtable();
            double mean_error_all = 0;
        //    System.IO.StreamWriter res1 = new System.IO.StreamWriter(@"D:\import_parameters\1\Export\result1.txt");
        //    System.IO.StreamWriter removed_images20 = new System.IO.StreamWriter(@"D:\import_parameters\1\Export\remooved_images.txt");
            for (int i = 0; i < removed_images_str.Count; i++)
            {
                //         removed_images20.WriteLine("" + removed_images_str[i]);
                List_of_removed_images.Add(removed_images_str[i]);
                
            }
            for (int iteration = 0; iteration < removed_images.Count; iteration++)
            {
                double mean_error = 0;
                int Points_involve;
                for (int i = 0; i < point_involve_iteration[iteration].Count ; i++)
                {
                    // mean_error = point_accuracy[point_involve_iteration[iteration][i]][0] + mean_error;
                    //  point_accuracy[point_involve_iteration[iteration][i]].Remove(point_accuracy[point_involve_iteration[iteration][i]][0]);
                    if(involve_points.ContainsKey("" + point_involve_iteration[iteration][i]))
                    {
                        involve_points["" + point_involve_iteration[iteration][i]] =Convert.ToInt32( involve_points["" + point_involve_iteration[iteration][i]]) + 1;
                    }
                    else
                    {
                        involve_points["" + point_involve_iteration[iteration][i]] = 0;
                    }
                    
                }
                foreach (DictionaryEntry de in involve_points)
                {
                    mean_error = point_accuracy[Convert.ToInt32( de.Key)][Convert.ToInt32(de.Value)] + mean_error;
                 //   point_accuracy[Convert.ToInt32(de.Key)].RemoveAt(0);
                }
            
                Points_involve = involve_points.Count;
                mean_error = mean_error / involve_points.Count;
                mean_error_all = mean_error;
                List_of_Increasing_Errors.Add(mean_error * 100 / point_complete_accuracy);
                //      res1.WriteLine("" + mean_error * 100 / point_complete_accuracy + " " + Points_involve);
            }
        //    res1.Close();
         //   removed_images20.Close();
        }
    }
}
