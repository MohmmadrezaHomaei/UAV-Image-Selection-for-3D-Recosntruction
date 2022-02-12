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
using System.IO.MemoryMappedFiles;
using System.Diagnostics;
using SVD_M;
using MathWorks.MATLAB.NET.Utility;
using MathWorks.MATLAB.NET.Arrays;
using System.Threading;
namespace theses
{
    class Grid
    {
        public static bool Intersect(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p0, Vector3 p00)
        {
            double Epsilon = 0.00000001;
            // Vectors from p1 to p2/p3 (edges)
            Vector3 e1, e2;

            Vector3 p, q, t;
            float det, invDet, u, v;


            //Find vectors for two edges sharing vertex/point p1
            e1 = p2 - p1;
            e2 = p3 - p1;

            // calculating determinant 
            p = Vector3.Cross(p0 - p00, e2);

            //Calculate determinat
            det = Vector3.Dot(e1, p);

            //if determinant is near zero, ray lies in plane of triangle otherwise not
            if (det > -Epsilon && det < Epsilon) { return false; }
            invDet = 1.0f / det;

            //calculate distance from p1 to ray origin
            t = p00 - p1;

            //Calculate u parameter
            u = Vector3.Dot(t, p) * invDet;

            //Check for ray hit
            if (u < 0 || u > 1) { return false; }

            //Prepare to test v parameter
            q = Vector3.Cross(t, e1);

            //Calculate v parameter
            v = Vector3.Dot(p0 - p00, q) * invDet;

            //Check for ray hit
            if (v < 0 || u + v > 1) { return false; }

            if ((Vector3.Dot(e2, q) * invDet) > Epsilon)
            {
                //ray does intersect
                return true;
            }
            return false;
        }

        public void threads_visibility_triangulation(int ii, int ie, int[,] triangle_int, double[,] triangle_coordinates, string[] image_name_Eterior, double[,] Exterior_Orientation_Eterior, double[] interior, out List<string> Visibility_str)
        {
            double x_image_maximim = (interior[10] / 2) * interior[11];
            double y_image_maximim = (interior[12] / 2) * interior[13];
            List<int>[] point_triangles = new List<int>[triangle_coordinates.GetLength(0)];
            for (int ij = 0; ij < triangle_coordinates.GetLength(0); ij++)
            {
                point_triangles[ij] = new List<int>();
            }
            for (int ij = 0; ij < triangle_int.GetLength(0); ij++)
            {
                point_triangles[triangle_int[ij, 0]].Add(ij);
                point_triangles[triangle_int[ij, 1]].Add(ij);
                point_triangles[triangle_int[ij, 2]].Add(ij);
            }
            List<int> test = new List<int>();
            double teta = 10;
            Visibility_str = new List<string>();
            double pi = Math.PI;
            for (int j2 = ii; j2 < ie; j2++)
            {

                double XL = Exterior_Orientation_Eterior[j2, 0];
                double YL = Exterior_Orientation_Eterior[j2, 1];
                double ZL = Exterior_Orientation_Eterior[j2, 2];
                double om = Exterior_Orientation_Eterior[j2, 3];
                double phi = Exterior_Orientation_Eterior[j2, 4];
                double k = Exterior_Orientation_Eterior[j2, 5];
                double ff = interior[0];


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
                double[,] triangle_coordinates_polar = new double[triangle_coordinates.GetLength(0), 4];
                //   List<string> image_points_visibility_str = new List<string>();
                int teta_shom = Convert.ToInt32(360 / teta);
                List<double>[,] teta_point = new List<double>[teta_shom, 3];
                for (int ij = 0; ij < teta_shom; ij++)
                {
                    teta_point[ij, 0] = new List<double>(); teta_point[ij, 1] = new List<double>(); teta_point[ij, 2] = new List<double>();
                }
                List<int> all_points = new List<int>();
                //System.IO.StreamWriter l1 = new System.IO.StreamWriter(@"E:\test\teta_test1.txt");
                //System.IO.StreamWriter l2 = new System.IO.StreamWriter(@"E:\test\teta_test2.txt");
                //System.IO.StreamWriter l3 = new System.IO.StreamWriter(@"E:\test\teta_test3.txt");
                //System.IO.StreamWriter l4 = new System.IO.StreamWriter(@"E:\test\teta_test4.txt");
                //System.IO.StreamWriter l5 = new System.IO.StreamWriter(@"E:\test\teta_test5.txt");
                //System.IO.StreamWriter l6 = new System.IO.StreamWriter(@"E:\test\teta_test6.txt");
                //System.IO.StreamWriter l7 = new System.IO.StreamWriter(@"E:\test\teta_test7.txt");
                //System.IO.StreamWriter l8 = new System.IO.StreamWriter(@"E:\test\teta_test8.txt");
                //System.IO.StreamWriter l9 = new System.IO.StreamWriter(@"E:\test\teta_test9.txt");
                //System.IO.StreamWriter l0 = new System.IO.StreamWriter(@"E:\test\teta_test0.txt");
                //System.IO.StreamWriter l10 = new System.IO.StreamWriter(@"E:\test\teta_test20.txt");
                //System.IO.StreamWriter l11 = new System.IO.StreamWriter(@"E:\test\teta_test16.txt");
                for (int i = 0; i < triangle_coordinates.GetLength(0); i++)
                {
                    double X = triangle_coordinates[i, 0];
                    double Y = triangle_coordinates[i, 1];
                    double Z = triangle_coordinates[i, 2];
                    double xx = -interior[0] * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                    double yy = -interior[0] * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                    double L = Math.Sqrt(((X - XL) * (X - XL)) + ((Y - YL) * (Y - YL)) + ((Z - ZL) * (Z - ZL)));
                    double ZZ11 = interior[0] * L / (Math.Sqrt((xx * xx) + (yy * yy) + (interior[0] * interior[0])));
                    double XX11 = -xx * ZZ11 / interior[0];
                    double YY11 = -yy * ZZ11 / interior[0];
                    triangle_coordinates_polar[i, 0] = (L * xx) / (Math.Sqrt((xx * xx) + (yy * yy) + (interior[0] * interior[0])));
                    triangle_coordinates_polar[i, 2] = -interior[0] * triangle_coordinates_polar[i, 0] / xx;
                    triangle_coordinates_polar[i, 1] = triangle_coordinates_polar[i, 2] * yy / (-interior[0]);
                    if (xx < x_image_maximim && xx > -x_image_maximim && yy < y_image_maximim && yy > -y_image_maximim)
                    {
                        double XX = triangle_coordinates_polar[i, 0], YY = triangle_coordinates_polar[i, 1],ZZ = triangle_coordinates_polar[i, 2];
                        double atan = Math.Atan2(triangle_coordinates_polar[i, 1], triangle_coordinates_polar[i, 0]) * 180 / pi;
                       
                        if (atan < 0)
                        {

                            atan = atan + 360;
                        }

                        double R = Math.Sqrt(((triangle_coordinates_polar[i, 0]) * (triangle_coordinates_polar[i, 0])) + ((triangle_coordinates_polar[i, 1]) * (triangle_coordinates_polar[i, 1])));
                        //    double grid_r = Math.Ceiling(Math.Sqrt(((X - XL) * (X - XL)) + ((Y - YL) * (Y - YL))) / r_interval);
                        int grid_teta = Convert.ToInt32(Math.Ceiling(atan / teta));
                        teta_point[grid_teta - 1, 0].Add(i);
                        teta_point[grid_teta - 1, 1].Add(R);
                        teta_point[grid_teta - 1, 2].Add(R/ triangle_coordinates_polar[i, 2]);
                    }
                }

                
                //    for (int yu = 0; yu < teta_point[0, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[0, 0][yu]);
                //        l0.WriteLine(triangle_coordinates[shom, 0]+" "+ triangle_coordinates[shom, 1]+" "+ triangle_coordinates[shom, 2]);
                //    }
                
               
                //    for (int yu = 0; yu < teta_point[17, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[17, 0][yu]);
                //        l1.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
                
                //    for (int yu = 0; yu < teta_point[2, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[2, 0][yu]);
                //        l2.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
               
                //    for (int yu = 0; yu < teta_point[3, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[3, 0][yu]);
                //        l3.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
               
                //    for (int yu = 0; yu < teta_point[4, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[4, 0][yu]);
                //        l4.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
               
                //    for (int yu = 0; yu < teta_point[5, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[5, 0][yu]);
                //        l5.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
               
                //    for (int yu = 0; yu < teta_point[6, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[6, 0][yu]);
                //        l6.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
               
                //    for (int yu = 0; yu < teta_point[7, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[7, 0][yu]);
                //        l7.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
              
                //    for (int yu = 0; yu < teta_point[8, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[8, 0][yu]);
                //        l8.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                
              
                //    for (int yu = 0; yu < teta_point[9, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[9, 0][yu]);
                //        l9.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
             
              
                //    for (int yu = 0; yu < teta_point[15, 0].Count; yu++)
                //    {
                //        int shom = Convert.ToInt32(teta_point[15, 0][yu]);
                //        l10.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //    }
                //for (int yu = 0; yu < teta_point[16, 0].Count; yu++)
                //{
                //    int shom = Convert.ToInt32(teta_point[16, 0][yu]);
                //    l11.WriteLine(triangle_coordinates[shom, 0] + " " + triangle_coordinates[shom, 1] + " " + triangle_coordinates[shom, 2]);
                //}
                //l0.Close(); l1.Close(); l2.Close(); l3.Close(); l4.Close();
                //l5.Close(); l6.Close(); l7.Close(); l8.Close(); l9.Close(); l10.Close(); l11.Close();
                System.IO.StreamWriter ll1l = new System.IO.StreamWriter(@"E:\test\point visible.txt");
                System.IO.StreamWriter ll1 = new System.IO.StreamWriter(@"E:\test\point visible1.txt");
                for (int i = 0; i < triangle_coordinates.GetLength(0); i++)
                {
                    
                    double X = triangle_coordinates[i, 0];
                    double Y = triangle_coordinates[i, 1];
                    double Z = triangle_coordinates[i, 2];
                    //if (X==1605383.818927 && Y==14726563.964230 && Z == 544.496502)
                    //{
                    //    MessageBox.Show("d");
                    //}
                    double xx = -interior[0] * (((m11 * (X - XL)) + (m12 * (Y - YL)) + (m13 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                    double yy = -interior[0] * (((m21 * (X - XL)) + (m22 * (Y - YL)) + (m23 * (Z - ZL))) / ((m31 * (X - XL)) + (m32 * (Y - YL)) + (m33 * (Z - ZL))));
                    double X_polar = 0;
                    double Y_polar = 0;
                    double Z_polar = 0;
                    if (xx < x_image_maximim && xx > -x_image_maximim && yy < y_image_maximim && yy > -y_image_maximim)
                    {
                        bool tri_ray_inter=false;
                        double L = Math.Sqrt(((X - XL) * (X - XL)) + ((Y - YL) * (Y - YL)) + ((Z - ZL) * (Z - ZL)));
                        X_polar = triangle_coordinates_polar[i, 0];

                        Z_polar = triangle_coordinates_polar[i, 2];
                        Y_polar = triangle_coordinates_polar[i, 1];
                        double atan = Math.Atan2(Y_polar, X_polar) * 180 / pi;
                        if (atan < 0)
                        {
                            atan = atan + 360;
                        }
                        double R = Math.Sqrt(((X_polar) * (X_polar)) + ((Y_polar) * (Y_polar)));
                        //    double grid_r = Math.Ceiling(Math.Sqrt(((X - XL) * (X - XL)) + ((Y - YL) * (Y - YL))) / r_interval);
                        int grid_teta = Convert.ToInt32(Math.Ceiling(atan / teta));
                        double atan_v = R/ Z_polar;
                        int [] i_tri = new int[teta_point[grid_teta - 1, 0].Count];
                        double[] R_tri = new double[teta_point[grid_teta - 1, 1].Count];
                        double[] atan_tri = new double[teta_point[grid_teta - 1, 2].Count];
                        for (int ji=0;ji< teta_point[grid_teta - 1, 0].Count; ji++)
                         {
                            i_tri[ji] =Convert.ToInt32( teta_point[grid_teta - 1, 0][ji]);
                            R_tri[ji] = teta_point[grid_teta - 1, 1][ji];
                            atan_tri[ji] = teta_point[grid_teta - 1, 2][ji];
                            if (R_tri[ji] < R)
                            {
                               // if (atan_tri[ji] > atan_v)
                               // {
                                    for (int ji2 = 0; ji2 < point_triangles[i_tri[ji]].Count; ji2++)
                                    {
                                        int g1 = i_tri[ji];
                                        int g2 = point_triangles[i_tri[ji]][ji2];
                                        int g3 = triangle_int[point_triangles[i_tri[ji]][ji2], 0];
                                        int g4 = triangle_int[point_triangles[i_tri[ji]][ji2], 1];
                                        int g5 = triangle_int[point_triangles[i_tri[ji]][ji2], 2];
                                  
                                        double g10 = triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 0], 0];
                                        Vector3 p1 = new Vector3();
                                        p1.X = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 0], 0]);
                                        p1.Y = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 0], 1]);
                                        p1.Z = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 0], 2] );
                                        Vector3 p2 = new Vector3();
                                        p2.X = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 1], 0]);
                                        p2.Y = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 1], 1]);
                                        p2.Z = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 1], 2] );
                                        Vector3 p3 = new Vector3();
                                        p3.X = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 2], 0]);
                                        p3.Y = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 2], 1]);
                                        p3.Z = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri[ji]][ji2], 2], 2] );
                                        Vector3 p00 = new Vector3();
                                        p00.X = Convert.ToSingle(X_polar);
                                        p00.Y = Convert.ToSingle(Y_polar);
                                        p00.Z = Convert.ToSingle(Z_polar);
                                        Vector3 p0 = new Vector3();
                                        p0.X = 0;
                                        p0.Y = 0;
                                        p0.Z = 0;
                                        tri_ray_inter = Intersect(p1, p2, p3, p0, p00);
                                        if (tri_ray_inter == true)
                                        {
                                            break;
                                        }
                                    }
                              //  }
                            }
                            if (tri_ray_inter == true)
                            {
                                ll1l.WriteLine("" + triangle_coordinates[i, 0] + " " + triangle_coordinates[i, 1] + " " + triangle_coordinates[i, 2]);
                                break;
                            }
                            if (ji == teta_point[grid_teta - 1, 0].Count-1)
                            {
                                ll1.WriteLine("" + triangle_coordinates[i, 0] + " " + triangle_coordinates[i, 1] + " " + triangle_coordinates[i, 2]);
                                Visibility_str.Add("" + i + " " + j2);
                            }
                        }
                   
                    }
                }
                ll1l.Close();
                ll1.Close();
            }
        }
        public void All_tie_points_geometric_visibility_saadat_all_triangulation(int[,] triangle_int, double[,] triangle_coordinates, string[] image_name_Eterior, double[,] Exterior_Orientation_Eterior, double[] interior, out List<int>[] visibility_point_camera, out List<int>[] visibility_camera_point)
        {
            List<string> Visibility_str0 = new List<string>();
            //   var s10 = Stopwatch.StartNew();
            //  threads_visibility(0, image_name_Eterior.Length, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str0);
            int odd = Convert.ToInt32(image_name_Eterior.Length / 8);
            List<string> Visibility_str1 = new List<string>();
            List<string> Visibility_str2 = new List<string>();
            List<string> Visibility_str3 = new List<string>();
            List<string> Visibility_str4 = new List<string>();
            List<string> Visibility_str5 = new List<string>();
            List<string> Visibility_str6 = new List<string>();
            List<string> Visibility_str7 = new List<string>();
            List<string> Visibility_str8 = new List<string>();
            threads_visibility_triangulation(0, image_name_Eterior.Length, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str1);
            //Thread t1 = new Thread(() => threads_visibility_triangulation(0, odd, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str1));
            //Thread t2 = new Thread(() => threads_visibility_triangulation(odd, odd * 2, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str2));
            //Thread t3 = new Thread(() => threads_visibility_triangulation(odd * 2, odd * 3, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str3));
            //Thread t4 = new Thread(() => threads_visibility_triangulation(odd * 3, odd * 4, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str4));
            //Thread t5 = new Thread(() => threads_visibility_triangulation(odd * 4, odd * 5, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str5));
            //Thread t6 = new Thread(() => threads_visibility_triangulation(odd * 5, odd * 6, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str6));
            //Thread t7 = new Thread(() => threads_visibility_triangulation(odd * 6, odd * 7, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str7));
            //Thread t8 = new Thread(() => threads_visibility_triangulation(odd * 7, image_name_Eterior.Length, triangle_int, triangle_coordinates, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str8));
            //t1.Start(); t2.Start(); t3.Start(); t4.Start(); t5.Start(); t6.Start(); t7.Start(); t8.Start();
            //t1.Join(); t2.Join(); t3.Join(); t4.Join(); t5.Join(); t6.Join(); t7.Join(); t8.Join();
            visibility_point_camera = new List<int>[triangle_coordinates.GetLength(0)];
            for (int ij = 0; ij < triangle_coordinates.GetLength(0); ij++)
            {
                visibility_point_camera[ij] = new List<int>();
            }
            visibility_camera_point = new List<int>[image_name_Eterior.GetLength(0)];
            for (int ij = 0; ij < image_name_Eterior.GetLength(0); ij++)
            {
                visibility_camera_point[ij] = new List<int>();
            }
            for (int i = 0; i < Visibility_str1.Count; i++)
            {
                string[] tokens = Visibility_str1[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
               
            }
            for (int i = 0; i < Visibility_str2.Count; i++)
            {
                string[] tokens = Visibility_str2[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
              
            }
            for (int i = 0; i < Visibility_str3.Count; i++)
            {
                string[] tokens = Visibility_str3[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
            
            }
            for (int i = 0; i < Visibility_str4.Count; i++)
            {
                string[] tokens = Visibility_str4[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
                
            }
            for (int i = 0; i < Visibility_str5.Count; i++)
            {
                string[] tokens = Visibility_str5[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
              
            }
            for (int i = 0; i < Visibility_str6.Count; i++)
            {
                string[] tokens = Visibility_str6[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
             
            }
            for (int i = 0; i < Visibility_str7.Count; i++)
            {
                string[] tokens = Visibility_str7[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
               
            }
            for (int i = 0; i < Visibility_str8.Count; i++)
            {
                string[] tokens = Visibility_str8[i].Split(new[] { " " }, StringSplitOptions.None);
                int point = Convert.ToInt32(tokens[0]);
                int camera = Convert.ToInt32(tokens[1]);
                visibility_point_camera[point].Add(camera);
                visibility_camera_point[camera].Add(point);
            }
    
            //s10.Stop();
            //MessageBox.Show("" + s10.ElapsedMilliseconds);

        }
        public void triangles_visibility(int[,] triangle_int, double[,] triangle_coordinates,double [,] normal, double[,] visible_point, out int[,] Visibility)
        {
            double[,] centroid_triangle_coordinates = new double[triangle_int.GetLength(0), 3];
            double pi = Math.PI;
            double min_distance = 20;
            double max_distance = 100;
            double landa_interval = 10;
            double phi_interval = 5;
            List<int>[] point_triangles = new List<int>[triangle_coordinates.GetLength(0)];
            for (int ij = 0; ij < triangle_coordinates.GetLength(0); ij++)
            {
                point_triangles[ij] = new List<int>();
            }
            for (int ij = 0; ij < triangle_int.GetLength(0); ij++)
            {

                point_triangles[triangle_int[ij, 0]].Add(ij);
                point_triangles[triangle_int[ij, 1]].Add(ij);
                point_triangles[triangle_int[ij, 2]].Add(ij);
            }
            for (int ij = 0; ij < triangle_coordinates.GetLength(0); ij++)
            {
                double phi = Math.Asin(normal[ij, 0]);
                double a = -normal[ij, 1] * Math.Cos(phi);
                double b = normal[ij, 2] * Math.Cos(phi);
            //    double tr = Math.Sin(phi);
                double c =1 -normal[ij, 0] * Math.Sin(phi);
                double om = 0;
                if ((((a * c) + (b * Math.Sqrt((a * a) + (b * b) - (c * c)))) / ((a * a) + (b * b)))<=1 && (((a * c) + (b * Math.Sqrt((a * a) + (b * b) - (c * c)))) / ((a * a) + (b * b))) >= -1)
                {
                    om = Math.Asin((((a * c) + (b * Math.Sqrt((a * a) + (b * b) - (c * c)))) / ((a * a) + (b * b))));
                }
                else
                {
                    om = Math.Asin((((a * c) - (b * Math.Sqrt((a * a) + (b * b) - (c * c)))) / ((a * a) + (b * b))));
                }
            //    double om =Math.Asin (((-4 * a * c) + (2 * b * Math.Sqrt((a * a) + (b * b) - (c * c)))) / (2 * ((a * a) + (b * b))));
                double a2 = -normal[ij, 0] * Math.Cos(phi);
                double b2 = normal[ij, 1] * Math.Cos(om);
                double c2 = -normal[ij, 1] * Math.Sin(phi)* Math.Sin(om);
                double d2 = normal[ij, 2] * Math.Sin(om);
                double e2 = normal[ij, 2] * Math.Sin(phi) * Math.Cos(om);
                double k = Math.Atan(-(b2 + d2) / (a2 + c2 + e2));
                double R11 = Math.Cos(phi) * Math.Cos(k);
                double R12 = (Math.Cos(om) * Math.Sin(k)) + (Math.Sin(om) * Math.Sin(phi) * Math.Cos(k));
                double R13 = (Math.Sin(om) * Math.Sin(k)) - (Math.Cos(om) * Math.Sin(phi) * Math.Cos(k));
                double R21 = -Math.Cos(phi) * Math.Sin(k);
                double R22 = (Math.Cos(om) * Math.Cos(k)) - (Math.Sin(om) * Math.Sin(phi) * Math.Sin(k));
                double R23 = (Math.Sin(om) * Math.Cos(k)) + (Math.Cos(om) * Math.Sin(phi) * Math.Sin(k));
                double R31 = Math.Sin(phi);
                double R32 =  -(Math.Sin(om) * Math.Cos(phi));
                double R33 =  (Math.Cos(om) * Math.Cos(phi) );
                double x0= -((triangle_coordinates[ij, 0] * R11) + (triangle_coordinates[ij, 1] * R12) + (triangle_coordinates[ij, 2] * R13));
                double y0 = -((triangle_coordinates[ij, 0] * R21) + (triangle_coordinates[ij, 1] * R22) + (triangle_coordinates[ij, 2] * R23));
                double z0 = -((triangle_coordinates[ij, 0] * R31) + (triangle_coordinates[ij, 1] * R32) + (triangle_coordinates[ij, 2] * R33));
                double teta = 1;
                List<double>[,] polar_triangles = new List<double>[(int) Math.Ceiling( 360/teta),2];
                for (int ijj = 0; ijj < polar_triangles.GetLength(0); ijj++)
                {
                    polar_triangles[ijj,0] = new List<double>(); polar_triangles[ijj, 1] = new List<double>();
                }
                double[,] triangle_coordinates_polar = new double[triangle_coordinates.GetLength(0),3];
                for (int ij2 = 0; ij2 < triangle_coordinates.GetLength(0); ij2++)
                {
                    double X = ((triangle_coordinates[ij2, 0] * R11) + (triangle_coordinates[ij2, 1] * R12) + (triangle_coordinates[ij2, 2] * R13)) + x0;
                    double Y = ((triangle_coordinates[ij2, 0] * R21) + (triangle_coordinates[ij2, 1] * R22) + (triangle_coordinates[ij2, 2] * R23)) + y0;
                    double Z = ((triangle_coordinates[ij2, 0] * R31) + (triangle_coordinates[ij2, 1] * R32) + (triangle_coordinates[ij2, 2] * R33)) + z0;
                    triangle_coordinates_polar[ij2, 0] = X;
                    triangle_coordinates_polar[ij2, 1] = Y;
                    triangle_coordinates_polar[ij2, 2] = Z;
                    if (Z > 0)
                    {
                        double atan = Math.Atan2(Y, X) * 180 / pi;

                        if (atan < 0)
                        {

                            atan = atan + 360;
                        }

                        double R = Math.Sqrt(((X) * (X)) + ((Y) * (Y)));
                        //    double grid_r = Math.Ceiling(Math.Sqrt(((X - XL) * (X - XL)) + ((Y - YL) * (Y - YL))) / r_interval);
                        int grid_teta = Convert.ToInt32(Math.Ceiling(atan / teta));
                        polar_triangles[grid_teta - 1, 0].Add(ij2);
                        polar_triangles[grid_teta - 1, 1].Add(R);
                    }
                }
                List<double> landa_s = new List<double>();
                List<double> phi_s = new List<double>();
                landa_s.Add(0);
                phi_s.Add(90);
                for (int i2 = 0; i2 < (int)Math.Ceiling(90 / phi_interval); i2++)
                {
                 double lan_cal2 =360 / (Math.Cos(i2 * phi_interval * pi / 180));
                    double lan_cal = 0;
                    for (int i = 0; i < (int)Math.Ceiling(360 / landa_interval); i++)
                    {
                        lan_cal = (i * landa_interval * lan_cal2 / 360);
                        if (lan_cal < 360)
                        {
                            
                            landa_s.Add(lan_cal);
                            phi_s.Add(i2 * phi_interval);
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                for (int i = 0; i < landa_s.Count; i++)
                {
                    double X_min = min_distance * Math.Sin(phi_s[i]) * Math.Cos(landa_s[i]);
                    double Y_min = min_distance * Math.Sin(phi_s[i]) * Math.Sin(landa_s[i]);
                    double Z_min = min_distance * Math.Cos(phi_s[i]);
                    double X_max = max_distance * Math.Sin(phi_s[i]) * Math.Cos(landa_s[i]);
                    double Y_max = max_distance * Math.Sin(phi_s[i]) * Math.Sin(landa_s[i]);
                    double Z_max = max_distance * Math.Cos(phi_s[i]);

                    int grid_teta = Convert.ToInt32(Math.Ceiling(landa_s[i] / teta));
                    bool test_dis = false;
                    List<int> i_in_angel = new List<int>();
                    for (int i2=0;i2< polar_triangles[grid_teta - 1, 0].Count; i2++)
                    {
                        int i_in_angel1 = Convert.ToInt32(polar_triangles[grid_teta - 1, 0][i2]);
                        i_in_angel.Add(i_in_angel1);

                        double distance_to_cam = Math.Sqrt(((triangle_coordinates_polar[i_in_angel1, 0] - X_min) * (triangle_coordinates_polar[i_in_angel1, 0] - X_min)) + ((triangle_coordinates_polar[i_in_angel1, 1] - Y_min) * (triangle_coordinates_polar[i_in_angel1, 1] - Y_min)) + ((triangle_coordinates_polar[i_in_angel1, 2] - Z_min) * (triangle_coordinates_polar[i_in_angel1, 2] - Z_min)));
                        if(distance_to_cam< min_distance)
                        {
                            break;
                        }
                       if(i2== polar_triangles[grid_teta - 1, 0].Count - 1)
                        {
                            test_dis = true;
                        }
                    }
                    if (test_dis == true)
                    {
                        double R_min = Math.Sqrt(((X_min) * (X_min)) + ((Y_min) * (Y_min)));
                        double R_max = Math.Sqrt(((X_max) * (X_max)) + ((Y_max) * (Y_max)));
                        for (int ji = 0; ji < polar_triangles[grid_teta - 1, 0].Count; ji++)
                        {
                            bool tri_ray_inter = false;
                            int i_tri = i_in_angel[ji];
                           double R_tri = polar_triangles[grid_teta - 1, 1][ji];
                          
                            if (R_tri < R_min)
                            {
                             
                                for (int ji2 = 0; ji2 < point_triangles[i_tri].Count; ji2++)
                                {
                                    int g1 = i_tri;
                                    int g2 = point_triangles[i_tri][ji2];
                                    int g3 = triangle_int[point_triangles[i_tri][ji2], 0];
                                    int g4 = triangle_int[point_triangles[i_tri][ji2], 1];
                                    int g5 = triangle_int[point_triangles[i_tri][ji2], 2];

                                    double g10 = triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 0], 0];
                                    Vector3 p1 = new Vector3();
                                    p1.X = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 0], 0]);
                                    p1.Y = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 0], 1]);
                                    p1.Z = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 0], 2]);
                                    Vector3 p2 = new Vector3();
                                    p2.X = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 1], 0]);
                                    p2.Y = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 1], 1]);
                                    p2.Z = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 1], 2]);
                                    Vector3 p3 = new Vector3();
                                    p3.X = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 2], 0]);
                                    p3.Y = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 2], 1]);
                                    p3.Z = Convert.ToSingle(triangle_coordinates_polar[triangle_int[point_triangles[i_tri][ji2], 2], 2]);
                                    Vector3 p0 = new Vector3();
                                    p0.X = Convert.ToSingle(X_min);
                                    p0.Y = Convert.ToSingle(Y_min);
                                    p0.Z = Convert.ToSingle(Z_min);
                                    Vector3 p00 = new Vector3();
                                    p00.X = 0;
                                    p00.Y = 0;
                                    p00.Z = 0;
                                    tri_ray_inter = Intersect(p1, p2, p3, p0, p00);
                                    if (tri_ray_inter == true)
                                    {
                                        break;
                                    }
                                }
                             
                            }
                            if (tri_ray_inter == true)
                            {
                           //     ll1l.WriteLine("" + triangle_coordinates[i, 0] + " " + triangle_coordinates[i, 1] + " " + triangle_coordinates[i, 2] + " " + visible_point[i, 3] + " " + visible_point[i, 4] + " " + visible_point[i, 5]);
                                break;
                            }
                            if (ji == polar_triangles[grid_teta - 1, 0].Count - 1)
                            {
                              //  ll1.WriteLine("" + triangle_coordinates[i, 0] + " " + triangle_coordinates[i, 1] + " " + triangle_coordinates[i, 2] + " " + visible_point[i, 3] + " " + visible_point[i, 4] + " " + visible_point[i, 5]);
                              //  Visibility_str.Add("" + i + " " + j2);
                            }
                        }
                    }

                }
                
                //double test = (normal[ij, 0] * R11) + (normal[ij, 1] * R12) + (normal[ij, 2] * R13);
                //double test1 = (normal[ij, 0] * R21) + (normal[ij, 1] * R22) + (normal[ij, 2] * R23);
                //double test2 = (normal[ij, 0] * R31) + (normal[ij, 1] * R32) + (normal[ij, 2] * R33);
            }

                List<string> Visibility_str = new List<string>();
            List<string> Visibility_str0 = new List<string>();
            //   var s10 = Stopwatch.StartNew();
            //  threads_visibility(0, image_name_Eterior.Length, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str0);
            //    int odd = Convert.ToInt32(image_name_Eterior.Length / 8);
            List<string> Visibility_str1 = new List<string>();
            List<string> Visibility_str2 = new List<string>();
            List<string> Visibility_str3 = new List<string>();
            List<string> Visibility_str4 = new List<string>();
            List<string> Visibility_str5 = new List<string>();
            List<string> Visibility_str6 = new List<string>();
            List<string> Visibility_str7 = new List<string>();
            List<string> Visibility_str8 = new List<string>();
            //       threads_visibility_triangulation(0, image_name_Eterior.Length, triangle_int, triangle_coordinates, visible_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str1);
            //Thread t1 = new Thread(() => threads_visibility(0, odd, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str1));
            //Thread t2 = new Thread(() => threads_visibility(odd, odd * 2, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str2));
            //Thread t3 = new Thread(() => threads_visibility(odd * 2, odd * 3, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str3));
            //Thread t4 = new Thread(() => threads_visibility(odd * 3, odd * 4, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str4));
            //Thread t5 = new Thread(() => threads_visibility(odd * 4, odd * 5, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str5));
            //Thread t6 = new Thread(() => threads_visibility(odd * 5, odd * 6, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str6));
            //Thread t7 = new Thread(() => threads_visibility(odd * 6, odd * 7, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str7));
            //Thread t8 = new Thread(() => threads_visibility(odd * 7, image_name_Eterior.Length, tie_point, image_name_Eterior, Exterior_Orientation_Eterior, interior, out Visibility_str8));
            //t1.Start(); t2.Start(); t3.Start(); t4.Start(); t5.Start(); t6.Start(); t7.Start(); t8.Start();
            //t1.Join(); t2.Join(); t3.Join(); t4.Join(); t5.Join(); t6.Join(); t7.Join(); t8.Join();
            for (int i = 0; i < Visibility_str1.Count; i++)
            {
                Visibility_str.Add(Visibility_str1[i]);
            }
            for (int i = 0; i < Visibility_str2.Count; i++)
            {
                Visibility_str.Add(Visibility_str2[i]);
            }
            for (int i = 0; i < Visibility_str3.Count; i++)
            {
                Visibility_str.Add(Visibility_str3[i]);
            }
            for (int i = 0; i < Visibility_str4.Count; i++)
            {
                Visibility_str.Add(Visibility_str4[i]);
            }
            for (int i = 0; i < Visibility_str5.Count; i++)
            {
                Visibility_str.Add(Visibility_str5[i]);
            }
            for (int i = 0; i < Visibility_str6.Count; i++)
            {
                Visibility_str.Add(Visibility_str6[i]);
            }
            for (int i = 0; i < Visibility_str7.Count; i++)
            {
                Visibility_str.Add(Visibility_str7[i]);
            }
            for (int i = 0; i < Visibility_str8.Count; i++)
            {
                Visibility_str.Add(Visibility_str8[i]);
            }
            string[] Visibility_str_array = Visibility_str.ToArray();
            Array.Sort(Visibility_str_array, StringComparer.InvariantCulture);
            Visibility = new int[Visibility_str_array.Length, 2];
            for (int i = 0; i < Visibility_str_array.Length; i++)
            {
                string[] tokens = Visibility_str_array[i].Split(new[] { " " }, StringSplitOptions.None);
                Visibility[i, 0] = Convert.ToInt32(tokens[0]); Visibility[i, 1] = Convert.ToInt32(tokens[1]);
            }
            //s10.Stop();
            //MessageBox.Show("" + s10.ElapsedMilliseconds);

        }

         public void hist(double [,] A,double [] interval, out int[,] B)
        {
            B = new int[10,3];
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < 10; j++)
                {
                    if (A[i, 0] < interval[0]*(j+1) && A[i, 0] > interval[0]*j)
                    {
                        B[j,0] = B[j,0] + 1;
                    }
                }
                for (int j = 0; j < 10; j++)
                {
                    if (A[i, 1] < interval[1] * (j + 1) && A[i, 1] > interval[1] * j)
                    {
                        B[j, 1] = B[j, 1] + 1;
                    }
                }
                for (int j = 0; j < 10; j++)
                {
                    if (A[i, 2] < interval[2] * (j + 1) && A[i, 2] > interval[2] * j)
                    {
                        B[j, 2] = B[j, 2] + 1;
                    }
                }

            }
        } 
        
    }
}
