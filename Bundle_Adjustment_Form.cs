using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using Accord.Math;
using System.Collections;
using System.Diagnostics;
using System.IO;
using Accord.Math.Decompositions;
using inv2;
using MathWorks.MATLAB.NET.Utility;
using MathWorks.MATLAB.NET.Arrays;
using System.Runtime.Serialization.Formatters.Binary;
using System.Data.SqlClient;
using save_m;

namespace theses
{
    public partial class Bundle_Adjustment_Form : Form
    {
        public Bundle_Adjustment_Form()
        {
            InitializeComponent();
            label14.Text = om.ToString();
            label13.Text = phi.ToString();
            label12.Text = k.ToString();
            label11.Text = xl.ToString();
            label9.Text = yl.ToString();
            label8.Text = zl.ToString();
        }
        public static  double[,] tie;
        public static List<int> count_camera;
        public static Hashtable hash_tie = new Hashtable();
        public static List<string> im_id2 = new List<string>();
        public static double[,] tie_points_coordinate; public double[,] Exterior_orientation;
        public static double[,] cov_exterior;
        public static double[,] tie_out;
        public static double[,] im_ob_sig;
       public static double[] Average_sigma;
        public static double[,] cam_observe;
        static string om="",phi = "", k = "", xl = "", yl = "", zl = "";

        private void button2_Click(object sender, EventArgs e)
        {
            
            
        }

        private void radioButton2_CheckedChanged(object sender, EventArgs e)
        {
            button1.Enabled = true;
           
        }

        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            button1.Enabled = true;
        }

        private void Bundle_Adjustment_Form_Load(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
         
           
            string[] image_name = ImportData.image_name;
            double[,] observation = ImportData.observation;
            double[] interior = ImportData.interior;
            List<double> X_export = ImportData.X_export;
            List<double> Y_export = ImportData.Y_export;
            List<double> Z_export = ImportData.Z_export;
            double[,] norm = ImportData.norm;
            string[] point_id = ImportData.point_id;
            double[,] observation2 = ImportData.observation2;
            string[] point_id_weight = ImportData.point_id_weight;
            double[,] GCP_weight = ImportData.GCP_weight;
            double[,] GPS_weight = ImportData.GPS_weight;
            Hashtable GPS = ImportData.GPS;
            List<string> im_id = ImportData.im_id;
            List<double> pixel_error = ImportData.pixel_error;
            Hashtable imageobserve = ImportData.imageobserve;
           string[] inv_observation_array = ImportData.inv_observation_array;
        double[,] tie_com = new double[ImportData.X_export.Count, 3];
           
            for (int i = 0; i < X_export.Count; i++)
            {
                tie_com[i, 0] = X_export[i];
                tie_com[i, 1] = Y_export[i];
                tie_com[i, 2] = Z_export[i];
            }
            Bundle_adjustment dd = new Bundle_adjustment();
            
            int atleast_camera_tie_point = 3;



            //   s1.Stop();
            //  MessageBox.Show("sina" + s1.ElapsedMilliseconds);

            try {
                
                    dd.Value_tie_point2(image_name, observation, inv_observation_array, im_id, interior, X_export, Y_export, Z_export, atleast_camera_tie_point,out cam_observe, out im_id2, out tie, out count_camera);
                  //  dd.Bundle_adjustment_fast_article(pixel_error,tie, point_id, observation2, point_id_weight, GCP_weight, image_name, observation, GPS, imageobserve, im_id2, interior, GPS_weight, out tie_points_coordinate, out Exterior_orientation, out cov_exterior, out tie_out, out Average_sigma);
                    dd.Bundle_adjustment_fast( tie, point_id, observation2, point_id_weight, GCP_weight, image_name, observation, GPS, cam_observe, im_id2, interior, GPS_weight, count_camera, out tie_points_coordinate, out Exterior_orientation, out cov_exterior, out tie_out, out Average_sigma);
                    im_ob_sig = cov_exterior;
                    om = Average_sigma[0].ToString(); phi = Average_sigma[1].ToString(); k = Average_sigma[2].ToString();
                    xl = Average_sigma[3].ToString(); yl = Average_sigma[4].ToString(); zl = Average_sigma[5].ToString();
                    label14.Text = om.ToString();
                    label13.Text = phi.ToString();
                    label12.Text = k.ToString();
                    label11.Text = xl.ToString();
                    label9.Text = yl.ToString();
                    label8.Text = zl.ToString();
                    (this.Owner as Form1).visibilityAnalysisToolStripMenuItem.Enabled = true;
          
                
            
            }
            catch
            {
                MessageBox.Show("Import files are wrong!");
            }
        }

        private void label7_Click(object sender, EventArgs e)
        {

        }
    }
}
