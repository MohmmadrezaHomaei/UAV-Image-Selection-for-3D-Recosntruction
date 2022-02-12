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
using System.Windows.Forms.DataVisualization.Charting;
namespace theses
{
    public partial class Optimum_selsection : Form
    {
        public Optimum_selsection()
        {
            InitializeComponent();

        }
       public static double[,] point_accuracy_complete;
      

        private void Optimum_selsection_Load(object sender, EventArgs e)
        {

        }

        static int[,] out_hist= new int [10,3] ;
        
        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

       public static List<string> List_of_removed_images = new List<string>();
       public static List<double> List_of_Increasing_Errors = new List<double>();

        private void button2_Click(object sender, EventArgs e)
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
            Hashtable imageobserve = ImportData.imageobserve;
            double[,] grid_coordinate = new double[2,2];
            Hashtable visibility_grid = new Hashtable ();
            Hashtable visibility_tie_com = new Hashtable() ;

            int[,] visibility_tie_com_int = new int[2,2];
            int[,] visibility = new int[2, 2];
            List<int>[] visibility_camera_point =Form1.visibility_camera_point;
            List<int>[] visibility_point_camera =Form1.visibility_point_camera;
        double[,] tri_coor = ImportData.tri_coor;
            double[,] normal = ImportData.normal;
            List<string> im_id2 = Bundle_Adjustment_Form.im_id2;
            double[,] im_ob_sig = Bundle_Adjustment_Form.im_ob_sig;
            List<double> pixel_error = ImportData.pixel_error;
            double[,] tie_com = new double[X_export.Count, 3];
            
            int n = 0;
         
            for (int i = 0; i < X_export.Count; i++)
            {
                tie_com[i, 0] = X_export[i];
                tie_com[i, 1] = Y_export[i];
                tie_com[i, 2] = Z_export[i];
            }
            Bundle_adjustment dd = new Bundle_adjustment();
            
            //   dd.Pre_analysis_image_mark_article(pixel_error,tie_com, visibility_tie_com, image_name, observation, interior, im_ob_sig, out point_accuracy_complete);
            dd.Pre_analysis_image_mark(tri_coor,normal, visibility_point_camera, visibility_camera_point, image_name, observation, interior, im_ob_sig, GPS_weight, out point_accuracy_complete, out List_of_removed_images,out List_of_Increasing_Errors);
            label2.Text = "" + List_of_Increasing_Errors[List_of_Increasing_Errors.Count - 1];
           label3.Text= "" + List_of_removed_images.Count;
            (this.Owner as Form1).dataGridView2.ColumnCount = 1;

            (this.Owner as Form1).dataGridView2.Columns[0].Name = "Image Name";
            for (int i = 0; i < List_of_removed_images.Count; i++)
            {
                (this.Owner as Form1).dataGridView2.Rows.Add(List_of_removed_images[i]);
                (this.Owner as Form1).dataGridView2.Rows[i].HeaderCell.ToolTipText = "Number of Removed Images = " + List_of_Increasing_Errors.Count;
            }
              (this.Owner as Form1).dataGridView2.Visible = true;

            (this.Owner as Form1).chart1.Visible = true;
            for (int i = 0; i < List_of_Increasing_Errors.Count; i++)
            {
                (this.Owner as Form1).chart1.Series[0].Points.AddXY(i + 1, List_of_Increasing_Errors[i]);
            }
            (this.Owner as Form1).chart1.ChartAreas["ChartArea1"].AxisX.MajorGrid.Enabled = false;
            (this.Owner as Form1).chart1.ChartAreas["ChartArea1"].AxisY.MajorGrid.Enabled = false;
            (this.Owner as Form1).chart1.ChartAreas["ChartArea1"].AxisX.Title = "Number of removed images";
            (this.Owner as Form1).chart1.ChartAreas["ChartArea1"].AxisY.Title = "Increasing Errors";
            (this.Owner as Form1).chart1.Series[0].IsValueShownAsLabel = false;
            (this.Owner as Form1).chart1.ChartAreas[0].AxisX.LabelStyle.Enabled = false;
            (this.Owner as Form1).toolTip1.SetToolTip((this.Owner as Form1).chart1, "Maximum Increasing Error = " + List_of_Increasing_Errors[List_of_Increasing_Errors.Count - 1] + "\n" + "Maximum Number of Removed Images = " + List_of_removed_images.Count);
         //   (this.Owner as Form1).toolTip1.SetToolTip((this.Owner as Form1).dataGridView2, "Number of Removed Images = " + List_of_Increasing_Errors.Count);
            (this.Owner as Form1).textBox1.Enabled = true;
            (this.Owner as Form1).textBox2.Enabled = true;
            (this.Owner as Form1).button1.Enabled = true;
        }
    }
}
