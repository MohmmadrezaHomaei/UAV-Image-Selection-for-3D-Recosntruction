using System;
using System.Collections.Generic;
using System.Collections;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace theses
{
    public partial class ImportData : Form
    {
        public ImportData()
        {
            InitializeComponent();
        }
        public static string[] image_name;
        public static double[,] observation;
        public static double[] interior;
        public static List<double> pixel_error = new List<double>();
        public static List<double> X_export = new List<double>();
        public static List<double> Y_export = new List<double>();
        public static List<double> Z_export = new List<double>();
        public static double[,] tie_point_sparse;
        public static double[,] tri_coor;
        public static int[,] tris;
       
        public static double[,] normal;
        public static double[,] norm;
        public static string[] point_id = new string[1];
        public static double[,] observation2 = new double[1, 1];
        public static string[] point_id_weight;
        public static double[,] GCP_weight = new double[1, 1];
        public static double[,] GPS_weight;
        public static List<string> el_im_ho = new List<string>();
        public static Hashtable GPS = new Hashtable();
        public static List<string> im_id = new List<string>();
        public static Hashtable imageobserve = new Hashtable();
        public static int[] bundle_button = new int[8];
        public static string[] eliminate_images;
        public static string[] inv_observation_array;
        static string sX_Averaged = "", sY_Averaged = "", sZ_Averaged = "", sX_maxd = "", sY_maxd = "", sZ_maxd = "";

        private void label9_Click(object sender, EventArgs e)
        {

        }

        private void label50_Click(object sender, EventArgs e)
        {

        }

        private void ImportData_Load(object sender, EventArgs e)
        {
            
        }
        public static double Epsilon = 0.00000001;
        private void button2_Click(object sender, EventArgs e)
        {


            Read_txt_file read = new Read_txt_file();
            FolderBrowserDialog folder = new FolderBrowserDialog();



            //read.Image_Observation2(Application.StartupPath + @"\\import_parameters\points.txt", out im_id, out imageobserve);

            if (folder.ShowDialog() == DialogResult.OK)
            {
                string folder_path = folder.SelectedPath;
                try
                {
                    read.triangle_mesh(folder_path + @"\\text_tri.txt", out tri_coor, out tris, out normal);
                    this.pictureBox10.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                }
                catch
                {

                }
                //.....................
                try
                {
                    read.pixel_errors(folder_path + @"\\pixel_errors.txt", out pixel_error);

                    this.pictureBox9.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");


                }
                catch
                {

                }
                //.....................
                try
                {
                    read.interior_parameters(folder_path + @"\\australispara.txt", out interior);

                    this.pictureBox1.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    bundle_button[0] = 1;

                }
                catch
                {

                }
                //...............
                try
                {
                    read.OM_PH_KA_XL_YL_ZL(folder_path + @"\\camera.txt", out image_name, out observation);

                    this.pictureBox2.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    label13.Text = image_name.Length.ToString();
                    bundle_button[1] = 1;

                }
                catch
                {

                }
                //..........................
                try
                {
                    read.tie_Point(folder_path + @"\\tie_complete.txt", out X_export, out Y_export, out Z_export, out norm);

                    this.pictureBox3.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    label14.Text = X_export.Count.ToString();
                    bundle_button[2] = 1;
                }
                catch
                {

                }
                //..................
                try
                {
                    read.GCP_Observation(folder_path + @"\\GCP.txt", out point_id, out observation2);
                    this.pictureBox4.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    bundle_button[3] = 1;
                }
                catch
                {

                }
                //..................
                try
                {
                    read.GCP_sigma(folder_path + @"\\GCP_sigma.txt", out point_id_weight, out GCP_weight);

                    this.pictureBox5.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    bundle_button[5] = 1;
                }
                catch
                {

                }
                //..............
                try
                {
                    read.GPS_Observation(folder_path + @"\\GPS.txt", out GPS);

                    this.pictureBox6.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    bundle_button[4] = 1;
                }
                catch
                {

                }
                //..............
                try
                {
                    read.GPS_sigma(folder_path + @"\\GPS_sigma.txt", out GPS_weight);

                    this.pictureBox7.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    bundle_button[6] = 1;
                }
                catch
                {

                }
                //..............
                try
                {
                    read.Image_Observation2(folder_path + @"\\points.txt", out im_id, out inv_observation_array);

                    this.pictureBox8.Image = Image.FromFile(Application.StartupPath + @"\\icons\images.jpg");
                    label15.Text = inv_observation_array.Length.ToString();
                    bundle_button[7] = 1;
                }
                catch
                {

                }
            }
            if (bundle_button[0] == 1 && bundle_button[1] == 1 && bundle_button[2] == 1 && bundle_button[7] == 1 && bundle_button[3] == 1 && bundle_button[4] == 1 && bundle_button[5] == 1 && bundle_button[6] == 1)
            {
                (this.Owner as Form1).bundleAdjustmentToolStripMenuItem.Enabled = true;
                
            }
            if (bundle_button[0] == 1 && bundle_button[1] == 1 && bundle_button[2] == 1 && bundle_button[7] == 1 && bundle_button[4] == 1 && bundle_button[6] == 1)
            {
                (this.Owner as Form1).bundleAdjustmentToolStripMenuItem.Enabled = true;
            }
            if (bundle_button[0] == 1 && bundle_button[1] == 1 && bundle_button[2] == 1 && bundle_button[7] == 1 && bundle_button[3] == 1 && bundle_button[5] == 1 && bundle_button[6] == 1)
            {
                (this.Owner as Form1).bundleAdjustmentToolStripMenuItem.Enabled = true;
            }
            //  image_name, out observation

            (this.Owner as Form1).dataGridView1.ColumnCount = 1;

            (this.Owner as Form1).dataGridView1.Columns[0].Name = "Image Name";
              for(int i=0;i< image_name.Length; i++)
            {
                (this.Owner as Form1).dataGridView1.Rows.Add(image_name[i]);
            }
              (this.Owner as Form1).dataGridView1.Visible = true;

        }
    }
}
