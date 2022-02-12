using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Collections;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using Accord.Math;
using System.Diagnostics;
using System.IO;
using System.Xml;
using Accord.Math.Decompositions;



namespace theses
{
    public partial class Form1 : Form
    {
        public string Camera_number;
       public static DataSet ds = new DataSet();
        public Hashtable dsgf = new Hashtable();
        DataColumn dc;
        public Form1()
        {
       
            InitializeComponent();
            
        }

        //public static string[] image_name;
        //public static double[,] observation;
        //public static double[] interior;
        //public static List<double> pixel_error = new List<double>();
        //public static List<double> X_export = new List<double>();
        //public static List<double> Y_export = new List<double>();
        //public static List<double> Z_export = new List<double>();
        //public static double[,] tie_point_sparse;
        //public static double[,] tri_coor;
        //public static int[,] tris;
        public static List<int>[] visibility_camera_point;
        public static List<int>[] visibility_point_camera;
     

       public static double Epsilon = 0.00000001;
          private void button1_Click(object sender, EventArgs e)
        {
          
    
        }
        private void cameraToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }
        private void fileToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }
        private void saveToolStripMenuItem_Click(object sender, EventArgs e)
        {
            SaveFileDialog saveFile = new SaveFileDialog();
            saveFile.Filter = "XML Files|*.xml";
            saveFile.Title = "Save a Xml File";
            saveFile.ShowDialog();
            if (saveFile.FileName != "")
            {
                FileStream fs =
                    (FileStream)saveFile.OpenFile();
                ds.WriteXml(fs);
            }



        }
        private void interiorParametersToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            Camera form = new Camera();
            form.Show();
            DataTable table = ds.Tables["australispara"];
           // MessageBox.Show("" + ds.Tables["australispara"].Rows[0][0] + "    " + ds.Tables["australispara"].Rows[1][0] + "    " + ds.Tables["australispara"].Rows[2][0]);
        }
        private void exteriorParametersToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            DataTable table = ds.Tables["camera"];
           //  MessageBox.Show("" + ds.Tables["camera"].Rows[0][5] + "    " + ds.Tables["camera"].Rows[11][13] + "    " + ds.Tables["camera"].Rows[10][4]);

        }
        private void imageObservationToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            DataTable table = ds.Tables["points"];
              MessageBox.Show("" + ds.Tables["points_id"].Rows[15][0] + "    " + ds.Tables["points_id"].Rows[154][0] + "    " + ds.Tables["points_id"].Rows[45][0] + "    " + ds.Tables["points_id"].Rows[450][0]);

        }
        private void tiePointsCoordinateToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            DataTable table = ds.Tables["XYZ__export"];
             MessageBox.Show("" + ds.Tables["XYZ__export"].Rows[1][0] + "    " + ds.Tables["XYZ__export"].Rows[10][1] + "    " + ds.Tables["XYZ__export"].Rows[8][0] + "    " + ds.Tables["XYZ__export"].Rows[4][1]);

        }
        private void gCPToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            DataTable table = ds.Tables["GCP"];
           //  MessageBox.Show("" + ds.Tables["GCP"].Rows[0][0]+"      "+ ds.Tables["GCP"].Rows[0][1] + "      " + ds.Tables["GCP"].Rows[0][2]);

        }
        private void exteriorAdjusmentToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }
        private void gPSToolStripMenuItem_Click(object sender, EventArgs e)
        {
            DataTable table = ds.Tables["GPS"];
             // MessageBox.Show("" + ds.Tables["GPS"].Rows[254][0]+"      "+ ds.Tables["GPS"].Rows[254][1] + "      " + ds.Tables["GPS"].Rows[306][0] + "      " + ds.Tables["GPS"].Rows[306][1]);

        }
        private void openToolStripMenuItem2_Click(object sender, EventArgs e)
        {
            OpenFileDialog open = new OpenFileDialog();
            open.Filter = "Image Files(*.xml)|*.xml";

            if (open.ShowDialog() == DialogResult.OK)
            {
                ds.ReadXml(open.FileName, XmlReadMode.InferSchema);
            }
            }
       
       
        private void bundleAdjustmentToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Bundle_Adjustment_Form baf = new Bundle_Adjustment_Form();
            baf.Owner = this;
            baf.Show();
            
        }
        private void Form1_Load(object sender, EventArgs e)
        {
        
            toolTip1.SetToolTip(panel1, "Maximum Increasing Error = " + 12.5 + "\n" + "Maximum Number of Removed Images = " + 22);
            //chart1.Series[0].Points.AddXY(2,5);
            //chart1.Series[0].Points.AddXY(2, 10);
            //chart1.Series[0].Points.AddXY(3, 11);
            //chart1.Series[0].Points.AddXY(4, 15);
            //chart1.Series[0].Points.AddXY(5, 16);
            //chart1.Series[0].Points.AddXY(6, 20);
            //chart1.ChartAreas["ChartArea1"].AxisX.MajorGrid.Enabled = false;
            //chart1.ChartAreas["ChartArea1"].AxisY.MajorGrid.Enabled = false;
            //chart1.ChartAreas["ChartArea1"].AxisX.Title = "Number of removed images";
            //chart1.ChartAreas["ChartArea1"].AxisY.Title = "Increasing Errors";
            //chart1.Series[0].IsValueShownAsLabel = false;
            //chart1.ChartAreas[0].AxisX.LabelStyle.Enabled = false;

        }


        private void menuStrip1_ItemClicked(object sender, ToolStripItemClickedEventArgs e)
        {

        }

        private void toolStripMenuItem1_Click(object sender, EventArgs e)
        {
            ImportData op = new ImportData();
            op.Owner = this;
            op.Show();
        }
        public static int image_report_index = new int();

        public void dataGridView1_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {
            string[] image_name = ImportData.image_name;
            double[,] observation = ImportData.observation;
            int row_index = e.RowIndex;
            int colmn_index = e.ColumnIndex;
            image_report_index = row_index;
           // string index_string= dataGridView1.Rows[row_index].Cells[colmn_index].Value.ToString();

            Image_report im_rep2 = new Image_report();
            im_rep2.Owner = this;
            im_rep2.Show();
           im_rep2.label14.Text = Math.Round(observation[row_index, 0],3).ToString();
            im_rep2.label13.Text = Math.Round(observation[row_index, 1],3).ToString();
            im_rep2.label12.Text = Math.Round(observation[row_index, 2],3).ToString();
            im_rep2.label11.Text = Math.Round((observation[row_index, 3]*180/Math.PI),3).ToString();
            im_rep2.label10.Text = Math.Round((observation[row_index, 4] * 180 / Math.PI),3).ToString();
            im_rep2.label9.Text = Math.Round((observation[row_index, 5] * 180 / Math.PI),3).ToString();
            double[,] Cov = Bundle_adjustment.cov_Exterior_orientation_Image;
            if (Cov != null)
            {
                im_rep2.label20.Text = Math.Sqrt(Cov[6 * row_index + 3, 6 * row_index + 3]).ToString();
                im_rep2.label19.Text = Math.Sqrt(Cov[6 * row_index + 4, 6 * row_index + 4]).ToString();
                im_rep2.label18.Text = Math.Sqrt(Cov[6 * row_index + 5, 6 * row_index + 3]).ToString();
                im_rep2.label17.Text = (Math.Sqrt(Cov[6 * row_index + 0, 6 * row_index + 0]) * 180 / Math.PI).ToString();
                im_rep2.label16.Text = (Math.Sqrt(Cov[6 * row_index + 1, 6 * row_index + 1]) * 180 / Math.PI).ToString();
                im_rep2.label15.Text = (Math.Sqrt(Cov[6 * row_index + 2, 6 * row_index + 2]) * 180 / Math.PI).ToString();
            }

        }

        private void visibilityAnalysisToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Grid gg = new Grid();
            int[,] tris = ImportData.tris;
            double[,] tri_coor = ImportData.tri_coor;
            string[] image_name = ImportData.image_name;
            double[,] observation = ImportData.observation;
            double[] interior = ImportData.interior;
            gg.All_tie_points_geometric_visibility_saadat_all_triangulation(tris, tri_coor, image_name, observation, interior, out visibility_point_camera, out visibility_camera_point);
            
            optimumImageSelectionToolStripMenuItem.Enabled = true;
        }
        private void optimumImageSelectionToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Optimum_selsection op = new Optimum_selsection();
            op.Owner = this;
            op.Show();
        }

        private void chart1_Click(object sender, EventArgs e)
        {
           
        }
     
        private void chart1_MouseMove(object sender, MouseEventArgs e)
        {
            
         
        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {
            List<string> List_of_removed_images = new List<string>();
            List<double> List_of_Increasing_Errors = new List<double>();
            List_of_Increasing_Errors = Optimum_selsection.List_of_Increasing_Errors;
            List_of_removed_images = Optimum_selsection.List_of_removed_images;
            if (textBox1.Text != "")
            {
                double desired_acccuracy = Convert.ToDouble(textBox1.Text);
          
                if (desired_acccuracy <= List_of_Increasing_Errors[List_of_Increasing_Errors.Count - 1]) 
                {
                    dataGridView2.Rows.Clear();
                    dataGridView2.Refresh();
                    int maximum_image = 0;
                    for (int i = 0; i < List_of_Increasing_Errors.Count; i++)
                    {
                        if (List_of_Increasing_Errors[i] > desired_acccuracy)
                        {
                            maximum_image = i;
                            break;
                        }
                    }
                    for (int i = 0; i < maximum_image; i++)
                    {
                        dataGridView2.Rows.Add(List_of_removed_images[i]);
                        dataGridView2.Rows[i].HeaderCell.ToolTipText = "Number of Removed Images = " + List_of_Increasing_Errors.Count;
                    }
                //    toolTip1.SetToolTip(dataGridView2, "Number of Removed Images = " + maximum_image);
                }
                else
                {
                   
                    List_of_Increasing_Errors = Optimum_selsection.List_of_Increasing_Errors;
                    List_of_removed_images = Optimum_selsection.List_of_removed_images;
                    dataGridView2.Rows.Clear();
                    dataGridView2.Refresh();

                    for (int i = 0; i < List_of_removed_images.Count; i++)
                    {
                        dataGridView2.Rows.Add(List_of_removed_images[i]);
                        dataGridView2.Rows[i].HeaderCell.ToolTipText = "Number of Removed Images = " + List_of_Increasing_Errors.Count;
                    }
                    //       toolTip1.SetToolTip(dataGridView2, "Number of Removed Images = " + List_of_removed_images.Count);
                }
            }
            else
            {
                
                List_of_Increasing_Errors = Optimum_selsection.List_of_Increasing_Errors;
                List_of_removed_images = Optimum_selsection.List_of_removed_images;
                dataGridView2.Rows.Clear();
                dataGridView2.Refresh();

                for (int i = 0; i < List_of_removed_images.Count; i++)
                {
                    dataGridView2.Rows.Add(List_of_removed_images[i]);
                    dataGridView2.Rows[i].HeaderCell.ToolTipText = "Number of Removed Images = " + List_of_Increasing_Errors.Count;
                }
         //       toolTip1.SetToolTip(dataGridView2, "Number of Removed Images = " + List_of_removed_images.Count);
            }
        }

        private void textBox2_TextChanged(object sender, EventArgs e)
        {
            if (textBox2.Text != "")
            {
                double maximum_image0 = Convert.ToDouble(textBox2.Text);
                int maximum_image = Convert.ToInt32(maximum_image0);
                List<string> List_of_removed_images = new List<string>();
                List<double> List_of_Increasing_Errors = new List<double>();
                List_of_Increasing_Errors = Optimum_selsection.List_of_Increasing_Errors;
                List_of_removed_images = Optimum_selsection.List_of_removed_images;
                dataGridView2.Rows.Clear();
                dataGridView2.Refresh();
                if (maximum_image < List_of_removed_images.Count)
                {
                    for (int i = 0; i < maximum_image; i++)
                    {
                        dataGridView2.Rows.Add(List_of_removed_images[i]);
                        dataGridView2.Rows[i].HeaderCell.ToolTipText = "Number of Removed Images = " + List_of_Increasing_Errors.Count;
                    }
                 //   toolTip1.SetToolTip(dataGridView2, "Number of Removed Images = " + maximum_image);
                }
            }
            else
            {
                List<string> List_of_removed_images = new List<string>();
                List<double> List_of_Increasing_Errors = new List<double>();
                List_of_Increasing_Errors = Optimum_selsection.List_of_Increasing_Errors;
                List_of_removed_images = Optimum_selsection.List_of_removed_images;
                dataGridView2.Rows.Clear();
                dataGridView2.Refresh();

                for (int i = 0; i < List_of_removed_images.Count; i++)
                {
                    dataGridView2.Rows.Add(List_of_removed_images[i]);
                    dataGridView2.Rows[i].HeaderCell.ToolTipText = "Number of Removed Images = " + List_of_Increasing_Errors.Count;
                }
              //  toolTip1.SetToolTip(dataGridView2, "Number of Removed Images = " + List_of_removed_images.Count);
            }
        }

        private void button1_Click_1(object sender, EventArgs e)
        {
            OpenFileDialog fdlg = new OpenFileDialog();
            fdlg.Title = "C# Corner Open File Dialog";
            fdlg.InitialDirectory = @"c:\";
            fdlg.Filter = "Agisoft files (*.psx*)|*.psx*|Agisoft files(*.psx*)|*.psx*"; ;
            fdlg.FilterIndex = 2;
            fdlg.RestoreDirectory = true;
            string agirout="";
            string filename="";
            if (fdlg.ShowDialog() == DialogResult.OK)
            {
                agirout = fdlg.FileName;
                filename = Path.GetFileName(agirout);
                agirout = agirout.Replace("\\"+ filename,"");
            }
            string zipFilePath = agirout + "\\align.files\\0\\chunk.zip";
            string extractionPath = agirout + "\\align.files\\0";
            using (Ionic.Zip.ZipFile zip = Ionic.Zip.ZipFile.Read(zipFilePath))
            {
                Directory.CreateDirectory(extractionPath);
                zip.ExtractAll(extractionPath, Ionic.Zip.ExtractExistingFileAction.OverwriteSilently);
            }
         
            //......................

            XmlDocument xdo = new XmlDocument();
            xdo.Load(agirout + "\\align.files\\0\\doc.xml");
            foreach (XmlNode ee in xdo.GetElementsByTagName("camera"))
            {
                        ee.Attributes["enabled"].Value = "true"; // FirstChild because the inner node is actually the inner text, yeah XmlNode is weird.
            }
            foreach (XmlNode ee in xdo.GetElementsByTagName("camera"))
            {
                for (int i = 0; i < dataGridView2.Rows.Count - 1; i++)
                {
                    if (ee.Attributes["label"].Value.Equals(dataGridView2.Rows[i].Cells[0].Value.ToString()))
                    {
                        ee.Attributes["enabled"].Value = "false"; // FirstChild because the inner node is actually the inner text, yeah XmlNode is weird.
                    }
                }
            }
            File.Delete(agirout + "\\align.files\\0\\doc.xml");
            xdo.Save(agirout + "\\align.files\\0\\doc.xml");
            //......................
    
            File.Delete(agirout + "\\align.files\\0\\chunk.zip");
            Ionic.Zip.ZipFile zip2 = new Ionic.Zip.ZipFile();
            
                zip2.AddFile(agirout + "\\align.files\\0\\doc.xml", "");
                FileInfo fi = new FileInfo(agirout + "\\align.files\\0\\doc.xml");
                System.IO.DirectoryInfo di = new System.IO.DirectoryInfo(agirout + "\\align.files\\0\\doc.xml");

                zip2.Save(agirout + "\\align.files\\0\\chunk.zip");
            File.Delete(agirout + "\\align.files\\0\\doc.xml");
        }
    }
}
