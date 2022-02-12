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
namespace theses
{
    class Read_txt_file
    {
        public void OM_PH_KA_XL_YL_ZL(string path, out string[] image_name, out double[,] Exterior_Orientation)
        {
            StreamReader txtfile2 = new StreamReader(path);
            double pis = Math.PI;
            string line;
            int count = 0;
            int count4 = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            image_name = new string[linnum-2];
            Exterior_Orientation = new double[linnum-2, 15];
            while ((line = txtfile.ReadLine()) != null)
            {
                if (count4 < 2)
                {
                    count4 = count4 + 1;
                }
                else {
                    double XL = 0, YL = 0, ZL = 0, OM = 0, PH = 0, KA = 0, r11 = 0, r12 = 0, r13 = 0, r21 = 0, r22 = 0, r23 = 0, r31 = 0, r32 = 0, r33 = 0;
                    string Image_Name = "";
                    string[] tokens = line.Split(new[] { "\t" }, StringSplitOptions.None);
                    string[] tokens2 = new string[16];
                    int o = 0;
                    for (int i = 0; i < tokens.Length; i++)
                    {
                        if (tokens[i] != "")
                        {
                            tokens2[o] = tokens[i];
                            o += 1;
                        }
                    }
         
                                Image_Name = tokens2[0];   
                                XL = Convert.ToDouble(tokens2[1]);
                                YL = Convert.ToDouble(tokens2[2]);
                                ZL = Convert.ToDouble(tokens2[3]);
                                OM = Convert.ToDouble(tokens2[4]);
                                PH = Convert.ToDouble(tokens2[5]);
                                KA = Convert.ToDouble(tokens2[6]);
                                r11 = Convert.ToDouble(tokens2[7]);
                                r12 = Convert.ToDouble(tokens2[8]);
                                r13 = Convert.ToDouble(tokens2[9]);
                                r21 = Convert.ToDouble(tokens2[10]);
                                r22 = Convert.ToDouble(tokens2[11]);
                                r23 = Convert.ToDouble(tokens2[12]);
                                r31 = Convert.ToDouble(tokens2[13]);
                                r32 = Convert.ToDouble(tokens2[14]);
                                r33 = Convert.ToDouble(tokens2[15]);
                
                    image_name[count] = Image_Name;
                    Exterior_Orientation[count, 0] = XL;
                    Exterior_Orientation[count, 1] = YL;
                    Exterior_Orientation[count, 2] = ZL;
                    Exterior_Orientation[count, 3] = OM * pis / 180;
                    Exterior_Orientation[count, 4] = PH * pis / 180;
                    Exterior_Orientation[count, 5] = KA * pis / 180;
                    Exterior_Orientation[count, 6] = r11;
                    Exterior_Orientation[count, 7] = r12;
                    Exterior_Orientation[count, 8] = r13;
                    Exterior_Orientation[count, 9] = r21;
                    Exterior_Orientation[count, 10] = r22;
                    Exterior_Orientation[count, 11] = r23;
                    Exterior_Orientation[count, 12] = r31;
                    Exterior_Orientation[count, 13] = r32;
                    Exterior_Orientation[count, 14] = r33;
                    count = count + 1;
                }
            }

        }
        public void interior_parameters(string path, out double[] interior_Orientation)
        {
            StreamReader txtfile2 = new StreamReader(path);
            interior_Orientation = new double[14];
            string line;
            int count = 0;
            int count4 = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            double f = 0, x0 = 0, y0 = 0, B1 = 0, B2 = 0, K1 = 0, P1 = 0, P2 = 0, K2 = 0, K3 = 0;
            while ((line = txtfile.ReadLine()) != null)
            {
                    count4 = count4 + 1;
                if (count4 ==15)
                {
                    char[] c = line.ToCharArray();
                    c = line.ToCharArray();
                    int count3 = 0;
                    while (c[count3] == ' ')
                    {
                        count3 = count3 + 1;
                    }
                    line = line.Substring(count3 + 2);
                    c = line.ToCharArray();
                    int count32 = 0;
                    while (c[count32] == ' ')
                    {
                        count32 = count32 + 1;
                    }
                    line = line.Substring(count32);
                    c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char1 = new char[len];

                    while (c[count2] != ' ')
                    {
                        image_name_char1[count2] = c[count2];
                        count2 = count2 + 1;
                    }
                    char[] image_name_char = new char[count2];
                    line = line.Substring(count2);
                    for (int jj = 0; jj < count2; jj++)
                    {
                        image_name_char[jj] = image_name_char1[jj];
                    }
                    interior_Orientation[10] = Convert.ToDouble(new string(image_name_char));
                    c = line.ToCharArray();
                    int count332 = 0;
                    while (c[count332] == ' ')
                    {
                        count332 = count332 + 1;
                    }
                    line = line.Substring(count332);
                    c = line.ToCharArray();
                    interior_Orientation[11] = Convert.ToDouble(new string(c))/1000;
                }


                if (count4 == 16)
                {
                    char[] c = line.ToCharArray();
                    c = line.ToCharArray();
                    int count3 = 0;
                    while (c[count3] == ' ')
                    {
                        count3 = count3 + 1;
                    }
                    line = line.Substring(count3 + 2);
                    c = line.ToCharArray();
                    int count32 = 0;
                    while (c[count32] == ' ')
                    {
                        count32 = count32 + 1;
                    }
                    line = line.Substring(count32);
                    c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char1 = new char[len];

                    while (c[count2] != ' ')
                    {
                        image_name_char1[count2] = c[count2];
                        count2 = count2 + 1;
                    }
                    char[] image_name_char = new char[count2];
                    line = line.Substring(count2);
                    for (int jj = 0; jj < count2; jj++)
                    {
                        image_name_char[jj] = image_name_char1[jj];
                    }
                    interior_Orientation[12] = Convert.ToDouble(new string(image_name_char));
                    c = line.ToCharArray();
                    int count332 = 0;
                    while (c[count332] == ' ')
                    {
                        count332 = count332 + 1;
                    }
                    line = line.Substring(count332);
                    c = line.ToCharArray();
                    interior_Orientation[13] = Convert.ToDouble(new string(c))/1000;
                }

                if (count4 > 19 && count4 <30)
                {
                    char[] c = line.ToCharArray();
                    c = line.ToCharArray();
                    int count3 = 0;
                    while (c[count3] == ' ')
                    {
                        count3 = count3 + 1;
                    }
                    line = line.Substring(count3+2);
                    c = line.ToCharArray();
                    int count32 = 0;
                    while (c[count32] == ' ')
                    {
                        count32 = count32 + 1;
                    }
                    line = line.Substring(count32);
                     c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char1 = new char[len];

                    while (c[count2] != ' ')
                    {
                        image_name_char1[count2] = c[count2];
                        count2 = count2 + 1;
                    }
                    char[] image_name_char = new char[count2];
                    for (int jj = 0; jj < count2; jj++)
                    {
                        image_name_char[jj] = image_name_char1[jj];
                    }
                    interior_Orientation[count4-20] = Convert.ToDouble(new string(image_name_char));
                    
                }
                
            }
            interior_Orientation[0] = interior_Orientation[0] / 1000 ;
            interior_Orientation[1] = interior_Orientation[1] / 1000;
            interior_Orientation[2] = interior_Orientation[2] / 1000;
            interior_Orientation[3] = interior_Orientation[3] / 1000;
            interior_Orientation[4] = interior_Orientation[4] / 1000;
            interior_Orientation[5] = interior_Orientation[5] / 1000;
            interior_Orientation[6] = interior_Orientation[6] / 1000;
            interior_Orientation[7] = interior_Orientation[7] / 1000;
            interior_Orientation[8] = interior_Orientation[8] / 1000;
            interior_Orientation[9] = interior_Orientation[9] / 1000;


        }
        public void pixel_errors(string path, out List<double> pixel_error)
        {
            StreamReader txtfile2 = new StreamReader(path);
            pixel_error = new List<double>();
            string line;

            StreamReader txtfile = new StreamReader(path);
         
            while ((line = txtfile.ReadLine()) != null)
            {
                pixel_error.Add(Convert.ToDouble(line));
            }
        }
        public void Eliminate_images_hosseini_naveh(string path, out List<string> Eliminate_images)
        {
            StreamReader txtfile2 = new StreamReader(path);
            Eliminate_images = new List<string>();
            string line;

            StreamReader txtfile = new StreamReader(path);

            while ((line = txtfile.ReadLine()) != null)
            {
                Eliminate_images.Add(Convert.ToString(line));
            }
        }
        public DataTable XYZ(string path)
        {
            StreamReader txtfile = new StreamReader(path);
            //..............................................................................................................
            DataTable txt_file_table = new DataTable("XYZ");
            DataColumn txt_file_table_colum;//dc
            DataRow txt_file_table_row;//dr

            // create Name column.
            txt_file_table_colum = new DataColumn();
            txt_file_table_colum.DataType = System.Type.GetType("System.Double");
            txt_file_table_colum.ColumnName = "X";
            txt_file_table_colum.Caption = "X";
            txt_file_table_colum.AutoIncrement = false;
            txt_file_table_colum.ReadOnly = false;
            txt_file_table_colum.Unique = false;

            // Add Name Column to the table.
            txt_file_table.Columns.Add(txt_file_table_colum);
            // Create Address column.
            txt_file_table_colum = new DataColumn();
            txt_file_table_colum.DataType = System.Type.GetType("System.Double");
            txt_file_table_colum.ColumnName = "Y";
            txt_file_table_colum.Caption = "Y";
            txt_file_table_colum.AutoIncrement = false;
            txt_file_table_colum.ReadOnly = false;
            txt_file_table_colum.Unique = false;

            // Add Name Column to the table.
            txt_file_table.Columns.Add(txt_file_table_colum);
            // Create Address column.
            txt_file_table_colum = new DataColumn();
            txt_file_table_colum.DataType = System.Type.GetType("System.Double");
            txt_file_table_colum.ColumnName = "Z";
            txt_file_table_colum.Caption = "Z";
            txt_file_table_colum.AutoIncrement = false;
            txt_file_table_colum.ReadOnly = false;
            txt_file_table_colum.Unique = false;

            // Add Name Column to the table.
            txt_file_table.Columns.Add(txt_file_table_colum);

            DataColumn[] PrimaryKeyColumns = new DataColumn[1];
            PrimaryKeyColumns[0] = txt_file_table.Columns["X"];
            txt_file_table.PrimaryKey = PrimaryKeyColumns;
            //..............................................................................................................
            string line;
            int count = 0;
            while ((line = txtfile.ReadLine()) != null)
            {
                txt_file_table_row = txt_file_table.NewRow();
                count = count + 1;
                for (int i = 1; i < 4; i++)
                {
                    char[] c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char = new char[len];
                    while (c[count2] != ' ')
                    {
                        image_name_char[count2] = c[count2];

                        count2 = count2 + 1;
                    }

                    line = line.Substring(count2 + 1);
                    switch (i)
                    {
                        case 1:
                            txt_file_table_row["X"] = new string(image_name_char);
                            break;
                        case 2:
                            txt_file_table_row["Y"] = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 3:
                            txt_file_table_row["Z"] = Convert.ToDouble(new string(image_name_char));
                            break;
                    }
                }
                try
                {
                    txt_file_table.Rows.Add(txt_file_table_row);
                }
                catch
                {
                }

            }

            return txt_file_table;

        }
       // public void Image_Observation2(string path, out List<string> point_id_list, out List<double> list_x, out List<double> list_y, out List<string> list_image, out List<int> list_image_num, out List<string> list_point)
          public void Image_Observation2(string path, out List<string> point_id_list, out string [] inv_observation_array)
        {
            var s1 = System.Diagnostics.Stopwatch.StartNew();
            StreamReader txtfile = new StreamReader(path);

            string line;
            string Image_Name = "";
            string point_id = "";
            Hashtable point_in = new Hashtable();
            int num = 0;
            point_id_list = new List<string>();
            List<double> list_x = new List<double>();
            List<double> list_y = new List<double>();
            StreamReader txtfile2 = new StreamReader(path);
            List<string> list_image = new List<string>();
            List<int> list_image_num = new List<int>();
            List<string>  list_point = new List<string>();
            int image_num_count = new int();
            List<string> inv_observation = new List<string>();
            int image_num = -1;
            string image_test = "";
            while ((line = txtfile2.ReadLine()) != null)
            {
                // create Name column.
                double x = 0, y = 0;
                string[] tokens = line.Split(new[] { " " }, StringSplitOptions.None);
                string[] tokens2 = new string[6];
                int o = 0;
                for (int i = 0; i < tokens.Length; i++)
                {
                    if (tokens[i] != "")
                    {
                        tokens2[o] = tokens[i];
                        o += 1;
                    }
                }
                if (image_test != tokens2[0]) image_num += 1;
                Image_Name = tokens2[0];
                image_test = Image_Name;
                point_id = tokens2[1];
           
                inv_observation.Add("" + point_id + " " + image_num + "  " + tokens2[2] + " " + tokens2[3]);



                //  list_point.Add(Image_Name + " " + point_id + " " + "y" + y);
                // Imgae_obesrvation.Add(Image_Name + " " + point_id + " " + "x", x / 1000);
                //  Imgae_obesrvation.Add(Image_Name + " " + point_id + " " + "y", y / 1000);


                if (!point_in.ContainsKey(point_id))
                {

                    point_id_list.Add(point_id);

                }
                point_in[point_id] = x / 1000;
            }
            list_image_num.Add(image_num_count);
            list_image_num.Remove(0);
           
            inv_observation_array = inv_observation.ToArray();
            Array.Sort(inv_observation_array, StringComparer.InvariantCulture);
            s1.Stop();
            System.IO.StreamWriter lS = new System.IO.StreamWriter(@"E:\test\invob.txt");
            for (int i = 0; i < inv_observation_array.Length; i++)
            {
                
                lS.WriteLine(inv_observation_array[i]);

            }
            lS.Close();

            // MessageBox.Show("sina" + s1.ElapsedMilliseconds);
        }
        public void Image_Observation(string path, out List<string> point_id_list, out Hashtable Imgae_obesrvation)
        {
            StreamReader txtfile = new StreamReader(path);
            Imgae_obesrvation = new Hashtable();
            string line;
            string Image_Name = "";
            string point_id = "";
            Hashtable point_in = new Hashtable();
            int num = 0;
            point_id_list = new List<string>();
            List<double> list_x = new List<double>();
            List<double> list_y = new List<double>();
            StreamReader txtfile2 = new StreamReader(path);
            List<string> list_all = new List<string>();

            while ((line = txtfile2.ReadLine()) != null)
            {
                // create Name column.
                double x = 0, y = 0;
                string[] tokens = line.Split(new[] { " " }, StringSplitOptions.None);
                string[] tokens2 = new string[6];
                int o = 0;
                for (int i = 0; i < tokens.Length; i++)
                {
                    if (tokens[i] != "")
                    {
                        tokens2[o] = tokens[i];
                        o += 1;
                    }
                }

                Image_Name = tokens2[0];

                point_id = tokens2[1];

                x = Convert.ToDouble(tokens2[2]);

                y = Convert.ToDouble(tokens2[3]);

              //  list_all.Add(Image_Name + " " + point_id + " " + "x" + x);
               // list_all.Add(Image_Name + " " + point_id + " " + "y" + y);
                Imgae_obesrvation.Add(Image_Name + " " + point_id + " " + "x", x / 1000);
                Imgae_obesrvation.Add(Image_Name + " " + point_id + " " + "y", y / 1000);
                if (!point_in.ContainsKey(point_id))
                {

                    point_id_list.Add(point_id);

                }
                point_in[point_id] = x / 1000;
            }
        }
        public void triangle_mesh(string path, out double[,] triangle_coordinates,out int[,] triangle_int,out double [,] triangles_normal)
        {
            var s1 = System.Diagnostics.Stopwatch.StartNew();
            StreamReader txtfile = new StreamReader(path);

            string line;
            int num = 0;

            StreamReader txtfile2 = new StreamReader(path);
            bool test = false;
            bool test2 = false;
            bool test3 = false;
            List<string> coordinates = new List<string>();
            List<string> triangles = new List<string>();
            List<string> normal = new List<string>();
            int num2 = 0;
            int num3 = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                num++;
                num2++;
                num3++;
                string[] tokens = line.Split(new[] { "  " }, StringSplitOptions.None);
                if (tokens.Length>4)
                {
                    if (tokens[4] == "coord Coordinate {")
                    {
                        test = true;
                        num = 0;
                 
                    }
                }
                if (test == true && num>2)
                {
                    string[] coordinate = tokens[6].Split(new[] { "," }, StringSplitOptions.None);

                    if (coordinate.Length == 5)
                    {
                        coordinates.Add(coordinate[0]);
                        coordinates.Add(coordinate[1]);
                        coordinates.Add(coordinate[2]);
                        coordinates.Add(coordinate[3]);
                    }
                    else
                    {
                        for(int ij=0;ij< coordinate.Length; ij++)
                        {
                            coordinates.Add(coordinate[ij]);
                        }
                        test = false;
                    }
                }
                if (tokens.Length == 5)
                {
                    if (tokens[4] == "normal Normal {")
                    {
                        test3 = true;
                        num3 = 0;
                    }

                }
                if (test3 == true && num3 > 2)
                {
                    string[] norm = tokens[6].Split(new[] { ", " }, StringSplitOptions.None);
                    if (norm.Length == 5)
                    {
                        normal.Add(norm[0]);
                        normal.Add(norm[1]);
                        normal.Add(norm[2]);
                        normal.Add(norm[3]);
                    }
                    else
                    {
                        for (int ij = 0; ij < norm.Length; ij++)
                        {
                            normal.Add(norm[ij]);
                        }
                      //  string[] triangle2 = norm[norm.Length - 1].Split(new[] { ", -1" }, StringSplitOptions.None);
                     //   normal.Add(triangle2[0]);
                        test3 = false;

                    }

                }
                if (tokens.Length == 5)
                {
                    if (tokens[4] == "coordIndex")
                    {
                        test2 = true;
                        num2 = 0;
                    }
                }
                if (test2 == true && num2 > 1)
                {
                    string[] triangle= tokens[5].Split(new[] { ", -1, " }, StringSplitOptions.None);
                    if (triangle.Length==7)
                    {
                        triangles.Add(triangle[0]);
                        triangles.Add(triangle[1]);
                        triangles.Add(triangle[2]);
                        triangles.Add(triangle[3]);
                        triangles.Add(triangle[4]);
                        triangles.Add(triangle[5]);
                    }
                    else
                    {
                        for (int ij = 0; ij < triangle.Length-1; ij++)
                        {
                            triangles.Add(triangle[ij]);
                        }
                        string[] triangle2 = triangle[triangle.Length-1].Split(new[] { ", -1" }, StringSplitOptions.None);
                        triangles.Add(triangle2[0]);
                        test2 = false;

                    }

                }
                }
          //  System.IO.StreamWriter ll1l = new System.IO.StreamWriter(@"E:\test\tri_coor.txt");
            triangle_coordinates = new double[coordinates.Count, 3];
            for(int i=0; i< coordinates.Count; i++)
            {
                string[] trian = coordinates[i].Split(new[] { " " }, StringSplitOptions.None);
                if (trian.Length == 3)
                {
                    triangle_coordinates[i, 0] = Convert.ToDouble( trian[0]);
                    triangle_coordinates[i, 1] = Convert.ToDouble(trian[1]);
                    triangle_coordinates[i, 2] = Convert.ToDouble(trian[2]);
                 //   ll1l.WriteLine("" + triangle_coordinates[i, 0] + " " + triangle_coordinates[i, 1] + " " + triangle_coordinates[i, 2]);
                }
                else
                {
                    triangle_coordinates[i, 0] = Convert.ToDouble(trian[1]);
                    triangle_coordinates[i, 1] = Convert.ToDouble(trian[2]);
                    triangle_coordinates[i, 2] = Convert.ToDouble(trian[3]);
                 //   ll1l.WriteLine("" + triangle_coordinates[i, 0] + " " + triangle_coordinates[i, 1] + " " + triangle_coordinates[i, 2]);
                }
               
            }
          //  ll1l.Close();
         //   System.IO.StreamWriter ll1 = new System.IO.StreamWriter(@"E:\test\tris.txt");
            triangle_int = new int[triangles.Count, 3];
            for (int i = 0; i < triangles.Count; i++)
            {
                string[] trian2 = triangles[i].Split(new[] { ", " }, StringSplitOptions.None);
                triangle_int[i, 0] = Convert.ToInt32(trian2[0]);
                triangle_int[i, 1] = Convert.ToInt32(trian2[1]);
                triangle_int[i, 2] = Convert.ToInt32(trian2[2]);
           //     ll1.WriteLine("" + triangle_int[i, 0] + " " + triangle_int[i, 1] + " " + triangle_int[i, 2]);
            }
            //   ll1.Close();

            triangles_normal = new double[normal.Count, 3];
            for (int i = 0; i < normal.Count; i++)
            {
                string[] trian2 = normal[i].Split(new[] { " " }, StringSplitOptions.None);
                triangles_normal[i, 0] = Convert.ToDouble(trian2[0]);
                triangles_normal[i, 1] = Convert.ToDouble(trian2[1]);
                triangles_normal[i, 2] = Convert.ToDouble(trian2[2]);
                //     ll1.WriteLine("" + triangle_int[i, 0] + " " + triangle_int[i, 1] + " " + triangle_int[i, 2]);
            }




            // MessageBox.Show("sina" + s1.ElapsedMilliseconds);
        }


        public void GPS_Observation(string path,out Hashtable GPS_observation)
        {
            StreamReader txtfile2 = new StreamReader(path);

            string line;
            int count4 = 0;
            int linnum = 0;
            int count5 = 0;

           GPS_observation = new Hashtable();
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);


            while ((line = txtfile.ReadLine()) != null && count5 < linnum - 3)
            {
                if (count4 < 2)
                {
                    count4 = count4 + 1;
                }
                else {
                    double X = 0, Y = 0, Z = 0;
                    string Image_Name = "";
                    string[] tokens = line.Split(new[] { "," }, StringSplitOptions.None);
                    string[] tokens2 = new string[4];
                    int o = 0;
                    for (int i = 0; i < 7; i++)
                    {
                        if (tokens[i] != "")
                        {
                            tokens2[o] = tokens[i];
                            o += 1;
                        }
                    }
             
                                Image_Name = tokens2[0];
                                
                                X = Convert.ToDouble(tokens2[1]);

                                Y = Convert.ToDouble(tokens2[2]);
                       
                                Z = Convert.ToDouble(tokens2[3]);

                    GPS_observation[Image_Name + " " + "X"] = X;
                    GPS_observation[Image_Name + " " + "Y"] = Y;
                    GPS_observation[Image_Name + " " + "Z"] = Z;
                    count5 = count5 + 1;
                }
            }
        }
        public void GCP_Observation(string path, out string[] point_name, out double[,] GCP_observayion)
        {
            StreamReader txtfile2 = new StreamReader(path);

            string line;
            int count = 0;
            int count4 = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            point_name = new string[linnum - 3];
            GCP_observayion = new double[linnum - 3, 3];
            while ((line = txtfile.ReadLine()) != null && count < linnum - 3)
            {
                if (count4 < 2)
                {
                    count4 = count4 + 1;
                }
                else {
                    double X = 0, Y = 0, Z = 0;
                    string Image_Name = "";
                    for (int i = 1; i < 5; i++)
                    {
                        char[] c = line.ToCharArray();
                        int len = c.Length;
                        int count2 = 0;
                        char[] image_name_char1 = new char[len];

                        while (c[count2] != ',')
                        {
                            image_name_char1[count2] = c[count2];
                            count2 = count2 + 1;
                        }
                        char[] image_name_char = new char[count2];
                        for (int jj = 0; jj < count2; jj++)
                        {
                            image_name_char[jj] = image_name_char1[jj];
                        }
                        line = line.Substring(count2 + 1);
                        c = line.ToCharArray();
                        switch (i)
                        {
                            case 1:
                                Image_Name = new string(image_name_char);
                                break;
                            case 2:
                                if (new string(image_name_char) != "")
                                {
                                    X = Convert.ToDouble(new string(image_name_char));
                                }
                                break;
                            case 3:
                                if (new string(image_name_char) != "")
                                {
                                    Y = Convert.ToDouble(new string(image_name_char));
                                }
                                break;
                            case 4:
                                if (new string(image_name_char) != "")
                                {
                                    Z = Convert.ToDouble(new string(image_name_char));
                                }
                                break;
                        }
                    }
                    if (X != 0)
                    {
                        point_name[count] = Image_Name.Replace(' ', '_');
                        GCP_observayion[count, 0] = X;
                        GCP_observayion[count, 1] = Y;
                        GCP_observayion[count, 2] = Z;
                        count = count + 1;
                    }
                }
            }
        }
        public void GCP_sigma(string path, out string[] point_name, out double[,] GCP_sigma)
        {
            StreamReader txtfile2 = new StreamReader(path);

            string line;
            int count = 0;
            int count4 = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            point_name = new string[linnum - 1];
            GCP_sigma = new double[linnum - 1, 3];
            while ((line = txtfile.ReadLine()) != null)
            {
                if (count4 < 1)
                {
                    count4 = count4 + 1;
                }
                else {
                    double X = 0, Y = 0, Z = 0;
                    string Image_Name = "";
                    for (int i = 1; i < 5; i++)
                    {
                        char[] c = line.ToCharArray();
                        int len = c.Length;
                        int count2 = 0;
                        char[] image_name_char1 = new char[len];
                        if (i < 4)
                        {
                            while (c[count2] != ' ')
                            {
                                image_name_char1[count2] = c[count2];
                                count2 = count2 + 1;
                            }
                            char[] image_name_char = new char[count2];
                            for (int jj = 0; jj < count2; jj++)
                            {
                                image_name_char[jj] = image_name_char1[jj];
                            }
                            line = line.Substring(count2 + 1);
                            c = line.ToCharArray();

                            switch (i)
                            {
                                case 1:
                                    Image_Name = new string(image_name_char);
                                    break;
                                case 2:
                                    if (new string(image_name_char) != "")
                                    {
                                        X = Convert.ToDouble(new string(image_name_char));
                                    }
                                    break;
                                case 3:
                                    if (new string(image_name_char) != "")
                                    {
                                        Y = Convert.ToDouble(new string(image_name_char));
                                    }
                                    break;
                             
                            }
                        }
                        else
                        {
                            Z = Convert.ToDouble(new string(c));
                        }
                    }
                    if (X != 0)
                    {
                        point_name[count] = Image_Name;
                        GCP_sigma[count, 0] = X;
                        GCP_sigma[count, 1] = Y;
                        GCP_sigma[count, 2] = Z;
                        count = count + 1;
                    }
                }
            }
        }
        public void GPS_sigma(string path, out double[,] GPS_sigma)
        {
            StreamReader txtfile2 = new StreamReader(path);

            string line;
            int count = 0;
            int count4 = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
           
            GPS_sigma = new double[linnum - 1, 6];
            while ((line = txtfile.ReadLine()) != null)
            {
                if (count4 < 1)
                {
                    count4 = count4 + 1;
                }
                else {
                    double sigom = 0,sigphi = 0,sigk = 0, sigX = 0, sigY = 0, sigZ = 0;
                   
                    for (int i = 1; i < 7; i++)
                    {
                        char[] c = line.ToCharArray();
                        int len = c.Length;
                        int count2 = 0;
                        char[] image_name_char1 = new char[len];
                        if (i < 6)
                        {
                            while (c[count2] != ' ')
                            {
                                image_name_char1[count2] = c[count2];
                                count2 = count2 + 1;
                            }
                            char[] image_name_char = new char[count2];
                            for (int jj = 0; jj < count2; jj++)
                            {
                                image_name_char[jj] = image_name_char1[jj];
                            }
                            line = line.Substring(count2 + 1);
                            c = line.ToCharArray();

                            switch (i)
                            {
                                case 1:
                                    sigom = Convert.ToDouble(new string(image_name_char));
                                    break;
                                case 2:
                                    if (new string(image_name_char) != "")
                                    {
                                        sigphi = Convert.ToDouble(new string(image_name_char));
                                    }
                                    break;
                                case 3:
                                    if (new string(image_name_char) != "")
                                    {
                                        sigk = Convert.ToDouble(new string(image_name_char));
                                    }
                                    break;
                                case 4:
                                    if (new string(image_name_char) != "")
                                    {
                                        sigX = Convert.ToDouble(new string(image_name_char));
                                    }
                                    break;
                                case 5:
                                    if (new string(image_name_char) != "")
                                    {
                                        sigY = Convert.ToDouble(new string(image_name_char));
                                    }
                                    break;

                            }
                        }
                        else
                        {
                            sigZ = Convert.ToDouble(new string(c));
                        }
                    }
                    if (sigX != 0)
                    {
                       
                        GPS_sigma[count, 0] = sigom;
                        GPS_sigma[count, 1] = sigphi;
                        GPS_sigma[count, 2] = sigk;
                        GPS_sigma[count, 3] = sigX;
                        GPS_sigma[count, 4] = sigY;
                        GPS_sigma[count, 5] = sigZ;
                        count = count + 1;
                    }
                }
            }
        }
        public void tie_Point(string path, out List<double> X, out List<double> Y, out List<double> Z,out double[,] norm)
        {
            StreamReader txtfile2 = new StreamReader(path);
            // Coulmn_name[0] = "X"; Coulmn_name[1] = "Y"; Coulmn_name[2] = "Z"; Coulmn_name[3] = "sigmaX"; Coulmn_name[4] = "sigmaY"; Coulmn_name[5] = "sigmaZ";
            string line;
            X = new List<double>();
            Y = new List<double>();
            Z = new List<double>();
            int count = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            norm = new double[linnum, 3];
            while ((line = txtfile.ReadLine()) != null)
            {
                double X1 = 0, Y1 = 0, Z1 = 0, n1=0, n2 = 0 , n3 = 0;
                string Image_Name = "";
                for (int i = 1; i < 10; i++)
                {
                    char[] c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char1 = new char[len];

                    while (c[count2] != ' ' && i < 9)
                    {
                        image_name_char1[count2] = c[count2];
                        count2 = count2 + 1;
                    }
                    char[] image_name_char = new char[count2];
                    for (int jj = 0; jj < count2; jj++)
                    {
                        image_name_char[jj] = image_name_char1[jj];
                    }
                    line = line.Substring(count2);
                    c = line.ToCharArray();
                    int count3 = 0;
                    while (c[count3] == ' ')
                    {
                        count3 = count3 + 1;
                    }
                    line = line.Substring(count3);
                    switch (i)
                    {
                        case 1:
                            X1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 2:
                            Y1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 3:
                            Z1 = Convert.ToDouble(new string(image_name_char));
                            break;

                    }
                }
                X.Add(X1);
                Y.Add(Y1);
                Z.Add(Z1);

                count = count + 1;

            }
        }
        public void tie_Point_sparse(string path, out double [,]tie_point_sparse)
        {
            StreamReader txtfile2 = new StreamReader(path);
            // Coulmn_name[0] = "X"; Coulmn_name[1] = "Y"; Coulmn_name[2] = "Z"; Coulmn_name[3] = "sigmaX"; Coulmn_name[4] = "sigmaY"; Coulmn_name[5] = "sigmaZ";
            string line;
       
            int count = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            tie_point_sparse = new double[linnum, 6];
            while ((line = txtfile.ReadLine()) != null)
            {
                double X1 = 0, Y1 = 0, Z1 = 0,c1 = 0, c2 = 0, c3 = 0;

                for (int i = 1; i < 10; i++)
                {
                    char[] c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char1 = new char[len];

                    while (c[count2] != ' ' && i < 9)
                    {
                        image_name_char1[count2] = c[count2];
                        count2 = count2 + 1;
                    }
                    char[] image_name_char = new char[count2];
                    for (int jj = 0; jj < count2; jj++)
                    {
                        image_name_char[jj] = image_name_char1[jj];
                    }
                    line = line.Substring(count2);
                    c = line.ToCharArray();
                    int count3 = 0;
                    while (c[count3] == ' ')
                    {
                        count3 = count3 + 1;
                    }
                    line = line.Substring(count3);
                    switch (i)
                    {
                        case 1:
                            X1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 2:
                            Y1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 3:
                            Z1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 4:
                            c1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 5:
                            c2 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 6:
                            c3 = Convert.ToDouble(new string(image_name_char));
                            break;

                    }
                }
                tie_point_sparse[count, 0] = X1;
                tie_point_sparse[count, 1] = Y1;
                tie_point_sparse[count, 2] = Z1;
                tie_point_sparse[count, 3] = c1;
                tie_point_sparse[count, 4] = c2;
                tie_point_sparse[count, 5] = c3;

                count = count + 1;

            }
        }
        public void tie_Point_sparse2(string path, out double[,] tie_point_sparse)
        {
            StreamReader txtfile2 = new StreamReader(path);
            // Coulmn_name[0] = "X"; Coulmn_name[1] = "Y"; Coulmn_name[2] = "Z"; Coulmn_name[3] = "sigmaX"; Coulmn_name[4] = "sigmaY"; Coulmn_name[5] = "sigmaZ";
            string line;

            int count = 0;
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            tie_point_sparse = new double[linnum, 6];
            while ((line = txtfile.ReadLine()) != null)
            {
                double X1 = 0, Y1 = 0, Z1 = 0, c1 = 0, c2 = 0, c3 = 0;

                for (int i = 1; i < 4; i++)
                {
                    char[] c = line.ToCharArray();
                    int len = c.Length;
                    int count2 = 0;
                    char[] image_name_char1 = new char[len];

                    while (c[count2] != ' ' && i < 3)
                    {
                        image_name_char1[count2] = c[count2];
                        count2 = count2 + 1;
                    }
                    char[] image_name_char = new char[count2];
                    for (int jj = 0; jj < count2; jj++)
                    {
                        image_name_char[jj] = image_name_char1[jj];
                    }
                    line = line.Substring(count2);
                    c = line.ToCharArray();
                    int count3 = 0;
                    while (c[count3] == ' ')
                    {
                        count3 = count3 + 1;
                    }
                    line = line.Substring(count3);
                    switch (i)
                    {
                        case 1:
                            X1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 2:
                            Y1 = Convert.ToDouble(new string(image_name_char));
                            break;
                        case 3:
                            Z1 = Convert.ToDouble(line);
                            break;

                    }
                }
                tie_point_sparse[count, 0] = X1;
                tie_point_sparse[count, 1] = Y1;
                tie_point_sparse[count, 2] = Z1;

                count = count + 1;

            }
        }

        public void OM_PH_KA_XL_YL_ZL_sigma(string path, out string[] image_name, out double[,] Exterior_Orientation)
        {
            StreamReader txtfile2 = new StreamReader(path);

            string line;
            int count = 0;

            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            image_name = new string[linnum];
            Exterior_Orientation = new double[linnum , 6];
            while ((line = txtfile.ReadLine()) != null)
            {

                double XL = 0, YL = 0, ZL = 0, OM = 0, PH = 0, KA = 0;
                    string Image_Name = "";
                    for (int i = 1; i < 8; i++)
                    {
                        char[] c = line.ToCharArray();
                        int len = c.Length;
                        int count2 = 0;
                        char[] image_name_char1 = new char[len];

                        while (c[count2] != ' ' && i != 7)
                        {
                            image_name_char1[count2] = c[count2];
                            count2 = count2 + 1;
                        }
                        char[] image_name_char = new char[count2];
                        for (int jj = 0; jj < count2; jj++)
                        {
                            image_name_char[jj] = image_name_char1[jj];
                        }
                        line = line.Substring(count2);
                        c = line.ToCharArray();
                        int count3 = 0;
                        while (c[count3] == ' ')
                        {
                            count3 = count3 + 1;
                        }
                        line = line.Substring(count3);

                        switch (i)
                        {
                            case 1:
                                Image_Name = new string(image_name_char);
                        //    MessageBox.Show("" + Image_Name);
                                break;
                            case 2:
                           
                            XL = Convert.ToDouble(new string(image_name_char));
                         //   MessageBox.Show("" + XL);
                            break;
                            case 3:
                                YL = Convert.ToDouble(new string(image_name_char));
                          //  MessageBox.Show("" + YL);
                            break;
                            case 4:
                                ZL = Convert.ToDouble(new string(image_name_char));
                         //   MessageBox.Show("" + ZL);

                            break;
                            case 5:
                                OM = Convert.ToDouble(new string(image_name_char));
                           /// MessageBox.Show("" + OM);
                            break;
                            case 6:

                                PH = Convert.ToDouble(new string(image_name_char));
                           // MessageBox.Show("" + PH);
                            break;
                            case 7:
                          
                            KA = Convert.ToDouble(line);
                           // MessageBox.Show("" + KA);
                            break;

                        }
                    }
                    image_name[count] = Image_Name;
                    Exterior_Orientation[count, 0] = XL;
                    Exterior_Orientation[count, 1] = YL;
                    Exterior_Orientation[count, 2] = ZL;
                    Exterior_Orientation[count, 3] = OM  ;
                    Exterior_Orientation[count, 4] = PH  ;
                    Exterior_Orientation[count, 5] = KA  ;

                    count = count + 1;
                
            }

        }
        public void read_eliminate_images(string path, out string[] image_name)
        {
            StreamReader txtfile2 = new StreamReader(path);

            string line;
           
            int linnum = 0;
            while ((line = txtfile2.ReadLine()) != null)
            {
                linnum = linnum + 1;
            }
            StreamReader txtfile = new StreamReader(path);
            image_name = new string[linnum ];
          
            int shom = 0;
            while ((line = txtfile.ReadLine()) != null)
            {
                

                        char[] c = line.ToCharArray();
                        int len = c.Length;
                      
                        char[] image_name_char1 = new char[len];

                    image_name[shom] = new string(c);
                shom = shom + 1;
            }
        }

    }
}