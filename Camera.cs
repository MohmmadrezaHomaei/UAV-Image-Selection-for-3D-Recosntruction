using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace theses
{
    public partial class Camera : Form
    {
        public Camera()
        {
            InitializeComponent();
        }

        private void Camera_Load(object sender, EventArgs e)
        {
            DataSet dataset = new DataSet();
            dataset = Form1.ds;
        }
    }
}
