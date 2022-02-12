namespace theses
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(Form1));
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea3 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend3 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series3 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Title title3 = new System.Windows.Forms.DataVisualization.Charting.Title();
            this.menuStrip1 = new System.Windows.Forms.MenuStrip();
            this.fileToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.openToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.openToolStripMenuItem2 = new System.Windows.Forms.ToolStripMenuItem();
            this.saveToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.saveAsToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.exiteToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.cameraToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.exportToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.interiorParametersToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.exteriorParametersToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.imageObservationToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.tiePointsCoordinateToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.gCPToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.gPSToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.exteriorAdjusmentToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.bundleAdjustmentToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.visibilityAnalysisToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.optimumImageSelectionToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.helpToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.panel1 = new System.Windows.Forms.Panel();
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.label1 = new System.Windows.Forms.Label();
            this.panel2 = new System.Windows.Forms.Panel();
            this.button1 = new System.Windows.Forms.Button();
            this.label6 = new System.Windows.Forms.Label();
            this.textBox2 = new System.Windows.Forms.TextBox();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.label5 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.panel3 = new System.Windows.Forms.Panel();
            this.dataGridView2 = new System.Windows.Forms.DataGridView();
            this.label3 = new System.Windows.Forms.Label();
            this.chart1 = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.label2 = new System.Windows.Forms.Label();
            this.toolTip1 = new System.Windows.Forms.ToolTip(this.components);
            this.menuStrip1.SuspendLayout();
            this.panel1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            this.panel2.SuspendLayout();
            this.panel3.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView2)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).BeginInit();
            this.SuspendLayout();
            // 
            // menuStrip1
            // 
            this.menuStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.fileToolStripMenuItem,
            this.cameraToolStripMenuItem,
            this.exteriorAdjusmentToolStripMenuItem,
            this.helpToolStripMenuItem});
            this.menuStrip1.Location = new System.Drawing.Point(0, 0);
            this.menuStrip1.Name = "menuStrip1";
            this.menuStrip1.Size = new System.Drawing.Size(905, 24);
            this.menuStrip1.TabIndex = 1;
            this.menuStrip1.Text = "menuStrip1";
            this.menuStrip1.ItemClicked += new System.Windows.Forms.ToolStripItemClickedEventHandler(this.menuStrip1_ItemClicked);
            // 
            // fileToolStripMenuItem
            // 
            this.fileToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.openToolStripMenuItem,
            this.openToolStripMenuItem2,
            this.saveToolStripMenuItem,
            this.saveAsToolStripMenuItem,
            this.exiteToolStripMenuItem});
            this.fileToolStripMenuItem.Name = "fileToolStripMenuItem";
            this.fileToolStripMenuItem.Size = new System.Drawing.Size(37, 20);
            this.fileToolStripMenuItem.Text = "File";
            this.fileToolStripMenuItem.Click += new System.EventHandler(this.fileToolStripMenuItem_Click);
            // 
            // openToolStripMenuItem
            // 
            this.openToolStripMenuItem.Image = ((System.Drawing.Image)(resources.GetObject("openToolStripMenuItem.Image")));
            this.openToolStripMenuItem.Name = "openToolStripMenuItem";
            this.openToolStripMenuItem.Size = new System.Drawing.Size(121, 22);
            this.openToolStripMenuItem.Text = "New File";
            // 
            // openToolStripMenuItem2
            // 
            this.openToolStripMenuItem2.Image = ((System.Drawing.Image)(resources.GetObject("openToolStripMenuItem2.Image")));
            this.openToolStripMenuItem2.Name = "openToolStripMenuItem2";
            this.openToolStripMenuItem2.Size = new System.Drawing.Size(121, 22);
            this.openToolStripMenuItem2.Text = "Open";
            this.openToolStripMenuItem2.Click += new System.EventHandler(this.openToolStripMenuItem2_Click);
            // 
            // saveToolStripMenuItem
            // 
            this.saveToolStripMenuItem.Image = ((System.Drawing.Image)(resources.GetObject("saveToolStripMenuItem.Image")));
            this.saveToolStripMenuItem.Name = "saveToolStripMenuItem";
            this.saveToolStripMenuItem.Size = new System.Drawing.Size(121, 22);
            this.saveToolStripMenuItem.Text = "Save";
            this.saveToolStripMenuItem.Click += new System.EventHandler(this.saveToolStripMenuItem_Click);
            // 
            // saveAsToolStripMenuItem
            // 
            this.saveAsToolStripMenuItem.Image = ((System.Drawing.Image)(resources.GetObject("saveAsToolStripMenuItem.Image")));
            this.saveAsToolStripMenuItem.Name = "saveAsToolStripMenuItem";
            this.saveAsToolStripMenuItem.Size = new System.Drawing.Size(121, 22);
            this.saveAsToolStripMenuItem.Text = "Save as...";
            // 
            // exiteToolStripMenuItem
            // 
            this.exiteToolStripMenuItem.Image = ((System.Drawing.Image)(resources.GetObject("exiteToolStripMenuItem.Image")));
            this.exiteToolStripMenuItem.Name = "exiteToolStripMenuItem";
            this.exiteToolStripMenuItem.Size = new System.Drawing.Size(121, 22);
            this.exiteToolStripMenuItem.Text = "Exit";
            // 
            // cameraToolStripMenuItem
            // 
            this.cameraToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.exportToolStripMenuItem});
            this.cameraToolStripMenuItem.Name = "cameraToolStripMenuItem";
            this.cameraToolStripMenuItem.Size = new System.Drawing.Size(48, 20);
            this.cameraToolStripMenuItem.Text = "Tools";
            this.cameraToolStripMenuItem.Click += new System.EventHandler(this.cameraToolStripMenuItem_Click);
            // 
            // exportToolStripMenuItem
            // 
            this.exportToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.interiorParametersToolStripMenuItem1,
            this.exteriorParametersToolStripMenuItem1,
            this.imageObservationToolStripMenuItem1,
            this.tiePointsCoordinateToolStripMenuItem1,
            this.gCPToolStripMenuItem1,
            this.gPSToolStripMenuItem});
            this.exportToolStripMenuItem.Name = "exportToolStripMenuItem";
            this.exportToolStripMenuItem.Size = new System.Drawing.Size(107, 22);
            this.exportToolStripMenuItem.Text = "Export";
            // 
            // interiorParametersToolStripMenuItem1
            // 
            this.interiorParametersToolStripMenuItem1.Name = "interiorParametersToolStripMenuItem1";
            this.interiorParametersToolStripMenuItem1.Size = new System.Drawing.Size(186, 22);
            this.interiorParametersToolStripMenuItem1.Text = "Interior Parameters";
            this.interiorParametersToolStripMenuItem1.Click += new System.EventHandler(this.interiorParametersToolStripMenuItem1_Click);
            // 
            // exteriorParametersToolStripMenuItem1
            // 
            this.exteriorParametersToolStripMenuItem1.Name = "exteriorParametersToolStripMenuItem1";
            this.exteriorParametersToolStripMenuItem1.Size = new System.Drawing.Size(186, 22);
            this.exteriorParametersToolStripMenuItem1.Text = "Exterior Parameters";
            this.exteriorParametersToolStripMenuItem1.Click += new System.EventHandler(this.exteriorParametersToolStripMenuItem1_Click);
            // 
            // imageObservationToolStripMenuItem1
            // 
            this.imageObservationToolStripMenuItem1.Name = "imageObservationToolStripMenuItem1";
            this.imageObservationToolStripMenuItem1.Size = new System.Drawing.Size(186, 22);
            this.imageObservationToolStripMenuItem1.Text = "Image Observation";
            this.imageObservationToolStripMenuItem1.Click += new System.EventHandler(this.imageObservationToolStripMenuItem1_Click);
            // 
            // tiePointsCoordinateToolStripMenuItem1
            // 
            this.tiePointsCoordinateToolStripMenuItem1.Name = "tiePointsCoordinateToolStripMenuItem1";
            this.tiePointsCoordinateToolStripMenuItem1.Size = new System.Drawing.Size(186, 22);
            this.tiePointsCoordinateToolStripMenuItem1.Text = "Tie points coordinate";
            this.tiePointsCoordinateToolStripMenuItem1.Click += new System.EventHandler(this.tiePointsCoordinateToolStripMenuItem1_Click);
            // 
            // gCPToolStripMenuItem1
            // 
            this.gCPToolStripMenuItem1.Name = "gCPToolStripMenuItem1";
            this.gCPToolStripMenuItem1.Size = new System.Drawing.Size(186, 22);
            this.gCPToolStripMenuItem1.Text = "GCP";
            this.gCPToolStripMenuItem1.Click += new System.EventHandler(this.gCPToolStripMenuItem1_Click);
            // 
            // gPSToolStripMenuItem
            // 
            this.gPSToolStripMenuItem.Name = "gPSToolStripMenuItem";
            this.gPSToolStripMenuItem.Size = new System.Drawing.Size(186, 22);
            this.gPSToolStripMenuItem.Text = "GPS";
            this.gPSToolStripMenuItem.Click += new System.EventHandler(this.gPSToolStripMenuItem_Click);
            // 
            // exteriorAdjusmentToolStripMenuItem
            // 
            this.exteriorAdjusmentToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripMenuItem1,
            this.bundleAdjustmentToolStripMenuItem,
            this.visibilityAnalysisToolStripMenuItem,
            this.optimumImageSelectionToolStripMenuItem});
            this.exteriorAdjusmentToolStripMenuItem.Name = "exteriorAdjusmentToolStripMenuItem";
            this.exteriorAdjusmentToolStripMenuItem.Size = new System.Drawing.Size(70, 20);
            this.exteriorAdjusmentToolStripMenuItem.Text = "Workflow";
            this.exteriorAdjusmentToolStripMenuItem.Click += new System.EventHandler(this.exteriorAdjusmentToolStripMenuItem_Click);
            // 
            // toolStripMenuItem1
            // 
            this.toolStripMenuItem1.Name = "toolStripMenuItem1";
            this.toolStripMenuItem1.Size = new System.Drawing.Size(167, 22);
            this.toolStripMenuItem1.Text = "Import Data";
            this.toolStripMenuItem1.Click += new System.EventHandler(this.toolStripMenuItem1_Click);
            // 
            // bundleAdjustmentToolStripMenuItem
            // 
            this.bundleAdjustmentToolStripMenuItem.Enabled = false;
            this.bundleAdjustmentToolStripMenuItem.Name = "bundleAdjustmentToolStripMenuItem";
            this.bundleAdjustmentToolStripMenuItem.Size = new System.Drawing.Size(167, 22);
            this.bundleAdjustmentToolStripMenuItem.Text = "Error Propagation";
            this.bundleAdjustmentToolStripMenuItem.Click += new System.EventHandler(this.bundleAdjustmentToolStripMenuItem_Click);
            // 
            // visibilityAnalysisToolStripMenuItem
            // 
            this.visibilityAnalysisToolStripMenuItem.Enabled = false;
            this.visibilityAnalysisToolStripMenuItem.Name = "visibilityAnalysisToolStripMenuItem";
            this.visibilityAnalysisToolStripMenuItem.Size = new System.Drawing.Size(167, 22);
            this.visibilityAnalysisToolStripMenuItem.Text = "Visibility Analysis";
            this.visibilityAnalysisToolStripMenuItem.Click += new System.EventHandler(this.visibilityAnalysisToolStripMenuItem_Click);
            // 
            // optimumImageSelectionToolStripMenuItem
            // 
            this.optimumImageSelectionToolStripMenuItem.Name = "optimumImageSelectionToolStripMenuItem";
            this.optimumImageSelectionToolStripMenuItem.Size = new System.Drawing.Size(167, 22);
            this.optimumImageSelectionToolStripMenuItem.Text = "Pre-Analysis";
            this.optimumImageSelectionToolStripMenuItem.Click += new System.EventHandler(this.optimumImageSelectionToolStripMenuItem_Click);
            // 
            // helpToolStripMenuItem
            // 
            this.helpToolStripMenuItem.Name = "helpToolStripMenuItem";
            this.helpToolStripMenuItem.Size = new System.Drawing.Size(44, 20);
            this.helpToolStripMenuItem.Text = "Help";
            // 
            // panel1
            // 
            this.panel1.BackColor = System.Drawing.SystemColors.ControlLightLight;
            this.panel1.Controls.Add(this.dataGridView1);
            this.panel1.Controls.Add(this.label1);
            this.panel1.Location = new System.Drawing.Point(27, 63);
            this.panel1.Name = "panel1";
            this.panel1.Size = new System.Drawing.Size(197, 329);
            this.panel1.TabIndex = 38;
            // 
            // dataGridView1
            // 
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.ImeMode = System.Windows.Forms.ImeMode.NoControl;
            this.dataGridView1.Location = new System.Drawing.Point(18, 26);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.ReadOnly = true;
            this.dataGridView1.Size = new System.Drawing.Size(158, 278);
            this.dataGridView1.TabIndex = 41;
            this.dataGridView1.Visible = false;
            this.dataGridView1.CellContentClick += new System.Windows.Forms.DataGridViewCellEventHandler(this.dataGridView1_CellContentClick);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(3, 0);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(68, 13);
            this.label1.TabIndex = 40;
            this.label1.Text = "Input Images";
            // 
            // panel2
            // 
            this.panel2.BackColor = System.Drawing.SystemColors.ControlLightLight;
            this.panel2.Controls.Add(this.button1);
            this.panel2.Controls.Add(this.label6);
            this.panel2.Controls.Add(this.textBox2);
            this.panel2.Controls.Add(this.textBox1);
            this.panel2.Controls.Add(this.label5);
            this.panel2.Controls.Add(this.label4);
            this.panel2.Controls.Add(this.panel3);
            this.panel2.Controls.Add(this.chart1);
            this.panel2.Controls.Add(this.label2);
            this.panel2.Location = new System.Drawing.Point(323, 63);
            this.panel2.Name = "panel2";
            this.panel2.Size = new System.Drawing.Size(530, 345);
            this.panel2.TabIndex = 42;
            // 
            // button1
            // 
            this.button1.Enabled = false;
            this.button1.Location = new System.Drawing.Point(183, 311);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(67, 23);
            this.button1.TabIndex = 47;
            this.button1.Text = "Browse";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click_1);
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(21, 316);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(88, 13);
            this.label6.TabIndex = 46;
            this.label6.Text = "Select Agisoft file";
            // 
            // textBox2
            // 
            this.textBox2.Enabled = false;
            this.textBox2.Location = new System.Drawing.Point(183, 270);
            this.textBox2.Name = "textBox2";
            this.textBox2.Size = new System.Drawing.Size(67, 20);
            this.textBox2.TabIndex = 45;
            this.textBox2.TextChanged += new System.EventHandler(this.textBox2_TextChanged);
            // 
            // textBox1
            // 
            this.textBox1.Enabled = false;
            this.textBox1.Location = new System.Drawing.Point(183, 232);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(67, 20);
            this.textBox1.TabIndex = 44;
            this.textBox1.TextChanged += new System.EventHandler(this.textBox1_TextChanged);
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(21, 277);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(140, 13);
            this.label5.TabIndex = 43;
            this.label5.Text = "Maximum Number of Images";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(21, 235);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(102, 13);
            this.label4.TabIndex = 42;
            this.label4.Text = "Desired Accuracy %";
            // 
            // panel3
            // 
            this.panel3.BackColor = System.Drawing.SystemColors.ControlLightLight;
            this.panel3.Controls.Add(this.dataGridView2);
            this.panel3.Controls.Add(this.label3);
            this.panel3.Location = new System.Drawing.Point(330, 16);
            this.panel3.Name = "panel3";
            this.panel3.Size = new System.Drawing.Size(197, 313);
            this.panel3.TabIndex = 42;
            // 
            // dataGridView2
            // 
            this.dataGridView2.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView2.ImeMode = System.Windows.Forms.ImeMode.NoControl;
            this.dataGridView2.Location = new System.Drawing.Point(18, 26);
            this.dataGridView2.Name = "dataGridView2";
            this.dataGridView2.ReadOnly = true;
            this.dataGridView2.Size = new System.Drawing.Size(158, 278);
            this.dataGridView2.TabIndex = 41;
            this.dataGridView2.Visible = false;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(3, 0);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(90, 13);
            this.label3.TabIndex = 40;
            this.label3.Text = "Removed Images";
            // 
            // chart1
            // 
            this.chart1.AccessibleName = "";
            chartArea3.Name = "ChartArea1";
            this.chart1.ChartAreas.Add(chartArea3);
            legend3.Enabled = false;
            legend3.Name = "Legend1";
            this.chart1.Legends.Add(legend3);
            this.chart1.Location = new System.Drawing.Point(13, 16);
            this.chart1.Name = "chart1";
            series3.ChartArea = "ChartArea1";
            series3.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Spline;
            series3.LabelBorderWidth = 5;
            series3.Legend = "Legend1";
            series3.Name = "Series1";
            this.chart1.Series.Add(series3);
            this.chart1.Size = new System.Drawing.Size(257, 186);
            this.chart1.TabIndex = 41;
            this.chart1.Tag = "";
            this.chart1.Text = "chart1";
            title3.Name = "Title1";
            title3.Text = "Average of Increasing Errors %";
            this.chart1.Titles.Add(title3);
            this.chart1.Visible = false;
            this.chart1.Click += new System.EventHandler(this.chart1_Click);
            this.chart1.MouseMove += new System.Windows.Forms.MouseEventHandler(this.chart1_MouseMove);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(237, 0);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(71, 13);
            this.label2.TabIndex = 40;
            this.label2.Text = "Data Analysis";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(905, 420);
            this.Controls.Add(this.panel2);
            this.Controls.Add(this.panel1);
            this.Controls.Add(this.menuStrip1);
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.Name = "Form1";
            this.Text = "UAV Image Selection";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.menuStrip1.ResumeLayout(false);
            this.menuStrip1.PerformLayout();
            this.panel1.ResumeLayout(false);
            this.panel1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.panel2.ResumeLayout(false);
            this.panel2.PerformLayout();
            this.panel3.ResumeLayout(false);
            this.panel3.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView2)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.MenuStrip menuStrip1;
        private System.Windows.Forms.ToolStripMenuItem fileToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem openToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem openToolStripMenuItem2;
        private System.Windows.Forms.ToolStripMenuItem saveToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem saveAsToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem exiteToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem cameraToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem exteriorAdjusmentToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem exportToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem interiorParametersToolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem exteriorParametersToolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem imageObservationToolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem tiePointsCoordinateToolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem gCPToolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem gPSToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem helpToolStripMenuItem;
        public System.Windows.Forms.ToolStripMenuItem visibilityAnalysisToolStripMenuItem;
        public System.Windows.Forms.ToolStripMenuItem optimumImageSelectionToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem1;
        public System.Windows.Forms.ToolStripMenuItem bundleAdjustmentToolStripMenuItem;
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.Label label1;
        public System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.Panel panel2;
        private System.Windows.Forms.Label label2;
        public System.Windows.Forms.DataGridView dataGridView2;
        private System.Windows.Forms.Label label3;
        public System.Windows.Forms.DataVisualization.Charting.Chart chart1;
        public System.Windows.Forms.ToolTip toolTip1;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label4;
        public System.Windows.Forms.Panel panel3;
        public System.Windows.Forms.TextBox textBox2;
        public System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.Label label6;
        public System.Windows.Forms.Button button1;
    }
}

